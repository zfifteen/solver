#include "linsolve/poisson_solver.hpp"
#include "solver/projection.hpp"

#include "operators/discrete_operators.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace solver {

namespace {

constexpr double kSmallNumber = 1.0e-30;

BoundaryFace lower_face(const Axis axis) {
  switch(axis) {
    case Axis::x:
      return BoundaryFace::x_min;
    case Axis::y:
      return BoundaryFace::y_min;
    case Axis::z:
      return BoundaryFace::z_min;
  }

  return BoundaryFace::x_min;
}

BoundaryFace upper_face(const Axis axis) {
  switch(axis) {
    case Axis::x:
      return BoundaryFace::x_max;
    case Axis::y:
      return BoundaryFace::y_max;
    case Axis::z:
      return BoundaryFace::z_max;
  }

  return BoundaryFace::x_max;
}

bool is_dirichlet_velocity_boundary(const BoundaryCondition& condition) {
  return condition.type == PhysicalBoundaryType::no_slip_wall ||
         condition.type == PhysicalBoundaryType::prescribed_velocity ||
         condition.type == PhysicalBoundaryType::symmetry;
}

struct GhostRelation {
  double coefficient = 0.0;
  double constant = 0.0;
};

double component_value(const BoundaryCondition& condition, const Axis axis) {
  return condition.velocity[static_cast<std::size_t>(axis_index(axis))];
}

GhostRelation tangential_ghost_relation(const BoundaryCondition& condition,
                                        const Axis component_axis) {
  switch(condition.type) {
    case PhysicalBoundaryType::no_slip_wall:
    case PhysicalBoundaryType::prescribed_velocity:
      return GhostRelation{
          .coefficient = -1.0,
          .constant = 2.0 * component_value(condition, component_axis),
      };
    case PhysicalBoundaryType::symmetry:
    case PhysicalBoundaryType::fixed_pressure:
      return GhostRelation{
          .coefficient = 1.0,
          .constant = 0.0,
      };
    case PhysicalBoundaryType::periodic:
      return GhostRelation{};
  }

  return GhostRelation{};
}

void require_same_layout(const StructuredField& left,
                         const StructuredField& right,
                         const char* operation_name) {
  if(!left.layout().same_shape_as(right.layout())) {
    throw std::invalid_argument(std::string(operation_name) + ": incompatible field layouts");
  }
}

void copy_range(StructuredField& destination,
                const IndexRange3D& destination_range,
                const StructuredField& source,
                const IndexRange3D& source_range) {
  const Extent3D destination_extent = destination_range.extent();
  const Extent3D source_extent = source_range.extent();
  if(destination_extent.nx != source_extent.nx || destination_extent.ny != source_extent.ny ||
     destination_extent.nz != source_extent.nz) {
    throw std::invalid_argument("copy_range expects equal extents");
  }

  for(int k = 0; k < destination_extent.nz; ++k) {
    for(int j = 0; j < destination_extent.ny; ++j) {
      for(int i = 0; i < destination_extent.nx; ++i) {
        destination(destination_range.i_begin + i, destination_range.j_begin + j,
                    destination_range.k_begin + k) =
            source(source_range.i_begin + i, source_range.j_begin + j, source_range.k_begin + k);
      }
    }
  }
}

void copy_active(const StructuredField& source, StructuredField& destination) {
  require_same_layout(source, destination, "copy_active");
  copy_range(destination, destination.layout().active_range(), source, source.layout().active_range());
}

double active_mean(const StructuredField& field) {
  const IndexRange3D active = field.layout().active_range();
  const std::size_t count = active.extent().cell_count();
  double sum = 0.0;

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        sum += field(i, j, k);
      }
    }
  }

  return sum / static_cast<double>(count);
}

double active_l2_norm(const StructuredField& field) {
  const IndexRange3D active = field.layout().active_range();
  const std::size_t count = active.extent().cell_count();
  double sum = 0.0;

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        sum += field(i, j, k) * field(i, j, k);
      }
    }
  }

  return std::sqrt(sum / static_cast<double>(count));
}

void axpy_active(StructuredField& destination, const StructuredField& source, const double scale) {
  require_same_layout(destination, source, "axpy_active");
  const IndexRange3D active = destination.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        destination(i, j, k) += scale * source(i, j, k);
      }
    }
  }
}

void copy_component(const StructuredField& source, StructuredField& destination) {
  require_same_layout(source, destination, "copy_component");

  const Extent3D storage = source.layout().storage_extent();
  for(int k = 0; k < storage.nz; ++k) {
    for(int j = 0; j < storage.ny; ++j) {
      for(int i = 0; i < storage.nx; ++i) {
        destination(i, j, k) = source(i, j, k);
      }
    }
  }
}

void scale_active(StructuredField& field, const double scale) {
  const IndexRange3D active = field.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        field(i, j, k) *= scale;
      }
    }
  }
}

bool pressure_axis_is_periodic(const PressureBoundarySet& conditions, const Axis axis) {
  const PressureBoundaryCondition& lower = conditions[lower_face(axis)];
  const PressureBoundaryCondition& upper = conditions[upper_face(axis)];
  if((lower.type == PressureBoundaryType::periodic) !=
     (upper.type == PressureBoundaryType::periodic)) {
    throw std::invalid_argument("periodic pressure boundaries must be paired on both sides");
  }
  return lower.type == PressureBoundaryType::periodic;
}

bool velocity_axis_is_periodic(const BoundaryConditionSet& conditions, const Axis axis) {
  const BoundaryCondition& lower = conditions[lower_face(axis)];
  const BoundaryCondition& upper = conditions[upper_face(axis)];
  if((lower.type == PhysicalBoundaryType::periodic) !=
     (upper.type == PhysicalBoundaryType::periodic)) {
    throw std::invalid_argument("periodic velocity boundaries must be paired on both sides");
  }
  return lower.type == PhysicalBoundaryType::periodic;
}

void apply_periodic_pair(StructuredField& field, const Axis axis, const bool duplicated_boundary_axis) {
  const IndexRange3D active = field.layout().active_range();
  const IndexRange3D lower_ghost = field.layout().ghost_range(lower_face(axis));
  const IndexRange3D upper_ghost = field.layout().ghost_range(upper_face(axis));
  copy_range(field, lower_ghost, field,
             duplicated_boundary_axis ? field.layout().boundary_active_range(upper_face(axis))
                                      : IndexRange3D{
                                            .i_begin = axis == Axis::x ? active.i_end - 1 : active.i_begin,
                                            .i_end = axis == Axis::x ? active.i_end : active.i_end,
                                            .j_begin = axis == Axis::y ? active.j_end - 1 : active.j_begin,
                                            .j_end = axis == Axis::y ? active.j_end : active.j_end,
                                            .k_begin = axis == Axis::z ? active.k_end - 1 : active.k_begin,
                                            .k_end = axis == Axis::z ? active.k_end : active.k_end,
                                        });
  copy_range(field, upper_ghost, field,
             duplicated_boundary_axis ? field.layout().boundary_active_range(lower_face(axis))
                                      : IndexRange3D{
                                            .i_begin = axis == Axis::x ? active.i_begin : active.i_begin,
                                            .i_end = axis == Axis::x ? active.i_begin + 1 : active.i_end,
                                            .j_begin = axis == Axis::y ? active.j_begin : active.j_begin,
                                            .j_end = axis == Axis::y ? active.j_begin + 1 : active.j_end,
                                            .k_begin = axis == Axis::z ? active.k_begin : active.k_begin,
                                            .k_end = axis == Axis::z ? active.k_begin + 1 : active.k_end,
                                        });
}

void apply_periodic_pair(FaceField& field, const Axis axis) {
  if(field.normal_axis() == axis) {
    const IndexRange3D lower_boundary = field.layout().boundary_active_range(lower_face(axis));
    const IndexRange3D upper_boundary = field.layout().boundary_active_range(upper_face(axis));
    const Extent3D extent = lower_boundary.extent();

    for(int k = 0; k < extent.nz; ++k) {
      for(int j = 0; j < extent.ny; ++j) {
        for(int i = 0; i < extent.nx; ++i) {
          const double shared = 0.5 * (field(lower_boundary.i_begin + i, lower_boundary.j_begin + j,
                                             lower_boundary.k_begin + k) +
                                       field(upper_boundary.i_begin + i, upper_boundary.j_begin + j,
                                             upper_boundary.k_begin + k));
          field(lower_boundary.i_begin + i, lower_boundary.j_begin + j, lower_boundary.k_begin + k) =
              shared;
          field(upper_boundary.i_begin + i, upper_boundary.j_begin + j, upper_boundary.k_begin + k) =
              shared;
        }
      }
    }

    const IndexRange3D lower_ghost = field.layout().ghost_range(lower_face(axis));
    const IndexRange3D upper_ghost = field.layout().ghost_range(upper_face(axis));
    const IndexRange3D source_lower = axis == Axis::x ? IndexRange3D{field.layout().active_range().i_end - 2,
                                                                     field.layout().active_range().i_end - 1,
                                                                     field.layout().active_range().j_begin,
                                                                     field.layout().active_range().j_end,
                                                                     field.layout().active_range().k_begin,
                                                                     field.layout().active_range().k_end}
                                                      : axis == Axis::y
                                                            ? IndexRange3D{field.layout().active_range().i_begin,
                                                                           field.layout().active_range().i_end,
                                                                           field.layout().active_range().j_end - 2,
                                                                           field.layout().active_range().j_end - 1,
                                                                           field.layout().active_range().k_begin,
                                                                           field.layout().active_range().k_end}
                                                            : IndexRange3D{field.layout().active_range().i_begin,
                                                                           field.layout().active_range().i_end,
                                                                           field.layout().active_range().j_begin,
                                                                           field.layout().active_range().j_end,
                                                                           field.layout().active_range().k_end - 2,
                                                                           field.layout().active_range().k_end - 1};
    const IndexRange3D source_upper = axis == Axis::x ? IndexRange3D{field.layout().active_range().i_begin + 1,
                                                                     field.layout().active_range().i_begin + 2,
                                                                     field.layout().active_range().j_begin,
                                                                     field.layout().active_range().j_end,
                                                                     field.layout().active_range().k_begin,
                                                                     field.layout().active_range().k_end}
                                                      : axis == Axis::y
                                                            ? IndexRange3D{field.layout().active_range().i_begin,
                                                                           field.layout().active_range().i_end,
                                                                           field.layout().active_range().j_begin + 1,
                                                                           field.layout().active_range().j_begin + 2,
                                                                           field.layout().active_range().k_begin,
                                                                           field.layout().active_range().k_end}
                                                            : IndexRange3D{field.layout().active_range().i_begin,
                                                                           field.layout().active_range().i_end,
                                                                           field.layout().active_range().j_begin,
                                                                           field.layout().active_range().j_end,
                                                                           field.layout().active_range().k_begin + 1,
                                                                           field.layout().active_range().k_begin + 2};
    copy_range(field, lower_ghost, field, source_lower);
    copy_range(field, upper_ghost, field, source_upper);
    return;
  }

  apply_periodic_pair(static_cast<StructuredField&>(field), axis, false);
}

void apply_pressure_face_boundary(PressureField& pressure,
                                  const BoundaryFace face,
                                  const PressureBoundaryCondition& condition) {
  const Axis axis = boundary_axis(face);
  const IndexRange3D boundary_active = pressure.layout().boundary_active_range(face);
  const IndexRange3D ghost_range = pressure.layout().ghost_range(face);
  const double spacing = pressure.layout().grid().spacing(axis);
  const Extent3D extent = boundary_active.extent();

  for(int k = 0; k < extent.nz; ++k) {
    for(int j = 0; j < extent.ny; ++j) {
      for(int i = 0; i < extent.nx; ++i) {
        const int bi = boundary_active.i_begin + i;
        const int bj = boundary_active.j_begin + j;
        const int bk = boundary_active.k_begin + k;
        const int gi = ghost_range.i_begin + i;
        const int gj = ghost_range.j_begin + j;
        const int gk = ghost_range.k_begin + k;
        const double active_value = pressure(bi, bj, bk);

        switch(condition.type) {
          case PressureBoundaryType::neumann:
            pressure(gi, gj, gk) =
                is_lower_boundary(face) ? active_value - condition.gradient * spacing
                                        : active_value + condition.gradient * spacing;
            break;
          case PressureBoundaryType::dirichlet:
            pressure(gi, gj, gk) = 2.0 * condition.value - active_value;
            break;
          case PressureBoundaryType::periodic:
            break;
        }
      }
    }
  }
}

void apply_total_pressure_gradient_face(PressureField& pressure_total,
                                        const BoundaryFace face,
                                        const FaceField& wall_normal_source) {
  const Axis axis = boundary_axis(face);
  const IndexRange3D boundary_active = pressure_total.layout().boundary_active_range(face);
  const IndexRange3D ghost_range = pressure_total.layout().ghost_range(face);
  const IndexRange3D wall_active = wall_normal_source.layout().active_range();
  const double spacing = pressure_total.layout().grid().spacing(axis);
  const bool lower = is_lower_boundary(face);
  const Extent3D extent = boundary_active.extent();

  for(int k = 0; k < extent.nz; ++k) {
    for(int j = 0; j < extent.ny; ++j) {
      for(int i = 0; i < extent.nx; ++i) {
        const int bi = boundary_active.i_begin + i;
        const int bj = boundary_active.j_begin + j;
        const int bk = boundary_active.k_begin + k;
        const int gi = ghost_range.i_begin + i;
        const int gj = ghost_range.j_begin + j;
        const int gk = ghost_range.k_begin + k;

        int wi = bi;
        int wj = bj;
        int wk = bk;
        switch(axis) {
          case Axis::x:
            wi = lower ? wall_active.i_begin : wall_active.i_end - 1;
            break;
          case Axis::y:
            wj = lower ? wall_active.j_begin : wall_active.j_end - 1;
            break;
          case Axis::z:
            wk = lower ? wall_active.k_begin : wall_active.k_end - 1;
            break;
        }

        const double gradient = wall_normal_source(wi, wj, wk);
        const double active_value = pressure_total(bi, bj, bk);
        pressure_total(gi, gj, gk) =
            lower ? active_value - gradient * spacing : active_value + gradient * spacing;
      }
    }
  }
}

void apply_face_field_boundary(FaceField& field,
                               const BoundaryFace face,
                               const BoundaryCondition& condition) {
  if(condition.type == PhysicalBoundaryType::periodic) {
    return;
  }

  const Axis boundary_axis_value = boundary_axis(face);
  const Axis component_axis = field.normal_axis();
  const IndexRange3D boundary_active = field.layout().boundary_active_range(face);
  const IndexRange3D ghost_range = field.layout().ghost_range(face);
  const Extent3D extent = boundary_active.extent();

  for(int k = 0; k < extent.nz; ++k) {
    for(int j = 0; j < extent.ny; ++j) {
      for(int i = 0; i < extent.nx; ++i) {
        const int bi = boundary_active.i_begin + i;
        const int bj = boundary_active.j_begin + j;
        const int bk = boundary_active.k_begin + k;
        const int gi = ghost_range.i_begin + i;
        const int gj = ghost_range.j_begin + j;
        const int gk = ghost_range.k_begin + k;

        if(component_axis == boundary_axis_value) {
          if(condition.type == PhysicalBoundaryType::fixed_pressure) {
            const int ii = is_lower_boundary(face) ? bi + (boundary_axis_value == Axis::x)
                                                   : bi - (boundary_axis_value == Axis::x);
            const int jj = is_lower_boundary(face) ? bj + (boundary_axis_value == Axis::y)
                                                   : bj - (boundary_axis_value == Axis::y);
            const int kk = is_lower_boundary(face) ? bk + (boundary_axis_value == Axis::z)
                                                   : bk - (boundary_axis_value == Axis::z);
            field(bi, bj, bk) = field(ii, jj, kk);
            field(gi, gj, gk) = field(bi, bj, bk);
            continue;
          }

          const double boundary_value =
              condition.type == PhysicalBoundaryType::symmetry ? 0.0
                                                               : component_value(condition, component_axis);
          const int ii = is_lower_boundary(face) ? bi + (boundary_axis_value == Axis::x)
                                                 : bi - (boundary_axis_value == Axis::x);
          const int jj = is_lower_boundary(face) ? bj + (boundary_axis_value == Axis::y)
                                                 : bj - (boundary_axis_value == Axis::y);
          const int kk = is_lower_boundary(face) ? bk + (boundary_axis_value == Axis::z)
                                                 : bk - (boundary_axis_value == Axis::z);

          field(bi, bj, bk) = boundary_value;
          field(gi, gj, gk) = 2.0 * boundary_value - field(ii, jj, kk);
          continue;
        }

        const int ii = bi;
        const int jj = bj;
        const int kk = bk;
        if(condition.type == PhysicalBoundaryType::no_slip_wall ||
           condition.type == PhysicalBoundaryType::prescribed_velocity) {
          const double wall_value = component_value(condition, component_axis);
          field(gi, gj, gk) = 2.0 * wall_value - field(ii, jj, kk);
        } else if(condition.type == PhysicalBoundaryType::symmetry ||
                  condition.type == PhysicalBoundaryType::fixed_pressure) {
          field(gi, gj, gk) = field(ii, jj, kk);
        }
      }
    }
  }
}

void apply_velocity_component_boundaries(const BoundaryConditionSet& boundary_conditions,
                                         FaceField& field) {
  for(const Axis axis : {Axis::x, Axis::y, Axis::z}) {
    if(velocity_axis_is_periodic(boundary_conditions, axis)) {
      apply_periodic_pair(field, axis);
      continue;
    }

    apply_face_field_boundary(field, lower_face(axis), boundary_conditions[lower_face(axis)]);
    apply_face_field_boundary(field, upper_face(axis), boundary_conditions[upper_face(axis)]);
  }
}

int sweep_begin_for_axis(const FaceField& field,
                         const Axis axis,
                         const BoundaryConditionSet& boundary_conditions) {
  if(field.normal_axis() != axis) {
    const IndexRange3D active = field.layout().active_range();
    switch(axis) {
      case Axis::x:
        return active.i_begin;
      case Axis::y:
        return active.j_begin;
      case Axis::z:
        return active.k_begin;
    }
  }

  const IndexRange3D active = field.layout().active_range();
  const bool lower_fixed = is_dirichlet_velocity_boundary(boundary_conditions[lower_face(axis)]);
  switch(axis) {
    case Axis::x:
      return active.i_begin + (lower_fixed ? 1 : 0);
    case Axis::y:
      return active.j_begin + (lower_fixed ? 1 : 0);
    case Axis::z:
      return active.k_begin + (lower_fixed ? 1 : 0);
  }

  return 0;
}

int sweep_end_for_axis(const FaceField& field,
                       const Axis axis,
                       const BoundaryConditionSet& boundary_conditions) {
  if(field.normal_axis() != axis) {
    const IndexRange3D active = field.layout().active_range();
    switch(axis) {
      case Axis::x:
        return active.i_end;
      case Axis::y:
        return active.j_end;
      case Axis::z:
        return active.k_end;
    }
  }

  const IndexRange3D active = field.layout().active_range();
  const bool upper_fixed = is_dirichlet_velocity_boundary(boundary_conditions[upper_face(axis)]);
  switch(axis) {
    case Axis::x:
      return active.i_end - (upper_fixed ? 1 : 0);
    case Axis::y:
      return active.j_end - (upper_fixed ? 1 : 0);
    case Axis::z:
      return active.k_end - (upper_fixed ? 1 : 0);
  }

  return 0;
}

void solve_tridiagonal(const std::vector<double>& lower,
                       const std::vector<double>& diagonal,
                       const std::vector<double>& upper,
                       const std::vector<double>& rhs,
                       std::vector<double>& solution) {
  const std::size_t n = diagonal.size();
  if(n == 0) {
    solution.clear();
    return;
  }

  std::vector<double> modified_upper(n, 0.0);
  std::vector<double> modified_rhs(n, 0.0);
  solution.assign(n, 0.0);

  double pivot = diagonal[0];
  if(std::abs(pivot) < kSmallNumber) {
    throw std::runtime_error("tridiagonal solve encountered a zero pivot");
  }

  modified_upper[0] = n > 1 ? upper[0] / pivot : 0.0;
  modified_rhs[0] = rhs[0] / pivot;

  for(std::size_t i = 1; i < n; ++i) {
    pivot = diagonal[i] - lower[i] * modified_upper[i - 1];
    if(std::abs(pivot) < kSmallNumber) {
      throw std::runtime_error("tridiagonal solve encountered a zero pivot");
    }

    modified_upper[i] = i + 1 < n ? upper[i] / pivot : 0.0;
    modified_rhs[i] = (rhs[i] - lower[i] * modified_rhs[i - 1]) / pivot;
  }

  solution[n - 1] = modified_rhs[n - 1];
  for(std::size_t i = n - 1; i-- > 0;) {
    solution[i] = modified_rhs[i] - modified_upper[i] * solution[i + 1];
  }
}

void solve_cyclic_tridiagonal(const std::vector<double>& lower,
                              const std::vector<double>& diagonal,
                              const std::vector<double>& upper,
                              const double lower_corner,
                              const double upper_corner,
                              const std::vector<double>& rhs,
                              std::vector<double>& solution) {
  const std::size_t n = diagonal.size();
  if(n == 0) {
    solution.clear();
    return;
  }

  if(n == 1) {
    const double pivot = diagonal[0] + lower_corner + upper_corner;
    if(std::abs(pivot) < kSmallNumber) {
      throw std::runtime_error("cyclic tridiagonal solve encountered a zero pivot");
    }
    solution.assign(1, rhs[0] / pivot);
    return;
  }

  std::vector<double> diagonal_modified = diagonal;
  const double gamma = -diagonal[0];
  if(std::abs(gamma) < kSmallNumber) {
    throw std::runtime_error("cyclic tridiagonal solve encountered a zero gamma");
  }

  diagonal_modified[0] = diagonal[0] - gamma;
  diagonal_modified[n - 1] = diagonal[n - 1] - lower_corner * upper_corner / gamma;

  std::vector<double> lower_modified = lower;
  std::vector<double> upper_modified = upper;
  lower_modified[0] = 0.0;
  upper_modified[n - 1] = 0.0;

  std::vector<double> x;
  solve_tridiagonal(lower_modified, diagonal_modified, upper_modified, rhs, x);

  std::vector<double> u(n, 0.0);
  u[0] = gamma;
  u[n - 1] = lower_corner;

  std::vector<double> z;
  solve_tridiagonal(lower_modified, diagonal_modified, upper_modified, u, z);

  const double denominator = 1.0 + z[0] + upper_corner * z[n - 1] / gamma;
  if(std::abs(denominator) < kSmallNumber) {
    throw std::runtime_error("cyclic tridiagonal solve encountered a singular Sherman-Morrison update");
  }

  const double factor = (x[0] + upper_corner * x[n - 1] / gamma) / denominator;
  solution.resize(n);
  for(std::size_t i = 0; i < n; ++i) {
    solution[i] = x[i] - factor * z[i];
  }
}

void solve_component_sweep_periodic(const FaceField& input,
                                    const Axis axis,
                                    const double alpha,
                                    FaceField& output,
                                    int& line_solves) {
  const Grid& grid = input.layout().grid();
  const double spacing = grid.spacing(axis);
  const double beta = alpha / (spacing * spacing);
  const IndexRange3D active = input.layout().active_range();

  const int begin = [&]() {
    switch(axis) {
      case Axis::x:
        return active.i_begin;
      case Axis::y:
        return active.j_begin;
      case Axis::z:
        return active.k_begin;
    }
    return 0;
  }();
  const int end = [&]() {
    const int active_end = [&]() {
      switch(axis) {
        case Axis::x:
          return active.i_end;
        case Axis::y:
          return active.j_end;
        case Axis::z:
          return active.k_end;
      }
      return 0;
    }();
    return input.normal_axis() == axis ? active_end - 1 : active_end;
  }();
  const int unknowns = end - begin;
  if(unknowns <= 0) {
    return;
  }

  std::vector<double> lower(static_cast<std::size_t>(unknowns), -beta);
  std::vector<double> diagonal(static_cast<std::size_t>(unknowns), 1.0 + 2.0 * beta);
  std::vector<double> upper(static_cast<std::size_t>(unknowns), -beta);
  std::vector<double> rhs(static_cast<std::size_t>(unknowns), 0.0);
  std::vector<double> solution;
  lower[0] = 0.0;
  upper[static_cast<std::size_t>(unknowns - 1)] = 0.0;

  switch(axis) {
    case Axis::x:
      for(int k = active.k_begin; k < active.k_end; ++k) {
        for(int j = active.j_begin; j < active.j_end; ++j) {
          for(int i = begin; i < end; ++i) {
            rhs[static_cast<std::size_t>(i - begin)] = input(i, j, k);
          }
          solve_cyclic_tridiagonal(lower, diagonal, upper, -beta, -beta, rhs, solution);
          for(int i = begin; i < end; ++i) {
            output(i, j, k) = solution[static_cast<std::size_t>(i - begin)];
          }
          if(input.normal_axis() == axis) {
            output(active.i_end - 1, j, k) = output(begin, j, k);
          }
          ++line_solves;
        }
      }
      return;
    case Axis::y:
      for(int k = active.k_begin; k < active.k_end; ++k) {
        for(int i = active.i_begin; i < active.i_end; ++i) {
          for(int j = begin; j < end; ++j) {
            rhs[static_cast<std::size_t>(j - begin)] = input(i, j, k);
          }
          solve_cyclic_tridiagonal(lower, diagonal, upper, -beta, -beta, rhs, solution);
          for(int j = begin; j < end; ++j) {
            output(i, j, k) = solution[static_cast<std::size_t>(j - begin)];
          }
          if(input.normal_axis() == axis) {
            output(i, active.j_end - 1, k) = output(i, begin, k);
          }
          ++line_solves;
        }
      }
      return;
    case Axis::z:
      for(int j = active.j_begin; j < active.j_end; ++j) {
        for(int i = active.i_begin; i < active.i_end; ++i) {
          for(int k = begin; k < end; ++k) {
            rhs[static_cast<std::size_t>(k - begin)] = input(i, j, k);
          }
          solve_cyclic_tridiagonal(lower, diagonal, upper, -beta, -beta, rhs, solution);
          for(int k = begin; k < end; ++k) {
            output(i, j, k) = solution[static_cast<std::size_t>(k - begin)];
          }
          if(input.normal_axis() == axis) {
            output(i, j, active.k_end - 1) = output(i, j, begin);
          }
          ++line_solves;
        }
      }
      return;
  }
}

void solve_component_sweep(const FaceField& input,
                           const Axis axis,
                           const BoundaryConditionSet& boundary_conditions,
                           const double alpha,
                           FaceField& output,
                           int& line_solves) {
  if(velocity_axis_is_periodic(boundary_conditions, axis)) {
    solve_component_sweep_periodic(input, axis, alpha, output, line_solves);
    return;
  }

  const Grid& grid = input.layout().grid();
  const double spacing = grid.spacing(axis);
  const double beta = alpha / (spacing * spacing);
  const IndexRange3D active = input.layout().active_range();
  const int begin = sweep_begin_for_axis(input, axis, boundary_conditions);
  const int end = sweep_end_for_axis(input, axis, boundary_conditions);
  const int unknowns = end - begin;

  if(unknowns <= 0) {
    return;
  }

  std::vector<double> lower(static_cast<std::size_t>(unknowns), -beta);
  std::vector<double> diagonal(static_cast<std::size_t>(unknowns), 1.0 + 2.0 * beta);
  std::vector<double> upper(static_cast<std::size_t>(unknowns), -beta);
  std::vector<double> rhs(static_cast<std::size_t>(unknowns), 0.0);
  std::vector<double> solution;
  lower[0] = 0.0;
  upper[static_cast<std::size_t>(unknowns - 1)] = 0.0;

  switch(axis) {
    case Axis::x:
      for(int k = active.k_begin; k < active.k_end; ++k) {
        for(int j = active.j_begin; j < active.j_end; ++j) {
          std::fill(diagonal.begin(), diagonal.end(), 1.0 + 2.0 * beta);
          for(int i = begin; i < end; ++i) {
            const int local = i - begin;
            rhs[static_cast<std::size_t>(local)] = input(i, j, k);
          }

          if(input.normal_axis() == axis && is_dirichlet_velocity_boundary(boundary_conditions[lower_face(axis)])) {
            rhs[0] += beta * component_value(boundary_conditions[lower_face(axis)], axis);
          } else {
            const GhostRelation relation =
                tangential_ghost_relation(boundary_conditions[lower_face(axis)], input.normal_axis());
            diagonal[0] -= beta * relation.coefficient;
            rhs[0] += beta * relation.constant;
          }

          if(input.normal_axis() == axis && is_dirichlet_velocity_boundary(boundary_conditions[upper_face(axis)])) {
            rhs[static_cast<std::size_t>(unknowns - 1)] +=
                beta * component_value(boundary_conditions[upper_face(axis)], axis);
          } else {
            const GhostRelation relation =
                tangential_ghost_relation(boundary_conditions[upper_face(axis)], input.normal_axis());
            diagonal[static_cast<std::size_t>(unknowns - 1)] -= beta * relation.coefficient;
            rhs[static_cast<std::size_t>(unknowns - 1)] += beta * relation.constant;
          }

          solve_tridiagonal(lower, diagonal, upper, rhs, solution);

          for(int i = begin; i < end; ++i) {
            output(i, j, k) = solution[static_cast<std::size_t>(i - begin)];
          }
          ++line_solves;
        }
      }
      return;
    case Axis::y:
      for(int k = active.k_begin; k < active.k_end; ++k) {
        for(int i = active.i_begin; i < active.i_end; ++i) {
          std::fill(diagonal.begin(), diagonal.end(), 1.0 + 2.0 * beta);
          for(int j = begin; j < end; ++j) {
            const int local = j - begin;
            rhs[static_cast<std::size_t>(local)] = input(i, j, k);
          }

          if(input.normal_axis() == axis && is_dirichlet_velocity_boundary(boundary_conditions[lower_face(axis)])) {
            rhs[0] += beta * component_value(boundary_conditions[lower_face(axis)], axis);
          } else {
            const GhostRelation relation =
                tangential_ghost_relation(boundary_conditions[lower_face(axis)], input.normal_axis());
            diagonal[0] -= beta * relation.coefficient;
            rhs[0] += beta * relation.constant;
          }

          if(input.normal_axis() == axis && is_dirichlet_velocity_boundary(boundary_conditions[upper_face(axis)])) {
            rhs[static_cast<std::size_t>(unknowns - 1)] +=
                beta * component_value(boundary_conditions[upper_face(axis)], axis);
          } else {
            const GhostRelation relation =
                tangential_ghost_relation(boundary_conditions[upper_face(axis)], input.normal_axis());
            diagonal[static_cast<std::size_t>(unknowns - 1)] -= beta * relation.coefficient;
            rhs[static_cast<std::size_t>(unknowns - 1)] += beta * relation.constant;
          }

          solve_tridiagonal(lower, diagonal, upper, rhs, solution);

          for(int j = begin; j < end; ++j) {
            output(i, j, k) = solution[static_cast<std::size_t>(j - begin)];
          }
          ++line_solves;
        }
      }
      return;
    case Axis::z:
      for(int j = active.j_begin; j < active.j_end; ++j) {
        for(int i = active.i_begin; i < active.i_end; ++i) {
          std::fill(diagonal.begin(), diagonal.end(), 1.0 + 2.0 * beta);
          for(int k = begin; k < end; ++k) {
            const int local = k - begin;
            rhs[static_cast<std::size_t>(local)] = input(i, j, k);
          }

          if(input.normal_axis() == axis && is_dirichlet_velocity_boundary(boundary_conditions[lower_face(axis)])) {
            rhs[0] += beta * component_value(boundary_conditions[lower_face(axis)], axis);
          } else {
            const GhostRelation relation =
                tangential_ghost_relation(boundary_conditions[lower_face(axis)], input.normal_axis());
            diagonal[0] -= beta * relation.coefficient;
            rhs[0] += beta * relation.constant;
          }

          if(input.normal_axis() == axis && is_dirichlet_velocity_boundary(boundary_conditions[upper_face(axis)])) {
            rhs[static_cast<std::size_t>(unknowns - 1)] +=
                beta * component_value(boundary_conditions[upper_face(axis)], axis);
          } else {
            const GhostRelation relation =
                tangential_ghost_relation(boundary_conditions[upper_face(axis)], input.normal_axis());
            diagonal[static_cast<std::size_t>(unknowns - 1)] -= beta * relation.coefficient;
            rhs[static_cast<std::size_t>(unknowns - 1)] += beta * relation.constant;
          }

          solve_tridiagonal(lower, diagonal, upper, rhs, solution);

          for(int k = begin; k < end; ++k) {
            output(i, j, k) = solution[static_cast<std::size_t>(k - begin)];
          }
          ++line_solves;
        }
      }
      return;
  }
}

void solve_component_adi(FaceField& field,
                         const BoundaryConditionSet& boundary_conditions,
                         const double alpha,
                         int& line_solves) {
  apply_velocity_component_boundaries(boundary_conditions, field);
  if(alpha == 0.0) {
    return;
  }

  FaceField working = field;
  FaceField next(field.normal_axis(), field.layout().grid());

  for(const Axis axis : {Axis::x, Axis::y, Axis::z}) {
    copy_component(working, next);
    solve_component_sweep(working, axis, boundary_conditions, alpha, next, line_solves);
    apply_velocity_component_boundaries(boundary_conditions, next);
    working = next;
  }

  field = working;
}

}  // namespace

BoundaryConditionSet BoundaryConditionSet::all(const PhysicalBoundaryType type) {
  BoundaryConditionSet conditions;
  for(BoundaryCondition& condition : conditions.faces) {
    condition.type = type;
  }
  return conditions;
}

BoundaryConditionSet BoundaryConditionSet::cavity() {
  return all(PhysicalBoundaryType::no_slip_wall);
}

std::string to_string(const PhysicalBoundaryType type) {
  switch(type) {
    case PhysicalBoundaryType::no_slip_wall:
      return "no_slip_wall";
    case PhysicalBoundaryType::prescribed_velocity:
      return "prescribed_velocity";
    case PhysicalBoundaryType::symmetry:
      return "symmetry";
    case PhysicalBoundaryType::fixed_pressure:
      return "fixed_pressure";
    case PhysicalBoundaryType::periodic:
      return "periodic";
  }

  return "unknown";
}

std::string to_string(const PressureBoundaryType type) {
  switch(type) {
    case PressureBoundaryType::neumann:
      return "neumann";
    case PressureBoundaryType::dirichlet:
      return "dirichlet";
    case PressureBoundaryType::periodic:
      return "periodic";
  }

  return "unknown";
}

PressureBoundarySet derive_pressure_correction_boundary_conditions(
    const BoundaryConditionSet& boundary_conditions) {
  PressureBoundarySet mapped;

  for(int face_index = 0; face_index < static_cast<int>(mapped.faces.size()); ++face_index) {
    const BoundaryCondition& source = boundary_conditions.faces[static_cast<std::size_t>(face_index)];
    PressureBoundaryCondition& destination = mapped.faces[static_cast<std::size_t>(face_index)];

    switch(source.type) {
      case PhysicalBoundaryType::no_slip_wall:
      case PhysicalBoundaryType::prescribed_velocity:
      case PhysicalBoundaryType::symmetry:
        destination.type = PressureBoundaryType::neumann;
        destination.gradient = 0.0;
        break;
      case PhysicalBoundaryType::fixed_pressure:
        destination.type = PressureBoundaryType::dirichlet;
        destination.value = 0.0;
        break;
      case PhysicalBoundaryType::periodic:
        destination.type = PressureBoundaryType::periodic;
        break;
    }
  }

  for(const Axis axis : {Axis::x, Axis::y, Axis::z}) {
    pressure_axis_is_periodic(mapped, axis);
  }

  return mapped;
}

void apply_velocity_boundary_conditions(const BoundaryConditionSet& boundary_conditions,
                                        VelocityField& velocity) {
  apply_velocity_component_boundaries(boundary_conditions, velocity.x);
  apply_velocity_component_boundaries(boundary_conditions, velocity.y);
  apply_velocity_component_boundaries(boundary_conditions, velocity.z);
}

void apply_pressure_boundary_conditions(const PressureBoundarySet& boundary_conditions,
                                        PressureField& pressure) {
  for(const Axis axis : {Axis::x, Axis::y, Axis::z}) {
    if(pressure_axis_is_periodic(boundary_conditions, axis)) {
      apply_periodic_pair(static_cast<StructuredField&>(pressure), axis, false);
      continue;
    }

    apply_pressure_face_boundary(pressure, lower_face(axis), boundary_conditions[lower_face(axis)]);
    apply_pressure_face_boundary(pressure, upper_face(axis), boundary_conditions[upper_face(axis)]);
  }
}

void apply_total_pressure_boundary_conditions(const BoundaryConditionSet& boundary_conditions,
                                              const VelocityField& normal_pressure_gradient_source,
                                              PressureField& pressure_total) {
  for(const Axis axis : {Axis::x, Axis::y, Axis::z}) {
    if(velocity_axis_is_periodic(boundary_conditions, axis)) {
      apply_periodic_pair(static_cast<StructuredField&>(pressure_total), axis, false);
      continue;
    }

    for(const BoundaryFace face : {lower_face(axis), upper_face(axis)}) {
      const BoundaryCondition& condition = boundary_conditions[face];
      switch(condition.type) {
        case PhysicalBoundaryType::no_slip_wall:
        case PhysicalBoundaryType::prescribed_velocity:
        case PhysicalBoundaryType::symmetry:
          switch(axis) {
            case Axis::x:
              apply_total_pressure_gradient_face(pressure_total, face, normal_pressure_gradient_source.x);
              break;
            case Axis::y:
              apply_total_pressure_gradient_face(pressure_total, face, normal_pressure_gradient_source.y);
              break;
            case Axis::z:
              apply_total_pressure_gradient_face(pressure_total, face, normal_pressure_gradient_source.z);
              break;
          }
          break;
        case PhysicalBoundaryType::fixed_pressure: {
          const PressureBoundaryCondition pressure_boundary{
              .type = PressureBoundaryType::dirichlet,
              .value = condition.pressure,
          };
          apply_pressure_face_boundary(pressure_total, face, pressure_boundary);
          break;
        }
        case PhysicalBoundaryType::periodic:
          break;
      }
    }
  }
}

HelmholtzDiagnostics solve_predictor_adi(const VelocityField& rhs,
                                         const double alpha,
                                         const BoundaryConditionSet& boundary_conditions,
                                         VelocityField& predicted_velocity) {
  if(alpha < 0.0) {
    throw std::invalid_argument("solve_predictor_adi expects a non-negative alpha");
  }

  require_same_layout(rhs.x, predicted_velocity.x, "solve_predictor_adi");
  require_same_layout(rhs.y, predicted_velocity.y, "solve_predictor_adi");
  require_same_layout(rhs.z, predicted_velocity.z, "solve_predictor_adi");

  predicted_velocity = rhs;

  HelmholtzDiagnostics diagnostics;
  solve_component_adi(predicted_velocity.x, boundary_conditions, alpha, diagnostics.line_solves);
  solve_component_adi(predicted_velocity.y, boundary_conditions, alpha, diagnostics.line_solves);
  solve_component_adi(predicted_velocity.z, boundary_conditions, alpha, diagnostics.line_solves);
  return diagnostics;
}

void build_pressure_rhs(const VelocityField& predicted_velocity,
                        const BoundaryConditionSet& boundary_conditions,
                        const ProjectionOptions& options,
                        ScalarField& rhs) {
  if(options.dt <= 0.0) {
    throw std::invalid_argument("build_pressure_rhs expects dt > 0");
  }
  if(options.density <= 0.0) {
    throw std::invalid_argument("build_pressure_rhs expects density > 0");
  }

  VelocityField working = predicted_velocity;
  apply_velocity_boundary_conditions(boundary_conditions, working);
  operators::compute_divergence(working, rhs);
  scale_active(rhs, options.density / options.dt);
}

void correct_velocity(const VelocityField& predicted_velocity,
                      const PressureField& pressure_correction,
                      const BoundaryConditionSet& boundary_conditions,
                      const ProjectionOptions& options,
                      VelocityField& corrected_velocity) {
  if(options.dt <= 0.0) {
    throw std::invalid_argument("correct_velocity expects dt > 0");
  }
  if(options.density <= 0.0) {
    throw std::invalid_argument("correct_velocity expects density > 0");
  }

  require_same_layout(predicted_velocity.x, corrected_velocity.x, "correct_velocity");
  require_same_layout(predicted_velocity.y, corrected_velocity.y, "correct_velocity");
  require_same_layout(predicted_velocity.z, corrected_velocity.z, "correct_velocity");

  PressureField pressure_correction_with_ghosts = pressure_correction;
  apply_pressure_boundary_conditions(
      derive_pressure_correction_boundary_conditions(boundary_conditions),
      pressure_correction_with_ghosts);

  VelocityField pressure_gradient{pressure_correction.layout().grid()};
  operators::compute_gradient(pressure_correction_with_ghosts, pressure_gradient);

  corrected_velocity = predicted_velocity;
  axpy_active(corrected_velocity.x, pressure_gradient.x, -options.dt / options.density);
  axpy_active(corrected_velocity.y, pressure_gradient.y, -options.dt / options.density);
  axpy_active(corrected_velocity.z, pressure_gradient.z, -options.dt / options.density);
  apply_velocity_boundary_conditions(boundary_conditions, corrected_velocity);
}

ProjectionDiagnostics project_velocity(const VelocityField& predicted_velocity,
                                       const BoundaryConditionSet& boundary_conditions,
                                       const ProjectionOptions& options,
                                       PressureField& pressure_correction,
                                       VelocityField& corrected_velocity,
                                       ScalarField* pressure_rhs) {
  ProjectionDiagnostics diagnostics;
  const Grid& grid = pressure_correction.layout().grid();
  const PressureBoundarySet pressure_boundary_conditions =
      derive_pressure_correction_boundary_conditions(boundary_conditions);
  ScalarField rhs{grid};
  ScalarField divergence{grid};
  VelocityField working = predicted_velocity;

  apply_velocity_boundary_conditions(boundary_conditions, working);
  operators::compute_divergence(working, divergence);
  diagnostics.divergence_l2_before = active_l2_norm(divergence);

  build_pressure_rhs(working, boundary_conditions, options, rhs);
  diagnostics.rhs_l2 = active_l2_norm(rhs);
  diagnostics.pressure_solve = linsolve::solve_pressure_poisson(
      rhs, pressure_boundary_conditions, options, pressure_correction);
  correct_velocity(
      working, pressure_correction, boundary_conditions, options, corrected_velocity);

  operators::compute_divergence(corrected_velocity, divergence);
  diagnostics.divergence_l2_after = active_l2_norm(divergence);
  diagnostics.pressure_mean = active_mean(pressure_correction);

  if(pressure_rhs != nullptr) {
    copy_active(rhs, *pressure_rhs);
  }

  return diagnostics;
}

}  // namespace solver
