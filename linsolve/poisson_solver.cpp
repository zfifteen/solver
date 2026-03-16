#include "linsolve/poisson_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace solver::linsolve {

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

void require_same_layout(const StructuredField& left,
                         const StructuredField& right,
                         const char* operation_name) {
  if(!left.layout().same_shape_as(right.layout())) {
    throw std::invalid_argument(std::string(operation_name) + ": incompatible field layouts");
  }
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

double subtract_active_mean(StructuredField& field) {
  const double mean = active_mean(field);
  const IndexRange3D active = field.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        field(i, j, k) -= mean;
      }
    }
  }

  return mean;
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

double dot_active(const StructuredField& left, const StructuredField& right) {
  require_same_layout(left, right, "dot_active");
  const IndexRange3D active = left.layout().active_range();
  double sum = 0.0;

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        sum += left(i, j, k) * right(i, j, k);
      }
    }
  }

  return sum;
}

void copy_active(const StructuredField& source, StructuredField& destination) {
  require_same_layout(source, destination, "copy_active");
  const IndexRange3D active = source.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        destination(i, j, k) = source(i, j, k);
      }
    }
  }
}

void fill_storage(StructuredField& field, const double value) {
  const Extent3D storage = field.layout().storage_extent();

  for(int k = 0; k < storage.nz; ++k) {
    for(int j = 0; j < storage.ny; ++j) {
      for(int i = 0; i < storage.nx; ++i) {
        field(i, j, k) = value;
      }
    }
  }
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

PressureBoundarySet homogeneous_boundary_conditions(const PressureBoundarySet& actual) {
  PressureBoundarySet homogeneous = actual;

  for(PressureBoundaryCondition& boundary : homogeneous.faces) {
    boundary.value = 0.0;
    boundary.gradient = 0.0;
  }

  return homogeneous;
}

bool boundary_axis_is_periodic(const PressureBoundarySet& conditions, const Axis axis) {
  const PressureBoundaryCondition& lower = conditions[lower_face(axis)];
  const PressureBoundaryCondition& upper = conditions[upper_face(axis)];
  if((lower.type == PressureBoundaryType::periodic) !=
     (upper.type == PressureBoundaryType::periodic)) {
    throw std::invalid_argument("periodic pressure boundaries must be paired on both sides");
  }
  return lower.type == PressureBoundaryType::periodic;
}

bool pressure_problem_has_null_space(const PressureBoundarySet& boundary_conditions) {
  for(const PressureBoundaryCondition& boundary : boundary_conditions.faces) {
    if(boundary.type == PressureBoundaryType::dirichlet) {
      return false;
    }
  }

  return true;
}

double neighbor_value(const PressureField& input,
                      const PressureBoundarySet& boundary_conditions,
                      const int i,
                      const int j,
                      const int k,
                      const Axis axis,
                      const bool upper_neighbor) {
  const IndexRange3D active = input.layout().active_range();
  const int coordinate = axis == Axis::x ? i : axis == Axis::y ? j : k;
  const int begin = axis == Axis::x ? active.i_begin : axis == Axis::y ? active.j_begin : active.k_begin;
  const int end = axis == Axis::x ? active.i_end : axis == Axis::y ? active.j_end : active.k_end;

  if((!upper_neighbor && coordinate > begin) || (upper_neighbor && coordinate < end - 1)) {
    switch(axis) {
      case Axis::x:
        return input(i + (upper_neighbor ? 1 : -1), j, k);
      case Axis::y:
        return input(i, j + (upper_neighbor ? 1 : -1), k);
      case Axis::z:
        return input(i, j, k + (upper_neighbor ? 1 : -1));
    }
  }

  const BoundaryFace face = upper_neighbor ? upper_face(axis) : lower_face(axis);
  const PressureBoundaryCondition& boundary = boundary_conditions[face];
  const double center = input(i, j, k);
  const double spacing = input.layout().grid().spacing(axis);

  switch(boundary.type) {
    case PressureBoundaryType::periodic:
      switch(axis) {
        case Axis::x:
          return input(upper_neighbor ? active.i_begin : active.i_end - 1, j, k);
        case Axis::y:
          return input(i, upper_neighbor ? active.j_begin : active.j_end - 1, k);
        case Axis::z:
          return input(i, j, upper_neighbor ? active.k_begin : active.k_end - 1);
      }
      break;
    case PressureBoundaryType::neumann:
      return upper_neighbor ? center + boundary.gradient * spacing
                            : center - boundary.gradient * spacing;
    case PressureBoundaryType::dirichlet:
      return 2.0 * boundary.value - center;
  }

  return center;
}

void apply_poisson_operator(const PressureField& input,
                            const PressureBoundarySet& boundary_conditions,
                            PressureField& output) {
  require_same_layout(input, output, "apply_poisson_operator");
  const IndexRange3D active = input.layout().active_range();
  const Grid& grid = input.layout().grid();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        const double center = input(i, j, k);
        const double left = neighbor_value(input, boundary_conditions, i, j, k, Axis::x, false);
        const double right = neighbor_value(input, boundary_conditions, i, j, k, Axis::x, true);
        const double lower = neighbor_value(input, boundary_conditions, i, j, k, Axis::y, false);
        const double upper = neighbor_value(input, boundary_conditions, i, j, k, Axis::y, true);
        const double back = neighbor_value(input, boundary_conditions, i, j, k, Axis::z, false);
        const double front = neighbor_value(input, boundary_conditions, i, j, k, Axis::z, true);

        output(i, j, k) =
            -(right - 2.0 * center + left) / (grid.dx * grid.dx) -
            (upper - 2.0 * center + lower) / (grid.dy * grid.dy) -
            (front - 2.0 * center + back) / (grid.dz * grid.dz);
      }
    }
  }
}

void compute_operator_diagonal(const PressureField& template_field,
                               const PressureBoundarySet& boundary_conditions,
                               PressureField& diagonal) {
  const IndexRange3D active = template_field.layout().active_range();
  const Grid& grid = template_field.layout().grid();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        double value = 0.0;

        for(const Axis axis : {Axis::x, Axis::y, Axis::z}) {
          const double inverse_spacing_squared = 1.0 / std::pow(grid.spacing(axis), 2.0);
          double coefficient = 2.0 * inverse_spacing_squared;
          const int coordinate = axis == Axis::x ? i : axis == Axis::y ? j : k;
          const int begin =
              axis == Axis::x ? active.i_begin : axis == Axis::y ? active.j_begin : active.k_begin;
          const int end =
              axis == Axis::x ? active.i_end : axis == Axis::y ? active.j_end : active.k_end;

          if(boundary_axis_is_periodic(boundary_conditions, axis)) {
            value += coefficient;
            continue;
          }

          if(coordinate == begin) {
            const PressureBoundaryType type = boundary_conditions[lower_face(axis)].type;
            coefficient += type == PressureBoundaryType::dirichlet ? inverse_spacing_squared
                                                                   : -inverse_spacing_squared;
          }
          if(coordinate == end - 1) {
            const PressureBoundaryType type = boundary_conditions[upper_face(axis)].type;
            coefficient += type == PressureBoundaryType::dirichlet ? inverse_spacing_squared
                                                                   : -inverse_spacing_squared;
          }

          value += coefficient;
        }

        diagonal(i, j, k) = value;
      }
    }
  }
}

Grid coarsen_grid(const Grid& fine) {
  const auto next_cells = [](const int cells) {
    if(cells > 1 && cells % 2 == 0 && cells > 2) {
      return cells / 2;
    }
    return cells;
  };
  const auto next_spacing = [](const int fine_cells, const double fine_spacing) {
    if(fine_cells > 1 && fine_cells % 2 == 0 && fine_cells > 2) {
      return 2.0 * fine_spacing;
    }
    return fine_spacing;
  };

  return Grid{next_cells(fine.nx),
              next_cells(fine.ny),
              next_cells(fine.nz),
              next_spacing(fine.nx, fine.dx),
              next_spacing(fine.ny, fine.dy),
              next_spacing(fine.nz, fine.dz),
              fine.ghost_layers};
}

bool can_coarsen(const Grid& grid) {
  const auto axis_can_coarsen = [](const int cells) {
    return cells > 1 && cells % 2 == 0 && cells > 2;
  };

  return axis_can_coarsen(grid.nx) || axis_can_coarsen(grid.ny) || axis_can_coarsen(grid.nz);
}

struct HierarchyLevel {
  Grid grid;
  PressureField correction;
  PressureField rhs;
  PressureField residual;
  PressureField scratch;
  PressureField diagonal;

  explicit HierarchyLevel(const Grid& grid_in)
      : grid(grid_in),
        correction(grid),
        rhs(grid),
        residual(grid),
        scratch(grid),
        diagonal(grid) {}
};

std::vector<HierarchyLevel> build_hierarchy(const Grid& finest_grid,
                                            const PressureBoundarySet& homogeneous_boundaries) {
  std::vector<HierarchyLevel> hierarchy;
  hierarchy.emplace_back(finest_grid);
  compute_operator_diagonal(hierarchy.back().correction, homogeneous_boundaries,
                            hierarchy.back().diagonal);

  while(can_coarsen(hierarchy.back().grid)) {
    const Grid coarse_grid = coarsen_grid(hierarchy.back().grid);
    if(coarse_grid.nx == hierarchy.back().grid.nx && coarse_grid.ny == hierarchy.back().grid.ny &&
       coarse_grid.nz == hierarchy.back().grid.nz) {
      break;
    }

    hierarchy.emplace_back(coarse_grid);
    compute_operator_diagonal(hierarchy.back().correction, homogeneous_boundaries,
                              hierarchy.back().diagonal);
  }

  return hierarchy;
}

int refinement_ratio(const int fine_cells, const int coarse_cells) {
  if(fine_cells == coarse_cells) {
    return 1;
  }
  if(fine_cells == 2 * coarse_cells) {
    return 2;
  }

  throw std::invalid_argument("refinement_ratio expects equal or factor-of-two grids");
}

void restrict_full_weighting(const PressureField& fine, PressureField& coarse) {
  const Extent3D fine_extent = fine.layout().active_extent();
  const Extent3D coarse_extent = coarse.layout().active_extent();
  const int ratio_x = refinement_ratio(fine_extent.nx, coarse_extent.nx);
  const int ratio_y = refinement_ratio(fine_extent.ny, coarse_extent.ny);
  const int ratio_z = refinement_ratio(fine_extent.nz, coarse_extent.nz);
  const double inverse_children = 1.0 / static_cast<double>(ratio_x * ratio_y * ratio_z);

  for(int kc = 0; kc < coarse_extent.nz; ++kc) {
    for(int jc = 0; jc < coarse_extent.ny; ++jc) {
      for(int ic = 0; ic < coarse_extent.nx; ++ic) {
        double sum = 0.0;
        for(int dz = 0; dz < ratio_z; ++dz) {
          for(int dy = 0; dy < ratio_y; ++dy) {
            for(int dx = 0; dx < ratio_x; ++dx) {
              const Index3D fine_index = fine.layout().unchecked_storage_index_from_active(
                  ic * ratio_x + dx, jc * ratio_y + dy, kc * ratio_z + dz);
              sum += fine(fine_index.i, fine_index.j, fine_index.k);
            }
          }
        }

        const Index3D coarse_index = coarse.layout().unchecked_storage_index_from_active(ic, jc, kc);
        coarse(coarse_index.i, coarse_index.j, coarse_index.k) = inverse_children * sum;
      }
    }
  }
}

void prolongate_and_add(const PressureField& coarse, PressureField& fine) {
  const Extent3D fine_extent = fine.layout().active_extent();
  const Extent3D coarse_extent = coarse.layout().active_extent();
  const int ratio_x = refinement_ratio(fine_extent.nx, coarse_extent.nx);
  const int ratio_y = refinement_ratio(fine_extent.ny, coarse_extent.ny);
  const int ratio_z = refinement_ratio(fine_extent.nz, coarse_extent.nz);

  for(int kc = 0; kc < coarse_extent.nz; ++kc) {
    for(int jc = 0; jc < coarse_extent.ny; ++jc) {
      for(int ic = 0; ic < coarse_extent.nx; ++ic) {
        const Index3D coarse_index = coarse.layout().unchecked_storage_index_from_active(ic, jc, kc);
        const double coarse_value = coarse(coarse_index.i, coarse_index.j, coarse_index.k);

        for(int dz = 0; dz < ratio_z; ++dz) {
          for(int dy = 0; dy < ratio_y; ++dy) {
            for(int dx = 0; dx < ratio_x; ++dx) {
              const Index3D fine_index = fine.layout().unchecked_storage_index_from_active(
                  ic * ratio_x + dx, jc * ratio_y + dy, kc * ratio_z + dz);
              fine(fine_index.i, fine_index.j, fine_index.k) += coarse_value;
            }
          }
        }
      }
    }
  }
}

std::vector<double> flatten_active(const PressureField& field) {
  const IndexRange3D active = field.layout().active_range();
  std::vector<double> flattened;
  flattened.reserve(active.extent().cell_count());

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        flattened.push_back(field(i, j, k));
      }
    }
  }

  return flattened;
}

void unflatten_active(const std::vector<double>& values, PressureField& field) {
  const IndexRange3D active = field.layout().active_range();
  std::size_t cursor = 0;

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        field(i, j, k) = values[cursor++];
      }
    }
  }
}

void gaussian_elimination(std::vector<double>& matrix, std::vector<double>& rhs) {
  const std::size_t n = rhs.size();

  for(std::size_t pivot = 0; pivot < n; ++pivot) {
    std::size_t best = pivot;
    double best_abs = std::abs(matrix[pivot * n + pivot]);
    for(std::size_t row = pivot + 1; row < n; ++row) {
      const double candidate = std::abs(matrix[row * n + pivot]);
      if(candidate > best_abs) {
        best = row;
        best_abs = candidate;
      }
    }

    if(best_abs < kSmallNumber) {
      throw std::runtime_error("direct coarse solve encountered a singular matrix");
    }

    if(best != pivot) {
      for(std::size_t column = 0; column < n; ++column) {
        std::swap(matrix[pivot * n + column], matrix[best * n + column]);
      }
      std::swap(rhs[pivot], rhs[best]);
    }

    const double inverse_pivot = 1.0 / matrix[pivot * n + pivot];
    for(std::size_t row = pivot + 1; row < n; ++row) {
      const double factor = matrix[row * n + pivot] * inverse_pivot;
      if(std::abs(factor) < kSmallNumber) {
        continue;
      }

      for(std::size_t column = pivot; column < n; ++column) {
        matrix[row * n + column] -= factor * matrix[pivot * n + column];
      }
      rhs[row] -= factor * rhs[pivot];
    }
  }

  for(std::size_t row = n; row-- > 0;) {
    double sum = rhs[row];
    for(std::size_t column = row + 1; column < n; ++column) {
      sum -= matrix[row * n + column] * rhs[column];
    }
    rhs[row] = sum / matrix[row * n + row];
  }
}

void direct_coarse_solve(HierarchyLevel& level,
                         const PressureBoundarySet& boundary_conditions,
                         const bool enforce_zero_mean) {
  const std::vector<double> rhs_vector = flatten_active(level.rhs);
  const std::size_t n = rhs_vector.size();
  std::vector<double> matrix(n * n, 0.0);
  std::vector<double> solution_rhs = rhs_vector;
  PressureField basis(level.grid);
  PressureField image(level.grid);

  for(std::size_t column = 0; column < n; ++column) {
    fill_storage(basis, 0.0);
    std::vector<double> basis_vector(n, 0.0);
    basis_vector[column] = 1.0;
    unflatten_active(basis_vector, basis);

    apply_poisson_operator(basis, boundary_conditions, image);
    const std::vector<double> image_vector = flatten_active(image);
    for(std::size_t row = 0; row < n; ++row) {
      matrix[row * n + column] = image_vector[row];
    }
  }

  if(enforce_zero_mean) {
    for(std::size_t column = 0; column < n; ++column) {
      matrix[(n - 1) * n + column] = 1.0;
    }
    std::fill(solution_rhs.begin() + static_cast<std::ptrdiff_t>(n - 1), solution_rhs.end(), 0.0);
    solution_rhs[n - 1] = 0.0;
  }

  gaussian_elimination(matrix, solution_rhs);
  fill_storage(level.correction, 0.0);
  unflatten_active(solution_rhs, level.correction);
  if(enforce_zero_mean) {
    subtract_active_mean(level.correction);
  }
}

void smooth_jacobi(HierarchyLevel& level,
                   const PressureBoundarySet& boundary_conditions,
                   const MultigridPolicy& policy,
                   const int iterations,
                   const bool enforce_zero_mean) {
  for(int iteration = 0; iteration < iterations; ++iteration) {
    apply_poisson_operator(level.correction, boundary_conditions, level.scratch);
    level.residual = level.rhs;
    axpy_active(level.residual, level.scratch, -1.0);

    const IndexRange3D active = level.correction.layout().active_range();
    for(int k = active.k_begin; k < active.k_end; ++k) {
      for(int j = active.j_begin; j < active.j_end; ++j) {
        for(int i = active.i_begin; i < active.i_end; ++i) {
          level.scratch(i, j, k) =
              level.correction(i, j, k) +
              policy.jacobi_omega * level.residual(i, j, k) / level.diagonal(i, j, k);
        }
      }
    }

    copy_active(level.scratch, level.correction);
    if(enforce_zero_mean) {
      subtract_active_mean(level.correction);
    }
  }
}

void v_cycle(std::vector<HierarchyLevel>& hierarchy,
             const std::size_t level_index,
             const PressureBoundarySet& boundary_conditions,
             const MultigridPolicy& policy,
             const bool enforce_zero_mean) {
  HierarchyLevel& level = hierarchy[level_index];
  fill_storage(level.correction, 0.0);

  if(level_index + 1 == hierarchy.size()) {
    direct_coarse_solve(level, boundary_conditions, enforce_zero_mean);
    return;
  }

  smooth_jacobi(level, boundary_conditions, policy, policy.pre_smoothing_steps, enforce_zero_mean);
  apply_poisson_operator(level.correction, boundary_conditions, level.scratch);
  level.residual = level.rhs;
  axpy_active(level.residual, level.scratch, -1.0);
  if(enforce_zero_mean) {
    subtract_active_mean(level.residual);
  }

  HierarchyLevel& coarse = hierarchy[level_index + 1];
  fill_storage(coarse.rhs, 0.0);
  fill_storage(coarse.correction, 0.0);
  restrict_full_weighting(level.residual, coarse.rhs);
  if(enforce_zero_mean) {
    subtract_active_mean(coarse.rhs);
  }

  v_cycle(hierarchy, level_index + 1, boundary_conditions, policy, enforce_zero_mean);
  prolongate_and_add(coarse.correction, level.correction);
  if(enforce_zero_mean) {
    subtract_active_mean(level.correction);
  }

  smooth_jacobi(level, boundary_conditions, policy, policy.post_smoothing_steps, enforce_zero_mean);
}

PressureField build_boundary_source(const Grid& grid, const PressureBoundarySet& boundary_conditions) {
  PressureField zero_pressure(grid);
  PressureField source(grid);
  fill_storage(zero_pressure, 0.0);
  fill_storage(source, 0.0);
  apply_poisson_operator(zero_pressure, boundary_conditions, source);
  return source;
}

}  // namespace

std::string to_string(const MultigridCycle cycle) {
  switch(cycle) {
    case MultigridCycle::v_cycle:
      return "v_cycle";
  }

  return "unknown";
}

std::string to_string(const MultigridSmoother smoother) {
  switch(smoother) {
    case MultigridSmoother::damped_jacobi:
      return "damped_jacobi";
  }

  return "unknown";
}

MultigridPolicy default_multigrid_policy() {
  return MultigridPolicy{};
}

void build_poisson_rhs_from_pressure(const PressureField& pressure,
                                     const PressureBoundarySet& boundary_conditions,
                                     ScalarField& rhs) {
  require_same_layout(pressure, rhs, "build_poisson_rhs_from_pressure");
  PressureField system_rhs(pressure.layout().grid());
  apply_poisson_operator(pressure, boundary_conditions, system_rhs);
  copy_active(system_rhs, rhs);
  scale_active(rhs, -1.0);
}

PoissonSolveDiagnostics solve_pressure_poisson(const ScalarField& rhs,
                                               const PressureBoundarySet& boundary_conditions,
                                               const ProjectionOptions& options,
                                               PressureField& pressure) {
  if(options.poisson_max_iterations <= 0) {
    throw std::invalid_argument("solve_pressure_poisson expects iterations > 0");
  }
  if(options.poisson_tolerance <= 0.0) {
    throw std::invalid_argument("solve_pressure_poisson expects tolerance > 0");
  }

  require_same_layout(rhs, pressure, "solve_pressure_poisson");

  const MultigridPolicy policy = default_multigrid_policy();
  const PressureBoundarySet homogeneous_boundaries = homogeneous_boundary_conditions(boundary_conditions);
  const bool enforce_zero_mean = pressure_problem_has_null_space(boundary_conditions);
  const Grid& grid = pressure.layout().grid();

  PressureField system_rhs(grid);
  PressureField residual(grid);
  PressureField search_direction(grid);
  PressureField operator_result(grid);
  PressureField preconditioned_residual(grid);
  PressureField boundary_source = build_boundary_source(grid, boundary_conditions);

  copy_active(rhs, system_rhs);
  scale_active(system_rhs, -1.0);
  axpy_active(system_rhs, boundary_source, -1.0);

  PoissonSolveDiagnostics diagnostics;
  diagnostics.zero_mean_enforced = enforce_zero_mean;
  diagnostics.pre_smoothing_steps = policy.pre_smoothing_steps;
  diagnostics.post_smoothing_steps = policy.post_smoothing_steps;

  std::vector<HierarchyLevel> hierarchy = build_hierarchy(grid, homogeneous_boundaries);
  diagnostics.multigrid_levels = static_cast<int>(hierarchy.size());
  diagnostics.coarse_unknowns = static_cast<int>(hierarchy.back().grid.nx * hierarchy.back().grid.ny *
                                                 hierarchy.back().grid.nz);

  if(enforce_zero_mean) {
    diagnostics.removed_rhs_mean = subtract_active_mean(system_rhs);
    subtract_active_mean(pressure);
  }

  apply_poisson_operator(pressure, homogeneous_boundaries, operator_result);
  residual = system_rhs;
  axpy_active(residual, operator_result, -1.0);
  if(enforce_zero_mean) {
    subtract_active_mean(residual);
  }

  diagnostics.initial_residual_l2 = active_l2_norm(residual);
  diagnostics.final_residual_l2 = diagnostics.initial_residual_l2;
  diagnostics.relative_residual =
      diagnostics.initial_residual_l2 > 0.0 ? 1.0 : 0.0;
  diagnostics.residual_history.push_back(diagnostics.initial_residual_l2);

  if(diagnostics.initial_residual_l2 <= options.poisson_tolerance) {
    diagnostics.converged = true;
    diagnostics.relative_residual = 0.0;
    return diagnostics;
  }

  hierarchy.front().rhs = residual;
  v_cycle(hierarchy, 0, homogeneous_boundaries, policy, enforce_zero_mean);
  preconditioned_residual = hierarchy.front().correction;
  if(enforce_zero_mean) {
    subtract_active_mean(preconditioned_residual);
  }
  search_direction = preconditioned_residual;

  double residual_dot_preconditioned = dot_active(residual, preconditioned_residual);
  for(int iteration = 0; iteration < options.poisson_max_iterations; ++iteration) {
    apply_poisson_operator(search_direction, homogeneous_boundaries, operator_result);
    if(enforce_zero_mean) {
      subtract_active_mean(operator_result);
    }

    const double search_dot = dot_active(search_direction, operator_result);
    if(std::abs(search_dot) < kSmallNumber) {
      break;
    }

    const double alpha = residual_dot_preconditioned / search_dot;
    axpy_active(pressure, search_direction, alpha);
    if(enforce_zero_mean) {
      subtract_active_mean(pressure);
    }

    axpy_active(residual, operator_result, -alpha);
    if(enforce_zero_mean) {
      subtract_active_mean(residual);
    }

    diagnostics.iterations = iteration + 1;
    diagnostics.final_residual_l2 = active_l2_norm(residual);
    diagnostics.relative_residual =
        diagnostics.initial_residual_l2 > 0.0
            ? diagnostics.final_residual_l2 / diagnostics.initial_residual_l2
            : 0.0;
    diagnostics.residual_history.push_back(diagnostics.final_residual_l2);
    if(diagnostics.relative_residual <= options.poisson_tolerance) {
      diagnostics.converged = true;
      break;
    }

    hierarchy.front().rhs = residual;
    v_cycle(hierarchy, 0, homogeneous_boundaries, policy, enforce_zero_mean);
    preconditioned_residual = hierarchy.front().correction;
    if(enforce_zero_mean) {
      subtract_active_mean(preconditioned_residual);
    }

    const double next_dot = dot_active(residual, preconditioned_residual);
    const double beta = next_dot / residual_dot_preconditioned;
    const IndexRange3D active = search_direction.layout().active_range();
    for(int k = active.k_begin; k < active.k_end; ++k) {
      for(int j = active.j_begin; j < active.j_end; ++j) {
        for(int i = active.i_begin; i < active.i_end; ++i) {
          search_direction(i, j, k) =
              preconditioned_residual(i, j, k) + beta * search_direction(i, j, k);
        }
      }
    }
    if(enforce_zero_mean) {
      subtract_active_mean(search_direction);
    }
    residual_dot_preconditioned = next_dot;
  }

  return diagnostics;
}

}  // namespace solver::linsolve
