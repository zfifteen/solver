#include "core/fields.hpp"
#include "core/grid.hpp"
#include "core/runtime.hpp"
#include "io/checkpoint.hpp"
#include "io/vtk_export.hpp"
#include "linsolve/poisson_solver.hpp"
#include "metal/taylor_green_backend.hpp"
#include "operators/discrete_operators.hpp"
#include "solver/channel_flow.hpp"
#include "solver/lid_driven_cavity.hpp"
#include "solver/momentum_terms.hpp"
#include "solver/operator_verification.hpp"
#include "solver/projection.hpp"
#include "solver/taylor_green.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

template <typename Fn>
void require_exception_contains(Fn&& fn,
                                const std::string& expected_fragment,
                                const std::string& message) {
  try {
    fn();
  } catch(const std::exception& exception) {
    require(std::string(exception.what()).find(expected_fragment) != std::string::npos,
            message + ": wrong exception message");
    return;
  }

  throw std::runtime_error(message);
}

double pi() {
  return std::acos(-1.0);
}

std::string source_path(const std::string& relative_path) {
  return std::string(SOLVER_SOURCE_DIR) + "/" + relative_path;
}

std::filesystem::path temp_path(const std::string& filename) {
  return std::filesystem::temp_directory_path() / filename;
}

std::vector<std::uint8_t> read_binary_file(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  require(input.is_open(), "read_binary_file could not open path: " + path.string());
  return std::vector<std::uint8_t>(std::istreambuf_iterator<char>(input),
                                   std::istreambuf_iterator<char>());
}

void write_binary_file(const std::filesystem::path& path,
                       const std::vector<std::uint8_t>& bytes) {
  std::ofstream output(path, std::ios::binary);
  require(output.is_open(), "write_binary_file could not open path: " + path.string());
  output.write(reinterpret_cast<const char*>(bytes.data()),
               static_cast<std::streamsize>(bytes.size()));
  require(output.good(), "write_binary_file failed for path: " + path.string());
}

std::uint64_t read_le_u64(const std::vector<std::uint8_t>& bytes,
                          const std::size_t offset) {
  require(offset + 8 <= bytes.size(), "read_le_u64 out of range");
  std::uint64_t value = 0;
  for(int shift = 0; shift < 64; shift += 8) {
    value |= static_cast<std::uint64_t>(bytes[offset + static_cast<std::size_t>(shift / 8)]) << shift;
  }
  return value;
}

void write_le_u32(std::vector<std::uint8_t>& bytes,
                  const std::size_t offset,
                  const std::uint32_t value) {
  require(offset + 4 <= bytes.size(), "write_le_u32 out of range");
  for(int shift = 0; shift < 32; shift += 8) {
    bytes[offset + static_cast<std::size_t>(shift / 8)] =
        static_cast<std::uint8_t>((value >> shift) & 0xffu);
  }
}

void write_le_u64(std::vector<std::uint8_t>& bytes,
                  const std::size_t offset,
                  const std::uint64_t value) {
  require(offset + 8 <= bytes.size(), "write_le_u64 out of range");
  for(int shift = 0; shift < 64; shift += 8) {
    bytes[offset + static_cast<std::size_t>(shift / 8)] =
        static_cast<std::uint8_t>((value >> shift) & 0xffu);
  }
}

std::uint64_t fnv1a_64(const std::uint8_t* data, const std::size_t size) {
  std::uint64_t hash = 14695981039346656037ull;
  for(std::size_t index = 0; index < size; ++index) {
    hash ^= static_cast<std::uint64_t>(data[index]);
    hash *= 1099511628211ull;
  }
  return hash;
}

struct CheckpointSectionView {
  std::size_t payload_offset = 0;
  std::size_t payload_size = 0;
};

CheckpointSectionView find_checkpoint_section(const std::vector<std::uint8_t>& bytes,
                                              const std::array<char, 4>& tag) {
  constexpr std::size_t kMagicBytes = 8;
  constexpr std::size_t kHeaderBytes = 4 + 8 + 8;
  std::size_t cursor = kMagicBytes + kHeaderBytes;

  while(cursor < bytes.size()) {
    require(cursor + 12 <= bytes.size(), "checkpoint section header truncated in test helper");
    const std::array<char, 4> current_tag{
        static_cast<char>(bytes[cursor]),
        static_cast<char>(bytes[cursor + 1]),
        static_cast<char>(bytes[cursor + 2]),
        static_cast<char>(bytes[cursor + 3]),
    };
    const std::uint64_t section_size = read_le_u64(bytes, cursor + 4);
    const std::size_t payload_offset = cursor + 12;
    require(payload_offset + section_size <= bytes.size(),
            "checkpoint section payload truncated in test helper");
    if(current_tag == tag) {
      return CheckpointSectionView{
          .payload_offset = payload_offset,
          .payload_size = static_cast<std::size_t>(section_size),
      };
    }
    cursor = payload_offset + static_cast<std::size_t>(section_size);
  }

  throw std::runtime_error("checkpoint section not found in test helper");
}

void refresh_checkpoint_payload_header(std::vector<std::uint8_t>& bytes) {
  constexpr std::size_t kPayloadOffset = 8 + 4 + 8 + 8;
  require(bytes.size() >= kPayloadOffset,
          "refresh_checkpoint_payload_header requires a full checkpoint header");
  const std::uint64_t payload_size = bytes.size() - kPayloadOffset;
  write_le_u64(bytes, 12, payload_size);
  const std::uint64_t checksum =
      fnv1a_64(bytes.data() + kPayloadOffset, static_cast<std::size_t>(payload_size));
  write_le_u64(bytes, 20, checksum);
}

bool metal_backend_available() {
  static const bool available = solver::metal::is_backend_available();
  return available;
}

double square(const double value) {
  return value * value;
}

template <typename Field, typename ValueFn>
void fill_storage(Field& field, ValueFn&& value_fn) {
  const solver::Extent3D storage = field.layout().storage_extent();

  for(int k = 0; k < storage.nz; ++k) {
    const double z = field.layout().coordinate_for_storage_index(solver::Axis::z, k);
    for(int j = 0; j < storage.ny; ++j) {
      const double y = field.layout().coordinate_for_storage_index(solver::Axis::y, j);
      for(int i = 0; i < storage.nx; ++i) {
        const double x = field.layout().coordinate_for_storage_index(solver::Axis::x, i);
        field(i, j, k) = value_fn(x, y, z);
      }
    }
  }
}

template <typename Field, typename ValueFn>
double active_l2_error(const Field& field, ValueFn&& exact_value) {
  const solver::IndexRange3D active = field.layout().active_range();

  double sum = 0.0;
  std::size_t count = 0;
  for(int k = active.k_begin; k < active.k_end; ++k) {
    const double z = field.layout().coordinate_for_storage_index(solver::Axis::z, k);
    for(int j = active.j_begin; j < active.j_end; ++j) {
      const double y = field.layout().coordinate_for_storage_index(solver::Axis::y, j);
      for(int i = active.i_begin; i < active.i_end; ++i) {
        const double x = field.layout().coordinate_for_storage_index(solver::Axis::x, i);
        sum += square(field(i, j, k) - exact_value(x, y, z));
        ++count;
      }
    }
  }

  return std::sqrt(sum / static_cast<double>(count));
}

double observed_order(const double coarse_error, const double fine_error) {
  return std::log(coarse_error / fine_error) / std::log(2.0);
}

double wrap_periodic(const double coordinate, const double period) {
  const double wrapped = std::fmod(coordinate, period);
  return wrapped < 0.0 ? wrapped + period : wrapped;
}

double kinetic_energy(const solver::VelocityField& velocity) {
  const solver::Grid& grid = velocity.x.layout().grid();
  const solver::IndexRange3D cells = solver::FieldLayout::cell_centered(grid).active_range();
  double energy_sum = 0.0;
  std::size_t count = 0;

  for(int k = cells.k_begin; k < cells.k_end; ++k) {
    for(int j = cells.j_begin; j < cells.j_end; ++j) {
      for(int i = cells.i_begin; i < cells.i_end; ++i) {
        const double u_center = 0.5 * (velocity.x(i, j, k) + velocity.x(i + 1, j, k));
        const double v_center = 0.5 * (velocity.y(i, j, k) + velocity.y(i, j + 1, k));
        const double w_center = 0.5 * (velocity.z(i, j, k) + velocity.z(i, j, k + 1));
        energy_sum += 0.5 * (square(u_center) + square(v_center) + square(w_center));
        ++count;
      }
    }
  }

  return energy_sum / static_cast<double>(count);
}

double active_min(const solver::StructuredField& field) {
  const solver::IndexRange3D active = field.layout().active_range();
  double value = field(active.i_begin, active.j_begin, active.k_begin);

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        value = std::min(value, field(i, j, k));
      }
    }
  }

  return value;
}

double active_max(const solver::StructuredField& field) {
  const solver::IndexRange3D active = field.layout().active_range();
  double value = field(active.i_begin, active.j_begin, active.k_begin);

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        value = std::max(value, field(i, j, k));
      }
    }
  }

  return value;
}

double max_abs_active(const solver::StructuredField& field) {
  const solver::IndexRange3D active = field.layout().active_range();
  double value = 0.0;

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        value = std::max(value, std::abs(field(i, j, k)));
      }
    }
  }

  return value;
}

double active_mean_value(const solver::StructuredField& field) {
  const solver::IndexRange3D active = field.layout().active_range();
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

template <typename Field>
double active_l2_difference(const Field& left, const Field& right) {
  const solver::IndexRange3D active = left.layout().active_range();
  double sum = 0.0;
  std::size_t count = 0;

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        sum += square(left(i, j, k) - right(i, j, k));
        ++count;
      }
    }
  }

  return std::sqrt(sum / static_cast<double>(count));
}

template <typename Field>
double active_l2_norm_value(const Field& field) {
  const solver::IndexRange3D active = field.layout().active_range();
  double sum = 0.0;
  std::size_t count = 0;

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        sum += square(field(i, j, k));
        ++count;
      }
    }
  }

  return std::sqrt(sum / static_cast<double>(count));
}

template <typename Field>
double relative_active_l2_difference(const Field& left, const Field& right) {
  const double denominator = active_l2_norm_value(right);
  require(denominator > 0.0, "relative_active_l2_difference requires non-zero reference norm");
  return active_l2_difference(left, right) / denominator;
}

bool storage_bitwise_equal(const solver::StructuredField& left, const solver::StructuredField& right) {
  const solver::Extent3D left_storage = left.layout().storage_extent();
  const solver::Extent3D right_storage = right.layout().storage_extent();
  if(left_storage.nx != right_storage.nx || left_storage.ny != right_storage.ny ||
     left_storage.nz != right_storage.nz) {
    return false;
  }

  for(int k = 0; k < left_storage.nz; ++k) {
    for(int j = 0; j < left_storage.ny; ++j) {
      for(int i = 0; i < left_storage.nx; ++i) {
        if(std::bit_cast<std::uint64_t>(left(i, j, k)) !=
           std::bit_cast<std::uint64_t>(right(i, j, k))) {
          return false;
        }
      }
    }
  }

  return true;
}

bool metrics_bitwise_equal(const solver::SimulationStepMetrics& left,
                           const solver::SimulationStepMetrics& right) {
  return left.step == right.step &&
         std::bit_cast<std::uint64_t>(left.time) == std::bit_cast<std::uint64_t>(right.time) &&
         std::bit_cast<std::uint64_t>(left.dt) == std::bit_cast<std::uint64_t>(right.dt) &&
         std::bit_cast<std::uint64_t>(left.max_cfl) ==
             std::bit_cast<std::uint64_t>(right.max_cfl) &&
         std::bit_cast<std::uint64_t>(left.max_velocity_change) ==
             std::bit_cast<std::uint64_t>(right.max_velocity_change) &&
         std::bit_cast<std::uint64_t>(left.divergence_l2) ==
             std::bit_cast<std::uint64_t>(right.divergence_l2) &&
         left.pressure_iterations == right.pressure_iterations &&
         std::bit_cast<std::uint64_t>(left.pressure_relative_residual) ==
             std::bit_cast<std::uint64_t>(right.pressure_relative_residual);
}

template <typename Field>
void subtract_active_mean_inplace(Field& field) {
  const double mean = active_mean_value(field);
  const solver::IndexRange3D active = field.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        field(i, j, k) -= mean;
      }
    }
  }
}

template <typename FromField, typename ToField>
void copy_active_values(const FromField& source, ToField& destination) {
  const solver::IndexRange3D active = source.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        destination(i, j, k) = source(i, j, k);
      }
    }
  }
}

template <typename DestinationField, typename SourceField>
void axpy_active(DestinationField& destination, const SourceField& source, const double scale) {
  const solver::IndexRange3D active = destination.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        destination(i, j, k) += scale * source(i, j, k);
      }
    }
  }
}

void axpy_velocity(solver::VelocityField& destination,
                   const solver::VelocityField& source,
                   const double scale) {
  axpy_active(destination.x, source.x, scale);
  axpy_active(destination.y, source.y, scale);
  axpy_active(destination.z, source.z, scale);
}

template <typename Field>
void scale_active(Field& field, const double scale) {
  const solver::IndexRange3D active = field.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        field(i, j, k) *= scale;
      }
    }
  }
}

void scale_velocity_active(solver::VelocityField& field, const double scale) {
  scale_active(field.x, scale);
  scale_active(field.y, scale);
  scale_active(field.z, scale);
}

bool is_fixed_velocity_boundary(const solver::BoundaryCondition& condition) {
  return condition.type == solver::PhysicalBoundaryType::no_slip_wall ||
         condition.type == solver::PhysicalBoundaryType::prescribed_velocity ||
         condition.type == solver::PhysicalBoundaryType::symmetry;
}

void copy_storage_values(const solver::StructuredField& source, solver::StructuredField& destination) {
  const solver::Extent3D storage = source.layout().storage_extent();
  for(int k = 0; k < storage.nz; ++k) {
    for(int j = 0; j < storage.ny; ++j) {
      for(int i = 0; i < storage.nx; ++i) {
        destination(i, j, k) = source(i, j, k);
      }
    }
  }
}

void apply_component_boundary_conditions(const solver::BoundaryConditionSet& boundary_conditions,
                                         solver::FaceField& field) {
  solver::VelocityField velocity{field.layout().grid()};
  velocity.fill(0.0);

  switch(field.normal_axis()) {
    case solver::Axis::x:
      copy_storage_values(field, velocity.x);
      break;
    case solver::Axis::y:
      copy_storage_values(field, velocity.y);
      break;
    case solver::Axis::z:
      copy_storage_values(field, velocity.z);
      break;
  }

  solver::apply_velocity_boundary_conditions(boundary_conditions, velocity);

  switch(field.normal_axis()) {
    case solver::Axis::x:
      copy_storage_values(velocity.x, field);
      break;
    case solver::Axis::y:
      copy_storage_values(velocity.y, field);
      break;
    case solver::Axis::z:
      copy_storage_values(velocity.z, field);
      break;
  }
}

std::vector<solver::Index3D> collect_factorized_unknowns(const solver::FaceField& field,
                                                         const solver::BoundaryConditionSet& boundary_conditions) {
  solver::IndexRange3D active = field.layout().active_range();

  if(field.normal_axis() == solver::Axis::x) {
    if(is_fixed_velocity_boundary(boundary_conditions[solver::BoundaryFace::x_min])) {
      ++active.i_begin;
    }
    if(is_fixed_velocity_boundary(boundary_conditions[solver::BoundaryFace::x_max])) {
      --active.i_end;
    }
  } else if(field.normal_axis() == solver::Axis::y) {
    if(is_fixed_velocity_boundary(boundary_conditions[solver::BoundaryFace::y_min])) {
      ++active.j_begin;
    }
    if(is_fixed_velocity_boundary(boundary_conditions[solver::BoundaryFace::y_max])) {
      --active.j_end;
    }
  } else {
    if(is_fixed_velocity_boundary(boundary_conditions[solver::BoundaryFace::z_min])) {
      ++active.k_begin;
    }
    if(is_fixed_velocity_boundary(boundary_conditions[solver::BoundaryFace::z_max])) {
      --active.k_end;
    }
  }

  std::vector<solver::Index3D> unknowns;
  unknowns.reserve(static_cast<std::size_t>(active.extent().cell_count()));
  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        unknowns.push_back(solver::Index3D{i, j, k});
      }
    }
  }
  return unknowns;
}

std::vector<double> extract_unknown_values(const solver::FaceField& field,
                                           const std::vector<solver::Index3D>& unknowns) {
  std::vector<double> values;
  values.reserve(unknowns.size());
  for(const solver::Index3D index : unknowns) {
    values.push_back(field(index.i, index.j, index.k));
  }
  return values;
}

double second_derivative_axis(const solver::StructuredField& input,
                              const int i,
                              const int j,
                              const int k,
                              const solver::Axis axis) {
  const double spacing = input.layout().grid().spacing(axis);
  const double inverse_spacing_squared = 1.0 / (spacing * spacing);

  switch(axis) {
    case solver::Axis::x:
      return (input(i + 1, j, k) - 2.0 * input(i, j, k) + input(i - 1, j, k)) *
             inverse_spacing_squared;
    case solver::Axis::y:
      return (input(i, j + 1, k) - 2.0 * input(i, j, k) + input(i, j - 1, k)) *
             inverse_spacing_squared;
    case solver::Axis::z:
      return (input(i, j, k + 1) - 2.0 * input(i, j, k) + input(i, j, k - 1)) *
             inverse_spacing_squared;
  }

  return 0.0;
}

void apply_directional_helmholtz_operator(const solver::FaceField& input,
                                          const double alpha,
                                          const solver::Axis axis,
                                          solver::FaceField& output) {
  copy_storage_values(input, output);
  const solver::IndexRange3D active = input.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        output(i, j, k) = input(i, j, k) - alpha * second_derivative_axis(input, i, j, k, axis);
      }
    }
  }
}

solver::FaceField apply_factorized_operator(const solver::FaceField& input,
                                            const solver::BoundaryConditionSet& boundary_conditions,
                                            const double alpha) {
  solver::FaceField stage0(input.normal_axis(), input.layout().grid());
  solver::FaceField stage1(input.normal_axis(), input.layout().grid());
  solver::FaceField stage2(input.normal_axis(), input.layout().grid());
  solver::FaceField stage3(input.normal_axis(), input.layout().grid());

  copy_storage_values(input, stage0);
  apply_component_boundary_conditions(boundary_conditions, stage0);
  apply_directional_helmholtz_operator(stage0, alpha, solver::Axis::x, stage1);
  apply_component_boundary_conditions(boundary_conditions, stage1);
  apply_directional_helmholtz_operator(stage1, alpha, solver::Axis::y, stage2);
  apply_component_boundary_conditions(boundary_conditions, stage2);
  apply_directional_helmholtz_operator(stage2, alpha, solver::Axis::z, stage3);
  apply_component_boundary_conditions(boundary_conditions, stage3);
  return stage3;
}

std::vector<double> solve_dense_linear_system(std::vector<std::vector<double>> matrix,
                                              std::vector<double> rhs) {
  const std::size_t n = rhs.size();
  require(matrix.size() == n, "dense solve expects a square matrix");

  for(std::size_t pivot = 0; pivot < n; ++pivot) {
    std::size_t best = pivot;
    double max_value = std::abs(matrix[pivot][pivot]);
    for(std::size_t row = pivot + 1; row < n; ++row) {
      const double candidate = std::abs(matrix[row][pivot]);
      if(candidate > max_value) {
        max_value = candidate;
        best = row;
      }
    }

    require(max_value > 1.0e-14, "dense solve encountered a singular matrix");
    if(best != pivot) {
      std::swap(matrix[pivot], matrix[best]);
      std::swap(rhs[pivot], rhs[best]);
    }

    const double diagonal = matrix[pivot][pivot];
    for(std::size_t column = pivot; column < n; ++column) {
      matrix[pivot][column] /= diagonal;
    }
    rhs[pivot] /= diagonal;

    for(std::size_t row = 0; row < n; ++row) {
      if(row == pivot) {
        continue;
      }

      const double factor = matrix[row][pivot];
      if(std::abs(factor) <= 1.0e-14) {
        continue;
      }
      for(std::size_t column = pivot; column < n; ++column) {
        matrix[row][column] -= factor * matrix[pivot][column];
      }
      rhs[row] -= factor * rhs[pivot];
    }
  }

  return rhs;
}

double mixed_second_derivative_2d(const solver::StructuredField& input,
                                  const int i,
                                  const int j,
                                  const int k) {
  const double dy = input.layout().grid().dy;
  const double dy_squared = dy * dy;
  const double dyy_im1 =
      (input(i - 1, j + 1, k) - 2.0 * input(i - 1, j, k) + input(i - 1, j - 1, k)) / dy_squared;
  const double dyy_i =
      (input(i, j + 1, k) - 2.0 * input(i, j, k) + input(i, j - 1, k)) / dy_squared;
  const double dyy_ip1 =
      (input(i + 1, j + 1, k) - 2.0 * input(i + 1, j, k) + input(i + 1, j - 1, k)) / dy_squared;

  const double dx = input.layout().grid().dx;
  const double inverse_dx_squared = 1.0 / (dx * dx);
  return (dyy_ip1 - 2.0 * dyy_i + dyy_im1) * inverse_dx_squared;
}

double velocity_max_abs(const solver::VelocityField& velocity) {
  return std::max({max_abs_active(velocity.x), max_abs_active(velocity.y), max_abs_active(velocity.z)});
}

void test_build_profile_is_locked() {
  const solver::BuildInfo build_info = solver::get_build_info();
  require(build_info.build_profile == "deterministic" || build_info.build_profile == "benchmark",
          "unexpected build profile");
}

void test_runtime_platform_is_supported() {
  const solver::BuildInfo build_info = solver::get_build_info();
  require(build_info.supported_runtime_platform, "expected Apple Silicon runtime");
}

void test_project_name_matches_orchard_flow() {
  const solver::BuildInfo build_info = solver::get_build_info();
  require(build_info.project_name == "Orchard Flow", "unexpected project name");
}

void test_banner_contains_profile() {
  const solver::BuildInfo build_info = solver::get_build_info();
  const std::string banner = solver::format_build_banner(build_info);
  require(banner.find("profile: " + build_info.build_profile) != std::string::npos,
          "banner missing build profile");
}

void test_grid_coordinates() {
  const solver::Grid grid{4, 3, 2, 0.5, 0.25, 0.125, 1};

  require(!grid.is_2d(), "expected 3D grid");
  require(grid.cells(solver::Axis::x) == 4, "wrong x cell count");
  require(grid.cell_center(solver::Axis::x, 0) == 0.25, "wrong x cell center");
  require(grid.cell_center(solver::Axis::y, 1) == 0.375, "wrong y cell center");
  require(grid.face_coordinate(solver::Axis::z, 2) == 0.25, "wrong z face coordinate");
}

void test_pressure_layout_indexing() {
  const solver::Grid grid{4, 3, 2, 0.5, 0.25, 0.125, 1};
  const solver::PressureField pressure{grid};
  const solver::FieldLayout& layout = pressure.layout();
  const solver::Extent3D active = layout.active_extent();
  const solver::Extent3D storage = layout.storage_extent();
  const solver::IndexRange3D active_range = layout.active_range();

  require(active.nx == 4 && active.ny == 3 && active.nz == 2, "wrong pressure active extent");
  require(storage.nx == 6 && storage.ny == 5 && storage.nz == 4, "wrong pressure storage extent");
  require(active_range.i_begin == 1 && active_range.i_end == 5, "wrong active i range");
  require(active_range.j_begin == 1 && active_range.j_end == 4, "wrong active j range");
  require(active_range.k_begin == 1 && active_range.k_end == 3, "wrong active k range");
  require(layout.is_unit_stride_i(), "expected i to be unit stride");

  const std::size_t expected =
      static_cast<std::size_t>(2) + static_cast<std::size_t>(storage.nx) *
                                       (static_cast<std::size_t>(3) +
                                        static_cast<std::size_t>(storage.ny) *
                                            static_cast<std::size_t>(1));
  require(layout.index(2, 3, 1) == expected, "flat index mapping does not match contract");
}

void test_ghost_cell_access_and_boundary_ranges() {
  const solver::Grid grid{4, 3, 2, 0.5, 0.25, 0.125, 1};
  solver::PressureField pressure{grid};
  pressure.fill(0.0);
  pressure.fill_ghost_layer(solver::BoundaryFace::x_min, 0, 42.0);

  const solver::IndexRange3D x_min_ghost = pressure.layout().ghost_range(solver::BoundaryFace::x_min);
  const solver::IndexRange3D x_max_boundary =
      pressure.layout().boundary_active_range(solver::BoundaryFace::x_max);

  for(int k = x_min_ghost.k_begin; k < x_min_ghost.k_end; ++k) {
    for(int j = x_min_ghost.j_begin; j < x_min_ghost.j_end; ++j) {
      for(int i = x_min_ghost.i_begin; i < x_min_ghost.i_end; ++i) {
        require(pressure(i, j, k) == 42.0, "ghost fill helper missed a cell");
      }
    }
  }

  require(x_max_boundary.i_begin == pressure.layout().active_range().i_end - 1,
          "wrong x_max boundary slab");
  require(x_max_boundary.i_end == pressure.layout().active_range().i_end,
          "wrong x_max boundary slab end");

  const solver::Index3D first_active = pressure.layout().storage_index_from_active(0, 0, 0);
  require(pressure(first_active.i, first_active.j, first_active.k) == 0.0,
          "ghost fill touched active storage");
}

void test_memory_layout_and_alignment() {
  const solver::Grid grid{4, 3, 2, 0.5, 0.25, 0.125, 1};
  solver::PressureField pressure{grid};
  const solver::IndexRange3D active_range = pressure.layout().active_range();

  require(pressure.is_aligned(), "pressure storage is not aligned");
  require(&pressure(active_range.i_begin + 1, active_range.j_begin, active_range.k_begin) ==
              &pressure(active_range.i_begin, active_range.j_begin, active_range.k_begin) + 1,
          "i-direction is not contiguous");
}

void test_cell_and_face_placement() {
  const solver::Grid grid{4, 3, 2, 0.5, 0.25, 0.125, 1};
  const solver::PressureField pressure{grid};
  const solver::VelocityField velocity{grid};

  require(pressure.layout().active_extent().nx == 4, "wrong pressure x extent");
  require(velocity.x.layout().active_extent().nx == 5, "wrong u extent");
  require(velocity.y.layout().active_extent().ny == 4, "wrong v extent");
  require(velocity.z.layout().active_extent().nz == 3, "wrong w extent");

  require(pressure.layout().coordinate_at_active_index(solver::Axis::x, 0) == 0.25,
          "wrong pressure x coordinate");
  require(velocity.x.layout().coordinate_at_active_index(solver::Axis::x, 0) == 0.0,
          "wrong u-face x coordinate");
  require(velocity.x.layout().coordinate_at_active_index(solver::Axis::y, 0) == 0.125,
          "wrong u-face y coordinate");
  require(velocity.y.layout().coordinate_at_active_index(solver::Axis::y, 0) == 0.0,
          "wrong v-face y coordinate");
  require(velocity.z.layout().coordinate_at_active_index(solver::Axis::z, 0) == 0.0,
          "wrong w-face z coordinate");
}

void test_double_precision_storage_contract() {
  static_assert(std::is_same_v<solver::StructuredField::value_type, double>,
                "structured fields must store doubles");
  static_assert(std::is_same_v<solver::PressureField::value_type, double>,
                "pressure field must store doubles");
  static_assert(std::is_same_v<solver::ScalarField::value_type, double>,
                "scalar field must store doubles");

  const solver::Grid grid{4, 3, 2, 0.5, 0.25, 0.125, 1};
  const solver::PressureField pressure{grid};

  require(sizeof(*pressure.data()) == sizeof(double), "unexpected pressure storage precision");
}

void test_advection_options_and_cfl_diagnostic() {
  const solver::AdvectionOptions options{};
  require(solver::to_string(options.scheme) == "tvd", "wrong default advection scheme label");
  require(solver::to_string(options.limiter) == "van_leer", "wrong default limiter label");
  require(solver::describe(options).find("scheme=tvd") != std::string::npos,
          "missing scheme description");
  require(solver::describe(options).find("limiter=van_leer") != std::string::npos,
          "missing limiter description");

  const solver::Grid grid{8, 4, 1, 0.25, 0.5, 1.0, 1};
  solver::VelocityField velocity{grid};
  fill_storage(velocity.x, [](double, double, double) { return 2.0; });
  fill_storage(velocity.y, [](double, double, double) { return -1.0; });
  fill_storage(velocity.z, [](double, double, double) { return 0.0; });

  const double dt = 0.05;
  const solver::CflDiagnostics diagnostics = solver::compute_advective_cfl(velocity, dt);
  require(std::abs(diagnostics.max_u - 2.0) < 1.0e-12, "wrong max u in CFL diagnostic");
  require(std::abs(diagnostics.max_v - 1.0) < 1.0e-12, "wrong max v in CFL diagnostic");
  require(std::abs(diagnostics.max_w) < 1.0e-12, "wrong max w in CFL diagnostic");
  require(std::abs(diagnostics.max_cfl - dt * (2.0 / 0.25 + 1.0 / 0.5)) < 1.0e-12,
          "wrong CFL value");
}

void test_projection_boundary_mapping() {
  solver::BoundaryConditionSet boundary_conditions{};
  boundary_conditions[solver::BoundaryFace::x_min].type = solver::PhysicalBoundaryType::no_slip_wall;
  boundary_conditions[solver::BoundaryFace::x_max].type =
      solver::PhysicalBoundaryType::prescribed_velocity;
  boundary_conditions[solver::BoundaryFace::y_min].type = solver::PhysicalBoundaryType::symmetry;
  boundary_conditions[solver::BoundaryFace::y_max].type = solver::PhysicalBoundaryType::fixed_pressure;
  boundary_conditions[solver::BoundaryFace::y_max].pressure = 1.25;
  boundary_conditions[solver::BoundaryFace::z_min].type = solver::PhysicalBoundaryType::periodic;
  boundary_conditions[solver::BoundaryFace::z_max].type = solver::PhysicalBoundaryType::periodic;

  const solver::PressureBoundarySet mapped =
      solver::derive_pressure_correction_boundary_conditions(boundary_conditions);

  require(mapped[solver::BoundaryFace::x_min].type == solver::PressureBoundaryType::neumann,
          "no-slip wall should map to pressure Neumann");
  require(mapped[solver::BoundaryFace::x_max].type == solver::PressureBoundaryType::neumann,
          "prescribed velocity should map to pressure Neumann");
  require(mapped[solver::BoundaryFace::y_min].type == solver::PressureBoundaryType::neumann,
          "symmetry should map to pressure Neumann");
  require(mapped[solver::BoundaryFace::y_max].type == solver::PressureBoundaryType::dirichlet,
          "fixed pressure should map to pressure Dirichlet");
  require(std::abs(mapped[solver::BoundaryFace::y_max].value) < 1.0e-12,
          "pressure correction Dirichlet boundary should be homogeneous");
  require(mapped[solver::BoundaryFace::z_min].type == solver::PressureBoundaryType::periodic,
          "periodic boundary should map to periodic pressure");
  require(mapped[solver::BoundaryFace::z_max].type == solver::PressureBoundaryType::periodic,
          "periodic boundary should map to periodic pressure");
}

void test_predictor_adi_preserves_quiescent_state() {
  const solver::Grid grid{12, 10, 1, 1.0 / 12.0, 1.0 / 10.0, 1.0, 1};
  const solver::BoundaryConditionSet boundary_conditions = solver::BoundaryConditionSet::cavity();

  solver::VelocityField rhs{grid};
  solver::VelocityField predicted{grid};
  rhs.fill(0.0);
  predicted.fill(1.0);

  const solver::HelmholtzDiagnostics diagnostics =
      solver::solve_predictor_adi(rhs, 0.05, boundary_conditions, predicted);

  require(diagnostics.line_solves > 0, "ADI predictor should perform deterministic line solves");
  require(velocity_max_abs(predicted) <= 1.0e-12, "quiescent predictor state should remain zero");
}

void test_predictor_adi_preserves_periodic_couette_state() {
  solver::ChannelFlowConfig config = solver::default_channel_flow_config();
  config.case_kind = solver::ChannelFlowCase::couette;
  config.nx = 16;
  config.ny = 12;
  const solver::Grid grid{config.nx,
                          config.ny,
                          1,
                          1.0 / static_cast<double>(config.nx),
                          1.0 / static_cast<double>(config.ny),
                          1.0,
                          1};
  const solver::BoundaryConditionSet boundary_conditions =
      solver::make_channel_flow_boundary_conditions(config);

  solver::VelocityField rhs{grid};
  solver::VelocityField predicted{grid};
  fill_storage(rhs.x, [&config](double, const double y, double) {
    return config.top_velocity * y;
  });
  fill_storage(rhs.y, [](double, double, double) {
    return 0.0;
  });
  fill_storage(rhs.z, [](double, double, double) {
    return 0.0;
  });
  solver::apply_velocity_boundary_conditions(boundary_conditions, rhs);
  predicted = rhs;

  const solver::HelmholtzDiagnostics diagnostics =
      solver::solve_predictor_adi(rhs, 0.01, boundary_conditions, predicted);

  require(diagnostics.line_solves > 0, "periodic Couette predictor should perform line solves");
  require(relative_active_l2_difference(predicted.x, rhs.x) <= 1.0e-12,
          "periodic ADI sweep should preserve the exact Couette streamwise state");
  require(velocity_max_abs(predicted) <= config.top_velocity + 1.0e-12,
          "periodic Couette predictor should remain bounded");
}

void test_total_pressure_boundary_conditions_general_bc() {
  const solver::Grid grid{4, 3, 1, 0.25, 1.0 / 3.0, 1.0, 1};
  solver::BoundaryConditionSet boundary_conditions{};
  boundary_conditions[solver::BoundaryFace::x_min].type =
      solver::PhysicalBoundaryType::fixed_pressure;
  boundary_conditions[solver::BoundaryFace::x_min].pressure = 2.5;
  boundary_conditions[solver::BoundaryFace::x_max].type =
      solver::PhysicalBoundaryType::fixed_pressure;
  boundary_conditions[solver::BoundaryFace::x_max].pressure = 0.75;
  boundary_conditions[solver::BoundaryFace::y_min].type =
      solver::PhysicalBoundaryType::no_slip_wall;
  boundary_conditions[solver::BoundaryFace::y_max].type =
      solver::PhysicalBoundaryType::prescribed_velocity;
  boundary_conditions[solver::BoundaryFace::y_max].velocity = {1.0, 0.0, 0.0};
  boundary_conditions[solver::BoundaryFace::z_min].type = solver::PhysicalBoundaryType::symmetry;
  boundary_conditions[solver::BoundaryFace::z_max].type = solver::PhysicalBoundaryType::symmetry;

  solver::PressureField pressure_total{grid};
  solver::VelocityField diffusion{grid};
  pressure_total.fill(-99.0);
  diffusion.fill(0.0);

  fill_storage(pressure_total, [](const double x, const double y, double) {
    return 1.0 + 2.0 * x - 3.0 * y;
  });

  const solver::IndexRange3D pressure_active = pressure_total.layout().active_range();
  const solver::IndexRange3D u_active = diffusion.x.layout().active_range();
  const solver::IndexRange3D v_active = diffusion.y.layout().active_range();

  for(int k = u_active.k_begin; k < u_active.k_end; ++k) {
    for(int j = u_active.j_begin; j < u_active.j_end; ++j) {
      diffusion.x(u_active.i_begin, j, k) = 0.5 + 0.1 * static_cast<double>(j);
      diffusion.x(u_active.i_end - 1, j, k) = -0.25 + 0.05 * static_cast<double>(j);
    }
  }
  for(int k = v_active.k_begin; k < v_active.k_end; ++k) {
    for(int i = v_active.i_begin; i < v_active.i_end; ++i) {
      diffusion.y(i, v_active.j_begin, k) = -0.4 + 0.08 * static_cast<double>(i);
      diffusion.y(i, v_active.j_end - 1, k) = 0.3 - 0.07 * static_cast<double>(i);
    }
  }

  solver::apply_total_pressure_boundary_conditions(boundary_conditions, diffusion, pressure_total);

  const solver::IndexRange3D x_min_ghost = pressure_total.layout().ghost_range(solver::BoundaryFace::x_min);
  const solver::IndexRange3D x_max_ghost = pressure_total.layout().ghost_range(solver::BoundaryFace::x_max);
  const solver::IndexRange3D y_min_ghost = pressure_total.layout().ghost_range(solver::BoundaryFace::y_min);
  const solver::IndexRange3D y_max_ghost = pressure_total.layout().ghost_range(solver::BoundaryFace::y_max);
  const solver::IndexRange3D z_min_ghost = pressure_total.layout().ghost_range(solver::BoundaryFace::z_min);
  const solver::IndexRange3D z_max_ghost = pressure_total.layout().ghost_range(solver::BoundaryFace::z_max);

  for(int k = pressure_active.k_begin; k < pressure_active.k_end; ++k) {
    for(int j = pressure_active.j_begin; j < pressure_active.j_end; ++j) {
      const double x_min_expected = 2.0 * boundary_conditions[solver::BoundaryFace::x_min].pressure -
                                    pressure_total(pressure_active.i_begin, j, k);
      const double x_max_expected = 2.0 * boundary_conditions[solver::BoundaryFace::x_max].pressure -
                                    pressure_total(pressure_active.i_end - 1, j, k);
      require(std::abs(pressure_total(x_min_ghost.i_begin, j, k) - x_min_expected) <= 1.0e-12,
              "fixed-pressure x_min total-pressure ghost fill is inconsistent");
      require(std::abs(pressure_total(x_max_ghost.i_begin, j, k) - x_max_expected) <= 1.0e-12,
              "fixed-pressure x_max total-pressure ghost fill is inconsistent");
    }
  }

  for(int k = pressure_active.k_begin; k < pressure_active.k_end; ++k) {
    for(int i = pressure_active.i_begin; i < pressure_active.i_end; ++i) {
      const double y_min_expected =
          pressure_total(i, pressure_active.j_begin, k) -
          diffusion.y(i, v_active.j_begin, k) * grid.dy;
      const double y_max_expected =
          pressure_total(i, pressure_active.j_end - 1, k) +
          diffusion.y(i, v_active.j_end - 1, k) * grid.dy;
      require(std::abs(pressure_total(i, y_min_ghost.j_begin, k) - y_min_expected) <= 1.0e-12,
              "Neumann y_min total-pressure ghost fill is inconsistent");
      require(std::abs(pressure_total(i, y_max_ghost.j_begin, k) - y_max_expected) <= 1.0e-12,
              "Neumann y_max total-pressure ghost fill is inconsistent");
    }
  }

  for(int j = pressure_active.j_begin; j < pressure_active.j_end; ++j) {
    for(int i = pressure_active.i_begin; i < pressure_active.i_end; ++i) {
      require(std::abs(pressure_total(i, j, z_min_ghost.k_begin) -
                       pressure_total(i, j, pressure_active.k_begin)) <= 1.0e-12,
              "z_min total-pressure ghost fill should be zero-gradient");
      require(std::abs(pressure_total(i, j, z_max_ghost.k_begin) -
                       pressure_total(i, j, pressure_active.k_end - 1)) <= 1.0e-12,
              "z_max total-pressure ghost fill should be zero-gradient");
    }
  }
}

void test_predictor_adi_matches_factorized_dense_reference() {
  const solver::Grid grid{4, 3, 1, 0.25, 1.0 / 3.0, 1.0, 1};
  solver::BoundaryConditionSet boundary_conditions =
      solver::BoundaryConditionSet::all(solver::PhysicalBoundaryType::no_slip_wall);
  boundary_conditions[solver::BoundaryFace::z_min].type = solver::PhysicalBoundaryType::symmetry;
  boundary_conditions[solver::BoundaryFace::z_max].type = solver::PhysicalBoundaryType::symmetry;

  solver::VelocityField rhs{grid};
  rhs.fill(0.0);
  fill_storage(rhs.x, [](const double x, const double y, double) {
    return 0.2 + 0.4 * x - 0.3 * y + 0.1 * x * y;
  });
  solver::apply_velocity_boundary_conditions(boundary_conditions, rhs);

  const std::vector<solver::Index3D> unknowns = collect_factorized_unknowns(rhs.x, boundary_conditions);
  require(!unknowns.empty(), "factorized reference test requires interior unknowns");

  const double alpha = 0.0375;
  std::vector<std::vector<double>> matrix(unknowns.size(), std::vector<double>(unknowns.size(), 0.0));
  for(std::size_t column = 0; column < unknowns.size(); ++column) {
    solver::FaceField basis(solver::Axis::x, grid);
    basis.fill(0.0);
    basis(unknowns[column].i, unknowns[column].j, unknowns[column].k) = 1.0;
    const solver::FaceField image = apply_factorized_operator(basis, boundary_conditions, alpha);
    const std::vector<double> values = extract_unknown_values(image, unknowns);
    for(std::size_t row = 0; row < unknowns.size(); ++row) {
      matrix[row][column] = values[row];
    }
  }

  const std::vector<double> rhs_vector = extract_unknown_values(rhs.x, unknowns);
  const std::vector<double> dense_solution = solve_dense_linear_system(matrix, rhs_vector);

  solver::VelocityField predicted{grid};
  predicted.fill(0.0);
  const solver::HelmholtzDiagnostics diagnostics =
      solver::solve_predictor_adi(rhs, alpha, boundary_conditions, predicted);
  const std::vector<double> adi_solution = extract_unknown_values(predicted.x, unknowns);

  require(diagnostics.line_solves > 0, "factorized predictor should perform line solves");
  double max_error = 0.0;
  for(std::size_t index = 0; index < dense_solution.size(); ++index) {
    max_error = std::max(max_error, std::abs(adi_solution[index] - dense_solution[index]));
  }
  require(max_error <= 1.0e-10,
          "ADI predictor drifted from the dense factorized reference solve: max_error=" +
              std::to_string(max_error));
}

void test_lid_driven_cavity_predictor_rhs_matches_manual_formula() {
  solver::LidDrivenCavityConfig config = solver::default_lid_driven_cavity_config();
  config.nx = 8;
  config.ny = 6;
  const solver::Grid grid{config.nx,
                          config.ny,
                          1,
                          1.0 / static_cast<double>(config.nx),
                          1.0 / static_cast<double>(config.ny),
                          1.0,
                          1};
  const solver::BoundaryConditionSet boundary_conditions =
      solver::make_lid_driven_cavity_boundary_conditions(config);
  const double viscosity = config.lid_velocity / config.reynolds;
  const double dt = 0.4 * std::min(grid.dx, grid.dy);
  const double alpha = 0.5 * viscosity * dt;

  solver::VelocityField current_velocity{grid};
  solver::VelocityField advection_previous{grid};
  solver::VelocityField advection_current{grid};
  solver::VelocityField diffusion{grid};
  solver::VelocityField pressure_gradient{grid};
  solver::VelocityField factorized_correction{grid};
  solver::VelocityField predictor_rhs{grid};
  solver::VelocityField expected_advection{grid};
  solver::VelocityField expected_diffusion{grid};
  solver::VelocityField expected_pressure_gradient{grid};
  solver::VelocityField expected_rhs{grid};
  solver::PressureField pressure_total{grid};

  fill_storage(current_velocity.x, [](const double x, const double y, double) {
    return std::sin(pi() * x) * std::sin(pi() * y);
  });
  fill_storage(current_velocity.y, [](const double x, const double y, double) {
    return 0.25 * std::cos(pi() * x) * std::sin(pi() * y);
  });
  fill_storage(current_velocity.z, [](double, double, double) {
    return 0.0;
  });
  fill_storage(advection_previous.x, [](const double x, const double y, double) {
    return 0.1 + x - 0.5 * y;
  });
  fill_storage(advection_previous.y, [](const double x, const double y, double) {
    return -0.2 + 0.3 * x + 0.4 * y;
  });
  fill_storage(advection_previous.z, [](double, double, double) {
    return 0.0;
  });
  fill_storage(pressure_total, [](const double x, const double y, double) {
    return 0.3 * x * x - 0.2 * y + 0.1 * x * y;
  });

  solver::apply_velocity_boundary_conditions(boundary_conditions, current_velocity);
  solver::compute_diffusion_term(current_velocity, viscosity, expected_diffusion);
  solver::apply_total_pressure_boundary_conditions(boundary_conditions, expected_diffusion,
                                                   pressure_total);
  solver::compute_advection_term(current_velocity, config.advection, expected_advection);
  solver::operators::compute_gradient(pressure_total, expected_pressure_gradient);
  expected_rhs = current_velocity;
  axpy_velocity(expected_rhs, expected_pressure_gradient, -dt);
  axpy_velocity(expected_rhs, expected_diffusion, 0.5 * dt);
  axpy_velocity(expected_rhs, expected_advection, -1.5 * dt);
  axpy_velocity(expected_rhs, advection_previous, 0.5 * dt);

  solver::VelocityField expected_factorized_correction{grid};
  expected_factorized_correction.fill(0.0);
  for(int k = current_velocity.x.layout().active_range().k_begin;
      k < current_velocity.x.layout().active_range().k_end;
      ++k) {
    for(int j = current_velocity.x.layout().active_range().j_begin;
        j < current_velocity.x.layout().active_range().j_end;
        ++j) {
      for(int i = current_velocity.x.layout().active_range().i_begin;
          i < current_velocity.x.layout().active_range().i_end;
          ++i) {
        expected_factorized_correction.x(i, j, k) = mixed_second_derivative_2d(current_velocity.x, i, j, k);
      }
    }
  }
  for(int k = current_velocity.y.layout().active_range().k_begin;
      k < current_velocity.y.layout().active_range().k_end;
      ++k) {
    for(int j = current_velocity.y.layout().active_range().j_begin;
        j < current_velocity.y.layout().active_range().j_end;
        ++j) {
      for(int i = current_velocity.y.layout().active_range().i_begin;
          i < current_velocity.y.layout().active_range().i_end;
          ++i) {
        expected_factorized_correction.y(i, j, k) = mixed_second_derivative_2d(current_velocity.y, i, j, k);
      }
    }
  }
  axpy_velocity(expected_rhs, expected_factorized_correction, alpha * alpha);

  solver::detail::assemble_lid_driven_cavity_predictor_rhs(current_velocity,
                                                           pressure_total,
                                                           &advection_previous,
                                                           config.advection,
                                                           viscosity,
                                                           dt,
                                                           advection_current,
                                                           diffusion,
                                                           pressure_gradient,
                                                           factorized_correction,
                                                           predictor_rhs);

  require(active_l2_difference(advection_current.x, expected_advection.x) <= 1.0e-12,
          "predictor assembly should return the current advection term");
  require(active_l2_difference(diffusion.x, expected_diffusion.x) <= 1.0e-12,
          "predictor assembly should return the current diffusion term");
  require(active_l2_difference(pressure_gradient.x, expected_pressure_gradient.x) <= 1.0e-12,
          "predictor assembly should use the total-pressure gradient");
  require(max_abs_active(expected_factorized_correction.x) > 1.0e-4,
          "predictor correction test needs a nontrivial factorization correction");
  require(active_l2_difference(predictor_rhs.x, expected_rhs.x) <= 1.0e-12,
          "predictor RHS x component missed the factorized correction term");
  require(active_l2_difference(predictor_rhs.y, expected_rhs.y) <= 1.0e-12,
          "predictor RHS y component missed the factorized correction term");
}

void test_poisson_mgpcg_discrete_dirichlet_recovery() {
  const int resolution = 32;
  const solver::Grid grid{resolution,
                          resolution,
                          1,
                          1.0 / static_cast<double>(resolution),
                          1.0 / static_cast<double>(resolution),
                          1.0,
                          1};
  solver::PressureBoundarySet boundary_conditions{};
  boundary_conditions[solver::BoundaryFace::x_min].type = solver::PressureBoundaryType::dirichlet;
  boundary_conditions[solver::BoundaryFace::x_max].type = solver::PressureBoundaryType::dirichlet;
  boundary_conditions[solver::BoundaryFace::y_min].type = solver::PressureBoundaryType::dirichlet;
  boundary_conditions[solver::BoundaryFace::y_max].type = solver::PressureBoundaryType::dirichlet;
  boundary_conditions[solver::BoundaryFace::z_min].type = solver::PressureBoundaryType::neumann;
  boundary_conditions[solver::BoundaryFace::z_max].type = solver::PressureBoundaryType::neumann;

  solver::PressureField exact_pressure{grid};
  solver::ScalarField rhs{grid};
  solver::PressureField solution{grid};
  fill_storage(exact_pressure, [](const double x, const double y, double) {
    return std::sin(pi() * x) * std::sin(pi() * y);
  });

  solver::linsolve::build_poisson_rhs_from_pressure(exact_pressure, boundary_conditions, rhs);
  solution.fill(0.0);

  const solver::ProjectionOptions options{
      .dt = 1.0,
      .density = 1.0,
      .poisson_max_iterations = 200,
      .poisson_tolerance = 1.0e-10,
  };
  const solver::PoissonSolveDiagnostics diagnostics =
      solver::linsolve::solve_pressure_poisson(rhs, boundary_conditions, options, solution);

  const double relative_error = relative_active_l2_difference(solution, exact_pressure);
  require(diagnostics.converged, "MGPCG Dirichlet solve should converge");
  require(diagnostics.relative_residual <= 1.0e-10, "MGPCG relative residual is too large");
  require(relative_error <= 1.0e-8, "MGPCG Dirichlet recovery error is too large");
  require(diagnostics.multigrid_levels >= 2, "MGPCG should build a multigrid hierarchy");
  require(diagnostics.coarse_unknowns > 0, "MGPCG should report the coarse-grid size");
  require(diagnostics.solver == "mgpcg", "unexpected pressure solver label");
  require(diagnostics.preconditioner == "geometric_multigrid", "unexpected preconditioner label");
  require(diagnostics.cycle == "v_cycle", "unexpected multigrid cycle label");
  require(diagnostics.smoother == "damped_jacobi", "unexpected smoother label");
  require(diagnostics.pre_smoothing_steps == 2, "unexpected pre-smoothing policy");
  require(diagnostics.post_smoothing_steps == 2, "unexpected post-smoothing policy");
}

void test_poisson_mgpcg_pure_neumann_zero_mean_recovery() {
  const int resolution = 24;
  const solver::Grid grid{resolution,
                          resolution,
                          1,
                          1.0 / static_cast<double>(resolution),
                          1.0 / static_cast<double>(resolution),
                          1.0,
                          1};
  solver::PressureBoundarySet boundary_conditions{};
  for(solver::PressureBoundaryCondition& boundary : boundary_conditions.faces) {
    boundary.type = solver::PressureBoundaryType::neumann;
  }

  solver::PressureField exact_pressure{grid};
  solver::ScalarField rhs{grid};
  solver::PressureField solution{grid};
  fill_storage(exact_pressure, [](const double x, const double y, double) {
    return std::cos(2.0 * pi() * x) * std::cos(2.0 * pi() * y);
  });
  subtract_active_mean_inplace(exact_pressure);

  solver::linsolve::build_poisson_rhs_from_pressure(exact_pressure, boundary_conditions, rhs);
  solution.fill(0.0);

  const solver::ProjectionOptions options{
      .dt = 1.0,
      .density = 1.0,
      .poisson_max_iterations = 300,
      .poisson_tolerance = 1.0e-10,
  };
  const solver::PoissonSolveDiagnostics diagnostics =
      solver::linsolve::solve_pressure_poisson(rhs, boundary_conditions, options, solution);

  const double relative_error = relative_active_l2_difference(solution, exact_pressure);
  require(diagnostics.converged, "MGPCG pure-Neumann solve should converge");
  require(diagnostics.zero_mean_enforced, "pure-Neumann solve should enforce zero mean");
  require(std::abs(active_mean_value(solution)) <= 1.0e-10,
          "pure-Neumann solve should preserve zero-mean pressure");
  require(diagnostics.relative_residual <= 1.0e-10, "pure-Neumann relative residual is too large");
  require(relative_error <= 1.0e-8, "pure-Neumann recovery error is too large");
}

struct ManufacturedErrors {
  double gradient;
  double divergence;
  double laplacian;
};

ManufacturedErrors run_manufactured_solution_case(const int resolution) {
  const solver::OperatorManufacturedSolutionResult result =
      solver::run_operator_manufactured_solution_case(resolution);
  return ManufacturedErrors{
      .gradient = result.gradient_error,
      .divergence = result.divergence_error,
      .laplacian = result.laplacian_error,
  };
}

void test_manufactured_solution_convergence() {
  const ManufacturedErrors coarse = run_manufactured_solution_case(16);
  const ManufacturedErrors fine = run_manufactured_solution_case(32);

  require(coarse.gradient > fine.gradient, "gradient error did not decrease on refinement");
  require(coarse.divergence > fine.divergence, "divergence error did not decrease on refinement");
  require(coarse.laplacian > fine.laplacian, "laplacian error did not decrease on refinement");

  require(observed_order(coarse.gradient, fine.gradient) >= 1.8,
          "gradient convergence order below second-order target");
  require(observed_order(coarse.divergence, fine.divergence) >= 1.8,
          "divergence convergence order below second-order target");
  require(observed_order(coarse.laplacian, fine.laplacian) >= 1.8,
          "laplacian convergence order below second-order target");
}

void test_diffusion_term_matches_scaled_laplacian() {
  const double domain_length = 2.0 * pi();
  const int resolution = 32;
  const double spacing = domain_length / static_cast<double>(resolution);
  const solver::Grid grid{resolution, resolution, 1, spacing, spacing, 1.0, 1};
  const double viscosity = 0.125;

  solver::VelocityField velocity{grid};
  solver::VelocityField diffusion{grid};

  fill_storage(velocity.x, [](const double x, const double y, double) {
    return -std::cos(x) * std::sin(y);
  });
  fill_storage(velocity.y, [](const double x, const double y, double) {
    return std::sin(x) * std::cos(y);
  });
  fill_storage(velocity.z, [](double, double, double) {
    return 0.0;
  });

  solver::compute_diffusion_term(velocity, viscosity, diffusion);

  const double u_error = active_l2_error(diffusion.x, [viscosity](const double x, const double y, double) {
    return 2.0 * viscosity * std::cos(x) * std::sin(y);
  });
  const double v_error = active_l2_error(diffusion.y, [viscosity](const double x, const double y, double) {
    return -2.0 * viscosity * std::sin(x) * std::cos(y);
  });

  require(u_error < 2.0e-2, "diffusion u error too large");
  require(v_error < 2.0e-2, "diffusion v error too large");
}

struct TaylorGreenStepErrors {
  double velocity_error;
  double energy_error;
};

TaylorGreenStepErrors run_taylor_green_step_case(const int resolution) {
  const double domain_length = 2.0 * pi();
  const double viscosity = 0.01;
  const double spacing = domain_length / static_cast<double>(resolution);
  const solver::Grid grid{resolution, resolution, 1, spacing, spacing, 1.0, 1};
  const double dt = 0.1 * spacing * spacing;
  const solver::AdvectionOptions options{};

  solver::VelocityField velocity{grid};
  solver::VelocityField advection{grid};
  solver::VelocityField diffusion{grid};
  solver::VelocityField pressure_gradient{grid};
  solver::VelocityField updated_velocity{grid};
  solver::PressureField pressure{grid};

  fill_storage(velocity.x, [](const double x, const double y, double) {
    return -std::cos(x) * std::sin(y);
  });
  fill_storage(velocity.y, [](const double x, const double y, double) {
    return std::sin(x) * std::cos(y);
  });
  fill_storage(velocity.z, [](double, double, double) {
    return 0.0;
  });
  fill_storage(pressure, [](const double x, const double y, double) {
    return -0.25 * (std::cos(2.0 * x) + std::cos(2.0 * y));
  });

  updated_velocity = velocity;

  solver::compute_advection_term(velocity, options, advection);
  solver::compute_diffusion_term(velocity, viscosity, diffusion);
  solver::operators::compute_gradient(pressure, pressure_gradient);

  axpy_velocity(updated_velocity, advection, -dt);
  axpy_velocity(updated_velocity, pressure_gradient, -dt);
  axpy_velocity(updated_velocity, diffusion, dt);

  const double decay_velocity = std::exp(-2.0 * viscosity * dt);
  const double exact_energy = 0.25 * std::exp(-4.0 * viscosity * dt);
  const double numerical_energy = kinetic_energy(updated_velocity);

  const double u_error = active_l2_error(updated_velocity.x, [decay_velocity](const double x,
                                                                              const double y,
                                                                              double) {
    return -std::cos(x) * std::sin(y) * decay_velocity;
  });
  const double v_error = active_l2_error(updated_velocity.y, [decay_velocity](const double x,
                                                                              const double y,
                                                                              double) {
    return std::sin(x) * std::cos(y) * decay_velocity;
  });

  return TaylorGreenStepErrors{
      .velocity_error = std::sqrt(0.5 * (square(u_error) + square(v_error))),
      .energy_error = std::abs(numerical_energy - exact_energy),
  };
}

void test_taylor_green_step_behavior() {
  const TaylorGreenStepErrors coarse = run_taylor_green_step_case(16);
  const TaylorGreenStepErrors fine = run_taylor_green_step_case(32);

  require(coarse.velocity_error > fine.velocity_error,
          "Taylor-Green velocity error did not decrease on refinement");
  require(coarse.energy_error > fine.energy_error,
          "Taylor-Green energy error did not decrease on refinement");
  require(observed_order(coarse.velocity_error, fine.velocity_error) >= 1.5,
          "Taylor-Green velocity error did not show expected convergence");
  require(observed_order(coarse.energy_error, fine.energy_error) >= 1.5,
          "Taylor-Green energy error did not show expected convergence");
}

void test_bounded_advection_regression_case() {
  const int resolution_x = 64;
  const int resolution_y = 16;
  const double length = 1.0;
  const double dx = length / static_cast<double>(resolution_x);
  const double dy = length / static_cast<double>(resolution_y);
  const solver::Grid grid{resolution_x, resolution_y, 1, dx, dy, 1.0, 1};

  solver::VelocityField velocity{grid};
  solver::VelocityField advection{grid};
  solver::VelocityField updated{grid};
  solver::AdvectionOptions options{};
  const double dt = 0.4 * dx;

  fill_storage(velocity.x, [](double, double, double) {
    return 1.0;
  });
  fill_storage(velocity.y, [](double, double, double) {
    return 0.0;
  });
  fill_storage(velocity.z, [length](const double x, double, double) {
    const double wrapped_x = wrap_periodic(x, length);
    return wrapped_x < 0.5 ? 1.0 : 0.0;
  });

  updated = velocity;
  solver::compute_advection_term(velocity, options, advection);
  axpy_active(updated.z, advection.z, -dt);

  require(active_min(updated.z) >= -1.0e-12, "bounded TVD update undershot the minimum");
  require(active_max(updated.z) <= 1.0 + 1.0e-12, "bounded TVD update overshot the maximum");
}

void test_static_projection_preserves_quiescent_fluid() {
  const solver::Grid grid{16, 12, 1, 1.0 / 16.0, 1.0 / 12.0, 1.0, 1};
  const solver::BoundaryConditionSet boundary_conditions = solver::BoundaryConditionSet::cavity();
  const solver::ProjectionOptions options{
      .dt = 0.05,
      .density = 1.0,
      .poisson_max_iterations = 2000,
      .poisson_tolerance = 1.0e-12,
  };

  solver::VelocityField predicted{grid};
  solver::VelocityField corrected{grid};
  solver::PressureField pressure{grid};
  solver::ScalarField pressure_rhs{grid};
  predicted.fill(0.0);
  corrected.fill(0.0);
  pressure.fill(0.0);

  const solver::ProjectionDiagnostics diagnostics =
      solver::project_velocity(predicted, boundary_conditions, options, pressure, corrected, &pressure_rhs);

  require(diagnostics.pressure_solve.converged, "zero-RHS pressure solve should converge immediately");
  require(diagnostics.pressure_solve.iterations == 0, "quiescent projection should not iterate");
  require(diagnostics.rhs_l2 <= 1.0e-14, "quiescent projection RHS should remain zero");
  require(diagnostics.divergence_l2_before <= 1.0e-14,
          "quiescent predicted field should already be divergence free");
  require(diagnostics.divergence_l2_after <= 1.0e-14,
          "quiescent corrected field should remain divergence free");
  require(std::abs(diagnostics.pressure_mean) <= 1.0e-14, "quiescent pressure should remain zero");
  require(velocity_max_abs(corrected) <= 1.0e-12, "quiescent corrected velocity should remain zero");
}

void test_pure_neumann_projection_recovers_zero_mean_pressure() {
  const int resolution = 24;
  const solver::Grid grid{resolution, resolution, 1,
                          1.0 / static_cast<double>(resolution),
                          1.0 / static_cast<double>(resolution),
                          1.0,
                          1};
  const solver::BoundaryConditionSet boundary_conditions =
      solver::BoundaryConditionSet::all(solver::PhysicalBoundaryType::symmetry);
  const solver::ProjectionOptions options{
      .dt = 0.025,
      .density = 1.0,
      .poisson_max_iterations = 4000,
      .poisson_tolerance = 1.0e-12,
  };
  const solver::PressureBoundarySet pressure_boundaries =
      solver::derive_pressure_correction_boundary_conditions(boundary_conditions);

  solver::PressureField expected_pressure{grid};
  fill_storage(expected_pressure, [](const double x, const double y, double) {
    const double sx = std::sin(pi() * x);
    const double sy = std::sin(pi() * y);
    return sx * sx * sy * sy;
  });
  solver::apply_pressure_boundary_conditions(pressure_boundaries, expected_pressure);

  solver::VelocityField predicted{grid};
  solver::operators::compute_gradient(expected_pressure, predicted);
  scale_velocity_active(predicted, options.dt);
  solver::apply_velocity_boundary_conditions(boundary_conditions, predicted);

  solver::PressureField pressure{grid};
  solver::VelocityField corrected{grid};
  pressure.fill(0.0);
  corrected.fill(0.0);

  const solver::ProjectionDiagnostics diagnostics =
      solver::project_velocity(predicted, boundary_conditions, options, pressure, corrected);

  const double expected_mean = active_mean_value(expected_pressure);
  const solver::IndexRange3D active = expected_pressure.layout().active_range();
  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        expected_pressure(i, j, k) -= expected_mean;
      }
    }
  }
  solver::apply_pressure_boundary_conditions(pressure_boundaries, expected_pressure);

  const double pressure_error = active_l2_difference(pressure, expected_pressure);
  require(diagnostics.pressure_solve.converged, "pure-Neumann pressure solve should converge");
  require(diagnostics.pressure_solve.zero_mean_enforced,
          "pure-Neumann pressure solve should enforce zero mean");
  require(std::abs(diagnostics.pressure_mean) <= 1.0e-10,
          "pure-Neumann solve should remove the pressure null space");
  require(pressure_error <= 1.0e-8, "pressure recovery error is too large");
  require(velocity_max_abs(corrected) <= 1.0e-8, "projection should remove the gradient field");
  require(diagnostics.divergence_l2_after <= 1.0e-10,
          "projection should leave a divergence-free corrected field");
}

void test_lid_driven_cavity_config_loader_and_bc_subset() {
  const solver::LidDrivenCavityConfig config = solver::load_lid_driven_cavity_config(
      source_path("benchmarks/lid_driven_cavity_smoke.cfg"));
  require(config.nx == 32 && config.ny == 32, "smoke cavity config should load the grid size");
  require(std::abs(config.reynolds - 100.0) < 1.0e-12, "smoke cavity Reynolds number mismatch");
  require(!config.validate_reference, "smoke config should skip reference validation");

  const std::string summary = solver::describe(config);
  require(summary.find("validate_reference=false") != std::string::npos,
          "config description should report validation mode");

  const solver::BoundaryConditionSet boundary_conditions =
      solver::make_lid_driven_cavity_boundary_conditions(config);
  require(boundary_conditions[solver::BoundaryFace::x_min].type ==
              solver::PhysicalBoundaryType::no_slip_wall,
          "cavity x_min should remain a no-slip wall");
  require(boundary_conditions[solver::BoundaryFace::x_max].type ==
              solver::PhysicalBoundaryType::no_slip_wall,
          "cavity x_max should remain a no-slip wall");
  require(boundary_conditions[solver::BoundaryFace::y_min].type ==
              solver::PhysicalBoundaryType::no_slip_wall,
          "cavity y_min should remain a no-slip wall");
  require(boundary_conditions[solver::BoundaryFace::y_max].type ==
              solver::PhysicalBoundaryType::prescribed_velocity,
          "cavity lid should be a prescribed velocity boundary");
  require(std::abs(boundary_conditions[solver::BoundaryFace::y_max].velocity[0] -
                   config.lid_velocity) <= 1.0e-12,
          "cavity lid speed mismatch");
  require(boundary_conditions[solver::BoundaryFace::z_min].type ==
              solver::PhysicalBoundaryType::symmetry,
          "2D cavity z_min should map to symmetry");
  require(boundary_conditions[solver::BoundaryFace::z_max].type ==
              solver::PhysicalBoundaryType::symmetry,
          "2D cavity z_max should map to symmetry");
}

void test_channel_flow_config_loader_and_boundary_conditions() {
  const solver::ChannelFlowConfig couette = solver::load_channel_flow_config(
      source_path("benchmarks/channel_couette_smoke.cfg"));
  require(couette.case_kind == solver::ChannelFlowCase::couette,
          "Couette config should parse the case kind");
  require(couette.nx == 32 && couette.ny == 32, "Couette smoke config should load the grid size");
  require(couette.validate_profile, "Couette smoke config should validate the profile");

  const solver::BoundaryConditionSet couette_bc =
      solver::make_channel_flow_boundary_conditions(couette);
  require(couette_bc[solver::BoundaryFace::x_min].type ==
              solver::PhysicalBoundaryType::periodic,
          "Couette x_min should be periodic");
  require(couette_bc[solver::BoundaryFace::x_max].type ==
              solver::PhysicalBoundaryType::periodic,
          "Couette x_max should be periodic");
  require(couette_bc[solver::BoundaryFace::y_min].type ==
              solver::PhysicalBoundaryType::no_slip_wall,
          "Couette lower wall should be no-slip");
  require(couette_bc[solver::BoundaryFace::y_max].type ==
              solver::PhysicalBoundaryType::prescribed_velocity,
          "Couette upper wall should prescribe the lid velocity");
  require(std::abs(couette_bc[solver::BoundaryFace::y_max].velocity[0] - couette.top_velocity) <=
              1.0e-12,
          "Couette upper-wall speed mismatch");

  const solver::ChannelFlowConfig poiseuille = solver::load_channel_flow_config(
      source_path("benchmarks/channel_poiseuille_smoke.cfg"));
  require(poiseuille.case_kind == solver::ChannelFlowCase::poiseuille,
          "Poiseuille config should parse the case kind");

  const solver::BoundaryConditionSet poiseuille_bc =
      solver::make_channel_flow_boundary_conditions(poiseuille);
  require(poiseuille_bc[solver::BoundaryFace::x_min].type ==
              solver::PhysicalBoundaryType::fixed_pressure,
          "Poiseuille x_min should be fixed pressure");
  require(poiseuille_bc[solver::BoundaryFace::x_max].type ==
              solver::PhysicalBoundaryType::fixed_pressure,
          "Poiseuille x_max should be fixed pressure");
  require(std::abs(poiseuille_bc[solver::BoundaryFace::x_min].pressure -
                   poiseuille.pressure_drop) <= 1.0e-12,
          "Poiseuille inlet pressure mismatch");
  require(std::abs(poiseuille_bc[solver::BoundaryFace::x_max].pressure) <= 1.0e-12,
          "Poiseuille outlet pressure should default to zero");
  require(poiseuille_bc[solver::BoundaryFace::y_max].type ==
              solver::PhysicalBoundaryType::no_slip_wall,
          "Poiseuille upper wall should remain no-slip");
}

void test_channel_flow_couette_profile_validation() {
  const solver::ChannelFlowConfig config = solver::load_channel_flow_config(
      source_path("benchmarks/channel_couette_smoke.cfg"));
  const solver::ChannelFlowResult result = solver::run_channel_flow(config);

  require(result.final_step.step == config.steps, "Couette run should consume the configured steps");
  require(result.final_step.dt > 0.0, "Couette timestep should be positive");
  require(result.final_step.max_cfl <= config.cfl_limit + 1.0e-12,
          "Couette CFL should respect the configured ceiling");
  require(result.validation.reference_dataset == "analytic_couette_profile",
          "Couette validation should report the analytic profile dataset");
  require(result.validation.relative_l2_error <= 5.0e-3,
          "Couette profile error exceeded the M7 threshold");
  require(result.validation.pass, "Couette smoke validation should pass");
}

void test_channel_flow_poiseuille_profile_validation() {
  const solver::ChannelFlowConfig config = solver::load_channel_flow_config(
      source_path("benchmarks/channel_poiseuille_smoke.cfg"));
  const solver::ChannelFlowResult result = solver::run_channel_flow(config);

  require(result.final_step.step == config.steps,
          "Poiseuille run should consume the configured steps");
  require(result.final_step.dt > 0.0, "Poiseuille timestep should be positive");
  require(result.final_step.max_cfl <= config.cfl_limit + 1.0e-12,
          "Poiseuille CFL should respect the configured ceiling");
  require(result.validation.reference_dataset == "analytic_poiseuille_profile",
          "Poiseuille validation should report the analytic profile dataset");
  require(result.validation.relative_l2_error <= 5.0e-3,
          "Poiseuille profile error exceeded the M7 threshold");
  require(result.validation.pass, "Poiseuille smoke validation should pass");
}

void test_channel_flow_rejects_nonconverged_pressure_projection() {
  solver::ChannelFlowConfig config = solver::load_channel_flow_config(
      source_path("benchmarks/channel_poiseuille_smoke.cfg"));
  config.poisson_max_iterations = 1;
  config.poisson_tolerance = 1.0e-30;

  bool rejected = false;
  try {
    static_cast<void>(solver::run_channel_flow(config));
  } catch(const std::exception& exception) {
    rejected = std::string(exception.what()).find("channel flow pressure solve did not converge") !=
               std::string::npos;
  }
  require(rejected, "channel flow should reject a non-converged pressure projection");
}

void test_taylor_green_config_loader_and_boundary_conditions() {
  const solver::TaylorGreenConfig config = solver::load_taylor_green_config(
      source_path("benchmarks/taylor_green_smoke.cfg"));
  require(config.nx == 48 && config.ny == 48 && config.nz == 1,
          "Taylor-Green smoke config should load the 2D grid size");
  require(std::abs(config.viscosity - 0.01) <= 1.0e-12,
          "Taylor-Green smoke viscosity mismatch");
  require(config.validate_energy, "Taylor-Green smoke config should validate the decay");
  require(config.backend == solver::ExecutionBackend::cpu,
          "Taylor-Green configs should default to the CPU backend");

  const solver::BoundaryConditionSet boundary_conditions =
      solver::make_taylor_green_boundary_conditions(config);
  require(boundary_conditions[solver::BoundaryFace::x_min].type ==
              solver::PhysicalBoundaryType::periodic,
          "Taylor-Green x_min should be periodic");
  require(boundary_conditions[solver::BoundaryFace::x_max].type ==
              solver::PhysicalBoundaryType::periodic,
          "Taylor-Green x_max should be periodic");
  require(boundary_conditions[solver::BoundaryFace::y_min].type ==
              solver::PhysicalBoundaryType::periodic,
          "Taylor-Green y_min should be periodic");
  require(boundary_conditions[solver::BoundaryFace::y_max].type ==
              solver::PhysicalBoundaryType::periodic,
          "Taylor-Green y_max should be periodic");
  require(boundary_conditions[solver::BoundaryFace::z_min].type ==
              solver::PhysicalBoundaryType::symmetry,
          "Taylor-Green z_min should map to symmetry");
  require(boundary_conditions[solver::BoundaryFace::z_max].type ==
              solver::PhysicalBoundaryType::symmetry,
          "Taylor-Green z_max should map to symmetry");
}

void test_taylor_green_3d_config_loader_and_boundary_conditions() {
  const solver::TaylorGreenConfig config = solver::load_taylor_green_config(
      source_path("benchmarks/taylor_green_3d_smoke.cfg"));
  require(config.nx == 40 && config.ny == 40 && config.nz == 40,
          "Taylor-Green 3D smoke config should load the full 3D grid size");
  require(config.validate_energy, "Taylor-Green 3D smoke config should validate the decay");
  require(config.backend == solver::ExecutionBackend::cpu,
          "Taylor-Green 3D configs should default to the CPU backend");

  const solver::BoundaryConditionSet boundary_conditions =
      solver::make_taylor_green_boundary_conditions(config);
  require(boundary_conditions[solver::BoundaryFace::x_min].type ==
              solver::PhysicalBoundaryType::periodic,
          "Taylor-Green 3D x_min should be periodic");
  require(boundary_conditions[solver::BoundaryFace::x_max].type ==
              solver::PhysicalBoundaryType::periodic,
          "Taylor-Green 3D x_max should be periodic");
  require(boundary_conditions[solver::BoundaryFace::y_min].type ==
              solver::PhysicalBoundaryType::periodic,
          "Taylor-Green 3D y_min should be periodic");
  require(boundary_conditions[solver::BoundaryFace::y_max].type ==
              solver::PhysicalBoundaryType::periodic,
          "Taylor-Green 3D y_max should be periodic");
  require(boundary_conditions[solver::BoundaryFace::z_min].type ==
              solver::PhysicalBoundaryType::periodic,
          "Taylor-Green 3D z_min should be periodic");
  require(boundary_conditions[solver::BoundaryFace::z_max].type ==
              solver::PhysicalBoundaryType::periodic,
          "Taylor-Green 3D z_max should be periodic");
}

void test_taylor_green_smoke_validation() {
  const solver::TaylorGreenConfig config = solver::load_taylor_green_config(
      source_path("benchmarks/taylor_green_smoke.cfg"));
  const solver::TaylorGreenResult result = solver::run_taylor_green(config);

  require(result.final_step.step > 0, "Taylor-Green smoke run should advance at least one step");
  require(std::abs(result.final_step.time - config.final_time) <= 1.0e-12,
          "Taylor-Green smoke run should land on the configured final time");
  require(result.final_step.max_cfl <= config.cfl_limit + 1.0e-12,
          "Taylor-Green CFL should respect the configured ceiling");
  require(result.final_step.divergence_l2 <= 1.0e-10,
          "Taylor-Green smoke run should remain divergence free");
  require(result.validation.reference_dataset == "analytic_taylor_green_decay_2d",
          "Taylor-Green validation should report the 2D analytic decay dataset");
  require(result.validation.normalized_energy_error <= 1.0e-2,
          "Taylor-Green smoke energy error exceeded the M9 threshold");
  require(result.validation.pass, "Taylor-Green smoke validation should pass");
}

void test_taylor_green_3d_smoke_validation() {
  const solver::TaylorGreenConfig config = solver::load_taylor_green_config(
      source_path("benchmarks/taylor_green_3d_smoke.cfg"));
  const solver::TaylorGreenResult result = solver::run_taylor_green(config);

  require(result.final_step.step > 0, "Taylor-Green 3D smoke run should advance at least one step");
  require(std::abs(result.final_step.time - config.final_time) <= 1.0e-12,
          "Taylor-Green 3D smoke run should land on the configured final time");
  require(result.final_step.max_cfl <= config.cfl_limit + 1.0e-12,
          "Taylor-Green 3D CFL should respect the configured ceiling");
  require(result.final_step.divergence_l2 <= 1.0e-10,
          "Taylor-Green 3D smoke run should remain divergence free");
  require(result.validation.reference_dataset == "analytic_taylor_green_decay_3d",
          "Taylor-Green 3D validation should report the 3D analytic decay dataset");
  require(result.validation.normalized_energy_error <= 1.0e-2,
          "Taylor-Green 3D smoke energy error exceeded the M12 threshold");
  require(result.validation.pass, "Taylor-Green 3D smoke validation should pass");
}

void test_taylor_green_rejects_nonconverged_pressure_projection() {
  solver::TaylorGreenConfig config = solver::load_taylor_green_config(
      source_path("benchmarks/taylor_green_smoke.cfg"));
  config.poisson_max_iterations = 1;
  config.poisson_tolerance = 1.0e-30;

  bool rejected = false;
  try {
    static_cast<void>(solver::run_taylor_green(config));
  } catch(const std::exception& exception) {
    rejected =
        std::string(exception.what()).find("Taylor-Green pressure solve did not converge") !=
        std::string::npos;
  }
  require(rejected, "Taylor-Green should reject a non-converged pressure projection");
}

void test_taylor_green_backend_parser_and_metal_rejection() {
  require(solver::parse_execution_backend("cpu") == solver::ExecutionBackend::cpu,
          "backend parser should accept cpu");
  require(solver::parse_execution_backend("metal") == solver::ExecutionBackend::metal,
          "backend parser should accept metal");

  bool rejected = false;
  try {
    solver::TaylorGreenConfig config = solver::load_taylor_green_config(
        source_path("benchmarks/taylor_green_smoke.cfg"));
    config.backend = solver::ExecutionBackend::metal;
    static_cast<void>(solver::run_taylor_green(config));
  } catch(const std::exception&) {
    rejected = true;
  }
  require(rejected, "metal backend should reject the 2D Taylor-Green path");
}

void test_taylor_green_cpu_vs_metal_small_3d() {
  if(!metal_backend_available()) {
    return;
  }

  solver::TaylorGreenConfig cpu_config = solver::load_taylor_green_config(
      source_path("benchmarks/taylor_green_3d_smoke.cfg"));
  cpu_config.nx = 32;
  cpu_config.ny = 32;
  cpu_config.nz = 32;
  cpu_config.final_time = 0.01;
  cpu_config.poisson_max_iterations = 160;
  cpu_config.backend = solver::ExecutionBackend::cpu;

  solver::TaylorGreenConfig metal_config = cpu_config;
  metal_config.backend = solver::ExecutionBackend::metal;

  solver::TaylorGreenState cpu_state = solver::initialize_taylor_green_state(cpu_config);
  const solver::TaylorGreenResult cpu_result = solver::run_taylor_green(cpu_config, &cpu_state);

  solver::TaylorGreenState metal_state = solver::initialize_taylor_green_state(metal_config);
  const solver::TaylorGreenResult metal_result =
      solver::run_taylor_green(metal_config, &metal_state);

  solver::PressureField cpu_pressure = cpu_state.pressure_total;
  solver::PressureField metal_pressure = metal_state.pressure_total;
  subtract_active_mean_inplace(cpu_pressure);
  subtract_active_mean_inplace(metal_pressure);

  require(metal_result.backend_used == solver::ExecutionBackend::metal,
          "metal run should report the metal backend");
  require(!metal_result.accelerator_name.empty(),
          "metal run should report the accelerator name");
  require(metal_result.final_step.divergence_l2 <= 1.0e-10,
          "metal Taylor-Green run should remain divergence controlled");
  require(metal_result.validation.pass, "metal Taylor-Green run should pass the analytic gate");

  const double velocity_x_error =
      relative_active_l2_difference(metal_state.velocity.x, cpu_state.velocity.x);
  const double velocity_y_error =
      relative_active_l2_difference(metal_state.velocity.y, cpu_state.velocity.y);
  const double velocity_z_error =
      active_l2_difference(metal_state.velocity.z, cpu_state.velocity.z);
  const double pressure_error = relative_active_l2_difference(metal_pressure, cpu_pressure);
  const double energy_error =
      std::abs(metal_result.final_kinetic_energy - cpu_result.final_kinetic_energy) /
      cpu_result.final_kinetic_energy;

  require(std::max({velocity_x_error, velocity_y_error, velocity_z_error}) <= 5.0e-4,
          "metal velocity field drifted too far from the CPU reference");
  require(pressure_error <= 5.0e-2,
          "metal pressure field drifted too far from the CPU reference: pressure_error=" +
              std::to_string(pressure_error));
  require(energy_error <= 1.0e-4,
          "metal kinetic energy drifted too far from the CPU reference");
}

void test_taylor_green_metal_vtk_export() {
  if(!metal_backend_available()) {
    return;
  }

  solver::TaylorGreenConfig config = solver::load_taylor_green_config(
      source_path("benchmarks/taylor_green_3d_smoke.cfg"));
  config.nx = 16;
  config.ny = 16;
  config.nz = 16;
  config.final_time = 0.005;
  config.validate_energy = false;
  config.backend = solver::ExecutionBackend::metal;

  solver::TaylorGreenState state = solver::initialize_taylor_green_state(config);
  static_cast<void>(solver::run_taylor_green(config, &state));

  const std::filesystem::path vtk_path = temp_path("solver_taylor_green_metal.vtk");
  solver::io::write_mac_fields_vtk(vtk_path.string(), state.velocity, state.pressure_total);

  std::ifstream input(vtk_path);
  const std::string contents((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  require(contents.find("# vtk DataFile Version 3.0") != std::string::npos,
          "metal VTK export should write the legacy VTK header");
  require(contents.find("VECTORS velocity double") != std::string::npos,
          "metal VTK export should contain the velocity vector field");
  require(contents.find("SCALARS pressure double 1") != std::string::npos,
          "metal VTK export should contain the pressure scalar field");

  std::error_code ignore_error;
  std::filesystem::remove(vtk_path, ignore_error);
}

void test_taylor_green_metal_cleanup_metadata() {
  if(!metal_backend_available()) {
    return;
  }

  solver::TaylorGreenConfig config = solver::load_taylor_green_config(
      source_path("benchmarks/taylor_green_3d_smoke.cfg"));
  config.nx = 16;
  config.ny = 16;
  config.nz = 16;
  config.final_time = 0.005;
  config.validate_energy = false;
  config.backend = solver::ExecutionBackend::metal;

  const solver::metal::TaylorGreenMetalRun raw_run = solver::metal::run_taylor_green(config);
  const solver::BoundaryConditionSet boundary_conditions =
      solver::make_taylor_green_boundary_conditions(config);
  const solver::ProjectionOptions projection_options{
      .dt = raw_run.state.metrics.dt,
      .density = 1.0,
      .poisson_max_iterations = config.poisson_max_iterations,
      .poisson_tolerance = config.poisson_tolerance,
  };
  solver::PressureField pressure_correction{raw_run.state.grid};
  solver::VelocityField diffusion{raw_run.state.grid};
  solver::VelocityField corrected{raw_run.state.grid};
  solver::ScalarField pressure_rhs{raw_run.state.grid};
  const solver::ProjectionDiagnostics cleanup = solver::project_velocity(
      raw_run.state.velocity,
      boundary_conditions,
      projection_options,
      pressure_correction,
      corrected,
      &pressure_rhs);

  solver::TaylorGreenState expected_state = raw_run.state;
  expected_state.velocity = corrected;
  axpy_active(expected_state.pressure_total, pressure_correction, 1.0);
  axpy_active(expected_state.pressure_total,
              pressure_rhs,
              -0.5 * config.viscosity * raw_run.state.metrics.dt);
  solver::compute_diffusion_term(corrected, config.viscosity, diffusion);
  solver::apply_total_pressure_boundary_conditions(
      boundary_conditions, diffusion, expected_state.pressure_total);

  solver::TaylorGreenState final_state = solver::initialize_taylor_green_state(config);
  const solver::TaylorGreenResult result = solver::run_taylor_green(config, &final_state);
  require(result.backend_elapsed_seconds > 0.0,
          "metal result should expose backend elapsed time");
  require(result.cleanup_elapsed_seconds > 0.0,
          "metal result should expose cleanup elapsed time");
  require(result.final_step.pressure_iterations == cleanup.pressure_solve.iterations,
          "metal final-step iterations should report the cleanup projection solve");
  require(std::abs(result.final_step.pressure_relative_residual -
                   cleanup.pressure_solve.relative_residual) <= 1.0e-12,
          "metal final-step residual should report the cleanup projection solve");
  require(relative_active_l2_difference(final_state.pressure_total, expected_state.pressure_total) <= 1.0e-12,
          "metal cleanup should publish the CPU-equivalent total-pressure update");
}

void test_taylor_green_metal_rejects_nonconverged_cleanup_projection() {
  if(!metal_backend_available()) {
    return;
  }

  solver::TaylorGreenConfig config = solver::load_taylor_green_config(
      source_path("benchmarks/taylor_green_3d_smoke.cfg"));
  config.nx = 16;
  config.ny = 16;
  config.nz = 16;
  config.final_time = 0.005;
  config.validate_energy = false;
  config.backend = solver::ExecutionBackend::metal;
  config.poisson_max_iterations = 20;
  config.poisson_tolerance = 1.0e-12;

  bool rejected = false;
  try {
    static_cast<void>(solver::run_taylor_green(config));
  } catch(const std::exception& exception) {
    rejected = std::string(exception.what()).find(
                   "Taylor-Green metal cleanup pressure solve did not converge") !=
               std::string::npos;
  }
  require(rejected, "metal Taylor-Green should reject a non-converged cleanup projection");
}

void test_taylor_green_metal_fails_fast_on_nonconverged_pressure_solve() {
  if(!metal_backend_available()) {
    return;
  }

  solver::TaylorGreenConfig config = solver::load_taylor_green_config(
      source_path("benchmarks/taylor_green_3d_smoke.cfg"));
  config.nx = 16;
  config.ny = 16;
  config.nz = 16;
  config.final_time = 0.005;
  config.validate_energy = false;
  config.backend = solver::ExecutionBackend::metal;
  config.poisson_max_iterations = 1;
  config.poisson_tolerance = 1.0e-30;

  bool rejected = false;
  try {
    static_cast<void>(solver::run_taylor_green(config));
  } catch(const std::exception& exception) {
    rejected = std::string(exception.what()).find("did not converge") != std::string::npos;
  }
  require(rejected, "metal Taylor-Green should fail fast on a non-converged internal solve");
}

void test_lid_driven_cavity_checkpoint_roundtrip_and_checksum() {
  solver::LidDrivenCavityConfig config = solver::load_lid_driven_cavity_config(
      source_path("benchmarks/lid_driven_cavity_smoke.cfg"));
  config.max_steps = 32;
  config.min_steps = 32;
  solver::LidDrivenCavityState state = solver::initialize_lid_driven_cavity_state(config);
  solver::run_lid_driven_cavity_steps(config, 6, state);

  const std::filesystem::path checkpoint_path = temp_path("solver_lid_roundtrip.chk");
  const std::filesystem::path corrupted_path = temp_path("solver_lid_roundtrip_corrupt.chk");
  solver::io::write_lid_driven_cavity_checkpoint(checkpoint_path.string(), config, state);
  const solver::io::LidDrivenCavityCheckpoint loaded =
      solver::io::load_lid_driven_cavity_checkpoint(checkpoint_path.string(), config);

  require(loaded.metadata.format_version == 1u, "checkpoint format version mismatch");
  require(loaded.metadata.endianness == "little", "checkpoint should record little-endian encoding");
  require(loaded.metadata.scalar_bytes == sizeof(double),
          "checkpoint precision metadata mismatch");
  require(loaded.metadata.checksum != 0u, "checkpoint checksum should be populated");
  require(loaded.metadata.build_hash != 0u, "checkpoint build hash should be populated");
  require(loaded.metadata.configuration_hash != 0u,
          "checkpoint configuration hash should be populated");
  require(loaded.state.has_previous_advection == state.has_previous_advection,
          "checkpoint should preserve AB2 history availability");
  require(metrics_bitwise_equal(loaded.state.metrics, state.metrics),
          "checkpoint should preserve the full step metrics bitwise");
  require(storage_bitwise_equal(loaded.state.velocity.x, state.velocity.x),
          "checkpoint should preserve velocity.x bitwise");
  require(storage_bitwise_equal(loaded.state.velocity.y, state.velocity.y),
          "checkpoint should preserve velocity.y bitwise");
  require(storage_bitwise_equal(loaded.state.velocity.z, state.velocity.z),
          "checkpoint should preserve velocity.z bitwise");
  require(storage_bitwise_equal(loaded.state.advection_previous.x, state.advection_previous.x),
          "checkpoint should preserve advection_previous.x bitwise");
  require(storage_bitwise_equal(loaded.state.advection_previous.y, state.advection_previous.y),
          "checkpoint should preserve advection_previous.y bitwise");
  require(storage_bitwise_equal(loaded.state.advection_previous.z, state.advection_previous.z),
          "checkpoint should preserve advection_previous.z bitwise");
  require(storage_bitwise_equal(loaded.state.pressure_total, state.pressure_total),
          "checkpoint should preserve pressure_total bitwise");

  std::ifstream input(checkpoint_path, std::ios::binary);
  std::vector<char> bytes((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  require(bytes.size() > 32, "checkpoint corruption test needs a nontrivial payload");
  bytes.back() = static_cast<char>(bytes.back() ^ 0x5a);
  std::ofstream output(corrupted_path, std::ios::binary);
  output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
  output.close();

  bool checksum_failed = false;
  try {
    static_cast<void>(
        solver::io::load_lid_driven_cavity_checkpoint(corrupted_path.string(), config));
  } catch(const std::exception&) {
    checksum_failed = true;
  }
  require(checksum_failed, "corrupted checkpoint should fail checksum verification");

  std::error_code ignore_error;
  std::filesystem::remove(checkpoint_path, ignore_error);
  std::filesystem::remove(corrupted_path, ignore_error);
}

void test_lid_driven_cavity_rejects_nonconverged_pressure_projection() {
  solver::LidDrivenCavityConfig config = solver::load_lid_driven_cavity_config(
      source_path("benchmarks/lid_driven_cavity_smoke.cfg"));
  config.poisson_max_iterations = 1;
  config.poisson_tolerance = 1.0e-30;

  bool rejected = false;
  try {
    static_cast<void>(solver::run_lid_driven_cavity(config));
  } catch(const std::exception& exception) {
    rejected = std::string(exception.what()).find(
                   "lid-driven cavity pressure solve did not converge") != std::string::npos;
  }
  require(rejected, "lid-driven cavity should reject a non-converged pressure projection");
}

void test_lid_driven_cavity_restart_is_bitwise_deterministic() {
  solver::LidDrivenCavityConfig config = solver::load_lid_driven_cavity_config(
      source_path("benchmarks/lid_driven_cavity_smoke.cfg"));
  config.nx = 24;
  config.ny = 24;
  config.max_steps = 200;
  config.min_steps = 200;
  config.steady_tolerance = 1.0e-30;
  config.validate_reference = false;

  solver::LidDrivenCavityState uninterrupted = solver::initialize_lid_driven_cavity_state(config);
  solver::run_lid_driven_cavity_steps(config, 200, uninterrupted);

  solver::LidDrivenCavityState split = solver::initialize_lid_driven_cavity_state(config);
  solver::run_lid_driven_cavity_steps(config, 100, split);

  const std::filesystem::path checkpoint_path = temp_path("solver_lid_restart.chk");
  solver::io::write_lid_driven_cavity_checkpoint(checkpoint_path.string(), config, split);
  split = solver::io::load_lid_driven_cavity_checkpoint(checkpoint_path.string(), config).state;
  solver::run_lid_driven_cavity_steps(config, 100, split);

  require(split.has_previous_advection == uninterrupted.has_previous_advection,
          "restart path drifted in AB2 history availability");
  require(metrics_bitwise_equal(split.metrics, uninterrupted.metrics),
          "restart path drifted in final step metrics");
  require(storage_bitwise_equal(split.velocity.x, uninterrupted.velocity.x),
          "restart path drifted in velocity.x");
  require(storage_bitwise_equal(split.velocity.y, uninterrupted.velocity.y),
          "restart path drifted in velocity.y");
  require(storage_bitwise_equal(split.velocity.z, uninterrupted.velocity.z),
          "restart path drifted in velocity.z");
  require(storage_bitwise_equal(split.advection_previous.x, uninterrupted.advection_previous.x),
          "restart path drifted in advection_previous.x");
  require(storage_bitwise_equal(split.advection_previous.y, uninterrupted.advection_previous.y),
          "restart path drifted in advection_previous.y");
  require(storage_bitwise_equal(split.advection_previous.z, uninterrupted.advection_previous.z),
          "restart path drifted in advection_previous.z");
  require(storage_bitwise_equal(split.pressure_total, uninterrupted.pressure_total),
          "restart path drifted in pressure_total");

  const solver::LidDrivenCavityResult split_result =
      solver::finalize_lid_driven_cavity_result(config, split);
  const solver::LidDrivenCavityResult uninterrupted_result =
      solver::finalize_lid_driven_cavity_result(config, uninterrupted);
  require(metrics_bitwise_equal(split_result.final_step, uninterrupted_result.final_step),
          "restart path drifted in finalized metrics");
  for(std::size_t index = 0; index < split_result.u_vertical_centerline.value.size(); ++index) {
    require(std::bit_cast<std::uint64_t>(split_result.u_vertical_centerline.value[index]) ==
                std::bit_cast<std::uint64_t>(uninterrupted_result.u_vertical_centerline.value[index]),
            "restart path drifted in finalized u centerline");
  }
  for(std::size_t index = 0; index < split_result.v_horizontal_centerline.value.size(); ++index) {
    require(std::bit_cast<std::uint64_t>(split_result.v_horizontal_centerline.value[index]) ==
                std::bit_cast<std::uint64_t>(uninterrupted_result.v_horizontal_centerline.value[index]),
            "restart path drifted in finalized v centerline");
  }

  std::error_code ignore_error;
  std::filesystem::remove(checkpoint_path, ignore_error);
}

void test_lid_driven_cavity_checkpoint_rejects_unsupported_version() {
  solver::LidDrivenCavityConfig config = solver::load_lid_driven_cavity_config(
      source_path("benchmarks/lid_driven_cavity_smoke.cfg"));
  solver::LidDrivenCavityState state = solver::initialize_lid_driven_cavity_state(config);
  solver::run_lid_driven_cavity_steps(config, 2, state);

  const std::filesystem::path checkpoint_path = temp_path("solver_lid_version.chk");
  solver::io::write_lid_driven_cavity_checkpoint(checkpoint_path.string(), config, state);
  std::vector<std::uint8_t> bytes = read_binary_file(checkpoint_path);
  write_le_u32(bytes, 8, 2u);
  write_binary_file(checkpoint_path, bytes);

  require_exception_contains(
      [&] { static_cast<void>(solver::io::load_lid_driven_cavity_checkpoint(checkpoint_path.string(), config)); },
      "unsupported checkpoint format version",
      "checkpoint loader should reject unsupported format versions");

  std::error_code ignore_error;
  std::filesystem::remove(checkpoint_path, ignore_error);
}

void test_lid_driven_cavity_checkpoint_rejects_build_hash_mismatch() {
  solver::LidDrivenCavityConfig config = solver::load_lid_driven_cavity_config(
      source_path("benchmarks/lid_driven_cavity_smoke.cfg"));
  solver::LidDrivenCavityState state = solver::initialize_lid_driven_cavity_state(config);
  solver::run_lid_driven_cavity_steps(config, 2, state);

  const std::filesystem::path checkpoint_path = temp_path("solver_lid_build_hash.chk");
  solver::io::write_lid_driven_cavity_checkpoint(checkpoint_path.string(), config, state);
  std::vector<std::uint8_t> bytes = read_binary_file(checkpoint_path);
  const CheckpointSectionView metadata = find_checkpoint_section(bytes, {'M', 'E', 'T', 'A'});
  const std::size_t build_hash_offset = metadata.payload_offset + 4;
  write_le_u64(bytes, build_hash_offset, read_le_u64(bytes, build_hash_offset) ^ 0x1ull);
  refresh_checkpoint_payload_header(bytes);
  write_binary_file(checkpoint_path, bytes);

  require_exception_contains(
      [&] { static_cast<void>(solver::io::load_lid_driven_cavity_checkpoint(checkpoint_path.string(), config)); },
      "checkpoint build hash does not match the current executable",
      "checkpoint loader should reject build-hash mismatches");

  std::error_code ignore_error;
  std::filesystem::remove(checkpoint_path, ignore_error);
}

void test_lid_driven_cavity_checkpoint_rejects_configuration_mismatch() {
  solver::LidDrivenCavityConfig written_config = solver::load_lid_driven_cavity_config(
      source_path("benchmarks/lid_driven_cavity_smoke.cfg"));
  solver::LidDrivenCavityState state = solver::initialize_lid_driven_cavity_state(written_config);
  solver::run_lid_driven_cavity_steps(written_config, 2, state);

  const std::filesystem::path checkpoint_path = temp_path("solver_lid_config_hash.chk");
  solver::io::write_lid_driven_cavity_checkpoint(checkpoint_path.string(), written_config, state);

  solver::LidDrivenCavityConfig expected_config = written_config;
  expected_config.poisson_max_iterations += 1;
  require_exception_contains(
      [&] { static_cast<void>(solver::io::load_lid_driven_cavity_checkpoint(checkpoint_path.string(), expected_config)); },
      "checkpoint configuration hash does not match the requested config",
      "checkpoint loader should reject configuration mismatches");

  std::error_code ignore_error;
  std::filesystem::remove(checkpoint_path, ignore_error);
}

void test_checkpoint_test_helpers_reject_invalid_io() {
  const std::filesystem::path missing_path = temp_path("solver_missing_checkpoint_helper.chk");
  require_exception_contains(
      [&] { static_cast<void>(read_binary_file(missing_path)); },
      "read_binary_file could not open path",
      "read_binary_file should reject missing paths");

  const std::filesystem::path unwritable_path =
      std::filesystem::path("/definitely/not/a/real/directory") / "solver_write_fail.chk";
  require_exception_contains(
      [&] { write_binary_file(unwritable_path, {0x01u, 0x02u}); },
      "write_binary_file could not open path",
      "write_binary_file should reject unwritable paths");
}

void test_checkpoint_test_helpers_reject_undersized_header() {
  std::vector<std::uint8_t> bytes(8u, 0u);
  require_exception_contains(
      [&] { refresh_checkpoint_payload_header(bytes); },
      "refresh_checkpoint_payload_header requires a full checkpoint header",
      "refresh_checkpoint_payload_header should reject undersized buffers");
}

void test_lid_driven_cavity_vtk_export() {
  solver::LidDrivenCavityConfig config = solver::load_lid_driven_cavity_config(
      source_path("benchmarks/lid_driven_cavity_smoke.cfg"));
  solver::LidDrivenCavityState state = solver::initialize_lid_driven_cavity_state(config);
  solver::run_lid_driven_cavity_steps(config, 2, state);

  const std::filesystem::path vtk_path = temp_path("solver_lid_state.vtk");
  solver::io::write_lid_driven_cavity_vtk(vtk_path.string(), state);

  std::ifstream input(vtk_path);
  const std::string contents((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  require(contents.find("# vtk DataFile Version 3.0") != std::string::npos,
          "VTK export should write the legacy VTK header");
  require(contents.find("DATASET STRUCTURED_POINTS") != std::string::npos,
          "VTK export should declare a structured-points dataset");
  require(contents.find("VECTORS velocity double") != std::string::npos,
          "VTK export should contain the velocity vector field");
  require(contents.find("SCALARS pressure double 1") != std::string::npos,
          "VTK export should contain the pressure scalar field");

  std::error_code ignore_error;
  std::filesystem::remove(vtk_path, ignore_error);
}

void test_lid_driven_cavity_smoke_run() {
  const solver::LidDrivenCavityConfig config = solver::load_lid_driven_cavity_config(
      source_path("benchmarks/lid_driven_cavity_smoke.cfg"));
  const solver::LidDrivenCavityResult result = solver::run_lid_driven_cavity(config);

  require(result.final_step.step == config.max_steps,
          "smoke cavity run should consume the fixed short-step budget");
  require(result.final_step.dt > 0.0, "smoke cavity timestep should be positive");
  require(result.final_step.max_cfl <= config.cfl_limit + 1.0e-12,
          "smoke cavity CFL should respect the configured ceiling");
  require(result.final_step.divergence_l2 <= 1.0e-10,
          "smoke cavity run should remain divergence controlled");
  require(result.final_step.pressure_iterations >= 0,
          "smoke cavity run should report pressure iterations");
  require(result.validation.reference_dataset.empty(),
          "smoke cavity run should skip benchmark validation");
  require(result.extrema.u_vertical_max > 0.0, "smoke cavity lid should induce positive u motion");
  require(result.extrema.u_vertical_min < 0.0,
          "smoke cavity recirculation should induce negative u motion");
}

void test_lid_driven_cavity_reference_validation_gate() {
  const solver::LidDrivenCavityReference reference = solver::re100_centerline_reference_envelope();
  const solver::LidDrivenCavityReferencePoint& u_top = reference.points[0];
  const solver::LidDrivenCavityReferencePoint& u_mid = reference.points[1];
  const solver::LidDrivenCavityReferencePoint& v_left = reference.points[2];
  const solver::LidDrivenCavityReferencePoint& v_right = reference.points[3];

  const double u_slope = (u_top.value - u_mid.value) / (u_top.coordinate - u_mid.coordinate);
  const double u_intercept = u_top.value - u_slope * u_top.coordinate;
  const double v_slope = (v_right.value - v_left.value) / (v_right.coordinate - v_left.coordinate);
  const double v_intercept = v_left.value - v_slope * v_left.coordinate;

  solver::LidDrivenCavityResult passing{};
  passing.u_vertical_centerline.coordinate = {0.20, 0.40, 0.60, 0.80, 0.99};
  for(const double coordinate : passing.u_vertical_centerline.coordinate) {
    passing.u_vertical_centerline.value.push_back(u_slope * coordinate + u_intercept);
  }
  passing.v_horizontal_centerline.coordinate = {0.10, 0.30, 0.55, 0.72, 0.92};
  for(const double coordinate : passing.v_horizontal_centerline.coordinate) {
    passing.v_horizontal_centerline.value.push_back(v_slope * coordinate + v_intercept);
  }
  passing.final_step.divergence_l2 = 5.0e-11;

  const solver::LidDrivenCavityValidation validation =
      solver::validate_lid_driven_cavity_re100(passing);
  require(validation.pass, "reference-matching cavity result should pass validation");
  require(validation.reference_dataset == reference.dataset, "wrong reference dataset label");
  require(validation.max_relative_error <= 1.0e-12,
          "validation should reproduce the named reference sample points");
  require(validation.points[0].label == u_top.label, "validation point labels drifted");
  require(std::abs(validation.points[3].sample_value - v_right.value) <= 1.0e-12,
          "validation should sample the horizontal-right reference point");

  solver::LidDrivenCavityResult failing = passing;
  for(double& value : failing.v_horizontal_centerline.value) {
    value *= 0.97;
  }
  const solver::LidDrivenCavityValidation failed_validation =
      solver::validate_lid_driven_cavity_re100(failing);
  require(!failed_validation.pass, "3 percent sample drift should fail the 2 percent gate");
  require(failed_validation.max_relative_error >= 0.03 - 1.0e-12,
          "validation should report the largest sample-point drift");
}

}  // namespace

int main() {
  try {
  test_build_profile_is_locked();
  test_runtime_platform_is_supported();
  test_project_name_matches_orchard_flow();
  test_banner_contains_profile();
    test_grid_coordinates();
    test_pressure_layout_indexing();
    test_ghost_cell_access_and_boundary_ranges();
    test_memory_layout_and_alignment();
    test_cell_and_face_placement();
    test_double_precision_storage_contract();
    test_advection_options_and_cfl_diagnostic();
    test_projection_boundary_mapping();
    test_predictor_adi_preserves_quiescent_state();
    test_predictor_adi_preserves_periodic_couette_state();
    test_total_pressure_boundary_conditions_general_bc();
    test_predictor_adi_matches_factorized_dense_reference();
    test_lid_driven_cavity_predictor_rhs_matches_manual_formula();
    test_poisson_mgpcg_discrete_dirichlet_recovery();
    test_poisson_mgpcg_pure_neumann_zero_mean_recovery();
    test_manufactured_solution_convergence();
    test_diffusion_term_matches_scaled_laplacian();
    test_taylor_green_step_behavior();
    test_bounded_advection_regression_case();
    test_static_projection_preserves_quiescent_fluid();
    test_pure_neumann_projection_recovers_zero_mean_pressure();
    test_lid_driven_cavity_config_loader_and_bc_subset();
    test_channel_flow_config_loader_and_boundary_conditions();
    test_channel_flow_couette_profile_validation();
    test_channel_flow_poiseuille_profile_validation();
    test_channel_flow_rejects_nonconverged_pressure_projection();
    test_taylor_green_config_loader_and_boundary_conditions();
    test_taylor_green_3d_config_loader_and_boundary_conditions();
    test_taylor_green_smoke_validation();
    test_taylor_green_3d_smoke_validation();
    test_taylor_green_rejects_nonconverged_pressure_projection();
    test_taylor_green_backend_parser_and_metal_rejection();
    test_taylor_green_cpu_vs_metal_small_3d();
    test_taylor_green_metal_vtk_export();
    test_taylor_green_metal_cleanup_metadata();
    test_taylor_green_metal_rejects_nonconverged_cleanup_projection();
    test_taylor_green_metal_fails_fast_on_nonconverged_pressure_solve();
    test_lid_driven_cavity_checkpoint_roundtrip_and_checksum();
    test_lid_driven_cavity_checkpoint_rejects_unsupported_version();
    test_lid_driven_cavity_checkpoint_rejects_build_hash_mismatch();
    test_lid_driven_cavity_checkpoint_rejects_configuration_mismatch();
    test_checkpoint_test_helpers_reject_invalid_io();
    test_checkpoint_test_helpers_reject_undersized_header();
    test_lid_driven_cavity_restart_is_bitwise_deterministic();
    test_lid_driven_cavity_vtk_export();
    test_lid_driven_cavity_rejects_nonconverged_pressure_projection();
    test_lid_driven_cavity_smoke_run();
    test_lid_driven_cavity_reference_validation_gate();
  } catch(const std::exception& exception) {
    std::cerr << "solver_tests failed: " << exception.what() << '\n';
    return 1;
  }

  std::cout << "solver_tests passed" << '\n';
  return 0;
}
