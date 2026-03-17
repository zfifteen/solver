#include "solver/operator_verification.hpp"

#include "core/fields.hpp"
#include "operators/discrete_operators.hpp"

#include <cmath>
#include <stdexcept>

namespace solver {

namespace {

double pi() {
  return std::acos(-1.0);
}

double square(const double value) {
  return value * value;
}

template <typename Field, typename ValueFn>
void fill_storage(Field& field, ValueFn&& value_fn) {
  const Extent3D storage = field.layout().storage_extent();

  for(int k = 0; k < storage.nz; ++k) {
    const double z = field.layout().coordinate_for_storage_index(Axis::z, k);
    for(int j = 0; j < storage.ny; ++j) {
      const double y = field.layout().coordinate_for_storage_index(Axis::y, j);
      for(int i = 0; i < storage.nx; ++i) {
        const double x = field.layout().coordinate_for_storage_index(Axis::x, i);
        field(i, j, k) = value_fn(x, y, z);
      }
    }
  }
}

template <typename Field, typename ValueFn>
double active_l2_error(const Field& field, ValueFn&& exact_value) {
  const IndexRange3D active = field.layout().active_range();

  double sum = 0.0;
  std::size_t count = 0;
  for(int k = active.k_begin; k < active.k_end; ++k) {
    const double z = field.layout().coordinate_for_storage_index(Axis::z, k);
    for(int j = active.j_begin; j < active.j_end; ++j) {
      const double y = field.layout().coordinate_for_storage_index(Axis::y, j);
      for(int i = active.i_begin; i < active.i_end; ++i) {
        const double x = field.layout().coordinate_for_storage_index(Axis::x, i);
        sum += square(field(i, j, k) - exact_value(x, y, z));
        ++count;
      }
    }
  }

  if(count == 0) {
    return 0.0;
  }

  return std::sqrt(sum / static_cast<double>(count));
}

}  // namespace

OperatorManufacturedSolutionResult run_operator_manufactured_solution_case(const int resolution) {
  if(resolution <= 0) {
    throw std::invalid_argument("run_operator_manufactured_solution_case expects resolution > 0");
  }

  const double domain_length = 2.0 * pi();
  const double spacing = domain_length / static_cast<double>(resolution);
  const Grid grid{resolution, resolution, 1, spacing, spacing, 1.0, 1};

  PressureField pressure{grid};
  VelocityField exact_velocity{grid};
  VelocityField gradient{grid};
  ScalarField divergence{grid};
  PressureField laplacian{grid};

  fill_storage(pressure, [](const double x, const double y, double) {
    return std::sin(x) * std::cos(y);
  });

  fill_storage(exact_velocity.x, [](const double x, const double y, double) {
    return std::sin(x) * std::cos(y);
  });
  fill_storage(exact_velocity.y, [](const double x, const double y, double) {
    return std::cos(x) * std::sin(y);
  });
  fill_storage(exact_velocity.z, [](double, double, double) {
    return 0.0;
  });

  operators::compute_gradient(pressure, gradient);
  operators::compute_divergence(exact_velocity, divergence);
  operators::compute_laplacian(pressure, laplacian);

  const double gradient_x_error = active_l2_error(gradient.x, [](const double x, const double y, double) {
    return std::cos(x) * std::cos(y);
  });
  const double gradient_y_error = active_l2_error(gradient.y, [](const double x, const double y, double) {
    return -std::sin(x) * std::sin(y);
  });

  return OperatorManufacturedSolutionResult{
      .resolution = resolution,
      .domain_length = domain_length,
      .gradient_error = std::sqrt(0.5 * (square(gradient_x_error) + square(gradient_y_error))),
      .divergence_error = active_l2_error(divergence, [](const double x, const double y, double) {
        return 2.0 * std::cos(x) * std::cos(y);
      }),
      .laplacian_error = active_l2_error(laplacian, [](const double x, const double y, double) {
        return -2.0 * std::sin(x) * std::cos(y);
      }),
  };
}

}  // namespace solver
