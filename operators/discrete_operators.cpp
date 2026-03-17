#include "operators/discrete_operators.hpp"

#include <stdexcept>
#include <string>

namespace solver::operators {

namespace {

void require_same_layout(const StructuredField& left,
                         const StructuredField& right,
                         const char* operation_name) {
  if(!left.layout().same_shape_as(right.layout())) {
    throw std::invalid_argument(std::string(operation_name) + ": incompatible field layouts");
  }
}

void require_cell_centered_pressure(const PressureField& pressure) {
  if(pressure.layout().location() != FieldLocation::cell_center) {
    throw std::invalid_argument("compute_gradient expects cell-centered pressure");
  }
}

double second_derivative(const StructuredField& input,
                         const int i,
                         const int j,
                         const int k,
                         const Axis axis) {
  const double spacing = input.layout().grid().spacing(axis);
  const double inverse_spacing_squared = 1.0 / (spacing * spacing);

  switch(axis) {
    case Axis::x:
      return (input(i + 1, j, k) - 2.0 * input(i, j, k) + input(i - 1, j, k)) *
             inverse_spacing_squared;
    case Axis::y:
      return (input(i, j + 1, k) - 2.0 * input(i, j, k) + input(i, j - 1, k)) *
             inverse_spacing_squared;
    case Axis::z:
      return (input(i, j, k + 1) - 2.0 * input(i, j, k) + input(i, j, k - 1)) *
             inverse_spacing_squared;
  }

  __builtin_unreachable();
}

void compute_component_laplacian(const StructuredField& input, StructuredField& output) {
  require_same_layout(input, output, "compute_laplacian");

  const IndexRange3D active = input.layout().active_range();
  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        output(i, j, k) = second_derivative(input, i, j, k, Axis::x) +
                          second_derivative(input, i, j, k, Axis::y) +
                          second_derivative(input, i, j, k, Axis::z);
      }
    }
  }
}

}  // namespace

void compute_gradient(const PressureField& pressure, VelocityField& gradient) {
  require_cell_centered_pressure(pressure);

  const IndexRange3D active_x = gradient.x.layout().active_range();
  for(int k = active_x.k_begin; k < active_x.k_end; ++k) {
    for(int j = active_x.j_begin; j < active_x.j_end; ++j) {
      for(int i = active_x.i_begin; i < active_x.i_end; ++i) {
        gradient.x(i, j, k) =
            (pressure(i, j, k) - pressure(i - 1, j, k)) / pressure.layout().grid().dx;
      }
    }
  }

  const IndexRange3D active_y = gradient.y.layout().active_range();
  for(int k = active_y.k_begin; k < active_y.k_end; ++k) {
    for(int j = active_y.j_begin; j < active_y.j_end; ++j) {
      for(int i = active_y.i_begin; i < active_y.i_end; ++i) {
        gradient.y(i, j, k) =
            (pressure(i, j, k) - pressure(i, j - 1, k)) / pressure.layout().grid().dy;
      }
    }
  }

  const IndexRange3D active_z = gradient.z.layout().active_range();
  for(int k = active_z.k_begin; k < active_z.k_end; ++k) {
    for(int j = active_z.j_begin; j < active_z.j_end; ++j) {
      for(int i = active_z.i_begin; i < active_z.i_end; ++i) {
        gradient.z(i, j, k) =
            (pressure(i, j, k) - pressure(i, j, k - 1)) / pressure.layout().grid().dz;
      }
    }
  }
}

void compute_divergence(const VelocityField& velocity, ScalarField& divergence) {
  const FieldLayout& layout = divergence.layout();
  const IndexRange3D active = layout.active_range();
  const Grid& grid = layout.grid();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        const double du_dx = (velocity.x(i + 1, j, k) - velocity.x(i, j, k)) / grid.dx;
        const double dv_dy = (velocity.y(i, j + 1, k) - velocity.y(i, j, k)) / grid.dy;
        const double dw_dz = (velocity.z(i, j, k + 1) - velocity.z(i, j, k)) / grid.dz;
        divergence(i, j, k) = du_dx + dv_dy + dw_dz;
      }
    }
  }
}

void compute_laplacian(const StructuredField& input, StructuredField& output) {
  compute_component_laplacian(input, output);
}

void compute_laplacian(const VelocityField& input, VelocityField& output) {
  compute_component_laplacian(input.x, output.x);
  compute_component_laplacian(input.y, output.y);
  compute_component_laplacian(input.z, output.z);
}

}  // namespace solver::operators
