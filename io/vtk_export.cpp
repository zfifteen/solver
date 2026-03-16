#include "io/vtk_export.hpp"

#include <fstream>
#include <stdexcept>

namespace solver::io {

namespace {

bool same_grid(const Grid& left, const Grid& right) {
  return left.nx == right.nx && left.ny == right.ny && left.nz == right.nz &&
         left.dx == right.dx && left.dy == right.dy && left.dz == right.dz &&
         left.ghost_layers == right.ghost_layers;
}

}  // namespace

void write_mac_fields_vtk(const std::string& path,
                          const VelocityField& velocity,
                          const PressureField& pressure_total) {
  const Grid& pressure_grid = pressure_total.layout().grid();
  if(!same_grid(velocity.x.layout().grid(), pressure_grid) ||
     !same_grid(velocity.y.layout().grid(), pressure_grid) ||
     !same_grid(velocity.z.layout().grid(), pressure_grid)) {
    throw std::invalid_argument("VTK export expects velocity and pressure on the same grid");
  }

  std::ofstream output(path);
  if(!output.is_open()) {
    throw std::runtime_error("unable to open VTK file for writing: " + path);
  }

  const IndexRange3D active = pressure_total.layout().active_range();
  const std::size_t points = active.extent().cell_count();

  output << "# vtk DataFile Version 3.0\n";
  output << "solver checkpoint state\n";
  output << "ASCII\n";
  output << "DATASET STRUCTURED_POINTS\n";
  output << "DIMENSIONS " << pressure_grid.nx << ' ' << pressure_grid.ny << ' '
         << pressure_grid.nz << '\n';
  output << "ORIGIN " << 0.5 * pressure_grid.dx << ' ' << 0.5 * pressure_grid.dy << ' '
         << 0.5 * pressure_grid.dz << '\n';
  output << "SPACING " << pressure_grid.dx << ' ' << pressure_grid.dy << ' ' << pressure_grid.dz
         << '\n';
  output << "POINT_DATA " << points << '\n';
  output << "VECTORS velocity double\n";

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        const double u = 0.5 * (velocity.x(i, j, k) + velocity.x(i + 1, j, k));
        const double v = 0.5 * (velocity.y(i, j, k) + velocity.y(i, j + 1, k));
        const double w = 0.5 * (velocity.z(i, j, k) + velocity.z(i, j, k + 1));
        output << u << ' ' << v << ' ' << w << '\n';
      }
    }
  }

  output << "SCALARS pressure double 1\n";
  output << "LOOKUP_TABLE default\n";
  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        output << pressure_total(i, j, k) << '\n';
      }
    }
  }

  if(!output.good()) {
    throw std::runtime_error("failed while writing VTK file: " + path);
  }
}

void write_lid_driven_cavity_vtk(const std::string& path,
                                 const LidDrivenCavityState& state) {
  write_mac_fields_vtk(path, state.velocity, state.pressure_total);
}

}  // namespace solver::io
