#pragma once

#include "solver/lid_driven_cavity.hpp"

#include <string>

namespace solver::io {

void write_mac_fields_vtk(const std::string& path,
                          const VelocityField& velocity,
                          const PressureField& pressure_total);

void write_lid_driven_cavity_vtk(const std::string& path,
                                 const LidDrivenCavityState& state);

}  // namespace solver::io
