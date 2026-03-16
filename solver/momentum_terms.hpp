#pragma once

#include "core/fields.hpp"

#include <string>

namespace solver {

enum class AdvectionScheme : int {
  tvd = 0,
  upwind = 1,
  central = 2,
};

enum class FluxLimiter : int {
  van_leer = 0,
};

struct AdvectionOptions {
  AdvectionScheme scheme = AdvectionScheme::tvd;
  FluxLimiter limiter = FluxLimiter::van_leer;
};

struct CflDiagnostics {
  double dt = 0.0;
  double max_cfl = 0.0;
  double max_u = 0.0;
  double max_v = 0.0;
  double max_w = 0.0;
};

std::string to_string(AdvectionScheme scheme);
std::string to_string(FluxLimiter limiter);
std::string describe(const AdvectionOptions& options);

void compute_advection_term(const VelocityField& velocity,
                            const AdvectionOptions& options,
                            VelocityField& advection);

void compute_diffusion_term(const VelocityField& velocity,
                            double viscosity,
                            VelocityField& diffusion);

CflDiagnostics compute_advective_cfl(const VelocityField& velocity, double dt);

}  // namespace solver
