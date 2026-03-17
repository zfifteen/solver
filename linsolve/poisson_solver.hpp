#pragma once

#include "solver/projection.hpp"

#include <string>

namespace solver::linsolve {

enum class MultigridCycle : int {
  v_cycle = 0,
};

enum class MultigridSmoother : int {
  damped_jacobi = 0,
};

struct MultigridPolicy {
  MultigridCycle cycle = MultigridCycle::v_cycle;
  MultigridSmoother smoother = MultigridSmoother::damped_jacobi;
  int pre_smoothing_steps = 2;
  int post_smoothing_steps = 2;
  double jacobi_omega = 2.0 / 3.0;
};

std::string to_string(MultigridCycle cycle);
std::string to_string(MultigridSmoother smoother);

[[nodiscard]] MultigridPolicy default_multigrid_policy();

void build_poisson_rhs_from_pressure(const PressureField& pressure,
                                     const PressureBoundarySet& boundary_conditions,
                                     ScalarField& rhs);

[[nodiscard]] PoissonSolveDiagnostics solve_pressure_poisson(const ScalarField& rhs,
                                               const PressureBoundarySet& boundary_conditions,
                                               const ProjectionOptions& options,
                                               PressureField& pressure);

}  // namespace solver::linsolve
