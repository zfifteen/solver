#pragma once

#include "bc/boundary_conditions.hpp"
#include "core/fields.hpp"

#include <string>
#include <vector>

namespace solver {

struct HelmholtzDiagnostics {
  int line_solves = 0;
};

struct ProjectionOptions {
  double dt = 0.0;
  double density = 1.0;
  int poisson_max_iterations = 4000;
  double poisson_tolerance = 1.0e-12;
};

struct PoissonSolveDiagnostics {
  int iterations = 0;
  double initial_residual_l2 = 0.0;
  double final_residual_l2 = 0.0;
  double relative_residual = 0.0;
  double removed_rhs_mean = 0.0;
  int multigrid_levels = 0;
  int coarse_unknowns = 0;
  int pre_smoothing_steps = 0;
  int post_smoothing_steps = 0;
  bool converged = false;
  bool zero_mean_enforced = false;
  std::string solver = "mgpcg";
  std::string preconditioner = "geometric_multigrid";
  std::string cycle = "v_cycle";
  std::string smoother = "damped_jacobi";
  std::vector<double> residual_history{};
};

struct ProjectionDiagnostics {
  double rhs_l2 = 0.0;
  double divergence_l2_before = 0.0;
  double divergence_l2_after = 0.0;
  double pressure_mean = 0.0;
  PoissonSolveDiagnostics pressure_solve{};
};

HelmholtzDiagnostics solve_predictor_adi(const VelocityField& rhs,
                                         double alpha,
                                         const BoundaryConditionSet& boundary_conditions,
                                         VelocityField& predicted_velocity);

void build_pressure_rhs(const VelocityField& predicted_velocity,
                        const BoundaryConditionSet& boundary_conditions,
                        const ProjectionOptions& options,
                        ScalarField& rhs);

void correct_velocity(const VelocityField& predicted_velocity,
                      const PressureField& pressure_correction,
                      const BoundaryConditionSet& boundary_conditions,
                      const ProjectionOptions& options,
                      VelocityField& corrected_velocity);

ProjectionDiagnostics project_velocity(const VelocityField& predicted_velocity,
                                       const BoundaryConditionSet& boundary_conditions,
                                       const ProjectionOptions& options,
                                       PressureField& pressure_correction,
                                       VelocityField& corrected_velocity,
                                       ScalarField* pressure_rhs = nullptr);

}  // namespace solver
