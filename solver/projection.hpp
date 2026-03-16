#pragma once

#include "core/fields.hpp"

#include <array>
#include <string>
#include <vector>

namespace solver {

enum class PhysicalBoundaryType : int {
  no_slip_wall = 0,
  prescribed_velocity = 1,
  symmetry = 2,
  fixed_pressure = 3,
  periodic = 4,
};

enum class PressureBoundaryType : int {
  neumann = 0,
  dirichlet = 1,
  periodic = 2,
};

struct BoundaryCondition {
  PhysicalBoundaryType type = PhysicalBoundaryType::no_slip_wall;
  std::array<double, 3> velocity{0.0, 0.0, 0.0};
  double pressure = 0.0;
  double pressure_gradient = 0.0;
};

struct BoundaryConditionSet {
  std::array<BoundaryCondition, 6> faces{};

  [[nodiscard]] static BoundaryConditionSet all(PhysicalBoundaryType type);
  [[nodiscard]] static BoundaryConditionSet cavity();

  [[nodiscard]] BoundaryCondition& operator[](BoundaryFace face) noexcept {
    return faces[static_cast<std::size_t>(face)];
  }

  [[nodiscard]] const BoundaryCondition& operator[](BoundaryFace face) const noexcept {
    return faces[static_cast<std::size_t>(face)];
  }
};

struct PressureBoundaryCondition {
  PressureBoundaryType type = PressureBoundaryType::neumann;
  double value = 0.0;
  double gradient = 0.0;
};

struct PressureBoundarySet {
  std::array<PressureBoundaryCondition, 6> faces{};

  [[nodiscard]] PressureBoundaryCondition& operator[](BoundaryFace face) noexcept {
    return faces[static_cast<std::size_t>(face)];
  }

  [[nodiscard]] const PressureBoundaryCondition& operator[](BoundaryFace face) const noexcept {
    return faces[static_cast<std::size_t>(face)];
  }
};

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

std::string to_string(PhysicalBoundaryType type);
std::string to_string(PressureBoundaryType type);

// Projection solves for the pressure-correction variable phi.
PressureBoundarySet derive_pressure_correction_boundary_conditions(
    const BoundaryConditionSet& boundary_conditions);

void apply_velocity_boundary_conditions(const BoundaryConditionSet& boundary_conditions,
                                        VelocityField& velocity);

void apply_pressure_boundary_conditions(const PressureBoundarySet& boundary_conditions,
                                        PressureField& pressure);

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
