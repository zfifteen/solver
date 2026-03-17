#pragma once

#include "solver/momentum_terms.hpp"
#include "solver/projection.hpp"

#include <string>

namespace solver {

enum class ExecutionBackend : int {
  cpu = 0,
  metal = 1,
};

std::string to_string(ExecutionBackend backend);
ExecutionBackend parse_execution_backend(const std::string& value);

struct TaylorGreenConfig {
  int nx = 128;
  int ny = 128;
  int nz = 1;
  double viscosity = 0.01;
  double cfl_limit = 0.5;
  double final_time = 0.5;
  int poisson_max_iterations = 200;
  double poisson_tolerance = 1.0e-10;
  bool validate_energy = true;
  ExecutionBackend backend = ExecutionBackend::cpu;
  AdvectionOptions advection{};
};

struct TaylorGreenStepMetrics {
  int step = 0;
  double time = 0.0;
  double dt = 0.0;
  double max_cfl = 0.0;
  double max_velocity_change = 0.0;
  double divergence_l2 = 0.0;
  double max_divergence_l2 = 0.0;
  int pressure_iterations = 0;
  double pressure_relative_residual = 0.0;
};

struct TaylorGreenState {
  Grid grid;
  VelocityField velocity;
  VelocityField advection_previous;
  PressureField pressure_total;
  TaylorGreenStepMetrics metrics{};
  bool has_previous_advection = false;

  explicit TaylorGreenState(const Grid& grid_in)
      : grid(grid_in),
        velocity(grid_in),
        advection_previous(grid_in),
        pressure_total(grid_in) {}
};

struct TaylorGreenValidation {
  std::string reference_dataset;
  double normalized_energy_error = 0.0;
  double velocity_relative_l2_error = 0.0;
  double exact_kinetic_energy = 0.0;
  bool pass = false;
};

struct TaylorGreenResult {
  TaylorGreenConfig config{};
  ExecutionBackend backend_used = ExecutionBackend::cpu;
  std::string accelerator_name;
  double backend_elapsed_seconds = 0.0;
  double cleanup_elapsed_seconds = 0.0;
  TaylorGreenStepMetrics final_step{};
  double initial_kinetic_energy = 0.0;
  double final_kinetic_energy = 0.0;
  TaylorGreenValidation validation{};
};

[[nodiscard]] TaylorGreenConfig default_taylor_green_config();
[[nodiscard]] TaylorGreenConfig load_taylor_green_config(const std::string& path);
[[nodiscard]] std::string describe(const TaylorGreenConfig& config);
[[nodiscard]] Grid make_taylor_green_grid(const TaylorGreenConfig& config);
[[nodiscard]] BoundaryConditionSet make_taylor_green_boundary_conditions(
    const TaylorGreenConfig& config);
[[nodiscard]] double taylor_green_dt(const TaylorGreenConfig& config);
[[nodiscard]] TaylorGreenState initialize_taylor_green_state(const TaylorGreenConfig& config);
void run_taylor_green_steps(const TaylorGreenConfig& config, int step_count, TaylorGreenState& state);
[[nodiscard]] TaylorGreenResult finalize_taylor_green_result(const TaylorGreenConfig& config,
                                                             const TaylorGreenState& state);
[[nodiscard]] TaylorGreenResult run_taylor_green(const TaylorGreenConfig& config);
[[nodiscard]] TaylorGreenResult run_taylor_green(const TaylorGreenConfig& config,
                                                 TaylorGreenState* final_state);

}  // namespace solver
