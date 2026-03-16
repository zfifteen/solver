#pragma once

#include "solver/momentum_terms.hpp"
#include "solver/projection.hpp"

#include <array>
#include <string>
#include <vector>

namespace solver {

struct LidDrivenCavityConfig {
  int nx = 128;
  int ny = 128;
  double reynolds = 100.0;
  double lid_velocity = 1.0;
  double cfl_limit = 0.5;
  int max_steps = 8000;
  int min_steps = 200;
  double steady_tolerance = 1.0e-7;
  int poisson_max_iterations = 200;
  double poisson_tolerance = 1.0e-10;
  bool validate_reference = true;
  AdvectionOptions advection{};
};

struct SimulationStepMetrics {
  int step = 0;
  double time = 0.0;
  double dt = 0.0;
  double max_cfl = 0.0;
  double max_velocity_change = 0.0;
  double divergence_l2 = 0.0;
  int pressure_iterations = 0;
  double pressure_relative_residual = 0.0;
};

struct CenterlineProfile {
  std::vector<double> coordinate;
  std::vector<double> value;
};

struct CenterlineExtrema {
  double u_vertical_max = 0.0;
  double u_vertical_min = 0.0;
  double v_horizontal_max = 0.0;
  double v_horizontal_min = 0.0;
};

enum class CenterlineSampleKind : int {
  u_vertical = 0,
  v_horizontal = 1,
};

struct LidDrivenCavityReferencePoint {
  std::string label;
  CenterlineSampleKind line = CenterlineSampleKind::u_vertical;
  double coordinate = 0.0;
  double value = 0.0;
  double lower_bound = 0.0;
  double upper_bound = 0.0;
};

struct LidDrivenCavityReference {
  std::string dataset;
  std::array<LidDrivenCavityReferencePoint, 4> points{};
};

struct LidDrivenCavityValidationPoint {
  std::string label;
  CenterlineSampleKind line = CenterlineSampleKind::u_vertical;
  double coordinate = 0.0;
  double reference_value = 0.0;
  double reference_lower_bound = 0.0;
  double reference_upper_bound = 0.0;
  double sample_value = 0.0;
  double relative_error = 0.0;
};

struct LidDrivenCavityValidation {
  std::string reference_dataset;
  std::array<LidDrivenCavityValidationPoint, 4> points{};
  double max_relative_error = 0.0;
  double divergence_l2 = 0.0;
  bool pass = false;
};

struct LidDrivenCavityResult {
  LidDrivenCavityConfig config{};
  SimulationStepMetrics final_step{};
  CenterlineProfile u_vertical_centerline{};
  CenterlineProfile v_horizontal_centerline{};
  CenterlineExtrema extrema{};
  LidDrivenCavityValidation validation{};
};

[[nodiscard]] LidDrivenCavityConfig default_lid_driven_cavity_config();
[[nodiscard]] LidDrivenCavityConfig load_lid_driven_cavity_config(const std::string& path);
[[nodiscard]] std::string describe(const LidDrivenCavityConfig& config);
[[nodiscard]] std::string to_string(CenterlineSampleKind line);
[[nodiscard]] BoundaryConditionSet make_lid_driven_cavity_boundary_conditions(
    const LidDrivenCavityConfig& config);
[[nodiscard]] LidDrivenCavityReference re100_centerline_reference_envelope();
[[nodiscard]] LidDrivenCavityValidation validate_lid_driven_cavity_re100(
    const LidDrivenCavityResult& result);
[[nodiscard]] LidDrivenCavityResult run_lid_driven_cavity(const LidDrivenCavityConfig& config);

namespace detail {

void apply_lid_driven_cavity_total_pressure_boundary_conditions(const VelocityField& diffusion,
                                                                PressureField& pressure_total);

void assemble_lid_driven_cavity_predictor_rhs(const VelocityField& current_velocity,
                                              const PressureField& pressure_total,
                                              const VelocityField* previous_advection,
                                              const AdvectionOptions& advection_options,
                                              double viscosity,
                                              double dt,
                                              VelocityField& advection_current,
                                              VelocityField& diffusion,
                                              VelocityField& pressure_gradient,
                                              VelocityField& factorized_correction,
                                              VelocityField& predictor_rhs);

}  // namespace detail

}  // namespace solver
