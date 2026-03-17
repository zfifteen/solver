#include "solver/lid_driven_cavity.hpp"

#include "operators/discrete_operators.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace solver {

namespace {

std::string trim(const std::string& value) {
  const std::size_t first = value.find_first_not_of(" \t\r\n");
  if(first == std::string::npos) {
    return "";
  }

  const std::size_t last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

bool parse_bool(const std::string& value) {
  if(value == "true" || value == "1" || value == "yes" || value == "on") {
    return true;
  }
  if(value == "false" || value == "0" || value == "no" || value == "off") {
    return false;
  }
  throw std::runtime_error("invalid boolean value: " + value);
}

void axpy_active(StructuredField& destination, const StructuredField& source, const double scale) {
  const IndexRange3D active = destination.layout().active_range();

  for(int k = active.k_begin; k < active.k_end; ++k) {
    for(int j = active.j_begin; j < active.j_end; ++j) {
      for(int i = active.i_begin; i < active.i_end; ++i) {
        destination(i, j, k) += scale * source(i, j, k);
      }
    }
  }
}

void axpy_velocity(VelocityField& destination, const VelocityField& source, const double scale) {
  axpy_active(destination.x, source.x, scale);
  axpy_active(destination.y, source.y, scale);
  axpy_active(destination.z, source.z, scale);
}

double max_velocity_change(const VelocityField& left, const VelocityField& right) {
  double change = 0.0;

  const auto accumulate = [&change](const StructuredField& lhs, const StructuredField& rhs) {
    const IndexRange3D active = lhs.layout().active_range();
    for(int k = active.k_begin; k < active.k_end; ++k) {
      for(int j = active.j_begin; j < active.j_end; ++j) {
        for(int i = active.i_begin; i < active.i_end; ++i) {
          change = std::max(change, std::abs(lhs(i, j, k) - rhs(i, j, k)));
        }
      }
    }
  };

  accumulate(left.x, right.x);
  accumulate(left.y, right.y);
  accumulate(left.z, right.z);
  return change;
}

void validate_config(const LidDrivenCavityConfig& config) {
  if(config.nx <= 0 || config.ny <= 0) {
    throw std::invalid_argument("lid-driven cavity config expects positive nx and ny");
  }
  if(config.reynolds <= 0.0) {
    throw std::invalid_argument("lid-driven cavity config expects reynolds > 0");
  }
  if(config.lid_velocity <= 0.0) {
    throw std::invalid_argument("lid-driven cavity config expects lid_velocity > 0");
  }
  if(config.cfl_limit <= 0.0) {
    throw std::invalid_argument("lid-driven cavity config expects cfl_limit > 0");
  }
  if(config.max_steps <= 0 || config.min_steps < 0) {
    throw std::invalid_argument("lid-driven cavity config expects non-negative step controls");
  }
  if(config.poisson_max_iterations <= 0 || config.poisson_tolerance <= 0.0) {
    throw std::invalid_argument("lid-driven cavity config expects positive Poisson controls");
  }
}

double fixed_dt(const LidDrivenCavityConfig& config, const Grid& grid) {
  const double denominator = config.lid_velocity / grid.dx + config.lid_velocity / grid.dy;
  return config.cfl_limit / denominator;
}

CenterlineProfile sample_vertical_centerline_u(const VelocityField& velocity) {
  const FieldLayout& layout = velocity.x.layout();
  const IndexRange3D active = layout.active_range();
  const int i_center = active.i_begin + layout.active_extent().nx / 2;

  CenterlineProfile profile;
  profile.coordinate.reserve(static_cast<std::size_t>(layout.active_extent().ny));
  profile.value.reserve(static_cast<std::size_t>(layout.active_extent().ny));

  for(int j = active.j_begin; j < active.j_end; ++j) {
    profile.coordinate.push_back(layout.coordinate_for_storage_index(Axis::y, j));
    profile.value.push_back(velocity.x(i_center, j, active.k_begin));
  }

  return profile;
}

CenterlineProfile sample_horizontal_centerline_v(const VelocityField& velocity) {
  const FieldLayout& layout = velocity.y.layout();
  const IndexRange3D active = layout.active_range();
  const int j_center = active.j_begin + layout.active_extent().ny / 2;

  CenterlineProfile profile;
  profile.coordinate.reserve(static_cast<std::size_t>(layout.active_extent().nx));
  profile.value.reserve(static_cast<std::size_t>(layout.active_extent().nx));

  for(int i = active.i_begin; i < active.i_end; ++i) {
    profile.coordinate.push_back(layout.coordinate_for_storage_index(Axis::x, i));
    profile.value.push_back(velocity.y(i, j_center, active.k_begin));
  }

  return profile;
}

CenterlineExtrema compute_centerline_extrema(const CenterlineProfile& u_profile,
                                             const CenterlineProfile& v_profile) {
  if(u_profile.value.size() < 3 || v_profile.value.size() < 3) {
    throw std::invalid_argument("centerline extrema require at least three samples per profile");
  }

  const auto u_begin = u_profile.value.begin() + 1;
  const auto u_end = u_profile.value.end() - 1;
  const auto v_begin = v_profile.value.begin() + 1;
  const auto v_end = v_profile.value.end() - 1;
  const auto [u_min_it, u_max_it] = std::minmax_element(u_begin, u_end);
  const auto [v_min_it, v_max_it] = std::minmax_element(v_begin, v_end);

  return CenterlineExtrema{
      .u_vertical_max = *u_max_it,
      .u_vertical_min = *u_min_it,
      .v_horizontal_max = *v_max_it,
      .v_horizontal_min = *v_min_it,
  };
}

double relative_error(const double value, const double reference) {
  if(std::abs(reference) == 0.0) {
    return std::abs(value);
  }
  return std::abs(value - reference) / std::abs(reference);
}

double relative_miss_against_envelope(const double value,
                                      const double lower_bound,
                                      const double upper_bound) {
  const double lower = std::min(lower_bound, upper_bound);
  const double upper = std::max(lower_bound, upper_bound);
  if(value >= lower && value <= upper) {
    return 0.0;
  }

  const double nearest = value < lower ? lower : upper;
  return relative_error(value, nearest);
}

double interpolate_profile_value(const CenterlineProfile& profile, const double coordinate) {
  if(profile.coordinate.size() != profile.value.size() || profile.coordinate.size() < 2) {
    throw std::invalid_argument("centerline profile interpolation requires matching sample arrays");
  }
  if(coordinate < profile.coordinate.front() || coordinate > profile.coordinate.back()) {
    throw std::out_of_range("reference coordinate falls outside the sampled centerline range");
  }

  const auto upper = std::lower_bound(profile.coordinate.begin(), profile.coordinate.end(), coordinate);
  if(upper == profile.coordinate.begin()) {
    return profile.value.front();
  }
  if(upper == profile.coordinate.end()) {
    return profile.value.back();
  }
  if(*upper == coordinate) {
    return profile.value[static_cast<std::size_t>(upper - profile.coordinate.begin())];
  }

  const std::size_t upper_index = static_cast<std::size_t>(upper - profile.coordinate.begin());
  const std::size_t lower_index = upper_index - 1;
  if(profile.coordinate.size() >= 4) {
    const std::size_t start =
        std::min(lower_index > 0 ? lower_index - 1 : 0, profile.coordinate.size() - 4);
    double interpolated = 0.0;
    for(std::size_t i = 0; i < 4; ++i) {
      const std::size_t node = start + i;
      double basis = 1.0;
      for(std::size_t j = 0; j < 4; ++j) {
        if(i == j) {
          continue;
        }
        const std::size_t other = start + j;
        basis *= (coordinate - profile.coordinate[other]) /
                 (profile.coordinate[node] - profile.coordinate[other]);
      }
      interpolated += basis * profile.value[node];
    }
    return interpolated;
  }

  const double lower_coordinate = profile.coordinate[lower_index];
  const double upper_coordinate = profile.coordinate[upper_index];
  const double weight = (coordinate - lower_coordinate) / (upper_coordinate - lower_coordinate);
  return (1.0 - weight) * profile.value[lower_index] + weight * profile.value[upper_index];
}

double mixed_second_derivative_2d(const StructuredField& input, const int i, const int j, const int k) {
  const double dy = input.layout().grid().dy;
  const double dy_squared = dy * dy;
  const double dyy_im1 =
      (input(i - 1, j + 1, k) - 2.0 * input(i - 1, j, k) + input(i - 1, j - 1, k)) / dy_squared;
  const double dyy_i =
      (input(i, j + 1, k) - 2.0 * input(i, j, k) + input(i, j - 1, k)) / dy_squared;
  const double dyy_ip1 =
      (input(i + 1, j + 1, k) - 2.0 * input(i + 1, j, k) + input(i + 1, j - 1, k)) / dy_squared;

  const double dx = input.layout().grid().dx;
  const double inverse_dx_squared = 1.0 / (dx * dx);
  return (dyy_ip1 - 2.0 * dyy_i + dyy_im1) * inverse_dx_squared;
}

void compute_factorized_correction_2d(const VelocityField& current_velocity,
                                      VelocityField& factorized_correction) {
  factorized_correction.fill(0.0);

  const auto fill_component = [](const StructuredField& input, StructuredField& output) {
    const IndexRange3D active = input.layout().active_range();
    for(int k = active.k_begin; k < active.k_end; ++k) {
      for(int j = active.j_begin; j < active.j_end; ++j) {
        for(int i = active.i_begin; i < active.i_end; ++i) {
          output(i, j, k) = mixed_second_derivative_2d(input, i, j, k);
        }
      }
    }
  };

  fill_component(current_velocity.x, factorized_correction.x);
  fill_component(current_velocity.y, factorized_correction.y);
  fill_component(current_velocity.z, factorized_correction.z);
}

bool same_grid(const Grid& left, const Grid& right) {
  return left.nx == right.nx && left.ny == right.ny && left.nz == right.nz &&
         left.dx == right.dx && left.dy == right.dy && left.dz == right.dz &&
         left.ghost_layers == right.ghost_layers;
}

bool steady_stop_reached(const LidDrivenCavityConfig& config,
                         const SimulationStepMetrics& metrics) {
  return metrics.step >= config.min_steps &&
         metrics.max_velocity_change <= config.steady_tolerance;
}

}  // namespace

LidDrivenCavityState::LidDrivenCavityState(const Grid& grid_in)
    : grid(grid_in),
      velocity(grid),
      advection_previous(grid),
      pressure_total(grid) {}

LidDrivenCavityConfig default_lid_driven_cavity_config() {
  return LidDrivenCavityConfig{};
}

LidDrivenCavityConfig load_lid_driven_cavity_config(const std::string& path) {
  std::ifstream input(path);
  if(!input.is_open()) {
    throw std::runtime_error("unable to open cavity config: " + path);
  }

  LidDrivenCavityConfig config = default_lid_driven_cavity_config();
  std::string line;
  int line_number = 0;
  while(std::getline(input, line)) {
    ++line_number;
    const std::size_t comment = line.find('#');
    if(comment != std::string::npos) {
      line = line.substr(0, comment);
    }

    line = trim(line);
    if(line.empty()) {
      continue;
    }

    const std::size_t separator = line.find('=');
    if(separator == std::string::npos) {
      throw std::runtime_error("invalid config line " + std::to_string(line_number) + ": " + line);
    }

    const std::string key = trim(line.substr(0, separator));
    const std::string value = trim(line.substr(separator + 1));

    if(key == "nx") {
      config.nx = std::stoi(value);
    } else if(key == "ny") {
      config.ny = std::stoi(value);
    } else if(key == "reynolds") {
      config.reynolds = std::stod(value);
    } else if(key == "lid_velocity") {
      config.lid_velocity = std::stod(value);
    } else if(key == "cfl_limit") {
      config.cfl_limit = std::stod(value);
    } else if(key == "max_steps") {
      config.max_steps = std::stoi(value);
    } else if(key == "min_steps") {
      config.min_steps = std::stoi(value);
    } else if(key == "steady_tolerance") {
      config.steady_tolerance = std::stod(value);
    } else if(key == "poisson_max_iterations") {
      config.poisson_max_iterations = std::stoi(value);
    } else if(key == "poisson_tolerance") {
      config.poisson_tolerance = std::stod(value);
    } else if(key == "validate_reference") {
      config.validate_reference = parse_bool(value);
    } else {
      throw std::runtime_error("unsupported cavity config key: " + key);
    }
  }

  validate_config(config);
  return config;
}

std::string describe(const LidDrivenCavityConfig& config) {
  std::ostringstream builder;
  builder << "nx=" << config.nx << ", ny=" << config.ny << ", Re=" << config.reynolds
          << ", lid_velocity=" << config.lid_velocity << ", cfl_limit=" << config.cfl_limit
          << ", max_steps=" << config.max_steps << ", min_steps=" << config.min_steps
          << ", steady_tolerance=" << config.steady_tolerance
          << ", poisson_max_iterations=" << config.poisson_max_iterations
          << ", poisson_tolerance=" << config.poisson_tolerance
          << ", validate_reference=" << (config.validate_reference ? "true" : "false")
          << ", advection=" << describe(config.advection);
  return builder.str();
}

std::string to_string(const CenterlineSampleKind line) {
  switch(line) {
    case CenterlineSampleKind::u_vertical:
      return "u_vertical";
    case CenterlineSampleKind::v_horizontal:
      return "v_horizontal";
  }

  __builtin_unreachable();
}

BoundaryConditionSet make_lid_driven_cavity_boundary_conditions(
    const LidDrivenCavityConfig& config) {
  BoundaryConditionSet boundary_conditions = BoundaryConditionSet::cavity();
  boundary_conditions[BoundaryFace::y_max].type = PhysicalBoundaryType::prescribed_velocity;
  boundary_conditions[BoundaryFace::y_max].velocity = {config.lid_velocity, 0.0, 0.0};
  boundary_conditions[BoundaryFace::z_min].type = PhysicalBoundaryType::symmetry;
  boundary_conditions[BoundaryFace::z_max].type = PhysicalBoundaryType::symmetry;
  return boundary_conditions;
}

Grid make_lid_driven_cavity_grid(const LidDrivenCavityConfig& config) {
  return Grid{config.nx,
              config.ny,
              1,
              1.0 / static_cast<double>(config.nx),
              1.0 / static_cast<double>(config.ny),
              1.0,
              1};
}

double lid_driven_cavity_viscosity(const LidDrivenCavityConfig& config) {
  return config.lid_velocity / config.reynolds;
}

double lid_driven_cavity_dt(const LidDrivenCavityConfig& config) {
  return fixed_dt(config, make_lid_driven_cavity_grid(config));
}

LidDrivenCavityState initialize_lid_driven_cavity_state(const LidDrivenCavityConfig& config) {
  validate_config(config);

  const Grid grid = make_lid_driven_cavity_grid(config);
  const BoundaryConditionSet boundary_conditions =
      make_lid_driven_cavity_boundary_conditions(config);
  const double viscosity = lid_driven_cavity_viscosity(config);

  LidDrivenCavityState state{grid};
  state.metrics.dt = lid_driven_cavity_dt(config);
  state.velocity.fill(0.0);
  state.advection_previous.fill(0.0);
  state.pressure_total.fill(0.0);
  apply_velocity_boundary_conditions(boundary_conditions, state.velocity);

  VelocityField diffusion{grid};
  compute_diffusion_term(state.velocity, viscosity, diffusion);
  apply_total_pressure_boundary_conditions(boundary_conditions, diffusion, state.pressure_total);
  return state;
}

LidDrivenCavityReference re100_centerline_reference_envelope() {
  return LidDrivenCavityReference{
      .dataset = "Re100_centerline_literature_envelope",
      .points = {{
          {
              .label = "u_vertical_y0.9766",
              .line = CenterlineSampleKind::u_vertical,
              .coordinate = 0.9766,
              .value = 0.84123,
              .lower_bound = 0.84123,
              .upper_bound = 0.84123,
          },
          {
              .label = "u_vertical_y0.4531",
              .line = CenterlineSampleKind::u_vertical,
              .coordinate = 0.4531,
              .value = -0.21090,
              .lower_bound = -0.21400,
              .upper_bound = -0.21060,
          },
          {
              .label = "v_horizontal_x0.2344",
              .line = CenterlineSampleKind::v_horizontal,
              .coordinate = 0.2344,
              .value = 0.17527,
              .lower_bound = 0.17527,
              .upper_bound = 0.17960,
          },
          {
              .label = "v_horizontal_x0.8047",
              .line = CenterlineSampleKind::v_horizontal,
              .coordinate = 0.8047,
              .value = -0.24533,
              .lower_bound = -0.25400,
              .upper_bound = -0.24533,
          },
      }},
  };
}

LidDrivenCavityValidation validate_lid_driven_cavity_re100(const LidDrivenCavityResult& result) {
  const LidDrivenCavityReference reference = re100_centerline_reference_envelope();
  LidDrivenCavityValidation validation{
      .reference_dataset = reference.dataset,
      .divergence_l2 = result.final_step.divergence_l2,
  };

  double max_error = 0.0;
  bool pass = result.final_step.divergence_l2 <= 1.0e-10;

  for(std::size_t point_index = 0; point_index < reference.points.size(); ++point_index) {
    const LidDrivenCavityReferencePoint& point = reference.points[point_index];
    const CenterlineProfile& profile = point.line == CenterlineSampleKind::u_vertical
                                           ? result.u_vertical_centerline
                                           : result.v_horizontal_centerline;
    const double sample_value = interpolate_profile_value(profile, point.coordinate);

    validation.points[point_index] = LidDrivenCavityValidationPoint{
        .label = point.label,
        .line = point.line,
        .coordinate = point.coordinate,
        .reference_value = point.value,
        .reference_lower_bound = point.lower_bound,
        .reference_upper_bound = point.upper_bound,
        .sample_value = sample_value,
        .relative_error = relative_miss_against_envelope(sample_value, point.lower_bound, point.upper_bound),
    };
    max_error = std::max(max_error, validation.points[point_index].relative_error);
    pass = pass && validation.points[point_index].relative_error <= 0.02;
  }

  validation.max_relative_error = max_error;
  validation.pass = pass;
  return validation;
}

namespace detail {

void assemble_lid_driven_cavity_predictor_rhs(const VelocityField& current_velocity,
                                              const PressureField& pressure_total,
                                              const VelocityField* previous_advection,
                                              const AdvectionOptions& advection_options,
                                              const double viscosity,
                                              const double dt,
                                              VelocityField& advection_current,
                                              VelocityField& diffusion,
                                              VelocityField& pressure_gradient,
                                              VelocityField& factorized_correction,
                                              VelocityField& predictor_rhs) {
  compute_advection_term(current_velocity, advection_options, advection_current);
  compute_diffusion_term(current_velocity, viscosity, diffusion);
  operators::compute_gradient(pressure_total, pressure_gradient);
  compute_factorized_correction_2d(current_velocity, factorized_correction);

  predictor_rhs = current_velocity;
  axpy_velocity(predictor_rhs, pressure_gradient, -dt);
  axpy_velocity(predictor_rhs, diffusion, 0.5 * dt);

  if(previous_advection != nullptr) {
    axpy_velocity(predictor_rhs, advection_current, -1.5 * dt);
    axpy_velocity(predictor_rhs, *previous_advection, 0.5 * dt);
  } else {
    axpy_velocity(predictor_rhs, advection_current, -dt);
  }

  const double alpha = 0.5 * viscosity * dt;
  axpy_velocity(predictor_rhs, factorized_correction, alpha * alpha);
}

}  // namespace detail

namespace {

void advance_lid_driven_cavity_impl(const LidDrivenCavityConfig& config,
                                    const int step_count,
                                    const bool stop_on_steady,
                                    LidDrivenCavityState& state) {
  validate_config(config);
  if(step_count < 0) {
    throw std::invalid_argument("run_lid_driven_cavity_steps expects a non-negative step count");
  }

  const Grid grid = make_lid_driven_cavity_grid(config);
  if(!same_grid(state.grid, grid)) {
    throw std::invalid_argument("lid-driven cavity state grid does not match the config");
  }

  const BoundaryConditionSet boundary_conditions = make_lid_driven_cavity_boundary_conditions(config);
  const double viscosity = lid_driven_cavity_viscosity(config);
  const double dt = lid_driven_cavity_dt(config);
  const double predictor_alpha = 0.5 * viscosity * dt;
  const ProjectionOptions projection_options{
      .dt = dt,
      .density = 1.0,
      .poisson_max_iterations = config.poisson_max_iterations,
      .poisson_tolerance = config.poisson_tolerance,
  };

  VelocityField velocity{grid};
  VelocityField advection_current{grid};
  VelocityField advection_previous{grid};
  VelocityField diffusion{grid};
  VelocityField pressure_gradient{grid};
  VelocityField factorized_correction{grid};
  VelocityField predictor_rhs{grid};
  VelocityField predicted{grid};
  VelocityField corrected{grid};
  PressureField pressure_correction{grid};
  ScalarField pressure_rhs{grid};
  ProjectionDiagnostics projection{};

  for(int step = 0; step < step_count; ++step) {
    VelocityField current = state.velocity;
    apply_velocity_boundary_conditions(boundary_conditions, current);

    detail::assemble_lid_driven_cavity_predictor_rhs(current,
                                                     state.pressure_total,
                                                     state.has_previous_advection
                                                         ? &state.advection_previous
                                                         : nullptr,
                                                     config.advection,
                                                     viscosity,
                                                     dt,
                                                     advection_current,
                                                     diffusion,
                                                     pressure_gradient,
                                                     factorized_correction,
                                                     predictor_rhs);

    predicted = predictor_rhs;
    solve_predictor_adi(predictor_rhs, predictor_alpha, boundary_conditions, predicted);

    pressure_correction.fill(0.0);
    projection = project_velocity(
        predicted,
        boundary_conditions,
        projection_options,
        pressure_correction,
        corrected,
        &pressure_rhs);

    axpy_active(state.pressure_total, pressure_correction, 1.0);
    axpy_active(state.pressure_total, pressure_rhs, -0.5 * viscosity * dt);
    compute_diffusion_term(corrected, viscosity, diffusion);
    apply_total_pressure_boundary_conditions(boundary_conditions, diffusion, state.pressure_total);

    const CflDiagnostics cfl = compute_advective_cfl(corrected, dt);
    const double delta = max_velocity_change(state.velocity, corrected);
    if(!std::isfinite(delta) || !std::isfinite(projection.divergence_l2_after)) {
      throw std::runtime_error("lid-driven cavity solve produced a non-finite state");
    }

    state.velocity = corrected;
    state.advection_previous = advection_current;
    state.has_previous_advection = true;

    state.metrics = SimulationStepMetrics{
        .step = state.metrics.step + 1,
        .time = state.metrics.time + dt,
        .dt = dt,
        .max_cfl = cfl.max_cfl,
        .max_velocity_change = delta,
        .divergence_l2 = projection.divergence_l2_after,
        .pressure_iterations = projection.pressure_solve.iterations,
        .pressure_relative_residual = projection.pressure_solve.relative_residual,
    };

    if(stop_on_steady && steady_stop_reached(config, state.metrics)) {
      break;
    }
  }
}

}  // namespace

void run_lid_driven_cavity_steps(const LidDrivenCavityConfig& config,
                                 const int step_count,
                                 LidDrivenCavityState& state) {
  advance_lid_driven_cavity_impl(config, step_count, false, state);
}

LidDrivenCavityResult finalize_lid_driven_cavity_result(const LidDrivenCavityConfig& config,
                                                        const LidDrivenCavityState& state) {
  validate_config(config);
  const Grid grid = make_lid_driven_cavity_grid(config);
  if(!same_grid(state.grid, grid)) {
    throw std::invalid_argument("lid-driven cavity state grid does not match the config");
  }

  const BoundaryConditionSet boundary_conditions = make_lid_driven_cavity_boundary_conditions(config);
  VelocityField velocity = state.velocity;
  apply_velocity_boundary_conditions(boundary_conditions, velocity);

  LidDrivenCavityResult result{
      .config = config,
      .final_step = state.metrics,
      .u_vertical_centerline = sample_vertical_centerline_u(velocity),
      .v_horizontal_centerline = sample_horizontal_centerline_v(velocity),
  };
  result.extrema = compute_centerline_extrema(result.u_vertical_centerline, result.v_horizontal_centerline);
  if(config.validate_reference && std::abs(config.reynolds - 100.0) < 1.0e-12) {
    result.validation = validate_lid_driven_cavity_re100(result);
  }
  return result;
}

LidDrivenCavityResult run_lid_driven_cavity(const LidDrivenCavityConfig& config) {
  LidDrivenCavityState state = initialize_lid_driven_cavity_state(config);
  advance_lid_driven_cavity_impl(config, config.max_steps - state.metrics.step, true, state);
  return finalize_lid_driven_cavity_result(config, state);
}

}  // namespace solver
