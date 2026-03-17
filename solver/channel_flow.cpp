#include "solver/channel_flow.hpp"

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

ChannelFlowCase parse_case_kind(const std::string& value) {
  if(value == "couette") {
    return ChannelFlowCase::couette;
  }
  if(value == "poiseuille") {
    return ChannelFlowCase::poiseuille;
  }
  throw std::runtime_error("unsupported channel flow case: " + value);
}

void validate_config(const ChannelFlowConfig& config) {
  if(config.nx <= 0 || config.ny <= 0) {
    throw std::invalid_argument("channel flow config expects positive nx and ny");
  }
  if(config.viscosity <= 0.0) {
    throw std::invalid_argument("channel flow config expects viscosity > 0");
  }
  if(config.cfl_limit <= 0.0) {
    throw std::invalid_argument("channel flow config expects cfl_limit > 0");
  }
  if(config.steps <= 0) {
    throw std::invalid_argument("channel flow config expects steps > 0");
  }
  if(config.poisson_max_iterations <= 0 || config.poisson_tolerance <= 0.0) {
    throw std::invalid_argument("channel flow config expects positive Poisson controls");
  }
  if(config.case_kind == ChannelFlowCase::couette && config.top_velocity <= 0.0) {
    throw std::invalid_argument("Couette config expects top_velocity > 0");
  }
  if(config.case_kind == ChannelFlowCase::poiseuille && config.pressure_drop <= 0.0) {
    throw std::invalid_argument("Poiseuille config expects pressure_drop > 0");
  }
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

double reference_max_speed(const ChannelFlowConfig& config) {
  switch(config.case_kind) {
    case ChannelFlowCase::couette:
      return config.top_velocity;
    case ChannelFlowCase::poiseuille:
      return config.pressure_drop / (8.0 * config.viscosity);
  }

  return 0.0;
}

double fixed_dt(const ChannelFlowConfig& config, const Grid& grid) {
  const double max_u = reference_max_speed(config);
  return config.cfl_limit / (max_u / grid.dx);
}

double exact_streamwise_velocity(const ChannelFlowConfig& config, const double y) {
  switch(config.case_kind) {
    case ChannelFlowCase::couette:
      return config.top_velocity * y;
    case ChannelFlowCase::poiseuille:
      return 0.5 * config.pressure_drop * y * (1.0 - y) / config.viscosity;
  }

  return 0.0;
}

double exact_total_pressure(const ChannelFlowConfig& config, const double x) {
  switch(config.case_kind) {
    case ChannelFlowCase::couette:
      return 0.0;
    case ChannelFlowCase::poiseuille:
      return config.pressure_drop * (1.0 - x);
  }

  return 0.0;
}

void initialize_steady_channel_state(const ChannelFlowConfig& config,
                                     const BoundaryConditionSet& boundary_conditions,
                                     VelocityField& velocity,
                                     PressureField& pressure_total) {
  const IndexRange3D u_active = velocity.x.layout().active_range();
  for(int k = u_active.k_begin; k < u_active.k_end; ++k) {
    for(int j = u_active.j_begin; j < u_active.j_end; ++j) {
      const double y = velocity.x.layout().coordinate_for_storage_index(Axis::y, j);
      const double streamwise = exact_streamwise_velocity(config, y);
      for(int i = u_active.i_begin; i < u_active.i_end; ++i) {
        velocity.x(i, j, k) = streamwise;
      }
    }
  }

  const IndexRange3D v_active = velocity.y.layout().active_range();
  for(int k = v_active.k_begin; k < v_active.k_end; ++k) {
    for(int j = v_active.j_begin; j < v_active.j_end; ++j) {
      for(int i = v_active.i_begin; i < v_active.i_end; ++i) {
        velocity.y(i, j, k) = 0.0;
      }
    }
  }

  const IndexRange3D w_active = velocity.z.layout().active_range();
  for(int k = w_active.k_begin; k < w_active.k_end; ++k) {
    for(int j = w_active.j_begin; j < w_active.j_end; ++j) {
      for(int i = w_active.i_begin; i < w_active.i_end; ++i) {
        velocity.z(i, j, k) = 0.0;
      }
    }
  }

  const IndexRange3D pressure_active = pressure_total.layout().active_range();
  for(int k = pressure_active.k_begin; k < pressure_active.k_end; ++k) {
    for(int j = pressure_active.j_begin; j < pressure_active.j_end; ++j) {
      for(int i = pressure_active.i_begin; i < pressure_active.i_end; ++i) {
        const double x = pressure_total.layout().coordinate_for_storage_index(Axis::x, i);
        pressure_total(i, j, k) = exact_total_pressure(config, x);
      }
    }
  }

  apply_velocity_boundary_conditions(boundary_conditions, velocity);
}

void assemble_channel_flow_predictor_rhs(const VelocityField& current_velocity,
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

WallNormalProfile sample_streamwise_profile(const VelocityField& velocity) {
  const FieldLayout& layout = velocity.x.layout();
  const IndexRange3D active = layout.active_range();
  const int count_x = active.i_end - active.i_begin;

  WallNormalProfile profile;
  profile.coordinate.reserve(static_cast<std::size_t>(layout.active_extent().ny));
  profile.value.reserve(static_cast<std::size_t>(layout.active_extent().ny));

  for(int j = active.j_begin; j < active.j_end; ++j) {
    double sum = 0.0;
    for(int i = active.i_begin; i < active.i_end; ++i) {
      sum += velocity.x(i, j, active.k_begin);
    }
    profile.coordinate.push_back(layout.coordinate_for_storage_index(Axis::y, j));
    profile.value.push_back(sum / static_cast<double>(count_x));
  }

  return profile;
}

double relative_l2_profile_error(const ChannelFlowConfig& config, const WallNormalProfile& profile) {
  if(profile.coordinate.size() != profile.value.size() || profile.coordinate.empty()) {
    throw std::invalid_argument("channel profile validation requires matching non-empty arrays");
  }

  double numerator = 0.0;
  double denominator = 0.0;
  for(std::size_t index = 0; index < profile.coordinate.size(); ++index) {
    const double reference = exact_streamwise_velocity(config, profile.coordinate[index]);
    const double difference = profile.value[index] - reference;
    numerator += difference * difference;
    denominator += reference * reference;
  }

  if(denominator == 0.0) {
    return std::sqrt(numerator);
  }

  return std::sqrt(numerator / denominator);
}

}  // namespace

ChannelFlowConfig default_channel_flow_config() {
  return ChannelFlowConfig{};
}

ChannelFlowConfig load_channel_flow_config(const std::string& path) {
  std::ifstream input(path);
  if(!input.is_open()) {
    throw std::runtime_error("unable to open channel flow config: " + path);
  }

  ChannelFlowConfig config = default_channel_flow_config();
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

    if(key == "case") {
      config.case_kind = parse_case_kind(value);
    } else if(key == "nx") {
      config.nx = std::stoi(value);
    } else if(key == "ny") {
      config.ny = std::stoi(value);
    } else if(key == "viscosity") {
      config.viscosity = std::stod(value);
    } else if(key == "top_velocity") {
      config.top_velocity = std::stod(value);
    } else if(key == "pressure_drop") {
      config.pressure_drop = std::stod(value);
    } else if(key == "cfl_limit") {
      config.cfl_limit = std::stod(value);
    } else if(key == "steps") {
      config.steps = std::stoi(value);
    } else if(key == "poisson_max_iterations") {
      config.poisson_max_iterations = std::stoi(value);
    } else if(key == "poisson_tolerance") {
      config.poisson_tolerance = std::stod(value);
    } else if(key == "validate_profile") {
      config.validate_profile = parse_bool(value);
    } else {
      throw std::runtime_error("unsupported channel flow config key: " + key);
    }
  }

  validate_config(config);
  return config;
}

std::string to_string(const ChannelFlowCase case_kind) {
  switch(case_kind) {
    case ChannelFlowCase::couette:
      return "couette";
    case ChannelFlowCase::poiseuille:
      return "poiseuille";
  }

  return "unknown";
}

std::string describe(const ChannelFlowConfig& config) {
  std::ostringstream builder;
  builder << "case=" << to_string(config.case_kind) << ", nx=" << config.nx << ", ny=" << config.ny
          << ", viscosity=" << config.viscosity << ", top_velocity=" << config.top_velocity
          << ", pressure_drop=" << config.pressure_drop << ", cfl_limit=" << config.cfl_limit
          << ", steps=" << config.steps
          << ", poisson_max_iterations=" << config.poisson_max_iterations
          << ", poisson_tolerance=" << config.poisson_tolerance
          << ", validate_profile=" << (config.validate_profile ? "true" : "false")
          << ", advection=" << describe(config.advection);
  return builder.str();
}

BoundaryConditionSet make_channel_flow_boundary_conditions(const ChannelFlowConfig& config) {
  BoundaryConditionSet boundary_conditions = BoundaryConditionSet::all(PhysicalBoundaryType::no_slip_wall);
  boundary_conditions[BoundaryFace::z_min].type = PhysicalBoundaryType::symmetry;
  boundary_conditions[BoundaryFace::z_max].type = PhysicalBoundaryType::symmetry;

  switch(config.case_kind) {
    case ChannelFlowCase::couette:
      boundary_conditions[BoundaryFace::x_min].type = PhysicalBoundaryType::periodic;
      boundary_conditions[BoundaryFace::x_max].type = PhysicalBoundaryType::periodic;
      boundary_conditions[BoundaryFace::y_max].type = PhysicalBoundaryType::prescribed_velocity;
      boundary_conditions[BoundaryFace::y_max].velocity = {config.top_velocity, 0.0, 0.0};
      break;
    case ChannelFlowCase::poiseuille:
      boundary_conditions[BoundaryFace::x_min].type = PhysicalBoundaryType::fixed_pressure;
      boundary_conditions[BoundaryFace::x_min].pressure = config.pressure_drop;
      boundary_conditions[BoundaryFace::x_max].type = PhysicalBoundaryType::fixed_pressure;
      boundary_conditions[BoundaryFace::x_max].pressure = 0.0;
      break;
  }

  return boundary_conditions;
}

ChannelFlowValidation validate_channel_flow_profile(const ChannelFlowResult& result) {
  const ChannelFlowConfig& config = result.config;
  const double relative_error = relative_l2_profile_error(config, result.streamwise_profile);
  ChannelFlowValidation validation{
      .reference_dataset = config.case_kind == ChannelFlowCase::couette ? "analytic_couette_profile"
                                                                        : "analytic_poiseuille_profile",
      .relative_l2_error = relative_error,
      .divergence_l2 = result.final_step.divergence_l2,
      .pass = relative_error <= 5.0e-3 && result.final_step.divergence_l2 <= 1.0e-10,
  };
  return validation;
}

ChannelFlowResult run_channel_flow(const ChannelFlowConfig& config) {
  validate_config(config);

  const Grid grid{config.nx,
                  config.ny,
                  1,
                  1.0 / static_cast<double>(config.nx),
                  1.0 / static_cast<double>(config.ny),
                  1.0,
                  1};
  const BoundaryConditionSet boundary_conditions = make_channel_flow_boundary_conditions(config);
  const double dt = fixed_dt(config, grid);
  const double predictor_alpha = 0.5 * config.viscosity * dt;
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
  PressureField pressure_total{grid};
  PressureField pressure_correction{grid};
  ScalarField pressure_rhs{grid};
  ProjectionDiagnostics projection{};
  bool has_previous_advection = false;
  double time = 0.0;
  ChannelFlowStepMetrics metrics{.dt = dt};

  velocity.fill(0.0);
  pressure_total.fill(0.0);
  pressure_correction.fill(0.0);
  initialize_steady_channel_state(config, boundary_conditions, velocity, pressure_total);
  compute_diffusion_term(velocity, config.viscosity, diffusion);
  apply_total_pressure_boundary_conditions(boundary_conditions, diffusion, pressure_total);

  for(int step = 0; step < config.steps; ++step) {
    VelocityField current = velocity;
    apply_velocity_boundary_conditions(boundary_conditions, current);

    assemble_channel_flow_predictor_rhs(current,
                                        pressure_total,
                                        has_previous_advection ? &advection_previous : nullptr,
                                        config.advection,
                                        config.viscosity,
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

    axpy_active(pressure_total, pressure_correction, 1.0);
    axpy_active(pressure_total, pressure_rhs, -0.5 * config.viscosity * dt);
    compute_diffusion_term(corrected, config.viscosity, diffusion);
    apply_total_pressure_boundary_conditions(boundary_conditions, diffusion, pressure_total);

    const CflDiagnostics cfl = compute_advective_cfl(corrected, dt);
    const double delta = max_velocity_change(velocity, corrected);
    if(!std::isfinite(delta) || !std::isfinite(projection.divergence_l2_after)) {
      throw std::runtime_error("channel flow solve produced a non-finite state");
    }

    velocity = corrected;
    advection_previous = advection_current;
    has_previous_advection = true;
    time += dt;

    metrics = ChannelFlowStepMetrics{
        .step = step + 1,
        .time = time,
        .dt = dt,
        .max_cfl = cfl.max_cfl,
        .max_velocity_change = delta,
        .divergence_l2 = projection.divergence_l2_after,
        .pressure_iterations = projection.pressure_solve.iterations,
        .pressure_relative_residual = projection.pressure_solve.relative_residual,
    };
  }

  apply_velocity_boundary_conditions(boundary_conditions, velocity);

  ChannelFlowResult result{
      .config = config,
      .final_step = metrics,
      .streamwise_profile = sample_streamwise_profile(velocity),
  };
  if(config.validate_profile) {
    result.validation = validate_channel_flow_profile(result);
  }
  return result;
}

}  // namespace solver
