#include "solver/taylor_green.hpp"

#include "metal/taylor_green_backend.hpp"
#include "operators/discrete_operators.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace solver {

namespace {

double pi() {
  return std::acos(-1.0);
}

double square(const double value) {
  return value * value;
}

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

ExecutionBackend parse_backend_value(const std::string& value) {
  if(value == "cpu") {
    return ExecutionBackend::cpu;
  }
  if(value == "metal") {
    return ExecutionBackend::metal;
  }
  throw std::runtime_error("invalid Taylor-Green backend value: " + value);
}

void validate_config(const TaylorGreenConfig& config) {
  if(config.nx <= 0 || config.ny <= 0 || config.nz <= 0) {
    throw std::invalid_argument("Taylor-Green config expects positive nx, ny, and nz");
  }
  if(config.viscosity <= 0.0) {
    throw std::invalid_argument("Taylor-Green config expects viscosity > 0");
  }
  if(config.cfl_limit <= 0.0) {
    throw std::invalid_argument("Taylor-Green config expects cfl_limit > 0");
  }
  if(config.final_time <= 0.0) {
    throw std::invalid_argument("Taylor-Green config expects final_time > 0");
  }
  if(config.poisson_max_iterations <= 0 || config.poisson_tolerance <= 0.0) {
    throw std::invalid_argument("Taylor-Green config expects positive Poisson controls");
  }
}

double domain_length() {
  return 2.0 * pi();
}

bool is_3d_case(const TaylorGreenConfig& config) {
  return config.nz > 1;
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

int axis_di(const Axis axis, const int offset) {
  return axis == Axis::x ? offset : 0;
}

int axis_dj(const Axis axis, const int offset) {
  return axis == Axis::y ? offset : 0;
}

int axis_dk(const Axis axis, const int offset) {
  return axis == Axis::z ? offset : 0;
}

double second_derivative_axis(const StructuredField& input,
                              const int i,
                              const int j,
                              const int k,
                              const Axis axis) {
  const double spacing = input.layout().grid().spacing(axis);
  const double inverse_spacing_squared = 1.0 / (spacing * spacing);
  return (input(i + axis_di(axis, 1), j + axis_dj(axis, 1), k + axis_dk(axis, 1)) -
          2.0 * input(i, j, k) +
          input(i + axis_di(axis, -1), j + axis_dj(axis, -1), k + axis_dk(axis, -1))) *
         inverse_spacing_squared;
}

double mixed_second_derivative(const StructuredField& input,
                               const int i,
                               const int j,
                               const int k,
                               const Axis outer_axis,
                               const Axis inner_axis) {
  const double spacing = input.layout().grid().spacing(outer_axis);
  const double inverse_spacing_squared = 1.0 / (spacing * spacing);
  return (second_derivative_axis(
              input,
              i + axis_di(outer_axis, 1),
              j + axis_dj(outer_axis, 1),
              k + axis_dk(outer_axis, 1),
              inner_axis) -
          2.0 * second_derivative_axis(input, i, j, k, inner_axis) +
          second_derivative_axis(
              input,
              i + axis_di(outer_axis, -1),
              j + axis_dj(outer_axis, -1),
              k + axis_dk(outer_axis, -1),
              inner_axis)) *
         inverse_spacing_squared;
}

double triple_second_derivative(const StructuredField& input,
                                const int i,
                                const int j,
                                const int k,
                                const Axis axis_a,
                                const Axis axis_b,
                                const Axis axis_c) {
  const double spacing = input.layout().grid().spacing(axis_a);
  const double inverse_spacing_squared = 1.0 / (spacing * spacing);
  return (mixed_second_derivative(
              input,
              i + axis_di(axis_a, 1),
              j + axis_dj(axis_a, 1),
              k + axis_dk(axis_a, 1),
              axis_b,
              axis_c) -
          2.0 * mixed_second_derivative(input, i, j, k, axis_b, axis_c) +
          mixed_second_derivative(
              input,
              i + axis_di(axis_a, -1),
              j + axis_dj(axis_a, -1),
              k + axis_dk(axis_a, -1),
              axis_b,
              axis_c)) *
         inverse_spacing_squared;
}

void compute_factorized_correction(const VelocityField& current_velocity,
                                   const double alpha,
                                   VelocityField& factorized_correction) {
  factorized_correction.fill(0.0);
  const bool three_dimensional = !current_velocity.x.layout().grid().is_2d();

  const auto fill_component = [alpha, three_dimensional](const StructuredField& input,
                                                         StructuredField& output) {
    const IndexRange3D active = input.layout().active_range();
    for(int k = active.k_begin; k < active.k_end; ++k) {
      for(int j = active.j_begin; j < active.j_end; ++j) {
        for(int i = active.i_begin; i < active.i_end; ++i) {
          double correction =
              mixed_second_derivative(input, i, j, k, Axis::x, Axis::y);
          if(three_dimensional) {
            correction += mixed_second_derivative(input, i, j, k, Axis::x, Axis::z);
            correction += mixed_second_derivative(input, i, j, k, Axis::y, Axis::z);
            correction -= alpha *
                          triple_second_derivative(input, i, j, k, Axis::x, Axis::y, Axis::z);
          }
          output(i, j, k) = correction;
        }
      }
    }
  };

  fill_component(current_velocity.x, factorized_correction.x);
  fill_component(current_velocity.y, factorized_correction.y);
  fill_component(current_velocity.z, factorized_correction.z);
}

void assemble_predictor_rhs(const VelocityField& current_velocity,
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
  const double alpha = 0.5 * viscosity * dt;
  compute_factorized_correction(current_velocity, alpha, factorized_correction);

  predictor_rhs = current_velocity;
  axpy_velocity(predictor_rhs, pressure_gradient, -dt);
  axpy_velocity(predictor_rhs, diffusion, 0.5 * dt);

  if(previous_advection != nullptr) {
    axpy_velocity(predictor_rhs, advection_current, -1.5 * dt);
    axpy_velocity(predictor_rhs, *previous_advection, 0.5 * dt);
  } else {
    axpy_velocity(predictor_rhs, advection_current, -dt);
  }

  axpy_velocity(predictor_rhs, factorized_correction, alpha * alpha);
}

double exact_decay(const TaylorGreenConfig& config, const double time) {
  const double rate = is_3d_case(config) ? 3.0 : 2.0;
  return std::exp(-rate * config.viscosity * time);
}

double exact_kinetic_energy(const TaylorGreenConfig& config, const double time) {
  const double initial_energy = is_3d_case(config) ? 0.125 : 0.25;
  const double rate = is_3d_case(config) ? 6.0 : 4.0;
  return initial_energy * std::exp(-rate * config.viscosity * time);
}

double kinetic_energy(const VelocityField& velocity) {
  const Grid& grid = velocity.x.layout().grid();
  const IndexRange3D cells = FieldLayout::cell_centered(grid).active_range();
  double energy_sum = 0.0;
  std::size_t count = 0;

  for(int k = cells.k_begin; k < cells.k_end; ++k) {
    for(int j = cells.j_begin; j < cells.j_end; ++j) {
      for(int i = cells.i_begin; i < cells.i_end; ++i) {
        const double u_center = 0.5 * (velocity.x(i, j, k) + velocity.x(i + 1, j, k));
        const double v_center = 0.5 * (velocity.y(i, j, k) + velocity.y(i, j + 1, k));
        const double w_center = 0.5 * (velocity.z(i, j, k) + velocity.z(i, j, k + 1));
        energy_sum += 0.5 * (square(u_center) + square(v_center) + square(w_center));
        ++count;
      }
    }
  }

  return energy_sum / static_cast<double>(count);
}

double velocity_relative_l2_error(const VelocityField& velocity,
                                  const TaylorGreenConfig& config,
                                  const double time) {
  const double decay = exact_decay(config, time);
  double numerator = 0.0;
  double denominator = 0.0;

  const auto accumulate = [&](const StructuredField& field, const auto& exact_value) {
    const IndexRange3D active = field.layout().active_range();
    for(int k = active.k_begin; k < active.k_end; ++k) {
      const double z = field.layout().coordinate_for_storage_index(Axis::z, k);
      for(int j = active.j_begin; j < active.j_end; ++j) {
        const double y = field.layout().coordinate_for_storage_index(Axis::y, j);
        for(int i = active.i_begin; i < active.i_end; ++i) {
          const double x = field.layout().coordinate_for_storage_index(Axis::x, i);
          const double reference = exact_value(x, y, z, decay);
          numerator += square(field(i, j, k) - reference);
          denominator += square(reference);
        }
      }
    }
  };

  if(is_3d_case(config)) {
    accumulate(velocity.x, [](const double x, const double y, const double z, const double decay_value) {
      return -std::cos(x) * std::sin(y) * std::cos(z) * decay_value;
    });
    accumulate(velocity.y, [](const double x, const double y, const double z, const double decay_value) {
      return std::sin(x) * std::cos(y) * std::cos(z) * decay_value;
    });
  } else {
    accumulate(velocity.x, [](const double x, const double y, double, const double decay_value) {
      return -std::cos(x) * std::sin(y) * decay_value;
    });
    accumulate(velocity.y, [](const double x, const double y, double, const double decay_value) {
      return std::sin(x) * std::cos(y) * decay_value;
    });
  }
  accumulate(velocity.z, [](double, double, double, double) {
    return 0.0;
  });

  return denominator > 0.0 ? std::sqrt(numerator / denominator) : std::sqrt(numerator);
}

void initialize_exact_state(const TaylorGreenConfig& config,
                            const BoundaryConditionSet& boundary_conditions,
                            VelocityField& velocity,
                            PressureField& pressure_total) {
  const double time = 0.0;
  const double decay = exact_decay(config, time);

  const IndexRange3D u_active = velocity.x.layout().active_range();
  for(int k = u_active.k_begin; k < u_active.k_end; ++k) {
    const double z = velocity.x.layout().coordinate_for_storage_index(Axis::z, k);
    for(int j = u_active.j_begin; j < u_active.j_end; ++j) {
      const double y = velocity.x.layout().coordinate_for_storage_index(Axis::y, j);
      for(int i = u_active.i_begin; i < u_active.i_end; ++i) {
        const double x = velocity.x.layout().coordinate_for_storage_index(Axis::x, i);
        const double z_factor = is_3d_case(config) ? std::cos(z) : 1.0;
        velocity.x(i, j, k) = -std::cos(x) * std::sin(y) * z_factor * decay;
      }
    }
  }

  const IndexRange3D v_active = velocity.y.layout().active_range();
  for(int k = v_active.k_begin; k < v_active.k_end; ++k) {
    const double z = velocity.y.layout().coordinate_for_storage_index(Axis::z, k);
    for(int j = v_active.j_begin; j < v_active.j_end; ++j) {
      const double y = velocity.y.layout().coordinate_for_storage_index(Axis::y, j);
      for(int i = v_active.i_begin; i < v_active.i_end; ++i) {
        const double x = velocity.y.layout().coordinate_for_storage_index(Axis::x, i);
        const double z_factor = is_3d_case(config) ? std::cos(z) : 1.0;
        velocity.y(i, j, k) = std::sin(x) * std::cos(y) * z_factor * decay;
      }
    }
  }

  velocity.z.fill(0.0);

  const IndexRange3D p_active = pressure_total.layout().active_range();
  const double pressure_decay =
      std::exp((is_3d_case(config) ? -6.0 : -4.0) * config.viscosity * time);
  for(int k = p_active.k_begin; k < p_active.k_end; ++k) {
    const double z = pressure_total.layout().coordinate_for_storage_index(Axis::z, k);
    for(int j = p_active.j_begin; j < p_active.j_end; ++j) {
      const double y = pressure_total.layout().coordinate_for_storage_index(Axis::y, j);
      for(int i = p_active.i_begin; i < p_active.i_end; ++i) {
        const double x = pressure_total.layout().coordinate_for_storage_index(Axis::x, i);
        if(is_3d_case(config)) {
          pressure_total(i, j, k) =
              -0.0625 * (std::cos(2.0 * x) + std::cos(2.0 * y)) *
              (std::cos(2.0 * z) + 2.0) * pressure_decay;
        } else {
          pressure_total(i, j, k) =
              -0.25 * (std::cos(2.0 * x) + std::cos(2.0 * y)) * pressure_decay;
        }
      }
    }
  }

  apply_velocity_boundary_conditions(boundary_conditions, velocity);
  VelocityField diffusion{velocity.x.layout().grid()};
  compute_diffusion_term(velocity, config.viscosity, diffusion);
  apply_total_pressure_boundary_conditions(boundary_conditions, diffusion, pressure_total);
}

double stable_dt_upper_bound(const Grid& grid, const TaylorGreenConfig& config) {
  const double inverse_sum = 1.0 / grid.dx + 1.0 / grid.dy +
                             (grid.is_2d() ? 0.0 : 1.0 / grid.dz);
  return config.cfl_limit / inverse_sum;
}

TaylorGreenValidation build_validation(const TaylorGreenConfig& config,
                                       const TaylorGreenStepMetrics& metrics,
                                       const VelocityField& velocity) {
  const double exact_energy = exact_kinetic_energy(config, metrics.time);
  const double energy_error =
      std::abs(kinetic_energy(velocity) - exact_energy) / exact_energy;
  return TaylorGreenValidation{
      .reference_dataset = is_3d_case(config) ? "analytic_taylor_green_decay_3d"
                                              : "analytic_taylor_green_decay_2d",
      .normalized_energy_error = energy_error,
      .velocity_relative_l2_error = velocity_relative_l2_error(velocity, config, metrics.time),
      .exact_kinetic_energy = exact_energy,
      .pass = energy_error <= 1.0e-2 && metrics.divergence_l2 <= 1.0e-10,
  };
}

}  // namespace

TaylorGreenConfig default_taylor_green_config() {
  return TaylorGreenConfig{};
}

TaylorGreenConfig load_taylor_green_config(const std::string& path) {
  std::ifstream input(path);
  if(!input.is_open()) {
    throw std::runtime_error("unable to open Taylor-Green config: " + path);
  }

  TaylorGreenConfig config = default_taylor_green_config();
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
    } else if(key == "nz") {
      config.nz = std::stoi(value);
    } else if(key == "viscosity") {
      config.viscosity = std::stod(value);
    } else if(key == "cfl_limit") {
      config.cfl_limit = std::stod(value);
    } else if(key == "final_time") {
      config.final_time = std::stod(value);
    } else if(key == "poisson_max_iterations") {
      config.poisson_max_iterations = std::stoi(value);
    } else if(key == "poisson_tolerance") {
      config.poisson_tolerance = std::stod(value);
    } else if(key == "validate_energy") {
      config.validate_energy = parse_bool(value);
    } else if(key == "backend") {
      config.backend = parse_backend_value(value);
    } else {
      throw std::runtime_error("unsupported Taylor-Green config key: " + key);
    }
  }

  validate_config(config);
  return config;
}

std::string describe(const TaylorGreenConfig& config) {
  std::ostringstream builder;
  builder << "nx=" << config.nx << ", ny=" << config.ny << ", nz=" << config.nz
          << ", viscosity=" << config.viscosity
          << ", cfl_limit=" << config.cfl_limit
          << ", final_time=" << config.final_time
          << ", poisson_max_iterations=" << config.poisson_max_iterations
          << ", poisson_tolerance=" << config.poisson_tolerance
          << ", validate_energy=" << (config.validate_energy ? "true" : "false")
          << ", backend=" << to_string(config.backend)
          << ", advection=" << describe(config.advection);
  return builder.str();
}

Grid make_taylor_green_grid(const TaylorGreenConfig& config) {
  validate_config(config);
  return Grid{
      config.nx,
      config.ny,
      config.nz,
      domain_length() / static_cast<double>(config.nx),
      domain_length() / static_cast<double>(config.ny),
      is_3d_case(config) ? domain_length() / static_cast<double>(config.nz) : 1.0,
      1,
  };
}

BoundaryConditionSet make_taylor_green_boundary_conditions(const TaylorGreenConfig& config) {
  BoundaryConditionSet boundary_conditions = BoundaryConditionSet::all(PhysicalBoundaryType::periodic);
  if(!is_3d_case(config)) {
    boundary_conditions[BoundaryFace::z_min].type = PhysicalBoundaryType::symmetry;
    boundary_conditions[BoundaryFace::z_max].type = PhysicalBoundaryType::symmetry;
  }
  return boundary_conditions;
}

double taylor_green_dt(const TaylorGreenConfig& config) {
  validate_config(config);
  const Grid grid = make_taylor_green_grid(config);
  const double dt_upper = stable_dt_upper_bound(grid, config);
  const int step_count = std::max(1, static_cast<int>(std::ceil(config.final_time / dt_upper)));
  return config.final_time / static_cast<double>(step_count);
}

TaylorGreenState initialize_taylor_green_state(const TaylorGreenConfig& config) {
  validate_config(config);
  const Grid grid = make_taylor_green_grid(config);
  TaylorGreenState state{grid};
  state.metrics.dt = taylor_green_dt(config);
  initialize_exact_state(config,
                         make_taylor_green_boundary_conditions(config),
                         state.velocity,
                         state.pressure_total);
  state.advection_previous.fill(0.0);
  return state;
}

void run_taylor_green_steps(const TaylorGreenConfig& config,
                            const int step_count,
                            TaylorGreenState& state) {
  validate_config(config);
  if(step_count < 0) {
    throw std::invalid_argument("run_taylor_green_steps expects a non-negative step count");
  }

  const Grid expected_grid = make_taylor_green_grid(config);
  if(state.grid.nx != expected_grid.nx || state.grid.ny != expected_grid.ny ||
     state.grid.nz != expected_grid.nz || std::abs(state.grid.dx - expected_grid.dx) > 0.0 ||
     std::abs(state.grid.dy - expected_grid.dy) > 0.0 ||
     std::abs(state.grid.dz - expected_grid.dz) > 0.0) {
    throw std::invalid_argument("Taylor-Green state grid does not match the config");
  }

  const BoundaryConditionSet boundary_conditions = make_taylor_green_boundary_conditions(config);
  const double dt = taylor_green_dt(config);
  const double predictor_alpha = 0.5 * config.viscosity * dt;
  const ProjectionOptions projection_options{
      .dt = dt,
      .density = 1.0,
      .poisson_max_iterations = config.poisson_max_iterations,
      .poisson_tolerance = config.poisson_tolerance,
  };

  VelocityField advection_current{state.grid};
  VelocityField diffusion{state.grid};
  VelocityField pressure_gradient{state.grid};
  VelocityField factorized_correction{state.grid};
  VelocityField predictor_rhs{state.grid};
  VelocityField predicted{state.grid};
  VelocityField corrected{state.grid};
  PressureField pressure_correction{state.grid};
  ScalarField pressure_rhs{state.grid};
  ProjectionDiagnostics projection{};

  for(int step = 0; step < step_count; ++step) {
    VelocityField current = state.velocity;
    apply_velocity_boundary_conditions(boundary_conditions, current);

    assemble_predictor_rhs(current,
                           state.pressure_total,
                           state.has_previous_advection ? &state.advection_previous : nullptr,
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

    axpy_active(state.pressure_total, pressure_correction, 1.0);
    axpy_active(state.pressure_total, pressure_rhs, -0.5 * config.viscosity * dt);
    compute_diffusion_term(corrected, config.viscosity, diffusion);
    apply_total_pressure_boundary_conditions(boundary_conditions, diffusion, state.pressure_total);

    const CflDiagnostics cfl = compute_advective_cfl(corrected, dt);
    const double delta = max_velocity_change(state.velocity, corrected);
    if(!std::isfinite(delta) || !std::isfinite(projection.divergence_l2_after)) {
      throw std::runtime_error("Taylor-Green solve produced a non-finite state");
    }

    state.velocity = corrected;
    state.advection_previous = advection_current;
    state.has_previous_advection = true;

    state.metrics = TaylorGreenStepMetrics{
        .step = state.metrics.step + 1,
        .time = state.metrics.time + dt,
        .dt = dt,
        .max_cfl = cfl.max_cfl,
        .max_velocity_change = delta,
        .divergence_l2 = projection.divergence_l2_after,
        .max_divergence_l2 = std::max(state.metrics.max_divergence_l2, projection.divergence_l2_after),
        .pressure_iterations = projection.pressure_solve.iterations,
        .pressure_relative_residual = projection.pressure_solve.relative_residual,
    };
  }
}

TaylorGreenResult finalize_taylor_green_result(const TaylorGreenConfig& config,
                                               const TaylorGreenState& state) {
  validate_config(config);
  TaylorGreenResult result{
      .config = config,
      .backend_used = config.backend,
      .final_step = state.metrics,
      .initial_kinetic_energy = exact_kinetic_energy(config, 0.0),
      .final_kinetic_energy = kinetic_energy(state.velocity),
  };
  if(config.validate_energy) {
    result.validation = build_validation(config, state.metrics, state.velocity);
  }
  return result;
}

std::string to_string(const ExecutionBackend backend) {
  switch(backend) {
    case ExecutionBackend::cpu:
      return "cpu";
    case ExecutionBackend::metal:
      return "metal";
  }

  return "unknown";
}

ExecutionBackend parse_execution_backend(const std::string& value) {
  return parse_backend_value(value);
}

TaylorGreenResult run_taylor_green(const TaylorGreenConfig& config,
                                   TaylorGreenState* final_state) {
  validate_config(config);
  if(config.backend == ExecutionBackend::metal) {
    metal::TaylorGreenMetalRun metal_run = metal::run_taylor_green(config);
    double cleanup_elapsed_seconds = 0.0;
    if(metal_run.state.metrics.step > 0) {
      const BoundaryConditionSet boundary_conditions =
          make_taylor_green_boundary_conditions(config);
      const ProjectionOptions projection_options{
          .dt = metal_run.state.metrics.dt,
          .density = 1.0,
          .poisson_max_iterations = config.poisson_max_iterations,
          .poisson_tolerance = config.poisson_tolerance,
      };
      PressureField pressure_correction{metal_run.state.grid};
      VelocityField corrected{metal_run.state.grid};
      pressure_correction.fill(0.0);
      corrected.fill(0.0);

      const auto cleanup_started = std::chrono::steady_clock::now();
      const ProjectionDiagnostics cleanup = project_velocity(
          metal_run.state.velocity,
          boundary_conditions,
          projection_options,
          pressure_correction,
          corrected);
      cleanup_elapsed_seconds =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - cleanup_started).count();
      metal_run.state.velocity = corrected;
      axpy_active(metal_run.state.pressure_total, pressure_correction, 1.0);
      metal_run.state.metrics.max_cfl =
          compute_advective_cfl(corrected, metal_run.state.metrics.dt).max_cfl;
      metal_run.state.metrics.divergence_l2 = cleanup.divergence_l2_after;
      metal_run.state.metrics.max_divergence_l2 =
          std::max(metal_run.state.metrics.max_divergence_l2, cleanup.divergence_l2_after);
      metal_run.state.metrics.pressure_iterations = cleanup.pressure_solve.iterations;
      metal_run.state.metrics.pressure_relative_residual =
          cleanup.pressure_solve.relative_residual;
    }

    TaylorGreenResult result = finalize_taylor_green_result(config, metal_run.state);
    result.backend_used = ExecutionBackend::metal;
    result.accelerator_name = metal_run.device_name;
    result.backend_elapsed_seconds = metal_run.elapsed_seconds;
    result.cleanup_elapsed_seconds = cleanup_elapsed_seconds;
    if(final_state != nullptr) {
      *final_state = metal_run.state;
    }
    return result;
  }

  TaylorGreenState state = initialize_taylor_green_state(config);
  run_taylor_green_steps(config, static_cast<int>(std::ceil(config.final_time / taylor_green_dt(config))), state);
  TaylorGreenResult result = finalize_taylor_green_result(config, state);
  result.backend_used = ExecutionBackend::cpu;
  result.accelerator_name = "cpu";
  if(final_state != nullptr) {
    *final_state = state;
  }
  return result;
}

TaylorGreenResult run_taylor_green(const TaylorGreenConfig& config) {
  return run_taylor_green(config, nullptr);
}

}  // namespace solver
