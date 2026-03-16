#include "io/checkpoint.hpp"
#include "io/vtk_export.hpp"
#include "solver/lid_driven_cavity.hpp"

#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

struct CliOptions {
  std::string config_path = "benchmarks/lid_driven_cavity_re100_128.cfg";
  std::optional<int> steps_override{};
  std::string restart_path{};
  std::string checkpoint_out{};
  std::string vtk_out{};
};

CliOptions parse_args(int argc, char** argv) {
  CliOptions options;
  bool config_consumed = false;

  for(int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    const auto require_value = [&](const std::string& flag) -> std::string {
      if(index + 1 >= argc) {
        throw std::runtime_error("missing value for " + flag);
      }
      return argv[++index];
    };

    if(argument == "--steps") {
      options.steps_override = std::stoi(require_value(argument));
    } else if(argument == "--restart") {
      options.restart_path = require_value(argument);
    } else if(argument == "--checkpoint-out") {
      options.checkpoint_out = require_value(argument);
    } else if(argument == "--vtk-out") {
      options.vtk_out = require_value(argument);
    } else if(!config_consumed && !argument.empty() && argument.rfind("--", 0) != 0) {
      options.config_path = argument;
      config_consumed = true;
    } else {
      throw std::runtime_error("unsupported argument: " + argument);
    }
  }

  return options;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const CliOptions options = parse_args(argc, argv);
    const solver::LidDrivenCavityConfig config =
        solver::load_lid_driven_cavity_config(options.config_path);
    solver::LidDrivenCavityState state =
        options.restart_path.empty()
            ? solver::initialize_lid_driven_cavity_state(config)
            : solver::io::load_lid_driven_cavity_checkpoint(options.restart_path, config).state;

    const bool partial_run = options.steps_override.has_value();
    if(partial_run) {
      solver::run_lid_driven_cavity_steps(config, *options.steps_override, state);
    } else {
      while(state.metrics.step < config.max_steps) {
        solver::run_lid_driven_cavity_steps(config, 1, state);
        if(state.metrics.step >= config.min_steps &&
           state.metrics.max_velocity_change <= config.steady_tolerance) {
          break;
        }
      }
    }

    if(!options.checkpoint_out.empty()) {
      solver::io::write_lid_driven_cavity_checkpoint(options.checkpoint_out, config, state);
    }
    if(!options.vtk_out.empty()) {
      solver::io::write_lid_driven_cavity_vtk(options.vtk_out, state);
    }

    const solver::LidDrivenCavityResult result =
        solver::finalize_lid_driven_cavity_result(config, state);

    std::cout << "config_path: " << options.config_path << '\n';
    std::cout << "config: " << solver::describe(config) << '\n';
    std::cout << "steps: " << result.final_step.step << '\n';
    std::cout << "time: " << result.final_step.time << '\n';
    std::cout << "dt: " << result.final_step.dt << '\n';
    std::cout << "max_cfl: " << result.final_step.max_cfl << '\n';
    std::cout << "max_velocity_change: " << result.final_step.max_velocity_change << '\n';
    std::cout << "divergence_l2: " << result.final_step.divergence_l2 << '\n';
    std::cout << "pressure_iterations: " << result.final_step.pressure_iterations << '\n';
    std::cout << "pressure_relative_residual: " << result.final_step.pressure_relative_residual << '\n';
    std::cout << "u_vertical_max: " << result.extrema.u_vertical_max << '\n';
    std::cout << "u_vertical_min: " << result.extrema.u_vertical_min << '\n';
    std::cout << "v_horizontal_max: " << result.extrema.v_horizontal_max << '\n';
    std::cout << "v_horizontal_min: " << result.extrema.v_horizontal_min << '\n';
    std::cout << "reference_dataset: " << result.validation.reference_dataset << '\n';
    std::cout << "max_reference_relative_error: " << result.validation.max_relative_error << '\n';
    for(const solver::LidDrivenCavityValidationPoint& point : result.validation.points) {
      if(point.label.empty()) {
        continue;
      }
      std::cout << "reference_point[" << point.label << "]: line=" << solver::to_string(point.line)
                << ", coordinate=" << point.coordinate << ", reference=" << point.reference_value
                << ", reference_lower=" << point.reference_lower_bound
                << ", reference_upper=" << point.reference_upper_bound
                << ", sample=" << point.sample_value
                << ", relative_error=" << point.relative_error << '\n';
    }
    if(!options.restart_path.empty()) {
      std::cout << "restart_path: " << options.restart_path << '\n';
    }
    if(!options.checkpoint_out.empty()) {
      std::cout << "checkpoint_out: " << options.checkpoint_out << '\n';
    }
    if(!options.vtk_out.empty()) {
      std::cout << "vtk_out: " << options.vtk_out << '\n';
    }

    if(partial_run || !config.validate_reference) {
      std::cout << "benchmark_status: skipped" << '\n';
      return 0;
    }

    if(result.validation.reference_dataset.empty()) {
      std::cout << "benchmark_status: unavailable" << '\n';
      return 1;
    }

    std::cout << "benchmark_status: " << (result.validation.pass ? "pass" : "fail") << '\n';
    return result.validation.pass ? 0 : 1;
  } catch(const std::exception& exception) {
    std::cerr << "solver_cavity failed: " << exception.what() << '\n';
    return 1;
  }
}
