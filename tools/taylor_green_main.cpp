#include "io/vtk_export.hpp"
#include "solver/taylor_green.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
  try {
    std::string config_path = "benchmarks/taylor_green_128.cfg";
    std::string vtk_out;
    std::string backend_override;

    bool config_consumed = false;
    for(int index = 1; index < argc; ++index) {
      const std::string argument = argv[index];
      if(argument == "--vtk-out") {
        if(index + 1 >= argc) {
          throw std::runtime_error("missing value for --vtk-out");
        }
        vtk_out = argv[++index];
      } else if(argument == "--backend") {
        if(index + 1 >= argc) {
          throw std::runtime_error("missing value for --backend");
        }
        backend_override = argv[++index];
      } else if(!config_consumed && !argument.empty() && argument.rfind("--", 0) != 0) {
        config_path = argument;
        config_consumed = true;
      } else {
        throw std::runtime_error("unsupported argument: " + argument);
      }
    }

    solver::TaylorGreenConfig config = solver::load_taylor_green_config(config_path);
    if(!backend_override.empty()) {
      config.backend = solver::parse_execution_backend(backend_override);
    }
    solver::TaylorGreenState state = solver::initialize_taylor_green_state(config);
    const solver::TaylorGreenResult result = solver::run_taylor_green(config, &state);
    if(!vtk_out.empty()) {
      solver::io::write_mac_fields_vtk(vtk_out, state.velocity, state.pressure_total);
    }

    std::cout << "config_path: " << config_path << '\n';
    std::cout << "config: " << solver::describe(config) << '\n';
    std::cout << "backend_used: " << solver::to_string(result.backend_used) << '\n';
    if(!result.accelerator_name.empty()) {
      std::cout << "accelerator_name: " << result.accelerator_name << '\n';
    }
    if(result.backend_elapsed_seconds > 0.0) {
      std::cout << "backend_elapsed_seconds: " << result.backend_elapsed_seconds << '\n';
    }
    if(result.cleanup_elapsed_seconds > 0.0) {
      std::cout << "cleanup_elapsed_seconds: " << result.cleanup_elapsed_seconds << '\n';
    }
    std::cout << "steps: " << result.final_step.step << '\n';
    std::cout << "time: " << result.final_step.time << '\n';
    std::cout << "dt: " << result.final_step.dt << '\n';
    std::cout << "max_cfl: " << result.final_step.max_cfl << '\n';
    std::cout << "max_velocity_change: " << result.final_step.max_velocity_change << '\n';
    std::cout << "divergence_l2: " << result.final_step.divergence_l2 << '\n';
    std::cout << "max_divergence_l2: " << result.final_step.max_divergence_l2 << '\n';
    std::cout << "pressure_iterations: " << result.final_step.pressure_iterations << '\n';
    std::cout << "pressure_relative_residual: " << result.final_step.pressure_relative_residual
              << '\n';
    std::cout << "initial_kinetic_energy: " << result.initial_kinetic_energy << '\n';
    std::cout << "final_kinetic_energy: " << result.final_kinetic_energy << '\n';
    std::cout << "reference_dataset: " << result.validation.reference_dataset << '\n';
    std::cout << "exact_kinetic_energy: " << result.validation.exact_kinetic_energy << '\n';
    std::cout << "normalized_energy_error: " << result.validation.normalized_energy_error << '\n';
    std::cout << "velocity_relative_l2_error: " << result.validation.velocity_relative_l2_error
              << '\n';
    if(!vtk_out.empty()) {
      std::cout << "vtk_out: " << vtk_out << '\n';
    }

    if(!config.validate_energy) {
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
    std::cerr << "solver_taylor_green failed: " << exception.what() << '\n';
    return 1;
  }
}
