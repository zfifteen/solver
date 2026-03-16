#include "solver/lid_driven_cavity.hpp"

#include <iostream>
#include <string>

int main(int argc, char** argv) {
  try {
    const std::string config_path =
        argc > 1 ? argv[1] : "benchmarks/lid_driven_cavity_re100_128.cfg";
    const solver::LidDrivenCavityConfig config = solver::load_lid_driven_cavity_config(config_path);
    const solver::LidDrivenCavityResult result = solver::run_lid_driven_cavity(config);

    std::cout << "config_path: " << config_path << '\n';
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

    if(!config.validate_reference) {
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
