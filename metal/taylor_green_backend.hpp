#pragma once

#include "solver/taylor_green.hpp"

#include <string>

namespace solver::metal {

struct TaylorGreenMetalRun {
  TaylorGreenState state;
  std::string device_name;
  double elapsed_seconds = 0.0;
};

[[nodiscard]] TaylorGreenMetalRun run_taylor_green(const TaylorGreenConfig& config);

}  // namespace solver::metal
