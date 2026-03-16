#pragma once

#include "solver/lid_driven_cavity.hpp"

#include <cstdint>
#include <string>

namespace solver::io {

struct CheckpointMetadata {
  std::uint32_t format_version = 0;
  std::string endianness = "little";
  std::uint32_t scalar_bytes = 0;
  std::uint64_t checksum = 0;
  std::uint64_t build_hash = 0;
  std::uint64_t configuration_hash = 0;
};

struct LidDrivenCavityCheckpoint {
  CheckpointMetadata metadata{};
  LidDrivenCavityState state;

  explicit LidDrivenCavityCheckpoint(const Grid& grid)
      : state(grid) {}
};

void write_lid_driven_cavity_checkpoint(const std::string& path,
                                        const LidDrivenCavityConfig& config,
                                        const LidDrivenCavityState& state);

[[nodiscard]] LidDrivenCavityCheckpoint load_lid_driven_cavity_checkpoint(
    const std::string& path,
    const LidDrivenCavityConfig& expected_config);

}  // namespace solver::io
