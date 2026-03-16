#pragma once

#include <string>

namespace solver {

struct BuildInfo {
  std::string project_name;
  std::string project_version;
  std::string compiler_id;
  std::string compiler_version;
  std::string build_profile;
  bool fast_math_enabled;
  bool supported_runtime_platform;
};

BuildInfo get_build_info();
std::string format_build_banner(const BuildInfo& build_info);

}  // namespace solver

