#include "core/runtime.hpp"

#include <sstream>

namespace solver {

namespace {

bool is_supported_runtime_platform() {
#if defined(__APPLE__) && (defined(__aarch64__) || defined(__arm64__))
  return true;
#else
  return false;
#endif
}

}  // namespace

BuildInfo get_build_info() {
  return BuildInfo{
      .project_name = SOLVER_PROJECT_NAME,
      .project_version = SOLVER_PROJECT_VERSION,
      .compiler_id = SOLVER_COMPILER_ID,
      .compiler_version = SOLVER_COMPILER_VERSION,
      .build_profile = SOLVER_BUILD_PROFILE,
      .fast_math_enabled = SOLVER_FAST_MATH_ENABLED != 0,
      .supported_runtime_platform = is_supported_runtime_platform(),
  };
}

std::string format_build_banner(const BuildInfo& build_info) {
  std::ostringstream stream;
  stream << build_info.project_name << " v" << build_info.project_version << '\n'
         << "compiler: " << build_info.compiler_id << ' ' << build_info.compiler_version
         << '\n'
         << "profile: " << build_info.build_profile << '\n'
         << "fast_math: " << (build_info.fast_math_enabled ? "enabled" : "disabled") << '\n'
         << "runtime_platform: "
         << (build_info.supported_runtime_platform ? "apple-silicon" : "unsupported");
  return stream.str();
}

}  // namespace solver

