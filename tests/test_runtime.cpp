#include "core/runtime.hpp"

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}

void test_build_profile_is_locked() {
  const solver::BuildInfo build_info = solver::get_build_info();
  require(build_info.build_profile == "deterministic" || build_info.build_profile == "benchmark",
          "unexpected build profile");
}

void test_runtime_platform_is_supported() {
  const solver::BuildInfo build_info = solver::get_build_info();
  require(build_info.supported_runtime_platform, "expected Apple Silicon runtime");
}

void test_banner_contains_profile() {
  const solver::BuildInfo build_info = solver::get_build_info();
  const std::string banner = solver::format_build_banner(build_info);
  require(banner.find("profile: " + build_info.build_profile) != std::string::npos,
          "banner missing build profile");
}

}  // namespace

int main() {
  try {
    test_build_profile_is_locked();
    test_runtime_platform_is_supported();
    test_banner_contains_profile();
  } catch(const std::exception& exception) {
    std::cerr << "solver_tests failed: " << exception.what() << '\n';
    return 1;
  }

  std::cout << "solver_tests passed" << '\n';
  return 0;
}

