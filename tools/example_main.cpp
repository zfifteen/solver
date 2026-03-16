#include "core/runtime.hpp"

#include <iostream>

int main() {
  const solver::BuildInfo build_info = solver::get_build_info();

  std::cout << solver::format_build_banner(build_info) << '\n';
  std::cout << "example_status: ok" << '\n';

  return build_info.supported_runtime_platform ? 0 : 1;
}

