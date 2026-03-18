cmake_minimum_required(VERSION 3.24)

set(ORCHARD_FLOW_DISPLAY_NAME "Orchard Flow")

project(orchard_flow VERSION 0.1.0 LANGUAGES CXX OBJCXX)

include(CheckCXXCompilerFlag)
include(CTest)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT APPLE)
  message(FATAL_ERROR "solver supports macOS on Apple Silicon only.")
endif()

if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64)$")
  message(FATAL_ERROR "solver supports Apple Silicon (arm64) only.")
endif()

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  message(FATAL_ERROR "solver requires Apple Clang for the locked Milestone 0 environment.")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

set(SOLVER_BUILD_PROFILE "deterministic" CACHE STRING
    "Locked build profile: deterministic or benchmark")
set_property(CACHE SOLVER_BUILD_PROFILE PROPERTY STRINGS deterministic benchmark)

if(NOT SOLVER_BUILD_PROFILE STREQUAL "deterministic" AND
   NOT SOLVER_BUILD_PROFILE STREQUAL "benchmark")
  message(FATAL_ERROR "SOLVER_BUILD_PROFILE must be 'deterministic' or 'benchmark'.")
endif()

set(SOLVER_FAST_MATH_ENABLED 0)
if(SOLVER_BUILD_PROFILE STREQUAL "benchmark")
  set(SOLVER_FAST_MATH_ENABLED 1)
endif()

check_cxx_compiler_flag("-march=armv8-a" SOLVER_HAS_MARCH_ARMV8_A)
if(SOLVER_FAST_MATH_ENABLED)
  check_cxx_compiler_flag("-ffast-math" SOLVER_HAS_FAST_MATH)
else()
  set(SOLVER_HAS_FAST_MATH OFF)
endif()

add_library(solver_build_options INTERFACE)
target_compile_features(solver_build_options INTERFACE cxx_std_20)
target_compile_definitions(solver_build_options INTERFACE
  SOLVER_PROJECT_NAME="${ORCHARD_FLOW_DISPLAY_NAME}"
  SOLVER_PROJECT_VERSION="${PROJECT_VERSION}"
  SOLVER_COMPILER_ID="${CMAKE_CXX_COMPILER_ID}"
  SOLVER_COMPILER_VERSION="${CMAKE_CXX_COMPILER_VERSION}"
  SOLVER_BUILD_PROFILE="${SOLVER_BUILD_PROFILE}"
  SOLVER_FAST_MATH_ENABLED=${SOLVER_FAST_MATH_ENABLED})
target_compile_options(solver_build_options INTERFACE
  -Wall
  -Wextra
  -Wpedantic)

if(SOLVER_HAS_MARCH_ARMV8_A)
  target_compile_options(solver_build_options INTERFACE -march=armv8-a)
endif()

if(SOLVER_FAST_MATH_ENABLED)
  if(SOLVER_HAS_FAST_MATH)
    target_compile_options(solver_build_options INTERFACE -ffast-math)
  else()
    message(FATAL_ERROR "Benchmark profile requested -ffast-math, but the compiler rejected it.")
  endif()
endif()

message(STATUS "Configuring ${ORCHARD_FLOW_DISPLAY_NAME} ${PROJECT_VERSION}")
message(STATUS "  Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Locked profile: ${SOLVER_BUILD_PROFILE}")
message(STATUS "  -march=armv8-a: ${SOLVER_HAS_MARCH_ARMV8_A}")
message(STATUS "  -ffast-math: ${SOLVER_HAS_FAST_MATH}")

add_subdirectory(core)
add_subdirectory(operators)
add_subdirectory(linsolve)
add_subdirectory(metal)
add_subdirectory(solver)
add_subdirectory(io)
add_subdirectory(tools)
add_subdirectory(tests)
