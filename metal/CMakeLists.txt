find_program(SOLVER_XCRUN_EXECUTABLE xcrun REQUIRED)
find_library(SOLVER_FOUNDATION_FRAMEWORK Foundation REQUIRED)
find_library(SOLVER_METAL_FRAMEWORK Metal REQUIRED)

set(SOLVER_METAL_USE_FLOAT 1)

set(SOLVER_TAYLOR_GREEN_METAL_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/taylor_green_backend.metal)
set(SOLVER_TAYLOR_GREEN_AIR ${CMAKE_CURRENT_BINARY_DIR}/taylor_green_backend.air)
set(SOLVER_TAYLOR_GREEN_METALLIB ${CMAKE_CURRENT_BINARY_DIR}/taylor_green_backend.metallib)

add_custom_command(
  OUTPUT ${SOLVER_TAYLOR_GREEN_AIR}
  COMMAND ${SOLVER_XCRUN_EXECUTABLE} metal
          -c
          -DSOLVER_METAL_USE_FLOAT=${SOLVER_METAL_USE_FLOAT}
          ${SOLVER_TAYLOR_GREEN_METAL_SOURCE}
          -o ${SOLVER_TAYLOR_GREEN_AIR}
  DEPENDS ${SOLVER_TAYLOR_GREEN_METAL_SOURCE}
  VERBATIM
)

add_custom_command(
  OUTPUT ${SOLVER_TAYLOR_GREEN_METALLIB}
  COMMAND ${SOLVER_XCRUN_EXECUTABLE} metallib
          ${SOLVER_TAYLOR_GREEN_AIR}
          -o ${SOLVER_TAYLOR_GREEN_METALLIB}
  DEPENDS ${SOLVER_TAYLOR_GREEN_AIR}
  VERBATIM
)

add_custom_target(solver_metal_shaders ALL
  DEPENDS ${SOLVER_TAYLOR_GREEN_METALLIB}
)

add_library(solver_metal STATIC
  taylor_green_backend.mm
)

add_dependencies(solver_metal solver_metal_shaders)

target_include_directories(solver_metal
  PUBLIC
    ${PROJECT_SOURCE_DIR}
)

target_compile_definitions(solver_metal
  PRIVATE
    SOLVER_METAL_LIBRARY_PATH="${SOLVER_TAYLOR_GREEN_METALLIB}"
    SOLVER_METAL_USE_FLOAT=${SOLVER_METAL_USE_FLOAT}
)

target_compile_options(solver_metal
  PRIVATE
    -fobjc-arc
)

target_link_libraries(solver_metal
  PUBLIC
    solver_core
    ${SOLVER_FOUNDATION_FRAMEWORK}
    ${SOLVER_METAL_FRAMEWORK}
)
