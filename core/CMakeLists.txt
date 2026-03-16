add_library(solver_core STATIC
  runtime.cpp
)

target_include_directories(solver_core
  PUBLIC
    ${PROJECT_SOURCE_DIR}
)

target_link_libraries(solver_core
  PUBLIC
    solver_build_options
)

