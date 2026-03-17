add_library(solver_momentum STATIC
  channel_flow.cpp
  lid_driven_cavity.cpp
  momentum_terms.cpp
  operator_verification.cpp
  projection.cpp
  taylor_green.cpp
)

target_include_directories(solver_momentum
  PUBLIC
    ${PROJECT_SOURCE_DIR}
)

target_link_libraries(solver_momentum
  PUBLIC
    solver_operators
    solver_linsolve
    solver_metal
)
