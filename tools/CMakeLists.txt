add_executable(solver_example
  example_main.cpp
)

target_link_libraries(solver_example
  PRIVATE
    solver_momentum
)

add_executable(solver_cavity
  lid_driven_cavity_main.cpp
)

target_link_libraries(solver_cavity
  PRIVATE
    solver_io
)

add_executable(solver_channel
  channel_flow_main.cpp
)

target_link_libraries(solver_channel
  PRIVATE
    solver_momentum
)
