add_executable(solver_example
  example_main.cpp
)

target_link_libraries(solver_example
  PRIVATE
    solver_core
)

