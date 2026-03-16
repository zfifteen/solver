add_executable(solver_tests
  test_runtime.cpp
)

target_link_libraries(solver_tests
  PRIVATE
    solver_core
)

add_test(NAME solver_tests COMMAND solver_tests)

