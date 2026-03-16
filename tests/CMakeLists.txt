add_executable(solver_tests
  test_runtime.cpp
)

target_link_libraries(solver_tests
  PRIVATE
    solver_io
)

target_compile_definitions(solver_tests
  PRIVATE
    SOLVER_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)

add_test(NAME solver_tests COMMAND solver_tests)
