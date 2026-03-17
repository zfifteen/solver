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

find_package(Python3 COMPONENTS Interpreter REQUIRED)

add_test(
  NAME validation_suite_tests
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_validation_suite.py
)
