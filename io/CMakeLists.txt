add_library(solver_io STATIC
  checkpoint.cpp
  vtk_export.cpp
)

target_include_directories(solver_io
  PUBLIC
    ${PROJECT_SOURCE_DIR}
)

target_link_libraries(solver_io
  PUBLIC
    solver_momentum
    solver_core
)
