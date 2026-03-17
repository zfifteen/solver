#pragma once

#include <stdexcept>

namespace solver {

enum class Axis : int {
  x = 0,
  y = 1,
  z = 2,
};

enum class BoundaryFace : int {
  x_min = 0,
  x_max = 1,
  y_min = 2,
  y_max = 3,
  z_min = 4,
  z_max = 5,
};

inline constexpr int axis_index(const Axis axis) noexcept {
  return static_cast<int>(axis);
}

inline constexpr Axis boundary_axis(const BoundaryFace face) noexcept {
  switch(face) {
    case BoundaryFace::x_min:
    case BoundaryFace::x_max:
      return Axis::x;
    case BoundaryFace::y_min:
    case BoundaryFace::y_max:
      return Axis::y;
    case BoundaryFace::z_min:
    case BoundaryFace::z_max:
      return Axis::z;
  }

  __builtin_unreachable();
}

inline constexpr bool is_lower_boundary(const BoundaryFace face) noexcept {
  switch(face) {
    case BoundaryFace::x_min:
    case BoundaryFace::y_min:
    case BoundaryFace::z_min:
      return true;
    case BoundaryFace::x_max:
    case BoundaryFace::y_max:
    case BoundaryFace::z_max:
      return false;
  }

  __builtin_unreachable();
}

struct Grid {
  int nx;
  int ny;
  int nz;
  double dx;
  double dy;
  double dz;
  int ghost_layers;

  Grid(const int nx_in,
       const int ny_in,
       const int nz_in,
       const double dx_in,
       const double dy_in,
       const double dz_in,
       const int ghost_layers_in = 1)
      : nx(nx_in),
        ny(ny_in),
        nz(nz_in),
        dx(dx_in),
        dy(dy_in),
        dz(dz_in),
        ghost_layers(ghost_layers_in) {
    if(nx <= 0 || ny <= 0 || nz <= 0) {
      throw std::invalid_argument("grid dimensions must be positive");
    }
    if(dx <= 0.0 || dy <= 0.0 || dz <= 0.0) {
      throw std::invalid_argument("grid spacing must be positive");
    }
    if(ghost_layers <= 0) {
      throw std::invalid_argument("grid must define at least one ghost layer");
    }
  }

  [[nodiscard]] constexpr bool is_2d() const noexcept {
    return nz == 1;
  }

  [[nodiscard]] constexpr int cells(const Axis axis) const noexcept {
    switch(axis) {
      case Axis::x:
        return nx;
      case Axis::y:
        return ny;
      case Axis::z:
        return nz;
    }

    __builtin_unreachable();
  }

  [[nodiscard]] constexpr double spacing(const Axis axis) const noexcept {
    switch(axis) {
      case Axis::x:
        return dx;
      case Axis::y:
        return dy;
      case Axis::z:
        return dz;
    }

    __builtin_unreachable();
  }

  [[nodiscard]] double cell_center(const Axis axis, const int cell_index) const {
    if(cell_index < 0 || cell_index >= cells(axis)) {
      throw std::out_of_range("cell index outside grid extent");
    }

    return (static_cast<double>(cell_index) + 0.5) * spacing(axis);
  }

  [[nodiscard]] double face_coordinate(const Axis axis, const int face_index) const {
    if(face_index < 0 || face_index > cells(axis)) {
      throw std::out_of_range("face index outside grid extent");
    }

    return static_cast<double>(face_index) * spacing(axis);
  }
};

}  // namespace solver

