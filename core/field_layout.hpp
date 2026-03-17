#pragma once

#include "core/grid.hpp"

#include <cstddef>
#include <stdexcept>

namespace solver {

enum class FieldLocation : int {
  cell_center = 0,
  face_x = 1,
  face_y = 2,
  face_z = 3,
};

struct Extent3D {
  int nx;
  int ny;
  int nz;

  [[nodiscard]] constexpr int count(const Axis axis) const noexcept {
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

  [[nodiscard]] constexpr std::size_t cell_count() const noexcept {
    return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
           static_cast<std::size_t>(nz);
  }
};

struct Index3D {
  int i;
  int j;
  int k;
};

struct IndexRange3D {
  int i_begin;
  int i_end;
  int j_begin;
  int j_end;
  int k_begin;
  int k_end;

  [[nodiscard]] constexpr Extent3D extent() const noexcept {
    return Extent3D{
        .nx = i_end - i_begin,
        .ny = j_end - j_begin,
        .nz = k_end - k_begin,
    };
  }

  [[nodiscard]] constexpr bool contains(const int i, const int j, const int k) const noexcept {
    return i >= i_begin && i < i_end && j >= j_begin && j < j_end && k >= k_begin &&
           k < k_end;
  }
};

class FieldLayout {
 public:
  [[nodiscard]] static FieldLayout cell_centered(const Grid& grid) {
    return FieldLayout(FieldLocation::cell_center, grid, Extent3D{grid.nx, grid.ny, grid.nz});
  }

  [[nodiscard]] static FieldLayout face_centered(const Axis axis, const Grid& grid) {
    switch(axis) {
      case Axis::x:
        return FieldLayout(FieldLocation::face_x, grid, Extent3D{grid.nx + 1, grid.ny, grid.nz});
      case Axis::y:
        return FieldLayout(FieldLocation::face_y, grid, Extent3D{grid.nx, grid.ny + 1, grid.nz});
      case Axis::z:
        return FieldLayout(FieldLocation::face_z, grid, Extent3D{grid.nx, grid.ny, grid.nz + 1});
    }

    __builtin_unreachable();
  }

  [[nodiscard]] constexpr FieldLocation location() const noexcept {
    return location_;
  }

  [[nodiscard]] constexpr const Grid& grid() const noexcept {
    return grid_;
  }

  [[nodiscard]] constexpr int ghost_layers() const noexcept {
    return ghost_layers_;
  }

  [[nodiscard]] constexpr Extent3D active_extent() const noexcept {
    return active_;
  }

  [[nodiscard]] constexpr Extent3D storage_extent() const noexcept {
    return storage_;
  }

  [[nodiscard]] constexpr std::size_t storage_size() const noexcept {
    return storage_.cell_count();
  }

  [[nodiscard]] constexpr std::size_t storage_stride_j() const noexcept {
    return static_cast<std::size_t>(storage_.nx);
  }

  [[nodiscard]] constexpr std::size_t storage_stride_k() const noexcept {
    return static_cast<std::size_t>(storage_.nx) * static_cast<std::size_t>(storage_.ny);
  }

  [[nodiscard]] constexpr bool same_shape_as(const FieldLayout& other) const noexcept {
    return location_ == other.location_ && ghost_layers_ == other.ghost_layers_ &&
           active_.nx == other.active_.nx && active_.ny == other.active_.ny &&
           active_.nz == other.active_.nz && storage_.nx == other.storage_.nx &&
           storage_.ny == other.storage_.ny && storage_.nz == other.storage_.nz;
  }

  [[nodiscard]] constexpr bool is_unit_stride_i() const noexcept {
    return true;
  }

  [[nodiscard]] constexpr IndexRange3D active_range() const noexcept {
    return IndexRange3D{
        .i_begin = ghost_layers_,
        .i_end = ghost_layers_ + active_.nx,
        .j_begin = ghost_layers_,
        .j_end = ghost_layers_ + active_.ny,
        .k_begin = ghost_layers_,
        .k_end = ghost_layers_ + active_.nz,
    };
  }

  [[nodiscard]] constexpr bool is_storage_index(const int i, const int j, const int k) const noexcept {
    return i >= 0 && i < storage_.nx && j >= 0 && j < storage_.ny && k >= 0 && k < storage_.nz;
  }

  [[nodiscard]] constexpr bool is_active_storage_index(const int i,
                                                       const int j,
                                                       const int k) const noexcept {
    return active_range().contains(i, j, k);
  }

  [[nodiscard]] Index3D storage_index_from_active(const int i, const int j, const int k) const {
    if(i < 0 || i >= active_.nx || j < 0 || j >= active_.ny || k < 0 || k >= active_.nz) {
      throw std::out_of_range("active index outside field extent");
    }

    return unchecked_storage_index_from_active(i, j, k);
  }

  [[nodiscard]] constexpr Index3D unchecked_storage_index_from_active(const int i,
                                                                      const int j,
                                                                      const int k) const noexcept {
    return Index3D{
        .i = i + ghost_layers_,
        .j = j + ghost_layers_,
        .k = k + ghost_layers_,
    };
  }

  [[nodiscard]] std::size_t index(const int i, const int j, const int k) const {
    if(!is_storage_index(i, j, k)) {
      throw std::out_of_range("storage index outside field extent");
    }

    return unchecked_index(i, j, k);
  }

  [[nodiscard]] constexpr std::size_t unchecked_index(const int i,
                                                      const int j,
                                                      const int k) const noexcept {
    return static_cast<std::size_t>(i) + storage_stride_j() * static_cast<std::size_t>(j) +
           storage_stride_k() * static_cast<std::size_t>(k);
  }

  [[nodiscard]] IndexRange3D boundary_active_range(const BoundaryFace face) const {
    IndexRange3D range = active_range();
    const Axis axis = boundary_axis(face);
    const int first = is_lower_boundary(face) ? ghost_layers_
                                              : ghost_layers_ + active_.count(axis) - 1;
    const int last = first + 1;

    set_axis_range(axis, first, last, range);
    return range;
  }

  [[nodiscard]] IndexRange3D ghost_range(const BoundaryFace face, const int layer = 0) const {
    if(layer < 0 || layer >= ghost_layers_) {
      throw std::out_of_range("ghost layer outside configured range");
    }

    IndexRange3D range = active_range();
    const Axis axis = boundary_axis(face);

    if(is_lower_boundary(face)) {
      const int begin = ghost_layers_ - 1 - layer;
      set_axis_range(axis, begin, begin + 1, range);
    } else {
      const int begin = ghost_layers_ + active_.count(axis) + layer;
      set_axis_range(axis, begin, begin + 1, range);
    }

    return range;
  }

  [[nodiscard]] double coordinate_at_active_index(const Axis axis, const int active_index) const {
    if(active_index < 0 || active_index >= active_.count(axis)) {
      throw std::out_of_range("active coordinate index outside field extent");
    }

    if(location_ != FieldLocation::cell_center && axis == normal_axis(location_)) {
      return grid_.face_coordinate(axis, active_index);
    }

    return grid_.cell_center(axis, active_index);
  }

  [[nodiscard]] double coordinate_for_storage_index(const Axis axis, const int storage_index) const {
    if(storage_index < 0 || storage_index >= storage_.count(axis)) {
      throw std::out_of_range("storage coordinate index outside field extent");
    }

    const int shifted_index = storage_index - ghost_layers_;
    if(location_ != FieldLocation::cell_center && axis == normal_axis(location_)) {
      return static_cast<double>(shifted_index) * grid_.spacing(axis);
    }

    return (static_cast<double>(shifted_index) + 0.5) * grid_.spacing(axis);
  }

 private:
  FieldLayout(const FieldLocation location, const Grid& grid, const Extent3D active)
      : active_(active),
        storage_(
            Extent3D{active.nx + 2 * grid.ghost_layers, active.ny + 2 * grid.ghost_layers,
                     active.nz + 2 * grid.ghost_layers}),
        grid_(grid),
        location_(location),
        ghost_layers_(grid.ghost_layers) {}

  [[nodiscard]] static constexpr Axis normal_axis(const FieldLocation location) noexcept {
    switch(location) {
      case FieldLocation::cell_center:
        return Axis::x;
      case FieldLocation::face_x:
        return Axis::x;
      case FieldLocation::face_y:
        return Axis::y;
      case FieldLocation::face_z:
        return Axis::z;
    }

    __builtin_unreachable();
  }

  static constexpr void set_axis_range(const Axis axis,
                                       const int begin,
                                       const int end,
                                       IndexRange3D& range) noexcept {
    switch(axis) {
      case Axis::x:
        range.i_begin = begin;
        range.i_end = end;
        return;
      case Axis::y:
        range.j_begin = begin;
        range.j_end = end;
        return;
      case Axis::z:
        range.k_begin = begin;
        range.k_end = end;
        return;
    }
  }

  Extent3D active_;
  Extent3D storage_;
  Grid grid_;
  FieldLocation location_;
  int ghost_layers_;
};

}  // namespace solver
