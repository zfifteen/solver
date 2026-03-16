#pragma once

#include "core/aligned_buffer.hpp"
#include "core/field_layout.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>

namespace solver {

class StructuredField {
 public:
  using value_type = double;
  static constexpr std::size_t storage_alignment = 64;

  explicit StructuredField(FieldLayout layout)
      : layout_(layout),
        data_(layout_.storage_size()) {}

  [[nodiscard]] const FieldLayout& layout() const noexcept {
    return layout_;
  }

  [[nodiscard]] double* data() noexcept {
    return data_.data();
  }

  [[nodiscard]] const double* data() const noexcept {
    return data_.data();
  }

  [[nodiscard]] std::size_t size() const noexcept {
    return data_.size();
  }

  [[nodiscard]] bool is_aligned() const noexcept {
    return reinterpret_cast<std::uintptr_t>(data()) % storage_alignment == 0;
  }

  void fill(const double value) {
    std::fill_n(data(), size(), value);
  }

  void fill_range(const IndexRange3D& range, const double value) {
    for(int k = range.k_begin; k < range.k_end; ++k) {
      for(int j = range.j_begin; j < range.j_end; ++j) {
        for(int i = range.i_begin; i < range.i_end; ++i) {
          (*this)(i, j, k) = value;
        }
      }
    }
  }

  void fill_ghost_layer(const BoundaryFace face, const int layer, const double value) {
    fill_range(layout_.ghost_range(face, layer), value);
  }

  [[nodiscard]] double& at(const int i, const int j, const int k) {
    return data_[layout_.index(i, j, k)];
  }

  [[nodiscard]] const double& at(const int i, const int j, const int k) const {
    return data_[layout_.index(i, j, k)];
  }

  [[nodiscard]] double& unchecked(const int i, const int j, const int k) noexcept {
    return data_[layout_.unchecked_index(i, j, k)];
  }

  [[nodiscard]] const double& unchecked(const int i, const int j, const int k) const noexcept {
    return data_[layout_.unchecked_index(i, j, k)];
  }

  [[nodiscard]] double* row_ptr(const int j, const int k) noexcept {
    return data() + layout_.unchecked_index(0, j, k);
  }

  [[nodiscard]] const double* row_ptr(const int j, const int k) const noexcept {
    return data() + layout_.unchecked_index(0, j, k);
  }

  [[nodiscard]] double& operator()(const int i, const int j, const int k) {
#if !defined(NDEBUG)
    assert(layout_.is_storage_index(i, j, k));
#endif
    return data_[layout_.unchecked_index(i, j, k)];
  }

  [[nodiscard]] const double& operator()(const int i, const int j, const int k) const {
#if !defined(NDEBUG)
    assert(layout_.is_storage_index(i, j, k));
#endif
    return data_[layout_.unchecked_index(i, j, k)];
  }

 protected:
  FieldLayout layout_;
  AlignedBuffer<double, storage_alignment> data_;
};

class ScalarField final : public StructuredField {
 public:
  explicit ScalarField(const Grid& grid)
      : StructuredField(FieldLayout::cell_centered(grid)) {}
};

class PressureField final : public StructuredField {
 public:
  explicit PressureField(const Grid& grid)
      : StructuredField(FieldLayout::cell_centered(grid)) {}
};

class FaceField final : public StructuredField {
 public:
  FaceField(const Axis normal_axis, const Grid& grid)
      : StructuredField(FieldLayout::face_centered(normal_axis, grid)),
        normal_axis_(normal_axis) {}

  [[nodiscard]] constexpr Axis normal_axis() const noexcept {
    return normal_axis_;
  }

 private:
  Axis normal_axis_;
};

struct VelocityField {
  FaceField x;
  FaceField y;
  FaceField z;

  explicit VelocityField(const Grid& grid)
      : x(Axis::x, grid),
        y(Axis::y, grid),
        z(Axis::z, grid) {}

  void fill(const double value) {
    x.fill(value);
    y.fill(value);
    z.fill(value);
  }
};

}  // namespace solver
