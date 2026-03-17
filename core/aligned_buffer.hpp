#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <utility>

namespace solver {

template <typename T, std::size_t Alignment>
class AlignedBuffer {
 public:
  using value_type = T;
  static constexpr std::size_t alignment = Alignment;

  static_assert(Alignment >= alignof(T), "alignment must satisfy the value type");
  static_assert(Alignment % alignof(T) == 0, "alignment must be a multiple of alignof(T)");

  AlignedBuffer() = default;

  explicit AlignedBuffer(const std::size_t size) {
    resize(size);
  }

  AlignedBuffer(const AlignedBuffer& other)
      : AlignedBuffer(other.size_) {
    std::copy_n(other.data_, size_, data_);
  }

  AlignedBuffer(AlignedBuffer&& other) noexcept
      : data_(std::exchange(other.data_, nullptr)),
        size_(std::exchange(other.size_, 0)) {}

  AlignedBuffer& operator=(const AlignedBuffer& other) {
    if(this == &other) {
      return *this;
    }

    AlignedBuffer copy(other);
    swap(copy);
    return *this;
  }

  AlignedBuffer& operator=(AlignedBuffer&& other) noexcept {
    if(this == &other) {
      return *this;
    }

    reset();
    data_ = std::exchange(other.data_, nullptr);
    size_ = std::exchange(other.size_, 0);
    return *this;
  }

  ~AlignedBuffer() {
    reset();
  }

  void resize(const std::size_t size) {
    reset();
    if(size == 0) {
      return;
    }

    if(size > std::numeric_limits<std::size_t>::max() / sizeof(T)) {
      throw std::length_error("AlignedBuffer::resize: requested size overflows byte count");
    }

    data_ = static_cast<T*>(::operator new[](size * sizeof(T), std::align_val_t{Alignment}));
    size_ = size;
    std::fill_n(data_, size_, T{});
  }

  void reset() noexcept {
    if(data_ != nullptr) {
      ::operator delete[](data_, std::align_val_t{Alignment});
      data_ = nullptr;
    }
    size_ = 0;
  }

  void swap(AlignedBuffer& other) noexcept {
    std::swap(data_, other.data_);
    std::swap(size_, other.size_);
  }

  [[nodiscard]] T* data() noexcept {
    return data_;
  }

  [[nodiscard]] const T* data() const noexcept {
    return data_;
  }

  [[nodiscard]] std::size_t size() const noexcept {
    return size_;
  }

  [[nodiscard]] T& operator[](const std::size_t index) noexcept {
    return data_[index];
  }

  [[nodiscard]] const T& operator[](const std::size_t index) const noexcept {
    return data_[index];
  }

 private:
  T* data_ = nullptr;
  std::size_t size_ = 0;
};

}  // namespace solver

