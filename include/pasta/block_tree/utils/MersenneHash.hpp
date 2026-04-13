/*******************************************************************************
 * This file is part of pasta::block_tree
 *
 * Copyright (C) 2022 Daniel Meyer
 * Copyright (C) 2023 Etienne Palanga
 *
 * pasta::block_tree is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * pasta::block_tree is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with pasta::block_tree.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <pasta/bit_vector/bit_vector.hpp>
#include <span>
#include <vector>

#ifdef BT_INSTRUMENT
#  include <iostream>
#endif

namespace pasta {

#ifdef BT_INSTRUMENT
static std::atomic_size_t mersenne_hash_comparisons = 0;
static std::atomic_size_t mersenne_hash_equals = 0;
static std::atomic_size_t mersenne_hash_collisions = 0;

void print_hash_data() {
  std::cerr << "comparisons: " << mersenne_hash_comparisons
            << ", equals: " << mersenne_hash_equals
            << ", collisions: " << mersenne_hash_collisions
            << ", percent equals: "
            << 100 * mersenne_hash_equals / ((double)mersenne_hash_comparisons)
            << ", percent collisions: "
            << 100 * mersenne_hash_collisions /
                   ((double)mersenne_hash_comparisons)
            << std::endl;
}
#endif

template <typename T>
class MersenneHash {
public:
  __extension__ typedef unsigned __int128 uint128_t;
  std::span<const T> text_;
  uint128_t hash_;
  uint32_t start_;
  uint32_t length_;
  MersenneHash(std::vector<T> const& text,
               const uint128_t hash,
               const uint64_t start,
               const uint64_t length)
      : text_(text),
        hash_(hash),
        start_(start),
        length_(length){};

  MersenneHash(const std::span<const T> text,
               const uint128_t hash,
               const uint64_t start,
               const uint64_t length)
      : text_(text),
        hash_(hash),
        start_(start),
        length_(length){};

  constexpr MersenneHash() : text_(), hash_(0), start_(0), length_(0){};

  constexpr MersenneHash(const MersenneHash& other) = default;
  constexpr MersenneHash(MersenneHash&& other) = default;

  MersenneHash& operator=(const MersenneHash& other) = default;
  MersenneHash& operator=(MersenneHash&& other) = default;

  bool operator==(const MersenneHash& other) const {
#ifdef BT_INSTRUMENT
    ++mersenne_hash_comparisons;
#endif
    // if (length_ != other.length_)
    // return false;
    // std::cout << static_cast<uint64_t>(hash_) << ", "
    //          << static_cast<uint64_t>(other.hash_) << std::endl;
    if (hash_ != other.hash_)
      return false;

    const bool is_same = memcmp(text_.data() + start_,
                                other.text_.data() + other.start_,
                                length_) == 0;

#ifdef BT_INSTRUMENT
    if (!is_same) {
      // The hash is the same but the substring isn't => collision
      ++mersenne_hash_collisions;
    } else {
      // The substrings are the same
      ++mersenne_hash_equals;
    }
#endif
    return is_same;
  };
};

template <>
class MersenneHash<bool> {
public:
  __extension__ typedef unsigned __int128 uint128_t;
  /// @brief The whole string in which the substring lies
  const pasta::BitVector* text_;
  uint128_t hash_;
  /// @brief The bit start-position of the hashed substring
  uint64_t start_;
  /// @brief The number of bits in the hashed substring
  uint64_t length_;

  MersenneHash(const pasta::BitVector& text,
               const uint128_t hash,
               const uint64_t start,
               const uint64_t length)
      : text_{&text},
        hash_{hash},
        start_{start},
        length_{length} {};

  constexpr MersenneHash() : text_(nullptr), hash_(0), start_(0), length_(0){};

  constexpr MersenneHash(const MersenneHash& other) = default;
  constexpr MersenneHash(MersenneHash&& other) = default;

  constexpr MersenneHash& operator=(const MersenneHash& other) = default;
  constexpr MersenneHash& operator=(MersenneHash&& other) = default;

  inline bool operator==(const MersenneHash& other) const {
#ifdef BT_INSTRUMENT
    ++mersenne_hash_comparisons;
#endif
    if (hash_ != other.hash_) {
      //std::cout << "hashes unequal" << std::endl;
      return false;
    }

    size_t pos = 0;
    bool is_same = true;
    for (size_t remaining = length_; remaining > 64;
         pos += 64, remaining -= 64) {
      if (slice_at(*text_, start_ + pos) !=
          slice_at(*other.text_, other.start_ + pos)) {
        is_same = false;
        break;
      }
    }

    if (!is_same) {
      for (size_t i = pos; i < length_; ++i) {
        if ((*text_)[start_ + i] != (*other.text_)[other.start_ + i]) {
          is_same = false;
          break;
        }
      }
    }

#ifdef BT_INSTRUMENT
    if (!is_same) {
      // The hash is the same but the substring isn't => collision
      ++mersenne_hash_collisions;
    } else {
      // The substrings are the same
      ++mersenne_hash_equals;
    }
#endif
    return is_same;
  };

  static uint64_t slice_at(const pasta::BitVector& bv, const size_t i) {
    const std::span<uint64_t> backing = bv.data();
    const uint8_t offset = i % 64;
    const size_t data_index = i / 64;
    // TODO Check if the right shift actually shifts in zeros
    const uint64_t r =
        backing[data_index] & (~static_cast<uint64_t>(0) << offset);
    if (offset > 0) {
      const uint64_t l = backing[data_index + 1] &
                         (~static_cast<uint64_t>(0) >> (64 - offset));
      return (l << (64 - offset)) | (r >> offset);
    }
    return r;
  }
};

} // namespace pasta

template <typename T>
struct std::hash<pasta::MersenneHash<T>> {
  typename pasta::MersenneHash<T>::uint128_t
  operator()(const pasta::MersenneHash<T>& hS) const {
    return hS.hash_;
  }
}; // namespace std

/******************************************************************************/
