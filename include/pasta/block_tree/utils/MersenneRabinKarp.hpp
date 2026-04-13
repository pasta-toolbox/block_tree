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

#include "pasta/block_tree/utils/MersenneHash.hpp"

namespace pasta {

__extension__ typedef unsigned __int128 uint128_t;

template <uint8_t exponent>
static consteval uint128_t mersenne_prime() {
  uint128_t res = 1;
  for (size_t i = 0; i < exponent; i++) {
    res <<= 1;
  }
  return res - 1;
}

///
/// @brief A Rabin-Karp rolling hasher.
///
/// @tparam T The type of the characters in the text.
/// @tparam size_type The type to use for indexing etc.
/// @tparam mersenne_exponent If using a mersenne prime 2^p-1, then this should
///   be p. If this is 0, a normal modulus operation will be used.
///   Additionally, if this is != 0, then the prime_ attribute will be ignored
///   and '(1 << mersenne_exponent) - 1' will be used instead.
///
template <class T, class size_type, uint8_t mersenne_exponent = 0>
class MersenneRabinKarp {
public:
  /// The text being hashed
  std::span<const T> text_;
  uint128_t sigma_;
  /// The start index of the currently hashed window
  uint64_t init_;
  /// The window size of this hasher
  uint64_t length_;
  /// A large prime used for modulus operations
  uint128_t prime_;
  /// The current hash value
  uint128_t hash_;
  uint128_t max_sigma_;

  /// @brief Construct a new Rabin Karp hasher.
  /// @param text The text to hash. (not just the window but the entire text)
  /// @param sigma The alphabet size.
  /// @param init The start index of the first hashed window in the text.
  /// @param length The window size.
  /// @param prime A large prime used for modulus operations
  ///   iff not using mersenne_exponent.
  constexpr MersenneRabinKarp(const std::span<const T> text,
                              const uint64_t sigma,
                              const uint64_t init,
                              const uint64_t length,
                              const uint128_t prime)
      : text_(text),
        sigma_(sigma),
        init_(init),
        length_(length),
        prime_(prime) {
    max_sigma_ = 1;
    uint128_t fp = 0;
    uint128_t sigma_c = 1;
    for (uint64_t i = init_; i < init_ + length_; i++) {
      fp = fp * sigma;
      fp = mersenneModulo(fp + text_[i]);
    }
    for (uint64_t i = 0; i < length_ - 1; i++) {
      sigma_c = mersenneModulo(sigma_c * sigma_);
    }
    hash_ = fp;
    max_sigma_ = sigma_c;
  };

  constexpr MersenneRabinKarp(const std::vector<T>& text,
                              const uint64_t sigma,
                              const uint64_t init,
                              const uint64_t length,
                              const uint128_t prime)
      : MersenneRabinKarp(std::span(text), sigma, init, length, prime) {}

  /// @brief Moves the hasher to the specified start index in the backing
  ///   vector.
  void restart(const uint64_t index) {
    if (index + length_ >= text_.size()) {
      return;
    }
    init_ = index;
    uint128_t fp = 0;
    for (uint64_t i = init_; i < init_ + length_; i++) {
      fp = fp * sigma_;
      fp = mersenneModulo(fp + text_[i]);
    }
    hash_ = fp;
  };

  constexpr uint128_t mersenneModulo(uint128_t k) const {
    if constexpr (mersenne_exponent == 0) {
      return k % prime_;
    } else {
      constexpr uint128_t MERSENNE = mersenne_prime<mersenne_exponent>();
      uint128_t i = (k & MERSENNE) + (k >> mersenne_exponent);
      i -= (i >= MERSENNE) * MERSENNE;
      return i;
    }
  };

  /// @brief Retrieves the hash value at the hasher's current position.
  /// @return A MersenneHash object representing the current hash value.
  constexpr MersenneHash<T> current_hash() const {
    return MersenneHash<T>(text_, hash_, init_, length_);
  }

  /// @brief Advances the hasher by one character.
  constexpr void next() {
    if (text_.size() <= init_ + length_) {
      return;
    }

    uint128_t fp = hash_;
    T out_char = text_[init_];
    T in_char = text_[init_ + length_];
    const uint128_t out_char_influence = mersenneModulo(out_char * max_sigma_);
    // Conditionally add the prime, of the out_char_influence is too large
    if constexpr (mersenne_exponent == 0) {
      fp += prime_ * (out_char_influence > hash_) - out_char_influence;
    } else {
      fp += mersenne_prime<mersenne_exponent>() * (out_char_influence > hash_) -
            out_char_influence;
    }
    fp *= sigma_;
    fp += in_char;
    fp = mersenneModulo(fp);
    hash_ = fp;
    init_++;
  };
};

///
/// @brief A Rabin-Karp rolling hasher for bitstrings.
///
/// @tparam size_type The type to use for indexing etc.
/// @tparam mersenne_exponent If using a mersenne prime 2^p-1, then this should
///   be p. If this is 0, a normal modulus operation will be used.
///   Additionally, if this is != 0, then the prime_ attribute will be ignored
///   and '(1 << mersenne_exponent) - 1' will be used instead.
///
template <typename size_type, uint8_t mersenne_exponent>
class MersenneRabinKarp<bool, size_type, mersenne_exponent> {
public:
  /// The text being hashed
  const pasta::BitVector& text_;
  uint64_t init_;
  /// The window size of this hasher
  uint64_t length_;
  /// A large prime used for modulus operations
  uint128_t prime_;
  /// The current hash value
  uint128_t hash_;
  uint128_t max_sigma_;

  /// @brief Construct a new Rabin Karp hasher over a bitstring.
  /// @param text The bitvector to hash. (not just the window but the entire
  ///   text)
  /// @param init The start bit-index of the first hashed window in the text.
  /// @param length The window size in bits.
  /// @param prime A large prime used for modulus operations
  ///   iff not using mersenne_exponent.
  constexpr MersenneRabinKarp(const pasta::BitVector& text,
                              const uint64_t init,
                              const uint64_t length,
                              const uint128_t prime)
      : text_(text),
        init_(init),
        length_(length),
        prime_(prime) {
    max_sigma_ = 1;
    uint128_t fp = 0;
    uint128_t sigma_c = 1;
    for (uint64_t i = init_; i < init_ + length_; i++) {
      fp = mersenneModulo(2 * fp + get_bit(i));
    }
    for (uint64_t i = 0; i < length_ - 1; i++) {
      sigma_c = mersenneModulo(2 * sigma_c);
    }
    hash_ = fp;
    max_sigma_ = sigma_c;
  }

  /// @brief Moves the hasher to the specified start index in the backing
  ///   vector.
  constexpr void restart(const uint64_t index) {
    if (index + length_ >= text_.size()) {
      return;
    }
    init_ = index;
    uint128_t fp = 0;
    for (uint64_t i = init_; i < init_ + length_; i++) {
      fp = fp * 2;
      fp = mersenneModulo(fp + get_bit(i));
    }
    hash_ = fp;
  };

  [[nodiscard]] constexpr uint128_t mersenneModulo(const uint128_t k) const {
    if constexpr (mersenne_exponent == 0) {
      return k % prime_;
    } else {
      constexpr uint128_t MERSENNE = mersenne_prime<mersenne_exponent>();
      uint128_t i = (k & MERSENNE) + (k >> mersenne_exponent);
      i -= (i >= MERSENNE) * MERSENNE;
      return i;
    }
  };

  /// @brief Retrieves the hash value at the hasher's current position.
  /// @return A MersenneHash object representing the current hash value.
  [[nodiscard]] constexpr MersenneHash<bool> current_hash() const {
    return {text_, hash_, init_, length_};
  }

  /// @brief Advances the hasher by one character.
  constexpr void next() {
    if (text_.size() <= init_ + length_) {
      return;
    }

    uint128_t fp = hash_;
    const bool out_char = out_bit();
    const bool in_char = in_bit();
    const uint128_t out_char_influence = out_char * max_sigma_;
    // Conditionally add the prime, of the out_char_influence is too large
    if constexpr (mersenne_exponent == 0) {
      fp += prime_ * (out_char_influence > hash_) - out_char_influence;
    } else {
      fp += mersenne_prime<mersenne_exponent>() * (out_char_influence > hash_) -
            out_char_influence;
    }
    fp *= 2;
    fp += in_char;
    fp = mersenneModulo(fp);
    hash_ = fp;
    init_++;
  };

private:
  [[nodiscard]] constexpr bool out_bit() const {
    return get_bit(init_);
  }

  [[nodiscard]] constexpr bool in_bit() const {
    const size_t idx = init_ + length_;
    return get_bit(idx);
  }

  [[nodiscard]] constexpr bool get_bit(const size_t bit_index) const {
    return text_[bit_index];
  }
};

} // namespace pasta

/******************************************************************************/
