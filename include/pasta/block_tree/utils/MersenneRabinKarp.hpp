/*******************************************************************************
 * This file is part of pasta::block_tree
 *
 * Copyright (C) 2022 Daniel Meyer
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

#include <iostream>

namespace pasta {

template <class T, class size_type> class MersenneRabinKarp {
  __extension__ typedef unsigned __int128 uint128_t;

public:
  std::vector<T> const &text_;
  uint128_t sigma_;
  uint64_t init_;
  uint64_t length_;
  uint128_t prime_;
  uint64_t hash_;
  uint128_t max_sigma_;

  MersenneRabinKarp(std::vector<T> const &text, uint64_t sigma, uint64_t init,
                    uint64_t length, uint128_t prime)
      : text_(text), sigma_(sigma), init_(init), length_(length),
        prime_(prime) {
    max_sigma_ = 1;
    uint128_t fp = 0;
    uint128_t sigma_c = 1;
    for (uint64_t i = init_; i < init_ + length_; i++) {
      fp = mersenneModulo(fp * sigma);
      fp = mersenneModulo(fp + text_[i]);
    }
    for (uint64_t i = 0; i < length_ - 1; i++) {
      sigma_c = mersenneModulo(sigma_c * sigma_);
    }
    hash_ = (uint64_t)(fp);
    max_sigma_ = (uint64_t)(sigma_c);
  };

  void restart(uint64_t index) {
    if (index + length_ >= text_.size()) {
      return;
    }
    init_ = index;
    max_sigma_ = 1;
    uint128_t fp = 0;
    uint128_t sigma_c = 1;
    for (uint64_t i = init_; i < init_ + length_; i++) {
      fp = mersenneModulo(fp * sigma_);
      fp = mersenneModulo(fp + text_[i]);
    }
    for (uint64_t i = 0; i < length_ - 1; i++) {
      sigma_c = mersenneModulo(sigma_c * sigma_);
    }
    hash_ = (uint64_t)(fp);
    max_sigma_ = (uint64_t)(sigma_c);
  };

  inline uint128_t mersenneModulo(uint128_t k) {
    return k % prime_;
    //        uint128_t i = (k & prime_) + (k >> power_);
    //        return (i >= prime_) ? i - prime_ : i;
  };

  void next() {
    if (text_.size() <= init_ + length_) {
      return;
    }

    uint128_t fp = hash_;
    T out_char = text_[init_];
    T in_char = text_[init_ + length_];
    uint128_t out_char_influence = out_char * max_sigma_;
    out_char_influence = mersenneModulo(out_char_influence);
    if (out_char_influence < hash_) {
      fp -= out_char_influence;
    } else {
      fp = prime_ - (out_char_influence - fp);
    }
    fp *= sigma_;
    fp += in_char;
    fp = mersenneModulo(fp);
    hash_ = (uint64_t)(fp);
    init_++;
  };
};

} // namespace pasta

/******************************************************************************/
