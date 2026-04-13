/*******************************************************************************
 * This file is part of pasta::block_tree
 *
 * Copyright (C) 2022 Daniel Meyer
 * Copyright (C) 2023 Etienne Palanga <etienne.palanga@tu-dortmund.de>
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

#include <bit>
#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/optimized_for.hpp>
#include <pasta/bit_vector/support/rank_select.hpp>
#include <sdsl/int_vector.hpp>
#include <thread>
#include <vector>

namespace pasta {

template <std::signed_integral size_type, uint8_t recursion_level = 0>
class RecursiveBitBlockTree {
public:
  constexpr static bool types_is_block_tree = recursion_level > 0;
  using IsInternalType =
      std::conditional_t<types_is_block_tree,
                         RecursiveBitBlockTree<size_type, recursion_level - 1>,
                         pasta::BitVector>;
  using IsInternalRankType =
      std::conditional_t<types_is_block_tree,
                         RecursiveBitBlockTree<size_type, recursion_level - 1>,
                         pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>>;

  /// If this is true, then the only levels of the tree start to be
  ///   included starting at the first level that contains a back block
  ///
  /// For example, if levels 0 to 5 do not contain any back blocks, then the
  /// tree will only contain levels 6 and below.
  bool CUT_FIRST_LEVELS = true;

  /// The arity of the tree
  size_type tau_;
  size_type max_leaf_length_;
  /// The arity of the tree's root
  size_type s_ = 1;
  size_type leaf_size = 0;
  size_type amount_of_leaves = 0;
  size_type num_bits_;
  bool rank_support = false;
  /// Recursively compress the bit vectors of the tree
  std::vector<IsInternalType*> block_tree_types_;
  std::vector<IsInternalRankType*> block_tree_types_rs_;
  /// For each level and each back block, contains the index of the
  /// block's source
  std::vector<sdsl::int_vector<>*> block_tree_pointers_;
  std::vector<sdsl::int_vector<>*> block_tree_offsets_;
  //    std::vector<sdsl::int_vector<>*> block_tree_encoded_;
  std::vector<int64_t> block_size_lvl_;
  std::vector<int64_t> block_per_lvl_;
  std::vector<uint8_t> leaves_;

  std::vector<uint8_t> compress_map_;
  std::vector<uint8_t> decompress_map_;
  sdsl::int_vector<> compressed_leaves_;

  /// @brief For each level and each block, contains the number of 1s up to (and
  /// including) the block.
  std::vector<sdsl::int_vector<>> one_ranks_;
  /// @brief For each level and each back block,
  ///   contains the number of 1s up to (and including) the pointed-to area of
  ///   the back-block.
  std::vector<sdsl::int_vector<>> pointer_prefix_one_counts_;

  ~RecursiveBitBlockTree() {
    for (const IsInternalType* b : this->block_tree_types_) {
      delete b;
    }
    // in any other case, block_tree_types_ and block_tree_types_rs_ point to
    // the same object (a recursive block tree), so we may only free them once
    if constexpr (recursion_level == 0) {
      for (const RankSelect<OptimizedFor::ONE_QUERIES>* rs :
           this->block_tree_types_rs_) {
        delete rs;
      }
    }
    for (auto& ptrs : this->block_tree_pointers_) {
      delete ptrs;
    }
    for (auto& offsets : this->block_tree_offsets_) {
      delete offsets;
    }
  }

  [[nodiscard]] size_t height() const {
    return block_tree_types_.size();
  }

  [[nodiscard]] size_t size() const {
    return num_bits_;
  }

  bool operator[](const size_type bit_index) const {
    return access(bit_index);
  }

  bool access(const size_type bit_index) const {
    const int64_t byte_index = bit_index / 8;
    const int64_t bit_offset = bit_index % 8;

    int64_t block_size = block_size_lvl_[0];
    int64_t block_index = byte_index / block_size;
    int64_t off = byte_index % block_size;
    for (size_t i = 0; i < height(); i++) {
      const auto& is_internal = *block_tree_types_[i];
      const auto& is_internal_rank = *block_tree_types_rs_[i];
      const auto& pointers = *block_tree_pointers_[i];
      const auto& offsets = *block_tree_offsets_[i];
      if (!is_internal[block_index]) {
        // If this block is not internal, go to its pointed-to block
        const size_t back_block_index = is_internal_rank.rank0(block_index);
        off = off + offsets[back_block_index];
        block_index = pointers[back_block_index];
        if (off >= block_size) {
          ++block_index;
          off -= block_size;
        }
      }
      block_size /= tau_;
      const int64_t child = off / block_size;
      off %= block_size;
      block_index = is_internal_rank.rank1(block_index) * tau_ + child;
    }
    const uint8_t byte =
        decompress_map_[compressed_leaves_[block_index * leaf_size + off]];
    return ((1 << bit_offset) & byte) != 0;
  };

private:
  template <bool one>
  [[nodiscard]] size_t find_initial_block(const size_t rank) const {
    const auto& top_one_ranks = one_ranks_[0];
    const size_t block_size = block_size_lvl_[0];
    size_t start = (rank - 1) / (block_size * 8);
    size_t end = top_one_ranks.size() - 1;
    while (start != end) {
      const size_t middle = start + (end - start) / 2;
      size_t current_rank;
      if constexpr (one) {
        current_rank = (middle == 0) ? 0 : top_one_ranks[middle - 1];
      } else {
        const size_t middle_bits = middle * block_size * 8;
        current_rank =
            (middle == 0) ? 0 : middle_bits - top_one_ranks[middle - 1];
      }
      if (current_rank < rank) {
        if (start + 1 == end) {
          size_t bits;
          if constexpr (one) {
            bits = top_one_ranks[middle];
          } else {
            bits = (middle + 1) * block_size * 8 - top_one_ranks[middle];
          }
          // If there is only one block left, it's either the current or the
          // next block
          if (bits < rank) {
            start = middle + 1;
          }
          break;
        }
        start = middle;
      } else {
        end = middle - 1;
      }
    }
    return start;
  }

public:
  [[nodiscard("select result discarded")]] size_t select1(size_t rank) const {
    const auto& top_is_internal = *block_tree_types_[0];
    const auto& top_is_internal_rank = *block_tree_types_rs_[0];
    const auto& top_pointers = *block_tree_pointers_[0];
    const auto& top_offsets = *block_tree_offsets_[0];
    const auto& top_one_ranks = one_ranks_[0];
    size_t block_size = block_size_lvl_[0];

    // Binary Search for the correct top level block containing the correct 1
    size_t current_block = find_initial_block<true>(rank);

    size_t pos = (current_block * block_size * 8) - 1;
    // ReSharper disable once CppDFAUnreachableCode
    rank -= (current_block == 0) ? 0 : top_one_ranks[current_block - 1];

    // If that block is a back block, we need to move to the back-pointed block
    if (!top_is_internal[current_block]) {
      const size_t back_block_index = top_is_internal_rank.rank0(current_block);
      current_block = top_pointers[back_block_index];
      const size_t offset = top_offsets[back_block_index];
      size_t rank_d =
          (current_block == 0) ?
              top_one_ranks[current_block] :
              top_one_ranks[current_block] - top_one_ranks[current_block - 1];
      rank_d -= pointer_prefix_one_counts_[0][back_block_index];
      if (rank > rank_d) {
        rank -= rank_d;
        pos += (block_size - offset) * 8;
        ++current_block;
      } else {
        rank += pointer_prefix_one_counts_[0][back_block_index];
        pos -= offset * 8;
      }
    }

    size_t level = 1;
    while (level < height()) {
      const auto& pointer_ranks = pointer_prefix_one_counts_[level];
      const auto& is_internal = *block_tree_types_[level];
      const auto& is_internal_rank = *block_tree_types_rs_[level];
      const auto& prev_is_internal_rank = *block_tree_types_rs_[level - 1];
      const auto& offsets = *block_tree_offsets_[level];
      const auto& pointers = *block_tree_pointers_[level];
      const auto& one_ranks = one_ranks_[level];

      current_block = prev_is_internal_rank.rank1(current_block) * tau_;
      block_size /= tau_;
      const size_t start_block = current_block;
      while (one_ranks[current_block] < rank) {
        ++current_block;
      }
      rank -= (current_block == start_block) ? 0 : one_ranks[current_block - 1];
      pos += (current_block - start_block) * block_size * 8;
      if (!is_internal[current_block]) {
        size_t back_block_index = is_internal_rank.rank0(current_block);
        current_block = pointers[back_block_index];
        const size_t offset = offsets[back_block_index];
        size_t rank_d =
            (current_block % tau_ == 0) ?
                one_ranks[current_block] :
                one_ranks[current_block] - one_ranks[current_block - 1];
        rank_d -= pointer_ranks[back_block_index];
        if (rank > rank_d) {
          rank -= rank_d;
          pos += (block_size - offset) * 8;
          ++current_block;
        } else {
          rank += pointer_ranks[back_block_index];
          pos -= offset * 8;
        }
      }
      ++level;
    }

    current_block =
        block_tree_types_rs_[level - 1]->rank1(current_block) * tau_;
    size_t byte_offset = 0;
    while (rank > 0) {
      const uint8_t byte =
          decompress_map_[compressed_leaves_[current_block * leaf_size +
                                             byte_offset]];
      const uint8_t num_ones = std::popcount(byte);
      if (rank > num_ones) {
        rank -= num_ones;
        pos += 8;
        ++byte_offset;
      } else {
        for (size_t bit = 0; bit < 8 && rank > 0; ++bit) {
          pos++;
          rank -= ((1 << bit) & byte) > 0;
        }
      }
    }
    return pos;
  }

  [[nodiscard("select result discarded")]] size_t select0(size_t rank) const {
    const auto& top_is_internal = *block_tree_types_[0];
    const auto& top_is_internal_rank = *block_tree_types_rs_[0];
    const auto& top_pointers = *block_tree_pointers_[0];
    const auto& top_offsets = *block_tree_offsets_[0];
    const auto& top_one_ranks = one_ranks_[0];

    const size_t top_block_size = block_size_lvl_[0];
    const auto top_zero_ranks = [&top_one_ranks,
                                 top_block_size](const size_t i) -> size_t {
      return (i + 1) * top_block_size * 8 - top_one_ranks[i];
    };

    // Binary Search for the correct top level block containing the correct 1
    size_t current_block = find_initial_block<false>(rank);
    const size_t top_block_bits = top_block_size * 8;

    size_t pos = (current_block * top_block_bits) - 1;
    // ReSharper disable once CppDFAUnreachableCode
    rank -= (current_block == 0) ? 0 : top_zero_ranks(current_block - 1);
    // If that block is a back block, we need to move to the back-pointed block
    if (!top_is_internal[current_block]) {
      const size_t back_block_index = top_is_internal_rank.rank0(current_block);
      // const size_t child_block_bits =
      // height() == 1 ? leaf_size * 8 : block_size_lvl_[1] * 8;
      current_block = top_pointers[back_block_index];
      const size_t offset = top_offsets[back_block_index];
      const size_t prefix_bits = offset * 8;
      size_t rank_d =
          (current_block == 0) ?
              top_zero_ranks(current_block) :
              top_zero_ranks(current_block) - top_zero_ranks(current_block - 1);
      rank_d -= prefix_bits - pointer_prefix_one_counts_[0][back_block_index];
      if (rank > rank_d) {
        rank -= rank_d;
        pos += (top_block_size - offset) * 8;
        ++current_block;
      } else {
        rank += prefix_bits - pointer_prefix_one_counts_[0][back_block_index];
        pos -= offset * 8;
      }
    }

    size_t block_size = block_size_lvl_[0];
    size_t level = 1;
    while (level < height()) {
      const auto& pointer_ranks = pointer_prefix_one_counts_[level];
      const auto& is_internal = *block_tree_types_[level];
      const auto& is_internal_rank = *block_tree_types_rs_[level];
      const auto& prev_is_internal_rank = *block_tree_types_rs_[level - 1];
      const auto& offsets = *block_tree_offsets_[level];
      const auto& pointers = *block_tree_pointers_[level];
      const auto& one_ranks = one_ranks_[level];

      current_block = prev_is_internal_rank.rank1(current_block) * tau_;
      block_size /= tau_;

      const auto zero_ranks =
          [&one_ranks, this, block_size](const size_t i) -> size_t {
        const size_t rnk = (i % this->tau_ + 1) * block_size * 8 - one_ranks[i];
        return rnk;
      };
      const size_t start_block = current_block;
      while (zero_ranks(current_block) < rank) {
        ++current_block;
      }
      rank -=
          (current_block == start_block) ? 0 : zero_ranks(current_block - 1);
      pos += (current_block - start_block) * block_size * 8;
      if (!is_internal[current_block]) {
        size_t back_block_index = is_internal_rank.rank0(current_block);
        current_block = pointers[back_block_index];
        const size_t offset = offsets[back_block_index];
        const size_t prefix_bits = offset * 8;
        size_t rank_d =
            (current_block % tau_ == 0) ?
                zero_ranks(current_block) :
                zero_ranks(current_block) - zero_ranks(current_block - 1);
        rank_d -= prefix_bits - pointer_ranks[back_block_index];
        if (rank > rank_d) {
          rank -= rank_d;
          pos += (block_size - offset) * 8;
          ++current_block;
        } else {
          rank += prefix_bits - pointer_ranks[back_block_index];
          pos -= offset * 8;
        }
      }
      ++level;
    }

    current_block =
        block_tree_types_rs_[level - 1]->rank1(current_block) * tau_;
    size_t byte_offset = 0;
    while (rank > 0) {
      const uint8_t byte =
          decompress_map_[compressed_leaves_[current_block * leaf_size +
                                             byte_offset]];
      const uint8_t num_zeros = 8 - std::popcount(byte);
      if (rank > num_zeros) {
        rank -= num_zeros;
        pos += 8;
        byte_offset++;
      } else {
        for (size_t bit = 0; bit < 8 && rank > 0; ++bit) {
          pos++;
          rank -= ((1 << bit) & byte) == 0;
        }
      }
    }
    return pos;
  }

  /// @brief Counts the number of 1-bits up to (and excluding) an index.
  [[nodiscard("rank result discarded")]] size_t
  rank1(const size_type bit_index) const {
    const size_t byte_index = bit_index / 8;
    const auto& top_is_internal = *block_tree_types_[0];
    const auto& top_is_internal_rank = *block_tree_types_rs_[0];
    const auto& top_pointers = *block_tree_pointers_[0];
    const auto& top_offsets = *block_tree_offsets_[0];
    size_t block_size = block_size_lvl_[0];
    size_t block_index = byte_index / block_size;
    size_t block_offset = byte_index % block_size;
    size_t rank = (block_index == 0) ? 0 : one_ranks_[0][block_index - 1];
    if (!top_is_internal[block_index]) {
      // If the top block is a back block, go to it and adjust the offset
      const size_t back_block_index = top_is_internal_rank.rank0(block_index);
      rank -= pointer_prefix_one_counts_[0][back_block_index];
      block_offset += top_offsets[back_block_index];
      block_index = top_pointers[back_block_index];
      if (block_offset >= block_size) {
        // If we're exceeding the pointed-to block's offset,
        // add the ones inside of it
        rank +=
            (block_index == 0) ?
                one_ranks_[0][block_index] :
                (one_ranks_[0][block_index] - one_ranks_[0][block_index - 1]);
        ++block_index;
        block_offset -= block_size;
      }
    }

    // Go down to the next level
    block_size /= tau_;
    // How many children are we 'skipping over'
    size_t child = block_offset / block_size;
    block_offset %= block_size;
    block_index = top_is_internal_rank.rank1(block_index) * tau_ + child;

    size_t level = 1;
    while (level < height()) {
      const auto& ranks = one_ranks_[level];
      const auto& pointer_ranks = pointer_prefix_one_counts_[level];
      const auto& is_internal = *block_tree_types_[level];
      const auto& is_internal_rank = *block_tree_types_rs_[level];
      rank += (child == 0) ? 0 : ranks[block_index - 1];
      // If this block is internal, just go to the correct child
      if (is_internal[block_index]) {
        block_size /= tau_;
        child = block_offset / block_size;
        block_offset %= block_size;
        block_index = is_internal_rank.rank1(block_index) * tau_ + child;
        level++;
        continue;
      }

      // If we have a back block, we need to go to the pointed-to block
      const size_t back_block_index = is_internal_rank.rank0(block_index);
      rank -= pointer_ranks[back_block_index];
      block_offset += (*block_tree_offsets_[level])[back_block_index];
      block_index = (*block_tree_pointers_[level])[back_block_index];
      child = block_index % tau_;

      if (block_offset >= block_size) {
        // If we're exceeding the pointed-to block's offset,
        // add the ones inside of it and go to the next block
        rank += (child == 0) ? ranks[block_index] :
                               (ranks[block_index] - ranks[block_index - 1]);
        ++block_index;
        child = block_index % tau_;
        block_offset -= block_size;
      }
      const size_t remove_prefix = (child == 0) ? 0 : ranks[block_index - 1];
      rank -= remove_prefix;
    }

    // Number of leaves that exist before the leaves of the current block
    const size_type prefix_leaves = block_index - child;
    for (size_t block = 0; block < child * leaf_size; block++) {
      const uint8_t byte =
          decompress_map_[compressed_leaves_[prefix_leaves * leaf_size +
                                             block]];
      rank += std::popcount(byte);
    }
    for (size_t block = 0; block < block_offset; block++) {
      const uint8_t byte =
          decompress_map_[compressed_leaves_[block_index * leaf_size + block]];
      rank += std::popcount(byte);
    }

    // Masks to remove bits from the last byte,
    // that aren't part of the ran query
    static constexpr std::array<uint8_t, 8> MASKS = {
        0b0000'0000,
        0b0000'0001,
        0b0000'0011,
        0b0000'0111,
        0b0000'1111,
        0b0001'1111,
        0b0011'1111,
        0b0111'1111,
    };
    rank += std::popcount<uint8_t>(
        decompress_map_[compressed_leaves_[block_index * leaf_size +
                                           block_offset]] &
        MASKS[bit_index % 8]);
    return rank;
  }

  /// @brief Counts the number of 0-bits up to (and excluding) an index.
  size_t rank0(const size_type bit_index) const {
    return bit_index - rank1(bit_index);
  }

  [[nodiscard]] size_t print_space_usage() const {
    size_t space_usage = sizeof(tau_) + sizeof(max_leaf_length_) + sizeof(s_) +
                         sizeof(leaf_size);

    size_t delta_size = 0;
    for (const auto* bt : block_tree_types_) {
      if constexpr (types_is_block_tree) {
        space_usage += bt->print_space_usage();
        delta_size += bt->print_space_usage();
      } else {
        space_usage += bt->size() / 8;
        delta_size += bt->size() / 8;
      }
    }
#ifdef BT_DBG
    std::cout << "bv size: " << delta_size << std::endl;
    delta_size = 0;
#endif
    if constexpr (recursion_level == 0) {
      for (const auto* rs : block_tree_types_rs_) {
        space_usage += rs->space_usage();
        delta_size += rs->space_usage();
      }
#ifdef BT_DBG
      std::cout << "rs size: " << delta_size << std::endl;
#endif
    }
    delta_size = 0;
    for (const auto iv : block_tree_pointers_) {
      space_usage += sdsl::size_in_bytes(*iv);
      delta_size += sdsl::size_in_bytes(*iv);
      ;
    }
#ifdef BT_DBG
    std::cout << "ptrs size: " << delta_size << std::endl;
    delta_size = 0;
#endif
    for (const auto iv : block_tree_offsets_) {
      space_usage += sdsl::size_in_bytes(*iv);
      delta_size += sdsl::size_in_bytes(*iv);
    }
#ifdef BT_DBG
    std::cout << "offs size: " << delta_size << std::endl;
#endif
    space_usage += block_size_lvl_.size() *
                   sizeof(typename decltype(block_size_lvl_)::value_type);
    space_usage += block_per_lvl_.size() *
                   sizeof(typename decltype(block_per_lvl_)::value_type);

    if (rank_support) {
      for (auto& rs : one_ranks_) {
        space_usage += sdsl::size_in_bytes(rs);
      }
      for (auto& rs : pointer_prefix_one_counts_) {
        space_usage += sdsl::size_in_bytes(rs);
      }
    }

    // space_usage += leaves_.size() * sizeof(uint8_t);
    space_usage += sdsl::size_in_bytes(compressed_leaves_);
#ifdef BT_DBG
    std::cout << "leaf size: " << sdsl::size_in_bytes(compressed_leaves_)
              << std::endl;
#endif
    space_usage += compress_map_.size();

    return space_usage;
  };

  void
  add_bit_rank_support(size_t threads = std::thread::hardware_concurrency()) {
    if (rank_support) {
      return;
    }
    rank_support = true;

    if constexpr (recursion_level == 0) {
      threads = 1;
    }

    // Resize rank information vectors
    one_ranks_.resize(height(), sdsl::int_vector<0>());
    for (uint64_t level = 0; level < height(); level++) {
      one_ranks_[level].resize(block_tree_types_[level]->size());
    }
    pointer_prefix_one_counts_.resize(height(), sdsl::int_vector<0>());
    for (uint64_t level = 0; level < height(); level++) {
      pointer_prefix_one_counts_[level].resize(
          block_tree_pointers_[level]->size());
    }

    for (size_t block = 0; block < block_tree_types_[0]->size(); block++) {
      bit_rank_block(0, block);
    }

    for (size_t block = 1; block < block_tree_types_[0]->size(); block++) {
      one_ranks_[0][block] += one_ranks_[0][block - 1];
    }

#pragma omp parallel for default(none) num_threads(threads)
    for (size_t level = 1; level < height(); level++) {
      size_type counter = tau_;
      size_t acc = 0;
      for (size_t block = 0; block < one_ranks_[level].size(); block++) {
        const size_type ones_in_block = one_ranks_[level][block];
        acc += ones_in_block;
        one_ranks_[level][block] = acc;
        --counter;
        if (counter == 0) {
          acc = 0;
          counter = tau_;
        }
      }
    }
    for (auto& prefix_one_counts : pointer_prefix_one_counts_) {
      sdsl::util::bit_compress(prefix_one_counts);
    }
    for (auto& ranks : one_ranks_) {
      sdsl::util::bit_compress(ranks);
    }
  }

protected:
  void compress_leaves() {
    // Holds a 1 on every char that exists
    compress_map_.resize(256, 0);
    decompress_map_.resize(256, 0);
    for (size_t i = 0; i < this->leaves_.size(); ++i) {
      compress_map_[this->leaves_[i]] = 1;
    }
    for (size_t c = 0, cur_val = 0; c < this->compress_map_.size(); ++c) {
      const size_t tmp = compress_map_[c];
      compress_map_[c] = cur_val;
      decompress_map_[cur_val] = c;
      cur_val += tmp;
    }

    compressed_leaves_.resize(this->leaves_.size());
    for (size_t i = 0; i < this->leaves_.size(); ++i) {
      compressed_leaves_[i] = compress_map_[this->leaves_[i]];
    }
    sdsl::util::bit_compress(this->compressed_leaves_);
    leaves_.resize(0);
    leaves_.shrink_to_fit();
  }
  /// @brief Calculate the number of leading zeros for a 32-bit integer.
  /// This value is capped at 31.
  static size_type leading_zeros(const int32_t val) {
    return __builtin_clz(static_cast<unsigned int>(val) | 1);
  }

  /// @brief Calculate the number of leading zeros for a 64-bit integer.
  /// This value is capped at 64.
  static size_type leading_zeros(const int64_t val) {
    return __builtin_clzll(static_cast<unsigned long long>(val) | 1);
  }

  ///
  /// @brief Determine the padding and minimum height and the size of the blocks
  /// on the top level of a block tree with s top-level blocks and an arity of
  /// tau with leaves also of size tau.
  ///
  /// The height is the number of levels in the tree.
  /// The padding is the number of characters that the top-level exceeds the
  /// text length. For example, if the result was that the top level consists of
  /// s = 5 blocks of size 30 and the text size being 80, then the padding would
  /// be (5 * 30) - 80 = 70.
  ///
  /// @param[out] padding The number of characters in the last block (of the
  ///   first level of the tree) that are empty.
  /// @param[in] text_length The number of characters in the input string.
  /// @param[out] height The number of levels in the tree.
  /// @param[out] blk_size The size of blocks on the first level of the tree.
  ///
  void calculate_padding(int64_t& padding,
                         int64_t text_length,
                         int64_t& height,
                         int64_t& blk_size) {
    // This is the number of characters occupied by a tree with s*tau^h levels
    // and leaves of size tau. At the start, we only have a tree with the first
    // level with s leaf blocks which each have size tau. If we insert another
    // level, the number of leaf blocks (and therefore the number of occupied
    // characters) increases by a factor of tau.
    int64_t tmp_padding = this->s_ * this->tau_;
    int64_t h = 1;
    // Size of the blocks on the current level (starting at the leaf level)
    blk_size = tau_;
    // While the tree does not cover the entire text, add a level
    while (tmp_padding < text_length) {
      tmp_padding *= this->tau_;
      blk_size *= this->tau_;
      h++;
    }
    // once the tree has enough levels to cover the entire text, we set the
    // tree's values
    height = h;
    // The padding is the number of excess characters that the block tree covers
    // over the length of the text.
    padding = tmp_padding - text_length;
  }

  size_type bit_rank_block(size_type level, size_type block_index) {
    const auto& is_internal = *block_tree_types_[level];
    const auto& is_internal_rank = *block_tree_types_rs_[level];
    if (static_cast<uint64_t>(block_index) >= is_internal.size()) {
      return 0;
    }

    size_type num_ones = 0;
    if (is_internal[block_index]) {
      const size_type internal_index = is_internal_rank.rank1(block_index);
      if (static_cast<uint64_t>(level) < height() - 1) {
        // If we are not on the last level recursively call
        for (size_type k = 0; k < tau_; ++k) {
          num_ones += bit_rank_block(level + 1, internal_index * tau_ + k);
        }
      } else {
        // If we are on the last level
        for (size_type k = 0; k < tau_; ++k) {
          num_ones += bit_rank_leaf(internal_index * tau_ + k, leaf_size);
        }
      }
    } else {
      const size_type back_block_index = is_internal_rank.rank0(block_index);
      const size_type ptr = (*block_tree_pointers_[level])[back_block_index];
      const size_type off = (*block_tree_offsets_[level])[back_block_index];
      size_type num_ones_parts = 0;
      num_ones += one_ranks_[level][ptr];
      if (off > 0) {
        num_ones_parts = part_bit_rank_block(level, ptr, off);
        const size_type num_ones_2nd_part =
            part_bit_rank_block(level, ptr + 1, off);
        num_ones -= num_ones_parts;
        num_ones += num_ones_2nd_part;
      }
      pointer_prefix_one_counts_[level][back_block_index] = num_ones_parts;
    }
    one_ranks_[level][block_index] = num_ones;
    return num_ones;
  }

  size_type part_bit_rank_block(const size_type level,
                                const size_type block_index,
                                const size_type chars_to_process) {
    const auto& is_internal = *block_tree_types_[level];
    const auto& is_internal_rank = *block_tree_types_rs_[level];
    if (static_cast<uint64_t>(block_index) >= is_internal.size()) {
      return 0;
    }

    size_type num_ones = 0;
    if (is_internal[block_index]) {
      const size_type internal_index = is_internal_rank.rank1(block_index);
      size_type k = 0;
      size_type processed_chars = 0;
      if (static_cast<size_t>(level) < height() - 1) {
        const size_type child_size = block_size_lvl_[level + 1];
        // We're not on the last level
        // iterate over the children as long as we don't exceed the limit
        for (k = 0;
             k < tau_ && processed_chars + child_size <= chars_to_process;
             ++k) {
          num_ones += one_ranks_[level + 1][internal_index * tau_ + k];
          processed_chars += child_size;
        }

        // If we still need to process more chars and they end inside the next
        // child, rank that part of the next child
        if (processed_chars != chars_to_process) {
          num_ones += part_bit_rank_block(level + 1,
                                          internal_index * tau_ + k,
                                          chars_to_process - processed_chars);
        }
      } else {
        // We're on the last level
        for (k = 0; k < tau_ && processed_chars + leaf_size <= chars_to_process;
             ++k) {
          num_ones += bit_rank_leaf(internal_index * tau_ + k, leaf_size);
          processed_chars += leaf_size;
        }

        if (processed_chars != chars_to_process) {
          num_ones += bit_rank_leaf(internal_index * tau_ + k,
                                    chars_to_process % leaf_size);
        }
      }
    } else {
      const size_type back_block_index = is_internal_rank.rank0(block_index);
      const size_type ptr = (*block_tree_pointers_[level])[back_block_index];
      const size_type off = (*block_tree_offsets_[level])[back_block_index];

      // If we need to process chars beyond this block, we need to
      if (chars_to_process + off >= block_size_lvl_[level]) {
        // Ones in the entire block this block points to
        num_ones += one_ranks_[level][ptr];
        // Ones that overflow into the next block
        num_ones += part_bit_rank_block(level,
                                        ptr + 1,
                                        chars_to_process + off -
                                            block_size_lvl_[level]);
        // Num ones in the pointed-to block *before* the pointed-to area
        num_ones -= pointer_prefix_one_counts_[level][back_block_index];
      } else {
        // Number of ones up to the cutoff point
        num_ones += part_bit_rank_block(level, ptr, chars_to_process + off);
        // Num ones in the pointed-to block *before* the pointed-to area
        num_ones -= pointer_prefix_one_counts_[level][back_block_index];
      }
    }
    return num_ones;
  }

  ///
  /// @brief Count ones in leaf block.
  ///
  /// @param leaf_index The index of the leaf block.
  /// @param max_char_index The maximum character index (exclusive) to
  /// consider. This is used for when this block is at the end of the string.
  /// @return The number of ones in this block.
  ///
  size_type bit_rank_leaf(size_type leaf_index, size_type max_char_index) {
    if (static_cast<uint64_t>(leaf_index * leaf_size) >=
        compressed_leaves_.size()) {
      return 0;
    }

    size_type result = 0;
    for (size_type i = 0; i < max_char_index; ++i) {
      const uint8_t compressed_byte =
          compressed_leaves_[leaf_index * leaf_size + i];
      const uint8_t byte = decompress_map_[compressed_byte];
      result += std::popcount(byte);
    }
    return result;
  }
};

template <std::signed_integral size_type>
using BitBlockTree = RecursiveBitBlockTree<size_type>;

} // namespace pasta

/******************************************************************************/
