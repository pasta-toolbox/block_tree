/*******************************************************************************
 * This file is part of pasta::block_tree
 *
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

#include "pasta/bit_vector/bit_vector.hpp"
#include "pasta/block_tree/block_tree.hpp"
#include "pasta/block_tree/utils/MersenneHash.hpp"
#include "pasta/block_tree/utils/MersenneRabinKarp.hpp"
#include "pasta/block_tree/utils/sharded_util.hpp"

#include <ankerl/unordered_dense.h>
#include <memory>
#include <omp.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

using Clock = std::chrono::high_resolution_clock;
using TimePoint = Clock::time_point;
using Duration = Clock::duration;

__extension__ typedef unsigned __int128 uint128_t;

namespace pasta {

template <std::integral input_type, std::signed_integral size_type>
class BlockTreeFPParPH : public BlockTree<input_type, size_type> {
  using BitVector = pasta::BitVector;
  using Rank = pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>;

  /// A concurrent hash map
  /*template <typename key_type,
            typename value_type,
            typename hash_type = std::hash<key_type>,
            size_t num_submaps = 6,
            typename mutex_type = phmap::NullMutex>
  using HashMap = phmap::parallel_flat_hash_map<
      key_type,
      value_type,
      hash_type,
      phmap::priv::hash_default_eq<key_type>,
      phmap::priv::Allocator<
          typename phmap::priv::Pair<const key_type, value_type>>,
      num_submaps,
      mutex_type>;*/

  template <typename key_type,
            typename value_type,
            typename hash_type = std::hash<key_type>>
  using HashMap = ankerl::unordered_dense::map<key_type, value_type, hash_type>;

  /// A rabin karp hasher preconfigured for the current template parameters
  using RabinKarp = MersenneRabinKarp<input_type,
                                      size_type,
                                      internal::sharded::PRIME_EXPONENT>;
  /// A rabin karp hash for the preconfigured rabin karp hasher
  using RabinKarpHash = MersenneHash<input_type>;

  /// A hash map with rabin karp hashes as keys
  template <typename value_type>
  using RabinKarpMap =
      HashMap<RabinKarpHash, value_type, std::hash<RabinKarpHash>>;

  using LevelData = internal::sharded::LevelData<size_type, Rank>;
  using PairOccurrences = internal::sharded::PairOccurrences<size_type>;
  using BlockOccurrences = internal::sharded::BlockOccurrences<size_type>;

#ifdef BT_DBG
public:
  size_t bp_hash_pairs_ns = 0;
  size_t bp_scan_pairs_ns = 0;
  size_t bp_markings_ns = 0;
  size_t bp_bitvec_ns = 0;

  size_t b_hash_blocks_ns = 0;
  size_t b_scan_blocks_ns = 0;
  size_t b_update_blocks_ns = 0;

private:
#endif

  void construct(const std::vector<input_type>& text, const size_t threads) {
    const size_type text_len = text.size();
    /// The number of characters a block tree with s top-level blocks and arity
    /// of strictly tau would exceed over the text size
    int64_t padding;
    /// The height of the tree
    int64_t tree_height;
    /// The size of the largest blocks (i.e. the top level blocks)
    int64_t top_block_size;

    this->calculate_padding(padding, text_len, tree_height, top_block_size);

    const bool is_padded = padding > 0;

    std::vector<LevelData> levels;

    // Prepare the top level
    levels.emplace_back(0, top_block_size, text_len / top_block_size);
    LevelData& top_level = levels.back();
    top_level.block_starts->reserve(
        internal::sharded::ceil_div(text_len, top_level.block_size));
    for (size_type i = 0; i < text_len; i += top_level.block_size) {
      top_level.block_starts->push_back(i);
    }
    top_level.block_size = top_block_size;
    top_level.num_blocks = top_level.block_starts->size();

#ifdef BT_DBG
    size_t pairs_ns = 0;
    size_t blocks_ns = 0;
    size_t generate_ns = 0;
#endif

    // Construct the pre-pruned tree level by level
    for (size_t level = 0; level < static_cast<size_t>(tree_height); level++) {
#ifdef BT_DBG
      std::cout << "level " << level << std::endl;

      TimePoint now = Clock::now();
#endif
      LevelData& current = levels.back();
      scan_block_pairs(text, current, is_padded, threads);
#ifdef BT_DBG
      std::cout << "scanned block pairs " << std::endl;
      pairs_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                      Clock::now() - now)
                      .count();
      now = Clock::now();
#endif
      scan_blocks(text, current, is_padded, threads);
#ifdef BT_DBG
      std::cout << "scanned blocks " << std::endl;
      blocks_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                       Clock::now() - now)
                       .count();
      now = Clock::now();

#endif
      // Generate the next level (if we're not at the last level)
      if (level < static_cast<size_t>(tree_height) - 1 &&
          levels.back().block_size > this->max_leaf_length_ * this->tau_) {
        levels.push_back(std::move(generate_next_level(text, current)));
#ifdef BT_DBG
        generate_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                           Clock::now() - now)
                           .count();
#endif
      } else {
        break;
      }
    }
#ifdef BT_DBG
    TimePoint now = Clock::now();

    std::cout << "pairs: " << (pairs_ns / 1'000'000)
              << "ms,\n\thash pairs: " << (bp_hash_pairs_ns / 1'000'000)
              << "ms,\n\tscan pairs: " << (bp_scan_pairs_ns / 1'000'000)
              << "ms,\n\tmarkings: " << (bp_markings_ns / 1'000'000)
              << "ms,\n\tbitvec: " << (bp_bitvec_ns / 1'000'000)
              << "ms,\nblocks: " << (blocks_ns / 1'000'000)
              << "ms,\n\thash blocks: " << (b_hash_blocks_ns / 1'000'000)
              << "ms,\n\tscan blocks: " << (b_scan_blocks_ns / 1'000'000)
              << "ms,\n\tupdate blocks: " << (b_update_blocks_ns / 1'000'000)
              << "ms,\ngenerate_next: " << (generate_ns / 1'000'000) << "ms,"
              << std::endl;
#endif
    prune(levels);
#ifdef BT_DBG
    size_t prune_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - now)
            .count();
    now = Clock::now();

    std::cout << "prune: " << (prune_ns / 1'000'000) << "ms," << std::endl;
#endif
    make_tree(text, levels, padding);
#ifdef BT_DBG
    size_t make_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - now)
            .count();
    std::cout << "make: " << (make_ns / 1'000'000) << "ms" << std::endl;
#endif
  }

  /// @brief Scan through the blocks pairwise in order to identify which blocks
  /// should
  ///   be replaced with back blocks.
  ///
  /// @param text The input string.
  /// @param level The data for the current level.
  ///
  /// @return The block start indices for the next level of the tree
  void scan_block_pairs(const std::vector<input_type>& text,
                        LevelData& level,
                        const bool is_padded,
                        const size_t threads) {
    if (level.num_blocks < 4) {
      level.is_internal = std::make_unique<BitVector>(level.num_blocks, true);
      level.is_internal_rank = std::make_unique<Rank>(*level.is_internal);
      return;
    }

    // A map containing hashed block pairs mapped to their indices of the
    // pairs' first block respectively
    RabinKarpMap<PairOccurrences> map(level.num_blocks);

#ifdef BT_DBG
    TimePoint now = Clock::now();
#  pragma omp parallel default(none) num_threads(threads)                      \
      shared(level, map, text, now, is_padded)
#else
#  pragma omp parallel default(none) num_threads(threads)                      \
      shared(level, map, text, is_padded)
#endif
    {
      const size_t block_size = level.block_size;
      const size_t pair_size = 2 * block_size;
      const size_t num_blocks = level.num_blocks;
      const size_t num_block_pairs = num_blocks - 1 - is_padded;
      const auto& block_starts = *level.block_starts;
#pragma omp single
      {
        RabinKarp rk(text,
                     internal::sharded::SIGMA,
                     block_starts[0],
                     pair_size,
                     internal::sharded::PRIME);
        for (size_t i = 0; i < num_block_pairs; ++i) {
          // If the next block is not adjacent, we cannot hash the pair starting
          // at the current block
          if (!level.next_is_adjacent(i)) {
            continue;
          }
          rk.restart(block_starts[i]);
          // Move the hasher to the current block pair
          RabinKarpHash hash = rk.current_hash();
          // Try to find the hash in the map, insert a new entry if it doesn't
          // exist, and add the current block to the entry
          auto ptr = map.find(hash);
          if (ptr == map.end()) {
            auto [insert_ptr, _] = map.insert({hash, PairOccurrences(i)});
            ptr = insert_ptr;
          }
          ptr->second.add_block_pair(i);
          ptr->second.update(i);
        }
      }
#pragma omp barrier
#ifdef BT_DBG
#  pragma omp single
      {
        bp_hash_pairs_ns +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() -
                                                                 now)
                .count();
        now = Clock::now();
      }
#endif

      // Hash every window and determine for all block pairs whether they have
      // previous occurrences.
      size_t segment_size = std::max<size_t>(
          1,
          internal::sharded::ceil_div(num_block_pairs, omp_get_num_threads()));
      const size_t thread_id = omp_get_thread_num();

      // Start and end index of the current thread's segment
      const auto start = thread_id * segment_size;
      const auto end =
          std::min<size_t>(num_block_pairs, (thread_id + 1) * segment_size);

      if (start < static_cast<size_t>(num_block_pairs)) {
        RabinKarp rk(text,
                     internal::sharded::SIGMA,
                     block_starts[start],
                     pair_size,
                     internal::sharded::PRIME);
        for (size_t i = start; i < end; ++i) {
          if (!level.next_is_adjacent(i) | !level.next_is_adjacent(i + 1)) {
            continue;
          }
          scan_windows_in_block_pair(rk, map, block_size, i);
        }
      }
    }

#ifdef BT_DBG
    bp_scan_pairs_ns +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - now)
            .count();
    now = Clock::now();
#endif

    // Set up the packed array holding the markings for each block.
    // Each mark is a 2-bit number.
    // The MSB is 1 iff the block and its successor have a prior occurrence.
    // The LSB is 1 iff the block and its predecessor have a prior occurrence.
    sdsl::int_vector<2> markings(level.num_blocks, 0);
    for (auto it = map.begin(); it != map.end(); ++it) {
      const PairOccurrences& pair_occs = it->second;
      for (const size_type occ : pair_occs.occurrences) {
        if (pair_occs.first_occ_block < occ) {
          markings[occ] = markings[occ] | 0b10;
          markings[occ + 1] = markings[occ + 1] | 0b01;
        }
      }
    }
#ifdef BT_DBG
    bp_markings_ns +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - now)
            .count();
    now = Clock::now();
#endif

    // Generate the bit vector indicating which blocks are internal
    level.is_internal = std::make_unique<BitVector>(level.num_blocks);

    auto& is_internal = *level.is_internal;
    is_internal[0] = true;
    is_internal[level.num_blocks - 1] = markings[level.num_blocks - 1] != 0b01;
    for (size_type i = 0; i < level.num_blocks - 1; ++i) {
      const bool block_is_internal = markings[i] != 0b11;
      is_internal[i] = block_is_internal;
    }
#ifdef BT_DBG
    bp_bitvec_ns +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - now)
            .count();
#endif
    level.is_internal_rank = std::make_unique<Rank>(*level.is_internal);
  }

  /// @brief Scan through the windows starting in a block and mark
  /// them
  ///   accordingly if they represent the earliest occurrence of some block
  ///   hash.
  ///
  /// The supplied `RabinKarp` hasher must be at the start of the block.
  /// @param rk A Rabin-Karp hasher whose state is at the start of the block.
  /// @param map The map containing the hashes of block pairs mapped to their
  ///   block indexes at which they occur.
  /// @param current_block_index The index of the block being currently hashed.
  static inline void
  scan_windows_in_block_pair(RabinKarp& rk,
                             RabinKarpMap<PairOccurrences>& map,
                             const size_t block_size,
                             const size_type current_block_index) {
    for (size_t offset = 0; offset < block_size; ++offset, rk.next()) {
      RabinKarpHash current_hash = rk.current_hash();
      // Find the hash of the current window among the hashed block pairs.
      auto found = map.find(current_hash);
      if (found == map.end()) {
        continue;
      }
      PairOccurrences& occurrences = found->second;
      occurrences.update(current_block_index);
    }
  }

  /// @brief Determine the positions for each block's earliest occurrence if
  /// there is any.
  ///
  /// @param s The input text
  /// @param level_data The data for the current level
  /// @param is_padded true, iff the last block of the level extends past the
  ///   end of the text
  void scan_blocks(const std::vector<input_type>& text,
                   LevelData& level_data,
                   const bool is_padded,
                   const size_t threads) {
    const size_t num_blocks = level_data.num_blocks;

    level_data.pointers = std::make_unique<std::vector<size_type>>(
        num_blocks,
        internal::sharded::NO_EARLIER_OCC);
    level_data.offsets =
        std::make_unique<std::vector<size_type>>(num_blocks, 0);
    level_data.counters =
        std::make_unique<std::vector<size_type>>(num_blocks, 0);

    if (num_blocks <= 2) {
      return;
    }

    // A map with hashed slices as keys, which map to a vector of links,
    // describing a link between a (potential) back block to their source block.
    // In addition to the vector, there is a boolean which denotes whether a
    // hash has already been processed
    RabinKarpMap<BlockOccurrences> links(num_blocks);

#ifdef BT_DBG
    TimePoint now = Clock::now();
#  pragma omp parallel default(none) num_threads(threads)                      \
      shared(level_data, text, links, now, is_padded, std::cout)
#else
#  pragma omp parallel default(none) num_threads(threads)                      \
      shared(level_data, text, links, is_padded)
#endif
    {
      const std::vector<size_type>& block_starts = *level_data.block_starts;
#pragma omp single
      for (size_type i = 0; i < level_data.num_blocks - is_padded; ++i) {
        const RabinKarp rk(text,
                           internal::sharded::SIGMA,
                           block_starts[i],
                           level_data.block_size,
                           internal::sharded::PRIME);
        const RabinKarpHash hash = rk.current_hash();
        auto ptr = links.find(hash);
        if (ptr == links.end()) {
          auto [insert_ptr, _] = links.emplace(hash, BlockOccurrences(i));
          ptr = insert_ptr;
        }

        ptr->second.add_block(i);
        ptr->second.update(i, 0);
      }
#pragma omp barrier

#ifdef BT_DBG
#  pragma omp single
      {
        b_hash_blocks_ns +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() -
                                                                 now)
                .count();
        now = Clock::now();
      }
#endif

      const size_t num_total_iterations = level_data.num_blocks - is_padded - 1;

      const size_t thread_id = omp_get_thread_num();
      const size_t segment_size =
          internal::sharded::ceil_div(num_total_iterations,
                                      omp_get_num_threads());
      const size_t start = thread_id * segment_size;
      const size_t end = std::min<size_t>(num_total_iterations,
                                          (thread_id + 1) * segment_size);

      // Hash every window and find the first occurrences for every block.
      if (start < block_starts.size() - is_padded) {
        RabinKarp rk(text,
                     internal::sharded::SIGMA,
                     block_starts[start],
                     level_data.block_size,
                     internal::sharded::PRIME);
        for (size_t i = start; i < end; ++i) {
          if (!level_data.next_is_adjacent(i)) {
            continue;
          }
          if (static_cast<int64_t>(rk.init_) != block_starts[i]) {
            rk.restart(block_starts[i]);
          }
          scan_windows_in_block(rk, links, level_data, i);
        }
      }
    }
#ifdef BT_DBG
    b_scan_blocks_ns +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - now)
            .count();
    now = Clock::now();
#endif

    // By this point, the map should contain the first occurrences of every
    // respective block's content. We then fill the pointers and offsets with
    // this data and increment counters accordingly
    for (auto it = links.cbegin(); it != links.cend(); ++it) {
      // The occurrences of all blocks with a given hash
      const BlockOccurrences& occs = it->second;
      auto first_occ = occs.first_occ.load();
      for (const size_type occ : occs.occurrences) {
        if (occ == first_occ.block ||
            (first_occ.offset > 0 && occ == first_occ.block + 1)) {
          continue;
        }
        (*level_data.pointers)[occ] = first_occ.block;
        (*level_data.offsets)[occ] = first_occ.offset;
        const bool is_back_block = !(*level_data.is_internal)[occ];
        (*level_data.counters)[first_occ.block] += 1;
        (*level_data.counters)[first_occ.block + 1] +=
            is_back_block && (first_occ.offset > 0);
      }
    }
#ifdef BT_DBG
    b_update_blocks_ns +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - now)
            .count();
#endif
  }
  /// @brief Scans through block-sized windows starting inside one block and
  ///     tries to find earlier occurrences of blocks. Non-internal blocks will
  ///     have their respective m_source_blocks and m_offsets entries populated.
  /// @param rk A Rabin-Karp hasher whose current state is at a block start.
  /// @param links A map whose keys are hashed blocks and the values
  ///   are all block indices of blocks matching the hash in ascending order.
  /// @param current_block_internal_index The index of the block which the
  ///   Rabin-Karp hasher is situated in only with respect to *internal blocks*
  ///   on the current level, disregarding back blocks.
  /// @param num_hashes The number of times the Rabin-Karp hasher should hash.
  static void scan_windows_in_block(RabinKarp& rk,
                                    RabinKarpMap<BlockOccurrences>& links,
                                    LevelData& level_data,
                                    const size_type current_block_index) {
    for (size_type offset = 0; offset < level_data.block_size;
         ++offset, rk.next()) {
      const RabinKarpHash hash = rk.current_hash();
      // Find all blocks in the multimap that match our hash
      auto found = links.find(hash);
      if (found == links.end()) {
        continue;
      }
      found->second.update(current_block_index, offset);
    }
  }

  /// @brief Generate the block size, number of block and block start indices
  /// for the next level.
  ///
  /// This depends on the current level's block size, number of blocks and
  /// is_internal bit vector.
  ///
  /// @param text The input text.
  /// @param level The level data of the previous level.
  /// @return The level data of the next level.
  [[nodiscard]] LevelData
  generate_next_level(const std::vector<input_type>& text,
                      const LevelData& level) const {
    const size_t block_size = level.block_size;
    const size_t num_blocks = level.num_blocks;
    const auto& is_internal = *level.is_internal;
    const size_t next_block_size = block_size / this->tau_;

    std::vector<size_type> new_block_starts;
    new_block_starts.reserve(num_blocks * this->tau_);
    for (size_t i = 0; i < num_blocks; ++i) {
      if (!is_internal[i]) {
        continue;
      }

      // We generate up to tau new blocks for each internal block,
      // excluding blocks that start past the end of the text
      const auto parent_block_start = (*level.block_starts)[i];
      for (size_t j = 0, current_block_start = parent_block_start;
           j < static_cast<size_t>(this->tau_) &&
           current_block_start < text.size();
           ++j, current_block_start += next_block_size) {
        new_block_starts.push_back(current_block_start);
      }
    }

    LevelData next_level(level.level_index + 1,
                         next_block_size,
                         new_block_starts.size());
    next_level.block_starts =
        std::make_unique<std::vector<size_type>>(std::move(new_block_starts));
    return next_level;
  }

  ///
  /// @brief Takes a vector of levels and fills the block tree fields with them.
  ///
  /// @param[in] levels A vector containing data for each level, with the first
  /// entry corresponding to the topmost level.
  ///
  void make_tree(const std::vector<input_type>& text,
                 std::vector<LevelData>& levels,
                 int64_t padding) {
    const bool is_padded = padding > 0;

    // Count the current number of internal blocks per level
    std::vector<size_type> new_num_internal(levels.size(), 0);
    for (size_t level = 0; level < levels.size(); level++) {
      for (size_t block = 0; block < levels[level].is_internal->size();
           block++) {
        if ((*levels[level].is_internal)[block]) {
          new_num_internal[level]++;
        }
      }
    }

    // Create first level
    bool found_back_block = levels[0].is_internal->size() >
                                static_cast<size_t>(new_num_internal[0]) ||
                            !this->CUT_FIRST_LEVELS;
    LevelData& top_level = levels.front();
    if (found_back_block) {
      const size_t n = top_level.num_blocks;
      const size_t num_internal = new_num_internal[0];
      auto pointers = new sdsl::int_vector<>(n - num_internal, 0);
      auto offsets = new sdsl::int_vector<>(n - num_internal, 0);
      size_t num_back_blocks = 0;
      for (size_t i = 0; i < n; i++) {
        // if a back block is found, add its pointer and offset
        if (!(*top_level.is_internal)[i]) {
          (*pointers)[num_back_blocks] = (*top_level.pointers)[i];
          (*offsets)[num_back_blocks] = (*top_level.offsets)[i];
          num_back_blocks++;
        }
      }
      sdsl::util::bit_compress(*pointers);
      sdsl::util::bit_compress(*offsets);
      this->block_tree_types_.push_back(top_level.is_internal.release());
      this->block_tree_types_rs_.push_back(
          new Rank(*this->block_tree_types_.back()));
      this->block_tree_pointers_.push_back(pointers);
      this->block_tree_offsets_.push_back(offsets);
      this->block_size_lvl_.push_back(top_level.block_size);
    }
    top_level.pointers.reset();
    top_level.offsets.reset();
    top_level.counters.reset();

    // Add level data to the tree
    for (size_t level_index = 1; level_index < levels.size(); level_index++) {
      LevelData& level = levels[level_index];
      LevelData& previous_level = levels[level_index - 1];
      found_back_block |= static_cast<size_t>(new_num_internal[level_index]) <
                          levels[level_index].is_internal->size();
      if (!found_back_block) {
        level.is_internal.reset();
        level.is_internal_rank.reset();
        level.pointers.reset();
        level.offsets.reset();
        level.counters.reset();
        previous_level.block_starts.reset();
        continue;
      }

      make_tree_level(levels,
                      new_num_internal,
                      level_index,
                      is_padded,
                      text.size());

      // We don't need these anymore
      if (level_index < levels.size() - 1) {
        level.is_internal.reset();
      }
      level.is_internal_rank.reset();
      level.pointers.reset();
      level.offsets.reset();
      level.counters.reset();
      previous_level.block_starts.reset();
    }

    this->leaf_size = levels.back().block_size / this->tau_;
    // Construct the leaf string
    int64_t leaf_count = 0;
    auto& last_is_internal = *levels.back().is_internal;
    std::vector<size_type>& last_block_starts = *levels.back().block_starts;
    for (size_t block = 0; block < last_is_internal.size(); block++) {
      if (!last_is_internal[block]) {
        continue;
      }
      const size_type block_start = last_block_starts[block];
      // For every leaf on the last level, we have tau leaf blocks
      leaf_count += this->tau_;
      // Iterate through all characters in this child and
      // add them to the leaf string
      for (size_t b = 0; b < static_cast<size_t>(this->leaf_size * this->tau_);
           b++) {
        if (static_cast<size_t>(block_start + b) < text.size()) {
          this->leaves_.push_back(text[block_start + b]);
        }
      }
    }
    this->amount_of_leaves = leaf_count;
    this->compress_leaves();
  }

  /// @brief Generates a level and adds the relevant data to the block tree.
  ///
  /// @param levels The vector of levels of the tree.
  /// @param level_index The index of the level to generate. This must be
  ///   strictly greater than 0.
  /// @param is_padded Whether there is padding in the last block of the tree
  void make_tree_level(std::vector<LevelData>& levels,
                       const std::vector<size_type>& new_num_internal,
                       const size_t level_index,
                       const bool is_padded,
                       const size_t text_len) {
    LevelData& previous_level = levels[level_index - 1];
    LevelData& level = levels[level_index];

    size_type new_size =
        (new_num_internal[level_index - 1] - is_padded) * this->tau_;
    // Determine the number of children the last block generated
    if (is_padded) {
      const size_type last_block_parent_start =
          previous_level.block_starts->back();
      const size_type block_size = level.block_size;
      new_size +=
          internal::sharded::ceil_div(text_len - last_block_parent_start,
                                      block_size);
    }
    previous_level.block_starts.reset();
    const size_type num_internal = new_num_internal[level_index];

    // Allocate new vectors for the tree
    auto* is_internal = new BitVector(new_size);
    auto* pointers = new sdsl::int_vector<>(new_size - num_internal, 0);
    auto* offsets = new sdsl::int_vector<>(new_size - num_internal, 0);

    // Number of non-pruned blocks before the current block
    size_type num_non_pruned = 0;
    // Number of back blocks before the current block
    size_type num_back_blocks = 0;
    // Number of pruned blocks before the current block
    size_type num_pruned = 0;

    // We will reuse the allocated memory of the pointers vector to store
    // the number of pruned blocks before the block.
    // The invariant is that all values up to i are overwritten while all
    // values starting after i will still be valid pointers
    // This contains the number of pruned blocks before the block i
    std::vector<size_type>& prefix_pruned_blocks = *level.pointers;
    for (size_type i = 0; i < level.num_blocks; i++) {
      const size_type ptr = (*level.pointers)[i];
      prefix_pruned_blocks[i] = num_pruned;

      // If the current block is not pruned, add it to the new tree
      if (ptr == internal::sharded::PRUNED) {
        num_pruned++;
        continue;
      }

      // Add it to the is_internal bit vector
      const bool block_is_internal = (*level.is_internal)[i];
      (*is_internal)[num_non_pruned] = block_is_internal;
      num_non_pruned++;

      if (block_is_internal) {
        continue;
      }

      // If it is a back block, add its pointer and offset
      const size_type offset = (*level.offsets)[i];

      (*pointers)[num_back_blocks] = ptr - prefix_pruned_blocks[ptr];
      (*offsets)[num_back_blocks] = offset;
      num_back_blocks++;
    }

    sdsl::util::bit_compress(*pointers);
    sdsl::util::bit_compress(*offsets);
    this->block_tree_types_.push_back(is_internal);
    this->block_tree_types_rs_.push_back(new Rank(*is_internal));
    this->block_tree_pointers_.push_back(pointers);
    this->block_tree_offsets_.push_back(offsets);
    this->block_size_lvl_.push_back(level.block_size);
  }

  /// @brief Prunes the tree of unnecessary nodes.
  /// @param levels The levels of the tre represented as a vector of levels.
  void prune(std::vector<LevelData>& levels) {
    // We need to traverse the block tree in post order,
    // handling children from right to left.
    for (int block_index = levels[0].num_blocks - 1; block_index >= 0;
         --block_index) {
      prune_block(levels, 0, block_index);
    }
  }

  /// @brief Prunes a block and its descendants of unnecessary internal nodes.
  /// @param levels The WIP levels of the tree.
  /// @param level_index The level of the block to prune.
  /// @param block_index The index of the block to prune.
  /// @return Whether this block is/stays internal after the pruning process
  bool prune_block(std::vector<LevelData>& levels,
                   const size_t level_index,
                   const size_t block_index) {
    LevelData& level = levels[level_index];
    BitVector& is_internal = *level.is_internal;

    // If the current block is a back block already, there is nothing to prune
    if (!is_internal[block_index]) {
      return false;
    }

    const size_type first_child =
        level.is_internal_rank->rank1(block_index) * this->tau_;

    bool has_internal_children = false;

    // On the last level, all blocks just have leaves as children,
    // none of which can be pointed to. So only recurse, if we are not on the
    // last level.
    if (level_index < levels.size() - 1) {
      const size_type last_child =
          std::min<size_type>(first_child + this->tau_ - 1,
                              levels[level_index + 1].is_internal->size() - 1);
      // Iterate through children in reverse
      for (size_type child = last_child; child >= first_child; --child) {
        has_internal_children |= prune_block(levels, level_index + 1, child);
      }
    }

    // If any of the children is internal, this block stays internal as well
    if (has_internal_children) {
      return true;
    }

    const size_type pointer = (*level.pointers)[block_index];
    const size_type offset = (*level.offsets)[block_index];
    const size_type counter = (*level.counters)[block_index];
    // If there is no earlier occurrence or there are blocks pointing to this,
    // then this must stay internal
    if (pointer == internal::sharded::NO_EARLIER_OCC || counter > 0) {
      return true;
    }

    // Now we know that there is an earlier occurrence,
    // and nothing is pointing here.
    // We will make this block here into a back block...
    is_internal[block_index] = false;
    (*level.counters)[pointer] += 1;
    (*level.counters)[pointer + 1] += offset > 0;

    if (level_index == levels.size() - 1) {
      return false;
    }

    // ...and mark the children as pruned
    LevelData& child_level = levels[level_index + 1];
    const size_type last_child =
        std::min<size_type>(first_child + this->tau_ - 1,
                            child_level.is_internal->size() - 1);

    for (size_type child = last_child; child >= first_child; --child) {
      const size_type child_pointer = (*child_level.pointers)[child];
      const size_type child_offset = (*child_level.offsets)[child];
      // Decrement the counter of where the child points
      (*child_level.counters)[child_pointer] -= 1;
      (*child_level.counters)[child_pointer + 1] -= child_offset > 0;
      // Mark the child as pruned
      (*child_level.pointers)[child] = internal::sharded::PRUNED;
    }

    return false;
  }

public:
  BlockTreeFPParPH(const std::vector<input_type>& text,
                   const size_t arity,
                   const size_t root_arity,
                   const size_t max_leaf_length,
                   const size_t threads) {
    const auto old = omp_get_max_threads();
    const auto old_dynamic = omp_get_dynamic();
    omp_set_dynamic(0);
    omp_set_num_threads(static_cast<int>(threads));
    this->tau_ = arity;
    this->s_ = root_arity;
    this->max_leaf_length_ = max_leaf_length;
    this->map_unique_chars(text);
    construct(text, threads);
    omp_set_dynamic(old_dynamic);
    omp_set_num_threads(old);
  }

  ~BlockTreeFPParPH() {
    for (auto& rank : this->block_tree_types_rs_) {
      delete rank;
    }
    for (auto& bv : this->block_tree_types_) {
      delete bv;
    }
    for (auto& ptrs : this->block_tree_pointers_) {
      delete ptrs;
    }
    for (auto& offsets : this->block_tree_offsets_) {
      delete offsets;
    }
  }

  /// @brief Validates that a back-pointer actually points to the same text
  /// content.
  /// @param text The input text.
  /// @param level_index The index of the current level.
  /// @param block_index The block index.
  /// @param block_start The start index of the block's content in the text.
  /// @param source_start The start index of the source block's content in the
  /// text.
  /// @param source_pointer The block index of the source block.
  /// @param source_offset The offset from which the block copies out of the
  /// source block.
  /// @param block_size The block size.
  /// @return `true`, iff the pointer is valid. false otherwise
  bool debug_validate_pointer(const std::vector<input_type>& text,
                              const size_type level_index,
                              const size_type block_index,
                              const size_type block_start,
                              const size_type source_start,
                              const size_type source_pointer,
                              const size_type source_offset,
                              const size_type block_size) const {
    if (source_start + block_size > block_start) {
      std::cerr << "source overlapping block on level " << level_index
                << ":\n\tBlock Start: " << block_start
                << "\n\tSource Start: " << source_start
                << "\n\tBlock Size: " << block_size
                << "\n\tBlock: " << block_index
                << "\n\tSource Block: " << source_pointer
                << "\n\tSource Offset: " << source_offset << std::endl;
      return false;
    }
    for (size_type i = 0; i < block_size; i++) {
      if (text[block_start + i] != text[source_start + i]) {
        std::cerr << "source block mismatch on level " << level_index << ": "
                  << "\n\tBlock Start: " << block_start
                  << "\n\tSource Start: " << source_start
                  << "\n\tBlock Size: " << block_size
                  << "\n\tBlock: " << block_index
                  << "\n\tSource Block: " << source_pointer
                  << "\n\tSource Offset: " << source_offset << std::endl;
        return false;
      };
    }
    return true;
  }
}; // namespace pasta

} // namespace pasta