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

#include "pasta/block_tree/rec_block_tree.hpp"
#include "pasta/block_tree/utils/MersenneHash.hpp"
#include "pasta/block_tree/utils/MersenneRabinKarp.hpp"

#include <ankerl/unordered_dense.h>
#include <iostream>

__extension__ typedef unsigned __int128 uint128_t;

namespace pasta {

template <typename key_type,
          typename value_type,
          typename hash_type = std::hash<key_type>>
using HashMap = ankerl::unordered_dense::map<key_type, value_type, hash_type>;

template <typename input_type, typename size_type>
class BlockTreeFP : public BlockTree<input_type, size_type> {
public:
  size_type const_size = 0;
  size_type sigma_ = 0;
  ///
  /// @brief Prune this block and all of its children.
  ///
  /// This will set the pointers for each node that is pruned to PRUNED (== -2)
  /// and will increase and decrease pointers for each new back block or pruned
  /// back block respectively.
  ///
  /// @param[in] counter For each level (starting at the topmost level) and
  ///   block, contains the number of back blocks pointing to the block
  /// @param[in] pointer For each level (starting at the topmost level) and
  ///   block, contains the earliest block index in which the content of this
  /// @param[in] offset For each level (starting at the topmost level) and
  ///   block, contains the character offset at which this block's content can
  /// @param[in] marked_tree For each level (starting at the topmost level) and
  ///   block, has a 1 if the block is internal, 0 if it is a back block.
  /// @param[out] pruned_tree The resulting tree after pruning, with the same
  ///   content as marked_tree (as far as I can see this is unused)
  /// @param[in] i The current level, 0 being the top level
  /// @param[in] j The current block being pruned
  /// @param[in] ranks Rank data structures for each bit vector in marked_tree.
  /// @return true, if an this block was and stays internal, false otherwise
  ///
  bool prune_block(
      std::vector<std::vector<size_type>>& counter,
      std::vector<std::vector<size_type>>& pointer,
      std::vector<std::vector<size_type>>& offset,
      std::vector<pasta::BitVector*>& marked_tree,
      std::vector<pasta::BitVector*>& pruned_tree,
      size_type i,
      size_type j,
      std::vector<pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>>& ranks) {
    // String leaf children can always be pruned,
    // since they are contained in the leaf string.
    // Fully padded children don't exist and can be ignored/sanity check
    // assumes short circuit evaluation as compiler behaviour
    if (static_cast<uint64_t>(i) >= marked_tree.size() ||
        static_cast<uint64_t>(j) >= marked_tree[i]->size()) {
      return false;
    }

    bool marked_children = false;
    // We already incremented counters for back blocks during construction
    // and now only need to consider internal blocks
    if ((*marked_tree[i])[j] == 1) {
      // inverse postorder dfs. We handle children from back to front.
      // We need the rank over the internal nodes of this level,
      // because only the internal nodes generated children on the
      // next level.
      size_type rank_blk = ranks[i].rank1(j);
      for (size_type k = this->tau_ - 1; k >= 0; k--) {
        marked_children |= prune_block(counter,
                                       pointer,
                                       offset,
                                       marked_tree,
                                       pruned_tree,
                                       i + 1,
                                       rank_blk * this->tau_ + k,
                                       ranks);
      }
      // Conditions to be pruned are:
      // - no internal children,
      // - no pointers pointing to me
      // - an earlier occurrence in the text
      if (!marked_children && counter[i][j] == 0 &&
          pointer[i][j] != NO_FORMER_OCC) {
        // This block is no longer internal
        (*marked_tree[i])[j] = 0;
        // Since this is a back block now, we need to increment the counters for
        // the blocks this back block now points to
        counter[i][pointer[i][j]]++;
        if (offset[i][j] > 0) {
          counter[i][pointer[i][j] + 1]++;
        }
        // We only remove children if we're not at the last level
        if (static_cast<uint64_t>(i + 1) < counter.size()) {
          // Remove all of its children by decrementing counters
          // and marking them as PRUNED
          for (size_type k = this->tau_ - 1; k >= 0; k--) {
            // If this child node actually exists,
            // decrement the counters of the blocks the child points to and set
            if (static_cast<uint64_t>(rank_blk * this->tau_) + k <
                counter[i + 1].size()) {
              auto ptr_child = pointer[i + 1][rank_blk * this->tau_ + k];
              counter[i + 1][ptr_child]--;
              if (offset[i + 1][rank_blk * this->tau_ + k] > 0) {
                counter[i + 1][ptr_child + 1]--;
              }
              // Set this child's pointer to PRUNED
              pointer[i + 1][rank_blk * this->tau_ + k] = PRUNED;
            }
          }
        }
      }
    }
    // If this node has internal children, has other blocks pointing to itself
    // or has no earlier occurrence, then it remains internal
    return marked_children || counter[i][j] > 0 ||
           pointer[i][j] == NO_FORMER_OCC;
  };

  ///
  /// @brief Prunes the block tree
  ///
  /// @param[in] counter For each level (starting at the topmost level) and
  ///   block, contains the number of back blocks pointing to the block
  /// @param[in] pointer For each level (starting at the topmost level) and
  ///   block, contains the earliest block index in which the content of this
  ///   block can be found.
  /// @param[in] offset For each level (starting at the topmost level) and
  ///   block, contains the character offset at which this block's content can
  ///   be found in the block the back-pointer points to.
  /// @param[in] marked_tree For each level (starting at the topmost level) and
  ///   block, has a 1 if the block is internal, 0 if it is a back block.
  /// @param[out] pruned_tree The resulting tree after pruning, with the same
  ///   content as marked_tree
  /// @return 0
  ///
  int32_t pruning_extended(std::vector<std::vector<size_type>>& counter,
                           std::vector<std::vector<size_type>>& pointer,
                           std::vector<std::vector<size_type>>& offset,
                           std::vector<pasta::BitVector*>& marked_tree,
                           std::vector<pasta::BitVector*>& pruned_tree) {
    std::vector<pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>> ranks;
    for (auto bv : marked_tree) {
      ranks.push_back(pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
    }
    auto& top_lvl = *marked_tree[0];
    /// Prune the blocks on the top level from back to front
    for (size_type j = top_lvl.size() - 1; j >= 0; j--) {
      prune_block(counter,
                  pointer,
                  offset,
                  marked_tree,
                  pruned_tree,
                  0,
                  j,
                  ranks);
    }
    return 0;
  }

  int32_t pruning_simple(std::vector<pasta::BitVector*>& first_pass_bv,
                         std::vector<std::vector<int64_t>>& blk_lvl,
                         std::vector<pasta::BitVector*>& bv_pass_2,
                         std::vector<std::vector<size_type>>& pass1_pointer,
                         std::vector<std::vector<size_type>>& pass1_offset,
                         std::vector<std::vector<size_type>>& pass2_pointer,
                         std::vector<std::vector<size_type>>& pass2_offset,
                         std::vector<size_type>& pass2_max_pointer,
                         std::vector<size_type>& pass2_max_offset,
                         std::vector<size_type>& pass2_ones,
                         int64_t& block_size) {
    for (int64_t i = first_pass_bv.size() - 1; i >= 0; i--) {
      auto* bv = new pasta::BitVector(blk_lvl[i].size(), 0);
      size_type marked_counter = 0;
      if (static_cast<uint64_t>(i) != first_pass_bv.size() - 1) {
        for (uint64_t j = 0; j < bv->size(); j++) {
          if ((*first_pass_bv[i])[j] == 1) {
            for (size_type k = 0; k < this->tau_; k++) {
              if ((*bv_pass_2[bv_pass_2.size() -
                              1])[marked_counter * this->tau_ + k] == 1) {
                (*bv)[j] = 1;
              }
            }
            marked_counter++;
          }
        }
      }
      auto pointers = std::vector<size_type>();
      auto offsets = std::vector<size_type>();
      size_type max_pointer = 0;
      for (size_type j = blk_lvl[i].size() - 1; j >= 0; j--) {
        if ((*bv)[j] == 1) {
          continue;
        }
        bool has_ptr = false;
        auto ptr = pass1_pointer[i][j];
        auto off = pass1_offset[i][j];
        if (ptr == -1 || ptr + std::min(off, (size_type)1) >= j) {
          (*bv)[j] = 1;
        } else {
          size_type b = ptr;
          size_type current_offset = off;
          (*bv)[b] = 1;
          if (current_offset != 0) {
            (*bv)[b + 1] = 1;
          }
          pointers.push_back(b);
          offsets.push_back(current_offset);
          if (b > max_pointer) {
            max_pointer = b;
          }
          has_ptr = true;
        }
        if (!has_ptr) {
          (*bv)[j] = 1;
        }
      }
      pass2_max_offset.push_back(block_size);
      block_size *= this->tau_;
      bv_pass_2.push_back(bv);

      size_type ones_per_pass2 = 0;
      for (uint64_t j = 0; j < bv->size(); j++) {
        ones_per_pass2 += (*bv)[j];
      }
      pass2_ones.push_back(ones_per_pass2);
      pass2_pointer.push_back(pointers);
      pass2_offset.push_back(offsets);
      pass2_max_pointer.push_back(max_pointer);
    }
    return 0;
  }

  int32_t init_extended(std::vector<input_type>& text) {
    static constexpr uint128_t kPrime = 2305843009213693951ULL;
    /// The number of characters a block tree with s top-level blocks and arity
    /// of strictly tau would exceed over the text size
    int64_t added_padding = 0;
    /// The height of the tree
    int64_t tree_max_height = 0;
    /// The size of the largest blocks (i.e. the top level blocks)
    int64_t max_blk_size = 0;
    /// For each level (starting at the top) contains the text start indices of
    /// each block on the level
    std::vector<std::vector<int64_t>> blk_lvl;
    /// For each level and each block, contains the first block index at which
    /// the content of this block appears
    std::vector<std::vector<size_type>> pass1_pointer;
    /// For each level and each block, contains the offset at which
    /// the content of this block appears in its back-pointed block
    std::vector<std::vector<size_type>> pass1_offset;
    /// For each level contains a bit vector containing a 1 for each block that
    /// is internal and a 0 for each back block
    std::vector<pasta::BitVector*> bv_marked;
    /// For every level and block counts how many back blocks are pointing to
    /// the block
    std::vector<std::vector<size_type>> counter;
    std::vector<size_type> pass2_ones;
    /// The block size for each level, starting at the top level
    std::vector<int64_t> block_size_lvl_temp;
    this->calculate_padding(added_padding,
                            text.size(),
                            tree_max_height,
                            max_blk_size);
    auto is_padded = added_padding > 0 ? 1 : 0;
    /// The current block size starting at the top level
    int64_t block_size = max_blk_size;

    /// The text start indices of each block on the current level
    std::vector<int64_t> block_text_inx;
    for (uint64_t i = 0; i < text.size(); i += block_size) {
      block_text_inx.push_back(i);
    }
    // if the blocks on the current level are already below the max leaf length,
    // we may not divide them further. So the entire block tree just consists of
    // the current level verbatim, no pointers
    if (block_size <= this->max_leaf_length_) {
      auto* bv = new pasta::BitVector(block_text_inx.size(), 1);
      this->block_tree_types_rs_.push_back(
          new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
      auto p0 = new sdsl::int_vector<>(0, 0);
      auto o0 = new sdsl::int_vector<>(0, 0);
      auto& ptr0 = *p0;
      auto& off0 = *o0;
      sdsl::util::bit_compress(ptr0);
      sdsl::util::bit_compress(off0);
      this->block_tree_types_.push_back(bv);
      this->block_tree_pointers_.push_back(p0);
      this->block_tree_offsets_.push_back(o0);
      this->block_size_lvl_.push_back(block_size);
      this->leaf_size = block_size / this->tau_;
      this->leaves_ = std::vector<input_type>(text.begin(), text.end());
      this->compress_leaves();
      return 0;
    }
    while (block_size > this->max_leaf_length_) {
      block_size_lvl_temp.push_back(block_size);
      /// Marks whether a block should be internal or not
      auto* bv = new pasta::BitVector(block_text_inx.size(), false);
      /// If left[i] == 1, then there is an earlier occurrence of
      /// block[i]block[i+1]
      auto left = pasta::BitVector(block_text_inx.size(), false);
      /// If right[i] == 1, then there is an earlier occurrence of
      /// block[i-1]block[i]
      auto right = pasta::BitVector(block_text_inx.size(), false);
      auto pair_size = 2 * block_size;
      // Check, whether the last block's end extends past the end of the text
      auto last_block_padded =
          static_cast<uint64_t>(block_text_inx[block_text_inx.size() - 1] +
                                block_size) != text.size() ?
              1 :
              0;
      // map block pair hashes to the text index of their occurrences
      // collecting duplicates in a vector TODO
      HashMap<MersenneHash<uint8_t>, std::vector<size_type>> pairs(0);
      // map block hashes to their *block* index,
      // collecting duplicates in a vector TODO
      HashMap<MersenneHash<uint8_t>, std::vector<size_type>> blocks =
          HashMap<MersenneHash<uint8_t>, std::vector<size_type>>();
      // iterate through all blocks on the current level, skipping over the last
      // block if it is padded
      for (uint64_t i = 0; i < block_text_inx.size() - last_block_padded; i++) {
        // Hash the current block and insert it into the block hash map
        auto index = block_text_inx[i];
        MersenneRabinKarp<input_type, size_type> rk_block =
            MersenneRabinKarp<input_type, size_type>(text,
                                                     sigma_,
                                                     index,
                                                     block_size,
                                                     kPrime);
        MersenneHash<input_type> mh_block =
            MersenneHash<input_type>(text, rk_block.hash_, index, block_size);
        blocks[mh_block].push_back(i);
      }
      // Back pointers, offsets and the incoming-pointer-counters which we need
      // for pruning later
      std::vector<size_type> pointers(block_text_inx.size(), NO_FORMER_OCC);
      std::vector<size_type> offsets(block_text_inx.size(), 0);
      std::vector<size_type> counters(block_text_inx.size(), 0);
      // If a block pair is larger than the whole text, then there is
      // nothing really to do on this level. There cannot be any back pointers
      if (static_cast<uint64_t>(pair_size) > text.size()) {
        block_size /= this->tau_;
        // Start indices of blocks for the next level
        std::vector<int64_t> new_blocks(0);
        for (uint64_t i = 0; i < block_text_inx.size(); i++) {
          // All of the current blocks are internal
          (*bv)[i] = 1;
          // We split all blocks on the current level into tau sub-blocks but
          // only create blocks that start before the end of the text
          for (size_type j = 0; j < this->tau_; j++) {
            if (static_cast<uint64_t>(block_text_inx[i] + (j * block_size)) <
                text.size()) {
              new_blocks.push_back(block_text_inx[i] + (j * block_size));
            }
          }
        }
        /// There are no pointers, offsets or counters on this level
        std::vector<size_type> p(block_text_inx.size(), -1);
        std::vector<size_type> o(block_text_inx.size(), 0);
        std::vector<size_type> c(block_text_inx.size(), 0);
        blk_lvl.push_back(block_text_inx);
        block_text_inx = new_blocks;
        bv_marked.push_back(bv);
        pass1_pointer.push_back(p);
        pass1_offset.push_back(o);
        counter.push_back(c);
        continue;
      }
      // Iterate through the block pairs and add them to the pair hash table
      // along with their index *if they are consecutive*
      for (uint64_t i = 0; i < block_text_inx.size() - 1; i++) {
        if (block_text_inx[i] + block_size == block_text_inx[i + 1] &&
            static_cast<uint64_t>(block_text_inx[i] + pair_size) <=
                text.size()) {
          auto index = block_text_inx[i];
          MersenneRabinKarp<input_type, size_type> rk_pair =
              MersenneRabinKarp<input_type, size_type>(text,
                                                       sigma_,
                                                       index,
                                                       pair_size,
                                                       kPrime);
          MersenneHash<input_type> mh_pair =
              MersenneHash<input_type>(text, rk_pair.hash_, index, pair_size);
          pairs[mh_pair].push_back(i);
        }
      }
      // Find the occurrences of all block pairs' contents
      MersenneRabinKarp<input_type, size_type> rk_pair_sw =
          MersenneRabinKarp<input_type, size_type>(text,
                                                   sigma_,
                                                   0,
                                                   pair_size,
                                                   kPrime);
      // Hash each window in the text of the size of a block pair
      // and see if it corresponds to an actual block pair
      for (uint64_t i = 0; i < text.size() - pair_size; i++) {
        MersenneHash<input_type> mh_sw =
            MersenneHash<input_type>(text, rk_pair_sw.hash_, i, pair_size);
        if (pairs.find(mh_sw) != pairs.end()) {
          // If the current hash corresponds to a hashed block pair,
          // we update for those pairs, that they have an earlier occurrence
          for (auto b : pairs[mh_sw]) {
            if (i != static_cast<uint64_t>(block_text_inx[b])) {
              left[b] = 1;
              right[b + 1] = 1;
            }
          }
          pairs.erase(mh_sw);
        }
        rk_pair_sw.next();
      }
      auto old_block_size = block_size;
      auto new_block_size = block_size / this->tau_;
      std::vector<int64_t> new_blocks(0);
      for (uint64_t i = 0; i < block_text_inx.size(); i++) {
        /// This is true <=> the current block is adjacent
        /// to its predecessor and successor
        bool surrounded =
            (i > 0 && i < block_text_inx.size() - 1) &&
            block_text_inx[i] + old_block_size == block_text_inx[i + 1] &&
            block_text_inx[i - 1] + old_block_size == block_text_inx[i];
        /// marked <=> not internal (ONLY HERE IT SEEMS)
        bool marked = false;
        // if the block is adjacent to its predecessor and successor, then it
        // must have a previous occurrence both as a left and right part of a
        // block pair in order to be a back block.
        // Otherwise, only consecutive neighboring blocks need to be considered.
        if (surrounded) {
          marked = left[i] && right[i];
        } else {
          marked = left[i] || right[i];
        }
        // If either the block is internal or extends past the text end (i.e. is
        // padded), then we create tau child nodes for this block
        if (!(marked) || static_cast<uint64_t>(block_text_inx[i] +
                                               old_block_size) >= text.size()) {
          // This block is internal
          (*bv)[i] = 1;
          for (size_type j = 0; j < this->tau_; j++) {
            if (static_cast<uint64_t>(block_text_inx[i] +
                                      (j * new_block_size)) < text.size()) {
              new_blocks.push_back(block_text_inx[i] + (j * new_block_size));
            }
          }
        }
      }
      MersenneRabinKarp<input_type, size_type> rk_first_occ =
          MersenneRabinKarp<input_type, size_type>(text,
                                                   sigma_,
                                                   block_text_inx[0],
                                                   block_size,
                                                   kPrime);
      // Identify the first occurrence for each block on this level
      for (int64_t i = 0; static_cast<uint64_t>(i) < block_text_inx.size() - 1;
           i++) {
        // This is true <=>
        // This is not the last block,
        // this block is adjacent to the next and
        // the next block is internal
        // We need this, because the hasher overlaps the next block as well.
        bool followed =
            (static_cast<uint64_t>(i) < block_text_inx.size() - 1) &&
            block_text_inx[i] + block_size == block_text_inx[i + 1] &&
            (*bv)[i + 1] == 1;
        // If this block is internal
        if ((*bv)[i] == 1) {
          // If the hasher's current position is currently not at the current
          // block index, move it there
          if (rk_first_occ.init_ != static_cast<uint64_t>(block_text_inx[i])) {
            rk_first_occ.restart(block_text_inx[i]);
          }
          if (followed) {
            // We iterate through every window that starts in this block
            // and ends before the end of the text.
            // j is the offset into the current block
            for (int64_t j = 0; j < block_size &&
                                static_cast<uint64_t>(block_text_inx[i] + j +
                                                      block_size) < text.size();
                 j++) {
              // Hash the window and try to find an earlier occurrence
              MersenneHash<input_type> mh_first_occ =
                  MersenneHash<input_type>(text,
                                           rk_first_occ.hash_,
                                           block_text_inx[i] + j,
                                           block_size);
              if (blocks.find(mh_first_occ) != blocks.end()) {
                for (auto b : blocks[mh_first_occ]) {
                  // The if the current block (b) were i, it would reference
                  // itself. If j > 0 then the occurrence overlaps the block i
                  // + 1. Therefore, in that case b must be a block *past* i+1
                  if (b > i && (j <= 0 || b > i + 1)) {
                    // If all is well, set the back pointer and offsets
                    pointers[b] = i;
                    offsets[b] = j;
                    if ((*bv)[b] == 0) {
                      // If b is a back block, then we have another block
                      // pointing to i
                      counters[i]++;
                      if (j > 0) {
                        // if the offset is greater than 0, the copied area
                        // extends into the next block
                        counters[i + 1]++;
                      }
                    }
                  }
                }
                blocks.erase(mh_first_occ);
              }
              rk_first_occ.next();
            }
          } else {
            // If the next block is not adjacent, we only hash once
            MersenneHash<input_type> mh_first_occ =
                MersenneHash<input_type>(text,
                                         rk_first_occ.hash_,
                                         block_text_inx[i],
                                         block_size);
            if (blocks.find(mh_first_occ) != blocks.end()) {
              for (auto b : blocks[mh_first_occ]) {
                if (b != i) {
                  pointers[b] = i;
                  offsets[b] = 0;
                }
              }
              blocks.erase(mh_first_occ);
            }
          }
        }
      }
      // Add the values calculated on this level
      pass1_pointer.push_back(pointers);
      pass1_offset.push_back(offsets);
      counter.push_back(counters);
      const_size += pointers.size() * sizeof(size_type) * 2;
      blk_lvl.push_back(block_text_inx);
      block_text_inx = new_blocks;
      block_size = new_block_size;
      bv_marked.push_back(bv);
    }
    // By this point, the first pass is done and we have an unpruned block tree
    this->leaf_size = block_size;
    block_size *= this->tau_;
    // Prune the tree. Doing so will replace the pointers of pruned nodes with
    // PRUNED
    pruning_extended(counter,
                     pass1_pointer,
                     pass1_offset,
                     bv_marked,
                     bv_marked);

    std::vector<size_type> ones_per_lvl(bv_marked.size(), 0);
    // count 1s in each lvl;
    for (uint64_t i = 0; i < bv_marked.size(); i++) {
      auto& current_lvl = *bv_marked[i];
      for (uint64_t j = 0; j < bv_marked[i]->size(); j++) {
        if (current_lvl[j]) {
          ones_per_lvl[i]++;
        }
      }
    }

    auto& top_level = *bv_marked[0];
    bool found_back_block =
        top_level.size() != static_cast<uint64_t>(ones_per_lvl[0]) ||
        bv_marked.size() == 1;
    // If there is a back block on the first level, add its values to the tree
    if (found_back_block || !this->CUT_FIRST_LEVELS) {
      this->block_tree_types_.push_back(&top_level);
      this->block_tree_types_rs_.push_back(
          new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(top_level));
      auto p0 = new sdsl::int_vector<>(top_level.size() - ones_per_lvl[0], 0);
      auto o0 = new sdsl::int_vector<>(top_level.size() - ones_per_lvl[0], 0);
      auto& ptr0 = *p0;
      auto& off0 = *o0;
      size_type c = 0;
      for (uint64_t j = 0; j < top_level.size(); j++) {
        if (!top_level[j]) {
          ptr0[c] = pass1_pointer[0][j];
          off0[c] = pass1_offset[0][j];
          c++;
        }
      }
      sdsl::util::bit_compress(ptr0);
      sdsl::util::bit_compress(off0);
      this->block_tree_pointers_.push_back(p0);
      this->block_tree_offsets_.push_back(o0);
      this->block_size_lvl_.push_back(block_size_lvl_temp[0]);
    } else {
      delete bv_marked[0];
    }

    for (uint64_t i = 1; i < bv_marked.size(); i++) {
      // If the previous level is padded, we need to be careful, since the last
      // block possibly does not generate exactly tau children
      size_type new_size = (ones_per_lvl[i - 1] - is_padded) * this->tau_;
      auto last_block_parent = blk_lvl[i - 1][blk_lvl[i - 1].size() - 1];
      auto lvl_block_size = block_size_lvl_temp[i];
      // Determine the number of children the last block generated
      if (is_padded) {
        for (uint64_t j = 0; j < static_cast<uint64_t>(this->tau_); j++) {
          if (last_block_parent + j * lvl_block_size < text.size()) {
            new_size++;
          }
        }
      }
      // Check if we have found a back block on the current level
      found_back_block |= new_size != ones_per_lvl[i];
      // If there is a back block, we add this level's data to the tree
      if (found_back_block || !this->CUT_FIRST_LEVELS ||
          i == bv_marked.size() - 1) {
        // is_internal
        auto bit_vector = new pasta::BitVector(new_size, 0);
        auto& bv_ref = *bit_vector;
        auto p = new sdsl::int_vector<>(bv_ref.size() - ones_per_lvl[i], 0);
        auto o = new sdsl::int_vector<>(bv_ref.size() - ones_per_lvl[i], 0);
        auto& ptr = *p;
        auto& off = *o;
        // Maps block index => number of pruned blocks before this block
        HashMap<size_type, size_type> blocks_skipped;
        auto& lvl_pass1 = *bv_marked[i];
        // Number of non-pruned blocks so far
        size_type c = 0;
        //
        size_type c_u = 0;
        for (uint64_t j = 0; j < lvl_pass1.size(); j++) {
          blocks_skipped[j] = j - c;
          // If the current block is not pruned, add it to the new tree
          if (pass1_pointer[i][j] != PRUNED) {
            // Add it to the is_internal bit vector
            bv_ref[c] = (bool)lvl_pass1[j];
            // If it is a back block, add its pointer and offset
            if (!lvl_pass1[j]) {
              // We need to ignore the pruned blocks
              ptr[c_u] =
                  pass1_pointer[i][j] - blocks_skipped[pass1_pointer[i][j]];
              off[c_u] = pass1_offset[i][j];
              c_u++;
            }
            c++;
          }
        }
        // Add the new data to the tree
        this->block_tree_types_.push_back(&bv_ref);
        this->block_tree_types_rs_.push_back(
            new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(bv_ref));
        sdsl::util::bit_compress(ptr);
        sdsl::util::bit_compress(off);
        this->block_tree_pointers_.push_back(p);
        this->block_tree_offsets_.push_back(o);
        this->block_size_lvl_.push_back(block_size_lvl_temp[i]);
      }
      // Delete the old bitvec since we don't need it anymore.
      // If this is the last level, we still need the bv for constructing the
      // leaves
      if (i < bv_marked.size() - 1) {
        delete bv_marked[i];
      }
    }

    // Construct the leaf string
    int64_t leaf_count = 0;
    auto& last_level = (*bv_marked[bv_marked.size() - 1]);
    for (uint64_t i = 0; i < last_level.size(); i++) {
      if (last_level[i] == 1) {
        // For every leaf on the last level, we have tau leaf blocks
        leaf_count += this->tau_;
        // Iterate through all characters in this child and add them to the leaf
        // string
        for (uint64_t j = 0;
             j < static_cast<uint64_t>(this->leaf_size * this->tau_);
             j++) {
          if (static_cast<uint64_t>(blk_lvl[blk_lvl.size() - 1][i] + j) <
              text.size()) {
            this->leaves_.push_back(text[blk_lvl[blk_lvl.size() - 1][i] + j]);
          }
        }
      }
    }
    delete &last_level;
    this->amount_of_leaves = leaf_count;
    this->compress_leaves();
    return 0;
  }

  int32_t init_simple(std::vector<input_type>& text) {
    static constexpr uint128_t kPrime = 2305843009213693951ULL;
    int64_t added_padding = 0;
    int64_t tree_max_height = 0;
    int64_t max_blk_size = 0;
    std::vector<std::vector<int64_t>> blk_lvl;
    std::vector<std::vector<size_type>> pass1_pointer;
    std::vector<std::vector<size_type>> pass1_offset;
    std::vector<pasta::BitVector*> bv_pass_1;
    std::vector<pasta::BitVector*> bv_pass_2;
    std::vector<std::vector<size_type>> pass2_pointer;
    std::vector<std::vector<size_type>> pass2_offset;
    std::vector<size_type> pass2_max_pointer;
    std::vector<size_type> pass2_max_offset;
    std::vector<size_type> pass2_ones;
    std::vector<int64_t> block_size_lvl_temp;
    this->calculate_padding(added_padding,
                            text.size(),
                            tree_max_height,
                            max_blk_size);
    auto is_padded = added_padding > 0 ? 1 : 0;
    int64_t block_size = max_blk_size;
    std::vector<int64_t> block_text_inx;
    for (uint64_t i = 0; i < text.size(); i += block_size) {
      block_text_inx.push_back(i);
    }
    if (block_size <= this->max_leaf_length_) {
      auto* bv = new pasta::BitVector(block_text_inx.size(), 1);
      this->block_tree_types_rs_.push_back(
          new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
      auto p0 = new sdsl::int_vector<>(0, 0);
      auto o0 = new sdsl::int_vector<>(0, 0);
      auto& ptr0 = *p0;
      auto& off0 = *o0;
      sdsl::util::bit_compress(ptr0);
      sdsl::util::bit_compress(off0);
      this->block_tree_types_.push_back(bv);
      this->block_tree_pointers_.push_back(p0);
      this->block_tree_offsets_.push_back(o0);
      this->block_size_lvl_.push_back(block_size);
      this->leaf_size = block_size / this->tau_;
      this->leaves_ = std::vector<input_type>(text.begin(), text.end());
      this->compress_leaves();
      return 0;
    }
    bool found_back_block = this->max_leaf_length_ * this->tau_ >= block_size;
    while (block_size > this->max_leaf_length_) {
      block_size_lvl_temp.push_back(block_size);
      auto* bv = new pasta::BitVector(block_text_inx.size(), false);
      auto left = pasta::BitVector(block_text_inx.size(), false);
      auto right = pasta::BitVector(block_text_inx.size(), false);
      auto pair_size = 2 * block_size;
      auto last_block_padded =
          static_cast<uint64_t>(block_text_inx[block_text_inx.size() - 1] +
                                block_size) != text.size() ?
              1 :
              0;
      HashMap<MersenneHash<uint8_t>, std::vector<size_type>> pairs(0);
      HashMap<MersenneHash<uint8_t>, std::vector<size_type>> blocks =
          HashMap<MersenneHash<uint8_t>, std::vector<size_type>>();
      for (uint64_t i = 0; i < block_text_inx.size() - last_block_padded; i++) {
        auto index = block_text_inx[i];
        MersenneRabinKarp<input_type, size_type> rk_block =
            MersenneRabinKarp<input_type, size_type>(text,
                                                     sigma_,
                                                     index,
                                                     block_size,
                                                     kPrime);
        MersenneHash<input_type> mh_block =
            MersenneHash<input_type>(text, rk_block.hash_, index, block_size);
        blocks[mh_block].push_back(i);
      }
      std::vector<size_type> pointers(block_text_inx.size(), -1);
      std::vector<size_type> offsets(block_text_inx.size(), 0);
      if (static_cast<uint64_t>(pair_size) > text.size()) {
        block_size /= this->tau_;
        std::vector<int64_t> new_blocks(0);
        for (uint64_t i = 0; i < block_text_inx.size(); i++) {
          (*bv)[i] = 1;
          for (size_type j = 0; j < this->tau_; j++) {
            if (static_cast<uint64_t>(block_text_inx[i] + (j * block_size)) <
                text.size()) {
              new_blocks.push_back(block_text_inx[i] + (j * block_size));
            }
          }
        }
        std::vector<size_type> p(block_text_inx.size(), -1);
        std::vector<size_type> o(block_text_inx.size(), 0);
        blk_lvl.push_back(block_text_inx);
        block_text_inx = new_blocks;
        bv_pass_1.push_back(bv);
        pass1_pointer.push_back(p);
        pass1_offset.push_back(o);
        continue;
      }
      for (uint64_t i = 0; i < block_text_inx.size() - 1; i++) {
        if (block_text_inx[i] + block_size == block_text_inx[i + 1] &&
            static_cast<uint64_t>(block_text_inx[i] + pair_size) <=
                text.size()) {
          auto index = block_text_inx[i];
          MersenneRabinKarp<input_type, size_type> rk_pair =
              MersenneRabinKarp<input_type, size_type>(text,
                                                       sigma_,
                                                       index,
                                                       pair_size,
                                                       kPrime);
          MersenneHash<input_type> mh_pair =
              MersenneHash<input_type>(text, rk_pair.hash_, index, pair_size);
          pairs[mh_pair].push_back(i);
        }
      }
      // find pairs
      MersenneRabinKarp<input_type, size_type> rk_pair_sw =
          MersenneRabinKarp<input_type, size_type>(text,
                                                   sigma_,
                                                   0,
                                                   pair_size,
                                                   kPrime);
      for (uint64_t i = 0; i < text.size() - pair_size; i++) {
        MersenneHash<input_type> mh_sw =
            MersenneHash<input_type>(text, rk_pair_sw.hash_, i, pair_size);
        if (pairs.find(mh_sw) != pairs.end()) {
          for (auto b : pairs[mh_sw]) {
            if (i != static_cast<uint64_t>(block_text_inx[b])) {
              left[b] = 1;
              right[b + 1] = 1;
            }
          }
          pairs.erase(mh_sw);
        }
        rk_pair_sw.next();
      }
      auto old_block_size = block_size;
      auto new_block_size = block_size / this->tau_;
      std::vector<int64_t> new_blocks(0);
      for (uint64_t i = 0; i < block_text_inx.size(); i++) {
        bool surrounded =
            (i > 0 && i < block_text_inx.size() - 1) &&
            block_text_inx[i] + old_block_size == block_text_inx[i + 1] &&
            block_text_inx[i - 1] + old_block_size == block_text_inx[i];
        bool marked = false;
        if (surrounded) {
          marked = left[i] && right[i];
        } else {
          marked = left[i] || right[i];
        }
        if (!(marked) || static_cast<uint64_t>(block_text_inx[i] +
                                               old_block_size) >= text.size()) {
          (*bv)[i] = 1;
          for (size_type j = 0; j < this->tau_; j++) {
            if (static_cast<uint64_t>(block_text_inx[i] +
                                      (j * new_block_size)) < text.size()) {
              new_blocks.push_back(block_text_inx[i] + (j * new_block_size));
            }
          }
        }
      }
      for (uint64_t i = 0; i < block_text_inx.size() - 1; i++) {
        MersenneRabinKarp<input_type, size_type> rk_first_occ =
            MersenneRabinKarp<input_type, size_type>(text,
                                                     sigma_,
                                                     block_text_inx[i],
                                                     block_size,
                                                     kPrime);
        bool followed =
            (i < block_text_inx.size() - 1) &&
            block_text_inx[i] + block_size == block_text_inx[i + 1] &&
            (*bv)[i + 1] == 1;
        if ((*bv)[i] == 1) {
          if (followed) {
            for (uint64_t j = 0;
                 j < static_cast<uint64_t>(block_size) &&
                 block_text_inx[i] + j + block_size < text.size();
                 j++) {
              MersenneHash<input_type> mh_first_occ =
                  MersenneHash<input_type>(text,
                                           rk_first_occ.hash_,
                                           block_text_inx[i] + j,
                                           block_size);
              if (blocks.find(mh_first_occ) != blocks.end()) {
                for (auto b : blocks[mh_first_occ]) {
                  if (static_cast<uint64_t>(b) != i) {
                    pointers[b] = i;
                    offsets[b] = j;
                  }
                }
                blocks.erase(mh_first_occ);
              }
              rk_first_occ.next();
            }
          } else {
            MersenneHash<input_type> mh_first_occ =
                MersenneHash<input_type>(text,
                                         rk_first_occ.hash_,
                                         block_text_inx[i],
                                         block_size);
            if (blocks.find(mh_first_occ) != blocks.end()) {
              for (auto b : blocks[mh_first_occ]) {
                if (static_cast<uint64_t>(b) != i) {
                  pointers[b] = i;
                  offsets[b] = 0;
                }
              }
              blocks.erase(mh_first_occ);
            }
          }
        }
      }
      pass1_pointer.push_back(pointers);
      pass1_offset.push_back(offsets);
      const_size += pointers.size() * sizeof(size_type) * 2;
      blk_lvl.push_back(block_text_inx);
      block_text_inx = new_blocks;
      block_size = new_block_size;
      bv_pass_1.push_back(bv);
    }
    this->leaf_size = block_size;
    block_size *= this->tau_;
    pruning_simple(bv_pass_1,
                   blk_lvl,
                   bv_pass_2,
                   pass1_pointer,
                   pass1_offset,
                   pass2_pointer,
                   pass2_offset,
                   pass2_max_pointer,
                   pass2_max_offset,
                   pass2_ones,
                   block_size);
    auto size = pass2_pointer[pass2_pointer.size() - 1].size();
    found_back_block |= size != 0;
    if (found_back_block || !this->CUT_FIRST_LEVELS) {
      this->block_tree_types_.push_back(bv_pass_2[bv_pass_2.size() - 1]);
      this->block_tree_types_rs_.push_back(
          new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(
              *bv_pass_2[bv_pass_2.size() - 1]));

      auto p1 = new sdsl::int_vector<>(
          size,
          0,
          (8 * sizeof(size_type)) -
              this->leading_zeros(
                  pass2_max_pointer[pass2_max_pointer.size() - 1]));
      auto o1 = new sdsl::int_vector<>(
          size,
          0,
          (8 * sizeof(size_type)) -
              this->leading_zeros(
                  pass2_max_offset[pass2_max_offset.size() - 1]));
      for (uint64_t i = 0; i < pass2_pointer[pass2_pointer.size() - 1].size();
           i++) {
        (*p1)[i] = pass2_pointer[pass2_pointer.size() - 1][size - 1 - i];
        (*o1)[i] = pass2_offset[pass2_offset.size() - 1][size - 1 - i];
      }
      this->block_tree_pointers_.push_back(p1);
      this->block_tree_offsets_.push_back(o1);
      this->block_size_lvl_.push_back(block_size_lvl_temp[0]);
    } else {
      delete bv_pass_2[bv_pass_2.size() - 1];
    }

    size_type level = 0;
    for (size_type i = bv_pass_2.size() - 2; i >= 0; i--) {
      level++;
      auto pass1_i = bv_pass_2.size() - 1 - i;
      auto pass1_parent = pass1_i - 1;
      auto new_size = this->tau_ * (pass2_ones[i + 1] - is_padded);
      // the last block only spawns block that contain text
      auto last_block_parent =
          blk_lvl[pass1_parent][blk_lvl[pass1_parent].size() - 1];
      auto lvl_block_size = block_size_lvl_temp[level];
      if (is_padded) {
        for (size_type j = 0; j < this->tau_; j++) {
          if (static_cast<uint64_t>(last_block_parent + j * lvl_block_size) <
              text.size()) {
            new_size++;
          }
        }
      }
      found_back_block |= new_size != pass2_ones[i];
      if (found_back_block || !this->CUT_FIRST_LEVELS) {
        auto* bit_vector = new pasta::BitVector(new_size, 0);
        auto pointer = std::vector<size_type>();
        auto offset = std::vector<size_type>();
        size_type pointer_saved = 0;
        size_type pointer_skipped = 0;
        HashMap<size_type, size_type> blocks_skipped;
        size_type skip = 0;
        size_type replace = 0;
        for (uint64_t j = 0; j < bv_pass_1[pass1_i - 1]->size(); j++) {
          if ((*bv_pass_1[pass1_i - 1])[j] == 1) {
            if ((*bv_pass_2[i + 1])[j] == 1) {
              for (size_type k = 0;
                   k < this->tau_ && replace * this->tau_ + k < new_size;
                   k++) {
                bool x = (*bv_pass_2[i])[(j - skip) * this->tau_ + k];
                auto skipper = pointer_skipped + pointer_saved;
                (*bit_vector)[replace * this->tau_ + k] = x;
                if (x == 0) {
                  size_type z =
                      pass2_pointer[i][pass2_pointer[i].size() - 1 - skipper];
                  size_type y = blocks_skipped[z];
                  pointer.push_back(y);
                  offset.push_back(
                      pass2_offset[i][pass2_offset[i].size() - 1 - skipper]);
                  pointer_saved++;
                }
                blocks_skipped[(j - skip) * this->tau_ + k] =
                    replace * this->tau_ + k;
              }

              replace++;
            } else {
              pointer_skipped += this->tau_;
            }
          } else {
            skip++;
          }
        }
        auto p = new sdsl::int_vector<>(
            pointer.size(),
            0,
            (8 * sizeof(size_type)) -
                this->leading_zeros(pass2_max_pointer[i]));
        auto o = new sdsl::int_vector<>(
            pointer.size(),
            0,
            (8 * sizeof(size_type)) - this->leading_zeros(pass2_max_offset[i]));
        for (uint64_t j = 0; j < pointer.size(); j++) {
          (*p)[j] = pointer[j];
          (*o)[j] = offset[j];
        }
        this->block_size_lvl_.push_back(block_size_lvl_temp[level]);
        this->block_tree_pointers_.push_back(p);
        this->block_tree_offsets_.push_back(o);
        this->block_tree_types_.push_back(bit_vector);
        this->block_tree_types_rs_.push_back(
            new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(
                *bit_vector));
      }
    }

    int64_t leaf_count = 0;
    for (uint64_t i = 0; i < bv_pass_2[0]->size(); i++) {
      if ((*bv_pass_2[0])[i] == 1) {
        leaf_count += this->tau_;
        for (int j = 0; j < this->leaf_size * this->tau_; j++) {
          if (static_cast<uint64_t>(blk_lvl[blk_lvl.size() - 1][i] + j) <
              text.size()) {
            this->leaves_.push_back(text[blk_lvl[blk_lvl.size() - 1][i] + j]);
          }
        }
      }
    }
    this->amount_of_leaves = leaf_count;
    this->compress_leaves();
    for (auto bv : bv_pass_1) {
      delete bv;
    }
    return 0;
  };

  BlockTreeFP(std::vector<input_type>& text,
              size_type tau,
              size_type max_leaf_length,
              size_type s,
              size_type sigma,
              bool cut_first_levels,
              bool extended_prune) {
    sigma_ = sigma;
    this->CUT_FIRST_LEVELS = cut_first_levels;
    this->map_unique_chars(text);
    this->tau_ = tau;
    this->max_leaf_length_ = max_leaf_length;
    this->s_ = s;
    if (extended_prune) {
      init_extended(text);
    } else {
      init_simple(text);
    }
  };

private:
  // magic number to indicate that a block is pruned
  const int PRUNED = -2;
  // magic number to indicate that a block has no occurrences to its left side
  const int NO_FORMER_OCC = -1;
};

template <typename input_type, typename size_type>
auto* make_block_tree_fp(std::vector<input_type>& input,
                         size_type const tau,
                         size_type const max_leaf_length) {
  return new BlockTreeFP<input_type, size_type>(input,
                                                tau,
                                                max_leaf_length,
                                                1,
                                                256,
                                                true,
                                                true);
}

} // namespace pasta

/******************************************************************************/
