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

#include "pasta/block_tree/block_tree.hpp"
#include "pasta/block_tree/utils/MersenneRabinKarp.hpp"
#include "pasta/block_tree/utils/MersenneHash.hpp"

__extension__ typedef unsigned __int128 uint128_t;

template<typename input_type, typename size_type>
class BlockTreeFP : public BlockTree<input_type, size_type> {
public:
    size_type const_size = 0;
    size_type sigma_ = 0;
    bool prune_block(std::vector<std::vector<size_type>> &counter, std::vector<std::vector<size_type>> &pointer,
                     std::vector<std::vector<size_type>> &offset, std::vector<pasta::BitVector *> &marked_tree,
                     std::vector<pasta::BitVector *> &pruned_tree, size_type i, size_type j,
                     std::vector<pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>> &ranks) {
        // string leaf children can always be pruned
        // fully padded children don't exist and can be ignored/ sanity check
        // assumes short circuit evaluation as compiler behaviour
      if (static_cast<uint64_t>(i) >= marked_tree.size() || static_cast<uint64_t>(j) >= marked_tree[i]->size()) {
            return false;
        }


        bool marked_children = false;
        // we already incremented counters for unmarked blocks during a previous step and now only need to consider marked blocks
        if ((*marked_tree[i])[j] == 1) {
            // inverse postorder dfs
            // prune_block returns true if a marked block stays marked false otherwise
            size_type rank_blk = ranks[i].rank1(j);
            for (size_type k = this->tau_ - 1; k >= 0; k--) {
                marked_children |= prune_block(counter, pointer, offset, marked_tree, pruned_tree, i + 1,
                                               rank_blk * this->tau_ + k, ranks);
            }
            // conditions to be pruned are no marked children, no pointers pointing to me and a former occurrence in S
            if (!marked_children && counter[i][j] == 0 && pointer[i][j] != NO_FORMER_OCC) {

                (*marked_tree[i])[j] = 0;
                counter[i][pointer[i][j]]++;
                if (offset[i][j] > 0) {
                    counter[i][pointer[i][j] + 1]++;
                }
                if (static_cast<uint64_t>(i + 1) < counter.size()) {
                    // remove all of its children by decrementing counters and marking them as PRUNED
                    for (size_type k = this->tau_ - 1; k >= 0; k--) {
		      if (static_cast<uint64_t>(rank_blk * this->tau_) + k < counter[i + 1].size()) {
                            auto ptr_child = pointer[i + 1][rank_blk * this->tau_ + k];
                            counter[i + 1][ptr_child]--;
                            if (offset[i + 1][rank_blk * this->tau_ + k] > 0) {
                                counter[i + 1][ptr_child + 1]--;
                            }
                            pointer[i + 1][rank_blk * this->tau_ + k] = PRUNED;
                        }

                    }
                }

            }
        }
        return marked_children || counter[i][j] > 0 || pointer[i][j] == NO_FORMER_OCC;
    };

    int32_t pruning_extended(std::vector<std::vector<size_type>> &counter, std::vector<std::vector<size_type>> &pointer,
                             std::vector<std::vector<size_type>> &offset, std::vector<pasta::BitVector *> &marked_tree,
                             std::vector<pasta::BitVector *> &pruned_tree) {
        std::vector<pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>> ranks;
        for (auto bv: marked_tree) {
            ranks.push_back(pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
        }
        auto &top_lvl = *marked_tree[0];
        for (size_type j = top_lvl.size() - 1; j >= 0; j--) {
            prune_block(counter, pointer, offset, marked_tree, pruned_tree, 0, j, ranks);
        }
        return 0;
    }

    int32_t pruning_simple(std::vector<pasta::BitVector *> &first_pass_bv, std::vector<std::vector<int64_t>> &blk_lvl,
                           std::vector<pasta::BitVector *> &bv_pass_2,
                           std::vector<std::vector<size_type>> &pass1_pointer,
                           std::vector<std::vector<size_type>> &pass1_offset,
                           std::vector<std::vector<size_type>> &pass2_pointer,
                           std::vector<std::vector<size_type>> &pass2_offset,
                           std::vector<size_type> &pass2_max_pointer, std::vector<size_type> &pass2_max_offset,
                           std::vector<size_type> &pass2_ones,
                           int64_t &block_size) {

        for (int64_t i = first_pass_bv.size() - 1; i >= 0; i--) {
            auto *bv = new pasta::BitVector(blk_lvl[i].size(), 0);
            size_type marked_counter = 0;
            if (static_cast<uint64_t>(i) != first_pass_bv.size() - 1) {
                for (uint64_t j = 0; j < bv->size(); j++) {
                    if ((*first_pass_bv[i])[j] == 1) {
                        for (size_type k = 0; k < this->tau_; k++) {
                            if ((*bv_pass_2[bv_pass_2.size() - 1])[marked_counter * this->tau_ + k] == 1) {
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

    int32_t init_extended(std::vector<input_type> &text) {
        static constexpr uint128_t kPrime = 2305843009213693951ULL;
        int64_t added_padding = 0;
        int64_t tree_max_height = 0;
        int64_t max_blk_size = 0;
        std::vector<std::vector<int64_t>> blk_lvl;
        std::vector<std::vector<size_type>> pass1_pointer;
        std::vector<std::vector<size_type>> pass1_offset;
        std::vector<pasta::BitVector*> bv_marked;
        std::vector<std::vector<size_type>> counter;
        std::vector<size_type> pass2_ones;
        std::vector<int64_t> block_size_lvl_temp;
        this->calculate_padding(added_padding, text.size(), tree_max_height, max_blk_size);
        auto is_padded = added_padding > 0 ? 1 : 0;
        int64_t block_size = max_blk_size;

        std::vector<int64_t> block_text_inx;
        for (uint64_t i = 0; i < text.size(); i += block_size) {
            block_text_inx.push_back(i);
        }
        if (block_size <= this->max_leaf_length_) {
            auto *bv = new pasta::BitVector(block_text_inx.size(), 1);
            this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
            auto p0 = new sdsl::int_vector<>(0, 0);
            auto o0 = new sdsl::int_vector<>(0, 0);
            auto &ptr0 = *p0;
            auto &off0 = *o0;
            sdsl::util::bit_compress(ptr0);
            sdsl::util::bit_compress(off0);
            this->block_tree_types_.push_back(bv);
            this->block_tree_pointers_.push_back(p0);
            this->block_tree_offsets_.push_back(o0);
            this->block_size_lvl_.push_back(block_size);
            this->leaf_size = block_size/this->tau_;
            this->leaves_ = std::vector<input_type>(text.begin(), text.end());
	    this->compress_leaves();
            return 0;
        }
        while (block_size > this->max_leaf_length_) {
            block_size_lvl_temp.push_back(block_size);
            auto *bv = new pasta::BitVector(block_text_inx.size(), false);
            auto left = pasta::BitVector(block_text_inx.size(), false);
            auto right = pasta::BitVector(block_text_inx.size(), false);
            auto pair_size = 2 * block_size;
            auto last_block_padded = static_cast<uint64_t>(block_text_inx[block_text_inx.size() - 1] + block_size) != text.size() ? 1 : 0;
            std::unordered_map<MersenneHash<uint8_t>, std::vector<size_type>> pairs(0);
            std::unordered_map<MersenneHash<uint8_t>, std::vector<size_type>> blocks = std::unordered_map<MersenneHash<uint8_t>, std::vector<size_type>>();
            for (uint64_t i = 0; i < block_text_inx.size() - last_block_padded; i++) {
                auto index = block_text_inx[i];
                MersenneRabinKarp<input_type, size_type> rk_block = MersenneRabinKarp<input_type, size_type>(text, sigma_,
                                                                                                             index,
                                                                                                             block_size,
                                                                                                             kPrime);
                MersenneHash<input_type> mh_block = MersenneHash<input_type>(text, rk_block.hash_, index, block_size);
                blocks[mh_block].push_back(i);
            }
            std::vector<size_type> pointers(block_text_inx.size(), -1);
            std::vector<size_type> offsets(block_text_inx.size(), 0);
            std::vector<size_type> counters(block_text_inx.size(), 0);
            if (static_cast<uint64_t>(pair_size) > text.size()) {
                block_size /= this->tau_;
                std::vector<int64_t> new_blocks(0);
                for (uint64_t i = 0; i < block_text_inx.size(); i++) {
                    (*bv)[i] = 1;
                    for (size_type j = 0; j < this->tau_; j++) {
		      if (static_cast<uint64_t>(block_text_inx[i] + (j * block_size)) < text.size()) {
                            new_blocks.push_back(block_text_inx[i] + (j * block_size));
                        }
                    }
                }
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
            for (uint64_t i = 0; i < block_text_inx.size() - 1; i++) {
                if (block_text_inx[i] + block_size == block_text_inx[i + 1] &&
                    static_cast<uint64_t>(block_text_inx[i] + pair_size) <= text.size()) {
                    auto index = block_text_inx[i];
                    MersenneRabinKarp<input_type, size_type> rk_pair = MersenneRabinKarp<input_type, size_type>(text,
                                                                                                                sigma_,
                                                                                                                index,
                                                                                                                pair_size,
                                                                                                                kPrime);
                    MersenneHash<input_type> mh_pair = MersenneHash<input_type>(text, rk_pair.hash_, index, pair_size);
                    pairs[mh_pair].push_back(i);
                }
            }
            // find pairs
            MersenneRabinKarp<input_type, size_type> rk_pair_sw = MersenneRabinKarp<input_type, size_type>(text, sigma_, 0,
                                                                                                           pair_size,
                                                                                                           kPrime);
            for (uint64_t i = 0; i < text.size() - pair_size; i++) {
                MersenneHash<input_type> mh_sw = MersenneHash<input_type>(text, rk_pair_sw.hash_, i, pair_size);
                if (pairs.find(mh_sw) != pairs.end()) {
                    for (auto b: pairs[mh_sw]) {
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
                bool surrounded = (i > 0 && i < block_text_inx.size() - 1) &&
                                  block_text_inx[i] + old_block_size == block_text_inx[i + 1] &&
                                  block_text_inx[i - 1] + old_block_size == block_text_inx[i];
                bool marked = false;
                if (surrounded) {
                    marked = left[i] && right[i];
                } else {
                    marked = left[i] || right[i];
                }
                if (!(marked) || static_cast<uint64_t>(block_text_inx[i] + old_block_size) >= text.size()) {
                    (*bv)[i] = 1;
                    for (size_type j = 0; j < this->tau_; j++) {
		      if (static_cast<uint64_t>(block_text_inx[i] + (j * new_block_size)) < text.size()) {
                            new_blocks.push_back(block_text_inx[i] + (j * new_block_size));
                        }
                    }
                }
            }
            MersenneRabinKarp<input_type, size_type> rk_first_occ = MersenneRabinKarp<input_type, size_type>(
                    text, sigma_, block_text_inx[0], block_size, kPrime);
            for (int64_t i = 0; static_cast<uint64_t>(i) < block_text_inx.size() - 1; i++) {
                bool followed =
		  (static_cast<uint64_t>(i) < block_text_inx.size() - 1) && block_text_inx[i] + block_size == block_text_inx[i + 1] &&
                        (*bv)[i + 1] == 1;
                if ((*bv)[i] == 1) {
		  if (rk_first_occ.init_ != static_cast<uint64_t>(block_text_inx[i])) {
                        rk_first_occ.restart(block_text_inx[i]);
                    }
                    if (followed) {
		      for (int64_t j = 0; j < block_size && static_cast<uint64_t>(block_text_inx[i] + j + block_size) < text.size(); j++) {
                            MersenneHash<input_type> mh_first_occ = MersenneHash<input_type>(text, rk_first_occ.hash_,
                                                                                             block_text_inx[i] + j,
                                                                                             block_size);
                            if (blocks.find(mh_first_occ) != blocks.end()) {
                                for (auto b: blocks[mh_first_occ]) {
                                    // b cant be i and if j>0 then b cant follow on i (j>0) -> b > i + 1 (a -> b <=> not a or b)
                                    if (b > i && (j <= 0 || b > i + 1)) {
                                        pointers[b] = i;
                                        offsets[b] = j;
                                        if ((*bv)[b] == 0) {
                                            counters[i]++;
                                            if (j > 0) {
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
                        MersenneHash<input_type> mh_first_occ = MersenneHash<input_type>(text, rk_first_occ.hash_,
                                                                                         block_text_inx[i], block_size);
                        if (blocks.find(mh_first_occ) != blocks.end()) {
                            for (auto b: blocks[mh_first_occ]) {
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
            pass1_pointer.push_back(pointers);
            pass1_offset.push_back(offsets);
            counter.push_back(counters);
            const_size += pointers.size() * sizeof(size_type) * 2;
            blk_lvl.push_back(block_text_inx);
            block_text_inx = new_blocks;
            block_size = new_block_size;
            bv_marked.push_back(bv);
        }
        this->leaf_size = block_size;
        block_size *= this->tau_;
        pruning_extended(counter, pass1_pointer, pass1_offset, bv_marked, bv_marked);

        std::vector<size_type> ones_per_lvl(bv_marked.size(), 0);
        // count 1s in each lvl;
        for (uint64_t i = 0; i < bv_marked.size(); i++) {
            auto &current_lvl = *bv_marked[i];
            for (uint64_t j = 0; j < bv_marked[i]->size(); j++) {
                if (current_lvl[j]) {
                    ones_per_lvl[i]++;
                }
            }
        }


        auto &top_level = *bv_marked[0];
        bool found_back_block = top_level.size() != static_cast<uint64_t>(ones_per_lvl[0]) || bv_marked.size() == 1;
        if (found_back_block || !this->CUT_FIRST_LEVELS) {
            this->block_tree_types_.push_back(&top_level);
            this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(top_level));
            auto p0 = new sdsl::int_vector<>(top_level.size() - ones_per_lvl[0], 0);
            auto o0 = new sdsl::int_vector<>(top_level.size() - ones_per_lvl[0], 0);
            auto &ptr0 = *p0;
            auto &off0 = *o0;
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

            size_type new_size = (ones_per_lvl[i - 1] - is_padded) * this->tau_;
            auto last_block_parent = blk_lvl[i - 1][blk_lvl[i - 1].size() - 1];
            auto lvl_block_size = block_size_lvl_temp[i];
            if (is_padded) {
	      for (uint64_t j = 0; j < static_cast<uint64_t>(this->tau_); j++) {
                    if (last_block_parent + j * lvl_block_size < text.size()) {
                        new_size++;
                    }
                }
            }
            found_back_block |= new_size != ones_per_lvl[i];
            if (found_back_block || !this->CUT_FIRST_LEVELS) {
                auto bit_vector = new pasta::BitVector(new_size, 0);
                auto &bv_ref = *bit_vector;
                auto p = new sdsl::int_vector<>(bv_ref.size() - ones_per_lvl[i], 0);
                auto o = new sdsl::int_vector<>(bv_ref.size() - ones_per_lvl[i], 0);
                auto &ptr = *p;
                auto &off = *o;
                std::unordered_map<size_type, size_type> blocks_skipped;
                auto &lvl_pass1 = *bv_marked[i];
                size_type c = 0;
                size_type c_u = 0;
                for (uint64_t j = 0; j < lvl_pass1.size(); j++) {
                    blocks_skipped[j] = j - c;
                    if (pass1_pointer[i][j] != -2) {
                        bv_ref[c] = (bool) lvl_pass1[j];
                        if (!lvl_pass1[j]) {
                            ptr[c_u] = pass1_pointer[i][j] - blocks_skipped[pass1_pointer[i][j]];
                            off[c_u] = pass1_offset[i][j];
                            c_u++;
                        }
                        c++;
                    }
                }
                this->block_tree_types_.push_back(&bv_ref);
                this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(bv_ref));
                sdsl::util::bit_compress(ptr);
                sdsl::util::bit_compress(off);
                this->block_tree_pointers_.push_back(p);
                this->block_tree_offsets_.push_back(o);
                this->block_size_lvl_.push_back(block_size_lvl_temp[i]);
            } else {
                delete bv_marked[i];
            }
        }

        int64_t leaf_count = 0;
        auto &last_level = (*bv_marked[bv_marked.size() - 1]);
        for (uint64_t i = 0; i < last_level.size(); i++) {
            if (last_level[i] == 1) {
                leaf_count += this->tau_;
                for (uint64_t j = 0; j < static_cast<uint64_t>(this->leaf_size * this->tau_); j++) {
		  if (static_cast<uint64_t>(blk_lvl[blk_lvl.size() - 1][i] + j) < text.size()) {
                        this->leaves_.push_back(text[blk_lvl[blk_lvl.size() - 1][i] + j]);
                    }
                }
            }
        }
        this->amount_of_leaves = leaf_count;
	this->compress_leaves();
        return 0;
    }

    int32_t init_simple(std::vector<input_type> &text) {
        static constexpr uint128_t kPrime = 2305843009213693951ULL;
        int64_t added_padding = 0;
        int64_t tree_max_height = 0;
        int64_t max_blk_size = 0;
        std::vector<std::vector<int64_t>> blk_lvl;
        std::vector<std::vector<size_type>> pass1_pointer;
        std::vector<std::vector<size_type>> pass1_offset;
        std::vector<pasta::BitVector *> bv_pass_1;
        std::vector<pasta::BitVector *> bv_pass_2;
        std::vector<std::vector<size_type>> pass2_pointer;
        std::vector<std::vector<size_type>> pass2_offset;
        std::vector<size_type> pass2_max_pointer;
        std::vector<size_type> pass2_max_offset;
        std::vector<size_type> pass2_ones;
        std::vector<int64_t> block_size_lvl_temp;
        this->calculate_padding(added_padding, text.size(), tree_max_height, max_blk_size);
        auto is_padded = added_padding > 0 ? 1 : 0;
        int64_t block_size = max_blk_size;
        std::vector<int64_t> block_text_inx;
        for (uint64_t i = 0; i < text.size(); i += block_size) {
            block_text_inx.push_back(i);
        }
        if (block_size <= this->max_leaf_length_) {
            auto *bv = new pasta::BitVector(block_text_inx.size(), 1);
            this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
            auto p0 = new sdsl::int_vector<>(0, 0);
            auto o0 = new sdsl::int_vector<>(0, 0);
            auto &ptr0 = *p0;
            auto &off0 = *o0;
            sdsl::util::bit_compress(ptr0);
            sdsl::util::bit_compress(off0);
            this->block_tree_types_.push_back(bv);
            this->block_tree_pointers_.push_back(p0);
            this->block_tree_offsets_.push_back(o0);
            this->block_size_lvl_.push_back(block_size);
            this->leaf_size = block_size/this->tau_;
            this->leaves_ = std::vector<input_type>(text.begin(), text.end());
	    this->compress_leaves();
            return 0;
        }
        bool found_back_block = this->max_leaf_length_ * this->tau_ >= block_size;
        while (block_size > this->max_leaf_length_) {
            block_size_lvl_temp.push_back(block_size);
            auto *bv = new pasta::BitVector(block_text_inx.size(), false);
            auto left = pasta::BitVector(block_text_inx.size(), false);
            auto right = pasta::BitVector(block_text_inx.size(), false);
            auto pair_size = 2 * block_size;
            auto last_block_padded = static_cast<uint64_t>(block_text_inx[block_text_inx.size() - 1] + block_size) != text.size() ? 1 : 0;
            std::unordered_map<MersenneHash<uint8_t>, std::vector<size_type>> pairs(0);
            std::unordered_map<MersenneHash<uint8_t>, std::vector<size_type>> blocks = std::unordered_map<MersenneHash<uint8_t>, std::vector<size_type>>();
            for (uint64_t i = 0; i < block_text_inx.size() - last_block_padded; i++) {
                auto index = block_text_inx[i];
                MersenneRabinKarp<input_type, size_type> rk_block = MersenneRabinKarp<input_type, size_type>(text, sigma_,
                                                                                                             index,
                                                                                                             block_size,
                                                                                                             kPrime);
                MersenneHash<input_type> mh_block = MersenneHash<input_type>(text, rk_block.hash_, index, block_size);
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
		      if (static_cast<uint64_t>(block_text_inx[i] + (j * block_size)) < text.size()) {
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
                    static_cast<uint64_t>(block_text_inx[i] + pair_size) <= text.size()) {
                    auto index = block_text_inx[i];
                    MersenneRabinKarp<input_type, size_type> rk_pair = MersenneRabinKarp<input_type, size_type>(text,
                                                                                                                sigma_,
                                                                                                                index,
                                                                                                                pair_size,
                                                                                                                kPrime);
                    MersenneHash<input_type> mh_pair = MersenneHash<input_type>(text, rk_pair.hash_, index, pair_size);
                    pairs[mh_pair].push_back(i);
                }
            }
            // find pairs
            MersenneRabinKarp<input_type, size_type> rk_pair_sw = MersenneRabinKarp<input_type, size_type>(text, sigma_, 0,
                                                                                                           pair_size, kPrime);
            for (uint64_t i = 0; i < text.size() - pair_size; i++) {
                MersenneHash<input_type> mh_sw = MersenneHash<input_type>(text, rk_pair_sw.hash_, i, pair_size);
                if (pairs.find(mh_sw) != pairs.end()) {
                    for (auto b: pairs[mh_sw]) {
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
                bool surrounded = (i > 0 && i < block_text_inx.size() - 1) &&
                                  block_text_inx[i] + old_block_size == block_text_inx[i + 1] &&
                                  block_text_inx[i - 1] + old_block_size == block_text_inx[i];
                bool marked = false;
                if (surrounded) {
                    marked = left[i] && right[i];
                } else {
                    marked = left[i] || right[i];
                }
                if (!(marked) || static_cast<uint64_t>(block_text_inx[i] + old_block_size) >= text.size()) {
                    (*bv)[i] = 1;
                    for (size_type j = 0; j < this->tau_; j++) {
		      if (static_cast<uint64_t>(block_text_inx[i] + (j * new_block_size)) < text.size()) {
                            new_blocks.push_back(block_text_inx[i] + (j * new_block_size));
                        }
                    }
                }
            }
            for (uint64_t i = 0; i < block_text_inx.size() - 1; i++) {
                MersenneRabinKarp<input_type, size_type> rk_first_occ = MersenneRabinKarp<input_type, size_type>(text,
                                                                                                                 sigma_,
                                                                                                                 block_text_inx[i],
                                                                                                                 block_size,
                                                                                                                 kPrime);
                bool followed =
                        (i < block_text_inx.size() - 1) && block_text_inx[i] + block_size == block_text_inx[i + 1] &&
                        (*bv)[i + 1] == 1;
                if ((*bv)[i] == 1) {
                    if (followed) {
		      for (uint64_t j = 0; j < static_cast<uint64_t>(block_size) && block_text_inx[i] + j + block_size < text.size(); j++) {
                            MersenneHash<input_type> mh_first_occ = MersenneHash<input_type>(text, rk_first_occ.hash_,
                                                                                             block_text_inx[i] + j,
                                                                                             block_size);
                            if (blocks.find(mh_first_occ) != blocks.end()) {
                                for (auto b: blocks[mh_first_occ]) {
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
                        MersenneHash<input_type> mh_first_occ = MersenneHash<input_type>(text, rk_first_occ.hash_,
                                                                                         block_text_inx[i], block_size);
                        if (blocks.find(mh_first_occ) != blocks.end()) {
                            for (auto b: blocks[mh_first_occ]) {
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
        pruning_simple(bv_pass_1, blk_lvl, bv_pass_2, pass1_pointer, pass1_offset, pass2_pointer,
                       pass2_offset, pass2_max_pointer, pass2_max_offset, pass2_ones, block_size);
        auto size = pass2_pointer[pass2_pointer.size() - 1].size();
        found_back_block |= size != 0;
        if (found_back_block || !this->CUT_FIRST_LEVELS) {
            this->block_tree_types_.push_back(bv_pass_2[bv_pass_2.size() - 1]);
            this->block_tree_types_rs_.push_back(
                    new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv_pass_2[bv_pass_2.size() - 1]));

            auto p1 = new sdsl::int_vector<>(size, 0, (8 * sizeof(size_type)) - this->leading_zeros(
                    pass2_max_pointer[pass2_max_pointer.size() - 1]));
            auto o1 = new sdsl::int_vector<>(size, 0, (8 * sizeof(size_type)) - this->leading_zeros(
                    pass2_max_offset[pass2_max_offset.size() - 1]));
            for (uint64_t i = 0; i < pass2_pointer[pass2_pointer.size() - 1].size(); i++) {
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
            auto last_block_parent = blk_lvl[pass1_parent][blk_lvl[pass1_parent].size() - 1];
            auto lvl_block_size = block_size_lvl_temp[level];
            if (is_padded) {
                for (size_type j = 0; j < this->tau_; j++) {

		  if (static_cast<uint64_t>(last_block_parent + j * lvl_block_size) < text.size()) {
                        new_size++;
                    }
                }
            }
            found_back_block |= new_size != pass2_ones[i];
            if (found_back_block || !this->CUT_FIRST_LEVELS) {
                auto *bit_vector = new pasta::BitVector(new_size, 0);
                auto pointer = std::vector<size_type>();
                auto offset = std::vector<size_type>();
                size_type pointer_saved = 0;
                size_type pointer_skipped = 0;
                std::unordered_map<size_type, size_type> blocks_skipped;
                size_type skip = 0;
                size_type replace = 0;
                for (uint64_t j = 0; j < bv_pass_1[pass1_i - 1]->size(); j++) {
                    if ((*bv_pass_1[pass1_i - 1])[j] == 1) {
                        if ((*bv_pass_2[i + 1])[j] == 1) {
                            for (size_type k = 0; k < this->tau_ && replace * this->tau_ + k < new_size; k++) {
                                bool x = (*bv_pass_2[i])[(j - skip) * this->tau_ + k];
                                auto skipper = pointer_skipped + pointer_saved;
                                (*bit_vector)[replace * this->tau_ + k] = x;
                                if (x == 0) {
                                    size_type z = pass2_pointer[i][pass2_pointer[i].size() - 1 - skipper];
                                    size_type y = blocks_skipped[z];
                                    pointer.push_back(y);
                                    offset.push_back(pass2_offset[i][pass2_offset[i].size() - 1 - skipper]);
                                    pointer_saved++;
                                }
                                blocks_skipped[(j - skip) * this->tau_ + k] = replace * this->tau_ + k;
                            }

                            replace++;
                        } else {
                            pointer_skipped += this->tau_;
                        }
                    } else {
                        skip++;
                    }
                }
                auto p = new sdsl::int_vector<>(pointer.size(), 0,
                                                (8 * sizeof(size_type)) - this->leading_zeros(pass2_max_pointer[i]));
                auto o = new sdsl::int_vector<>(pointer.size(), 0,
                                                (8 * sizeof(size_type)) - this->leading_zeros(pass2_max_offset[i]));
                for (uint64_t j = 0; j < pointer.size(); j++) {
                    (*p)[j] = pointer[j];
                    (*o)[j] = offset[j];
                }
                this->block_size_lvl_.push_back(block_size_lvl_temp[level]);
                this->block_tree_pointers_.push_back(p);
                this->block_tree_offsets_.push_back(o);
                this->block_tree_types_.push_back(bit_vector);
                this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bit_vector));
            }
        }

        int64_t leaf_count = 0;
        for (uint64_t i = 0; i < bv_pass_2[0]->size(); i++) {
            if ((*bv_pass_2[0])[i] == 1) {
                leaf_count += this->tau_;
                for (int j = 0; j < this->leaf_size * this->tau_; j++) {
		  if (static_cast<uint64_t>(blk_lvl[blk_lvl.size() - 1][i] + j) < text.size()) {
                        this->leaves_.push_back(text[blk_lvl[blk_lvl.size() - 1][i] + j]);
                    }
                }
            }
        }
        this->amount_of_leaves = leaf_count;
	this->compress_leaves();
        for (auto bv: bv_pass_1) {
            delete bv;
        }
        return 0;
    };

    BlockTreeFP(std::vector<input_type> &text, size_type tau, size_type max_leaf_length, size_type s, size_type sigma,
                           bool cut_first_levels, bool extended_prune) {
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

    ~BlockTreeFP() {
        for(auto& bt_t: this->block_tree_types_) {
            delete bt_t;
        }
        for(auto& bt_rs: this->block_tree_types_rs_) {
            delete bt_rs;
        }
        for(auto& bt_p: this->block_tree_pointers_) {
            delete bt_p;
        }
        for(auto& bt_o: this->block_tree_offsets_) {
            delete bt_o;
        }

    };
private:
    // magic number to indicate that a block is pruned
    const int PRUNED = -2;
    // magic number to indicate that a block has no occurrences to its left side
    const int NO_FORMER_OCC = -1;
};

template <typename input_type, typename size_type>
auto* make_block_tree_fp(std::vector<input_type>& input, size_type const tau, size_type const max_leaf_length) {
  return new BlockTreeFP<input_type, size_type>(input, tau, max_leaf_length, 1, 256, true, true);
}

/******************************************************************************/
