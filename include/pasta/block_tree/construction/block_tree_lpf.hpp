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
#include "pasta/block_tree/utils/lpf_array.hpp"

template<typename input_type, typename size_type>
class BlockTreeLPF: public BlockTree<input_type,size_type> {
public:
    BlockTreeLPF(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, size_type s, bool mark, bool cut_first_level, bool dp) {
        this->CUT_FIRST_LEVELS = cut_first_level;
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        this->s_ = s;
        std::vector<size_type> lpf(text.size());
        std::vector<size_type> lpf_ptr(text.size());
        lpf_array_stack(text, lpf, lpf_ptr);
        if (dp) init_dp(text, lpf, lpf_ptr, mark);
        else init(text, lpf, lpf_ptr,  mark);
    };
    bool prune_block(std::vector<std::vector<size_type>>& counter, std::vector<std::vector<size_type>>& pointer, std::vector<std::vector<size_type>>& offset, std::vector<pasta::BitVector*>& marked_tree, std::vector<pasta::BitVector*>& pruned_tree, size_type i, size_type j, std::vector<pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>>& ranks) {
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
                marked_children |= prune_block(counter, pointer, offset, marked_tree, pruned_tree, i+1, rank_blk * this->tau_ + k,ranks);
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
		      if (static_cast<uint64_t>(rank_blk * this->tau_ + k) < counter[i + 1].size()) {
                            auto ptr_child = pointer[i + 1][rank_blk * this->tau_ + k];
                            counter[i + 1][ptr_child]--;
                            if (offset[i + 1][rank_blk * this->tau_ + k] > 0) {
                                counter[i + 1][ptr_child + 1]--;
                            }
                            pointer[i+1][rank_blk * this->tau_ + k] = PRUNED;
                        }

                    }
                }

            }
        }
        return marked_children || counter[i][j] > 0 || pointer[i][j] == NO_FORMER_OCC;
    }
    int32_t pruning_extended(std::vector<std::vector<size_type>>& counter, std::vector<std::vector<size_type>>& pointer, std::vector<std::vector<size_type>>& offset, std::vector<pasta::BitVector*>& marked_tree, std::vector<pasta::BitVector*>& pruned_tree) {
        std::vector<pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>> ranks;
        for (auto bv: marked_tree) {
            ranks.emplace_back(*bv);
        }
        auto& top_lvl = *marked_tree[0];
        for (size_type j = top_lvl.size() - 1; j >= 0; j--) {
            prune_block(counter, pointer, offset, marked_tree, pruned_tree, 0, j,ranks);
        }

        return 0;
    }
    int32_t init_dp(std::vector<input_type>& text, std::vector<size_type>& lpf, std::vector<size_type>& prevOcc, bool mark) {
        // bv_marked contains all the marked
        std::vector<pasta::BitVector*> bv_marked;
        std::vector<pasta::BitVector*> bv_pruned;
        std::vector<std::vector<size_type>> counter;
        std::vector<std::vector<size_type>> pass2_pointer;
        std::vector<std::vector<size_type>> pass2_offset;
        std::vector<int64_t> pass2_ones;
        std::vector<std::vector<int64_t>> blk_lvl;
        std::vector<int64_t> block_size_lvl_temp;
        int64_t added_padding = 0;
        int64_t tree_max_height = 0;
        int64_t max_blk_size = 0;
        this->calculate_padding(added_padding, text.size(), tree_max_height, max_blk_size);
        int64_t is_padded = (added_padding > 0) ? 1 : 0;
        int64_t block_size = max_blk_size;
        std::vector<int64_t> block_text_inx;
        for (uint64_t i = 0; i < text.size(); i+= block_size) {
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
            return 0;
        }
        for (uint64_t i = 0; i < prevOcc.size(); i++) {
            size_type p = prevOcc[i];
            if (lpf[i] >= block_size && (lpf[p] >= block_size || lpf[i] <= lpf[p])) {
                prevOcc[i] = prevOcc[p];
            }
        }
        // Firstly we build the theory structure or even just generate an all internal blocks tree
        while (block_size > this->max_leaf_length_) {
            // if marking enabled we initialize everything to 0 and mark else is init_simple to 1
            auto lvl_bv = new pasta::BitVector(block_text_inx.size(),!mark);
            auto& bv = *lvl_bv;
            if (mark) {
                mark_blocks(bv, lpf, block_text_inx, block_size);
            }
            auto pointers = std::vector<size_type>(block_text_inx.size(), -1);
            auto offsets = std::vector<size_type>(block_text_inx.size(), -1);
            auto counter_lvl = std::vector<size_type>(block_text_inx.size(), 0);

            for (uint64_t z = 0; z < block_text_inx.size(); z++) {
                size_type first_ind = block_text_inx[z];
                size_type ind = first_ind;
                if (lpf[ind] >= block_size) {
                    size_type ptr = prevOcc[ind];
                    if (block_size + ptr - 1 < first_ind && lpf[ptr] < block_size) {
                        size_type b = this->find_next_smallest_index_binary_search(ptr, block_text_inx);
                        size_type current_offset = ptr % block_size;
                        if (!bv[z]) {
                            counter_lvl[b]++;
                            if (current_offset > 0) {
                                counter_lvl[b + 1]++;
                            }
                        }
                        pointers[z] = b;
                        offsets[z] = current_offset;
                    }
                }
            }
            block_size_lvl_temp.push_back(block_size);
            blk_lvl.push_back(block_text_inx);
            block_size = block_size/ this->tau_;
            for (uint64_t b = 0; b < block_text_inx.size(); b++) {
	      for (size_type j = 0; j < block_size * this->tau_ && static_cast<uint64_t>(block_text_inx[b] + j) < lpf.size(); j++) {
                    size_type i = block_text_inx[b] + j;
                    size_type p = prevOcc[i];
                    if ((lpf[i] >= block_size && lpf[p] >= block_size) || ((p != -1) && lpf[i] <= lpf[p])) {
                        prevOcc[i] = prevOcc[p];
                    }
                }
            }
            std::vector<int64_t> block_text_inx_new;
            generate_next_level(block_text_inx,block_text_inx_new, bv, text.size(), block_size);
            block_text_inx = block_text_inx_new;
            counter.push_back(counter_lvl);
            pass2_pointer.push_back(pointers);
            pass2_offset.push_back(offsets);
            bv_marked.push_back(&bv);
        }

        this->leaf_size = block_size;
        block_size *= this->tau_;
        pruning_extended(counter, pass2_pointer, pass2_offset, bv_marked, bv_marked);
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
        bool found_back_block = top_level.size() != static_cast<uint64_t>(ones_per_lvl[0]);
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
                    ptr0[c] = pass2_pointer[0][j];
                    off0[c] = pass2_offset[0][j];
                    c++;
                }
            }
            sdsl::util::bit_compress(ptr0);
            sdsl::util::bit_compress(off0);
            this->block_tree_pointers_.push_back(p0);
            this->block_tree_offsets_.push_back(o0);
            this->block_size_lvl_.push_back(block_size_lvl_temp[0]);
        }
        for (uint64_t i = 1; i < bv_marked.size(); i++) {
            size_type new_size = (ones_per_lvl[i - 1] - is_padded) * this->tau_;
            auto last_block_parent = blk_lvl[i - 1][blk_lvl[i - 1].size() - 1];
            auto lvl_block_size = block_size_lvl_temp[i];
            if (is_padded) {
                for (size_type j = 0; j < this->tau_; j++) {
		  if (static_cast<uint64_t>(last_block_parent + j * lvl_block_size) < text.size()) {
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
                    if (pass2_pointer[i][j] != -2) {
                        bv_ref[c] = (bool) lvl_pass1[j];
                        if (!lvl_pass1[j]) {
                            ptr[c_u] = pass2_pointer[i][j] - blocks_skipped[pass2_pointer[i][j]];
                            off[c_u] = pass2_offset[i][j];
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
            }
        }
        int64_t leaf_count = 0;
        auto& last_level = (*bv_marked[bv_marked.size() - 1]);
        for (uint64_t i = 0; i < last_level.size(); i++) {
            if (last_level[i] == 1) {
                leaf_count += this->tau_;
                for (int j = 0; j < this->leaf_size * this->tau_; j++) {
		  if (static_cast<uint64_t>(blk_lvl[blk_lvl.size() -1][i] + j) < text.size()) {
                        this->leaves_.push_back(text[blk_lvl[blk_lvl.size() -1][i] + j]);
                    }
                }
            }
        }
        this->amount_of_leaves = leaf_count;
        return 0;
    };
    int32_t init(std::vector<input_type>& text, std::vector<size_type>& lpf, std::vector<size_type>& prevOcc, bool mark) {
        // bv_marked contains all the marked
        std::vector<pasta::BitVector*> bv_marked;
        std::vector<pasta::BitVector*> bv_pruned;
        std::vector<std::vector<size_type>> counter;
        std::vector<std::vector<size_type>> pass2_pointer;
        std::vector<std::vector<size_type>> pass2_offset;
        std::vector<size_type> pass2_max_pointer;
        std::vector<size_type> pass2_max_offset;
        std::vector<int64_t> pass2_ones;
        std::vector<std::vector<int64_t>> blk_lvl;
        std::vector<int64_t> block_size_lvl_temp;
        int64_t added_padding = 0;
        int64_t tree_max_height = 0;
        int64_t max_blk_size = 0;
        this->calculate_padding(added_padding, text.size(), tree_max_height, max_blk_size);
        auto is_padded = added_padding > 0 ? 1 : 0;
        int64_t block_size = max_blk_size;
        std::vector<int64_t> block_text_inx;

        for (uint64_t i = 0; i < text.size(); i+= block_size) {
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
            return 0;
        }
        for (size_type i = 0; static_cast<uint64_t>(i) < prevOcc.size(); i++) {
            size_type p = prevOcc[i];
            if (lpf[i] >= block_size && (lpf[p] >= block_size || lpf[i] <= lpf[p])) {
                prevOcc[i] = prevOcc[p];
            }
        }

        // Firstly we build the theory structure or even just generate an all internal blocks tree
        while (block_size > this->max_leaf_length_) {
            // if marking enabled we initialize everything to 0 and mark else is init_simple to 1
            auto lvl_bv = new pasta::BitVector(block_text_inx.size(),!mark);
            auto& bv = *lvl_bv;
            if (mark) {
                mark_blocks(bv, lpf, block_text_inx, block_size);
            }
            auto pointers = std::vector<size_type>(block_text_inx.size(), -1);
            auto offsets = std::vector<size_type>(block_text_inx.size(), -1);
            auto counter_lvl = std::vector<size_type>(block_text_inx.size(), 0);

            for (size_type i = 0; static_cast<uint64_t>(i) < block_text_inx.size(); i++) {
                size_type first_ind = block_text_inx[i];
                size_type ind = first_ind;
                [[maybe_unused]] bool has_ptr = false;
                while (lpf[ind] >= block_size) {
                    size_type ptr = prevOcc[ind];
                    if (block_size + prevOcc[ind] - 1 >= first_ind) {
                        ind = prevOcc[ind];
                    } else if (lpf[prevOcc[ind]] >= block_size) {
                        ind = prevOcc[ind];
                    } else {
                        // If we don't mark the blocks and instead build the whole tree we can skip the search
                        // div should be cheap, because we need to % anyway
                        size_type b = this->find_next_smallest_index_binary_search(ptr, block_text_inx);
                        size_type current_offset = ptr % block_size;
                        if (!bv[i]) {
                            counter_lvl[b]++;
                            if (current_offset > 0) {
                                counter_lvl[b + 1]++;
                            }
                        }
                        pointers[i] = b;
                        offsets[i] = current_offset;
                        break;
                    }
                }

            }

            block_size_lvl_temp.push_back(block_size);
            blk_lvl.push_back(block_text_inx);
            block_size = block_size/ this->tau_;
//            for (size_type b = 0; b < block_text_inx.size(); b++) {
//                for (size_type j = 0; j < block_size * this->tau_ && block_text_inx[b] + j < text.size(); j++) {
//                    size_type i = block_text_inx[b] + j;
//                    size_type p = prevOcc[i];
//                    if ((lpf[i] >= block_size && lpf[p] >= block_size) || ((p != -1) && lpf[i] <= lpf[p])) {
//                        prevOcc[i] = prevOcc[p];
//                    }
//                }
//            }

            std::vector<int64_t> block_text_inx_new;
            generate_next_level(block_text_inx,block_text_inx_new, bv, text.size(), block_size);
            block_text_inx = block_text_inx_new;
            counter.push_back(counter_lvl);
            pass2_pointer.push_back(pointers);
            pass2_offset.push_back(offsets);
            bv_marked.push_back(&bv);
        }

        this->leaf_size = block_size;
        block_size *= this->tau_;
        pruning_extended(counter, pass2_pointer, pass2_offset, bv_marked, bv_marked);
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
        bool found_back_block = top_level.size() != static_cast<uint64_t>(ones_per_lvl[0]);
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
                    ptr0[c] = pass2_pointer[0][j];
                    off0[c] = pass2_offset[0][j];
                    c++;
                }
            }
            sdsl::util::bit_compress(ptr0);
            sdsl::util::bit_compress(off0);
            this->block_tree_pointers_.push_back(p0);
            this->block_tree_offsets_.push_back(o0);
            this->block_size_lvl_.push_back(block_size_lvl_temp[0]);
        }
        for(uint64_t i = 1; i < bv_marked.size(); i++) {
            size_type new_size = (ones_per_lvl[i - 1] - is_padded) * this->tau_;
            auto last_block_parent = blk_lvl[i - 1][blk_lvl[i - 1].size() - 1];
            auto lvl_block_size = block_size_lvl_temp[i];
            if (is_padded) {
                for (size_type j = 0; j < this->tau_; j++) {
		  if (static_cast<uint64_t>(last_block_parent + j * lvl_block_size) < text.size()) {
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
                [[maybe_unused]] size_type pointer_saved = 0;
                [[maybe_unused]] size_type pointer_skipped = 0;
                std::unordered_map<size_type, size_type> blocks_skipped;
                [[maybe_unused]] size_type skip = 0;
                [[maybe_unused]] size_type replace = 0;
                [[maybe_unused]] auto &lvl_above_pass1 = *bv_marked[i - 1];
                auto &lvl_pass1 = *bv_marked[i];
                size_type c = 0;
                size_type c_u = 0;
                for (uint64_t j = 0; j < lvl_pass1.size(); j++) {
                    blocks_skipped[j] = j - c;
                    if (pass2_pointer[i][j] != -2) {
                        bv_ref[c] = (bool) lvl_pass1[j];
                        if (!lvl_pass1[j]) {
                            ptr[c_u] = pass2_pointer[i][j] - blocks_skipped[pass2_pointer[i][j]];
                            off[c_u] = pass2_offset[i][j];
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
            }
        }

        int64_t leaf_count = 0;
        auto& last_level = (*bv_marked[bv_marked.size() - 1]);
        for (uint64_t i = 0; i < last_level.size(); i++) {
            if (last_level[i] == 1) {
                leaf_count += this->tau_;
                for (int j = 0; j < this->leaf_size * this->tau_; j++) {
		  if (static_cast<uint64_t>(blk_lvl[blk_lvl.size() -1][i] + j) < text.size()) {
                        this->leaves_.push_back(text[blk_lvl[blk_lvl.size() -1][i] + j]);
                    }
                }
            }
        }
        this->amount_of_leaves = leaf_count;

        return 0;
    }
    BlockTreeLPF(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, bool mark, bool cut_first_level, bool dp) {
        this->CUT_FIRST_LEVELS = cut_first_level;
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        std::vector<size_type> lpf(text.size());
        std::vector<size_type> lpf_ptr(text.size());
        std::vector<size_type> lz;
        lpf_array_stack(text, lpf, lpf_ptr);
        calculate_lz_factor(this->s_,lpf, lz);
        if (dp) init_dp(text, lpf, lpf_ptr,  mark);
        else init(text, lpf, lpf_ptr, mark);
    };
    BlockTreeLPF(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, std::vector<size_type>& lz, bool mark, bool cut_first_level) {
        this->CUT_FIRST_LEVELS = cut_first_level;
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        this->s_ = lz.size();
        init(text, lpf, lpf_ptr, mark);
    };
    BlockTreeLPF(std::vector<input_type>& text, size_type tau, size_type max_leaf_length,size_type s, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, [[maybe_unused]] std::vector<size_type>& lz, bool mark, bool cut_first_level) {
        this->CUT_FIRST_LEVELS = cut_first_level;
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        this->s_ = s;
        init(text, lpf, lpf_ptr, mark);
    };
    ~BlockTreeLPF() {
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
    size_type generate_next_level(std::vector<int64_t> &old_level, std::vector<int64_t> &new_level, pasta::BitVector& bv, int64_t N, int64_t block_size) {
        for (uint64_t i = 0; i < bv.size(); i++) {
            if (bv[i] == 1) {
                for (size_type j = 0; j < this->tau_ ; j++) {
                    if (old_level[i] + (j * block_size) < N) {
                        new_level.push_back(old_level[i] + (j * block_size));
                    }
                }
            }
        }
        return 0;
    }

    size_type mark_blocks(pasta::BitVector& bv, std::vector<size_type>& lpf, std::vector<int64_t>& block_text_inx, int64_t block_size) {
        bv[0] = 1;
        for (uint64_t i = 1; i < block_text_inx.size() - 1; i++) {
            if (block_text_inx[i - 1] + block_size == block_text_inx[i] && block_text_inx[i] + block_size == block_text_inx[i+1] &&
                (lpf[block_text_inx[i - 1]] < 2 * block_size || lpf[block_text_inx[i]] < 2 * block_size)) {
                bv[i] = 1;
            }
        }
        if (bv.size() > 1) {
            size_type last = bv.size() - 1;
            bv[last] = block_text_inx[last - 1] + block_size == block_text_inx[last] && lpf[block_text_inx[last - 1]] < 2 * block_size;
        }
        return 0;
        return 0;
    }
};

template <typename input_type, typename size_type>
auto* make_block_tree_lpf(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, bool set_s_to_z) {
  std::vector<size_type> lpf(text.size());
  std::vector<size_type> lpf_ptr(text.size());
  std::vector<size_type> lz;
  size_type lzn = 0;
  lpf_array(text, lpf, lpf_ptr);
  calculate_lz_factor(lzn, lpf, lz);
  
  return new BlockTreeLPF<input_type, size_type>(text, tau, max_leaf_length, (set_s_to_z ? lzn : 1), lpf, lpf_ptr, lz, false, true);
}

template <typename input_type, typename size_type>
auto* make_block_tree_lpf_parallel(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, bool set_s_to_z, size_t threads) {
  std::vector<size_type> lpf(text.size());
  std::vector<size_type> lpf_ptr(text.size());
  std::vector<size_type> lz;
  size_type lzn = 0;
  lpf_array_ansv(text, lpf, lpf_ptr, threads);
  calculate_lz_factor(lzn, lpf, lz);
  
  return new BlockTreeLPF<input_type, size_type>(text, tau, max_leaf_length, (set_s_to_z ? lzn : 1), lpf, lpf_ptr, lz, false, true);
}

/******************************************************************************/
