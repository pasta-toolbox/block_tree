//
// Created by daniel on 23.08.22.
//
#include "bv_blocktree.hpp"
#include <lpfarray.hpp>
#ifndef BLOCK_TREE_BV_BLOCKTREE_LZ_HPP
#define BLOCK_TREE_BV_BLOCKTREE_LZ_HPP
template<typename input_type, typename size_type>
class BV_BlockTree_lpf_pruned: public BV_Block_Tree<input_type,size_type> {
public:
    int32_t pruning_extended(std::vector<std::vector<size_type>>& counter, std::vector<std::vector<size_type>>& pointer, std::vector<std::vector<size_type>>& offset, std::vector<pasta::BitVector*>& marked_tree, std::vector<pasta::BitVector*>& pruned_tree) {
        std::vector<pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>> ranks;
        for (auto bv: marked_tree) {
            ranks.push_back(pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
        }
        auto& lowest_lvl = *marked_tree[marked_tree.size() - 1];
        auto &lowest_lvl_counter = counter[counter.size() - 1];
        auto &lowest_lvl_pointer = pointer[counter.size() - 1];
        auto &lowest_lvl_offset = offset[counter.size() - 1];
        for (size_type i = lowest_lvl.size() - 1; i >= 0; i--) {
            if (lowest_lvl[i]) {
                if (lowest_lvl_counter[i] == 0 && lowest_lvl_pointer[i] != -1) {
                    lowest_lvl_counter[lowest_lvl_pointer[i]]++;
                    if (lowest_lvl_offset[i] >0) lowest_lvl_counter[lowest_lvl_pointer[i] + 1]++;
                    lowest_lvl[i] = 0;
                }
            }
        }
        for (size_type i = marked_tree.size() - 2; i >= 0; i--) {
            auto &current_lvl_bv = *marked_tree[i];
            auto &current_lvl_counter = counter[i];
            auto &current_lvl_pointer = pointer[i];
            auto &current_lvl_offset = offset[i];
            auto &below_lvl_bv = *marked_tree[i + 1];
            auto &below_lvl_counter = counter[i + 1];
            auto &below_lvl_pointer = pointer[i + 1];
            auto &below_lvl_offset = offset[i + 1];
            size_type first_child = ((below_lvl_bv.size() - 1)/ this->tau_) * this->tau_;
            size_type last_child = below_lvl_bv.size() - 1;
            for (size_type j = current_lvl_bv.size() - 1; j >= 0; j--) {
                if (current_lvl_bv[j]) {
                    bool has_marked_children = false;
                    size_type child = last_child;
                    while (child >= first_child) {
                        bool child_not_replaceable = false;
                        if (i + 2 < marked_tree.size() && below_lvl_bv[child]) {
                            auto below_rank = pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(below_lvl_bv);
                            size_type x = below_rank.rank1(child) * this->tau_;
                            auto &grandchild_lvl_bv = *marked_tree[i + 2];
                            for (size_type k = 0; k < this->tau_ && x + k < grandchild_lvl_bv.size(); k++) {
                                child_not_replaceable |= grandchild_lvl_bv[x + k];
                            }
                        }
                        below_lvl_bv[child] = (below_lvl_counter[child] > 0 || below_lvl_pointer[child] == -1) || child_not_replaceable;
                        has_marked_children = has_marked_children || below_lvl_bv[child];
                        child--;
                    }
                    if (current_lvl_counter[j] == 0 && current_lvl_pointer[j] > -1) {
                        if (!has_marked_children) {
                            current_lvl_bv[j] = 0;
                            current_lvl_counter[current_lvl_pointer[j]]++;
                            if (current_lvl_offset[j] > 0) {
                                current_lvl_counter[current_lvl_pointer[j] + 1]++;
                            }
                            size_type pointer_child = last_child;
                            while (pointer_child >= first_child) {
                                below_lvl_counter[below_lvl_pointer[pointer_child]]--;
                                if (below_lvl_offset[pointer_child] > 0) {
                                    below_lvl_counter[below_lvl_pointer[pointer_child] + 1]--;
                                }
                                below_lvl_pointer[pointer_child] = -2;
                                pointer_child--;
                            }

                            //delete pointer_child pointer
                        }
                    }
                    last_child = first_child - 1;
                    first_child -= this->tau_;
                }
            }
        }

        return 0;
    }
    int32_t init(std::vector<input_type>& text, std::vector<size_type>& lpf, std::vector<size_type>& prevOcc, std::vector<size_type>& lz, bool mark) {
        // bv_marked contains all the marked
        std::vector<pasta::BitVector*> bv_marked;
        std::vector<pasta::BitVector*> bv_pruned;
        std::vector<std::vector<size_type>> counter;
        std::vector<std::vector<size_type>> pass2_pointer;
        std::vector<std::vector<size_type>> pass2_offset;
        std::vector<size_type> pass2_max_pointer;
        std::vector<size_type> pass2_max_offset;
        std::vector<size_type> pass2_ones;
        std::vector<std::vector<int64_t>> blk_lvl;

        int64_t added_padding = 0;
        int64_t tree_max_height = 0;
        int64_t max_blk_size = 0;
        this->calculate_padding(added_padding, text.size(), tree_max_height, max_blk_size);
        auto is_padded = added_padding > 0 ? 1 : 0;
        int64_t block_size = max_blk_size;
        std::vector<int64_t> block_text_inx;

        for (int64_t i = 0; i < text.size(); i+= block_size) {
            block_text_inx.push_back(i);
        }
        for(size_type i = 0; i < prevOcc.size(); i++) {
            if (lpf[i] <= lpf[prevOcc[i]] && prevOcc[i] == i - 1) {
                prevOcc[i] = prevOcc[prevOcc[i]];
            }
        }
        // Firstly we build the theory structure or even just generate an all internal blocks tree
        auto t01 = std::chrono::high_resolution_clock::now();
        while (block_size > this->max_leaf_length_) {
            // if marking enabled we initialize everything to 0 and mark else is init to 1
            auto lvl_bv = new pasta::BitVector(block_text_inx.size(),!mark);
            auto& bv = *lvl_bv;
            if (mark) {
                mark_blocks(bv, lz, block_text_inx, block_size);
            }
            auto pointers = std::vector<size_type>(block_text_inx.size(), -1);
            auto offsets = std::vector<size_type>(block_text_inx.size(), -1);
            auto counter_lvl = std::vector<size_type>(block_text_inx.size(), 0);
            for (size_type i = 0; i < block_text_inx.size(); i++) {
                size_type first_ind = block_text_inx[i];
                size_type ind = first_ind;
                bool has_ptr = false;
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
            this->block_size_lvl_.push_back(block_size);
            blk_lvl.push_back(block_text_inx);
            block_size = block_size/ this->tau_;
            std::vector<int64_t> block_text_inx_new;
            generate_next_level(block_text_inx,block_text_inx_new, bv, text.size(), block_size);
            block_text_inx = block_text_inx_new;
            counter.push_back(counter_lvl);
            pass2_pointer.push_back(pointers);
            pass2_offset.push_back(offsets);
            bv_marked.push_back(&bv);
        }
        auto t02 = std::chrono::high_resolution_clock::now();
        this->leaf_size = block_size;
        block_size *= this->tau_;
        pruning_extended(counter, pass2_pointer, pass2_offset, bv_marked);
        auto t03 = std::chrono::high_resolution_clock::now();
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << std::endl;
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t03 - t02).count() << std::endl;
        std::vector<size_type> ones_per_lvl(bv_marked.size(), 0);
        // count 1s in each lvl;
        for (size_type i = 0; i < bv_marked.size(); i++) {
            auto& current_lvl = *bv_marked[i];
            for (size_type j = 0; j < bv_marked[i]->size(); j++) {
                if (current_lvl[j]) {
                    ones_per_lvl[i]++;
                }
                if (current_lvl[j] && pass2_pointer[i][j] == -2) {
                    std::cout << i << ":" << j << std::endl;
                }
            }
        }
        auto neg2c = 0;
        for (size_type i = 0; i < pass2_pointer.size(); i++) {
            for (size_type j = 0; j < bv_marked[i]->size(); j++) {
                neg2c += (-2 == pass2_pointer[i][j]);
            }
            std::cout << neg2c << std::endl;
        }
        std::cout << "Remove " << neg2c << std::endl;
        for (auto v: ones_per_lvl) {
            std::cout << v << std::endl;
        }
        std::cout << "done" << std::endl;
        for (auto v: bv_marked) {
            std::cout << v->size() << std::endl;
        }
        std::cout << "done4" << std::endl;
//        std::cout << *bv_marked[0] << std::endl;


        auto& top_level = *bv_marked[0];
        this->block_tree_types_.push_back(&top_level);
        this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(top_level));
        auto p0 = new sdsl::int_vector<>(top_level.size() - ones_per_lvl[0],0);
        auto o0 = new sdsl::int_vector<>(top_level.size() - ones_per_lvl[0],0);
        auto& ptr0 = *p0;
        auto& off0 = *o0;
        size_type c = 0;
        for (size_type j = 0; j < top_level.size(); j++) {
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
        for(size_type i = 1; i < bv_marked.size(); i++) {
            size_type new_size = (ones_per_lvl[i - 1] - is_padded) * this->tau_;
            auto last_block_parent = blk_lvl[i - 1][blk_lvl[i - 1].size() - 1];
            auto lvl_block_size = this->block_size_lvl_[i];
            if (is_padded) {
                for (size_type j = 0; j < this->tau_; j++) {
                    if (last_block_parent + j * lvl_block_size < text.size()) {
                        new_size++;
                    }
                }
            }

            auto bit_vector = new pasta::BitVector(new_size,0);
            auto& bv_ref = *bit_vector;
            auto p = new sdsl::int_vector<>(bv_ref.size() - ones_per_lvl[i],0);
            auto o = new sdsl::int_vector<>(bv_ref.size() - ones_per_lvl[i],0);
            auto& ptr = *p;
            auto& off = *o;
            size_type pointer_saved = 0;
            size_type pointer_skipped = 0;
            std::unordered_map<size_type, size_type> blocks_skipped;
            size_type skip = 0;
            size_type replace = 0;
            auto& lvl_above_pass1 = *bv_marked[i - 1];
            auto& lvl_pass1 = *bv_marked[i];
            size_type c = 0;
            size_type c_u = 0;
            std::cout << "Made it til here" << std::endl;
            for (size_type j = 0; j < lvl_pass1.size(); j++) {
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
            std::cout << "Made it til here2" << std::endl;
            this->block_tree_types_.push_back(&bv_ref);
            this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(bv_ref));
            sdsl::util::bit_compress(ptr);
            sdsl::util::bit_compress(off);
            this->block_tree_pointers_.push_back(p);
            this->block_tree_offsets_.push_back(o);
        }
        int64_t leaf_count = 0;
        auto& last_level = (*bv_marked[bv_marked.size() - 1]);
        for (size_type i = 0; i < last_level.size(); i++) {
            if (last_level[i] == 1) {
                leaf_count += this->tau_;
                for (int j = 0; j < this->leaf_size * this->tau_; j++) {
                    if (blk_lvl[blk_lvl.size() -1][i] + j < text.size()) {
                        this->leaves_.push_back(text[blk_lvl[blk_lvl.size() -1][i] + j]);
                    }
                }
            }
        }
        this->amount_of_leaves = leaf_count;
//        for (size_type i = bv_marked.size() - 2; i >= 0; i--) {
//            auto bv = pasta::BitVector(blk_lvl[i].size(),0);
//            // mark nodes that have children
//            size_type marked_counter = 0;
//            for (size_type j = 0; j < bv.size(); j++) {
//                auto& lvl_bv = (*bv_marked[i]);
//                if (lvl_bv[j] == 1) {
//
//                    for (size_type k = 0; k < this->tau_; k++) {
//                        auto& lvl_below = (*bv_pruned[bv_pruned.size() - 1]);
//                        if (lvl_below[marked_counter * this->tau_ + k] == 1) {
//                            bv[j] = 1;
//                        }
//                    }
//                    marked_counter++;
//                }
//            }
//            size_type max_pointer = 0;
//            size_type max_offset = 0;
//            auto pointers = std::vector<size_type>();
//            auto offsets = std::vector<size_type>();
//            for (size_type j = blk_lvl[i].size() - 1; j >= 0; j--) {
//                if (bv[j] == 1) {
//                    continue;
//                }
//                size_type first_ind = blk_lvl[i][j];
//                size_type ind = first_ind;
//                bool has_ptr = false;
//                while (lpf[ind] >= block_size) {
//                    size_type ptr = prevOcc[ind];
//                    if (block_size + prevOcc[ind] - 1 >= first_ind) {
//                        ind = prevOcc[ind];
//                    } else if (lpf[prevOcc[ind]] >= block_size) {
//                        ind = prevOcc[ind];
//                    } else {
//                        size_type b = this->find_next_smallest_index_binary_search(ptr, blk_lvl[i]);
//
//                        size_type current_offset = ptr % block_size;
//                        bv[b] = 1;
//                        if (current_offset != 0) {
//                            bv[b + 1] = 1;
//                        }
//                        pointers.push_back(b);
//                        offsets.push_back(current_offset);
//                        if (b > max_pointer) {
//                            max_pointer = b;
//                        }
//                        if (current_offset > max_offset) {
//                            max_offset = current_offset;
//                        }
//                        has_ptr = true;
//                        break;
//                    }
//                }
//                if (!has_ptr) {
//                    bv[j] = 1;
//                }
////                if (j == 340) {
////                    std::cout<< i << " " << (*tree_bv)[j] << std::endl;
////                }
//            }
//
//            block_size *= this->tau_;
//            bv_pruned.push_back(&bv);
////            if (i == bv_marked.size() -1) {
////                std::cout << *tree_bv << std::endl;
////            }
//            size_type ones_per_pass2 = 0;
//            for (size_type j = 0; j < bv.size(); j++) {
//                ones_per_pass2 += bv[j];
//            }
//            pass2_ones.push_back(ones_per_pass2);
//            pass2_pointer.push_back(pointers);
//            pass2_offset.push_back(offsets);
//            pass2_max_offset.push_back(max_offset);
//            pass2_max_pointer.push_back(max_pointer);
//
//        }
        //determine lvl final size
        // prepare first level;
//        auto& lvl = (*bv_pruned[bv_pruned.size() - 1]);
//        this->block_tree_types_.push_back(&lvl);
//        this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(lvl));
//        auto size = pass2_pointer[pass2_pointer.size() -1].size();
//        auto p1 = new sdsl::int_vector<>(size,0,(8 * sizeof(size_type)) -  this->leading_zeros(pass2_max_pointer[pass2_max_pointer.size() -1]));
//        auto o1 = new sdsl::int_vector<>(size,0,(8 * sizeof(size_type)) -  this->leading_zeros(pass2_max_offset[pass2_max_offset.size() -1]));
//        auto& ptr1 = *p1;
//        auto& off1 = *o1;
//        for(size_type i = 0; i < pass2_pointer[pass2_pointer.size() - 1].size(); i++) {
//            ptr1[i] = pass2_pointer[pass2_pointer.size() - 1][size- 1 -i];
//            off1[i] = pass2_offset[pass2_offset.size() - 1][size - 1 -i];
//        }
//        this->block_tree_pointers_.push_back(p1);
//        this->block_tree_offsets_.push_back(o1);
//        for (size_type i = bv_pruned.size() - 2; i >= 0; i--) {
//            auto pass1_i = bv_pruned.size() - 1 - i;
//            auto pass1_parent = pass1_i - 1;
//            auto new_size = this->tau_ * (pass2_ones[i + 1] - is_padded);
//            // the last block only spawns block that contain text
//            auto last_block_parent = blk_lvl[pass1_parent][blk_lvl[pass1_parent].size() - 1];
//            auto lvl_block_size = this->block_size_lvl_[pass1_i];
//            if (is_padded) {
//                for (size_type j = 0; j < this->tau_; j++) {
//                    if (last_block_parent + j * lvl_block_size < text.size()) {
//                        new_size++;
//                    }
//                }
//            }
//
//            auto bit_vector = new pasta::BitVector(new_size,0);
//            auto& bv_ref = *bit_vector;
//            auto pointer = std::vector<size_type>();
//            auto offset = std::vector<size_type>();
//            size_type pointer_saved = 0;
//            size_type pointer_skipped = 0;
//            std::unordered_map<size_type, size_type> blocks_skipped;
//            size_type skip = 0;
//            size_type replace = 0;
//            auto& lvl_below_pass1 = *bv_marked[pass1_i - 1];
//            auto& lvl_below_pass2 = *bv_pruned[i + 1];
//            auto& lvl_pass2 = *bv_pruned[i];
//            for (size_type j = 0; j < lvl_below_pass1.size(); j++) {
//                if (lvl_below_pass1[j] == 1) {
//                    if (lvl_below_pass2[j] == 1) {
//                        for (size_type k = 0; k < this->tau_ && replace * this->tau_ + k < new_size; k++) {
//                            bool x = lvl_pass2[(j - skip) * this->tau_ + k];
//                            auto skipper = pointer_skipped + pointer_saved;
//                            bv_ref[replace * this->tau_ + k] = x;
//                            if (x == 0) {
//                                size_type z = pass2_pointer[i][pass2_pointer[i].size() - 1 - skipper];
//                                size_type y = blocks_skipped[z];
//                                pointer.push_back(y);
//                                offset.push_back(pass2_offset[i][pass2_offset[i].size() - 1 - skipper]);
//                                pointer_saved++;
//                            }
//                            blocks_skipped[(j - skip) * this->tau_ + k] = replace * this->tau_ + k;
//                        }
//
//                        replace++;
//                    } else {
//                        pointer_skipped+= this->tau_;
//                    }
//                } else {
//                    skip++;
//                }
//            }
//            auto p = new sdsl::int_vector<>(pointer.size(),0,(8 * sizeof(size_type)) -  this->leading_zeros(pass2_max_pointer[i]));
//            auto o = new sdsl::int_vector<>(pointer.size(),0,(8 * sizeof(size_type)) -  this->leading_zeros(pass2_max_offset[i]));
//            auto& ptr = *p;
//            auto& off = *o;
//            for(size_type j = 0; j < pointer.size(); j++) {
//                ptr[j] = pointer[j];
//                off[j] = offset[j];
//            }
//            this->block_tree_pointers_.push_back(p);
//            this->block_tree_offsets_.push_back(o);
//            this->block_tree_types_.push_back(bit_vector);
////            std::cout << (*bv_marked[pass1_i - 1]) << std::endl;
////            std::cout << (*bv_pruned[i + 1]) << std::endl;
////            std::cout << *bit_vector << std::endl;
//
//            this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(bv_ref));
//        }
//        int64_t leaf_count = 0;
//        auto& last_level = (*bv_pruned[0]);
//        for (size_type i = 0; i < last_level.size(); i++) {
//            if (last_level[i] == 1) {
//                leaf_count += this->tau_;
//                for (int j = 0; j < this->leaf_size * this->tau_; j++) {
//                    if (blk_lvl[blk_lvl.size() -1][i] + j < text.size()) {
//                        this->leaves_.push_back(text[blk_lvl[blk_lvl.size() -1][i] + j]);
//                    }
//                }
//            }
//        }
//        this->amount_of_leaves = leaf_count;
        // final pass
//        for (auto a: bv_marked) {
//            std::cout << *a << std::endl << std::endl;
//        }
//        std::cout << "alles kllar" << std::endl;
//        for (auto a: bv_pruned) {
//            std::cout << *a << std::endl<< std::endl;
//        }

        return 0;
    };

    BV_BlockTree_lpf_pruned(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, bool mark) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        std::vector<size_type> lpf(text.size());
        std::vector<size_type> lpf_ptr(text.size());
        std::vector<size_type> lz;
        lpf_array(text, lpf, lpf_ptr);
        calculate_lz_factor(this->s_,lpf, lz);
        init(text, lpf, lpf_ptr, lz, mark);
    };
    BV_BlockTree_lpf_pruned(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, std::vector<size_type>& lz, bool mark) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        this->s_ = lz.size();
        init(text, lpf, lpf_ptr, lz, mark);
    };
    BV_BlockTree_lpf_pruned(std::vector<input_type>& text, size_type tau, size_type max_leaf_length,size_type s, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, std::vector<size_type>& lz, bool mark) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        this->s_ = s;
        init(text, lpf, lpf_ptr, lz, mark);
    };
    ~BV_BlockTree_lpf_pruned() {
    };
private:
    size_type generate_next_level(std::vector<int64_t> &old_level, std::vector<int64_t> &new_level, pasta::BitVector& bv, int64_t N, int64_t block_size) {
//        std::cout << block_size << " this size " << " " << old_level.size() << " " << new_level.size() << std::endl;
        for (size_type i = 0; i < bv.size(); i++) {
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

    size_type mark_blocks(pasta::BitVector& bv, std::vector<size_type>& lz, std::vector<int64_t>& block_text_inx, int64_t block_size) {

        size_type j = 0;
        for (size_type i = 0; i < lz.size() - 1; i++) {
            size_type f = lz[i];
            while (j < block_text_inx.size() - 1 && block_text_inx[j + 1] <= f) {
                j++;
            }

            bv[j] = 1;

            if (j > 0 && block_text_inx[j - 1] + block_size == block_text_inx[j]) {
                bv[j - 1] = 1;
            }

            if (j + 1 < bv.size() && block_text_inx[j] + block_size == block_text_inx[j + 1]) {
                bv[j + 1] = 1;
            }
        }
        return 0;
    }
};
#endif //BLOCK_TREE_BV_BLOCKTREE_LZ_HPP

