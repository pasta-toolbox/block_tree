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
    int32_t init([[maybe_unused]] std::vector<input_type>& text, std::vector<size_type>& lpf, std::vector<size_type>& prevOcc, std::vector<size_type>& lz, bool mark) {
        // bv_pass_1 contains all the marked
        std::vector<pasta::BitVector*> bv_pass_1;
        std::vector<pasta::BitVector*> bv_pass_2;
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
                        // If we dont mark the blocks and instead build the whole tree we can skip the search
                        // div should be cheap, because we need to % anyway
                        size_type b = (mark) ? this->find_next_smallest_index_binary_search(ptr, block_text_inx) : ptr/block_size;
                        size_type current_offset = ptr % block_size;
                        counter_lvl[b]++;
                        if (current_offset != 0) {
                            counter_lvl[b + 1]++;
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
            bv_pass_1.push_back(&bv);

            for (auto c: counter_lvl) {
                std::cout << c << ",";
            }
            std::cout << std::endl;
        }
        this->leaf_size = block_size;
        block_size *= this->tau_;
        std::vector<size_type> pos_lvl(bv_pass_1.size(),0);
        for (size_type i = 0; i < bv_pass_1.size(); i++) {

            pos_lvl[i] = bv_pass_1[i]->size() - 1;
        }
        auto& first_lvl_bv = (*bv_pass_1[0]);
        for (size_type i = first_lvl_bv.size() - 1; i >= 0; i++) {
            if (first_lvl_bv[i]) {
                for (size_type j = bv_pass_1.size() - 2; j > 0; j++) {
                    auto& lvl_bv = (*bv_pass_1[j]);
                    if (lvl_bv[i]) {
                        for (size_type k = this->tau_ - 1; k >= 0; k++) {
                            std::cout << j << " " << pos_lvl[j + 1] - k;
                        }
                        pos_lvl[j] -= this->tau_;
                    }
                }
            }
        }
        for (size_type i = bv_pass_1.size() - 2; i >= 0; i--) {
            auto bv = pasta::BitVector(blk_lvl[i].size(),0);
            // mark nodes that have children
            size_type marked_counter = 0;
            for (size_type j = 0; j < bv.size(); j++) {
                auto& lvl_bv = (*bv_pass_1[i]);
                if (lvl_bv[j] == 1) {
                    if (counter[i][j] == 0) {

                    }
                    for (size_type k = 0; k < this->tau_; k++) {
                        auto& lvl_below = (*bv_pass_2[bv_pass_2.size() - 1]);
                        if (lvl_below[marked_counter * this->tau_ + k] == 1) {
                            bv[j] = 1;
                        }
                    }
                    marked_counter++;
                }
            }
            size_type max_pointer = 0;
            size_type max_offset = 0;
            auto pointers = std::vector<size_type>();
            auto offsets = std::vector<size_type>();
            for (size_type j = blk_lvl[i].size() - 1; j >= 0; j--) {
                if (bv[j] == 1) {
                    continue;
                }
                size_type first_ind = blk_lvl[i][j];
                size_type ind = first_ind;
                bool has_ptr = false;
                while (lpf[ind] >= block_size) {
                    size_type ptr = prevOcc[ind];
                    if (block_size + prevOcc[ind] - 1 >= first_ind) {
                        ind = prevOcc[ind];
                    } else if (lpf[prevOcc[ind]] >= block_size) {
                        ind = prevOcc[ind];
                    } else {
                        size_type b = this->find_next_smallest_index_binary_search(ptr, blk_lvl[i]);

                        size_type current_offset = ptr % block_size;
                        bv[b] = 1;
                        if (current_offset != 0) {
                            bv[b + 1] = 1;
                        }
                        pointers.push_back(b);
                        offsets.push_back(current_offset);
                        if (b > max_pointer) {
                            max_pointer = b;
                        }
                        if (current_offset > max_offset) {
                            max_offset = current_offset;
                        }
                        has_ptr = true;
                        break;
                    }
                }
                if (!has_ptr) {
                    bv[j] = 1;
                }
//                if (j == 340) {
//                    std::cout<< i << " " << (*bv)[j] << std::endl;
//                }
            }

            block_size *= this->tau_;
            bv_pass_2.push_back(&bv);
//            if (i == bv_pass_1.size() -1) {
//                std::cout << *bv << std::endl;
//            }
            size_type ones_per_pass2 = 0;
            for (size_type j = 0; j < bv.size(); j++) {
                ones_per_pass2 += bv[j];
            }
            pass2_ones.push_back(ones_per_pass2);
            pass2_pointer.push_back(pointers);
            pass2_offset.push_back(offsets);
            pass2_max_offset.push_back(max_offset);
            pass2_max_pointer.push_back(max_pointer);

        }
        //determine lvl final size
        // prepare first level;
        auto& lvl = (*bv_pass_2[bv_pass_2.size() - 1]);
        this->block_tree_types_.push_back(&lvl);
        this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(lvl));
        auto size = pass2_pointer[pass2_pointer.size() -1].size();
        auto p1 = new sdsl::int_vector<>(size,0,(8 * sizeof(size_type)) -  this->leading_zeros(pass2_max_pointer[pass2_max_pointer.size() -1]));
        auto o1 = new sdsl::int_vector<>(size,0,(8 * sizeof(size_type)) -  this->leading_zeros(pass2_max_offset[pass2_max_offset.size() -1]));
        auto& ptr1 = *p1;
        auto& off1 = *o1;
        for(size_type i = 0; i < pass2_pointer[pass2_pointer.size() - 1].size(); i++) {
            ptr1[i] = pass2_pointer[pass2_pointer.size() - 1][size- 1 -i];
            off1[i] = pass2_offset[pass2_offset.size() - 1][size - 1 -i];
        }
        this->block_tree_pointers_.push_back(p1);
        this->block_tree_offsets_.push_back(o1);
        for (size_type i = bv_pass_2.size() - 2; i >= 0; i--) {
            auto pass1_i = bv_pass_2.size() - 1 - i;
            auto pass1_parent = pass1_i - 1;
            auto new_size = this->tau_ * (pass2_ones[i + 1] - is_padded);
            // the last block only spawns block that contain text
            auto last_block_parent = blk_lvl[pass1_parent][blk_lvl[pass1_parent].size() - 1];
            auto lvl_block_size = this->block_size_lvl_[pass1_i];
            if (is_padded) {
                for (size_type j = 0; j < this->tau_; j++) {
                    if (last_block_parent + j * lvl_block_size < text.size()) {
                        new_size++;
                    }
                }
            }

            auto bit_vector = new pasta::BitVector(new_size,0);
            auto& bv_ref = *bit_vector;
            auto pointer = std::vector<size_type>();
            auto offset = std::vector<size_type>();
            size_type pointer_saved = 0;
            size_type pointer_skipped = 0;
            std::unordered_map<size_type, size_type> blocks_skipped;
            size_type skip = 0;
            size_type replace = 0;
            auto& lvl_below_pass1 = *bv_pass_1[pass1_i - 1];
            auto& lvl_below_pass2 = *bv_pass_2[i + 1];
            auto& lvl_pass2 = *bv_pass_2[i];
            for (size_type j = 0; j < lvl_below_pass1.size(); j++) {
                if (lvl_below_pass1[j] == 1) {
                    if (lvl_below_pass2[j] == 1) {
                        for (size_type k = 0; k < this->tau_ && replace * this->tau_ + k < new_size; k++) {
                            bool x = lvl_pass2[(j - skip) * this->tau_ + k];
                            auto skipper = pointer_skipped + pointer_saved;
                            bv_ref[replace * this->tau_ + k] = x;
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
                        pointer_skipped+= this->tau_;
                    }
                } else {
                    skip++;
                }
            }
            auto p = new sdsl::int_vector<>(pointer.size(),0,(8 * sizeof(size_type)) -  this->leading_zeros(pass2_max_pointer[i]));
            auto o = new sdsl::int_vector<>(pointer.size(),0,(8 * sizeof(size_type)) -  this->leading_zeros(pass2_max_offset[i]));
            auto& ptr = *p;
            auto& off = *o;
            for(size_type j = 0; j < pointer.size(); j++) {
                ptr[j] = pointer[j];
                off[j] = offset[j];
            }
            this->block_tree_pointers_.push_back(p);
            this->block_tree_offsets_.push_back(o);
            this->block_tree_types_.push_back(bit_vector);
//            std::cout << (*bv_pass_1[pass1_i - 1]) << std::endl;
//            std::cout << (*bv_pass_2[i + 1]) << std::endl;
//            std::cout << *bit_vector << std::endl;

            this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(bv_ref));
        }
        int64_t leaf_count = 0;
        auto& last_level = (*bv_pass_2[0]);
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
        // final pass
//        for (auto a: bv_pass_1) {
//            std::cout << *a << std::endl << std::endl;
//        }
//        std::cout << "alles kllar" << std::endl;
//        for (auto a: bv_pass_2) {
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
        for (size_type i = 0; i < lz.size(); i++) {
            while (block_text_inx[j] + block_size < lz[i]) {
                j++;
            }

            bv[j] = 1;
//            if (j > 0) {
//                std::cout << block_text_inx[j-1] << " " << block_text_inx[j] << " " << block_size << std::endl;
//            }
//            && block_text_inx[j-1] + block_size == block_text_inx[j])
            if (j > 0 && block_text_inx[j-1] + block_size == block_text_inx[j]){
                bv[j - 1] = 1;
            }
//            && block_text_inx[j] + block_size == block_text_inx[j + 1]
            if (j + 1 < bv.size() && block_text_inx[j] + block_size == block_text_inx[j + 1]) {
                bv[j + 1] = 1;
            }
        }
        return 0;
    }
};
#endif //BLOCK_TREE_BV_BLOCKTREE_LZ_HPP

