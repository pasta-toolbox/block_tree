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
    int32_t init(std::vector<input_type>& text, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, std::vector<size_type>& lz) {
        std::vector<pasta::BitVector*> bv_pass_1;
        std::vector<pasta::BitVector*> bv_pass_2;
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
        for(size_type i = 0; i < lpf_ptr.size(); i++) {
            if (lpf[i] <= lpf[lpf_ptr[i]] && lpf_ptr[i] == i-1) {
                lpf_ptr[i] = lpf_ptr[lpf_ptr[i]];
            }
        }
        do {
            auto* bv = new pasta::BitVector(block_text_inx.size(),0);
            mark_blocks(bv, lz, block_text_inx, block_size);
            this->block_size_lvl_.push_back(block_size);
            blk_lvl.push_back(block_text_inx);
            block_size = block_size/ this->tau_;
            std::vector<int64_t> block_text_inx_new;
            generate_next_level(block_text_inx,block_text_inx_new, bv, text.size(), block_size);
            block_text_inx = block_text_inx_new;
            bv_pass_1.push_back(bv);
//            if (block_size <= this->max_leaf_length_ ) {
//                std::cout << *bv << std::endl << std::endl;
//            }
        } while (block_size > this->max_leaf_length_);
        this->leaf_size = block_size;
        block_size *= this->tau_;
        for (size_type i = bv_pass_1.size() - 1; i >= 0; i--) {
            auto* bv = new pasta::BitVector(blk_lvl[i].size(),0);
            size_type marked_counter = 0;
            if (i != bv_pass_1.size() - 1) {
                for (size_type j = 0; j < bv->size(); j++) {
                    if ((*bv_pass_1[i])[j] == 1) {
                        for (size_type k = 0; k < this->tau_; k++) {
                            if ((*bv_pass_2[bv_pass_2.size() - 1])[marked_counter * this->tau_ + k] == 1) {
                                (*bv)[j] = 1;
                            }
                        }
                        marked_counter++;
                    }
                }
            }
            size_type max_pointer = 0;
            size_type max_offset = 0;
            auto pointers = std::vector<size_type>();
            auto offsets = std::vector<size_type>();
            for (size_type j = blk_lvl[i].size() - 1; j >= 0; j--) {
                if ((*bv)[j] == 1) {
                    continue;
                }
                size_type first_ind = blk_lvl[i][j];
                size_type ind = first_ind;
                bool has_ptr = false;
                while (lpf[ind] >= block_size) {
                    size_type ptr = lpf_ptr[ind];
                    if (block_size + lpf_ptr[ind] - 1 >= first_ind) {
                        ind = lpf_ptr[ind];
                    } else if (lpf[lpf_ptr[ind]] >= block_size) {
                        ind = lpf_ptr[ind];
                    } else {
                        size_type b = this->find_next_smallest_index_binary_search(ptr, blk_lvl[i]);
//                        if (b == 340) {
//                            std::cout<< i << " marks 340: " << j << std::endl;
//                        }
                        size_type current_offset = ptr % block_size;
                        (*bv)[b] = 1;
                        if (current_offset != 0) {
                            (*bv)[b + 1] = 1;
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
                    (*bv)[j] =1;
                }
//                if (j == 340) {
//                    std::cout<< i << " " << (*bv)[j] << std::endl;
//                }
            }

            block_size *= this->tau_;
            bv_pass_2.push_back(bv);
//            if (i == bv_pass_1.size() -1) {
//                std::cout << *bv << std::endl;
//            }
            size_type ones_per_pass2 = 0;
            for (size_type j = 0; j < bv->size(); j++) {
                ones_per_pass2 += (*bv)[j];
            }
            pass2_ones.push_back(ones_per_pass2);
            pass2_pointer.push_back(pointers);
            pass2_offset.push_back(offsets);
            pass2_max_offset.push_back(max_offset);
            pass2_max_pointer.push_back(max_pointer);
        }
        //determine lvl final size
        // prepare first level;
        this->block_tree_types_.push_back(bv_pass_2[bv_pass_2.size() - 1]);
        this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv_pass_2[bv_pass_2.size() - 1]));
        auto size = pass2_pointer[pass2_pointer.size() -1].size();
        auto p1 = new sdsl::int_vector<>(size,0,(8 * sizeof(size_type)) -  this->leading_zeros(pass2_max_pointer[pass2_max_pointer.size() -1]));
        auto o1 = new sdsl::int_vector<>(size,0,(8 * sizeof(size_type)) -  this->leading_zeros(pass2_max_offset[pass2_max_offset.size() -1]));
        for(size_type i = 0; i < pass2_pointer[pass2_pointer.size() - 1].size(); i++) {
            (*p1)[i] = pass2_pointer[pass2_pointer.size() - 1][size- 1 -i];
            (*o1)[i] = pass2_offset[pass2_offset.size() - 1][size - 1 -i];
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
//            std::cout << "block size revi " << this->block_size_lvl_[pass1_i] << " " << this->block_size_lvl_[pass1_i + 1] << std::endl;
//            std::cout << "new size " << new_size << " " << is_padded << " " << last_block_parent << std::endl;
            if (is_padded) {
                for (size_type j = 0; j < this->tau_; j++) {
//                    std::cout << last_block_parent + j * lvl_block_size << " " << text.size() << std::endl;
                    if (last_block_parent + j * lvl_block_size < text.size()) {
                        new_size++;
                    }
                }
            }
//            std::cout << "new size " << new_size << " " << is_padded << " " << last_block_parent << std::endl;
            auto* bit_vector = new pasta::BitVector(new_size,0);
            auto pointer = std::vector<size_type>();
            auto offset = std::vector<size_type>();
            size_type pointer_saved = 0;
            size_type pointer_skipped = 0;
            std::unordered_map<size_type, size_type> blocks_skipped;
            size_type skip = 0;
            size_type replace = 0;
            for (size_type j = 0; j < bv_pass_1[pass1_i - 1]->size(); j++) {
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
                        pointer_skipped+= this->tau_;
                    }
                } else {
                    skip++;
                }
            }
//            std::cout << i << " skipped " << skip << " pointer skipped " << pointer_skipped << " pointer saved " << pointer_saved << " rep " << replace << std::endl;
//            for (size_type j = 0; j < bv_pass_2[i + 1]->size(); j++) {
//                if ((*bv_pass_2[i + 1])[j] == 1) {
//                    for (size_type k = 0; k < this->tau_; k++) {
//                        bool x = (*bv_pass_2[i])[higher_lvl_ones * this->tau_ + k];
////                        std::cout << j * this->tau_ + k << " "<< parent * this->tau_ + k << " " << (*bv_pass_2[i]).size() << std::endl;
//                        (*bit_vector)[higher_lvl_ones * this->tau_ + k] = x;
//                    }
//                    higher_lvl_ones++;
//                }
//            }
//            std::cout << *bit_vector << std::endl;
            auto p = new sdsl::int_vector<>(pointer.size(),0,(8 * sizeof(size_type)) -  this->leading_zeros(pass2_max_pointer[i]));
            auto o = new sdsl::int_vector<>(pointer.size(),0,(8 * sizeof(size_type)) -  this->leading_zeros(pass2_max_offset[i]));
            for(size_type j = 0; j < pointer.size(); j++) {
                (*p)[j] = pointer[j];
                (*o)[j] = offset[j];
            }

            this->block_tree_pointers_.push_back(p);
            this->block_tree_offsets_.push_back(o);
            this->block_tree_types_.push_back(bit_vector);
            this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bit_vector));
        }
        int64_t leaf_count = 0;
        for (size_type i = 0; i < bv_pass_2[0]->size(); i++) {
            if ((*bv_pass_2[0])[i] == 1) {
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
        for (auto bv : bv_pass_1) {
            delete bv;
        }
        for (size_type i = 0; i < bv_pass_2.size() - 1; i++) {
            delete bv_pass_2[i];
        }
        return 0;
    };
    BV_BlockTree_lpf_pruned(std::vector<input_type>& text, size_type tau, size_type max_leaf_length) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        std::vector<size_type> lpf(text.size());
        std::vector<size_type> lpf_ptr(text.size());
        std::vector<size_type> lz;
        lpf_array(text, lpf, lpf_ptr);
        calculate_lz_factor(this->s_,lpf, lz);
        init(text, lpf, lpf_ptr, lz);
    };
    BV_BlockTree_lpf_pruned(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, std::vector<size_type>& lz) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        this->s_ = lz.size();
        init(text, lpf, lpf_ptr, lz);
    };
    BV_BlockTree_lpf_pruned(std::vector<input_type>& text, size_type tau, size_type max_leaf_length,size_type s, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, std::vector<size_type>& lz) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        this->s_ = s;
        init(text, lpf, lpf_ptr, lz);
    };
    ~BV_BlockTree_lpf_pruned() {
        for (auto a : this->block_tree_types_) {
            delete a;
        }
        for (auto a: this->block_tree_types_rs_) {
            delete a;
        }
        for (auto a: this->block_tree_pointers_) {
            delete a;
        }
        for (auto a: this->block_tree_offsets_) {
            delete a;
        }

    };
private:
    size_type generate_next_level(std::vector<int64_t> &old_level, std::vector<int64_t> &new_level, pasta::BitVector* bv, int64_t N, int64_t block_size) {
//        std::cout << block_size << " this size " << " " << old_level.size() << " " << new_level.size() << std::endl;
        for (size_type i = 0; i < bv->size(); i++) {
            if ((*bv)[i] == 1) {
                for (size_type j = 0; j < this->tau_ ; j++) {
                    if (old_level[i] + (j * block_size) < N) {
                        new_level.push_back(old_level[i] + (j * block_size));
                    }
                }
            }
        }
        return 0;
    }

    size_type mark_blocks(pasta::BitVector* bv, std::vector<size_type>& lz, std::vector<int64_t>& block_text_inx, int64_t block_size) {
        int blocks_marked = 0;
        size_type j = 0;
        for (size_type i = 0; i < lz.size(); i++) {

            size_type f = lz[i];

            while (j < block_text_inx.size() -1 && block_text_inx[j + 1] <= f) {
                j++;
            }
            (*bv)[j] = 1;
//            if (j > 0) {
//                std::cout << block_text_inx[j-1] << " " << block_text_inx[j] << " " << block_size << std::endl;
//            }
            if (j > 0 && block_text_inx[j-1] + block_size == block_text_inx[j]) {
                (*bv)[j - 1] = 1;
            }

            if (j + 1< bv->size() && block_text_inx[j] + block_size == block_text_inx[j + 1]) {
                (*bv)[j + 1] = 1;
            }
        }
        return 0;
    }
};
#endif //BLOCK_TREE_BV_BLOCKTREE_LZ_HPP

