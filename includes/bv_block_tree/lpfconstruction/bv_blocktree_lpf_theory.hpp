//
// Created by daniel on 03.09.22.
//
#include "bv_blocktree.hpp"
#ifndef BLOCK_TREE_BV_BLOCKTREE_THEORY_H
#define BLOCK_TREE_BV_BLOCKTREE_THEORY_H

template<typename input_type, typename size_type>
class BV_BlockTree_lpf_theory: public BV_Block_Tree<input_type,size_type> {
public:
    int32_t init_dp(std::vector<input_type>& text, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr) {
        int64_t added_padding = 0;
        int64_t tree_max_height = 0;
        int64_t max_blk_size = 0;
        std::vector<int64_t> block_size_lvl_temp;
        this->calculate_padding(added_padding, text.size(), tree_max_height, max_blk_size);
        int64_t block_size = max_blk_size;
        std::vector<int64_t> block_text_inx;
        for (int64_t i = 0; i < text.size(); i+= block_size) {
            block_text_inx.push_back(i);
        }
        for (size_type i = 0; i < lpf_ptr.size(); i++) {
            size_type p = lpf_ptr[i];
            if (lpf[i] >= block_size && (lpf[p] >= block_size || lpf[i] <= lpf[p])) {
                lpf_ptr[i] = lpf_ptr[p];
            }
        }
        bool found_back_block = false;
        size_type level = 0;
        while (block_size > this->max_leaf_length_) {

            size_type max_pointer = 0;
            size_type max_offset = 0;
            block_size_lvl_temp.push_back(block_size);
            auto bit_vector = new pasta::BitVector(block_text_inx.size(),false);
            auto& bv = *bit_vector;
            auto pointers = std::vector<size_type>();
            auto offsets = std::vector<size_type>();
            mark_blocks_percise(bv, lpf, block_text_inx, block_size);

            size_type counter = 0;
            for (int i = 0; i < bv.size(); i++) {
                if (bv[i] == 0) counter++;
            }

            for (size_type i = 0; i < block_text_inx.size(); i++) {
                if (bv[i] == 0) {
                    size_type first_ind = block_text_inx[i];
                    size_type ind = first_ind;
                    if (lpf[ind] >= block_size) {
                        size_type ptr = lpf_ptr[ind];
                        if (block_size + ptr - 1 < first_ind && lpf[ptr] < block_size) {
                            size_type b = this->find_next_smallest_index_binary_search(ptr, block_text_inx);
                            size_type current_offset = ptr % block_size;
                            pointers.push_back(b);
                            offsets.push_back(current_offset);
                        }
                    }

                }
            }
            block_size /= this->tau_;
            for (size_type b = 0; b < block_text_inx.size(); b++) {
                for (size_type j = 0; j < block_size * this->tau_ && block_text_inx[b] + j < text.size(); j++) {
                    size_type i = block_text_inx[b] + j;
                    size_type p = lpf_ptr[i];
                    if ((lpf[i] >= block_size && lpf[p] >= block_size) || ((p != -1) && lpf[i] <= lpf[p])) {
                        lpf_ptr[i] = lpf_ptr[p];
                    }
                }
            }
            std::vector<int64_t> block_text_inx_new;
            size_type ones_in_lvl = 0;
            for (size_type i = 0; i < bv.size(); i++) {
                if (bv[i] == 1) {
                    ones_in_lvl++;
                    for (size_type j = 0; j < this->tau_ ; j++) {
                        if (block_text_inx[i] + (j * block_size) < text.size()) {
                            block_text_inx_new.push_back(block_text_inx[i] + (j * block_size));
                        }
                    }
                }
            }
            block_text_inx = block_text_inx_new;
            found_back_block |= ones_in_lvl != bv.size();
            if (found_back_block || !this->CUT_FIRST_LEVELS) {
                this->block_tree_types_.push_back(bit_vector);
                this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(bv));
                auto p = new sdsl::int_vector<>(pointers.size(), 0);
                auto o = new sdsl::int_vector<>(offsets.size(),0);
                auto& ptr = *p;
                auto& off = *o;
                for(int j = 0; j < pointers.size(); j++) {
                    auto pointer = pointers[j];
                    auto offset = offsets[j];
                    ptr[j] = pointer;
                    off[j] = offset;
                }
                sdsl::util::bit_compress(ptr);
                sdsl::util::bit_compress(off);
                this->block_tree_pointers_.push_back(p);
                this->block_tree_offsets_.push_back(o);
                this->block_size_lvl_.push_back(block_size * this->tau_);
            }
            level++;

        }
        this->leaf_size = block_size;
        this->amount_of_leaves = block_text_inx.size();
        for (size_type ptr: block_text_inx) {
            for (int i = 0; i < block_size; i++) {
                if (ptr + i < text.size()) {
                    this->leaves_.push_back(text[ptr + i]);
                }
            }
        }
        return 0;
    };
    int32_t init(std::vector<input_type>& text, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr) {
        int64_t added_padding = 0;
        int64_t tree_max_height = 0;
        int64_t max_blk_size = 0;
        std::vector<int64_t> block_size_lvl_temp;
        this->calculate_padding(added_padding, text.size(), tree_max_height, max_blk_size);
        int64_t block_size = max_blk_size;
        std::vector<int64_t> block_text_inx;
        for (int64_t i = 0; i < text.size(); i+= block_size) {
            block_text_inx.push_back(i);
        }
        for (size_type i = 0; i < lpf_ptr.size(); i++) {
            size_type p = lpf_ptr[i];
            if ((lpf[i] >= block_size && lpf[p] >= block_size) || ((p != -1) && lpf[i] <= lpf[p])) {
                lpf_ptr[i] = lpf_ptr[p];
            }
        }
        bool found_back_block = false;
        size_type level = 0;
        while (block_size > this->max_leaf_length_) {

            size_type max_pointer = 0;
            size_type max_offset = 0;
            block_size_lvl_temp.push_back(block_size);
            auto bit_vector = new pasta::BitVector(block_text_inx.size(),false);
            auto& bv = *bit_vector;
            auto pointers = std::vector<size_type>();
            auto offsets = std::vector<size_type>();
            mark_blocks_percise(bv, lpf, block_text_inx, block_size);

            size_type counter = 0;
            for (int i = 0; i < bv.size(); i++) {
                if (bv[i] == 0) counter++;
            }

            for (size_type i = 0; i < block_text_inx.size(); i++) {
                if (bv[i] == 0) {
                    bool set = false;
                    size_type first_ind = block_text_inx[i];
                    size_type ind = first_ind;

                    while (lpf[ind] >= block_size) {

                        size_type ptr = lpf_ptr[ind];

                        if (block_size + ptr - 1 >= first_ind) {
                            ind = lpf_ptr[ind];
                        } else if (lpf[ptr] >= block_size) {
                            ind = lpf_ptr[ind];
                        } else {

                            size_type b = this->find_next_smallest_index_binary_search(ptr, block_text_inx);
                            size_type current_offset = ptr % block_size;
                            pointers.push_back(b);
                            offsets.push_back(current_offset);
                            if (b > max_pointer) {
                                max_pointer = b;
                            }
                            if (current_offset > max_offset) {
                                max_offset = current_offset;
                            }
                            set = true;
                            break;
                        }
                    }

                }
            }
            block_size /= this->tau_;
//            for (size_type b = 0; b < block_text_inx.size(); b++) {
//                for (size_type j = 0; j < block_size * this->tau_; j++) {
//                    size_type i = block_text_inx[b] + j;
//                    size_type p = lpf_ptr[i];
//                    if ((lpf[i] >= block_size && lpf[p] >= block_size) || ((p != -1) && lpf[i] <= lpf[p])) {
//                        lpf_ptr[i] = lpf_ptr[p];
//                    }
//                }
//            }
            std::vector<int64_t> block_text_inx_new;
            size_type ones_in_lvl = 0;
            for (size_type i = 0; i < bv.size(); i++) {
                if (bv[i] == 1) {
                    ones_in_lvl++;
                    for (size_type j = 0; j < this->tau_ ; j++) {
                        if (block_text_inx[i] + (j * block_size) < text.size()) {
                            block_text_inx_new.push_back(block_text_inx[i] + (j * block_size));
                        }
                    }
                }
            }
            block_text_inx = block_text_inx_new;
            found_back_block |= ones_in_lvl != bv.size();
            if (found_back_block || !this->CUT_FIRST_LEVELS) {
                this->block_tree_types_.push_back(bit_vector);
                this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(bv));
                auto p = new sdsl::int_vector<>(pointers.size(), 0);
                auto o = new sdsl::int_vector<>(offsets.size(),0);
                auto& ptr = *p;
                auto& off = *o;
                for(int j = 0; j < pointers.size(); j++) {
                    auto pointer = pointers[j];
                    auto offset = offsets[j];
                    ptr[j] = pointer;
                    off[j] = offset;
                }
                sdsl::util::bit_compress(ptr);
                sdsl::util::bit_compress(off);
                this->block_tree_pointers_.push_back(p);
                this->block_tree_offsets_.push_back(o);
                this->block_size_lvl_.push_back(block_size * this->tau_);
            }
            level++;

        }
        this->leaf_size = block_size;
        this->amount_of_leaves = block_text_inx.size();
        for (size_type ptr: block_text_inx) {
            for (int i = 0; i < block_size; i++) {
                if (ptr + i < text.size()) {
                    this->leaves_.push_back(text[ptr + i]);
                }
            }
        }
        return 0;
    };
    BV_BlockTree_lpf_theory(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, bool dp) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        // first we create lpf and lpf_ptr arrays;
        std::vector<size_type> lpf(text.size());
        std::vector<size_type> lpf_ptr(text.size());
        std::vector<size_type> lz;
        lpf_array_stack(text, lpf, lpf_ptr);
        calculate_lz_factor(this->s_,lpf, lz);
        if (dp) init_dp(text, lpf, lpf_ptr);
        else init(text, lpf, lpf_ptr);
    };
    BV_BlockTree_lpf_theory(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, size_type s, bool dp) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        this->s_ = s;
        // first we create lpf and lpf_ptr arrays;
        std::vector<size_type> lpf(text.size());
        std::vector<size_type> lpf_ptr(text.size());
        lpf_array_stack(text, lpf, lpf_ptr);
        if (dp) init_dp(text, lpf, lpf_ptr);
        else init(text, lpf, lpf_ptr);
    };
    BV_BlockTree_lpf_theory(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, std::vector<size_type>& lz) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        // first we create lpf and lpf_ptr arrays;
        this->s_ = lz.size();
        init_dp(text, lpf, lpf_ptr);
    };
    BV_BlockTree_lpf_theory(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, size_type s, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, std::vector<size_type>& lz) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        // first we create lpf and lpf_ptr arrays;
        this->s_ = s;
        init_dp(text, lpf, lpf_ptr);
    };
    ~BV_BlockTree_lpf_theory() {
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
    size_type generate_next_level(std::vector<size_type> &old_level, std::vector<size_type> &new_level, pasta::BitVector* bv, size_type N, size_type block_size) {

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


        size_type mark_blocks(pasta::BitVector& bv, std::vector<size_type>& lz, std::vector<int64_t>& block_text_inx, int64_t block_size) {
            size_type j = 0;
            for (size_type i = 0; i < lz.size() - 1; i++) {
                size_type f = lz[i];
                while (j < block_text_inx.size() - 1 && block_text_inx[j + 1] <= f) {
                    j++;
                }
                bv[j] = 1;

                if (j > 0 && block_text_inx[j - 1] + block_size == block_text_inx[j] && f != block_text_inx[j] + block_size) {
                    bv[j - 1] = 1;
                }

                if (j + 1 < bv.size() && block_text_inx[j] + block_size == block_text_inx[j + 1] && f != block_text_inx[j]) {
                    bv[j + 1] = 1;
                }
            }
            return 0;

    }
    size_type mark_blocks_percise(pasta::BitVector& bv, std::vector<size_type>& lpf, std::vector<int64_t>& block_text_inx, int64_t block_size) {
        bv[0] = 1;
        for (size_type i = 1; i < block_text_inx.size() - 1; i++) {
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
    }
};
#endif //BLOCK_TREE_BV_BLOCKTREE_THEORY_H
