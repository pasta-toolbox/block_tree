//
// Created by daniel on 03.09.22.
//
#include "bv_blocktree.hpp"
#ifndef BLOCK_TREE_BV_BLOCKTREE_THEORY_H
#define BLOCK_TREE_BV_BLOCKTREE_THEORY_H

template<typename input_type, typename size_type>
class BV_BlockTree_lpf_theory: public BV_Block_Tree<input_type,size_type> {
public:
    int32_t init(std::vector<input_type>& text, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, std::vector<size_type>& lz) {
        int64_t added_padding = 0;
        int64_t tree_max_height = 0;
        int64_t max_blk_size = 0;
        this->calculate_padding(added_padding, text.size(), tree_max_height, max_blk_size);
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

        while (block_size > this->max_leaf_length_) {
            std::cout << block_size << std::endl;
            size_type max_pointer = 0;
            size_type max_offset = 0;
            this->block_size_lvl_.push_back(block_size);
            this->block_per_lvl_.push_back(block_text_inx.size());
            auto bit_vector = new pasta::BitVector(block_text_inx.size(),0);
            auto& bv = *bit_vector;
            auto pointers = std::vector<size_type>();
            auto offsets = std::vector<size_type>();
            mark_blocks(bv, lz, block_text_inx, block_size);
            for (size_type i = 0; i < block_text_inx.size(); i++) {
                if (bv[i] == 1) {
                    continue;
                }
                size_type first_ind = block_text_inx[i];
                size_type ind = first_ind;
                while (lpf[ind] >= block_size) {
                    size_type ptr = lpf_ptr[ind];
                    if (block_size + lpf_ptr[ind] - 1 >= first_ind) {
                        ind = lpf_ptr[ind];
                    } else if (lpf[lpf_ptr[ind]] >= block_size) {
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
                        break;
                    }
                }
            }
            block_size /= this->tau_;
            std::vector<int64_t> block_text_inx_new;
            for (size_type i = 0; i < bv.size(); i++) {
                if (bv[i] == 1) {
                    for (size_type j = 0; j < this->tau_ ; j++) {
                        if (block_text_inx[i] + (j * block_size) < text.size()) {
                            block_text_inx_new.push_back(block_text_inx[i] + (j * block_size));
                        }
                    }
                }
            }
            this->block_tree_types_.push_back(bit_vector);
            this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(bv));
            size_type ptr_size = (8 * sizeof(max_pointer)) - this->leading_zeros(max_pointer);
            size_type off_size = (8 * sizeof(max_offset)) - this->leading_zeros(max_offset);
            auto p = new sdsl::int_vector(pointers.size(), 0,ptr_size);
            auto o = new sdsl::int_vector(offsets.size(),0,off_size);
            auto& ptr = *p;
            auto& off = *o;
            for(int j = 0; j < pointers.size(); j++) {
                auto pointer = pointers[j];
                auto offset = offsets[j];
                ptr[j] = pointer;
                off[j] = offset;
            }
            this->block_tree_pointers_.push_back(p);
            this->block_tree_offsets_.push_back(o);
            block_text_inx = block_text_inx_new;

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
    BV_BlockTree_lpf_theory(std::vector<input_type>& text, size_type tau, size_type max_leaf_length) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        // first we create lpf and lpf_ptr arrays;
        std::vector<size_type> lpf(text.size());
        std::vector<size_type> lpf_ptr(text.size());
        std::vector<size_type> lz;
        lpf_array(text, lpf, lpf_ptr);
        calculate_lz_factor(this->s_,lpf, lz);
        init(text, lpf, lpf_ptr, lz);
    };
    BV_BlockTree_lpf_theory(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, std::vector<size_type>& lz) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        // first we create lpf and lpf_ptr arrays;
        this->s_ = lz.size();
        init(text, lpf, lpf_ptr, lz);
    };
    BV_BlockTree_lpf_theory(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, size_type s, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, std::vector<size_type>& lz) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        // first we create lpf and lpf_ptr arrays;
        this->s_ = s;
        init(text, lpf, lpf_ptr, lz);
    };
    ~BV_BlockTree_lpf_theory() {

    };
private:
    size_type generate_next_level(std::vector<size_type> &old_level, std::vector<size_type> &new_level, pasta::BitVector* bv, size_type N, size_type block_size) {
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

    size_type mark_blocks(pasta::BitVector& bv, std::vector<size_type>& lz, std::vector<int64_t>& block_text_inx, int64_t block_size) {
        int blocks_marked = 0;
        size_type j = 0;
        for (size_type i = 0; i < lz.size(); i++) {

            size_type f = lz[i];

            while (j < block_text_inx.size() -1 && block_text_inx[j + 1] <= f) {
                j++;
            }
            bv[j] = 1;
//            if (j > 0) {
//                std::cout << block_text_inx[j-1] << " " << block_text_inx[j] << " " << block_size << std::endl;
//            }
            if (j > 0 && block_text_inx[j-1] + block_size == block_text_inx[j]) {
                bv[j - 1] = 1;
            }

            if (j + 1< bv.size() && block_text_inx[j] + block_size == block_text_inx[j + 1]) {
                bv[j + 1] = 1;
            }
        }
        return 0;
    }
};
#endif //BLOCK_TREE_BV_BLOCKTREE_THEORY_H
