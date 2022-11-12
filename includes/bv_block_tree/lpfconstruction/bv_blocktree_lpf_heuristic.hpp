#include "bit_vector.hpp"
#include "support/find_l2_flat_with.hpp"
#include "support/flat_rank.hpp"
#include "support/flat_rank_select.hpp"
#include "support/optimized_for.hpp"
#include "support/rank.hpp"
#include "support/rank_select.hpp"
#include "support/wide_rank.hpp"
#include "support/wide_rank_select.hpp"
#include <chrono>
#include "int_vector.hpp"
#include "bv_blocktree.hpp"

#ifndef BLOCK_TREE_BV_BLOCKTREE_LPF_H
#define BLOCK_TREE_BV_BLOCKTREE_LPF_H
template<typename input_type, typename size_type>
class BV_BlockTree_lpf_heuristic: public BV_Block_Tree<input_type,size_type> {
public:

    BV_BlockTree_lpf_heuristic(std::vector<input_type>& text, size_type tau, size_type max_leaf_length) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        // first we create lpf and lpfptr arrays;
        std::vector<size_type> lpf(text.size());
        std::vector<size_type> lpf_ptr(text.size());
        std::vector<size_type> lz;
        lpf_array(text, lpf, lpf_ptr);
        calculate_lz_factor(this->s_,lpf, lz);
        init(text, lpf, lpf_ptr, lz);
    };
    BV_BlockTree_lpf_heuristic(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, std::vector<size_type>& lz) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        // first we create lpf and lpfptr arrays;
        this->s_ = lz.size();
        init(text, lpf, lpf_ptr, lz);
    };
    BV_BlockTree_lpf_heuristic(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, size_type s, std::vector<size_type>& lpf, std::vector<size_type>& lpf_ptr, std::vector<size_type>& lz) {
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        // first we create lpf and lpfptr arrays;;
        this->s_ = s;
        init(text, lpf, lpf_ptr, lz);
    };
    ~BV_BlockTree_lpf_heuristic() {
    };

private:
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
//        manipulate lpf_ptr
        for(size_type i = 0; i < lpf_ptr.size(); i++) {
            if (lpf[i] <= lpf[lpf_ptr[i]] && lpf_ptr[i] == i-1) {
                lpf_ptr[i] = lpf_ptr[lpf_ptr[i]];
            }
        }
        bool found_back_block = false;
        size_type level = 0;
        do {

//            std::cout << "Blocksize: " << block_size << " Blocks in Level: " << block_text_inx.size();

            this->block_per_lvl_.push_back(block_text_inx.size());
            auto bit_vector = new pasta::BitVector(block_text_inx.size(),false);
            auto& bv = *bit_vector;
            std::vector<size_type> pointers = std::vector<size_type>();
            std::vector<size_type> offsets = std::vector<size_type>();
            for(int64_t block_in_lvl =(int64_t)block_text_inx.size() - 1; block_in_lvl >= 0; block_in_lvl--) {
                if (bv[block_in_lvl] == 1) {
                    continue;
                }
                int64_t initial_index = block_text_inx[block_in_lvl];
                int64_t ind = initial_index;
                int64_t current_pointer;
                int64_t current_offset;
                bool has_ptr = false;
                while (lpf[ind] >= block_size) {
                    int64_t ptr = lpf_ptr[ind];
                    if (block_size + lpf_ptr[ind] - 1 >= initial_index) {
                        ind = lpf_ptr[ind];
                    } else if (lpf[lpf_ptr[ind]] >= block_size ) {
                        ind = lpf_ptr[ind];

                    } else {
                        int64_t b = this->find_next_smallest_index_binary_search(ptr, block_text_inx);
                        if (ptr <= block_text_inx[b] + block_size - 1 && ((block_text_inx[b] + block_size ==  block_text_inx[b+1]) || ptr % block_size == 0)) {
                            current_pointer = b;
                            current_offset = ptr % block_size;
                            has_ptr = true;
                        }
                        if (has_ptr) {
                            bv[current_pointer] = true;
                            if (ptr != block_text_inx[b]) {
                                bv[current_pointer + 1] = true;
                            }
                            pointers.push_back(current_pointer);
                            offsets.push_back(current_offset);
                        }
                        break;
                    }
                }
                if (!has_ptr) {
                    bv[block_in_lvl] = true;
                }
            }
            std::vector<int64_t> block_text_inx_new(0,0);
            block_size = block_size / this->tau_;
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
            found_back_block |= ones_in_lvl != bit_vector->size();
            if (found_back_block) {
                this->block_tree_types_.push_back(bit_vector);
                this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(bv));
                auto p = new sdsl::int_vector(pointers.size(), 0, 64);
                auto o = new sdsl::int_vector(offsets.size(), 0, 64);
                auto &ptr = *p;
                auto &off = *o;
                for (int64_t j = 0; j < pointers.size(); j++) {
                    auto pointer = pointers[pointers.size() - 1 - j];
                    auto offset = offsets[pointers.size() - 1 - j];
                    ptr[j] = pointer;
                    off[j] = offset;
                }
                this->block_size_lvl_.push_back(block_size * this->tau_);
                this->block_tree_pointers_.push_back(p);
                this->block_tree_offsets_.push_back(o);
            }
            block_text_inx = block_text_inx_new;
        } while (block_size > this->max_leaf_length_);
        this->leaf_size = block_size;
        this->amount_of_leaves = block_text_inx.size();
        for (size_type ptr: block_text_inx) {
            for (int64_t i = 0; i < block_size; i++) {
                if (ptr + i < text.size()) {
                    this->leaves_.push_back(text[ptr + i]);
                }
            }
        }
        return 0;
    }
};
//                        for (int j = 1; j < block_text_inx.size(); j++) {
//                            total_cmprs++;
//                            if (ptr < block_text_inx[j]) {
//                                if (r== j) {
//                                    hits++;
//                                } else {
//                                    misses++;
//                                }
//                                //std::cout << "found a block " << i << " pointing to Block " << j-1 <<  std::endl;
//                                (*bv)[j - 1] = 1;
//                                if (ptr % block_size != 0) {
//                                    (*bv)[j] = 1;
//                                }
//                                int64_t p = j - 1;
//                                pointers->push_back(p);
//                                offsets->push_back(ptr % block_size);
//                                break;
//                            }
//                        }
//        if (added_padding > 0) {
//            size_type init_block_size = (added_padding + text.size()) / this->s_;
//            this->s_ = text.size() / init_block_size;
//            if (text.size() / init_block_size != 0) {
//                this->s_++;
//            }
//            size_type padding =  this->s_ * init_block_size - text.size();
//            std::cout << "psdd " << padding << std::endl;
//            std::vector<size_type> lpf_padding(padding);
//            std::vector<size_type> lpf_ptr_padding(padding);
//            lpf_padding[0] = 0;
//            lpf_ptr_padding[0] = 0;
//            for (size_type i = 1; i < lpf_padding.size(); i++) {
//                lpf_padding[i] = lpf_padding.size() - i;
//            }
//            for (size_type i = 1; i < lpf_ptr_padding.size(); i++) {
//                lpf_ptr_padding[i] = lpf.size() + i - 1;
//            }
//            std::cout << this->s_ << std::endl;
//            lpf.insert(lpf.end(), lpf_padding.begin(), lpf_padding.end());
//            lpf_ptr.insert(lpf_ptr.end(), lpf_ptr_padding.begin(), lpf_ptr_padding.end());
//        }
//        std::cout << "after Padding: " <<  lpf.size() << std::endl;
//        int blocks = this->s_;
//            if (this->block_tree_types_.size() == 0) {
//                std::vector<size_type> top_level(this->u_chars_,0);
//                std::vector<std::vector<size_type>> temp(this->u_chars_);
//                size_type b = 0;
//                for (size_type t_index = 0; t_index < text.size() && b < block_text_inx.size(); t_index++) {
//
//                    if (block_text_inx[b] == t_index) {
//                        for (int i = 0; i < this->chars_index_.size(); i++) {
//                            auto c = this->chars_index_[i];
//                            if (c != -1) {
//                                temp[c].push_back(top_level[c]);
//                            }
//                        }
//                        b++;
//                    }
//                    top_level[this->chars_index_[text[t_index]]]++;
//                }
//                std::cout << " " << temp.size() << std::endl;
//                for (auto v: temp) {
//                    std::vector<size_type>* cv = new std::vector<size_type>;
//                    for (auto val: v) {
//                        cv->push_back(val);
//                    }
//                    this->top_level_c_ranks_.push_back(cv);
//                }
//            } else {
//                std::vector<std::vector<size_type>> ranks_per_b_c = std::vector<std::vector<size_type>>(this->u_chars_, std::vector<size_type>(block_text_inx.size()));
//                std::vector<std::vector<size_type>*>* lvl = new std::vector<std::vector<size_type>*>;
//                for (size_type i = 0; i < block_text_inx.size(); i++) {
//                    std::vector<size_type> ranks(this->u_chars_,0);
//                    size_type block_start = block_text_inx[i];
//                    for (int j = 0; j < block_size && block_start + j < text.size(); j++) {
//                        ranks[this->chars_index_[text[block_start + j]]]++;
//                    }
//                    for (auto c: this->chars_index_) {
//                        if (c != -1) {
//                            ranks_per_b_c[c].push_back(ranks[c]);
//                        }
//                    }
//                }
//                for (auto v: ranks_per_b_c) {
//                    std::vector<size_type>* cv = new std::vector<size_type>;
//                    for (auto val: v) {
//                        cv->push_back(val);
//                    }
//                    lvl->push_back(cv);
//                }
//                this->c_ranks_.push_back(lvl);
//            }
//                        int64_t b = this->find_next_smallest_index_binary_search(ptr, block_text_inx);
//                        if (b != -1 && ptr <= block_text_inx[b] + block_size - 1 && ((block_text_inx[b] + block_size ==  block_text_inx[b+1]) || ptr % block_size == 0)) {
//                            current_pointer = b;
//                            current_offset = ptr % block_size;
//                            has_ptr = true;
//                        }
//                        size_type b = this->find_next_smallest_index_binary_search(ptr, block_text_inx);
//                        if (ptr >= block_text_inx[b] && ptr <= block_text_inx[b] + block_size - 1 && ((block_text_inx[b] + block_size ==  block_text_inx[b+1]) || ptr % block_size == 0)) {
//                            current_pointer = b;
//                            current_offset = ptr % block_size;
//                            has_ptr = true;
//                        }
//                        int64_t b = find_next_smallest_index_block_tree(ptr);
//                        if (b != -1 && ptr <= block_text_inx[b] + block_size - 1 && ((block_text_inx[b] + block_size ==  block_text_inx[b+1]) || ptr % block_size == 0)) {
//                            current_pointer = b;
//                            current_offset = ptr % block_size;
//                            has_ptr = true;
//                        }
//uint64_t rank(uint8_t c, int64_t index) {
//    auto c_id = this->chars_index_[c];
////        std::cout << this->chars_index_[c] << std::endl;
//    int64_t block_size = this->block_size_lvl_[0];
//    int64_t blk_pointer = index / block_size;
//    int64_t off = index % block_size;
//    uint64_t rank = (*this->top_level_c_ranks_[c_id])[blk_pointer];
////        std::cout << "Zeichen " << index << " steht im " << blk_pointer << " Block an Stelle " << off << " vor dem Block stehen " << rank << " cahrs"<< std::endl;
//    int64_t child = 0;
//    for (int i = 0; i < this->block_tree_types_.size(); i++) {
////            if ((*this->block_tree_types_[i])[blk_pointer] == 0) {
////                std::cout << "Im " << i << "-ten Level steht es in einem BackBlock ";
////            } else {
////                std::cout << "Im " << i << "-ten Level steht es in einem Internen Block" << std::endl;
////            }
//        if ((*this->block_tree_types_[i])[blk_pointer] == 0) {
//            int64_t blk = this->block_tree_types_rs_[i]->rank0(blk_pointer);
//
//            int64_t to = off + (*this->block_tree_offsets_[i])[blk];
//            blk_pointer = (*this->block_tree_pointers_[i])[blk];
//            if (to >= block_size) {
//                blk_pointer++;
//            }
//            off = to % block_size;
//            std::cout << "und zeigt auf "<<  blk_pointer << " stelle " << off << std::endl;
//        }
//        size_type rank_blk = this->block_tree_types_rs_[i]->rank1(blk_pointer);
//        std::cout << "Es ist der " << rank_blk << "-te internal im " << i << "-ten level" << std::endl;
//        blk_pointer = rank_blk * this->tau_;
//        block_size /= this->tau_;
//        child = off / block_size;
//        off = off % block_size;
//        blk_pointer += child;
//    }
//    for (int i = 0; i < off; i++) {
//        if (this->leaves_[blk_pointer * this->leaf_size + i] == c) {
//            rank++;
//        }
//    }
//    return rank;
//}
#endif //BLOCK_TREE_BV_BLOCKTREE_LPF_H
