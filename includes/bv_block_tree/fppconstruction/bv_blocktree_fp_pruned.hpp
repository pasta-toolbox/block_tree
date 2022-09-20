//
// Created by daniel on 07.09.22.
//
#include <bv_blocktree.hpp>
#include <MersenneRabinKarp.h>
#include <MersenneHash.h>
#ifndef BLOCK_TREE_BV_BLOCKTREE_FP_PRUNED_H
#define BLOCK_TREE_BV_BLOCKTREE_FP_PRUNED_H
__extension__ typedef unsigned __int128 uint128_t;
#endif //BLOCK_TREE_BV_BLOCKTREE_FP_PRUNED_H
template<typename input_type, typename size_type>
class BV_BlockTree_fp_pruned: public BV_Block_Tree<input_type,size_type> {
public:
    int32_t init(std::vector<input_type>& text) {
        std::cout << "HI" << std::endl;
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
        std::cout << "HI1" << std::endl;
        this->calculate_padding(added_padding, text.size(),tree_max_height, max_blk_size);
        std::cout << "HI2" << std::endl;
        auto is_padded = added_padding > 0 ? 1 : 0;
        std::cout << "HI3" << std::endl;
        int64_t block_size = max_blk_size;
        std::cout << "HI4" << std::endl;
        std::vector<int64_t> block_text_inx;
        std::cout << "HI5 " << block_size << std::endl;
        for (int64_t i = 0; i < text.size(); i+= block_size) {
            block_text_inx.push_back(i);
        }
        std::cout << "HI6" << std::endl;
        while (block_size > this->max_leaf_length_) {
            std::cout << block_size << std::endl;
            this->block_size_lvl_.push_back(block_size);
            this->block_per_lvl_.push_back(block_text_inx.size());
            auto* bv = new pasta::BitVector(block_text_inx.size(),false);
            auto left = pasta::BitVector(block_text_inx.size(),false);
            auto right = pasta::BitVector(block_text_inx.size(),false);
            auto pair_size = 2 * block_size;
            auto last_block_padded = block_text_inx[block_text_inx.size() - 1] + block_size != text.size() ? 1 : 0;
            std::unordered_map<MersenneHash<uint8_t>, std::vector<size_type>> pairs(0);
            std::unordered_map<MersenneHash<uint8_t>, std::vector<size_type>> blocks = std::unordered_map<MersenneHash<uint8_t>, std::vector<size_type>>();
            for (size_type i = 0; i < block_text_inx.size() - last_block_padded; i++) {
                auto index = block_text_inx[i];
                MersenneRabinKarp<input_type, size_type> rk_block = MersenneRabinKarp<input_type, size_type>(text, 256, index, block_size);
                MersenneHash<input_type> mh_block = MersenneHash<input_type>(text, rk_block.hash_, index,block_size);
                blocks[mh_block].push_back(i);
            }
            std::vector<size_type> pointers(block_text_inx.size(),-1);
            std::vector<size_type> offsets(block_text_inx.size(),0);
            if (pair_size > text.size()) {
                block_size /= this->tau_;
                std::vector<int64_t> new_blocks(0);
                for (size_type i = 0; i < block_text_inx.size(); i++) {
                    (*bv)[i] = 1;
                    for (size_type j = 0; j < this->tau_ ; j++) {
                        if (block_text_inx[i] + (j * block_size) < text.size()) {
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
            for (size_type i = 0; i < block_text_inx.size() - 1; i++) {
                if (block_text_inx[i] + block_size == block_text_inx[i+1] && block_text_inx[i] + pair_size <= text.size()) {
                    auto index = block_text_inx[i];
                    MersenneRabinKarp<input_type, size_type> rk_pair = MersenneRabinKarp<input_type, size_type>(text, 256, index, pair_size);
                    MersenneHash<input_type> mh_pair = MersenneHash<input_type>(text, rk_pair.hash_, index,pair_size);
                    pairs[mh_pair].push_back(i);
                }
            }
            // find pairs
            MersenneRabinKarp<input_type, size_type> rk_pair_sw = MersenneRabinKarp<input_type, size_type>(text, 256, 0, pair_size);
            for (size_type i = 0; i < text.size() - pair_size; i++) {
                MersenneHash<input_type> mh_sw = MersenneHash<input_type>(text, rk_pair_sw.hash_, i, pair_size);
                if (pairs.find(mh_sw) != pairs.end()) {
                    for (auto b: pairs[mh_sw]) {
                        if (i != block_text_inx[b]) {
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
            for (size_type i = 0; i < block_text_inx.size(); i++) {
                bool surrounded = (i > 0 && i < block_text_inx.size() - 1) && block_text_inx[i] + old_block_size == block_text_inx[i+1] && block_text_inx[i - 1] + old_block_size == block_text_inx[i];
                bool marked = false;
                if (surrounded) {
                    marked = left[i] && right[i];
                } else {
                    marked = left[i] || right[i];
                }
                if (!(marked) || block_text_inx[i] + old_block_size >= text.size()) {
                    (*bv)[i] = 1;
                    for (size_type j = 0; j < this->tau_ ; j++) {
                        if (block_text_inx[i] + (j * new_block_size) < text.size()) {
                            new_blocks.push_back(block_text_inx[i] + (j * new_block_size));
                        }
                    }
                }
            }
            for (size_type i = 0; i < block_text_inx.size() - 1; i++) {
                MersenneRabinKarp<input_type, size_type> rk_first_occ = MersenneRabinKarp<input_type, size_type>(text, 256, block_text_inx[i], block_size);
                bool followed = (i < block_text_inx.size() - 1) && block_text_inx[i] + block_size == block_text_inx[i+1] && (*bv)[i+1]==1;
                if ((*bv)[i] == 1) {
                    if (followed) {
                        for (size_type j = 0; j < block_size && block_text_inx[i] + j + block_size < text.size(); j++) {
                            MersenneHash<input_type> mh_first_occ = MersenneHash<input_type>(text, rk_first_occ.hash_, block_text_inx[i] + j,block_size);
                            if (blocks.find(mh_first_occ) != blocks.end()) {
                                for (auto b: blocks[mh_first_occ]) {
//                                    if (i == 340) {
//                                        std::cout << block_size << " " << b << std::endl;
//                                    }
                                    if (b != i) {
                                        pointers[b] = i;
                                        offsets[b] = j;
                                    }
                                }
                                blocks.erase(mh_first_occ);
                            }
                            rk_first_occ.next();
                        }
                    } else {
                        MersenneHash<input_type> mh_first_occ = MersenneHash<input_type>(text, rk_first_occ.hash_, block_text_inx[i], block_size);
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
            blk_lvl.push_back(block_text_inx);
            block_text_inx = new_blocks;
            block_size = new_block_size;
            bv_pass_1.push_back(bv);
//            if (block_size <= this->max_leaf_length_ ) {
//                std::cout << *bv << std::endl << std::endl;
//            }

        }

        this->leaf_size = block_size;
        block_size *= this->tau_;
        std::cout << "step 1 done" << std::endl;
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
//                if (ptr == 340) {
//                    std::cout<< i << " " << (*bv)[j]  << " marks 340: " << j << " " << off << std::endl;
//                }
//                if (j == 2318) {
//                    std::cout<< i << " " << (*bv)[j]  << " marks 340: " << j << " " << ptr << std::endl;
//                }
                if (ptr == -1 || ptr + std::min(off, 1) >= j ) {
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
//                if (j == 340) {
//                    std::cout<< i << " " << (*bv)[j]  << std::endl;
//                }
            }
            pass2_max_offset.push_back(block_size);
            block_size *= this->tau_;
            bv_pass_2.push_back(bv);
//            if (i == bv_pass_1.size() - 1) {
//                std::cout << *bv << std::endl;
//            }
            size_type ones_per_pass2 = 0;
            for (size_type j = 0; j < bv->size(); j++) {
                ones_per_pass2 += (*bv)[j];
            }
            pass2_ones.push_back(ones_per_pass2);
            pass2_pointer.push_back(pointers);
            pass2_offset.push_back(offsets);
            pass2_max_pointer.push_back(max_pointer);
        }
        std::cout << "step 2 done" << std::endl;
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
        std::cout << "step 2 done" << std::endl;
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
        std::cout << "leaf done" << std::endl;
        for (auto bv : bv_pass_1) {
            delete bv;
        }
        for (size_type i = 0; i < bv_pass_2.size() - 1; i++) {
            delete bv_pass_2[i];
        }
        return 0;
    };
    BV_BlockTree_fp_pruned(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, size_type s) {
        this->map_unique_chars(text);
    this->tau_ = tau;
    this->max_leaf_length_ = max_leaf_length;
    this->s_ = s;
    std::cout << this->tau_ << " " << this->s_ << std::endl;
    init(text);
    };
    ~BV_BlockTree_fp_pruned() {
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
        for (auto a: this->top_level_c_ranks_) {
            delete a;
        }
        for (auto a: this->c_ranks_) {
            for (auto b: *a) {
                delete b;
            }
            delete a;
        }
        for (auto a: this->g_ranks_) {
            for (auto b: *a) {
                delete b;
            }
            delete a;
        }

    };
};