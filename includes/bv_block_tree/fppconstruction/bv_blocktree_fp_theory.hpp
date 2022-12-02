//
// Created by daniel on 05.09.22.
//
#include <bv_blocktree.hpp>
#include <MersenneHash.h>
#include <MersenneRabinKarp.h>
#include <limits.h>

#ifndef BLOCK_TREE_BV_BLOCKTREE_FP_THEORY_H
#define BLOCK_TREE_BV_BLOCKTREE_FP_THEORY_H
__extension__ typedef unsigned __int128 uint128_t;

template<typename input_type, typename size_type>
class BV_BlockTree_fp_theory : public BV_Block_Tree<input_type, size_type> {
public:
    size_type sigma_ = 0;

    int32_t init(std::vector<input_type> &text) {
        static constexpr uint128_t kPrime = 2305843009213693951ULL;
        int64_t added_padding = 0;
        int64_t tree_max_height = 0;
        int64_t max_blk_size = 0;
        this->calculate_padding(added_padding, text.size(), tree_max_height, max_blk_size);
        int64_t block_size = max_blk_size;
        std::vector<int64_t> block_text_inx;
        for (int64_t i = 0; i < text.size(); i += block_size) {
            block_text_inx.push_back(i);
        }
        bool found_back_block = false;
        while (block_size > this->max_leaf_length_) {
            auto *bv = new pasta::BitVector(block_text_inx.size(), 0);
            auto left = pasta::BitVector(block_text_inx.size(), 0);
            auto right = pasta::BitVector(block_text_inx.size(), 0);
            auto pair_size = 2 * block_size;
            std::unordered_map<MersenneHash<uint8_t>, std::vector<size_type>> pairs(0);
            if (pair_size > text.size()) {
                block_size /= this->tau_;
                std::vector<int64_t> new_blocks(0);
                for (size_type i = 0; i < block_text_inx.size(); i++) {
                    (*bv)[i] = 1;
                    for (size_type j = 0; j < this->tau_; j++) {
                        if (block_text_inx[i] + (j * block_size) < text.size()) {
                            new_blocks.push_back(block_text_inx[i] + (j * block_size));
                        }
                    }
                }
                block_text_inx = new_blocks;
                if (!this->CUT_FIRST_LEVELS) {
                    this->block_tree_types_.push_back(bv);
                    this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
                    auto *pointers = new sdsl::int_vector(0, 0, 1);
                    auto *offsets = new sdsl::int_vector(0, 0, 1);
                    this->block_tree_pointers_.push_back(pointers);
                    this->block_tree_offsets_.push_back(offsets);
                    this->block_size_lvl_.push_back(block_size * this->tau_);
                } else {
                    delete bv;
                }
                continue;
            }
            for (size_type i = 0; i < block_text_inx.size() - 1; i++) {
                if (block_text_inx[i] + block_size == block_text_inx[i + 1] &&
                    block_text_inx[i] + pair_size <= text.size()) {
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
            MersenneRabinKarp<input_type, size_type> rk_pair_sw = MersenneRabinKarp<input_type, size_type>(text, sigma_,
                                                                                                           0,
                                                                                                           pair_size,
                                                                                                           kPrime);
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
            size_type ones_per_level = 0;
            for (int64_t i = 0; i < block_text_inx.size(); i++) {
                bool surrounded = (i > 0 && i < block_text_inx.size() - 1) &&
                                  block_text_inx[i] + old_block_size == block_text_inx[i + 1] &&
                                  block_text_inx[i - 1] + old_block_size == block_text_inx[i];
                bool marked = false;
                if (surrounded) {
                    marked = left[i] && right[i];
                } else {
                    marked = left[i] || right[i];
                }
                if (!(marked)) {
                    (*bv)[i] = 1;
                    ones_per_level++;
                    for (size_type j = 0; j < this->tau_; j++) {
                        if (block_text_inx[i] + (j * new_block_size) < text.size()) {
                            new_blocks.push_back(block_text_inx[i] + (j * new_block_size));
                        }
                    }
                }
            }
            std::unordered_map<MersenneHash<uint8_t>, std::vector<size_type>> blocks(0);
            std::unordered_map<size_type, size_type> pointer(0);
            std::unordered_map<size_type, size_type> offset(0);

            for (size_type i = 0; i < block_text_inx.size(); i++) {
                if ((*bv)[i] == 0) {
                    auto index = block_text_inx[i];
                    MersenneRabinKarp<input_type, size_type> rk_block = MersenneRabinKarp<input_type, size_type>(text,
                                                                                                                 sigma_,
                                                                                                                 index,
                                                                                                                 block_size,
                                                                                                                 kPrime);
                    MersenneHash<input_type> mh_block = MersenneHash<input_type>(text, rk_block.hash_, index,
                                                                                 block_size);
                    blocks[mh_block].push_back(i);
                }
            }
//            if (blocks.empty()) {
//                block_size /= this->tau_;
//                block_text_inx = new_blocks;
//                continue;
//            }

//            for (size_type i = 0; i < text.size() - block_size; i++) {
//                MersenneHash<input_type> mh_first_occ = MersenneHash<input_type>(text, rk_first_occ.hash_, i,block_size);
//                if (pairs.find(mh_first_occ) != pairs.end()) {
//                    for (auto b: pairs[mh_first_occ]) {
//                        pointer[b] = i;
//                        if (i > max_pointer) {
//                            max_pointer = i;
//                        }
//                        offset[b] = i % block_size;
//                    }
//                    pairs.erase(mh_first_occ);
//                }
//                rk_first_occ.next();
//            }
            if (!blocks.empty()) {
                MersenneRabinKarp<input_type, size_type> rk_first_occ = MersenneRabinKarp<input_type, size_type>(
                        text, sigma_, block_text_inx[0], block_size, kPrime);
                for (size_type i = 0; i < block_text_inx.size() - 1; i++) {
                    bool followed = (i < block_text_inx.size() - 1) &&
                                    block_text_inx[i] + old_block_size == block_text_inx[i + 1] && (*bv)[i + 1] == 1;
                    if ((*bv)[i] == 1) {
                        if (followed) {
                            if (rk_first_occ.init_ != block_text_inx[i]) {
                                rk_first_occ.restart(block_text_inx[i]);
                            }
                            for (size_type j = 0; j < block_size && block_text_inx[i] + j < text.size(); j++) {
                                MersenneHash<input_type> mh_first_occ = MersenneHash<input_type>(text,
                                                                                                 rk_first_occ.hash_,
                                                                                                 block_text_inx[i] + j,
                                                                                                 block_size);
                                if (blocks.find(mh_first_occ) != blocks.end()) {
                                    for (auto b: blocks[mh_first_occ]) {
                                        pointer[b] = i;
                                        offset[b] = j;
                                    }
                                    blocks.erase(mh_first_occ);
                                }
                                rk_first_occ.next();

                            }
                        } else {
                            if (rk_first_occ.init_ != block_text_inx[i]) {
                                rk_first_occ.restart(block_text_inx[i]);
                            }
                            MersenneHash<input_type> mh_first_occ = MersenneHash<input_type>(text, rk_first_occ.hash_,
                                                                                             block_text_inx[i],
                                                                                             block_size);
                            if (blocks.find(mh_first_occ) != blocks.end()) {
                                for (auto b: blocks[mh_first_occ]) {
                                    pointer[b] = i;
                                    offset[b] = 0;
                                }
                                blocks.erase(mh_first_occ);
                            }
                        }
                    }
                }
            }
            auto *pointers = new sdsl::int_vector(pointer.size(), 0, 64);
            auto *offsets = new sdsl::int_vector(offset.size(), 0, 64);
            auto j = 0;
            for (size_type i = 0; i < block_text_inx.size(); i++) {
                if ((*bv)[i] == 0) {
                    if (pointer.find(i) != pointer.end()) {
                        (*pointers)[j] = pointer[i];
                    }
                    (*offsets)[j] = (size_type) offset[i];
                    j++;
                }
            }

            found_back_block |= ones_per_level != bv->size();
            if (found_back_block || !this->CUT_FIRST_LEVELS) {
                this->block_tree_types_.push_back(bv);
                this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
                sdsl::util::bit_compress(*pointers);
                sdsl::util::bit_compress(*offsets);
                this->block_tree_pointers_.push_back(pointers);
                this->block_tree_offsets_.push_back(offsets);
                this->block_size_lvl_.push_back(block_size);
            } else {
                delete bv;
                delete pointers;
                delete offsets;
            }
            block_text_inx = new_blocks;
            block_size = new_block_size;
//            mark_blocks(bv, text, block_text_inx, block_size);
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
    }

    ~BV_BlockTree_fp_theory() {
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

    BV_BlockTree_fp_theory(std::vector<input_type> &text, size_type tau, size_type max_leaf_length, size_type s,
                           size_type sigma, bool cut_first_levels) {
        this->CUT_FIRST_LEVELS = cut_first_levels;
        sigma_ = sigma;
        this->map_unique_chars(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        this->s_ = s;
        // first we create lpf and lpf_ptr arrays;
        init(text);
    };
};

#endif //BLOCK_TREE_BV_BLOCKTREE_FP_THEORY_H
