//
// Created by daniel on 05.09.22.
//
#include <bv_blocktree.hpp>
#include <MersenneHash.h>
#include <MersenneRabinKarp.h>
#ifndef BLOCK_TREE_BV_BLOCKTREE_FP_THEORY_H
#define BLOCK_TREE_BV_BLOCKTREE_FP_THEORY_H
__extension__ typedef unsigned __int128 uint128_t;
template<typename input_type, typename size_type>
class BV_BlockTree_fp_theory: public BV_Block_Tree<input_type,size_type> {
public:
    int32_t init(std::vector<input_type>& text) {
        static constexpr uint128_t kPrime = 2305843009213693951ULL;
        size_type added_padding = 0;
        size_type tree_max_height = 0;
        size_type max_blk_size = 0;
        this->calculate_padding(added_padding, text.size(),tree_max_height, max_blk_size);
        size_type block_size = max_blk_size;
        std::vector<size_type> block_text_inx;
        for (int i = 0; i < text.size(); i+= block_size) {
            block_text_inx.push_back(i);
        }
        while (block_size > this->max_leaf_length_) {
//            std::cout << block_size << " " << block_text_inx.size() << std::endl;

            this->block_size_lvl_.push_back(block_size);
            this->block_per_lvl_.push_back(block_text_inx.size());
            auto* bv = new pasta::BitVector(block_text_inx.size(),0);
            auto left = pasta::BitVector(block_text_inx.size(),0);
            auto right = pasta::BitVector(block_text_inx.size(),0);
            auto pair_size = 2 * block_size;
            std::unordered_map<MersenneHash<uint8_t>, std::vector<size_type>> pairs(0);
            if (pair_size > text.size()) {
                block_size /= this->tau_;
                std::vector<size_type> new_blocks(0);
                for (size_type i = 0; i < block_text_inx.size(); i++) {
                    (*bv)[i] = 1;
                    for (size_type j = 0; j < this->tau_ ; j++) {
                        if (block_text_inx[i] + (j * block_size) < text.size()) {
                            new_blocks.push_back(block_text_inx[i] + (j * block_size));
                        }
                    }
                }
                block_text_inx = new_blocks;
                this->block_tree_types_.push_back(bv);
                this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
                auto *pointers = new  sdsl::int_vector(0,0, 1);
                auto *offsets = new sdsl::int_vector(0,0, 1);
                this->block_tree_pointers_.push_back(pointers);
                this->block_tree_offsets_.push_back(offsets);
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
//            std::cout << "pairs size " <<pairs.size() << std::endl;
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
//            std::cout << "pairs size " <<pairs.size() << std::endl;
            auto old_block_size = block_size;
            auto new_block_size = block_size / this->tau_;
            std::vector<size_type> new_blocks(0);
            for (size_type i = 0; i < block_text_inx.size(); i++) {
                bool surrounded = (i > 0 && i < block_text_inx.size() - 1) && block_text_inx[i] + old_block_size == block_text_inx[i+1] && block_text_inx[i - 1] + old_block_size == block_text_inx[i];
                bool marked = false;
                if (surrounded) {
                    marked = left[i] && right[i];
                } else {
                    marked = left[i] || right[i];
                }
                if (!(marked)) {
                    (*bv)[i] = 1;
                    for (size_type j = 0; j < this->tau_ ; j++) {
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
                    MersenneRabinKarp<input_type, size_type> rk_block = MersenneRabinKarp<input_type, size_type>(text, 256, index, block_size);
                    MersenneHash<input_type> mh_block = MersenneHash<input_type>(text, rk_block.hash_, index,block_size);
                    blocks[mh_block].push_back(i);
                }
            }
//            std::cout << "hallo " << blocks.size() << std::endl;
            if (blocks.empty()) {
                block_size /= this->tau_;
                block_text_inx = new_blocks;
                this->block_tree_types_.push_back(bv);
                this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
                auto *pointers = new  sdsl::int_vector(0,0, 1);
                auto *offsets = new sdsl::int_vector(0,0, 1);
                this->block_tree_pointers_.push_back(pointers);
                this->block_tree_offsets_.push_back(offsets);
                continue;
            }
            size_type max_pointer = 0;
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

            for (size_type i = 0; i < block_text_inx.size() - 1; i++) {
                MersenneRabinKarp<input_type, size_type> rk_first_occ = MersenneRabinKarp<input_type, size_type>(text, 256, block_text_inx[i], block_size);
                bool followed = (i < block_text_inx.size() - 1) && block_text_inx[i] + old_block_size == block_text_inx[i+1] && (*bv)[i+1]==1;
                if ((*bv)[i] == 1) {
                    if (followed) {
                        for (size_type j = 0; j < block_size && block_text_inx[i] + j < text.size(); j++) {
                            MersenneHash<input_type> mh_first_occ = MersenneHash<input_type>(text, rk_first_occ.hash_, block_text_inx[i] + j,block_size);
                            if (blocks.find(mh_first_occ) != blocks.end()) {
                                for (auto b: blocks[mh_first_occ]) {
                                    pointer[b] = i;
                                    if (i > max_pointer) {
                                        max_pointer = i;
                                    }
                                    offset[b] = j;
                                }
                                blocks.erase(mh_first_occ);
                            }
                            rk_first_occ.next();

                        }
                    } else {
                        MersenneHash<input_type> mh_first_occ = MersenneHash<input_type>(text, rk_first_occ.hash_, block_text_inx[i], block_size);
                        if (blocks.find(mh_first_occ) != blocks.end()) {
                            for (auto b: blocks[mh_first_occ]) {
//                                if (b == i) {
//                                    std::cout << b << " " << i << std::endl;
//                                }
                                pointer[b] = i;
                                if (i > max_pointer) {
                                    max_pointer = i;
                                }
                                offset[b] = 0;
                            }
                            blocks.erase(mh_first_occ);
                        }
                    }
                }
            }

            auto pointer_width = 8 * sizeof(max_pointer) - this->leading_zeros(max_pointer);
//            std::cout << pointer_width << " so much space "  << max_pointer << std::endl;
            auto offset_width = 8 * sizeof(block_size) - this->leading_zeros(block_size - 1);
            auto *pointers = new  sdsl::int_vector(pointer.size(),0, pointer_width);
            auto *offsets = new sdsl::int_vector(offset.size(),0, offset_width);
            auto j = 0;
            for (size_type i = 0; i < block_text_inx.size(); i++) {
                if ((*bv)[i] == 0) {
                    if (pointer.find(i) != pointer.end()) {
                        (*pointers)[j] = pointer[i];
                    }
                    (*offsets)[j] = (size_type)offset[i];
                    j++;
                }
            }
            block_text_inx = new_blocks;
            block_size = new_block_size;
            this->block_tree_types_.push_back(bv);
            this->block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
            this->block_tree_pointers_.push_back(pointers);
            this->block_tree_offsets_.push_back(offsets);
//            mark_blocks(bv, text, block_text_inx, block_size);
        }
        this->leaf_size = block_size;
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
    BV_BlockTree_fp_theory(std::vector<input_type>& text, size_type tau, size_type max_leaf_length, size_type s) {
        this->map_unique_charas(text);
        this->tau_ = tau;
        this->max_leaf_length_ = max_leaf_length;
        this->s_ = s;
        // first we create lpf and lpf_ptr arrays;
        init(text);
    };
};
#endif //BLOCK_TREE_BV_BLOCKTREE_FP_THEORY_H
