
#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/find_l2_flat_with.hpp>
#include <pasta/bit_vector/support/flat_rank.hpp>
#include <pasta/bit_vector/support/flat_rank_select.hpp>
#include <pasta/bit_vector/support/optimized_for.hpp>
#include <pasta/bit_vector/support/rank.hpp>
#include <pasta/bit_vector/support/rank_select.hpp>
#include <pasta/bit_vector/support/wide_rank.hpp>
#include <pasta/bit_vector/support/wide_rank_select.hpp>
#include <chrono>
#ifndef BLOCK_TREE_BV_BLOCKTREE_LPF_H
#define BLOCK_TREE_BV_BLOCKTREE_LPF_H
class BV_BlockTree_lpf_64 {
public:
    int64_t tau_;
    int64_t max_leaf_length_;
    int64_t s_;
    std::vector<pasta::BitVector*> block_tree_types_;
    std::vector<pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>*> block_tree_types_rs_;
    std::vector<std::vector<int64_t>*> block_tree_pointers_;
    std::vector<std::vector<int64_t>*> block_tree_offsets_;
    std::vector<int64_t> block_size_lvl_;
    std::vector<int64_t> leaf_pointer_;
    std::vector<u_int8_t> leaves_;
    BV_BlockTree_lpf_64(std::vector<uint8_t>& text, int32_t tau, int32_t max_leaf_length, int64_t s) {
        s_ = s;
        tau_ = tau;
        max_leaf_length_ = max_leaf_length;
        // first we create lpf and lpfptr arrays;
        std::vector<int64_t> lpf(text.size());
        std::vector<int64_t> lpf_ptr(text.size());
        lpf_array64(text, lpf, lpf_ptr);
        int64_t z = 0;
        int64_t added_padding = 0;
        calculate_lz_factor(s_,lpf);
        s_ += 2;
        calculate_padding(added_padding, text.size());
        if (added_padding > 0) {
            int64_t init_block_size = (added_padding + text.size()) / s_;
            s_ = std::ceil(text.size() / init_block_size);
            int64_t padding = init_block_size - (text.size() % s_);
            std::vector<int64_t> lpf_padding(padding);
            std::vector<int64_t> lpf_ptr_padding(padding);
            lpf_padding[0] = 0;
            lpf_ptr_padding[0] = 0;
            for (int64_t i = 1; i < lpf_padding.size(); i++) {
                lpf_padding[i] = lpf_padding.size() - i;
            }
            for (int64_t i = 1; i < lpf_ptr_padding.size(); i++) {
                lpf_ptr_padding[i] = lpf.size() + i - 1;
            }
            lpf.insert(lpf.end(), lpf_padding.begin(), lpf_padding.end());
            lpf_ptr.insert(lpf_ptr.end(), lpf_ptr_padding.begin(), lpf_ptr_padding.end());
        }
        std::cout << "after Padding: " <<  lpf.size() << std::endl;
        int blocks = s_;
        int block_size = lpf.size() / blocks;
        std::vector<int64_t> block_text_inx(blocks);
        for (int i = 0; i < blocks; i++) {
            block_text_inx[i] = (i * block_size);
        }
        //manipulate lpf_ptr
        for(int64_t i = 0; i < lpf_ptr.size(); i++) {
            if (lpf[i] <= lpf[lpf_ptr[i]]) {
                lpf_ptr[i] = lpf_ptr[lpf_ptr[i]];
            }
        }
        auto t01 = std::chrono::high_resolution_clock::now();
        int64_t total_chains = 0;
        int64_t total_cmprs = 0;
        int64_t total_blocks = 0;
        int64_t hits = 0;
        int64_t misses = 0;
        while (block_size > max_leaf_length_) {
            std::cout << "Blockssize: " << block_size << " #Blocks: " << block_text_inx.size();
            block_size_lvl_.push_back(block_size);
            pasta::BitVector* bv = new pasta::BitVector(block_text_inx.size(),0);
            std::vector<int64_t>* pointers = new std::vector<int64_t>();
            std::vector<int64_t>* offsets = new std::vector<int64_t>();
            std::vector<int64_t> block_text_inx_new;
            for(int i = block_text_inx.size() - 1; i >=  0; i--) {
                total_blocks++;
                // if block is already marked we don't have to check for replacements
                if ((*bv)[i] == 1) {
                    continue;
                }
                // ind is initially block i's starting pos in the text
                int64_t ind = block_text_inx[i];
                int64_t ur_ind = ind;
                bool has_ptr = false;
                // we look until we find leftmost occ of [ind..ind+blocksize]
                while (lpf[ind] >= block_size) {
                    if (lpf[lpf_ptr[ind]] >= block_size) {
                        lpf_ptr[ind] = lpf_ptr[lpf_ptr[ind]];
//                        lpf_ptr[ur_ind] = lpf_ptr[lpf_ptr[ind]];
                    }
                    if (block_size + lpf_ptr[ind] - 1 >= ur_ind || lpf[lpf_ptr[ind]] >= block_size) {
//                        std::cout << ind << " I" << text[ind] << " ?" << std::endl;
                        // our found occs overlapps with the block
                        // we will try if we find another one

                        ind = lpf_ptr[ind];
                        total_chains++;
                    } else {
                        // we found the leftmost occ
                        has_ptr = true;
                        // determine which block(s) it points at
                        int64_t ptr = lpf_ptr[ind];
                        lpf_ptr[ur_ind] = ptr;
                        int64_t block_in_lvl = ptr / block_size_lvl_[0];
                        int64_t off_in_lvl = off_in_lvl % block_size_lvl_[0];
                        int64_t child = off_in_lvl / (block_size_lvl_[0]/tau_);
                        for (int j  = 1; j < block_size_lvl_.size() - 1; j++) {
                            block_in_lvl = (int64_t)(*block_tree_types_rs_[j]).rank1(block_in_lvl ) ;
                            block_in_lvl += child;
                            off_in_lvl = off_in_lvl % block_size_lvl_[j];
                            child = off_in_lvl / block_size_lvl_[j+1];
                        }
                        int64_t l = 0;
                        int64_t r = i;
                        while (l < r) {
                            total_cmprs++;
                            int64_t m = std::floor((l+r)/2);
                            if (ptr < block_text_inx[m]) {
                                r = m;
                            } else {
                                l = m + 1;
                            }
                        }
                        if (block_in_lvl + child == r-1) {
                            hits++;
                        } else {
                            misses++;
                        }
                        (*bv)[r - 1] = 1;
                        if (ptr % block_size != 0) {
                            (*bv)[r] = 1;
                        }
                        int64_t p = r - 1;
                        pointers->push_back(p);
                        offsets->push_back(ptr % block_size);
                        break;
                    }
                }
                // if there is no block its pointing at it is its leftmost occ
                if (!has_ptr) {
                    (*bv)[i] = 1;
                }
            }
            block_size = block_size / tau_;
            for (int i = 0; i < bv->size(); i++) {
                if ((*bv)[i] == 1) {
                    for (int j = 0; j < tau_; j++) {
                        block_text_inx_new.push_back(block_text_inx[i] + j * block_size);

                    }
                }
            }
            block_tree_types_.push_back(bv);
            block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
            block_tree_pointers_.push_back(pointers);
            block_tree_offsets_.push_back(offsets);
            block_text_inx = block_text_inx_new;

            auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t01);
            t01 = std::chrono::high_resolution_clock::now();
            std::cout << " Time: " << ms_int.count() << std::endl;
        }
        std::cout << "Total chains: " << total_chains << " Total cmprs: " << total_cmprs <<  " by total blocks: " << total_blocks << std::endl;
        std::cout << "hits " << hits << " misses " << misses << std::endl;
        leaf_pointer_ = block_text_inx;
        for (auto ptr: leaf_pointer_) {
            for (int i = 0; i < block_size; i++) {
                if (ptr + i < text.size()) {
                    leaves_.push_back(text[ptr + i]);
                } else {
                    leaves_.push_back(0);
                }

            }
        }
    };
    uint8_t access(int64_t index) {

        int64_t blk_pointer = index / block_size_lvl_[0];
        int64_t off = index % block_size_lvl_[0];
        std::cout << "access: " << index << " in " << blk_pointer<< " " <<off << std::endl;
        for (int i = 0; i < block_tree_types_.size() - 1; i++) {
            if ((*block_tree_types_[i])[blk_pointer] == 1) {
                int rank_blk = block_tree_types_rs_[i]->rank1(blk_pointer);
                std::cout << i << "-lvl: we are the " << rank_blk << " marked block" <<std::endl;
                if (off < block_size_lvl_[i]/tau_) {
                    std::cout << "the first child which is block " << rank_blk * 2 << std::endl;
                    blk_pointer = rank_blk * 2;
                } else {
                    std::cout << "the 2nd child which is block " << rank_blk * 2 + 1<< " off is" << off % block_size_lvl_[i]/2 << std::endl;
                    blk_pointer = rank_blk * 2 + 1;
                    off = off % block_size_lvl_[i]/2;
                }
            } else {
                int rank_blk = block_tree_types_rs_[i]->rank1(blk_pointer);
                off = off + (*block_tree_offsets_[i])[block_tree_offsets_[i]->size() -1 + blk_pointer];
                blk_pointer = (*block_tree_pointers_[i])[block_tree_offsets_[i]->size() -1 + blk_pointer];
                if (off >= block_size_lvl_[i]) {
                    off = off % block_size_lvl_[i];
                    blk_pointer++;
                }
                if (off < block_size_lvl_[i]/2) {
                    blk_pointer = block_tree_types_rs_[i]->rank1(blk_pointer) * 2;
                } else {
                    blk_pointer = block_tree_types_rs_[i]->rank1(blk_pointer) * 2 + 1;
                    off = off % block_size_lvl_[i]/2;
                }
            }
        }

        int rank_blk = block_tree_types_rs_[block_tree_types_rs_.size() -1]->rank1(blk_pointer);
        std::cout << rank_blk << std::endl;
        return leaves_[leaf_pointer_[rank_blk] + off];
    }
private:
    void calculate_lz_factor(int64_t &z, std::vector<int64_t> &lpf) {
        int64_t temp_z = 0;
        int64_t lz = 0;
        while (lz < lpf.size()) {
            lz = lz + std::max(1LL, lpf[lz + 1]);
            temp_z++;
        }
        z = temp_z;
        std::cout << "Given Text has " << z << " LZ-factors"<< std::endl;
    }
    void calculate_padding(int64_t &padding, int64_t text_length) {
        int64_t tmp_padding = s_ * tau_;
        int64_t h = 1;
        while (tmp_padding < text_length) {
            tmp_padding *= tau_;
            h++;
        }
        padding = tmp_padding - text_length;
        std::cout << "Padding: " << padding << " h: " << h << " SIZE: " << tmp_padding <<  std::endl;
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
#endif //BLOCK_TREE_BV_BLOCKTREE_LPF_H
