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
    int64_t leaf_size = 0;
    std::vector<pasta::BitVector*> block_tree_types_;
    std::vector<pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>*> block_tree_types_rs_;
    std::vector<std::vector<int64_t>> block_tree_pointers_;
    std::vector<std::vector<int64_t>> block_tree_offsets_;
    std::vector<int64_t> block_size_lvl_;
    std::vector<int64_t> block_per_lvl_;
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
        int64_t new_blocks = 0;
        while (block_size > max_leaf_length_) {
            block_size_lvl_.push_back(block_size);
            block_per_lvl_.push_back(block_text_inx.size());
            std::cout << "Blockssize: " << block_size << " #Blocks: " << block_text_inx.size();
            pasta::BitVector* bv = new pasta::BitVector(block_text_inx.size(),0);
            std::vector<int64_t> pointers = std::vector<int64_t>();
            std::vector<int64_t> offsets = std::vector<int64_t>();

            for(int i = block_text_inx.size() - 1; i >=  0; i--) {
                if ((*bv)[i] == 1) {
                    continue;
                }
                int64_t ind = block_text_inx[i];
                int64_t ur_ind = ind;
                bool has_ptr = false;
                while (lpf[ind] >= block_size) {
                    if (lpf[lpf_ptr[ind]] >= block_size) {
                        lpf_ptr[ind] = lpf_ptr[lpf_ptr[ind]];
                    }
                    if (block_size + lpf_ptr[ind] - 1 >= ur_ind || lpf[lpf_ptr[ind]] >= block_size) {
                        ind = lpf_ptr[ind];
                        total_chains++;
                    } else {
                        has_ptr = true;
                        int64_t ptr = lpf_ptr[ind];
                        lpf_ptr[ur_ind] = ptr;
                        int64_t l = 0;
                        int64_t r = i;
//                        int id = ptr;
//                        int bs = block_size_lvl_[0];
//                        int block = id/bs;
//                        int off = id % bs;
//                        for (int j = 1; j < block_size_lvl_.size(); j++) {
//                            block = block_tree_types_rs_[j - 1]->rank1(block) * tau_ + off/block_size_lvl_[j];
//                            off %= block_size_lvl_[j];
//                        }
                        while (l < r) {
                            total_cmprs++;
                            int64_t m = std::floor((l+r)/2);
                            if (ptr < block_text_inx[m]) {
                                r = m;
                            } else {
                                l = m + 1;
                            }
                        }
//                        if (block_text_inx[r] < ptr || ptr > block_text_inx[r] + block_size || block_text_inx[r + 1] < ptr) {
//                            std::cout << "FEHLER" << std::endl;
//                        }
//                        if (block + 1 == r) {
//                            hits++;
//                        } else {
//                            std::cout << std::endl << block_text_inx[block + 1] << " " << block_text_inx[r] << " " << ptr << std::endl;
//                            misses++;
//                        }
//                        std::cout << "From " << i << " to " << r  - 1<< std::endl;
                        (*bv)[r - 1] = 1;
                        if (ptr % block_size != 0) {
                            (*bv)[r] = 1;
//                            std::cout << "From " << i << " to " << r << std::endl;
                        }
                        int64_t p = r - 1;
                        pointers.push_back(p);
                        offsets.push_back(ptr % block_size);
                        break;
                    }
                }
                // if there is no block its pointing at it, is its leftmost occ
                if (!has_ptr) {
                    (*bv)[i] = 1;
                }
            }
            std::vector<int64_t> block_text_inx_new(0,0);
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
            std::vector<int64_t> p(pointers.size());
            std::vector<int64_t> o(offsets.size());
            for(int j = 0; j < p.size(); j++) {
                p[j] = pointers[pointers.size() - 1 - j];
                o[j] = offsets[pointers.size() - 1 - j];
            }
            block_tree_pointers_.push_back(p);
            block_tree_offsets_.push_back(o);
            block_text_inx = block_text_inx_new;
            auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t01);
            t01 = std::chrono::high_resolution_clock::now();
            std::cout << " Time: " << ms_int.count() << " H/M: "<< hits << " " << misses << std::endl;
        }
        std::cout << "Total chains: " << total_chains << " Total cmprs: " << total_cmprs <<  " by total blocks: " << total_blocks << std::endl;
        std::cout << "hits " << hits << " misses " << misses << std::endl;
        leaf_size = block_size;
        for (int64_t ptr: block_text_inx) {
            for (int i = 0; i < block_size; i++) {
                if (ptr + i < text.size()) {
//                    std::cout << "blk " << block_size << " i " << i << " ptr " << ptr <<  std::endl;
                    leaves_.push_back(text[ptr + i]);
                } else {
                    leaves_.push_back(0);
                }

            }
        }
        for (int i = 0; i < block_tree_types_.size(); i++) {
            print_bv(i);
            for(int j = 0; j < block_tree_pointers_[i].size(); j++) {
                std::cout << block_tree_pointers_[i][j] <<",";
            }
            std::cout << std::endl;
            for(int j = 0; j < block_tree_offsets_[i].size(); j++) {
                std::cout << block_tree_offsets_[i][j] <<",";
            }
            std::cout << std::endl;
        }
        print_leafs();
    };
    uint8_t access(int64_t index) {
        int64_t block_size = block_size_lvl_[0];
        // select the top level block and it's offset
        int64_t blk_pointer = index / block_size;
        int64_t off = index % block_size;
        int64_t child = 0;
//        std::cout << "access: " << index << " in " << blk_pointer<< " " <<off << std::endl;
        for (int i = 0; i < block_tree_types_.size(); i++) {
            if ((*block_tree_types_[i])[blk_pointer] == 0) {
                int64_t blk = block_tree_types_rs_[i]->rank0(blk_pointer);
                int64_t to = off + block_tree_offsets_[i][blk];
                std::cout << "1: " <<blk_pointer << std::endl;
                blk_pointer = block_tree_pointers_[i][blk];
                if (to >= block_size) {
                    blk_pointer++;
                }
                off = to % block_size;
                std::cout << "2: " << blk_pointer << std::endl;
            }
            int64_t rank_blk = block_tree_types_rs_[i]->rank1(blk_pointer);
            blk_pointer = rank_blk * tau_;
            block_size /= tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer += child;
//            std::cout << "Blocksize  " << block_size_lvl_[i] << ": " << blk_pointer;
//            std::cout << " we find our " << index/(block_size_lvl_[i]/tau_) << "th child @ ";


//            std::cout << "ptr " << blk_pointer << " rnk " << rank_blk << std::endl;

//            std::cout << blk_pointer <<"/" << off << std::endl;
        }
//        std::cout << leaf_pointer_[blk_pointer] << std::endl;
        return leaves_[blk_pointer * leaf_size + off];
    }
private:
    void calculate_lz_factor(int64_t &z, std::vector<int64_t> &lpf) {
        int64_t temp_z = 0;
        int64_t lz = 0;
        while (lz < lpf.size()) 
	{
	    if (lpf[lz +1] > 1) {
		    lz = lz + lpf[lz + 1];
		} else {
			lz++;
		}
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
    int32_t print_bv(int i) {
        for (int j = 0; j < block_per_lvl_[i]; j++) {
            std::cout << (*block_tree_types_[i])[j];
        }
        std::cout << std::endl;
        return 0;
    }
    int32_t print_leafs() {
        for (int j = 0; j <leaves_.size(); j++) {
            std::cout << leaves_[j];
        }
        std::cout << std::endl;
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
#endif //BLOCK_TREE_BV_BLOCKTREE_LPF_H
