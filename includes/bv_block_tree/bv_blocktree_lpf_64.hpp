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
#include <sdsl/int_vector.hpp>

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
    std::vector<sdsl::int_vector<>*> block_tree_pointers_;
    std::vector<sdsl::int_vector<>*> block_tree_offsets_;
    std::vector<sdsl::int_vector<>*> first_lvl_ranks;
    std::vector<int64_t> block_size_lvl_;
    std::vector<int64_t> block_per_lvl_;
    std::vector<u_int8_t> leaves_;
    std::vector<int32_t> chars_index_;
    int32_t u_chars_;
    std::vector<std::vector<uint32_t>*> top_level_c_ranks_;
    std::vector<std::vector<std::vector<uint32_t>*>*> c_ranks_;
    std::vector<std::vector<std::vector<uint32_t>*>*> g_ranks_;
    int64_t find_next_smallest_index_ls(int64_t i, std::vector<int64_t>& pVector) {
        int64_t b = 0;
        while (b < pVector.size() && i >= pVector[b]) {
            b++;
        }
        return b - 1;
    }
    int64_t find_next_smallest_index_bt(int64_t index) {
        int64_t block_size = block_size_lvl_[0];
        int64_t blk_pointer = index / block_size;
        int64_t off = index % block_size;
        int64_t child = 0;
        for (int i = 0; i < block_tree_types_.size(); i++) {
            if ((*block_tree_types_[i])[blk_pointer] == 0) {
                return -1;
            }
            if (off > 0 && (*block_tree_types_[i])[blk_pointer + 1] == 0) {
                return -1;
            }
            int64_t rank_blk = block_tree_types_rs_[i]->rank1(blk_pointer);

//            std::cout << "Es ist der " << rank_blk << "-te internal im " << i << "-ten level" << std::endl;
            blk_pointer = rank_blk * tau_;
            block_size /= tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer += child;
        }
        return blk_pointer;
    }
    int64_t find_next_smallest_index_bs(int64_t i, std::vector<int64_t>& pVector) {
        int64_t l = 0;
        int64_t r = pVector.size();
        while (l < r) {
            int64_t m = std::floor((l+r)/2);
            if (i < pVector[m]) {
                r = m;
            } else {
                l = m + 1;
            }
        }
        return r -1 ;
    }

    BV_BlockTree_lpf_64(std::vector<uint8_t>& text, int32_t tau, int32_t max_leaf_length, int64_t s) {
        chars_index_ = std::vector<int32_t>(256, -1);
        u_chars_ = 0;
        uint8_t i = 0;
        for (auto a: text) {
            if (chars_index_[a] == -1 ) {
                chars_index_[a] = i;
                i++;
            }
        }
        u_chars_ = i;
        for (auto a: chars_index_) {
            if (a != -1) {
                std::cout << a << std::endl;
            }
        }
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
        calculate_padding(added_padding, text.size());
        if (added_padding > 0) {
            int64_t init_block_size = (added_padding + text.size()) / s_;
            s_ = text.size() / init_block_size;
            if (text.size() / init_block_size != 0) {
                s_++;
            }
            int64_t padding =  s_ * init_block_size - text.size();
            std::cout << "psdd " << padding << std::endl;
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
            std::cout << s_ << std::endl;
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
//        for(int64_t i = 0; i < lpf_ptr.size(); i++) {
//            if (lpf[i] <= lpf[lpf_ptr[i]] && lpf_ptr[i] == i-1) {
//                lpf_ptr[i] = lpf_ptr[lpf_ptr[i]];
//            }
//        }
        auto t01 = std::chrono::high_resolution_clock::now();
        int64_t total_chains = 0;
        int64_t total_cmprs = 0;
        int64_t total_blocks = 0;
        int64_t hits = 0;
        int64_t misses = 0;
        int64_t new_blocks = 0;
        while (block_size > max_leaf_length_) {
            uint64_t max_pointer = 0;
            uint64_t max_offset = 0;
            std::cout << "Blocksize: " << block_size << " Blocks in Level: " << block_text_inx.size();
            block_size_lvl_.push_back(block_size);
            block_per_lvl_.push_back(block_text_inx.size());
            pasta::BitVector* bv = new pasta::BitVector(block_text_inx.size(),0);
            std::vector<int64_t> pointers = std::vector<int64_t>();
            std::vector<int64_t> offsets = std::vector<int64_t>();
            for(int block_in_lvl = block_text_inx.size() - 1; block_in_lvl >= 0; block_in_lvl--) {
                if ((*bv)[block_in_lvl] == 1) {
                    continue;
                }
                int64_t initial_index = block_text_inx[block_in_lvl];
                int64_t ind = initial_index;
                int64_t current_pointer = 0;
                int64_t current_offset = 0;
                bool has_ptr = false;
                while (lpf[ind] >= block_size) {
                    int64_t ptr = lpf_ptr[ind];
                    if (block_size + lpf_ptr[ind] - 1 >= initial_index) {
                        ind = lpf_ptr[ind];
                    } else if (lpf[lpf_ptr[ind]] >= block_size ) {
                        ind = lpf_ptr[ind];
//                        int64_t b = find_next_smallest_index_bt(ptr);
//                        if (b != -1 && ptr <= block_text_inx[b] + block_size - 1 && ((block_text_inx[b] + block_size ==  block_text_inx[b+1]) || ptr % block_size == 0)) {
//                            current_pointer = b;
//                            current_offset = ptr % block_size;
//                            has_ptr = true;
//                        }
                        int64_t b = find_next_smallest_index_bs(ptr, block_text_inx);
                        if (ptr <= block_text_inx[b] + block_size - 1 && (block_text_inx[b] + block_size ==  block_text_inx[b+1])) {
                            current_pointer = b;
                            current_offset = ptr % block_size;
                            has_ptr = true;
                        }
                    } else {
//                        int64_t b = find_next_smallest_index_bt(ptr);
//                        if (b != -1 && ptr <= block_text_inx[b] + block_size - 1 && ((block_text_inx[b] + block_size ==  block_text_inx[b+1]) || ptr % block_size == 0)) {
//                            current_pointer = b;
//                            current_offset = ptr % block_size;
//                            has_ptr = true;
//                        }
                        int64_t b = find_next_smallest_index_bs(ptr, block_text_inx);
                        if (ptr <= block_text_inx[b] + block_size - 1 && ((block_text_inx[b] + block_size ==  block_text_inx[b+1]) || ptr % block_size == 0)) {
                            current_pointer = b;
                            current_offset = ptr % block_size;
                            has_ptr = true;
                        }
                        if (has_ptr) {
                            (*bv)[current_pointer] = 1;
                            if (ptr != block_text_inx[b]) {
                                (*bv)[current_pointer + 1] = 1;
                            }
                            pointers.push_back(current_pointer);
                            offsets.push_back(current_offset);

                            if (current_pointer > max_pointer) {
                                max_pointer = current_pointer;
                            }
                            if (current_offset > max_offset) {
                                max_offset = current_offset;
                            }
                        }
                        break;
                    }
                }
                // if there is no block its pointing at it, is its leftmost occ
                if (!has_ptr) {
                    (*bv)[block_in_lvl] = 1;
                }
            }
            if (block_tree_types_.size() == 0) {
                std::vector<uint64_t> top_level(u_chars_,0);
                std::vector<std::vector<uint64_t>> temp(u_chars_);
                int64_t b = 0;
                for (int64_t t_index = 0; t_index < text.size() && b < block_text_inx.size(); t_index++) {

                    if (block_text_inx[b] == t_index) {
                        for (int i = 0; i < chars_index_.size(); i++) {
                            auto c = chars_index_[i];
                            if (c != -1) {
                                temp[c].push_back(top_level[c]);
                            }
                        }
                        b++;
                    }
                    top_level[chars_index_[text[t_index]]]++;
                }
                std::cout << " " << temp.size() << std::endl;
                for (auto v: temp) {
                    std::vector<uint32_t>* cv = new std::vector<uint32_t>;
                    for (auto val: v) {
                        cv->push_back(val);
                    }
                    top_level_c_ranks_.push_back(cv);
                }
            } else {
                std::vector<std::vector<uint64_t>> ranks_per_b_c = std::vector<std::vector<uint64_t>>(u_chars_, std::vector<uint64_t>(block_text_inx.size()));
                std::vector<std::vector<uint32_t>*>* lvl = new std::vector<std::vector<uint32_t>*>;
                for (int64_t i = 0; i < block_text_inx.size(); i++) {
                    std::vector<uint64_t> ranks(u_chars_,0);
                    int64_t block_start = block_text_inx[i];
                    for (int j = 0; j < block_size && block_start + j < text.size(); j++) {
                        ranks[chars_index_[text[block_start + j]]]++;
                    }
                    for (auto c: chars_index_) {
                        if (c != -1) {
                            ranks_per_b_c[c].push_back(ranks[c]);
                        }
                    }
                }
                for (auto v: ranks_per_b_c) {
                    std::vector<uint32_t>* cv = new std::vector<uint32_t>;
                    for (auto val: v) {
                        cv->push_back(val);
                    }
                    lvl->push_back(cv);
                }
                c_ranks_.push_back(lvl);
            }
            std::vector<int64_t> block_text_inx_new(0,0);
            block_size = block_size / tau_;
            for (int i = 0; i < bv->size(); i++) {
                if ((*bv)[i] == 1) {
                    for (int j = 0; j < tau_; j++) {
                        block_text_inx_new.push_back(block_text_inx[i] + (j * block_size));
                    }
                }
            }

            block_tree_types_.push_back(bv);
            block_tree_types_rs_.push_back(new pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>(*bv));
            int ptr_size = 64 - __builtin_clzll(max_pointer | 1);
            int off_size = 64 - __builtin_clzll(max_offset | 1);
            sdsl::int_vector<>*p = new sdsl::int_vector(pointers.size(), 0,ptr_size);
            sdsl::int_vector<>*o = new sdsl::int_vector(offsets.size(),0,off_size);
            std::vector<std::vector<uint64_t>> g_temp(u_chars_);
            for(int j = 0; j < pointers.size(); j++) {
                auto pointer = pointers[pointers.size() - 1 - j];
                auto offset = offsets[pointers.size() - 1 - j];
                (*p)[j] = pointer;
                (*o)[j] = offset;
                std::vector<uint64_t> g_counter(u_chars_);
                for (int i = 0; i < offset; i++) {
                    g_counter[chars_index_[text[block_text_inx[pointer] + i]]]++;
                }
                for (int i = 0; i < g_counter.size(); i++) {
                    g_temp[i].push_back(g_counter[i]);
                }
            }
            std::vector<std::vector<uint32_t>*>* lvl = new std::vector<std::vector<uint32_t>*>;
            for (auto v: g_temp) {
                std::vector<uint32_t>* cv = new std::vector<uint32_t>;
                for (auto val: v) {
                    cv->push_back(val);
                }
                lvl->push_back(cv);
            }
            g_ranks_.push_back(lvl);
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
                    leaves_.push_back(text[ptr + i]);
                }
            }
        }
//        for (int i = 0; i < block_tree_types_.size(); i++) {
//            print_bv(i);
//            for(int j = 0; j < block_tree_pointers_[i].size(); j++) {
//                std::cout << block_tree_pointers_[i][j] <<",";
//            }
//            std::cout << std::endl;
//            for(int j = 0; j < block_tree_offsets_[i].size(); j++) {
//                std::cout << block_tree_offsets_[i][j] <<",";
//            }
//            std::cout << std::endl;
//        }
//        print_leafs();
        std::cout << leaves_.size() << " so much space" << std::endl;
    };
    ~BV_BlockTree_lpf_64() {
        for (auto a : block_tree_types_) {
            delete a;
        }
        for (auto a: block_tree_types_rs_) {
            delete a;
        }
        for (auto a: block_tree_pointers_) {
            delete a;
        }
        for (auto a: block_tree_offsets_) {
            delete a;
        }
        for (auto a: top_level_c_ranks_) {
            delete a;
        }
        for (auto a: c_ranks_) {
            for (auto b: *a) {
                delete b;
            }
            delete a;
        }
        for (auto a: g_ranks_) {
            for (auto b: *a) {
                delete b;
            }
            delete a;
        }

    };
    uint8_t access(int64_t index) {

        int64_t block_size = block_size_lvl_[0];
        int64_t blk_pointer = index / block_size;
        int64_t off = index % block_size;
//        std::cout << "Zeichen " << index << " steht im " << blk_pointer << " Block an Stelle " << off << std::endl;
        int64_t child = 0;
        for (int i = 0; i < block_tree_types_.size(); i++) {
//            if ((*block_tree_types_[i])[blk_pointer] == 0) {
//                std::cout << "Im " << i << "-ten Level steht es in einem BackBlock ";
//            } else {
//                std::cout << "Im " << i << "-ten Level steht es in einem Internen Block" << std::endl;
//            }
            if ((*block_tree_types_[i])[blk_pointer] == 0) {
                int64_t blk = block_tree_types_rs_[i]->rank0(blk_pointer);

                int64_t to = off + (*block_tree_offsets_[i])[blk];
                blk_pointer = (*block_tree_pointers_[i])[blk];
                if (to >= block_size) {
                    blk_pointer++;
                }
                off = to % block_size;
//                std::cout << "und zeigt auf "<<  blk_pointer << " stelle " << off << std::endl;
            }
            int64_t rank_blk = block_tree_types_rs_[i]->rank1(blk_pointer);
//            std::cout << "Es ist der " << rank_blk << "-te internal im " << i << "-ten level" << std::endl;
            blk_pointer = rank_blk * tau_;
            block_size /= tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer += child;
        }
        return leaves_[blk_pointer * leaf_size + off];
    }
    uint64_t rank(uint8_t c, int64_t index) {
        auto c_id = chars_index_[c];
        std::cout << chars_index_[c] << std::endl;
        int64_t block_size = block_size_lvl_[0];
        int64_t blk_pointer = index / block_size;
        int64_t off = index % block_size;
        uint64_t rank = (*top_level_c_ranks_[c_id])[blk_pointer];
        std::cout << "Zeichen " << index << " steht im " << blk_pointer << " Block an Stelle " << off << " vor dem Block stehen " << rank << " cahrs"<< std::endl;
        int64_t child = 0;
        for (int i = 0; i < block_tree_types_.size(); i++) {
            if ((*block_tree_types_[i])[blk_pointer] == 0) {
                std::cout << "Im " << i << "-ten Level steht es in einem BackBlock ";
            } else {
                std::cout << "Im " << i << "-ten Level steht es in einem Internen Block" << std::endl;
            }
            if ((*block_tree_types_[i])[blk_pointer] == 0) {
                int64_t blk = block_tree_types_rs_[i]->rank0(blk_pointer);

                int64_t to = off + (*block_tree_offsets_[i])[blk];
                blk_pointer = (*block_tree_pointers_[i])[blk];
                if (to >= block_size) {
                    blk_pointer++;
                }
                off = to % block_size;
                std::cout << "und zeigt auf "<<  blk_pointer << " stelle " << off << std::endl;
            }
            int64_t rank_blk = block_tree_types_rs_[i]->rank1(blk_pointer);
            std::cout << "Es ist der " << rank_blk << "-te internal im " << i << "-ten level" << std::endl;
            blk_pointer = rank_blk * tau_;
            block_size /= tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer += child;
        }
        for (int i = 0; i < off; i++) {
            if (leaves_[blk_pointer * leaf_size + i] == c) {
                rank++;
            }
        }
        return rank;
    }
private:
    void calculate_lz_factor(int64_t &z, std::vector<int64_t> &lpf) {
        int64_t temp_z = 0;
        int64_t lz = 0;
        while (lz < lpf.size() - 1)
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
