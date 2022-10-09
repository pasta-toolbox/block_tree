
#include <iostream>
#include <vector>
#include <bit_vector.hpp>
#include <support/find_l2_flat_with.hpp>
#include <support/flat_rank.hpp>
#include <support/flat_rank_select.hpp>
#include <support/optimized_for.hpp>
#include <support/rank.hpp>
#include <support/rank_select.hpp>
#include <support/wide_rank.hpp>
#include <support/wide_rank_select.hpp>
#ifndef BLOCK_TREE_BV_BLOCKTREE_HPP
#define BLOCK_TREE_BV_BLOCKTREE_HPP

template<typename input_type, typename size_type>
class BV_Block_Tree {
public:
    size_type tau_;
    size_type max_leaf_length_;
    size_type s_ = 1;
    size_type leaf_size = 0;
    size_type amount_of_leaves = 0;
    std::vector<pasta::BitVector*> block_tree_types_;
    std::vector<pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>*> block_tree_types_rs_;
    std::vector<sdsl::int_vector<>*> block_tree_pointers_;
    std::vector<sdsl::int_vector<>*> block_tree_offsets_;
    std::vector<int64_t> block_size_lvl_;
    std::vector<int64_t> block_per_lvl_;
    std::vector<input_type> leaves_;
    std::unordered_map<input_type, size_type> chars_index_;
    std::vector<input_type> chars_;
    size_type u_chars_;
    std::vector<std::vector<int64_t>> top_level_c_ranks_;
    std::vector<std::vector<sdsl::int_vector<>>> c_ranks_;
    std::vector<std::vector<sdsl::int_vector<>>> pointer_c_ranks_;

    int64_t access(size_type index) {
        size_type block_size = block_size_lvl_[0];
        size_type blk_pointer = index / block_size;
        size_type off = index % block_size;
        size_type child = 0;
        for (size_type i = 0; i < block_tree_types_.size(); i++) {
            if ((*block_tree_types_[i])[blk_pointer] == 0) {
                size_type blk = block_tree_types_rs_[i]->rank0(blk_pointer);
                size_type to = off + (*block_tree_offsets_[i])[blk];
                blk_pointer = (*block_tree_pointers_[i])[blk];
                if (to >= block_size) {
                    blk_pointer++;
                }
                off = to % block_size;
            }
            size_type rank_blk = block_tree_types_rs_[i]->rank1(blk_pointer);
            blk_pointer = rank_blk * tau_;
            block_size /= tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer += child;
        }
        return leaves_[blk_pointer * leaf_size + off];
    };
    int64_t rank(input_type c, size_type index) {
        int64_t c_index = chars_index_[c];
        size_type block_size = block_size_lvl_[0];
        size_type blk_pointer = index / block_size;
        size_type off = index % block_size;
        int64_t rank = top_level_c_ranks_[c_index][blk_pointer];
        size_type child = 0;
        if ((*block_tree_types_[0])[blk_pointer] == 0) {
            size_type blk = block_tree_types_rs_[0]->rank0(blk_pointer);
            size_type to = off + (*block_tree_offsets_[0])[blk];
        }
        for (size_type i = 0; i < block_tree_types_.size(); i++) {
            if ((*block_tree_types_[i])[blk_pointer]) {

            } else {
                size_type blk = block_tree_types_rs_[i]->rank0(blk_pointer);
                size_type to = off + (*block_tree_offsets_[i])[blk];
                blk_pointer = (*block_tree_pointers_[i])[blk];
                if (to >= block_size) {
                    blk_pointer++;
                }
                off = to % block_size;
            }
            size_type rank_blk = block_tree_types_rs_[i]->rank1(blk_pointer);
            blk_pointer = rank_blk * tau_;
            block_size /= tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer += child;
        }
        return rank;
    };
    size_type print_space_usage() {
        size_type space_usage = sizeof(tau_) + sizeof(max_leaf_length_) + sizeof(s_) + sizeof(leaf_size);
        for (auto bv: block_tree_types_) {
            space_usage += bv->size()/8;
        }
        for (auto rs: block_tree_types_rs_) {
            space_usage += rs->space_usage();
        }
        for (auto iv: block_tree_pointers_) {
            space_usage += (iv->size() * iv->width())/8;
        }
        for (auto iv: block_tree_offsets_) {
            space_usage += (iv->size() * iv->width())/8;
        }
        for (auto v: block_size_lvl_) {
            space_usage += sizeof(v);
        }
        for (auto v: block_per_lvl_) {
            space_usage += sizeof(v);
        }
        space_usage += leaves_.size() * sizeof(input_type);
//        for (auto v: top_level_c_ranks_) {
//            space_usage += v->size() * sizeof(size_type);
//        }
//        for (auto v: c_ranks_) {
//            for (auto iv: *v) {
//                space_usage += iv->size() * sizeof(size_type);
//            }
//        }
//        for (auto v: g_ranks_) {
//            for (auto iv: *v) {
//                space_usage += iv->size() * sizeof(size_type);
//            }
//        }
//        std::cout << "we use " << space_usage << "Bytes" << std::endl;
        return space_usage;
    };
    int32_t add_rank_support() {
        c_ranks_.resize(chars_.size(), std::vector<sdsl::int_vector<>>());
        pointer_c_ranks_.resize(chars_.size(), std::vector<sdsl::int_vector<>>());
        top_level_c_ranks_.resize(chars_.size(), std::vector<int64_t>(block_per_lvl_[0],0));
        for (int i = 0; i < c_ranks_.size(); i++) {
            c_ranks_[i].resize(block_tree_types_.size(), sdsl::int_vector<>());
            for (int j = 0; j < c_ranks_[i].size(); j++) {
                c_ranks_[i][j].resize(block_tree_types_[j]->size());
            }
        }
        for (int i = 0; i < pointer_c_ranks_.size(); i++) {
            pointer_c_ranks_[i].resize(block_tree_pointers_.size(), sdsl::int_vector<>());
            for (int j = 0; j < pointer_c_ranks_[i].size(); j++) {
                pointer_c_ranks_[i][j].resize(block_tree_pointers_[j]->size());
            }
        }
        for (auto c: chars_) {
            for (size_type i = 0; i < block_tree_types_[0]->size(); i++) {
                rank_block(c, 0, i);
            }
            std::cout << "In here?" << std::endl;
            size_type temp1 = c_ranks_[chars_index_[c]][0][0];
            size_type temp2 = 0;
            for (size_type i = 0; i < block_tree_types_[0]->size() - 1; i++) {
                temp2 = c_ranks_[chars_index_[c]][0][i + 1];
                c_ranks_[chars_index_[c]][0][i + 1] += c_ranks_[chars_index_[c]][0][i];
                c_ranks_[chars_index_[c]][0][i] -= temp1;
                temp1 = temp2;
            }
            for (size_type i = 1; i < block_tree_types_.size(); i++) {
                size_type j = 0;
                std::cout << i << std::endl;
                while (j < block_tree_types_[i]->size() - tau_) {
                    temp1 = 0;
                    temp2 = 0;
                    for (size_type k = 0; k < tau_ - 1; k++) {

                        temp1 = c_ranks_[chars_index_[c]][i][j + k];
                        temp2 = c_ranks_[chars_index_[c]][i][j + k + 1];
                        c_ranks_[chars_index_[c]][i][j + k + 1] += c_ranks_[chars_index_[c]][i][j + k];
                        c_ranks_[chars_index_[c]][i][j + k] -= temp1;
                        temp1 = temp2;
                    }
                    j += tau_;
                }
                std::cout << " test " << std::endl;
                temp1 = 0;
                temp2 = 0;
                while (j < block_tree_types_[i]->size() - 1) {
                    temp1 = c_ranks_[chars_index_[c]][i][j];
                    temp2 = c_ranks_[chars_index_[c]][i][j + 1];
                    c_ranks_[chars_index_[c]][i][j + 1] += c_ranks_[chars_index_[c]][i][j];
                    c_ranks_[chars_index_[c]][i][j] -= temp1;
                    temp1 = temp2;
                    j++;
                }
            }
            std::cout << "Not in here?" << std::endl;
        }
        for (auto c: chars_) {
            top_level_c_ranks_[chars_index_[c]][0] = 0;
            for (size_type j = 1; j < c_ranks_[chars_index_[c]][0].size(); j++) {
                if (c_ranks_[chars_index_[c]][0][j - 1] < 0) {
                    std::cout << c << " " << c_ranks_[chars_index_[c]][0][j - 1] << std::endl;
                }
                std::cout << c << " " << j << " " << top_level_c_ranks_[chars_index_[c]][j - 1] << " " << c_ranks_[chars_index_[c]][0][j - 1] << std::endl;
                top_level_c_ranks_[chars_index_[c]][j] = top_level_c_ranks_[chars_index_[c]][j - 1] + c_ranks_[chars_index_[c]][0][j - 1];
                std::cout << top_level_c_ranks_[chars_index_[c]][j] << std::endl;
//                top_level_c_ranks_[chars_index_[c]][j] += (size_type) c_ranks_[chars_index_[c]][0][j - 1];
            }
        }
        std::vector<std::vector<size_type>> leaf_ranks;
//        c_ranks_.resize(u_chars_, )
//        for (auto c: chars_) {
//            size_type rank_c = 0;
//            for (size_type i = 0; i < block_tree_types_[0]->size(); i++) {
//                rank_c += rank_block(c,0,i);
//            }
//            std::cout << c << ": " << rank_c << std::endl;
//            for (int i = 0; i < c_ranks_[chars_index_[c]].size(); i++) {
//                for (int j = 0; j < c_ranks_[chars_index_[c]][i].size(); j++) {
//                    std::cout << c_ranks_[chars_index_[c]][i][j] << " ";
//                }
//                std::cout << std::endl;
//            }
//            for (int i = 0; i < pointer_c_ranks_[chars_index_[c]].size(); i++) {
//                for (int j = 0; j < pointer_c_ranks_[chars_index_[c]][i].size(); j++) {
//                    std::cout << pointer_c_ranks_[chars_index_[c]][i][j] << " ";
//                }
//                std::cout << std::endl;
//            }
//        }
//        for (int64_t i = 0; i < amount_of_leaves; i++) {
//            leaf_ranks[chars_index_[c]][i]  = rank_leaf(c,i, leaf_size);
//        }
//        std::vector<std::vector<std::vector<size_type>>> ranks;
//        std::vector<std::vector<std::vector<size_type>>> pointer_ranks;
//        for (int i = block_tree_types_.size() - 1; i >= 0; i++) {
//            for (int j = 0; j < block_tree_types_[i]->size(); j++) {
//                size_type child = j % tau_;
//                std::vector<size_type> accumulator = std::vector<size_type>(u_chars_, 0);
//                for (int k = 0; k < tau_; k++) {
////                    for (auto c: chars_) {
////                        ranks[i][j][chars_index_[c]]
////                    }
//                }
//
//                // ranks as child
//                if (!(*block_tree_types_[i])[j]) {
//                    // ranks internal
//                }
//            }
//        }

        return 0;
    }
    inline size_type leading_zeros(int32_t val) {
        return __builtin_clz(static_cast<unsigned int>(val) | 1);
    }
    inline size_type leading_zeros(int64_t val) {
        return __builtin_clzll(static_cast<unsigned long long>(val) | 1);
    }
    void calculate_padding(int64_t &padding, int64_t text_length, int64_t &height, int64_t &blk_size) {
        int64_t tmp_padding = this->s_ * this->tau_;
        int64_t h = 1;
        blk_size = tau_;
        while (tmp_padding < text_length) {
            std::cout << tmp_padding << std::endl;
            tmp_padding *= this->tau_;
            blk_size *= this->tau_;
            h++;
        }
        height = h;
        padding = tmp_padding - text_length;
//        std::cout << "Padding: " << padding << " h: " << h << " SIZE: " << tmp_padding << " BLK_SIZE: " <<  blk_size <<   std::endl;
    }
    size_type rank_block(input_type c, size_type i, size_type j) {
        size_type rank_c = 0;
        if ((*block_tree_types_[i])[j] == 1) {
            if (i != block_tree_types_.size() - 1) {
                size_type rank_blk = block_tree_types_rs_[i]->rank1(j);
                for (size_type k = 0; k < tau_; k++) {
                    rank_c += rank_block(c, i + 1, rank_blk * tau_ + k);
                }
            } else {
                size_type rank_blk = block_tree_types_rs_[i]->rank1(j);
                for (size_type k = 0; k < tau_; k++) {
                    rank_c += rank_leaf(c, rank_blk * tau_ + k, leaf_size);
                }
            }
        } else {
            size_type rank_0 = block_tree_types_rs_[i]->rank0(j);
            size_type ptr = (*block_tree_pointers_[i])[rank_0];
            size_type off = (*block_tree_offsets_[i])[rank_0];
            size_type rank_g = 0;
            rank_c += c_ranks_[chars_index_[c]][i][ptr];
            if (off != 0) {
                rank_g = part_rank_block(c,i, ptr, off - 1);
                size_type rank_2nd = part_rank_block(c,i,ptr + 1, off - 1);
                if (c == 'N')
                    std::cout << i << ":" << j << "ptr " << ptr << " rank_g: " << rank_g << " rank_2nd: " << rank_2nd <<  " rankc: " << c_ranks_[chars_index_[c]][i][ptr] <<std::endl;
                if (rank_c < rank_g) {
                    std::cout << "wtf" << std::endl;
                }
                rank_c -= rank_g;
                rank_c += rank_2nd;
            }
            pointer_c_ranks_[chars_index_[c]][i][rank_0] = rank_g;
        }
        c_ranks_[chars_index_[c]][i][j] = rank_c;
        return rank_c;
    }
    size_type part_rank_block(input_type c, size_type i, size_type j, size_type g) {
        bool debug = false;
        size_type rank_c = 0;
        size_type blk_size = block_size_lvl_[i];
        if ((*block_tree_types_[i])[j] == 1) {
            if (i != block_tree_types_.size() - 1) {
                size_type rank_blk = block_tree_types_rs_[i]->rank1(j);
                size_type k = 0;
                size_type k_sum = 0;
                for (k = 0; k < tau_ && k_sum + block_size_lvl_[i + 1] <= g + 1; k++) {
                    rank_c += c_ranks_[chars_index_[c]][i + 1][rank_blk * tau_ + k];
                    k_sum+=block_size_lvl_[i + 1];
                    if (debug) {
                        std::cout << k << " " << rank_c << std::endl;
                    }
                }
                if (debug) {
                    std::cout << i + 1 << " " << rank_blk * tau_ + k << " " << g % block_size_lvl_[i + 1] << std::endl;
                    std::cout << i + 1 << " " << c_ranks_[chars_index_[c]][i + 1].size() << std::endl;
                }
                if (k_sum != g + 1) {
                    rank_c += part_rank_block(c, i + 1, rank_blk * tau_ + k, g % block_size_lvl_[i + 1]);
                }
                if (debug) {
                    std::cout << k << " " << rank_c << std::endl;
                }
            } else {
                size_type rank_blk = block_tree_types_rs_[i]->rank1(j);
                size_type k = 0;
                size_type k_sum = 0;
                for (k = 0; k < tau_ && k_sum + leaf_size <= g + 1; k++) {
                    rank_c += rank_leaf(c, rank_blk * tau_ + k, leaf_size);
                    k_sum += leaf_size;
                }
                if (k_sum != g + 1) {
                    rank_c += rank_leaf(c, rank_blk * tau_ + k, g % leaf_size);
                }
            }
        } else {
            size_type rank_0 = block_tree_types_rs_[i]->rank0(j);
            size_type ptr = (*block_tree_pointers_[i])[rank_0];
            size_type off = (*block_tree_offsets_[i])[rank_0];
            if (g + off > block_size_lvl_[i]) {
                rank_c += c_ranks_[chars_index_[c]][i][ptr] - pointer_c_ranks_[chars_index_[c]][i][rank_0] +
                          part_rank_block(c, i, ptr + 1, g - off - 1);
            } else {
                rank_c += part_rank_block(c, i, ptr, g + off - 1) - pointer_c_ranks_[chars_index_[c]][i][rank_0];
            }
        }
        return rank_c;
    }
    size_type rank_leaf(input_type c, size_type leaf_index,size_type i) {

        if (leaf_index * leaf_size >= leaves_.size()) {
            return 0;
        }
//        size_type x = leaves_.size() - leaf_index * this->tau_;
//        i = std::min(i, x);
        size_type result = 0;
        for (size_type ind = 0; ind < i; ind++) {
            if (leaves_[leaf_index * leaf_size + ind] == c) {
                result++;
            }
        }
        return result;
    }

    size_type map_unique_chars(std::vector<input_type> &text) {
        this->u_chars_ = 0;
        input_type i = 0;
        for (auto a: text) {
            if (chars_index_.find(a) == chars_index_.end()) {
                chars_index_[a] = i;
                i++;
                chars_.push_back(a);
            }
        }
        this->u_chars_ = i;
        return 0;
    };
    size_type find_next_smallest_index_binary_search(size_type i, std::vector<int64_t>& pVector) {
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
    };
    int64_t find_next_smallest_index_linear_scan(size_type i, std::vector<size_type>& pVector) {
        int64_t b = 0;
        while (b < pVector.size() && i >= pVector[b]) {
            b++;
        }
        return b - 1;
    };
    size_type find_next_smallest_index_block_tree(size_type index) {
        size_type block_size = this->block_size_lvl_[0];
        size_type blk_pointer = index / block_size;
        size_type off = index % block_size;
        size_type child = 0;
        for (size_type i = 0; i < this->block_tree_types_.size(); i++) {
            if ((*this->block_tree_types_[i])[blk_pointer] == 0) {
                return -1;
            }
            if (off > 0 && (*this->block_tree_types_[i])[blk_pointer + 1] == 0) {
                return -1;
            }
            size_type rank_blk = this->block_tree_types_rs_[i]->rank1(blk_pointer);
            blk_pointer = rank_blk * this->tau_;
            block_size /= this->tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer += child;
        }
        return blk_pointer;
    };
};
#endif //BLOCK_TREE_BV_BLOCKTREE_HPP
