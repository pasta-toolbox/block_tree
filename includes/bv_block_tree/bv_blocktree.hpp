
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
        bool narrate = true;
        bool debug = c == 't';
        if (narrate) {
            std::cout << "Rank query on " << index << " for " << c << std::endl;
        }
        int64_t c_index = chars_index_[c];
        size_type block_size = block_size_lvl_[0];
        size_type blk_pointer = index / block_size;
        size_type off = index % block_size;
        int64_t rank = (blk_pointer == 0) ? 0 :c_ranks_[c_index][0][blk_pointer - 1];
        size_type child = 0;
        if (narrate) {
            std::cout << "On the first level we are in block " << blk_pointer << " with offset " << off << std::endl;
            std::cout << "Before our first block we have " << rank << " " << c << "s" << std::endl;
        }
        if ((*block_tree_types_[0])[blk_pointer]) {
            std::cout << "marked block on the first level so we just select the proper child" << std::endl;
            block_size /= tau_;
            child = off / block_size;
            blk_pointer = block_tree_types_rs_[0]->rank1(blk_pointer) * tau_ + child;
            off = off % block_size;
        } else {

        }
        // we first calculate the 
        size_type i = 1;
        while (i < block_tree_types_.size()) {
            std::cout << "On the " << i << "th level we look at block " << blk_pointer << std::endl;
            if ((*block_tree_types_[i])[blk_pointer]) {
                std::cout << index << ":" << i << " " << rank << std::endl;
                rank += (child == 0) ? 0 :c_ranks_[c_index][i][blk_pointer - 1];
                size_type rank_blk = block_tree_types_rs_[i]->rank1(blk_pointer);
                blk_pointer = rank_blk * tau_;
                block_size /= tau_;
                child = off / block_size;
                off = off % block_size;
                blk_pointer += child;
                i++;
            } else {
                std::cout << "seems like we have a unmarked block we need to add our prefix " << blk_pointer << " " << c_ranks_[c_index][i].size() << std::endl;
                rank += (child == 0) ? 0 : c_ranks_[c_index][i][blk_pointer - 1];
                std::cout << "then we need to translate the query first we need - all cs b4 the offset" << std::endl;
                size_type blk = block_tree_types_rs_[i]->rank0(blk_pointer);
                rank -= pointer_c_ranks_[c_index][i][blk];
                size_type to = off + (*block_tree_offsets_[i])[blk];
                blk_pointer = (*block_tree_pointers_[i])[blk];
                std::cout << "then we need to decide which block were translateing to" << blk_pointer << std::endl;
                if (to >= block_size) {
                    auto adder = (child == 0) ? c_ranks_[c_index][i][blk_pointer] : c_ranks_[c_index][i][blk_pointer] - c_ranks_[c_index][i][blk_pointer - 1];
                    rank += adder;
                    blk_pointer++;
                }
                off = to % block_size;
                child = blk_pointer / tau_;
            }
        }
        std::cout << 122121221 << std::endl;
//        for (size_type i = 0; i < block_tree_types_.size(); i++) {
//            if ((*block_tree_types_[i])[blk_pointer]) {
//                rank += (child == 0) ? 0 :c_ranks_[c_index][0][blk_pointer - 1];
//                size_type rank_blk = block_tree_types_rs_[i]->rank1(blk_pointer);
//                blk_pointer = rank_blk * tau_;
//                block_size /= tau_;
//                child = off / block_size;
//                off = off % block_size;
//            } else {
//                size_type blk = block_tree_types_rs_[i]->rank0(blk_pointer);
//                size_type to = off + (*block_tree_offsets_[i])[blk];
//                blk_pointer = (*block_tree_pointers_[i])[blk];
//                if (to >= block_size) {
////                    if (i != 0 && blk_pointer % 2 == 1) {
////                        std::cout << "modcheck" <<  blk_pointer << std::endl;
////                    }
//                    if (c_ranks_[c_index][i][blk_pointer + 1] < c_ranks_[c_index][i][blk_pointer]) {
//                        std::cout << "WTFFFFFFFFFFFFFFFFF "<< to << std::endl;
//                    }
//                    if (debug) {
//                        std::cout << i << " 324:" <<blk_pointer%2<< " "  << c_ranks_[c_index][i][blk_pointer + 1] << " " << c_ranks_[c_index][i][blk_pointer] << std::endl;
//                    }
//                    rank = rank + c_ranks_[c_index][i][blk_pointer + 1] - c_ranks_[c_index][i][blk_pointer];
//                    blk_pointer++;
//                }
//                if (debug) {
//                    std::cout << i << " " << blk_pointer <<  " " << pointer_c_ranks_[c_index][i][blk] << " " << rank << std::endl;
//                }
//                rank -= pointer_c_ranks_[c_index][i][blk];
//                if (debug) {
//                    std::cout << i << " " << blk_pointer <<  " " << pointer_c_ranks_[c_index][i][blk] << " " << rank << std::endl;
//                }
//                off = to % block_size;
//            }
//            size_type rank_blk = block_tree_types_rs_[i]->rank1(blk_pointer);
//            blk_pointer = rank_blk * tau_;
//            block_size /= tau_;
//            child = off / block_size;
//            off = off % block_size;
//            blk_pointer += child;
//        }
//        if (debug) {
//            std::cout << "rank: " << rank << std::endl;
//        }
        std::cout << rank << " 12321321" << std::endl;
        size_type prefix_leaves = blk_pointer - child;
        for (int j = 0; j < child * leaf_size; j++) {
            if ((leaves_)[prefix_leaves + j] == c) rank++;
        }
        for (int j = 0; j <= off; j++) {
            if ((leaves_)[blk_pointer + j] == c) rank++;
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
        std::cout << "Made it until here" << std::endl;
        c_ranks_.resize(chars_.size(), std::vector<sdsl::int_vector<>>());
        pointer_c_ranks_.resize(chars_.size(), std::vector<sdsl::int_vector<>>());
        std::cout << "Made it until here" << std::endl;
        for (int i = 0; i < c_ranks_.size(); i++) {
            c_ranks_[i].resize(block_tree_types_.size(), sdsl::int_vector<>());
            for (int j = 0; j < c_ranks_[i].size(); j++) {
                c_ranks_[i][j].resize(block_tree_types_[j]->size());
            }
        }
        std::cout << "Made it until here" << std::endl;
        for (int i = 0; i < pointer_c_ranks_.size(); i++) {
            pointer_c_ranks_[i].resize(block_tree_pointers_.size(), sdsl::int_vector<>());
            for (int j = 0; j < pointer_c_ranks_[i].size(); j++) {
                pointer_c_ranks_[i][j].resize(block_tree_pointers_[j]->size());
            }
        }
        std::cout << "Made it until here" << std::endl;
        for (auto c: chars_) {
            for (size_type i = 0; i < block_tree_types_[0]->size(); i++) {
                rank_block(c, 0, i);
            }
            std::cout << c << " " << "rank dones" << std::endl;
            size_type temp1 = 0;
            size_type temp2 = 0;
            for (size_type i = 1; i < block_tree_types_[0]->size(); i++) {
                c_ranks_[chars_index_[c]][0][i] += c_ranks_[chars_index_[c]][0][i - 1];
            }
            for (size_type i = 1; i < block_tree_types_.size(); i++) {
                size_type j = 0;
                while (j < block_tree_types_[i]->size() - tau_) {
                    temp1 = c_ranks_[chars_index_[c]][i][j];
                    temp2 = 0;
                    for (size_type k = 0; k < tau_ - 1; k++) {
                        temp1 = c_ranks_[chars_index_[c]][i][j + k];
                        temp2 = c_ranks_[chars_index_[c]][i][j + k + 1];
                        if (c == 't' && i == 2) {
                            std::cout << "dieses " << temp1 << " " << temp2 << std::endl;
                        }
                        c_ranks_[chars_index_[c]][i][j + k + 1] += c_ranks_[chars_index_[c]][i][j + k];
                        temp1 = temp2;
                    };
                    j += tau_;
                }
                temp1 = c_ranks_[chars_index_[c]][i][j];
                temp2 = 0;
                while (j < block_tree_types_[i]->size() - 1) {
                    temp1 = c_ranks_[chars_index_[c]][i][j];
                    temp2 = c_ranks_[chars_index_[c]][i][j + 1];
                    c_ranks_[chars_index_[c]][i][j + 1] += c_ranks_[chars_index_[c]][i][j];
                    temp1 = temp2;
                    j++;
                }
            }
        }
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
        if (j >= block_tree_types_[i]->size()) {
            return 0;
        }
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
                rank_c -= rank_g;
                rank_c += rank_2nd;
            }
            pointer_c_ranks_[chars_index_[c]][i][rank_0] = rank_g;
        }
        c_ranks_[chars_index_[c]][i][j] = rank_c;
        return rank_c;
    }
    size_type part_rank_block(input_type c, size_type i, size_type j, size_type g) {
        if (j >= block_tree_types_[i]->size()) {
            return 0;
        }
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
