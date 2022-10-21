#include <iostream>
#include "libsais.h"
#include <lpfarray.hpp>
#include <vector>
#include "bit_vector.hpp"
#include "cmdline_parser.hpp"
#include <fstream>
#include <sstream>
#include <lpfconstruction/bv_blocktree_lpf_heuristic.hpp>
#include <lpfconstruction/bv_blocktree_lpf_pruned.hpp>
#include <lpfconstruction/bv_blocktree_lpf_theory.hpp>
#include <fppconstruction/bv_blocktree_fp_theory.hpp>
#include <fppconstruction/bv_blocktree_fp_pruned.hpp>
#include <chrono>
#include <type_traits>
#include <iostream>
#include <unordered_set>
int main(int argc, char* argv[]) {
    tlx::CmdlineParser cp;
    // add a byte size argument which the user can enter like '1gi'
    uint64_t a_size = 0;
    cp.add_bytes('s', "size", a_size,
                 "Number of bytes to process.");

    // process command line
    if (!cp.process(argc, argv))
        return -1; // some error occurred and help was always written to user.

//    std::cout << "Command line parsed okay." << std::endl;
    std::cout <<  std::endl  << "Run with "<< a_size << " Bytes" << std::endl;
    std::string test(a_size, ' ');
//    std::ifstream t("/home/daniel/blocktree-experiments/data/Escherichia_Coli");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/english");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/einstein.de.txt");
    std::ifstream t("/home/daniel/blocktree-experiments/data/einstein.en.txt");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/influenza");
//    std::ifstream t("/Users/daniel/Downloads/einstein.en.txt");
    std::stringstream buffer;
    t.read(&test[0], a_size);
//    test = "NNBOBOTWNNBOBIOOTBSHTFNEBOBOTWNEBOBOTWNEBOBIOOTBSHTFNSBOBOTW";
//   test = "ABCDABCDEFGHABCDDE12";
    std::vector<uint8_t> vec(test.begin(), test.end());
    auto t01 = std::chrono::high_resolution_clock::now();

//    BV_BlockTree_fp_theory<uint8_t, int32_t> *fp_bt = new BV_BlockTree_fp_theory<uint8_t, int32_t>(vec,2,1,1);

//    BV_BlockTree_fp_pruned<uint8_t, int32_t> *fp_bt = new BV_BlockTree_fp_pruned<uint8_t, int32_t>(vec, 16, 4, 1);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << "time " << ms_int.count() << std::endl;
    std::vector<int32_t> lpf(vec.size());
    std::vector<int32_t> lpf_ptr(vec.size());
    std::vector<int32_t> lz;
    int32_t lzn = 0;
    lpf_array(vec, lpf, lpf_ptr);
    calculate_lz_factor(lzn,lpf, lz);
    BV_BlockTree_fp_pruned<uint8_t, int32_t>*  fp_bt = new BV_BlockTree_fp_pruned<uint8_t, int32_t>(vec, 2, 8,1);
    BV_BlockTree_lpf_pruned<uint8_t, int32_t>*  lpf_bt = new BV_BlockTree_lpf_pruned<uint8_t, int32_t>(vec, 2, 8,lpf, lpf_ptr, lz);
    std::cout << "Errors space " << lpf_bt->print_space_usage() << std::endl;
    int  j = 0;
    for (int i = 0; i < vec.size(); i++) {
        auto x = lpf_bt->access(i);
        if (x != vec[i]) {

            j++;
        }

    }
    std::cout << "Errors " << j << std::endl;
    std::cout << "lpf" << std::endl;
    lpf_bt->add_rank_support();
    fp_bt->add_rank_support();
    for (auto c: lpf_bt->chars_) {
        std::cout << c << ":";
        for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_[c]].size(); i++) {
            std::cout << (int) lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][i].width() << "/" << (int) lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][i].size() << " ";
        }
        std::cout << std::endl;
    }
    for (auto c: fp_bt->chars_) {
        std::cout << c << ":";
        for (int i = 0; i < fp_bt->c_ranks_[fp_bt->chars_index_[c]].size(); i++) {
            std::cout << (int) fp_bt->c_ranks_[fp_bt->chars_index_[c]][i].width() << "/" << (int) fp_bt->c_ranks_[fp_bt->chars_index_[c]][i].size() << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Errors space " << lpf_bt->print_space_usage() << std::endl;
//    for (int i = 0; i < fp_bt->block_tree_types_.size(); i++) {
//        std::cout << i << " lvl/bv_s/bv_rs/pointer(n,w,s)/pointer(n,w,s) ";
//        std::cout << fp_bt->block_tree_types_[i]->size() << "/";
//        std::cout << fp_bt->block_tree_types_rs_[i]->space_usage() << "/";
//        std::cout << "(" << fp_bt->block_tree_pointers_[i]->size() << "," << (int) fp_bt->block_tree_pointers_[i]->width() << "," << fp_bt->block_tree_pointers_[i]->bit_size() << ")/";
//        std::cout << "(" << fp_bt->block_tree_offsets_[i]->size() << "," << (int) fp_bt->block_tree_offsets_[i]->width() << "," << fp_bt->block_tree_offsets_[i]->bit_size() << ")" << std::endl;
//    }

    for (int i = 0; i < lpf_bt->block_tree_types_.size(); i++) {
        std::cout << i << " lvl/bv_s/bv_rs/pointer(n,w,s)/pointer(n,w,s) ";
        std::cout << lpf_bt->block_tree_types_[i]->size() << "/";
        std::cout << lpf_bt->block_tree_types_rs_[i]->space_usage() << "/";
        std::cout << "(" << lpf_bt->block_tree_pointers_[i]->size() << "," << (int) lpf_bt->block_tree_pointers_[i]->width() << "," << lpf_bt->block_tree_pointers_[i]->bit_size() << ")/";
        std::cout << "(" << lpf_bt->block_tree_offsets_[i]->size() << "," << (int) lpf_bt->block_tree_offsets_[i]->width() << "," << lpf_bt->block_tree_offsets_[i]->bit_size() << ")" << std::endl;
    }
    for (int i = 0; i < fp_bt->block_tree_types_.size(); i++) {
        std::cout << i << " lvl/bv_s/bv_rs/pointer(n,w,s)/pointer(n,w,s) ";
        std::cout << fp_bt->block_tree_types_[i]->size() << "/";
        std::cout << fp_bt->block_tree_types_rs_[i]->space_usage() << "/";
        std::cout << "(" << fp_bt->block_tree_pointers_[i]->size() << "," << (int) fp_bt->block_tree_pointers_[i]->width() << "," << fp_bt->block_tree_pointers_[i]->bit_size() << ")/";
        std::cout << "(" << fp_bt->block_tree_offsets_[i]->size() << "," << (int) fp_bt->block_tree_offsets_[i]->width() << "," << fp_bt->block_tree_offsets_[i]->bit_size() << ")" << std::endl;
    }
    for (int i = 0; i < lpf_bt->block_tree_types_.size(); i++) {
        std::vector<int> state = std::vector<int>(lpf_bt->block_tree_types_[i]->size(), 2);
        for (int k = 0; k < lpf_bt->block_tree_types_[i]->size(); k++) {
            if ((*lpf_bt->block_tree_types_[i])[k] == 0) {
                state[k] = 0;
            }
        }
        auto result = std::vector<int>(4,0);
        for (auto s: state) {
            result[s]++;
        }
        std::cout << "result " << i << ":";
        for (auto r: result) {
            std::cout << r << " ";
        }
        std::cout << std::endl;
        for (int k = 0; k < lpf_bt->block_tree_pointers_[i]->size(); k++) {
            int ptr = (*lpf_bt->block_tree_pointers_[i])[k];
            int off = (*lpf_bt->block_tree_offsets_[i])[k];
            if (state[ptr] == 2) {
                state[ptr] = 1;
            }
            if (state[ptr] == 0) {
                state[ptr] = 3;
            }
            if (off > 0) {
                if (state[ptr + 1] == 2) {
                    state[ptr + 1] = 1;
                }
                if (state[ptr + 1] == 0) {
                    state[ptr + 1] = 3;
                }
            }


        }
        result = std::vector<int>(4,0);
        for (auto s: state) {
            result[s]++;
        }
        std::cout << "result " << i << ":";
        for (auto r: result) {
            std::cout << r << " ";
        }
        std::cout << std::endl;
    }
//    std::cout << *lpf_bt->block_tree_types_[7] << std::endl;
//    std::cout << *lpf_bt->block_tree_types_[8] << std::endl;
//    std::cout << *lpf_bt->block_tree_types_[9] << std::endl;
//    std::cout << *lpf_bt->block_tree_types_[10] << std::endl;
//    for (int i = 0; i < fp_bt->block_tree_types_.size(); i++) {
//        for (int j = 0; j < lpf_bt->block_tree_pointers_[i]->size() && j < fp_bt->block_tree_pointers_[i]->size(); j++)
//        {
//            if ((*fp_bt->block_tree_pointers_[i])[j] != (*lpf_bt->block_tree_pointers_[i])[j]) {
//                std::cout << i << " " << j << " " << (*fp_bt->block_tree_pointers_[i])[j] << " " << (*lpf_bt->block_tree_pointers_[i])[j] << std::endl;
//            }
//        }
//    }
//    for (int i = 0; i < fp_bt->block_tree_types_.size(); i++) {
//        for (int j = 0; j < lpf_bt->block_tree_offsets_[i]->size() && j < fp_bt->block_tree_offsets_[i]->size(); j++)
//        {
//            if ((*fp_bt->block_tree_offsets_[i])[j] != (*lpf_bt->block_tree_offsets_[i])[j]) {
//                std::cout << i << " " << j << " " << (*fp_bt->block_tree_offsets_[i])[j] << " " << (*lpf_bt->block_tree_offsets_[i])[j] << std::endl;
//            }
//        }
//    }
//    for (auto bv: lpf_bt->block_tree_types_) {
//        std::cout << *bv << std::endl;
//    }
//    for (auto bv: fp_bt->block_tree_types_) {
//        std::cout << *bv << std::endl;
//    }

    j = 0;

//    for (auto c: lpf_bt->chars_) {
//        int count = 0;
//        for (int i = 0; i < vec.size(); i++) {
//            if (vec[i] == c) {
//                count++;
//            }
//            auto x = lpf_bt->rank(c, i);
//            if (x != count) {
//                std::cout << c << ":" << i << " " << x << " " << count << std::endl;
//                j++;
//            }
//
//
//        }
//        if (c=='e') std::cout << count <<" e's" << std::endl;
//        std::cout << c << " Errors " << j << std::endl;
//    }
    std::cout << "Errors " << j << std::endl;
//    for (int i = 0; i < (*lpf_bt->block_tree_pointers_[1]).size(); i++) {
//        if ((*lpf_bt->block_tree_pointers_[1])[i] % lpf_bt->tau_ == lpf_bt->tau_ - 1) {
//            std::cout << "warum " << i << (*lpf_bt->block_tree_pointers_[1])[i] <<  std::endl;
//        }
//    }
//    for (auto c: lpf_bt->chars_) {
//        auto d = lpf_bt->chars_index_[c];
//        for (int i = 0; i < lpf_bt->c_ranks_[d][1].size(); i++) {
//            if (i % 2 == 1) {
//                if (lpf_bt->c_ranks_[d][1][i - 1] > lpf_bt->c_ranks_[d][1][i]) {
//                    std::cout << "help" << std::endl;
//                }
//            }
//        }
//    }
//    std::cout << "Pointers " << std::endl;
//    for (int i = 0; i < lpf_bt->block_tree_pointers_.size(); i++) {
//        for (int j = 0; j < lpf_bt->block_tree_pointers_[i]->size(); j++) {
//            std::cout << (*lpf_bt->block_tree_pointers_[i])[j] << " ";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << "Offsets " << std::endl;
//    for (int i = 0; i < lpf_bt->block_tree_offsets_.size(); i++) {
//        for (int j = 0; j < lpf_bt->block_tree_offsets_[i]->size(); j++) {
//            std::cout << (*lpf_bt->block_tree_offsets_[i])[j] << " ";
//        }
//        std::cout << std::endl;
//    }
//    for (int i = 0; i < vec.size(); i++) {
//        if (fp_bt->access(i) != vec[i]) {
//            j++;
////            std::cout << i << std::endl;
//        }
////        std::cout << fp_bt->access(i) << " " <<  vec[i] << std::endl;
//    }
//    std::vector<int> error = {206,207,268,269,772,773,802,803,830,831,1014,1015,1022,1023};
//    for (auto q: error) {
//        auto x = fp_bt->access(q);
//        if (x != vec[q]) {
//            std::cout << q <<" " << x <<  " " << vec[q] <<  std::endl;
//        }
//    }
//    for (int i = 0; i < lpf_bt->block_tree_types_.size(); i++) {
//        for (int j = 0; j < (*lpf_bt->block_tree_types_[i]).size(); j++) {
//            if ((bool)(*lpf_bt->block_tree_types_[i])[j] != (bool)(*fp_bt->block_tree_types_[i])[j]) {
//                std::cout << i << " " << j << " lpf " << (*lpf_bt->block_tree_types_[i])[j] << " fp " << (*fp_bt->block_tree_types_[i])[j] << std::endl;
//            }
//        }
//    }



//    for (auto bv: lpf_bt->block_tree_types_) {
//        std::cout << *bv << std::endl;
//    }
//    std::vector<int32_t> lpf(test.size());
//    std::vector<int32_t> lpf_ptr(test.size());
//    int64_t lz = 0;
//    int64_t k = 0;
//    while (lz < test.size()) {
//        lz = lz + std::max(1LL, lpf[lz + 1]);
//        k++;
//    }
//    std::cout << "lz " << lz.size() << std::endl;
//    delete lpf_bt;
//    delete fp_bt;
//
//    for (int i = 0; i < lpf_bt->block_size_lvl_.size(); i++) {
//        std::cout << i << ":" << lpf_bt->block_size_lvl_[i] << std::endl;
//    }
//    auto d = lpf_bt->chars_index_['e'];
//    std::cout << "char_ranks "  << 'e' << std::endl;
//    for (int i = 0;i < lpf_bt->c_ranks_[d].size(); i++) {
//        for (int j = 0;j < lpf_bt->c_ranks_[d][i].size(); j++) {
//            std::cout << lpf_bt->c_ranks_[d][i][j] << " ";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << "pointer_ranks "  << 'e' << std::endl;
//    for (int i = 0;i < lpf_bt->pointer_c_ranks_[d].size(); i++) {
//        for (int j = 0;j < lpf_bt->pointer_c_ranks_[d][i].size(); j++) {
//            std::cout << lpf_bt->pointer_c_ranks_[d][i][j] << " ";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << test << std::endl;
//    for (auto c: lpf_bt->leaves_) std::cout << c;
//    std::cout << std::endl;
//    std::cout << test[1743] << std::endl;
//    std::cout << test[1744] << std::endl;
    std::unordered_map<char, int64_t> hist;
    std::unordered_set<int> characters;
    std::random_device rnd_device;
    std::mt19937 mersenne_engine(rnd_device());
    std::vector<int> access_queries_;
    std::vector<int> select_queries_;
    std::vector<uint8_t> select_c_;
    std::uniform_int_distribution<uint64_t> dist(0, test.size() - 1);
    for (int i =  0; i < test.size(); i++) {
        hist[test[i]] = hist[test[i]] + 1;
    }
    for (size_t i = 0; i < 1000000; ++i) {
        access_queries_.push_back(dist(mersenne_engine));
    }
    for (size_t i = 0; i < 1000000; ++i) {
        uint8_t x = 0;
        int xsum = 0;
        while (xsum + hist[x] < access_queries_[i]) {
            xsum += hist[x];
            x++;
        }
        select_c_.push_back(x);
        select_queries_.push_back(access_queries_[i] - xsum);
    }
    std::cout << "Starting Queries" << "\n";
    size_t result = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (auto const& query : access_queries_) {
        result += lpf_bt->access(query);
    }
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "Starting Rank Queries" << "\n";
    auto start2 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < select_c_.size(); i++) {
        result += lpf_bt->rank(select_c_[i], access_queries_[i]);
    }
    auto elapsed2 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start2).count();
    double nanosec_per_access = elapsed / static_cast<double>(access_queries_.size());
    double nanosec_per_rank = elapsed2 / static_cast<double>(access_queries_.size());
    std::cout << "#access queries " << access_queries_.size() << std::endl;
    std::cout << nanosec_per_access << " ns per query\n";
    std::cout << "#rank queries " << access_queries_.size() << std::endl;
    std::cout << nanosec_per_rank << " ns per query\n";
    delete lpf_bt;
return 0;
}
