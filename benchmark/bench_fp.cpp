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
//    std::ifstream t("/home/daniel/blocktree-experiments/data/world_leaders");

//    std::ifstream t("/home/daniel/blocktree-experiments/data/influenza");
//    std::ifstream t("/Users/daniel/Downloads/einstein.en.txt");
    std::stringstream buffer;
    t.read(&test[0], a_size);
//    test = "NNBOBOTWNNBOBIOOTBSHTFNEBOBOTWNEBOBOTWNEBOBIOOTBSHTFNSBOBOTW";
//    test = "abbaabbbaaabab";
//test = "aabaaaaaaa";
//   test = "ABCDABCDEFGHABCDDE12";
    std::vector<uint8_t> vec(test.begin(), test.end());


//    BV_BlockTree_fp_theory<uint8_t, int32_t> *fp_bt = new BV_BlockTree_fp_theory<uint8_t, int32_t>(vec,2,1,1);

//    BV_BlockTree_fp_pruned<uint8_t, int32_t> *fp_bt = new BV_BlockTree_fp_pruned<uint8_t, int32_t>(vec, 16, 4, 1);



//    std::vector<int64_t> lpf2(vec.size());
//    std::vector<int64_t> lpf_ptr2(vec.size());
//    std::vector<int64_t> lz2;

    int64_t numberOfLzPhrases = 0;
    auto t0a = std::chrono::high_resolution_clock::now();
//    lpf_array_ansv(vec, lpfArray, prevOcc);

//    lpf_array_stack(vec, lpf2, lpf_ptr2);
    auto t0b = std::chrono::high_resolution_clock::now();
    auto ms_intstack = std::chrono::duration_cast<std::chrono::milliseconds>(t0b - t0a);
    std::cout << "LPF TIME stack " << ms_intstack.count() << std::endl;


    auto t01 = std::chrono::high_resolution_clock::now();
//    lpf_array_ansv(vec, lpfArray, prevOcc);
    auto t02 = std::chrono::high_resolution_clock::now();

//    for (int i = 0; i < lpfArray.size(); i++) {
//        std::cout << lpfArray[i] << " ";
//    }
//    std::cout << std::endl;
//
//
//    for (int i = 0; i < prevOcc.size(); i++) {
//        std::cout << prevOcc[i] << " ";
//    }

    std::cout << std::endl;
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << "LPF TIME " << ms_int.count() << std::endl;
    bool DONT_CUT_FIRST_LEVELS = false;
    bool CUT_FIRST_LEVELS = true;
    bool EXTENDED_PRUNE = true;
    bool SIMPLE_PRUNE = false;
    int tau = 16;
    int mls = 16;
    int64_t s = 1;
    std::cout << "#LZ PHRASES " << s << std::endl;
//    BV_BlockTree_lpf_heuristic<uint8_t, int32_t>* fp_bt = new BV_BlockTree_lpf_heuristic<uint8_t, int32_t>(vec, 2, 1);
//    BV_BlockTree_lpf_theory<uint8_t, int32_t>* fp_bt = new BV_BlockTree_lpf_theory<uint8_t, int32_t>(vec, 2, 1,1,lpfArray, prevOcc, lzPhrases);
    auto t03 = std::chrono::high_resolution_clock::now();
    auto* fpTheory_bt = new BV_BlockTree_fp_theory<uint8_t, int32_t>(vec, tau, mls,s,256, CUT_FIRST_LEVELS);
    auto t04 = std::chrono::high_resolution_clock::now();
    auto ms_int_fp_theo = std::chrono::duration_cast<std::chrono::milliseconds>(t04 - t03);
    std::cout << "FP THEORY TIME " << ms_int_fp_theo.count() << std::endl;
    auto* fpPruned_bt = new BV_BlockTree_fp_pruned<uint8_t, int32_t>(vec, tau, mls,s,256,CUT_FIRST_LEVELS,
                                                                     EXTENDED_PRUNE);
    auto t05 = std::chrono::high_resolution_clock::now();
    auto ms_int_fp_prune = std::chrono::duration_cast<std::chrono::milliseconds>(t05 - t04);
    std::cout << "FP PRUNE TIME " << ms_int_fp_prune.count() << std::endl;
    auto t0x = std::chrono::high_resolution_clock::now();
    auto* fpPruned_simple_bt = new BV_BlockTree_fp_pruned<uint8_t, int32_t>(vec, tau, mls,s,256, CUT_FIRST_LEVELS,
                                                                            SIMPLE_PRUNE);
    auto t0y = std::chrono::high_resolution_clock::now();
    auto ms_int_fp_prune_simple = std::chrono::duration_cast<std::chrono::milliseconds>(t0y - t0x);
    std::cout << "FP PRUNE SIMPLE TIME " << ms_int_fp_prune_simple.count() << std::endl;
    auto* lpfHeuristic_bt = new BV_BlockTree_lpf_heuristic<uint8_t, int32_t>(vec, tau, mls);
    auto t06 = std::chrono::high_resolution_clock::now();
    auto ms_int_lpf_heu = std::chrono::duration_cast<std::chrono::milliseconds>(t06 - t0y);
    std::cout << "LPF Heuristic TIME " << ms_int_lpf_heu.count() << std::endl;
    auto* lpfTheory_bt = new BV_BlockTree_lpf_theory<uint8_t, int32_t>(vec, tau, mls, false);
    auto t07 = std::chrono::high_resolution_clock::now();
    auto ms_int_lpf_theo = std::chrono::duration_cast<std::chrono::milliseconds>(t07 - t06);
    std::cout << "LPF THEORY TIME " << ms_int_lpf_theo.count()<< std::endl;
    auto* lpfTheory_bt_dp = new BV_BlockTree_lpf_theory<uint8_t, int32_t>(vec, tau, mls, true);
    auto t0aa = std::chrono::high_resolution_clock::now();
    auto ms_int_lpf_theo_dp = std::chrono::duration_cast<std::chrono::milliseconds>(t0aa - t07);
    std::cout << "LPF THEORY DP TIME " << ms_int_lpf_theo_dp.count()<< std::endl;
    auto* lpfPruned_bt = new BV_BlockTree_lpf_pruned<uint8_t, int32_t>(vec, tau, mls,
                                                                       true, CUT_FIRST_LEVELS, true);
    auto t08 = std::chrono::high_resolution_clock::now();
    auto ms_int_lpf_prune = std::chrono::duration_cast<std::chrono::milliseconds>(t08 - t0aa);
    std::cout << "LPF PRUNED TIME " << ms_int_lpf_prune.count() << std::endl;
    auto* lpfPruned_bt_dp = new BV_BlockTree_lpf_pruned<uint8_t, int32_t>(vec, tau, mls,
                                                                       true, CUT_FIRST_LEVELS, false);
    auto t09 = std::chrono::high_resolution_clock::now();
    auto ms_int_lpf_prune_dp = std::chrono::duration_cast<std::chrono::milliseconds>(t09 - t08);
    std::cout << "LPF PRUNED DP TIME " << ms_int_lpf_prune_dp.count() << std::endl;
    std::cout << "time " << ms_int.count() << std::endl;
    std::cout << "fp theory space " << fpTheory_bt->print_space_usage() << " " << fpTheory_bt->block_tree_types_.size() <<  std::endl;
    std::cout << "fp pruned space " << fpPruned_bt->print_space_usage() << " " << fpPruned_bt->block_tree_types_.size() <<   std::endl;
    std::cout << "fp pruned simple space " << fpPruned_simple_bt->print_space_usage() << " " << fpPruned_bt->block_tree_types_.size() <<   std::endl;
    std::cout << "lpfArray heuristic space " << lpfHeuristic_bt->print_space_usage() << " " << lpfHeuristic_bt->block_tree_types_.size() <<   std::endl;
    std::cout << "lpfArray theory space " << lpfTheory_bt->print_space_usage() << " " << lpfTheory_bt->block_tree_types_.size() <<   std::endl;
    std::cout << "lpfArray theory dp space " << lpfTheory_bt_dp->print_space_usage() << " " << lpfTheory_bt_dp->block_tree_types_.size() <<   std::endl;
    std::cout << "lpfArray pruned space " << lpfPruned_bt->print_space_usage() << " " << lpfPruned_bt->block_tree_types_.size() <<  std::endl;
    std::cout << "lpfArray pruned dp space " << lpfPruned_bt_dp->print_space_usage() << " " << lpfPruned_bt_dp->block_tree_types_.size() <<  std::endl;
//    BV_BlockTree_lpf_pruned<uint8_t, int32_t>*  lpf_bt2 = new BV_BlockTree_lpf_pruned<uint8_t, int32_t>(vec, 2, 1, 15,lpfArray, prevOcc, lzPhrases, false);
    std::cout << "lpfArray" << std::endl;
//    BV_BlockTree_lpf_pruned<uint8_t, int32_t>*  lpf_bt = new BV_BlockTree_lpf_pruned<uint8_t, int32_t>(vec, 2, 1, 15,lpfArray, prevOcc, lzPhrases, true);

//    for (int i = 0; i < lpf_bt->block_tree_types_.size(); i++) {
//        std::cout << i << " lvl/bv_s/bv_rs/pointer(n,w,s)/pointer(n,w,s) ";
//        std::cout << lpf_bt->block_tree_types_[i]->size() << "/";
//        std::cout << lpf_bt->block_tree_types_rs_[i]->space_usage() << "/";
//        std::cout << "(" << lpf_bt->block_tree_pointers_[i]->size() << "," << (int) lpf_bt->block_tree_pointers_[i]->width() << "," << lpf_bt->block_tree_pointers_[i]->bit_size() << ")/";
//        std::cout << "(" << lpf_bt->block_tree_offsets_[i]->size() << "," << (int) lpf_bt->block_tree_offsets_[i]->width() << "," << lpf_bt->block_tree_offsets_[i]->bit_size() << ")" << std::endl;
//    }
//    for (int i = 0; i < lpf_bt2->block_tree_types_.size(); i++) {
//        std::cout << i << " lvl/bv_s/bv_rs/pointer(n,w,s)/pointer(n,w,s) ";
//        std::cout << lpf_bt2->block_tr ee_types_[i]->size() << "/";
//        std::cout << lpf_bt2->block_tree_types_rs_[i]->space_usage() << "/";
//        std::cout << "(" << lpf_bt2->block_tree_pointers_[i]->size() << "," << (int) lpf_bt2->block_tree_pointers_[i]->width() << "," << lpf_bt2->block_tree_pointers_[i]->bit_size() << ")/";
//        std::cout << "(" << lpf_bt2->block_tree_offsets_[i]->size() << "," << (int) lpf_bt2->block_tree_offsets_[i]->width() << "," << lpf_bt2->block_tree_offsets_[i]->bit_size() << ")" << std::endl;
//    }
//
//    for (int i = 0; i < fp_bt->block_tree_types_.size(); i++) {
//        std::cout << i << " lvl/bv_s/bv_rs/pointer(n,w,s)/pointer(n,w,s) ";
//        std::cout << fp_bt->block_tree_types_[i]->size() << "/";
//        std::cout << fp_bt->block_tree_types_rs_[i]->space_usage() << "/";
//        std::cout << "(" << fp_bt->block_tree_pointers_[i]->size() << "," << (int) fp_bt->block_tree_pointers_[i]->width() << "," << fp_bt->block_tree_pointers_[i]->bit_size() << ")/";
//        std::cout << "(" << fp_bt->block_tree_offsets_[i]->size() << "," << (int) fp_bt->block_tree_offsets_[i]->width() << "," << fp_bt->block_tree_offsets_[i]->bit_size() << ")" << std::endl;
//    }
    int  j = 0;
    bool error_mode = false;
    for (int i = 0; i < vec.size(); i++) {
        auto x = lpfPruned_bt->access(i);
        if (x != vec[i]) {
            if (!error_mode) {
//            std::cout << i << std::endl;
                error_mode = true;
            }
//            std::cout << i << std::endl;
            j++;
        }

    }
    std::cout << "Errors " << j << std::endl;


    for (int i = 0; i < vec.size(); i++) {
        auto x = lpfTheory_bt->access(i);
        if (x != vec[i]) {
            j++;
        }

    }
    int32_t threads = 6;
    std::cout << "Errors " << j << std::endl;
    for (int i = 0; i < vec.size(); i++) {
        auto x = lpfPruned_bt_dp->access(i);
        if (x != vec[i]) {
            j++;
        }

    }
    std::cout << "Errors " << j << std::endl;


    for (int i = 0; i < vec.size(); i++) {
        auto x = lpfHeuristic_bt->access(i);
        if (x != vec[i]) {
            j++;
        }

    }

    std::cout << "Errors " << j << std::endl;
    auto t0r1 = std::chrono::high_resolution_clock::now();
    fpPruned_bt->add_rank_support_omp(threads);
    auto t0r2 = std::chrono::high_resolution_clock::now();
    auto ms_int_rs = std::chrono::duration_cast<std::chrono::milliseconds>(t0r2 - t0r1);
    std::cout << "RS TIME " << threads << ": " <<  ms_int_rs.count() << std::endl;
    auto t0r3 = std::chrono::high_resolution_clock::now();
    lpfPruned_bt->add_rank_support();
    auto t0r4 = std::chrono::high_resolution_clock::now();
    auto ms_int_rs2 = std::chrono::duration_cast<std::chrono::milliseconds>(t0r4  - t0r3);
    std::cout << "RS TIME " << 1 << ": " <<  ms_int_rs2.count() << std::endl;

//    for (auto c: lpf_bt->chars_) {
//        std::cout << c << ":";
//        for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_[c]].size(); i++) {
//            std::cout << (int) lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][i].width() << "/" << (int) lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][i].size() << " ";
//        }
//        std::cout << std::endl;
//    }
    std::cout << "Errors space " << fpPruned_bt->print_space_usage() << std::endl;
//    for (int i = 0; i < fp_bt->block_tree_types_.size(); i++) {
//        std::cout << i << " lvl/bv_s/bv_rs/pointer(n,w,s)/pointer(n,w,s) ";
//        std::cout << fp_bt->block_tree_types_[i]->size() << "/";
//        std::cout << fp_bt->block_tree_types_rs_[i]->space_usage() << "/";
//        std::cout << "(" << fp_bt->block_tree_pointers_[i]->size() << "," << (int) fp_bt->block_tree_pointers_[i]->width() << "," << fp_bt->block_tree_pointers_[i]->bit_size() << ")/";
//        std::cout << "(" << fp_bt->block_tree_offsets_[i]->size() << "," << (int) fp_bt->block_tree_offsets_[i]->width() << "," << fp_bt->block_tree_offsets_[i]->bit_size() << ")" << std::endl;
//    }

//    for (int i = 0; i < lpf_bt->block_tree_types_.size(); i++) {
//        std::vector<int> state = std::vector<int>(lpf_bt->block_tree_types_[i]->size(), 2);
//        for (int k = 0; k < lpf_bt->block_tree_types_[i]->size(); k++) {
//            if ((*lpf_bt->block_tree_types_[i])[k] == 0) {
//                state[k] = 0;
//            }
//        }
//        auto result = std::vector<int>(4,0);
//        for (auto s: state) {
//            result[s]++;
//        }
//        std::cout << "result " << i << ":";
//        for (auto r: result) {
//            std::cout << r << " ";
//        }
//        std::cout << std::endl;
//        for (int k = 0; k < lpf_bt->block_tree_pointers_[i]->size(); k++) {
//            int ptr = (*lpf_bt->block_tree_pointers_[i])[k];
//            int off = (*lpf_bt->block_tree_offsets_[i])[k];
//            if (state[ptr] == 2) {
//                state[ptr] = 1;
//            }
//            if (state[ptr] == 0) {
//                state[ptr] = 3;
//            }
//            if (off > 0) {
//                if (state[ptr + 1] == 2) {
//                    state[ptr + 1] = 1;
//                }
//                if (state[ptr + 1] == 0) {
//                    state[ptr + 1] = 3;
//                }
//            }
//
//
//        }
//        result = std::vector<int>(4,0);
//        for (auto s: state) {
//            result[s]++;
//        }
//        std::cout << "result " << i << ":";
//        for (auto r: result) {
//            std::cout << r << " ";
//        }
//        std::cout << std::endl;
//    }
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
//
//    j = 0;
//
//    for (auto c: lpfPruned_bt->chars_) {
//        int count = 0;
//        for (int i = 0; i < vec.size(); i++) {
//            if (vec[i] == c) {
//                count++;
//            }
//            auto x = lpfPruned_bt->rank(c, i);
//            if (x != count) {
//                std::cout << c << ":" << i << " " << x << " " << count << std::endl;
//                j++;
//            }
//
//
//        }
//        if (c=='e') std::cout << count <<" e's" << std::endl;
//        std::cout << c << "  rank Errors " << j << std::endl;
//    }
//    std::cout << "Errors " << j << std::endl;
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
//                std::cout << i << " " << j << " lpfArray " << (*lpf_bt->block_tree_types_[i])[j] << " fp " << (*fp_bt->block_tree_types_[i])[j] << std::endl;
//            }
//        }
//    }



//    for (auto bv: lpf_bt->block_tree_types_) {
//        std::cout << *bv << std::endl;
//    }
//    std::vector<int32_t> lpfArray(test.size());
//    std::vector<int32_t> prevOcc(test.size());
//    int64_t lzPhrases = 0;
//    int64_t k = 0;
//    while (lzPhrases < test.size()) {
//        lzPhrases = lzPhrases + std::max(1LL, lpfArray[lzPhrases + 1]);
//        k++;
//    }
//    std::cout << "lzPhrases " << lzPhrases.size() << std::endl;
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
        result += lpfPruned_bt->access(query);
    }
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "Starting Encoded Queries" << "\n";
    result = 0;
    auto start3 = std::chrono::high_resolution_clock::now();
    for (auto const& query : access_queries_) {
        result += lpfPruned_bt->access(query);
    }
    auto elapsed3 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start3).count();
    std::cout << "Starting Rank Queries" << "\n";
    auto start2 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < select_c_.size(); i++) {
        result += lpfPruned_bt->rank(select_c_[i], access_queries_[i]);
    }
    auto elapsed2 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start2).count();
    auto start4 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < select_c_.size(); i++) {
            result += lpfPruned_bt->rank_base(select_c_[i], access_queries_[i]);
    }

    auto elapsed4 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start4).count();
    auto start5 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < select_c_.size(); i++) {
        result += lpfPruned_bt->select(select_c_[i], select_queries_[i]);
    }

    auto elapsed5 = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start5).count();
    for (int i = 0; i < select_c_.size(); i++) {
        if (lpfPruned_bt->rank_base(select_c_[i], access_queries_[i]) != lpfPruned_bt->rank(select_c_[i], access_queries_[i])) std::cout << i << std::endl;
    }
    double nanosec_per_access = elapsed / static_cast<double>(access_queries_.size());
    double nanosec_per_access2 = elapsed3 / static_cast<double>(access_queries_.size());
    double nanosec_per_rank = elapsed2 / static_cast<double>(access_queries_.size());
    double nanosec_per_rank2 = elapsed4 / static_cast<double>(access_queries_.size());
    double nanosec_per_select = elapsed5 / static_cast<double>(access_queries_.size());
    std::cout << "#access queries " << access_queries_.size() << std::endl;
    std::cout << nanosec_per_access << " ns per query\n";
    std::cout << "#access encoded queries " << access_queries_.size() << std::endl;
    std::cout << nanosec_per_access2 << " ns per query\n";
    std::cout << "#rank queries " << access_queries_.size() << std::endl;
    std::cout << nanosec_per_rank << " ns per query\n";
    std::cout << "#rank queries base " << access_queries_.size() << std::endl;
    std::cout << nanosec_per_rank2 << " ns per query\n";
    std::cout << "#select queries " << access_queries_.size() << std::endl;
    std::cout << nanosec_per_select << " ns per query\n";
//    for (int i = 0; i < lpf_bt->block_tree_types_.size(); i++) {
//        for (int j = 0; j < (*lpf_bt->block_tree_types_[i]).size(); j++) {
//            if ((bool)((*lpf_bt->block_tree_types_[i])[j]) != (bool)((*fp_bt->block_tree_types_[i])[j])) {
//                std::cout << i << " " << j << " " << (*lpf_bt->block_tree_types_[i])[j] << " " << (*fp_bt->block_tree_types_[i])[j] << std::endl;
//            }
//        }
//    }
    delete fpTheory_bt;
    delete fpPruned_bt;
    delete fpPruned_simple_bt;
    delete lpfPruned_bt;
    delete lpfPruned_bt_dp;
    delete lpfHeuristic_bt;
    delete lpfTheory_bt;
    delete lpfTheory_bt_dp;
return 0;
}
