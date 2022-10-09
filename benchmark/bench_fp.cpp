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
#include <malloc_count.h>
#include <type_traits>
#include <iostream>
#include "malloc_count.h"
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
//    std::ifstream t("/home/daniel/blocktree-experiments/data/english.1024MB");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/einstein.de.txt");
    std::ifstream t("/home/daniel/blocktree-experiments/data/einstein.en.txt");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/influenza");
//    std::ifstream t("/Users/daniel/Downloads/einstein.en.txt");
    std::stringstream buffer;
    t.read(&test[0], a_size);
    test = "NNBOBOTWNNBOBIOOTBSHTFNEBOBOTWNEBOBOTWNEBOBIOOTBSHTFNSBOBOTW";
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
    BV_BlockTree_lpf_heuristic<uint8_t, int32_t>*  lpf_bt = new BV_BlockTree_lpf_heuristic<uint8_t, int32_t>(vec, 2, 1, 15, lpf, lpf_ptr, lz);
    std::cout << "lpf" << std::endl;
    lpf_bt->add_rank_support();
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
    int  j = 0;
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
    std::cout << j << " Errors " << std::endl;
    for (int i = 0; i < test.size(); i++) {
        std::cout << i << " " << lpf_bt->rank('N', i) << std::endl;
    }

    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_['N']][0].size(); i++) {
        std::cout << lpf_bt->c_ranks_[lpf_bt->chars_index_['N']][0][i] << std::endl;
    }
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_['N']][1].size(); i++) {
        std::cout << lpf_bt->c_ranks_[lpf_bt->chars_index_['N']][1][i] << std::endl;
    }
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
    delete lpf_bt;
return 0;
}
