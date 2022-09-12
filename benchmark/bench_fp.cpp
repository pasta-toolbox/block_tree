#include <iostream>
#include "libsais.h"
#include "lpfarray.hpp"
#include <vector>
#include "bit_vector.hpp"
#include "cmdline_parser.hpp"
#include <fstream>
#include <sstream>
#include "lpfconstruction/bv_blocktree_lpf_heuristic.hpp"
#include "lpfconstruction/bv_blocktree_lpf_pruned.hpp"
#include "lpfconstruction/bv_blocktree_lpf_theory.hpp"
#include "fppconstruction/bv_blocktree_fp_theory.hpp"
#include "fppconstruction/bv_blocktree_fp_pruned.hpp"
#include <chrono>
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
    std::cout << std::endl << "Run with " << a_size << " Bytes" << std::endl;
    std::string test(a_size, ' ');
//    std::ifstream t("/home/daniel/blocktree-experiments/data/Escherichia_Coli");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/english.1024MB");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/einstein.de.txt");
    std::ifstream t("/home/daniel/blocktree-experiments/data/einstein.en.txt");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/influenza");
//    std::ifstream t("/Users/daniel/Downloads/einstein.en.txt");
    std::stringstream buffer;
    t.read(&test[0], a_size);
//    std::string test = "NNBOBOTWNNBOBIOOTBSHTFNEBOBOTWNEBOBOTWNEBOBIOOTBSHTFNSBOBOTW";
//   test = "ABCDABCDEFGHABCDDE12";
    std::vector <uint8_t> vec(test.begin(), test.end());
    auto t01 = std::chrono::high_resolution_clock::now();
    std::cout << vec.size() << std::endl;
//    BV_BlockTree_fp_theory<uint8_t, int32_t> *bt = new BV_BlockTree_fp_theory<uint8_t, int32_t>(vec,2,1,1);

    BV_BlockTree_fp_pruned<uint8_t, int32_t> *bt2 = new BV_BlockTree_fp_pruned<uint8_t, int32_t>(vec,2,1,1);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << "time " << ms_int.count() << std::endl;
    int  j = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (bt2->access(i)  != vec[i]) {
            std::cout << i <<" " << bt2->access(i) <<  " " << vec[i] <<  std::endl;
            j++;
        }
    }
    std::cout << j << " Errors " << std::endl;
//    for (int i = 0; i < bt->block_tree_types_.size(); i++) {
//        std::cout << i << " lvl/bv_s/bv_rs/pointer(n,w,s)/pointer(n,w,s) ";
//        std::cout << bt->block_tree_types_[i]->size() << "/";
//        std::cout << bt->block_tree_types_rs_[i]->space_usage() << "/";
//        std::cout << "("<<bt->block_tree_pointers_[i]->size()<<","<<(int) bt->block_tree_pointers_[i]->width()<< ","<<bt->block_tree_pointers_[i]->bit_size()<<")/";
//        std::cout << "("<<bt->block_tree_offsets_[i]->size()<<","<< (int) bt->block_tree_offsets_[i]->width()<< ","<<bt->block_tree_offsets_[i]->bit_size()<<")"<<std::endl;
//    }
//    for (auto bv: bt->block_tree_types_) {
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
//    delete bt;
    delete bt2;
return 0;
}
