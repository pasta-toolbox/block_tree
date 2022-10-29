#include <iostream>
#include <lpfarray.hpp>
#include <vector>
#include "bit_vector.hpp"
#include "cmdline_parser.hpp"
#include <lpfconstruction/bv_blocktree_lpf_heuristic.hpp>
#include <lpfconstruction/bv_blocktree_lpf_pruned.hpp>
#include <lpfconstruction/bv_blocktree_lpf_theory.hpp>
#include <fppconstruction/bv_blocktree_fp_theory.hpp>
#include <fppconstruction/bv_blocktree_fp_pruned.hpp>
#include <chrono>
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
    BV_BlockTree_lpf_pruned<uint8_t, int32_t>*  lpf_bt = new BV_BlockTree_lpf_pruned<uint8_t, int32_t>(vec, 2, 1, 15, lpf, lpf_ptr, lz, true);
    lpf_bt->add_rank_support();

    std::cout << "select " << lpf_bt->select('N', 0) << std::endl;
    std::cout << "select " << lpf_bt->select('N', 1) << std::endl;
    std::cout << "select " << lpf_bt->select('N', 2) << std::endl;
    std::cout << "select " << lpf_bt->select('N', 3) << std::endl;
    std::cout << "select " << lpf_bt->select('N', 4) << std::endl;
    std::cout << "select " << lpf_bt->select('N', 5) << std::endl;
    std::cout << "select " << lpf_bt->select('N', 6) << std::endl;
    std::cout << "select " << lpf_bt->select('N', 7) << std::endl;
    std::cout << "select " << lpf_bt->select('N', 8) << std::endl;
    std::cout << "select " << lpf_bt->select('N', 9) << std::endl;
    std::cout << "select " << lpf_bt->select('N', 10) << std::endl;
    std::cout << "select " << lpf_bt->select('N', 11) << std::endl;

    delete lpf_bt;
    return 0;
}