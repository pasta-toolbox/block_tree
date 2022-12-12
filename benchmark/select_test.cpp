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
    std::ifstream t("/home/daniel/blocktree-experiments/data/pitches");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/world_leaders");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/einstein.de.txt");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/pitches");
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

    std::vector<int32_t> lpf(vec.size());
    std::vector<int32_t> lpf_ptr(vec.size());
    lpf_array_stack(vec,lpf,lpf_ptr);
    std::vector<int32_t> lz;

    int32_t lzn = 0;
    calculate_lz_factor(lzn,lpf,lz);
    BV_BlockTree_fp_theory<uint8_t, int32_t>*  lpf_bt = new BV_BlockTree_fp_theory<uint8_t, int32_t>(vec, 16, 16,lzn, 256,
                                                                                                       true);
    std::cout << "time " << lpf_bt->block_tree_types_.size() << std::endl;
    lpf_bt->add_rank_support();
    std::cout << "time " << lpf_bt->block_tree_types_.size() << std::endl;

        int j = 0;
        for (auto c: lpf_bt->chars_) {
        int count = 0;

        for (int i = 0; i < vec.size(); i++) {
            if (vec[i] == c) {

                count++;
                auto x = lpf_bt->select(c, count);
                if (x != i) {
                    std::cout << (int)c << ":" << i << " " << x << " " << count << std::endl;
                    j++;
                }
            }



        }
        std::cout << c << "  select Errors " << j << std::endl;
    }
    out:
    delete lpf_bt;
    return 0;
}