#include <iostream>
#include <libsais.h>
#include <lpfarray.hpp>
#include <vector>
#include <bit_vector.hpp>
#include <tlx/cmdline_parser.hpp>
#include <fstream>
#include <sstream>
#include <bv_blocktree_lpf_64.hpp>
#include <chrono>
#include <type_traits>
#include <iostream>

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
//    std::ifstream t("/home/daniel/blocktree-experiments/data/influenza");
    std::ifstream t("/Users/daniel/Downloads/einstein.en.txt");
    std::stringstream buffer;
    t.read(&test[0], a_size);
//    test = "NNBOBOTWNNBOBIOOTBSHTFNEBOBOTWNEBOBOTWNEBOBIOOTBSHTFNSBOBOTW";
//   test = "ABCDABCDEFGHABCDDE12";
    std::vector<uint8_t> vec(test.begin(), test.end());
    std::vector<int64_t> lpf(test.size());
    std::vector<int64_t> lpf_ptr(test.size());
//    int64_t lz = 0;
//    int64_t k = 0;
//    while (lz < test.size()) {
//        lz = lz + std::max(1LL, lpf[lz + 1]);
//        k++;
//    }
//    std::cout << "lz " << lz.size() << std::endl;
    auto t01 = std::chrono::high_resolution_clock::now();
    BV_BlockTree_lpf_64 bt(vec, 2, 1, 15);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << "time " <<  ms_int.count() << std::endl;
    std::vector<uint8_t> access_vec(0);
    int j = 0;
    for (int i = 0; i < test.size(); i++) {
        auto x= bt.access(i);
        access_vec.push_back(x);
        if (static_cast<int>(x) != static_cast<int>(vec[i])) {
            std::cout <<i << " " <<  static_cast<int>(x) << " " << static_cast<int>(vec[i]) << std::endl;
            j++;
        }
    }
        std::cout << j << std::endl;
//    std::cout << bt.access(6910) << " " << test[6910] << std::endl;
//    std::cout << bt.access(6911) << " " << test[6911] << std::endl;
////    std::cout << bt.access(412) << " " << test[412] << std::endl;
////    std::cout << bt.access(784) << " " << test[784] << std::endl;
////    std::cout << bt.access(785) << " " << test[785] << std::endl;
    return 0;
}
