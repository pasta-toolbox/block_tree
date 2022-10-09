#include <tlx/die.hpp>
#include <vector>
#include <lpfconstruction/bv_blocktree_lpf_heuristic.hpp>
#include <lpfarray.hpp>
#include <chrono>
#include "gtest/gtest.h"

TEST(BlockTreeTest, paper_example_test_access) {
    std::string test = "NNBOBOTWNNBOBIOOTBSHTFNEBOBOTWNEBOBOTWNEBOBIOOTBSHTFNSBOBOTW";
    std::vector<uint8_t> vec(test.begin(), test.end());
    std::vector<int32_t> lpf(vec.size());
    std::vector<int32_t> lpf_ptr(vec.size());
    std::vector<int32_t> lz;
    int32_t lzn = 0;
    lpf_array(vec, lpf, lpf_ptr);
    calculate_lz_factor(lzn,lpf, lz);
    BV_BlockTree_lpf_heuristic<uint8_t, int32_t>*  lpf_bt = new BV_BlockTree_lpf_heuristic<uint8_t, int32_t>(vec, 2, 1, 15, lpf, lpf_ptr, lz);

    // test access
    for (int i = 0; i < vec.size(); i++) {
        ASSERT_EQ(lpf_bt->access(i), vec[i]);
    }
}

TEST(BlockTreeTest, paper_example_test_rank) {
    std::string test = "NNBOBOTWNNBOBIOOTBSHTFNEBOBOTWNEBOBOTWNEBOBIOOTBSHTFNSBOBOTW";
    std::vector<uint8_t> vec(test.begin(), test.end());
    std::vector<int32_t> lpf(vec.size());
    std::vector<int32_t> lpf_ptr(vec.size());
    std::vector<int32_t> lz;
    int32_t lzn = 0;
    lpf_array(vec, lpf, lpf_ptr);
    calculate_lz_factor(lzn,lpf, lz);
    BV_BlockTree_lpf_heuristic<uint8_t, int32_t>*  lpf_bt = new BV_BlockTree_lpf_heuristic<uint8_t, int32_t>(vec, 2, 1, 15, lpf, lpf_ptr, lz);
    lpf_bt->add_rank_support();
    std::vector<uint8_t> chars = {'N', 'B', 'O', 'T', 'W' , 'I', 'S', 'H', 'F', 'E'};
    // test rank
    for (auto c: chars) {
        int count = 0;
        for (int i = 0; i < vec.size(); i++) {
            if (vec[i] == c) {
                count++;
            }
            ASSERT_EQ(lpf_bt->rank(c, i), count) << " rank and count differ at index " << i;
        }
    }
}
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}