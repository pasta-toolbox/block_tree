#include <vector>
#include <lpfconstruction/bv_blocktree_lpf_heuristic.hpp>
#include <lpfarray.hpp>
#include <chrono>
#include "gtest/gtest.h"
class btTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::string str = "NNBOBOTWNNBOBIOOTBSHTFNEBOBOTWNEBOBOTWNEBOBIOOTBSHTFNSBOBOTW";
        vec = std::vector<uint8_t>(str.begin(), str.end());
        std::vector<int32_t> lpf(vec.size());
        std::vector<int32_t> lpf_ptr(vec.size());
        std::vector<int32_t> lz;
        int32_t lzn = 0;
        lpf_array(vec, lpf, lpf_ptr);
        calculate_lz_factor(lzn,lpf, lz);
        lpf_bt = new BV_BlockTree_lpf_heuristic<uint8_t, int32_t>(vec, 2, 1, 15, lpf, lpf_ptr, lz);
        lpf_bt->add_rank_support();
        ranks_first_lvl['N'] = {0,2,2,4,4,4,5,5,6,6,7,7,7,7,8};
        ranks_first_lvl['W'] = {0,0,1,1,1,1,1,1,2,2,3,3,3,3,3};
        ranks_first_lvl['O'] = {0,1,2,3,5,5,5,7,7,9,9,10,12,12,13};
        ranks_first_lvl['T'] = {0,0,1,1,1,2,3,3,4,4,5,5,6,7,7};
        ranks_first_lvl['S'] = {0,0,0,0,0,1,1,1,1,1,1,1,1,2,3};
        ranks_first_lvl['I'] = {0,0,0,0,1,1,1,1,1,1,1,2,2,2,2};
        ranks_first_lvl['E'] = {0,0,0,0,0,0,1,1,2,2,3,3,3,3,3};
        ranks_first_lvl['B'] = {0,1,2,3,4,5,5,7,7,9,9,11,12,12,13};
        ranks_first_lvl['H'] = {0,0,0,0,0,1,1,1,1,1,1,1,1,2,2};
        ranks_first_lvl['F'] = {0,0,0,0,0,0,1,1,1,1,1,1,1,2,2};
        ranks_pointer_first_lvl['N'] = {2,2,0,2,0,0,0};
        ranks_pointer_first_lvl['W'] = {0,0,0,0,0,0,0};
        ranks_pointer_first_lvl['O'] = {0,0,0,0,0,0,0};
        ranks_pointer_first_lvl['T'] = {0,0,0,0,0,1,0};
        ranks_pointer_first_lvl['S'] = {0,0,0,0,0,0,0};
        ranks_pointer_first_lvl['I'] = {0,0,0,0,1,0,0};
        ranks_pointer_first_lvl['E'] = {0,0,0,0,0,0,0};
        ranks_pointer_first_lvl['B'] = {0,0,0,0,1,1,0};
        ranks_pointer_first_lvl['H'] = {0,0,0,0,0,0,0};
        ranks_pointer_first_lvl['F'] = {0,0,0,0,0,0,0};

        ranks_second_lvl['N'] = {0,2,0,0,0,2,0,0,0,0,0,0,0,0,0,1};
        ranks_second_lvl['W'] = {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0};
        ranks_second_lvl['O'] = {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0};
        ranks_second_lvl['T'] = {0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0};
        ranks_second_lvl['S'] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
        ranks_second_lvl['I'] = {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0};
        ranks_second_lvl['E'] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        ranks_second_lvl['B'] = {0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0};
        ranks_second_lvl['H'] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        ranks_second_lvl['F'] = {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0};
    }

    void TearDown() override {delete lpf_bt;}
    BV_BlockTree_lpf_heuristic<uint8_t, int32_t>*  lpf_bt{};
    std::vector<uint8_t> vec;
    std::unordered_map<uint8_t, std::vector<int>> ranks_first_lvl;
    std::unordered_map<uint8_t, std::vector<int>> ranks_second_lvl;
    std::unordered_map<uint8_t, std::vector<int>> ranks_pointer_first_lvl;
    std::unordered_map<uint8_t, std::vector<int>> ranks_pointer_second_lvl;
};
TEST_F(btTest, accessTest) {
    // test access
    for (int i = 0; i < vec.size(); i++) {
        ASSERT_EQ(lpf_bt->access(i), vec[i]);
    }
}
TEST_F(btTest, firstLevelPrefixRanksN) {
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_['N']][0].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_['N']][0][i], ranks_first_lvl['N'][i]);
    }
}
TEST_F(btTest, firstLevelPrefixRanksW) {
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_['W']][0].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_['W']][0][i], ranks_first_lvl['W'][i]);
    }
}
TEST_F(btTest, firstLevelPrefixRanksO) {
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_['O']][0].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_['O']][0][i], ranks_first_lvl['O'][i]);
    }
}
TEST_F(btTest, firstLevelPrefixRanksT) {
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_['T']][0].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_['T']][0][i], ranks_first_lvl['T'][i]);
    }
}
TEST_F(btTest, firstLevelPrefixRanksS) {
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_['S']][0].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_['S']][0][i], ranks_first_lvl['S'][i]);
    }
}
TEST_F(btTest, firstLevelPrefixRanksI) {
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_['I']][0].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_['I']][0][i], ranks_first_lvl['I'][i]);
    }
}
TEST_F(btTest, firstLevelPrefixRanksE) {
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_['E']][0].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_['E']][0][i], ranks_first_lvl['E'][i]);
    }
}
TEST_F(btTest, firstLevelPrefixRanksB) {
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_['B']][0].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_['B']][0][i], ranks_first_lvl['B'][i]);
    }
}
TEST_F(btTest, firstLevelPrefixRanksH) {
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_['H']][0].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_['H']][0][i], ranks_first_lvl['H'][i]);
    }
}
TEST_F(btTest, firstLevelPrefixRanksF) {
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_['F']][0].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_['F']][0][i], ranks_first_lvl['F'][i]);
    }
}
TEST_F(btTest, sndLevelPrefixRanksN) {
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_['N']][1].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_['N']][1][i], ranks_second_lvl['N'][i]);
    }
}
TEST_F(btTest, sndLevelPrefixRanksW) {
    uint8_t c = 'W';
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1][i], ranks_second_lvl[c][i]);
    }
}
TEST_F(btTest, sndLevelPrefixRanksO) {
    uint8_t c = 'O';
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1][i], ranks_second_lvl[c][i]);
    }
}
TEST_F(btTest, sndLevelPrefixRanksT) {
    uint8_t c = 'T';
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1][i], ranks_second_lvl[c][i]);
    }
}
TEST_F(btTest, sndLevelPrefixRanksS) {
    uint8_t c = 'S';
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1][i], ranks_second_lvl[c][i]);
    }
}
TEST_F(btTest, sndLevelPrefixRanksI) {
    uint8_t c = 'I';
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1][i], ranks_second_lvl[c][i]);
    }
}
TEST_F(btTest, sndLevelPrefixRanksE) {
    uint8_t c = 'E';
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1][i], ranks_second_lvl[c][i]);
    }
}
TEST_F(btTest, sndLevelPrefixRanksB) {
    uint8_t c = 'B';
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1][i], ranks_second_lvl[c][i]);
    }
}
TEST_F(btTest, sndLevelPrefixRanksH) {
    uint8_t c = 'H';
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1][i], ranks_second_lvl[c][i]);
    }
}
TEST_F(btTest, sndLevelPrefixRanksF) {
    uint8_t c = 'F';
    // test access
    for (int i = 0; i < lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->c_ranks_[lpf_bt->chars_index_[c]][1][i], ranks_second_lvl[c][i]);
    }
}
TEST_F(btTest, firstLevelPointerRanksN) {
    uint8_t c = 'N';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i], ranks_pointer_first_lvl[c][i]);
    }
}
TEST_F(btTest, firstLevelPointerRanksW) {
    uint8_t c = 'W';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i], ranks_pointer_first_lvl[c][i]);
    }
}
TEST_F(btTest, firstLevelPointerRanksO) {
    uint8_t c = 'O';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i], ranks_pointer_first_lvl[c][i]);
    }
}
TEST_F(btTest, firstLevelPointerRanksT) {
    uint8_t c = 'T';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i], ranks_pointer_first_lvl[c][i]);
    }
}
TEST_F(btTest, firstLevelPointerRanksS) {
    uint8_t c = 'S';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i], ranks_pointer_first_lvl[c][i]);
    }
}
TEST_F(btTest, firstLevelPointerRanksI) {
    uint8_t c = 'I';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i], ranks_pointer_first_lvl[c][i]);
    }
}
TEST_F(btTest, firstLevelPointerRanksE) {
    uint8_t c = 'E';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i], ranks_pointer_first_lvl[c][i]);
    }
}
TEST_F(btTest, firstLevelPointerRanksB) {
    uint8_t c = 'B';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i], ranks_pointer_first_lvl[c][i]);
    }
}
TEST_F(btTest, firstLevelPointerRanksH) {
    uint8_t c = 'H';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i], ranks_pointer_first_lvl[c][i]);
    }
}
TEST_F(btTest, firstLevelPointerRanksF) {
    uint8_t c = 'F';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i], ranks_pointer_first_lvl[c][i]);
    }
}
TEST_F(btTest, sndLevelPointerRanksN) {
    uint8_t c = 'N';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1][i], 0) << "at " << i << " " << lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i];
    }
}
TEST_F(btTest, sndLevelPointerRanksW) {
    uint8_t c = 'W';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1][i], 0) << "at " << i << " " << lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i];
    }
}
TEST_F(btTest, sndLevelPointerRanksO) {
    uint8_t c = 'O';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1][i], 0) << "at " << i << " " << lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i];
    }
}
TEST_F(btTest, sndLevelPointerRanksT) {
    uint8_t c = 'T';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1][i], 0) << "at " << i << " " << lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i];
    }
}
TEST_F(btTest, sndLevelPointerRanksS) {
    uint8_t c = 'S';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1][i], 0) << "at " << i << " " << lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i];
    }
}
TEST_F(btTest, sndLevelPointerRanksI) {
    uint8_t c = 'I';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1][i], 0) << "at " << i << " " << lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i];
    }
}
TEST_F(btTest, sndLevelPointerRanksE) {
    uint8_t c = 'E';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1][i], 0) << "at " << i << " " << lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i];
    }
}
TEST_F(btTest, sndLevelPointerRanksB) {
    uint8_t c = 'B';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1][i], 0) << "at " << i << " " << lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i];
    }
}
TEST_F(btTest, sndLevelPointerRanksH) {
    uint8_t c = 'H';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1][i], 0) << "at " << i << " " << lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i];
    }
}
TEST_F(btTest, sndLevelPointerRanksF) {
    uint8_t c = 'F';
    // test access
    for (int i = 0; i < lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1].size(); i++) {
        ASSERT_EQ(lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][1][i], 0) << "at " << i << " " << lpf_bt->pointer_c_ranks_[lpf_bt->chars_index_[c]][0][i];
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
    auto*  lpf_bt = new BV_BlockTree_lpf_heuristic<uint8_t, int32_t>(vec, 2, 1, 15, lpf, lpf_ptr, lz);
    lpf_bt->add_rank_support();
    std::vector<uint8_t> chars = {'N', 'B', 'O', 'T', 'W' , 'I', 'S', 'H', 'F', 'E'};
    // test rank
    for (auto c: chars) {
        int count = 0;
        for (int i = 0; i < vec.size(); i++) {
            if (vec[i] == c) {
                count++;
            }
            ASSERT_EQ(lpf_bt->rank(c, i), count) << " rank and count differ at index " << i << " for " << c;

        }
    }
}
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}