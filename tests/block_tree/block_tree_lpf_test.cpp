/*******************************************************************************
 * This file is part of pasta::block_tree
 *
 * Copyright (C) 2022 Daniel Meyer
 * Copyright (C) 2023 Florian Kurpicz <florian@kurpicz.org>
 *
 * pasta::block_tree is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * pasta::block_tree is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with pasta::block_tree.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <random>
#include <vector>

#include <gtest/gtest.h>

#include <pasta/block_tree/construction/block_tree_lpf.hpp>
#include <pasta/block_tree/utils/lpf_array.hpp>

class BlockTreeLPFTest : public ::testing::Test {

protected:

  std::vector<uint8_t> text;
  std::vector<uint8_t> gappy_alphabet_text;

  pasta::BlockTreeLPF<uint8_t, int32_t>* bt;
  pasta::BlockTreeLPF<uint8_t, int32_t>* gappy_alphabet_bt;
  
  void SetUp() override {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint8_t> dist(0, 15);

    size_t const string_length = 100000;
    text.resize(string_length);
	gappy_alphabet_text.resize(string_length);
    for (size_t i = 0; i < text.size(); ++i) {
      text[i] = dist(gen);
	  gappy_alphabet_text[i] = 2*dist(gen);
    }

    bt = pasta::make_block_tree_lpf<uint8_t, int32_t>(text, 2, 1, true);
    bt->add_rank_support();
	gappy_alphabet_bt = pasta::make_block_tree_lpf<uint8_t, int32_t>(gappy_alphabet_text, 2, 1, true);
    gappy_alphabet_bt->add_rank_support();
  }

  void TearDown() override {
    delete bt;
	delete gappy_alphabet_bt;
  }

};

TEST_F(BlockTreeLPFTest, access) {
  for (size_t i = 0; i < text.size(); ++i) {
    ASSERT_EQ(bt->access(i), text[i]);
  }
}

TEST_F(BlockTreeLPFTest, access_gappy_alphabet) {
  for (size_t i = 0; i < text.size(); ++i) {
    ASSERT_EQ(gappy_alphabet_bt->access(i), gappy_alphabet_text[i]);
  }
}

TEST_F(BlockTreeLPFTest, rank) {
  std::array<size_t, 256> hist = {0};

  for (size_t i = 0; i < text.size() - 1; ++i) {
    ++hist[text[i]];
    ASSERT_EQ(bt->rank(text[i], i), hist[text[i]]);
  }
}

TEST_F(BlockTreeLPFTest, select) {
  std::array<size_t, 256> hist = {0};

  for (size_t i = 0; i < text.size() - 1; ++i) {
    ++hist[text[i]];
    ASSERT_EQ(bt->select(text[i], hist[text[i]]), i);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

/******************************************************************************/
