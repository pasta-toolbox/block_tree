/*******************************************************************************
 * This file is part of pasta::block_tree
 *
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

#include <cstdint>
#include <iostream>
#include <random>

#include <pasta/block_tree/construction/block_tree_lpf.hpp>

int32_t main()
{

  // Generate some random text
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<uint8_t> dist(0, 15);

  size_t const string_length = 100000;
  std::vector<uint8_t> text; 
  text.resize(string_length);
  for (size_t i = 0; i < text.size(); ++i) {
    text[i] = dist(gen);
  }

  // Build the block tree for the random text
  auto* bt = pasta::make_block_tree_lpf<uint8_t, int32_t>(text, 2, 1, true);

  // Use the block tree to access individual characters of the text
  std::cout << "# Acces" << "\n";
  for (size_t i = 0; i < 100; ++i) {
    std::cout << bt->access(i) << ", ";
  }
  std::cout << "\n";

  // Add additional rank and select support
  bt->add_rank_support();

  // Get the rank of the first character for the first 10 characters
  std::cout << "# Rank" << "\n";
  for (size_t i = 0; i < 15; ++i) {
    std::cout << bt->rank(0, i) << ", ";
  }
  std::cout << "\n";

  // Get the positions of the first 5 in the text
  std::cout << "# Select" << "\n";
  std::cout << bt->select(5, 1) << "\n";

  // Clean-up
  delete bt;
  
  return 0;
}



/******************************************************************************/
