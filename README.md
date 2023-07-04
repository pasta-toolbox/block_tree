# pasta::block_tree

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![pasta::block_tree CI badge](https://github.com/pasta-toolbox/block_tree/actions/workflows/ctest.yml/badge.svg)](https://github.com/pasta-toolbox/block_tree/actions/workflows/ctest.yml)

This header-only library contains two highly configurable block tree implementations and is largely based on Daniel Meyer's [code](https://github.com/uqdwq/block_tree) developed during his Master's thesis.

If you use this code in a scientific context, please cite our paper.
```bibtex
@inproceedings{KoepplKM2023LPFBlockTree,
  author    = {Dominik Köppl and
               Florian Kurpicz and
               Daniel Meyer},
  title     = {Faster Block Tree Construction},
  booktitle = {Accepted at {ESA}},
  year      = {2023}
}
```

## Content
This repository contains a generic [interface](include/pasta/block_tree/block_tree.hpp) for block trees and two construction algorithms.
1. A [construction algorithm](include/pasta/block_tree/construction/block_tree_fp.hpp) based on Karp-Rabin fingerprints, which is basically a re-implementation of the original block tree implementation described by [Belazzougui et al. J. Comput. Syst. Sci. ’21](https://doi.org/10.1016/j.jcss.2020.11.002).
2. A [construction algorithm](include/pasta/block_tree/construction/block_tree_lpf.hpp) based on the longest previous factor array, which can be up to an order of magnitude faster than the first algorithm, see our [benchmarks](#markdown-header-benchmarks-and-tests-title) below.

### Easy to Use

This library has only a single dependency in the form of the Succinct Data Structure Library ([SDSL](https://github.com/simongog/sdsl-lite)).[^removeSDSL]
Please make sure you have this library installed, such that CMake can find it.
You can see in our [example](examples/block_tree_construction.cpp) below, that the block tree is easy to compute and easy to use.

[^removeSDSL]: We are working hard on removing this dependency.

```cpp
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
```

Note that this software is currently in early development and the interface might change during the following releases.

### Benchmarks and Tests
There exists an easy to use [benchmark repository][], which helps to compare the implementations contained in this repository with other block tree implementations.
The script reproduces the experiments conducted in our paper (including plotting the data).
Below, you can find some of the figures we present in the paper.

![Sequential construction time plot](https://raw.githubusercontent.com/pasta-toolbox/block_tree/main/docs/images/construction_time_repetitive_wo_rank_select_v0.1.0.png)

![Parallel construction time plot](https://raw.githubusercontent.com/pasta-toolbox/block_tree/main/docs/images/parallel_construction_time_repetitive_wo_rank_select_v0.1.0.png)

[benchmark repository]: https://github.com/pasta-toolbox/block_tree_experiments