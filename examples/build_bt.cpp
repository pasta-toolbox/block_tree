/*******************************************************************************
 * This file is part of pasta::block_tree
 *
 * Copyright (C) 2023 Etienne Palanga <etienne.palanga@tu-dortmund.de>
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

#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <pasta/bit_vector/support/flat_rank_select.hpp>
#include <syncstream>

#define REC_PAR_SHARDED
#define REC_DENSE_BIT

using SizeType = int64_t;

#if defined REC_BIT || defined REC_DENSE_BIT || defined REC_PAR_SHARDED
constexpr size_t RECURSION_LEVELS = 0;
#else
constexpr size_t RECURSION_LEVELS = 0;
#endif

#if defined REC_DENSE_BIT
#  include <pasta/block_tree/construction/rec_dense_bit_block_tree_sharded.hpp>
using BBT =
    pasta::RecursiveDenseBitBlockTreeSharded<SizeType, RECURSION_LEVELS>;
#  define BIT_ALGO_NAME "rec_dense_bit"
#elif defined REC_BIT
#  include <pasta/block_tree/construction/rec_bit_block_tree_sharded.hpp>
using BBT = pasta::RecursiveBitBlockTreeSharded<SizeType, RECURSION_LEVELS>;
#  define BIT_ALGO_NAME "rec_bit"
#endif

#ifdef FP
#  include <pasta/block_tree/construction/block_tree_fp.hpp>
std::unique_ptr<pasta::BlockTreeFP<uint8_t, SizeType>>
make_bt(std::vector<uint8_t>& text,
        const size_t arity,
        const size_t leaf_length,
        const size_t,
        const size_t) {
  ;
  return std::unique_ptr<pasta::BlockTreeFP<uint8_t, SizeType>>(
      pasta::make_block_tree_fp<uint8_t, SizeType>(text, arity, leaf_length));
}
#  define ALGO_NAME "fp"
#elif defined FP2
#  include <pasta/block_tree/construction/block_tree_fp2_seq.hpp>
std::unique_ptr<pasta::BlockTreeFP2<uint8_t, SizeType>>
make_bt(std::vector<uint8_t>& text,
        const size_t arity,
        const size_t leaf_length,
        const size_t,
        const size_t) {
  ;
  return std::make_unique<pasta::BlockTreeFP2<uint8_t, SizeType>>(text,
                                                                  arity,
                                                                  1,
                                                                  leaf_length);
}
#  define ALGO_NAME "fp2"
#elif defined LPF
#  include <pasta/block_tree/construction/block_tree_lpf.hpp>
std::unique_ptr<pasta::BlockTreeLPF<uint8_t, SizeType>>
make_bt(std::vector<uint8_t>& text,
        const size_t arity,
        const size_t leaf_length,
        const size_t threads,
        const size_t) {
  ;
  return std::unique_ptr<pasta::BlockTreeLPF<uint8_t, SizeType>>(
      pasta::make_block_tree_lpf_parallel<uint8_t, SizeType>(text,
                                                             arity,
                                                             leaf_length,
                                                             true,
                                                             threads));
}
#  define ALGO_NAME "lpf"
#elif defined PAR_SHARDED
#  include <pasta/block_tree/construction/block_tree_fp_par_sharded.hpp>
std::unique_ptr<pasta::BlockTreeFPParSharded<uint8_t, SizeType>>
make_bt(std::vector<uint8_t>& text,
        const size_t arity,
        const size_t leaf_length,
        const size_t threads,
        const size_t) {
  ;
  return std::make_unique<pasta::BlockTreeFPParSharded<uint8_t, SizeType>>(
      text,
      arity,
      1,
      leaf_length,
      threads);
}
#  define ALGO_NAME "shard"
#elif defined PAR_SHARDED_SYNC
#  include <pasta/block_tree/construction/block_tree_fp_par_sync_sharded.hpp>
std::unique_ptr<pasta::BlockTreeFPParShardedSync<uint8_t, SizeType>>
make_bt(std::vector<uint8_t>& text,
        const size_t arity,
        const size_t leaf_length,
        const size_t threads,
        const size_t queue_size) {
  ;
  return std::make_unique<pasta::BlockTreeFPParShardedSync<uint8_t, SizeType>>(
      text,
      arity,
      1,
      leaf_length,
      threads,
      queue_size);
}
#  define ALGO_NAME "shard_sync"
#elif defined PAR_SHARDED_SYNC_SMALL
#  include <pasta/block_tree/construction/rec_block_tree_sharded.hpp>
std::unique_ptr<pasta::RecursiveBlockTreeSharded<uint8_t, SizeType, 0>>
make_bt(std::vector<uint8_t>& text,
        const size_t arity,
        const size_t leaf_length,
        const size_t threads,
        const size_t queue_size) {
  return std::make_unique<
      pasta::RecursiveBlockTreeSharded<uint8_t, SizeType, 0>>(text,
                                                              arity,
                                                              1,
                                                              leaf_length,
                                                              threads,
                                                              queue_size);
}
#  define ALGO_NAME "shard_sync_small"
#elif defined REC_PAR_SHARDED
#  include <pasta/block_tree/construction/rec_block_tree_sharded.hpp>
std::unique_ptr<
    pasta::RecursiveBlockTreeSharded<uint8_t, int64_t, RECURSION_LEVELS>>
make_bt(std::vector<uint8_t>& text,
        const size_t arity,
        const size_t leaf_length,
        const size_t threads,
        const size_t queue_size) {
  return std::make_unique<
      pasta::RecursiveBlockTreeSharded<uint8_t, int64_t, RECURSION_LEVELS>>(
      text,
      arity,
      1,
      leaf_length,
      threads,
      queue_size);
}
#  define ALGO_NAME "rec_shard"
#elif defined PAR_PHMAP
#  include <pasta/block_tree/construction/block_tree_fp_par_phmap.hpp>
std::unique_ptr<pasta::BlockTreeFPParPH<uint8_t, SizeType>>
make_bt(std::vector<uint8_t>& text,
        const size_t arity,
        const size_t leaf_length,
        const size_t threads,
        const size_t) {
  ;
  return std::make_unique<pasta::BlockTreeFPParPH<uint8_t, SizeType>>(
      text,
      arity,
      1,
      leaf_length,
      threads);
}
#  define ALGO_NAME "par_map"
#elif defined PAR_PARLAY
#  include <pasta/block_tree/construction/block_tree_fp_par_parlay.hpp>
std::unique_ptr<pasta::BlockTreeFPParParlay<uint8_t, SizeType>>
make_bt(std::vector<uint8_t>& text,
        const size_t arity,
        const size_t leaf_length,
        const size_t threads,
        const size_t) {
  ;
  return std::make_unique<pasta::BlockTreeFPParParlay<uint8_t, SizeType>>(
      text,
      arity,
      1,
      leaf_length,
      threads);
}
#  define ALGO_NAME "par_parlay"
#endif

#include <sstream>
#include <string>
#include <tlx/cmdline_parser.hpp>

using Clock = std::chrono::high_resolution_clock;
using TimePoint = Clock::time_point;
using Duration = Clock::duration;

int main(int argc, char** argv) {
  using namespace pasta;

  tlx::CmdlineParser cp;
  cp.set_description(
      "Build a block tree for a given input text or bit vector.");

  std::string file;
  cp.add_param_string("file",
                      file,
                      "The path to the file which to build a block tree from");

  size_t arity = 0;
  cp.add_param_size_t("arity", arity, "The arity of the block tree");
  size_t leaf_length = 0;
  cp.add_param_size_t(
      "leaf",
      leaf_length,
      "The maximum number of characters saved verbatim per leaf block.");

  size_t threads = 1;
  cp.add_size_t('t',
                "threads",
                threads,
                "The number of threads to use for parallel algorithms (ignored "
                "for sequential algorithms)");
  size_t queue_size = 1024;
  cp.add_size_t(
      'q',
      "queue_size",
      queue_size,
      "The size of each thread's queue used for sharded hash map algorithms");

  bool make_bv = false;
  cp.add_bool('b',
              "bitvec",
              make_bv,
              "Whether to interpret the input as a bitvector. If \"-o\" is not "
              "used, each byte of the input file will be interpreted as 8 bits "
              "of the bit vector respectively. In each byte, the least "
              "significant bit is index 0.");

  std::string one_chars;
  cp.add_string(
      'o',
      "one_chars",
      one_chars,
      "A string to be used with the \"-b\" flag. If used, each character of "
      "the input represents one bit of the bit vector. It will be a 1 if this "
      "parameter string contains the respective character, 0 otherwise.");

  bool verify = false;
  cp.add_bool(
      'v',
      "verify",
      verify,
      "Verify whether all queries on the produced block tree are correct.");

  if (!cp.process(argc, argv)) {
    return 1;
  }

#ifdef BT_DBG
  std::cout << "building block tree with parameters:"
            << "\narity: " << arity << "\nmax leaf length: " << leaf_length
            << "\nusing " << threads << " threads" << std::endl;
#endif

  if (!std::filesystem::exists(file)) {
    std::cerr << "File does not exist" << std::endl;
    exit(1);
  }

  std::unique_ptr<pasta::BitVector> bv;
  std::vector<uint8_t> text;
  {
    const size_t input_size =
        std::ifstream(file, std::ios::binary | std::ios::ate).tellg();
    std::ifstream t(file);
    if (make_bv) {
      if (one_chars.empty()) {
        // Interpret each character as 8 bits
        bv = std::make_unique<pasta::BitVector>(input_size * 8);
        std::span<std::byte> bytes = std::as_writable_bytes(bv->data());
        for (size_t i = 0; i < input_size; ++i) {
          char next_byte;
          t >> next_byte;
          bytes[i] = std::byte{static_cast<uint8_t>(next_byte)};
        }
      } else {
        // Interpret each character as a bit
        bv = std::make_unique<pasta::BitVector>(input_size);
        std::array<bool, 256> is_one{};
        for (char c : one_chars) {
          is_one[static_cast<uint8_t>(c)] = true;
        }
        for (size_t i = 0; i < input_size; ++i) {
          char next_byte;
          t >> next_byte;
          (*bv)[i] = is_one[static_cast<uint8_t>(next_byte)];
        }
      }
    } else {
      std::stringstream buffer;
      buffer << t.rdbuf();
      std::string input = buffer.str();
      text = std::vector<uint8_t>(input.begin(), input.end());
    }
  }

  std::cout << "RESULT"
            << " file=" << std::filesystem::path(file).filename().string()
            << " threads=" << threads << " arity=" << arity
            << " leaf_length=" << leaf_length;
  if (make_bv) {
    std::cout << " input_size=" << bv->size() / 8;
  } else {
    std::cout << " input_size=" << text.size();
  }
  TimePoint now = Clock::now();

  if (make_bv) {
    std::cout << " algo=" << BIT_ALGO_NAME;
    // Make bit vector block tree

    /*
        auto bt = std::make_unique<
            RecursiveBitBlockTreeSharded<SizeType, RECURSION_LEVELS>>(*bv,
                                                                     arity,
                                                                     1,
                                                                     leaf_length,
                                                                     threads,
                                                                     queue_size);
    */
    auto bt =
        std::make_unique<BBT>(*bv, arity, 1, leaf_length, threads, queue_size);

#ifdef BT_DBG
/*
    size_t cnt = 0;
    for (const auto& b : *bt->leaf_bits_) {
      if (b) {
        cnt++;
      }
    }
    std::cout << "num ones: " << cnt << "/" << bt->leaf_bits_->size() << " ("
              << static_cast<double>(cnt) * 100 / bt->leaf_bits_->size() << "%)"
              << std::endl;
              */
#endif
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                       Clock::now() - now)
                       .count();
    const size_t no_rs_space = bt->print_space_usage();
    bt->add_bit_rank_support();
    auto elapsed_rs = std::chrono::duration_cast<std::chrono::milliseconds>(
                          Clock::now() - now)
                          .count();
    const size_t rs_space = bt->print_space_usage();
    std::cout << " rec=" << RECURSION_LEVELS;
    std::cout << " time=" << elapsed << " space=" << no_rs_space
              << " time_rs=" << elapsed_rs << " space_rs=" << rs_space;
    std::cout << std::endl;

    if (!verify) {
      return 0;
    }

    std::cerr << "Start verification...\n";

    FlatRankSelect<> frs(*bv);

#if defined BT_INSTRUMENT && defined BT_DBG
    pasta::print_hash_data();
#endif
    std::cerr << "Access queries... " << std::flush;
#pragma omp parallel for
    for (size_t i = 0; i < bv->size(); ++i) {
      const bool c = bt->access(i);
      if (c != (*bv)[i]) {
        std::osyncstream(std::cerr)
            << "\nAccess error at position " << i
            << "\nExpected: " << std::boolalpha << (*bv)[i] << "\nActual: " << c
            << std::noboolalpha << std::endl;
        exit(1);
      }
    }
    std::cerr << "successful\n";
    std::cerr << "Rank 1 queries... " << std::flush;
#pragma omp parallel for
    for (size_t i = 0; i < bv->size(); i++) {
      const size_t bt_rank = bt->rank1(i);
      const size_t bv_rank = frs.rank1(i);

      if (bv_rank != bt_rank) {
        std::osyncstream(std::cerr)
            << "\nRank one error at position " << i << "\nExpected: " << bv_rank
            << "\nActual: " << bt_rank << std::endl;
        throw std::runtime_error("oof");
      }
    }

    std::cerr << "successful\n";
    std::cerr << "Rank 0 queries... " << std::flush;
#pragma omp parallel for
    for (size_t i = 0; i < bv->size(); i++) {
      const size_t bt_rank = bt->rank0(i);
      const size_t bv_rank = frs.rank0(i);

      if (bv_rank != bt_rank) {
        std::osyncstream(std::cerr) << "\nRank zero error at position " << i
                                    << "\nExpected: " << bv_rank
                                    << "\nActual: " << bt_rank << std::endl;
        throw std::runtime_error("oof");
      }
    }

    const size_t num_zeros = frs.rank0(bv->size());
    const size_t num_ones = frs.rank1(bv->size());

    std::cerr << "successful\n";
    std::cerr << "Select 1 queries... " << std::flush;
#pragma omp parallel for
    for (size_t i = 1; i <= num_ones; i++) {
      const size_t bv_rank = frs.select1(i);
      const size_t bt_rank = bt->select1(i);
      if (bv_rank != bt_rank) {
        std::osyncstream(std::cerr) << "\nSelect one error at position " << i
                                    << "\nExpected: " << bv_rank
                                    << "\nActual: " << bt_rank << std::endl;
        throw std::runtime_error("oof");
      }
    }

    std::cerr << "successful\n";
    std::cerr << "Select 0 queries... " << std::flush;
#pragma omp parallel for
    for (size_t i = 1; i <= num_zeros; i++) {
      const size_t bv_rank = frs.select0(i);
      const size_t bt_rank = bt->select0(i);
      if (bv_rank != bt_rank) {
        std::osyncstream(std::cerr) << "\nSelect zero error at position " << i
                                    << "\nExpected: " << bv_rank
                                    << "\nActual: " << bt_rank << std::endl;
        throw std::runtime_error("oof");
      }
    }
    std::cerr << "successful" << std::endl;
  } else {
    std::cout << " algo=" << ALGO_NAME;
    // Make text block tree
    auto bt = make_bt(text, arity, leaf_length, threads, queue_size);
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                       Clock::now() - now)
                       .count();

    std::cout << " rec=" << RECURSION_LEVELS;
    std::cout << " time=" << elapsed << " space=" << bt->print_space_usage();

    if (!verify) {
      std::cout << std::endl;
      return 0;
    }

    // std::cerr << "Start verification...\n";
#if defined BT_INSTRUMENT && defined BT_DBG
    pasta::print_hash_data();
#endif

    // std::cerr << "Access queries... " << std::flush;
#pragma omp parallel for
    for (size_t i = 0; i < text.size(); ++i) {
      const auto c = bt->access(i);
      if (c != text[i]) {
        #pragma omp critical
        std::cout << " verification=failed" << std::endl;
        std::exit(-1);
      }
    }
  }
  std::cout << " verification=passed" << std::endl;

  return 0;
}

/******************************************************************************/
