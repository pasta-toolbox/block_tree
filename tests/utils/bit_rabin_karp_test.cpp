
#include "pasta/block_tree/utils/MersenneHash.hpp"

#include <bitset>
#include <cstddef>
#include <cstdint>
#include <gtest/gtest.h>
#include <pasta/block_tree/utils/MersenneRabinKarp.hpp>
#include <random>

class BitRabinKarpTest : public ::testing::Test {
protected:

  static constexpr size_t unique_len = 100000;
  pasta::BitVector bv;

  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_int_distribution<uint8_t> dist(10);
    bv.resize(unique_len * 2 + 3);
    for (size_t i = 0; i < unique_len; ++i) {
      bv[i] = dist(gen) % 2 == 0;
    }
    // Offset it by some amount to check that even misaligned bit sequences
    // correctly match
    bv[unique_len] = false;
    bv[unique_len + 1] = false;
    bv[unique_len + 2] = false;
    for (size_t i = 0; i < unique_len; ++i) {
      bv[unique_len + 3 + i] = static_cast<bool>(bv[i]);
    }
  }

public:
  bool
  compare_ranges(const size_t s1, const size_t s2, const size_t len) const {
    for (size_t i = 0; i < len; ++i) {
      if (get_bit(s1 + i) != get_bit(s2 + i)) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] bool get_bit(const size_t bit_index) const {
    return bv[bit_index];
  }
};

TEST_F(BitRabinKarpTest, test_slice) {
  for (size_t i = 0; i < bv.size() - 1; i++) {
    uint64_t manual_slice = 0;
    for (size_t j = 0; j < 64; j++) {
      manual_slice |= (static_cast<uint64_t>(bv[i + j]) << j);
    }

    uint64_t slice = pasta::MersenneHash<bool>::slice_at(bv, i);
    ASSERT_EQ(manual_slice, slice)
        << " bit index " << i
        << "\nwith manual slice: " << std::bitset<64>(manual_slice)
        << "\nwith hash slice:   " << std::bitset<64>(slice);
  }
}

TEST_F(BitRabinKarpTest, test_eq) {
  constexpr std::array<size_t, 4> sizes = {13, 24, 59, 1220};

  for (const size_t size : sizes) {
    pasta::MersenneRabinKarp<bool, int32_t> rk1(bv, 0, size, (1ULL << 61) - 1);
    pasta::MersenneRabinKarp<bool, int32_t> rk2(bv,
                                                unique_len + 3,
                                                size,
                                                (1ULL << 61) - 1);
    for (size_t i = 0; i < unique_len - size; ++i) {
      const auto h1 = rk1.current_hash();
      const auto h2 = rk2.current_hash();
      ASSERT_TRUE(h1 == h2) << "offset: " << i << " for half len " << unique_len
                            << " and size " << size;
      rk1.next();
      rk2.next();
    }
  }
}

TEST_F(BitRabinKarpTest, test_rnd) {
  constexpr std::array<size_t, 4> sizes = {13, 24, 59, 1220};

  for (const size_t size : sizes) {
    pasta::MersenneRabinKarp<bool, int32_t> rk1(bv, 0, size, (1ULL << 61) - 1);
    pasta::MersenneRabinKarp<bool, int32_t> rk2(bv,
                                                unique_len,
                                                size,
                                                (1ULL << 61) - 1);
    size_t offset = 0;
    for (size_t i = 0; i < unique_len - size; ++i) {
      const auto h1 = rk1.current_hash();
      const auto h2 = rk2.current_hash();
      ASSERT_EQ((h1 == h2), compare_ranges(offset, unique_len + offset, size))
          << "error at offset " << offset << " for half_len " << unique_len
          << " and window size " << size;
      rk1.next();
      rk2.next();
      offset++;
    }
  }
}
