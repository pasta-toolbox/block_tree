#pragma once

#include "pasta/bit_vector/bit_vector.hpp"
#include "pasta/block_tree/utils/MersenneHash.hpp"
#include "pasta/block_tree/utils/MersenneRabinKarp.hpp"

#include <cstddef>
#include <cstdint>
#include <list>
#include <memory>
#include <vector>

/// @brief Utilities for the construction algorithms using the sharded hash map
namespace pasta {

namespace internal::sharded {

__extension__ typedef unsigned __int128 uint128_t;

/// @brief A marker for a block that has no earlier occurrence
inline constexpr int64_t NO_EARLIER_OCC = -1;

/// @brief A marker for a block that has been pruned
inline constexpr int64_t PRUNED = -2;

/// @brief For some block size (in bytes) i, return the number of trailing
///   zeros in a 64 bit integer when zeroing out characters that are not part
///   of the block.
constexpr static uint64_t MASK_TRAILING_ZEROS[9] =
    {64, 56, 48, 40, 32, 24, 16, 8, 0};

/// @brief Masks used for the identity hash. These depend on endianness
constexpr static std::array<uint64_t, 9> masks() {
  if constexpr (std::endian::native == std::endian::big) {
    return {0,
            static_cast<uint64_t>(~0) << MASK_TRAILING_ZEROS[1],
            static_cast<uint64_t>(~0) << MASK_TRAILING_ZEROS[2],
            static_cast<uint64_t>(~0) << MASK_TRAILING_ZEROS[3],
            static_cast<uint64_t>(~0) << MASK_TRAILING_ZEROS[4],
            static_cast<uint64_t>(~0) << MASK_TRAILING_ZEROS[5],
            static_cast<uint64_t>(~0) << MASK_TRAILING_ZEROS[6],
            static_cast<uint64_t>(~0) << MASK_TRAILING_ZEROS[7],
            static_cast<uint64_t>(~0) << MASK_TRAILING_ZEROS[8]};
  } else {
    return {0,
            static_cast<uint64_t>(~0) >> MASK_TRAILING_ZEROS[1],
            static_cast<uint64_t>(~0) >> MASK_TRAILING_ZEROS[2],
            static_cast<uint64_t>(~0) >> MASK_TRAILING_ZEROS[3],
            static_cast<uint64_t>(~0) >> MASK_TRAILING_ZEROS[4],
            static_cast<uint64_t>(~0) >> MASK_TRAILING_ZEROS[5],
            static_cast<uint64_t>(~0) >> MASK_TRAILING_ZEROS[6],
            static_cast<uint64_t>(~0) >> MASK_TRAILING_ZEROS[7],
            static_cast<uint64_t>(~0) >> MASK_TRAILING_ZEROS[8]};
  }
}

/// @brief Masks for identity hashes for a block size i (in bytes)
inline constexpr std::array<uint64_t, 9> HASH_MASKS = masks();

/// @brief Base of the polynomial used for the Rabin-Karp hasher
inline constexpr size_t SIGMA = 256;

/// @brief The exponent of the mersenne prime used for the Rabin-Karp hasher
inline constexpr uint8_t PRIME_EXPONENT = 107;

/// @brief A mersenne prime used for the Rabin-Karp hasher
inline constexpr uint128_t PRIME = pasta::mersenne_prime<PRIME_EXPONENT>();

/// @brief Determine whether to use a Rabin-Karp hash for hashing text windows
/// or just use the block's content itself as a hash, stored in an integer.
enum class UseHash {
  /// @brief Use a Rabin-Karp hash
  RABIN_KARP,
  /// @brief Use the block's content as a hash
  IDENTITY
};

/// @brief Contains data about a block tree level under construction
template <std::signed_integral size_type, typename RankSelect>
struct LevelData {
  /// @brief Contains a 1 for each internal block (= block with children)
  ///   and a 0 for each block that has a back pointer
  std::unique_ptr<pasta::BitVector> is_internal;
  /// @brief Rank data structure for is_internal
  std::unique_ptr<RankSelect> is_internal_rank;
  /// @brief The block from which a back block is copying
  std::unique_ptr<std::vector<size_type>> pointers;
  /// @brief The offset into the block from which the back block is copying
  std::unique_ptr<std::vector<size_type>> offsets;
  /// @brief The number of back blocks pointing to the block
  std::unique_ptr<std::vector<size_type>> counters;
  /// @brief Block start indices
  std::unique_ptr<std::vector<size_type>> block_starts;
  /// @brief The block size on this level
  int64_t block_size;
  /// @brief The index of the current level.
  /// First level is 0, second level is 1 etc.
  int64_t level_index;
  /// @brief The number of blocks on the current level
  int64_t num_blocks;

  LevelData(const int64_t level_index_,
            const int64_t block_size_,
            const int64_t num_blocks_)
      : is_internal(nullptr),
        is_internal_rank(nullptr),
        pointers(new std::vector<size_type>()),
        offsets(new std::vector<size_type>()),
        counters(new std::vector<size_type>()),
        block_starts(new std::vector<size_type>()),
        block_size(block_size_),
        level_index(level_index_),
        num_blocks(num_blocks_) {}

  /// @brief Checks whether a block is adjacent in the text
  ///   to its successor on this level
  [[nodiscard]] bool next_is_adjacent(size_t i) const {
    return (*block_starts)[i] + block_size == (*block_starts)[i + 1];
  }
};

/// @brief Contains data about the occurrences of a hashed block pair
template <std::signed_integral size_type>
struct PairOccurrences {
  /// @brief The first block in the text in which the content appears
  size_type first_occ_block;
  /// @brief A list of block indices in which the content of the hashed block
  /// pair appears
  ///
  /// We're using an std::list here instead of an std::vector, since the
  /// reallocation upon insertion lead to issues during parallel access, when
  /// another thread tries to access the vector during reallocation.
  std::vector<size_type> occurrences;

  /// @brief Initialize the occurrences of a hashed block pair.
  ///
  /// Note, that this only sets the first occurrence to the given block index,
  /// but does not add it to the occurrences list.
  /// @param first_occ_block_ The block index of the pair's first block.
  inline explicit PairOccurrences(size_type first_occ_block_)
      : first_occ_block(first_occ_block_),
        occurrences() {}

  PairOccurrences(PairOccurrences&&) noexcept = default;
  PairOccurrences& operator=(PairOccurrences&&) = default;

  /// @brief Add a block index to the occurrences.
  /// @param block_index The block index to add to the occurrences.
  void add_block_pair(size_type block_index) {
    occurrences.push_back(block_index);
  }

  /// @brief If the given block index is an earlier occurrence, update it
  /// @param block_index The block index of an occurrence
  void update(size_type block_index) {
    first_occ_block = std::min<size_type>(first_occ_block, block_index);
  }
};

/// @brief Contains data about the occurrences of a hashed block
template <std::signed_integral size_type>
struct BlockOccurrences {
  /// @brief Represents the first occurrence of a block
  struct [[gnu::packed]] FirstOccurrence {
    /// @brief Block index of the first occurrence of the block's content
    int64_t block : 40;
    /// @brief The offset into the block at which that first occurrence occurs
    int32_t offset : 24;

    inline FirstOccurrence(int64_t first_occ_block_,
                           int32_t first_occ_offset_)
        : block(first_occ_block_),
          offset(first_occ_offset_) {}
  };

  static_assert(std::atomic<FirstOccurrence>::is_always_lock_free,
                "first occurrence must be able to be atomically updated");
  static_assert(sizeof(FirstOccurrence) == 8,
                "should be size of computer word");

  // @brief The block index and offset of the first occurrence of this block's
  //   content
  std::atomic<FirstOccurrence> first_occ;

  /// @brief A list of block indices in which the content of the hashed block
  ///   occurs
  std::vector<size_type> occurrences;

  /// @brief Initialize the occurrences of a hashed block.
  ///
  /// Note, that this only sets the first occurrence to the given block index,
  /// but does not add it to the occurrences list.
  /// @param first_occ_block_ The block index of the block's first occurrence.
  explicit BlockOccurrences(size_type first_occ_block_)
      : first_occ({first_occ_block_, 0}),
        occurrences() {}

  BlockOccurrences(const BlockOccurrences& other)
      : first_occ(other.first_occ.load()),
        occurrences(other.occurrences) {}

  BlockOccurrences(BlockOccurrences&& other) noexcept
      : first_occ(other.first_occ.load()),
        occurrences(std::move(other.occurrences)) {}

  ~BlockOccurrences() = default;

  BlockOccurrences& operator=(BlockOccurrences&& other) noexcept {
    first_occ = other.first_occ.load();
    occurrences = std::move(other.occurrences);
    return *this;
  }

  /// @brief Add a block index to the occurrences.
  /// @param block_index The block index to add to the occurrences.
  void add_block(size_type block_index) {
    occurrences.push_back(block_index);
  }

  /// @brief If the given block index and offset are an earlier occurrence,
  ///   update them
  /// @param block_index The block index of an occurrence
  /// @param block_offset The offset of that occurrence
  void update(size_type block_index, size_type block_offset) {
    FirstOccurrence prev_first_occ = this->first_occ.load();
    FirstOccurrence set(block_index, block_offset);
    while (block_index < prev_first_occ.block &&
           !first_occ.compare_exchange_weak(prev_first_occ, set)) {
    }
  }
};

/// @brief An update function for the sharded hash map that updates the
///   occurrences of a hashed block pair
template <std::integral input_type, std::signed_integral size_type>
struct UpdatePairOccurrences {
  /// @brief The block index to add to the occurrences
  using InputValue = size_type;
  /// @brief Update the occurrences of a hashed block pair by adding the new
  ///   block index and updating the first occurrence if needed
  /// @param occurrences A reference to the occurrences in the map
  /// @param input_value The new block index to add to the occurrences
  inline static void update(const MersenneHash<input_type>&,
                            PairOccurrences<size_type>& occurrences,
                            InputValue&& input_value) {
    occurrences.add_block_pair(input_value);
    occurrences.update(input_value);
  }

  /// @brief Initialize the occurrences of a hashed block pair
  /// @param input_value The block index of the pair's first block
  /// @return The initialized occurrences only containing the given block pair
  inline static PairOccurrences<size_type> init(const MersenneHash<input_type>&,
                                                InputValue&& input_value) {
    PairOccurrences occurrences(input_value);
    occurrences.add_block_pair(input_value);
    occurrences.update(input_value);
    return occurrences;
  }
};

/// @brief An update function for the sharded hash map that updates the
///   occurrences of a hashed block
template <std::integral input_type, std::signed_integral size_type>
struct UpdateBlockOccurrences {
  /// @brief A pair of the block index
  ///   and offset of the first occurrence of a block
  using InputValue = std::pair<size_type, size_type>;

  /// @brief Update the occurrences of a hashed block by adding the new
  ///   block index and offset and updating the first occurrence if needed
  /// @param occurrences A reference to the occurrences in the map
  /// @param input_value The new block index and offset to add to the
  ///   occurrences
  inline static void update(const MersenneHash<input_type>&,
                            BlockOccurrences<size_type>& occurrences,
                            InputValue&& input_value) {
    occurrences.add_block(input_value.first);
    occurrences.update(input_value.first, input_value.second);
  }

  /// @brief Initialize the occurrences of a hashed block.
  /// @param input_value A pair of the block index and offset of one of the
  ///   block's occurrences
  /// @return The initialized occurrences only containing the given block
  inline static BlockOccurrences<size_type>
  init(const MersenneHash<input_type>&, InputValue&& input_value) {
    BlockOccurrences occurrences(input_value.first);
    occurrences.add_block(input_value.first);
    occurrences.update(input_value.first, input_value.second);
    return occurrences;
  }
};

/// @brief A mixing functions to provide better avalanching to intermediate hash
/// values.
///
/// https://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html
constexpr uint64_t mix_select(uint64_t key) {
  key ^= (key >> 31);
  key *= 0x7fb5d329728ea185;
  key ^= (key >> 27);
  key *= 0x81dadef4bc2dd44d;
  key ^= (key >> 33);
  return key;
}

/// @brief Returns the ceiling of x / y for x > 0;
///
/// https://stackoverflow.com/questions/2745074/fast-ceiling-of-an-integer-division-in-c-c
size_t ceil_div(std::integral auto x, std::integral auto y) {
  return 1 + (static_cast<size_t>(x) - 1) / static_cast<size_t>(y);
}

} // namespace internal::sharded

///
/// @brief An update function for sharded maps which on update just overwrites
/// the value.
///
/// @tparam K The key type saved in the hash map.
/// @tparam V The value type saved in the hash map.
///
template <typename K, typename V>
struct Overwrite {
  using InputValue = V;

  inline static void update(K&, V& value, V&& input_value) {
    value = input_value;
  }

  inline static V init(K&, V&& input_value) {
    return input_value;
  }
};

///
/// @brief An update function for sharded maps which upon update does nothing
/// besides inserting the value if it doesn't exist.
///
/// @tparam K The key type saved in the hash map.
/// @tparam V The value type saved in the hash map.
///
template <typename K, typename V>
struct Keep {
  using InputValue = V;
  inline static void update(K&, V&, V&&) {}

  inline static V init(K&, V&& input_value) {
    return input_value;
  }
};

} // namespace pasta
