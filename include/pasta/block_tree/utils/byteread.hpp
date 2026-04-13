#pragma once

#include "concepts.hpp"

namespace pasta {

/// @brief Swaps the bytes in an integer.
template <std::integral Int>
Int byteswap(Int& i) {
  switch (sizeof(Int)) {
    case 16:
      __bswap_16(i);
      break;
    case 32:
      __bswap_32(i);
      break;
    case 64:
      __bswap_64(i);
      break;
    default: {
    }
  }
  return i;
}

/// @brief Copies an integer from a pointer using the architecture's native
/// endianness.
/// @param ptr The pointer to copy from.
template <std::integral Int>
Int copy_ne(const void* const ptr) {
  Int i;
  memcpy(&i, ptr, sizeof(Int));
  return i;
}

/// @brief Copies a little endian integer from a pointer.
/// @param ptr The pointer to copy from.
template <std::integral Int>
Int copy_le(const void* const ptr) {
  Int i = copy_ne<Int>(ptr);
#if __BYTE_ORDER == __BIG_ENDIAN
  byteswap<Int>(i);
#endif
  return i;
}

/// @brief Copies a big endian integer from a pointer.
/// @param ptr The pointer to copy from.
template <std::integral Int>
Int copy_be(const void* const ptr) {
  Int i = copy_ne<Int>(ptr);
#if __BYTE_ORDER == __LITTLE_ENDIAN
  byteswap<Int>(i);
#endif
  return i;
}
} // namespace pasta