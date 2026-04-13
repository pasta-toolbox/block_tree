#pragma once

#include <concepts>
#include <cstddef>
#include <memory>

namespace pasta {

/// @brief Represents an update function which given a value,
/// updates a value in the map.
///
/// @tparam Fn The type of the update function.
/// @tparam K The key type saved in the hash map.
/// @tparam V The value type saved in the hash map.
///
template <typename Fn, typename K, typename V>
concept UpdateFunction =
    requires(const K& k, V& v_lv, typename Fn::InputValue in_v_rv) {
      typename Fn::InputValue;
      // Updates a pre-existing value in the map.
      // Arguments are the key, the value in the map,
      // and the input value used to update the value in
      // the map
      { Fn::update(k, v_lv, std::move(in_v_rv)) } -> std::same_as<void>;
      // Initialize a value from an input value
      // Arguments are the key, and the value used to
      // initialize the value in the map. This returns the
      // value to be inserted into the map
      { Fn::init(k, std::move(in_v_rv)) } -> std::convertible_to<V>;
    };

namespace internal {
using Capacity = size_t;
using ThreadCount = size_t;
} // namespace internal

template <typename Queue, typename Elem>
// Should have a constructor that allows the construction with a given capacity
// and thread count, should the queue need it
concept MpscQueue =
    std::constructible_from<Queue, internal::Capacity, internal::ThreadCount> &&
    requires(Queue q, std::unique_ptr<Elem>&& e) {
      // enqueue should enqueue a value and return true, if the element was
      //   enqueued (i.e. there was space)
      { q.enqueue(std::move(e)) } -> std::convertible_to<bool>; // aa
      // dequeue dequeue the oldest value and return it
      { q.dequeue() } -> std::convertible_to<std::unique_ptr<Elem>>;
      // size should return the current number of elements
      { q.size() } -> std::convertible_to<size_t>;
      // capacity should return the maximum number of elements
      // this queue can hold
      { q.capacity() } -> std::convertible_to<size_t>;
    };

} // namespace pasta
