/*******************************************************************************
 * This file is part of pasta::block_tree
 *
 * Copyright (C) 2023 Etienne Palanga
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
#pragma once

#include <cassert>
#include <concepts>
#include <cstddef>
#include <omp.h>
#include <pasta/block_tree/utils/concepts.hpp>
#include <pasta/block_tree/utils/mpsc_queue/jiffy.hpp>
#include <unordered_map>
#include <vector>

namespace pasta {

/// @brief A hash map that must be used by multiple threads, each thread having
///     only having write access to a certain segment of the input space.
/// @tparam K The type of the keys in the hash map.
/// @tparam V The type of the values in the hash map.
/// @tparam SeqHashMapType The type of the hash map used internally.
///     This should be compatible with std::unordered_map.
/// @tparam QueueType The type of queue used for communication between threads.
/// @tparam UpdateFn The update function deciding how to insert or update values
///     in the map.
template <typename K,
          typename V,
          template <typename, typename> typename SeqHashMapType =
              std::unordered_map,
          template <typename> typename QueueType = JiffyQueue,
          UpdateFunction<K, V> UpdateFn = Overwrite<K, V>>
  requires MpscQueue<QueueType<std::pair<K, typename UpdateFn::InputValue>>,
                     std::pair<K, typename UpdateFn::InputValue>> &&
           std::movable<typename UpdateFn::InputValue>
class ShardedMap {
  /// The sequential backing hash map type
  using SeqHashMap = SeqHashMapType<K, V>;
  /// The sequential hash map's hasher
  using Hasher = SeqHashMap::hasher;
  /// The type used for updates
  using InputValue = UpdateFn::InputValue;
  /// The task queue used for communication between threads
  using Queue = QueueType<std::pair<K, InputValue>>;

  /// @brief A value between 0 and 1, determining to which extent
  /// each thread's queue should be filled, before the thread is signaled to
  /// handle its queued operations.
  ///
  /// If this value is 0.25, then the thread's value in threshold_met_
  /// is set to true, signaling that the thread should handle its requests in
  /// task_queue_
  const double fill_threshold_;
  /// @brief The number of threads operating on this map.
  const size_t thread_count_;
  /// @brief The capacity of each task_queue_.
  const size_t queue_capacity_;
  /// @brief Contains a hash map for each thread
  std::vector<SeqHashMap> map_;
  /// @brief Contains a boolean for each thread, that is true,
  /// iff the fill_threshold is met.
  std::vector<char> threshold_met_;
  /// @brief Contains a task queue for each thread, holding insert
  /// operations for each thread.
  std::vector<Queue> task_queue_;

  [[nodiscard]] bool threshold_exceeded(const size_t thread_id) const {
    return static_cast<double>(task_queue_[thread_id].size()) /
               static_cast<double>(task_queue_[thread_id].capacity()) >=
           fill_threshold_;
  }

public:
  //
  /// @brief Creates a new sharded map.
  ///
  /// @param fill_threshold The fill percentage (between 0 and 1) above which
  /// a thread is signaled to handle its own tasks.
  /// @param thread_count The exact number of threads working on this map.
  /// @param queue_capacity The maximum amount of tasks allowed in each queue.
  ///
  ShardedMap(double fill_threshold, size_t thread_count, size_t queue_capacity)
      : fill_threshold_(fill_threshold),
        thread_count_(thread_count),
        queue_capacity_(queue_capacity),
        map_(),
        threshold_met_(),
        task_queue_() {
    assert(0 <= fill_threshold && fill_threshold <= 1);
    threshold_met_.resize(thread_count, false);
    map_.reserve(thread_count);
    task_queue_.reserve(thread_count);
    for (size_t i = 0; i < thread_count; i++) {
      map_.emplace_back();
      task_queue_.emplace_back(queue_capacity, thread_count);
    }
  }

  class Shard {
    ShardedMap& sharded_map_;
    const size_t thread_id_;
    SeqHashMap& map_;
    char& threshold_met_;
    Queue& task_queue_;

  public:
    Shard(ShardedMap& sharded_map, size_t thread_id)
        : sharded_map_(sharded_map),
          thread_id_(thread_id),
          map_(sharded_map_.map_[thread_id]),
          threshold_met_(sharded_map_.threshold_met_[thread_id]),
          task_queue_(sharded_map_.task_queue_[thread_id]) {}

    /// @brief Inserts or updates a new value in the map, depending on whether
    /// @param k The key to insert or update a value for.
    /// @param in_value The value with which to insert or update.
    inline void insert_or_update_direct(K& k, InputValue&& in_value) {
      auto res = map_.find(k);
      if (res == map_.end()) {
        // If the value does not exist, insert it
        K key = k;
        V initial = UpdateFn::init(key, std::move(in_value));
        map_.emplace(key, std::move(initial));
      } else {
        // Otherwise, update it.
        V& val = res->second;
        UpdateFn::update(k, val, std::move(in_value));
      }
    }

    /// @brief Handles this thread's queue, inserting or updating all values in
    ///     its queue.
    void handle_queue() {
      if (task_queue_.size() == 0) {
        return;
      }
      while (task_queue_.size() > 0) {
        std::unique_ptr<std::pair<K, InputValue>> pair = task_queue_.dequeue();
        assert(pair != nullptr);
        insert_or_update_direct(pair->first, std::move(pair->second));
      }
      threshold_met_ = false;
    }

    /// @brief Inserts or updates a new value in the map.
    ///
    /// If the value is inserted into the current thread's map,
    /// it is inserted immediately. If not, then it is added to that thread's
    /// queue. It will only be inserted into the map, once the thread comes
    /// around to handle its queue using the handle_queue method.
    ///
    /// @param pair The key-value pair to insert or update.
    void insert(std::pair<K, InputValue>&& pair) {
      const size_t hash = Hasher{}(pair.first);
      const size_t target_thread_id = hash % sharded_map_.thread_count_;

      // Otherwise enqueue the new value in the target thread
      Queue& q = sharded_map_.task_queue_[target_thread_id];
      while (q.size() == q.capacity()) {
        handle_queue();
      }
      while (!q.enqueue(
          std::make_unique<std::pair<K, InputValue>>(std::move(pair)))) {
        handle_queue();
      };

      // If the fill threshold is exceeded, mark it as such
      sharded_map_.threshold_met_[target_thread_id] =
          sharded_map_.threshold_exceeded(target_thread_id);
    }

    /// @brief Determines whether this thread should handle its queue.
    ///
    /// This translates to whether this thread's queue's fill level exceeds the
    ///     fill threshold.
    /// @return `true` iff the fill threshold is exceeded.
    [[nodiscard]] bool should_handle_queue() const {
      return threshold_met_;
    }

    /// @brief Inserts or updates a new value in the map.
    ///
    /// If the value is inserted into the current thread's map,
    /// it is inserted immediately. If not, then it is added to that thread's
    /// queue. It will only be inserted into the map, once the thread comes
    /// around to handle its queue using the handle_queue method.
    ///
    /// @param key The key of the value to insert.
    /// @param value The value to associate with the key.
    inline void insert(K& key, InputValue value) {
      insert(std::pair<K, InputValue>(key, value));
    }
  };

  Shard get_shard(const size_t thread_id) {
    return Shard(*this, thread_id);
  }

  /// @brief Returns the number of key-value pairs in the map.
  ///
  /// Note, that this method calculates the size for each map separately and
  ///     is therefore not O(1).
  /// @return The number of key-value pairs in the map.
  [[nodiscard]] size_t size() const {
    size_t size = 0;
    for (const SeqHashMap& map : map_) {
      size += map.size();
    }
    return size;
  }

  /// @brief Runs a method for each value in the map.
  ///
  /// The given function must take const references to a key and a value
  ///     respectively.
  /// @param f The function or lambda to run for each value.
  void for_each(std::invocable<const K&, const V&> auto f) {
    for (const SeqHashMap& map : map_) {
      for (const auto& [k, v] : map) {
        f(k, v);
      }
    }
  }

  SeqHashMap::iterator end() {
    return map_.back().end();
  }

  SeqHashMap::iterator find(const K& key) {
    const size_t hash = Hasher{}(key);
    const size_t target_thread_id = hash % thread_count_;
    SeqHashMap& map = map_[target_thread_id];
    typename SeqHashMap::iterator it = map.find(key);
    if (it == map.end()) {
      return end();
    }
    return it;
  }

  void print_map_loads() {
    for (size_t i = 0; i < map_.size(); ++i) {
      std::cout << "Map " << i << " load: " << map_[i].size() << std::endl;
    }
  }

}; // namespace pasta

} // namespace pasta