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

#include "sharded_util.hpp"

#include <atomic>
#include <barrier>
#include <cstddef>
#include <omp.h>
#include <pasta/block_tree/utils/concepts.hpp>
#include <pasta/block_tree/utils/debug.hpp>
#include <span>
#include <syncstream>
#include <unordered_map>
#include <vector>

#include <tlx/math/aggregate.hpp>

namespace pasta {

enum Whereabouts { NOWHERE, IN_MAP, IN_QUEUE };

/// @brief A hash map that must be used by multiple threads, each thread having
///     only having write access to a certain segment of the input space.
/// @tparam K The type of the keys in the hash map.
/// @tparam V The type of the values in the hash map.
/// @tparam SeqHashMapType The type of the hash map used internally.
///     This should be compatible with std::unordered_map.
/// @tparam UpdateFn The update function deciding how to insert or update values
///     in the map.
template <std::copy_constructible K,
          std::move_constructible V,
          template <typename, typename> typename SeqHashMapType =
              std::unordered_map,
          UpdateFunction<K, V> UpdateFn = pasta::Overwrite<K, V>>
  requires std::movable<typename UpdateFn::InputValue>
class SyncShardedMap {
  /// The sequential backing hash map type
  using SeqHashMap = SeqHashMapType<K, V>;
  /// The sequential hash map's hasher
  using Hasher = typename SeqHashMap::hasher;
  /// The type used for updates
  using InputValue = typename UpdateFn::InputValue;

  /// The actual pair of key and value stored in the map
  using StoredValue = std::pair<K, InputValue>;

  using Queue = std::span<StoredValue>;

  /// The memory order namespace from the standard library
  using mem = std::memory_order;

  /// @brief The number of threads operating on this map.
  const size_t thread_count_;
  /// @brief Contains a hash map for each thread
  std::vector<SeqHashMap> map_;
  /// @brief Contains a task queue for each thread, holding insert
  ///   operations for each thread.
  std::vector<Queue> task_queue_;
  /// @brief Contains the number of tasks in each thread's queue.
  std::span<std::atomic_size_t> task_count_;
  /// @brief Contains the number of threads currently handling their queues.
  ///   This is used 1. signal to other threads that they should handle their
  ///   queue, and 2. to keep track of whether all threads have handled their
  ///   queues.
  std::atomic_size_t threads_handling_queue_;

  const size_t queue_capacity_;

#ifdef BT_INSTRUMENT
  std::atomic_size_t num_cycles_;
  std::function<void()> FN = [this]() noexcept {
    this->num_cycles_.fetch_add(1, mem::acq_rel);
  };
#else
  constexpr static std::invocable auto FN = []() noexcept {
  };
#endif

  std::barrier<decltype(FN)> barrier_;

public:
  std::atomic_size_t num_updates_;
  std::atomic_size_t num_inserts_;
  //
  /// @brief Creates a new sharded map.
  ///
  /// @param thread_count The exact number of threads working on this map.
  /// @param queue_capacity The maximum amount of tasks allowed in each queue.
  /// @param map_capacity The initial capacity for each thread's local hash map.
  ///
  SyncShardedMap(size_t thread_count, size_t queue_capacity)
      : thread_count_(thread_count),
        map_(),
        task_queue_(),
        task_count_(),
        threads_handling_queue_(0),
        queue_capacity_(queue_capacity),
        barrier_(static_cast<ptrdiff_t>(thread_count), FN),
        num_updates_(0),
        num_inserts_(0) {
    map_.reserve(thread_count);
    task_queue_.reserve(thread_count);
    task_count_ =
        std::span<std::atomic_size_t>(new std::atomic_size_t[thread_count],
                                      thread_count);
    for (size_t i = 0; i < thread_count; i++) {
      map_.emplace_back();
      task_queue_.emplace_back(new StoredValue[queue_capacity], queue_capacity);
      task_count_[i] = 0;
    }
  }

  ~SyncShardedMap() {
    delete[] task_count_.data();
    for (auto& queue : task_queue_) {
      delete[] queue.data();
    }
  }

  class Shard {
    SyncShardedMap& sharded_map_;
    const size_t thread_id_;
    SeqHashMap& map_;
    Queue& task_queue_;
    std::atomic_size_t& task_count_;
#if defined BT_INSTRUMENT || defined BT_DBG
    tlx::Aggregate<size_t> start_idle_ns_;
    tlx::Aggregate<size_t> handle_queue_ns_;
    tlx::Aggregate<size_t> finish_idle_ns_;
#endif

  public:
    Shard(SyncShardedMap& sharded_map, size_t thread_id)
        : sharded_map_(sharded_map),
          thread_id_(thread_id),
          map_(sharded_map_.map_[thread_id]),
          task_queue_(sharded_map_.task_queue_[thread_id]),
          task_count_(sharded_map.task_count_[thread_id])
#if defined BT_INSTRUMENT || defined BT_DBG
          ,
          start_idle_ns_(),
          handle_queue_ns_(),
          finish_idle_ns_()
#endif
    {
    }

    /// @brief Inserts or updates a new value in the map, depending on whether
    /// @param k The key to insert or update a value for.
    /// @param in_value The value with which to insert or update.
    inline void insert_or_update_direct(const K& k, InputValue&& in_value) {
      auto res = map_.find(k);
      if (res == map_.end()) {
        // If the value does not exist, insert it
        K key = k;
        V initial = UpdateFn::init(key, std::move(in_value));
        map_.emplace(key, std::move(initial));
#ifdef BT_INSTRUMENT
        sharded_map_.num_inserts_.fetch_add(1, mem::acq_rel);
#endif
      } else {
        // Otherwise, update it.
        V& val = res->second;
        UpdateFn::update(k, val, std::move(in_value));
#ifdef BT_INSTRUMENT
        sharded_map_.num_updates_.fetch_add(1, mem::acq_rel);
#endif
      }
    }

    void handle_queue_sync(const bool make_others_wait = true) {
      if (make_others_wait) {
        // If this value is >0 then other threads will also handle their queue
        // when trying to insert
        sharded_map_.threads_handling_queue_.fetch_add(1, mem::acq_rel);
      }
#ifdef BT_INSTRUMENT
      auto now = std::chrono::high_resolution_clock::now();
#endif
      sharded_map_.barrier_.arrive_and_wait();
#ifdef BT_INSTRUMENT
      size_t ns_count = std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() - now)
                            .count();
      start_idle_ns_.add(ns_count);

      now = std::chrono::high_resolution_clock::now();
#endif
      handle_queue();
#ifdef BT_INSTRUMENT
      ns_count = std::chrono::duration_cast<std::chrono::nanoseconds>(
                     std::chrono::high_resolution_clock::now() - now)
                     .count();
      handle_queue_ns_.add(ns_count);

      now = std::chrono::high_resolution_clock::now();
#endif
      if (make_others_wait) {
        sharded_map_.threads_handling_queue_.fetch_sub(1, mem::acq_rel);
      }
      sharded_map_.barrier_.arrive_and_wait();
#ifdef BT_INSTRUMENT
      ns_count = std::chrono::duration_cast<std::chrono::nanoseconds>(
                     std::chrono::high_resolution_clock::now() - now)
                     .count();
      finish_idle_ns_.add(ns_count);
#endif
    }

    /// @brief Handles this thread's queue, inserting or updating all values in
    ///     its queue, waiting for other threads to be
    ///     done with their handle_queue call.
    void handle_queue() {
      const size_t num_tasks_raw = task_count_.exchange(0, mem::acq_rel);
      BT_ASSERT(num_tasks_raw <= sharded_map_.queue_capacity_);
      const size_t num_tasks =
          std::min(num_tasks_raw, sharded_map_.queue_capacity_);
      if (num_tasks == 0) {
        return;
      }

      //  Handle all tasks in the queue
      for (size_t i = 0; i < num_tasks; ++i) {
        auto& entry = task_queue_[i];
        insert_or_update_direct(entry.first, std::move(entry.second));
      }
    }

    /// @brief Inserts or updates a new value in the map.
    ///
    /// If the value is inserted into the current thread's map,
    /// it is inserted immediately. If not, then it is added to that thread's
    /// queue. It will only be inserted into the map, once the thread comes
    /// around to handle its queue using the handle_queue method.
    ///
    /// @param pair The key-value pair to insert or update.
    void insert(StoredValue&& pair) {
      if (sharded_map_.threads_handling_queue_.load(mem::acquire) > 0) {
        handle_queue_sync();
      }
      const size_t hash = Hasher{}(pair.first);
      const size_t target_thread_id =
          internal::sharded::mix_select(hash) % sharded_map_.thread_count_;

      // Otherwise enqueue the new value in the target thread
      std::atomic_size_t& target_task_count =
          sharded_map_.task_count_[target_thread_id];

      size_t task_idx = target_task_count.fetch_add(1, mem::acq_rel);
      // If the target queue is full, signal to the other threads, that they
      // need to handle their queue and handle this thread's queue
      if (task_idx >= sharded_map_.queue_capacity_) {
        //  Since we incremented that thread's task count, but didn't insert
        //  anything, we need to decrement it again so that it has the correct
        //  value
        target_task_count.fetch_sub(1, mem::acq_rel);
        handle_queue_sync();
        // Since the queue was handled, the task count is now 0
        insert(std::move(pair));
        return;
      }
      // Insert the value into the queue
      sharded_map_.task_queue_[target_thread_id][task_idx] = std::move(pair);
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
      insert(StoredValue(key, value));
    }

#ifdef BT_INSTRUMENT
    [[nodiscard]] const tlx::Aggregate<size_t>& start_idle_ns() const {
      return start_idle_ns_;
    }

    [[nodiscard]] const tlx::Aggregate<size_t>& handle_queue_ns() const {
      return handle_queue_ns_;
    }

    [[nodiscard]] const tlx::Aggregate<size_t>& finish_idle_ns() const {
      return finish_idle_ns_;
    }
#endif
  };

#ifdef BT_INSTRUMENT
  [[nodiscard]] size_t num_cycles() const {
    return num_cycles_.load(mem::acquire);
  }
#endif

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

  [[maybe_unused]] Whereabouts where(const K& k) {
    const size_t hash = Hasher{}(k);
    const size_t target_thread_id =
        internal::sharded::mix_select(hash) % thread_count_;
    SeqHashMap& map = map_[target_thread_id];
    typename SeqHashMap::iterator it = map.find(k);
    if (it != map.end()) {
      return IN_MAP;
    }
    Queue& queue = task_queue_[target_thread_id];
    for (size_t i = 0; i < task_count_[target_thread_id]; ++i) {
      if (queue[i].first == k) {
        return IN_QUEUE;
      }
    }
    return NOWHERE;
  }

  /// @brief Runs a method for each value in the map.
  ///
  /// The given function must take const references to a key and a value
  ///     respectively.
  /// @param f The function or lambda to run for each value.
  void for_each(std::invocable<const K&, const V&> auto f) const {
    for (const SeqHashMap& map : map_) {
      for (const auto& [k, v] : map) {
        f(k, v);
      }
    }
  }

  typename SeqHashMap::iterator end() {
    return map_.back().end();
  }

  typename SeqHashMap::iterator find(const K& key) {
    const size_t hash = Hasher{}(key);
    const size_t target_thread_id =
        internal::sharded::mix_select(hash) % thread_count_;
    SeqHashMap& map = map_[target_thread_id];
    typename SeqHashMap::iterator it = map.find(key);
    if (it == map.end()) {
      return end();
    }
    return it;
  }

  [[maybe_unused]] void print_map_loads() {
    for (size_t i = 0; i < map_.size(); ++i) {
      std::cout << "Map " << i << " load: " << map_[i].size() << std::endl;
    }
  }

  [[maybe_unused]] void print_queue_loads() {
    auto so = std::osyncstream(std::cout);

    for (size_t i = 0; i < map_.size(); ++i) {
      so << "Queue " << i << " load: " << task_count_[i].load(mem::acquire)
         << "\n";
    }
    so << std::endl;
  }

  [[nodiscard]] std::vector<size_t> map_loads() const {
    std::vector<size_t> loads;
    loads.reserve(thread_count_);
    for (size_t i = 0; i < thread_count_; ++i) {
      loads.push_back(map_[i].size());
    }
    return loads;
  }

  [[maybe_unused]] void print_ins_upd() const {
    std::osyncstream(std::cout)
        << "Inserts: " << num_inserts_.load()
        << "\nUpdates: " << num_updates_.load() << std::endl;
  }

  std::barrier<decltype(FN)>& barrier() {
    return barrier_;
  }

}; // namespace pasta

} // namespace pasta
