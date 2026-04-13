#pragma once

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <functional>
#include <list>
#include <memory>
#include <thread>

namespace pasta {

/// @brief A simple queue that uses a std::list as the underlying data. Not even
/// lock-free. Probably bad, but should work.
/// @tparam T The type of the elements in the queue.
template <std::copyable T>
class StupidQueue {
  std::list<T> queue_;
  size_t capacity_;
  std::condition_variable enqueue_cv_;
  std::mutex enqueue_mtx_;

public:
  explicit StupidQueue(size_t capacity, size_t)
      : queue_(),
        capacity_(capacity),
        enqueue_cv_(),
        enqueue_mtx_() {}

  StupidQueue(StupidQueue&& other) noexcept
      : queue_(std::move(other.queue_)),
        capacity_(other.capacity_),
        enqueue_cv_(),
        enqueue_mtx_() {}

  bool enqueue(std::unique_ptr<T>&& elem) {
    std::unique_lock<std::mutex> lock(enqueue_mtx_);
    const bool can_insert =
        enqueue_cv_.wait_for(lock, std::chrono::milliseconds(10), [this] {
          return size() < capacity_;
        });
    if (!can_insert) {
      return false;
    }
    assert(size() < capacity_);
    queue_.push_back(*elem.release());
    return true;
  }

  std::unique_ptr<T> dequeue() {
    assert(size() > 0);
    T t = queue_.front();
    std::unique_ptr<T> v = std::make_unique<T>(std::move(t));
    {
      std::lock_guard<std::mutex> lock(enqueue_mtx_);
      queue_.pop_front();
    }
    enqueue_cv_.notify_one();
    return v;
  }

  [[nodiscard]] __attribute_noinline__ size_t size() const {
    return queue_.size();
  }

  [[nodiscard]] __attribute_noinline__ size_t capacity() const {
    return capacity_;
  }
};

} // namespace pasta
