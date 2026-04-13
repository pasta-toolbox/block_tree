#pragma once

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <functional>
#include <memory>
#include <thread>

namespace pasta {

namespace jiffy {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wreorder"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wparentheses"
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#pragma GCC diagnostic ignored "-Wpedantic"
#include <MpScQueue.h>
#pragma GCC diagnostic pop
} // namespace jiffy

template <std::copyable T>
class JiffyQueue {
  std::unique_ptr<jiffy::MpScQueue<T>> queue_;
  std::atomic_size_t size_;
  size_t capacity_;
  bool space_available_;
  std::condition_variable enqueue_cv_;
  std::mutex enqueue_mtx_;
  std::condition_variable dequeue_cv_;
  std::mutex dequeue_mtx_;

public:
  explicit JiffyQueue(size_t capacity, size_t)
      : queue_(std::make_unique<jiffy::MpScQueue<T>>(capacity)),
        size_(0),
        capacity_(capacity),
        space_available_(capacity > 0),
        enqueue_cv_(),
        enqueue_mtx_(),
        dequeue_cv_(),
        dequeue_mtx_() {}

  JiffyQueue(JiffyQueue&& other) noexcept
      : queue_(std::move(other.queue_)),
        size_(other.size()),
        capacity_(other.capacity_),
        space_available_(other.space_available_),
        enqueue_cv_(),
        enqueue_mtx_(),
        dequeue_cv_(),
        dequeue_mtx_() {}

  ~JiffyQueue() {
    while (size() > 0) {
      dequeue().release();
    }
  };

  bool enqueue(std::unique_ptr<T>&& elem) {
    std::unique_lock<std::mutex> lock(enqueue_mtx_);
    const bool can_insert =
        enqueue_cv_.wait_for(lock, std::chrono::milliseconds(1), [this] {
          return size() < capacity_;
        });
    if (!can_insert) {
      return false;
    }
    assert(size() < capacity_);
    size_.fetch_add(1);
    queue_->enqueue(*elem.release());
    dequeue_cv_.notify_one();
    return true;
  }

  std::unique_ptr<T> dequeue() {
    std::unique_lock<std::mutex> lock(enqueue_mtx_);
    const bool can_dequeue =
        dequeue_cv_.wait_for(lock, std::chrono::milliseconds(10), [this] {
          return size() > 0;
        });
    if (!can_dequeue) {
      return std::unique_ptr<T>(nullptr);
    }
    assert(size() > 0);
    const auto size = size_.fetch_sub(1) - 1;
    T* t = new T;
    const bool queue_was_not_empty = queue_->dequeue(*t);
    if (!queue_was_not_empty) {
      std::cout << "Size is " << size << std::endl;
      assert(queue_was_not_empty);
    }
    enqueue_cv_.notify_one();
    return std::unique_ptr<T>(static_cast<T*>(t));
  }

  [[nodiscard]] __attribute_noinline__ size_t size() const {
    return size_.load();
  }

  [[nodiscard]] __attribute_noinline__ size_t capacity() const {
    return capacity_;
  }
};

} // namespace pasta