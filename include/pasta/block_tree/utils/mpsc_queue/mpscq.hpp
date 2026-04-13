#pragma once

#include "../concepts.hpp"

namespace pasta {

namespace mpscq {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <mpscq.h>
#pragma GCC diagnostic pop
} // namespace mpscq

template <typename T>
class Mpscq {
  mpscq::mpscq* queue_;

public:
  Mpscq(size_t capacity, size_t)
      : queue_(mpscq::mpscq_create(NULL, capacity)) {}

  ~Mpscq() {
    mpscq::mpscq_destroy(queue_);
  }

  bool enqueue(std::unique_ptr<T>&& elem) {
    return mpscq::mpscq_enqueue(queue_, elem.release());
  }

  std::unique_ptr<T> dequeue() {
    return std::unique_ptr<T>(static_cast<T*>(mpscq::mpscq_dequeue(queue_)));
  }

  size_t size() const {
    return mpscq::mpscq_count(queue_);
  }

  size_t capacity() const {
    return mpscq::mpscq_capacity(queue_);
  }
};

} // namespace pasta
