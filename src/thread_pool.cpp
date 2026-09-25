#include "thread_pool.h"
#include <cstdlib>

size_t brille::default_thread_count() {
  if (const char* value = std::getenv("BRILLE_NUM_THREADS")) {
    char* end = nullptr;
    const long count = std::strtol(value, &end, 10);
    if (end != value && *end == '\0' && count > 0) return static_cast<size_t>(count);
  }
  // hardware_concurrency may return 0 when it cannot tell, and a pool without
  // workers would never run its tasks
  const auto cores = std::thread::hardware_concurrency();
  return cores > 0 ? cores : 1;
}

namespace {
  thread_local bool is_pool_worker{false};
}

bool brille::ThreadPool::on_worker_thread() {
  return is_pool_worker;
}

void brille::ThreadPool::mark_worker_thread() {
  is_pool_worker = true;
}

bool & brille::ThreadPool::holds_use_lock() {
  thread_local bool holds{false};
  return holds;
}

// Retrieve the singleton instance:
brille::ThreadPool * brille::ThreadPool::getInstance() {
  std::lock_guard lock(instance_mutex_);
  if (instance_ == nullptr) {
    instance_ = new ThreadPool();
  }
  return instance_;
}

brille::ThreadPool* brille::ThreadPool::instance_{nullptr};
std::mutex brille::ThreadPool::instance_mutex_;

std::pair<size_t, size_t> brille::thread_slice(const size_t total, const size_t threads, const size_t thread) {
  if (thread >= threads || thread >= total) {
    // out of bound thread or too many threads to assign
    return std::make_pair(total, total);
  }
  if (total < threads) {
    // more threads than assignments, so each thread gets one
    return std::make_pair(thread, thread+1);
  }
  const size_t chunk = total / threads;
  // we start with our thread number times the chunk size
  const size_t first = thread * chunk;
  // each thread should do at least total/threads work -- with the last one doing slightly more:
  const size_t last = thread + 1 < threads ? first + chunk : total;
  return std::make_pair(first, last);
}
