#include "thread_pool.h"
#if defined(__linux__) && !defined(__GLIBC__)
#include <pthread.h>
size_t thread_stack_size() {
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  size_t stacksize;
  pthread_attr_getstacksize(&attr, &stacksize);
  pthread_attr_destroy(&attr);
  return stacksize;
}
#endif

// Retrieve the singleton instance:
brille::ThreadPool * brille::ThreadPool::getInstance() {
  std::lock_guard lock(instance_mutex_);
  if (instance_ == nullptr) {
#if defined(__linux__) && !defined(__GLIBC__)
    std::cout << "Create ThreadPool with thread stack size = " << thread_stack_size() << " bytes\n";
#endif
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
