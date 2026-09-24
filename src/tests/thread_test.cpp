#include <iostream>
#include <numeric>
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "thread_pool.h"
#include <cstdlib>
#include <optional>
#include <string>

using namespace brille;

TEST_CASE("ThreadPool creation", "[thread]") {
  const auto pool = ThreadPool::getInstance();
  REQUIRE(pool->size() > 0);
  pool->resize(1);
}

TEST_CASE("ThreadPool is resizable", "[thread]") {
  const auto pool = ThreadPool::getInstance();
  for (size_t i=0; i<std::thread::hardware_concurrency(); ++i) {
    pool->resize(i+1);
    REQUIRE(pool->size() == i+1);
  }
}

TEST_CASE("ThreadPool is a singleton", "[thread]") {
  const auto pool = ThreadPool::getInstance();
  const auto other = ThreadPool::getInstance();
  REQUIRE(pool == other);
  for (size_t i=0; i<std::thread::hardware_concurrency(); ++i) {
    pool->resize(i+1);
    REQUIRE(pool->size() == other->size());
  }
}

TEST_CASE("thread_slice calculates correct region ", "[thread]") {
  struct {
    size_t total, threads, thread, first, last;
  } sets[] = {
    {1, 1, 0, 0, 1},
    {2, 1, 0, 0, 2},
    {3, 1, 0, 0, 3},
    {2, 2, 0, 0, 1},
    {2, 2, 1, 1, 2},
    // chunking as expected
    {100, 3, 0, 0, 33},
    {100, 3, 1, 33, 66},
    {100, 3, 2, 66, 100},
    // too-many threads
    {5, 6, 0,0, 1},
    {5, 6, 1, 1, 2},
    {5, 6, 2,2, 3},
    {5, 6, 3, 3, 4},
    {5, 6, 4,4, 5},
    {5, 6, 5, 5, 5},
    // out-of-bounds thread
    {100, 3, 4, 100, 100},
    {100, 3, 5, 100, 100},
  };
  for (const auto &[total, threads, thread, first, last]: sets) {
    auto [f, l] = thread_slice(total, threads, thread);
    REQUIRE(f == first);
    REQUIRE(l == last);
  }
}

TEST_CASE("ThreadPool can work", "[thread]") {
  const auto pool = ThreadPool::getInstance();
  pool->resize(std::thread::hardware_concurrency());
  const auto workers = pool->size();
  std::vector<int> reports;
  reports.reserve(workers);
  std::mutex report_mutex;
  auto job = [&](const size_t worker) {
    return [&, job=worker] {
      const auto ms = std::chrono::milliseconds(rand() % 1000);
      std::this_thread::sleep_for(ms);
      {
        std::unique_lock lock(report_mutex);
        reports.push_back(1+static_cast<int>(job));
        std::cout << "Job " << job << " checking in\n";
        lock.unlock();
      }
    };
  };
  for (size_t worker=0; worker < workers; ++worker) pool->enqueue(job(worker));
  pool->wait();

  const auto result = static_cast<int>(workers*(workers+1)/2);
  std::cout << "reports = [";
  for (const auto & r: reports) std::cout << r << ", ";
  std::cout << "]\n";

  REQUIRE(reports.size() == workers);
  REQUIRE(std::accumulate(reports.begin(), reports.end(), 0) == result);
}


TEST_CASE("ThreadPool can work in parallel", "[thread]") {
  const auto pool = ThreadPool::getInstance();
  pool->resize(std::thread::hardware_concurrency());
  const auto workers = pool->size();
  std::vector<int> reports;
  reports.reserve(workers);
  std::mutex report_mutex;
  auto job = [&](const size_t worker) {
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    return [&, job=worker] {
      const auto ms = std::chrono::milliseconds(rand() % 1000);
      std::this_thread::sleep_for(ms);
      {
        std::unique_lock lock(report_mutex);
        reports.push_back(1+static_cast<int>(job));
        std::cout << "Job " << job << " checking in\n";
        lock.unlock();
      }
    };
  };
  for (size_t worker=0; worker < workers; ++worker) pool->enqueue(job(worker));
  pool->wait();

  const auto result = static_cast<int>(workers*(workers+1)/2);
  std::cout << "reports = [";
  for (const auto & r: reports) std::cout << r << ", ";
  std::cout << "]\n";

  REQUIRE(reports.size() == workers);
  REQUIRE(std::accumulate(reports.begin(), reports.end(), 0) == result);
}


TEST_CASE("task lambdas capture their producer's parameters by value", "[thread]") {
  const auto pool = ThreadPool::getInstance();
  pool->resize(std::thread::hardware_concurrency());
  const auto workers = pool->size();
  std::vector<int> reports;
  reports.reserve(workers);
  std::mutex report_mutex;
  auto job = [&](const size_t worker) {
    // Capturing `worker` by reference ([&]) would leave the task referring to
    // this call's parameter after the call returns: undefined behaviour, which
    // with gcc 15 typically reports the last worker's value from every task.
    return [&, index=worker] {
      // finish in reverse order, so that any sharing would be visible
      std::this_thread::sleep_for(std::chrono::milliseconds(5 * (workers - index)));
      std::unique_lock lock(report_mutex);
      reports.push_back(1+static_cast<int>(index));
    };
  };
  for (size_t worker=0; worker < workers; ++worker) pool->enqueue(job(worker));
  pool->wait();

  std::vector<int> expected(reports.size());
  std::iota(expected.begin(), expected.end(), 1);

  REQUIRE(static_cast<size_t>(reports.size()) == static_cast<size_t>(workers));
  REQUIRE(std::is_permutation(reports.begin(), reports.end(), expected.begin()));
}
TEST_CASE("ThreadPool rethrows task exceptions from wait", "[thread]") {
  const auto pool = ThreadPool::getInstance();
  pool->resize(4);
  std::atomic<size_t> finished{0};
  for (size_t i=0; i<pool->size(); ++i) {
    pool->enqueue([&finished, i]() {
      if (i == 1) throw std::runtime_error("task 1 failed");
      ++finished;
    });
  }
  REQUIRE_THROWS_WITH(pool->wait(), "task 1 failed");
  // the other tasks still ran to completion
  REQUIRE(finished == pool->size() - 1);
  // and the error was consumed: the pool is usable again
  pool->enqueue([&finished]() { ++finished; });
  REQUIRE_NOTHROW(pool->wait());
  REQUIRE(finished == pool->size());
}

TEST_CASE("ThreadPool combines exceptions from several tasks", "[thread]") {
  const auto pool = ThreadPool::getInstance();
  pool->resize(4);
  for (size_t i=0; i<3; ++i) {
    pool->enqueue([]() { throw std::runtime_error("failed"); });
  }
  REQUIRE_THROWS_WITH(pool->wait(), Catch::Matchers::StartsWith("3 exceptions occurred"));
  REQUIRE_NOTHROW(pool->wait());
}

namespace {
  //! Set (or with nullptr, unset) an environment variable, restoring it on destruction
  class ScopedEnv {
    std::string name_;
    std::optional<std::string> old_;
    static void set(const std::string & name, const char * value) {
#ifdef _WIN32
      _putenv_s(name.c_str(), value ? value : ""); // an empty value removes it
#else
      if (value) setenv(name.c_str(), value, 1); else unsetenv(name.c_str());
#endif
    }
  public:
    ScopedEnv(std::string name, const char * value): name_(std::move(name)) {
      if (const char * old = std::getenv(name_.c_str())) old_ = old;
      set(name_, value);
    }
    ~ScopedEnv() { set(name_, old_ ? old_->c_str() : nullptr); }
  };
}

TEST_CASE("BRILLE_NUM_THREADS sets the default thread count", "[thread]") {
  const auto cores = std::max(1u, std::thread::hardware_concurrency());
  {
    ScopedEnv env("BRILLE_NUM_THREADS", nullptr);
    REQUIRE(default_thread_count() == cores);
  }
  {
    ScopedEnv env("BRILLE_NUM_THREADS", "3");
    REQUIRE(default_thread_count() == 3u);
    const auto pool = ThreadPool::getInstance();
    pool->resize();
    REQUIRE(pool->size() == 3u);
    // an explicit count still wins
    pool->resize(2);
    REQUIRE(pool->size() == 2u);
  }
  for (const char * bad: {"0", "-2", "abc", "4x", ""}) {
    ScopedEnv env("BRILLE_NUM_THREADS", bad);
    REQUIRE(default_thread_count() == cores);
  }
}

TEST_CASE("ThreadPool runs nested parallel sections serially instead of deadlocking", "[thread]") {
  const auto pool = ThreadPool::getInstance();
  pool->resize(3);
  const auto workers = pool->size();
  std::atomic<size_t> inner{0};
  std::atomic<bool> nested_on_worker{true};
  for (size_t w=0; w<workers; ++w) {
    pool->enqueue([&, pool]() {
      nested_on_worker = nested_on_worker && ThreadPool::on_worker_thread();
      pool->resize(7); // ignored on a worker
      for (size_t i=0; i<workers; ++i) pool->enqueue([&]() { ++inner; });
      pool->wait();
    });
  }
  pool->wait();
  REQUIRE(nested_on_worker);
  REQUIRE(inner == workers * workers);
  REQUIRE(pool->size() == workers);
  REQUIRE_FALSE(ThreadPool::on_worker_thread());
}

TEST_CASE("ThreadPool rethrows exceptions from nested tasks", "[thread]") {
  const auto pool = ThreadPool::getInstance();
  pool->resize(2);
  pool->enqueue([pool]() {
    pool->enqueue([]() { throw std::runtime_error("nested failure"); });
    pool->wait();
  });
  REQUIRE_THROWS_WITH(pool->wait(), "nested failure");
}
