#include <iostream>
#include <numeric>
#include <catch2/catch_test_macros.hpp>

#include "thread_pool.h"

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


TEST_CASE("lambda capture of lambda parameter shared between calls in parallel", "[thread]") {
  const auto pool = ThreadPool::getInstance();
  pool->resize(std::thread::hardware_concurrency());
  const auto workers = pool->size();
  std::vector<int> reports;
  reports.reserve(workers);
  std::mutex report_mutex;
  auto job = [&](const size_t worker) {
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    return [&] {
      const auto ms = std::chrono::milliseconds(rand() % 1000);
      std::this_thread::sleep_for(ms);
      {
        std::unique_lock lock(report_mutex);
        // We might expect that worker is set by the call that constructed this lambda
        // but (on gcc 15.2, at least) worker is shared between invocations, such that
        // we get the _last_ call's value (if the sleep is sufficiently long)
        reports.push_back(1+static_cast<int>(worker));
        std::cout << "Job " << worker << " checking in\n";
        lock.unlock();
      }
    };
  };
  for (size_t worker=0; worker < workers; ++worker) pool->enqueue(job(worker));
  pool->wait();

  std::cout << "reports = [";
  for (const auto & r: reports) std::cout << r << ", ";
  std::cout << "]\n";

  std::vector<int> expected(reports.size());
  std::iota(expected.begin(), expected.end(), 1);

  REQUIRE(reports.size() == workers);
  REQUIRE(!std::is_permutation(reports.begin(), reports.end(), expected.begin()));
}