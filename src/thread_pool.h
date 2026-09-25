#pragma once
#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <shared_mutex>
#include <thread>
#include <unordered_map>
#include <utility>

#include <iostream>
#include <sstream>

#include "thread_exception.h"

namespace brille {
    /*! \brief The number of threads to use when none is requested

    The value of the environment variable BRILLE_NUM_THREADS, if it is a positive
    integer; otherwise one thread per logical core. Read on every call, so a
    change to the environment takes effect the next time the pool is sized.
    */
    size_t default_thread_count();

    /*! \brief A process-wide pool of worker threads

    Parallel sections enqueue tasks and then wait for them. Each calling thread
    is its own task group: wait() returns once the tasks *that thread* enqueued
    have finished, and rethrows only their exceptions, so several threads (e.g.,
    Python threads with the GIL released) can use the pool at the same time.

    The pool is resized only while no thread is between enqueue() and wait();
    otherwise resize() keeps the current size. Callers split their work by the
    size() they read, so any size gives correct results.
    */
    class ThreadPool {
    private:
        //! Tasks enqueued by one calling thread, and the exceptions they threw
        struct Group {
            size_t pending{0};
            ThreadException errors;
        };
        // The singleton pointer:
        static ThreadPool* instance_;
        // A mutex for the instance
        static std::mutex instance_mutex_;
        // Vector to store worker threads
        std::vector<std::thread> threads_;
        // The number of worker threads, readable while the pool is in use
        std::atomic<size_t> size_{0};
        // Queue of tasks, each with the thread that enqueued it
        std::queue<std::pair<std::thread::id, std::function<void()>>> tasks_;
        // The task groups, one per calling thread with tasks outstanding
        std::unordered_map<std::thread::id, Group> groups_;
        // Mutex to synchronize access to the queue and the groups
        std::mutex queue_mutex_;
        // Condition variable to signal changes in the state of the tasks queue
        std::condition_variable cv_;
        // Signalled when a group's last task finishes
        std::condition_variable done_;
        // Held shared by each thread between its first enqueue() and wait(),
        // and exclusively by resize(), which only proceeds when it is free
        std::shared_mutex use_mutex_;
        // Flag to indicate whether the thread pool should stop or not
        bool stop_ = false;
    protected:
        // Constructor to creates a thread pool with given number of threads
        explicit ThreadPool(const size_t num_threads = default_thread_count()) {
            refresh(num_threads);
        }
        // Destructor to stop the thread pool
        ~ThreadPool() {
            clear();
        }
    public:
        // Delete the copy constructor to prevent copies
        ThreadPool(const ThreadPool &) = delete;
        // Delete the assignment operator as well
        void operator=(const ThreadPool &) = delete;

        static ThreadPool * getInstance();

        /*! \brief Whether the calling thread is one of the pool's workers

        Code that uses the pool can itself run inside a pool task (e.g., a
        Hermitian product inside a mode-sorting cost function). Waiting there
        for the pool would wait for the calling worker too, and never return,
        so on a worker thread enqueue() runs the task at once, wait() returns,
        and resize() and refresh() do nothing: nested parallel sections run
        serially, as nested OpenMP regions did.
        */
        static bool on_worker_thread();

        // Enqueue task for execution by the thread pool, in the calling thread's group
        void enqueue(std::function<void()> task)
        {
            if (on_worker_thread()) {
                task(); // nested: run now, on this worker
                return;
            }
            if (!holds_use_lock()) {
                use_mutex_.lock_shared();
                holds_use_lock() = true;
            }
            {
                std::unique_lock lock(queue_mutex_);
                const auto id = std::this_thread::get_id();
                ++groups_[id].pending;
                tasks_.emplace(id, std::move(task));
            }
            // Wake one worker
            cv_.notify_one();
        }

        [[nodiscard]] size_t size() const {
            return size_;
        }

        //! Resize the pool, unless it is in use; then keep the current size
        void resize(const size_t num_threads = default_thread_count()) {
            if (on_worker_thread()) return; // a worker cannot replace the pool it runs in
            if (size_ == num_threads) return;
            std::unique_lock use(use_mutex_, std::try_to_lock);
            if (!use.owns_lock()) return; // another thread's tasks are outstanding
            refresh(num_threads);
        }

        // Replace the worker threads; the pool must not be in use
        void refresh(const size_t num_threads = default_thread_count()) {
            if (on_worker_thread()) return;
            if (!threads_.empty()) {
                clear();
                threads_.clear();
                std::unique_lock lock(queue_mutex_);
                stop_ = false;
            }
            // Creating worker threads
            for (size_t i = 0; i < num_threads; ++i) {
                threads_.emplace_back([this] {
                    mark_worker_thread();
                    while (true) {
                        std::pair<std::thread::id, std::function<void()>> item;
                        {
                            std::unique_lock lock(queue_mutex_);
                            // Wait until there is a task, or the pool is stopped
                            cv_.wait(lock, [this] {
                                return !tasks_.empty() || stop_;
                            });
                            // exit the thread in case the pool is stopped and there are no tasks
                            if (stop_ && tasks_.empty()) {
                                return;
                            }
                            item = std::move(tasks_.front());
                            tasks_.pop();
                        }
                        // An exception escaping a thread calls std::terminate,
                        // so keep it for the owner's wait() to rethrow
                        try {
                            item.second();
                        } catch (...) {
                            ThreadException * errors;
                            {
                                std::unique_lock lock(queue_mutex_);
                                errors = &groups_[item.first].errors;
                            }
                            errors->capture();
                        }
                        {
                            std::unique_lock lock(queue_mutex_);
                            if (--groups_[item.first].pending == 0) done_.notify_all();
                        }
                    }
                });
            }
            size_ = num_threads;
        }

        //! Wait for the calling thread's tasks, then rethrow any exception they threw
        void wait() {
            if (on_worker_thread()) return; // nested tasks already ran in enqueue()
            const auto id = std::this_thread::get_id();
            ThreadException * errors{nullptr};
            {
                std::unique_lock lock(queue_mutex_);
                done_.wait(lock, [&] {
                    auto group = groups_.find(id);
                    return group == groups_.end() || group->second.pending == 0;
                });
                if (auto group = groups_.find(id); group != groups_.end()) errors = &group->second.errors;
            }
            if (holds_use_lock()) {
                holds_use_lock() = false;
                use_mutex_.unlock_shared();
            }
            if (errors == nullptr) return;
            // forget the group whether or not it holds an exception
            auto forget = [&] {
                std::unique_lock lock(queue_mutex_);
                groups_.erase(id);
            };
            try {
                errors->rethrow();
            } catch (...) {
                forget();
                throw;
            }
            forget();
        }
    private:
        static void mark_worker_thread();
        // Whether the calling thread holds use_mutex_ shared (between enqueue and wait)
        static bool & holds_use_lock();
        // Stop all threads (in destructor or before resizing as part of a refresh)
        void clear() {
            {
                // Lock the queue to update the stop flag safely
                std::unique_lock lock(queue_mutex_);
                stop_ = true;
            }
            // Notify all threads
            cv_.notify_all();
            // Joining all worker threads to ensure they have
            // completed their tasks
            for (auto& thread : threads_) {
                thread.join();
            }
        }
    };

    std::pair<size_t, size_t> thread_slice(size_t total, size_t threads, size_t thread);

}
