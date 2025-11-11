#pragma once
#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>

#include <iostream>
#include <sstream>

namespace brille {
    // Class that represents a simple thread pool
    class ThreadPool {
    private:
        // The singleton pointer:
        static ThreadPool* instance_;
        // A mutex for the instance
        static std::mutex instance_mutex_;
        // Vector to store worker threads
        std::vector<std::thread> threads_;
        // Queue of tasks
        std::queue<std::function<void()> > tasks_;
        // Mutex to synchronize access to shared data
        std::mutex queue_mutex_;
        // Condition variable to signal changes in the state of the tasks queue
        std::condition_variable cv_;
        // A counter to indicate whether any threads are waiting
        std::atomic_uint64_t wait_count_{0};
        // std::size_t wait_count_{0};
        std::mutex wait_count_mutex_;
        std::condition_variable wait_condition_;
        // Flag to indicate whether the thread pool should stop or not
        bool stop_ = false;
    protected:
        // Constructor to creates a thread pool with given number of threads
        explicit ThreadPool(const size_t num_threads = std::thread::hardware_concurrency()) {
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

        // Enqueue task for execution by the thread pool
        void enqueue(std::function<void()> task)
        {
            {
                std::unique_lock lock(queue_mutex_);
                tasks_.emplace(std::move(task));
            }
            // Wake one worker
            cv_.notify_one();
            // Wake any waiters due to the queue state change
            wait_condition_.notify_all();
        }

        [[nodiscard]] size_t size() const {
            return threads_.size();
        }

        void resize(const size_t num_threads = std::thread::hardware_concurrency()) {
            if (threads_.size() != num_threads) {
                refresh(num_threads);
            }
        }

        // Resize the pool to a specified number of threads
        void refresh(const size_t num_threads = std::thread::hardware_concurrency()) {
            if (!threads_.empty()) {
                clear();
                threads_.clear();
                std::unique_lock lock(queue_mutex_); // there are no threads now. Is this lock really necessary?
                stop_ = false;
            }
            // Creating worker threads
            for (size_t i = 0; i < num_threads; ++i) {
                threads_.emplace_back([this] {
                    while (true) {
                        std::function<void()> task;
                        // The reason for putting the below code
                        // here is to unlock the queue before
                        // executing the task so that other
                        // threads can perform enqueue tasks
                        {
                            // Locking the queue so that data
                            // can be shared safely
                            std::unique_lock lock(queue_mutex_);
                            ++wait_count_;
                            // Notify any waiters due to the wait count change
                            wait_condition_.notify_all();

                            // Waiting until there is a task to
                            // execute or the pool is stopped
                            cv_.wait(lock, [this] {
                                return !tasks_.empty() || stop_;
                            });
                            --wait_count_;
                            // Notify any waiters due to the wait count change
                            wait_condition_.notify_all();

                            // exit the thread in case the pool
                            // is stopped and there are no tasks
                            if (stop_ && tasks_.empty()) {
                                return;
                            }
                            // Get the next task from the queue
                            task = std::move(tasks_.front());
                            tasks_.pop();
                            // std::stringstream s;
                            // s << "A thread woke up to do a job! ";
                            // s << wait_count_ << " waiting, " << tasks_.size() << " jobs remain";
                            // std::cout << s.str() << std::endl;

                            // Notify any waiters due to the queue change
                            wait_condition_.notify_all();

                        }
                        task();
                        wait_condition_.notify_all();
                    }
                });
            }
        }

        // Wait for all threads to finish their work
        void wait() {
            std::unique_lock lock(queue_mutex_);
            // std::cout << wait_count_ << " waiting " << " and " << tasks_.size() << " tasks" << std::endl;
            wait_condition_.wait(lock, [this] {
                // std::cout << wait_count_ << " of " << threads_.size() << " waiting and " << tasks_.size() << " tasks" << std::endl;
                return tasks_.empty() && wait_count_ >= threads_.size();
            });
            // std::cout << wait_count_ << " waiting " << " and " << tasks_.size() << " tasks" << std::endl;
            // lock.unlock();
        }
    private:
        // Stop all threads (in destructor or before resizing as part of a refresh)
        void clear() {
            {
                // Lock the queue to update the stop flag safely
                std::unique_lock lock(queue_mutex_);
                stop_ = true;
            }
            // Notify all threads
            cv_.notify_all();
            // Notify waiter so they re-evaluate
            wait_condition_.notify_all();
            // Joining all worker threads to ensure they have
            // completed their tasks
            for (auto& thread : threads_) {
                thread.join();
            }
        }
    };

    std::pair<size_t, size_t> thread_slice(size_t total, size_t threads, size_t thread);

}