#ifndef BRILLE_THREAD_EXCEPTION_H
#define BRILLE_THREAD_EXCEPTION_H
#include <exception>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <utility>

/*! \brief Collect exceptions thrown on worker threads, to rethrow on the calling thread

An exception that escapes a thread's function calls std::terminate, so worker
code catches everything and hands it to `capture`. After the workers finish,
the calling thread calls `rethrow`, which throws the captured exception (or a
std::runtime_error combining all messages, if more than one was captured) and
resets the object so that it can be reused.
*/
class ThreadException {
  size_t count_=0;
  std::exception_ptr ptr_=nullptr;
  std::mutex         lock_;
  static std::string message(const std::exception_ptr & ptr) {
    try {
      if (ptr) std::rethrow_exception(ptr);
    } catch (const std::exception & e){
      return e.what();
    } catch (...) {
      return "unknown exception";
    }
    return "";
  }
public:
  void rethrow(){
    std::exception_ptr ptr;
    size_t count;
    {
      std::unique_lock<std::mutex> guard(lock_);
      ptr = std::exchange(ptr_, nullptr);
      count = std::exchange(count_, 0);
    }
    if (count > 1) {
      std::ostringstream oss;
      oss << count << " exceptions occurred:\n" << message(ptr);
      throw std::runtime_error(oss.str());
    }
    if (ptr) std::rethrow_exception(ptr);
  }
  //! Call only from within a catch block
  void capture() {
    std::unique_lock<std::mutex> guard(lock_);
    if (count_++) {
      std::ostringstream oss;
      oss << message(ptr_) << "\n" << message(std::current_exception());
      ptr_ = std::make_exception_ptr(std::runtime_error(oss.str()));
    } else {
      ptr_ = std::current_exception();
    }
  }
  template <typename Function, typename... Parameters>
  void run(Function f, Parameters... params)
  {
    try {
      f(params...);
    } catch (...) {
      capture();
    }
  }
};

#endif // BRILLE_THREAD_EXCEPTION_H
