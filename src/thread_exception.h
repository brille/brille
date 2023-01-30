#ifndef BRILLE_THREAD_EXCEPTION_H
#define BRILLE_THREAD_EXCEPTION_H
#include <exception>
#include <mutex>
#include <sstream>

class ThreadException {
  size_t count_=0;
  std::exception_ptr ptr_=nullptr;
  std::mutex         lock_;
public:
  void rethrow(){
    if (count_>1) {
      std::ostringstream oss;
      oss << count_ << " exceptions occurred:\n";
      try {
        if (ptr_) std::rethrow_exception(ptr_);
      } catch (const std::exception & e){
        oss << e.what();
      }
      ptr_ = std::make_exception_ptr(std::runtime_error(oss.str()));
    }
    if (ptr_) std::rethrow_exception(ptr_);
  }
  void capture() {
    std::unique_lock<std::mutex> guard(lock_);
    if (count_++) {
      std::ostringstream oss;
      try {
        if (ptr_) std::rethrow_exception(ptr_);
      } catch (const std::exception & e){
        oss << e.what();
      }
      try {
        std::rethrow_exception(std::current_exception());
      } catch (const std::exception & e){
        oss << "\n" << e.what();
      }
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
