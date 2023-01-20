/* This file is part of brille.

Copyright © 2019-2022 Greg Tucker <gregory.tucker@ess.eu>

brille is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

brille is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with brille. If not, see <https://www.gnu.org/licenses/>.            */
#ifndef BRILLE_DEBUG_HPP_
#define BRILLE_DEBUG_HPP_
/*! \file
    \author Greg Tucker
    \brief Defines utilities for simple logging output to stdout
*/
#include <stdio.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <array>
#include <utility>
#include <vector>
#include <complex>
#include <chrono>
//#define VERBOSE_DEBUG
//#define DEBUG // comment-out for no debugging output
//#define PROFILING

#define PRINTING_PRECISION 16
namespace brille {
  /*! \brief Absolute value of signed or unsigned values

  \param x A value
  \returns |x|, the absolute value of x

  `std::abs()` is not guaranteed to be defined for unsigned integers.
  This template function returns its input for unsigned `T` or non-scalar `T`
  (like the element reference into a std::vector<bool> which is *not* a `bool`!)
  and otherwise returns `std::abs(x)` for signed scalar T.
  */
  template<typename T>
  inline T abs(const T & x) {
      if constexpr (std::is_unsigned_v<T>||!std::is_scalar_v<T>){
          return x;
      } else {
          return std::abs(x);
      }
  }

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
  template<class T>
  inline std::enable_if_t<!std::is_unsigned_v<T>, bool>
  is_negative(const T x){ return x < 0; }
  template<class T>
  inline std::enable_if_t<std::is_unsigned_v<T>, bool>
  is_negative(const T){ return false; }
#endif
#if defined(DOXYGEN_SHOULD_SKIP_THIS)
  /*! \brief Return the if a scalar is negative

  \param x A scalar value
  \return The equivalent of `x < 0`

  Some compilers complain about the comparison `x<0` for unsigned integer x.
  This function has been overloaded to avoid that comparison through template
  substitution failure.
  */
  template<class T> bool is_negative(const T x);
#endif

/*! \brief Determine the width of the current terminal window

\return The terminal width in characters or 2¹⁵ if the terminal width is zero.
*/
int terminal_width();
/*! \brief Determine the height of the current terminal window

\return The terminal height in lines or 2¹⁵ if the terminal height is zero.
*/
int terminal_height();

/*! \brief A utility structure to identify containers in templates

The `value` type is true for `std::vector` and `std::array` and otherwise false.
*/
template <typename T> struct is_container {
  enum { value = false };
};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template <typename T> struct is_container<std::vector<T>> {
  enum { value = true };
};
template <typename T, size_t N> struct is_container<std::array<T, N>> {
  enum { value = true };
};
template<bool C, typename T> using enable_if_t = typename std::enable_if<C,T>::type;
#endif

/*! \brief Construct a string representation of non-container input */
template<typename T, typename=typename std::enable_if<!is_container<T>::value>::type>
std::string my_to_string(const T x, const size_t width=0){
  std::ostringstream streamobj;
  size_t w{width};
  if constexpr (!std::is_integral<T>::value){
    streamobj << std::fixed;
    streamobj << std::setprecision(PRINTING_PRECISION);
    if (w > PRINTING_PRECISION) w -= 1 + PRINTING_PRECISION; // account for the decimal mark and precision
  }
  // char may or may not be signed, depending on the system
  if constexpr (std::is_base_of<char,T>::value || (std::is_integral<T>::value && std::is_unsigned<T>::value) ){
    if (w) streamobj << std::setw(w);
    streamobj << x;
  } else {
    if (w) streamobj << std::setw(w-1); // -1 to account for the sign
    streamobj << (brille::is_negative(x) ? "-" : " ") << brille::abs(x);
  }
  return streamobj.str();
}
/*! \brief Construct a string representation of complex non-container input */
template<typename T, typename=typename std::enable_if<!is_container<T>::value>::type>
std::string my_to_string(const std::complex<T> x, const size_t width=0){
  T r = std::real(x), i=std::imag(x);
  std::ostringstream streamobj;
  size_t w{width};
  if (!std::is_integral<T>::value){
    streamobj << std::fixed;
    streamobj << std::setprecision(PRINTING_PRECISION);
    if (w > 1 + 2 * PRINTING_PRECISION) w -= 2 + 2 * PRINTING_PRECISION;
  }
  if (!std::is_integral<T>::value || std::is_signed<T>::value){
    if (w>3) streamobj << std::setw(w-3); // -3 for -±i
    streamobj << (brille::is_negative(r) ? "-" : " ") << brille::abs(r);
    streamobj << (brille::is_negative(i) ? "-i" : "+i") << brille::abs(i);
  } else {
    if (w>2) streamobj << std::setw(w-2); // -2 for +i
    streamobj << r << "+i" << i;
  }
  return streamobj.str();
}
/*! \brief Construct a string representation of a vector of containers */
template<typename T, template<class> class C,
        typename=typename std::enable_if<!is_container<T>::value>::type,
        typename=typename std::enable_if<is_container<C<T>>::value>::type>
std::string my_to_string(const std::vector<C<T>>& v, const size_t){
  std::string s;
  for (C<T> x: v) s += my_to_string(x) + "\n";
  return s;
}
/*! \brief Construct a string representation of a container input */
template<typename T, typename=typename std::enable_if<is_container<T>::value>::type>
std::string my_to_string(const T & a, const size_t w=0){
  std::string s;
  for (auto x: a) s += my_to_string(x, w);
  return s;
}
template<class T, class R> std::string my_to_string(const std::pair<T,R>& p, const size_t w=0){
  size_t h = (w > 2u) ? (w - 1u) / 2u : 0u;  // ensure we do not overflow to h = 2⁶³ -- GitHub issue #79
  std::string s = my_to_string(p.first, h) + "—" + my_to_string(p.second, h);
  return s;
}

template<typename T>
[[nodiscard]] std::string list_to_string(const std::vector<T>& v){
  std::ostringstream s;
  s << "[";
  bool o{};
  for (const auto & x: v) s << (o ? ", " : (o=1, "")) << my_to_string(x);
  s << "]";
  return s.str();
}
template<typename T, size_t N>
[[nodiscard]] std::string list_to_string(const std::array<T,N>& v){
  std::ostringstream s;
  s << "[";
  bool o{};
  for (const auto & x: v) s << (o ? ", " : (o=1, "")) << my_to_string(x);
  s << "]";
  return s.str();
}
template<typename T>
[[nodiscard]] std::string lists_to_string(const std::vector<std::vector<T>>& vv){
  std::ostringstream s;
  s << "[";
  bool o{};
  for (const auto & v: vv) s << (o ? ",\n " : (o=1, "")) << list_to_string(v);
  s << "]";
  return s.str();
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  /*! The length of the string representation of a scalar */
  template<typename T>
  static enable_if_t<!is_container<T>::value, size_t>
  max_element_length(const T& v){
    return my_to_string(v).size();
  }
  /* The longest string representation of the elements of a container */
  template<typename T>
  static enable_if_t<is_container<T>::value, size_t>
  max_element_length(const T& v){
    size_t l=0;
    for (auto x: v){
      size_t t = max_element_length(x);
      if (t > l) l = t;
    }
    return l;
  }
#else
  /*! \brief Find the maximum string-representation length of a scalar or element

  \param soc A single scalar value or a container of scalars
  \return The length of the longest string representation of the scalar(s)
  */
  template<class T>
  static size_t max_element_length(const T& soc);
#endif

/*! \brief Construct a timestamp string

\param time The `time_point` to produce a string timestamp from
\return A string of the form `"[YYYY-MM-DD HH:MM:SS:mmm]"`
*/
template<typename T>
std::string time_to_string(std::chrono::time_point<T> time) {
  using namespace std;
  using namespace std::chrono;

  time_t curr_time = T::to_time_t(time);
  char sRep[100];
  strftime(sRep, sizeof(sRep), "%Y-%m-%d %H:%M:%S", localtime(&curr_time));

  typename T::duration since_epoch = time.time_since_epoch();
  seconds s = duration_cast<seconds>(since_epoch);
  since_epoch -= s;
  milliseconds milli = duration_cast<milliseconds>(since_epoch);

  stringstream buffer;
  buffer << '[' << sRep << ':' << setfill('0') << setw(3) << milli.count() << "] ";
  return buffer.str();
}

/*! \brief A class for simple standard output logging

Provides variable input argument methods for displaying simple logging messages.
\note A single namespace-wide object is used in conjunction with vararg macros.
\see info_update, info_update_if, debug_update_if, debug_update, verbose_update
\see verbose_update_if, profile_update, profile_update_if
*/
class DebugPrinter{
  std::string last_function; // replace this with the stack?
  bool _silenced;
  bool emit_datetime;
  size_t _before;
public:
  DebugPrinter(std::string s)
  : last_function(std::move(s)), _silenced(false), emit_datetime(false), _before(0)
  {
    #if defined(PROFILING)
    emit_datetime = true;
    #endif
  };
  [[nodiscard]] bool silenced() const {return _silenced;}
  bool silenced(bool slc) {_silenced = slc; return _silenced;}
  bool silence() {_silenced = true; return _silenced;}
  bool unsilence() {_silenced = false; return !_silenced;}
  [[nodiscard]] bool datetime() const {return emit_datetime;}
  bool datetime(bool edt) {emit_datetime = edt; return emit_datetime;}
  template<typename... L> void print(const std::string& fnc, L... l){
    if (!_silenced){
      size_t torem{0};
      if (last_function.compare(fnc)){
        last_function = fnc;
        std::cout << std::endl << fnc << std::endl;
      }
      if (emit_datetime){
        // auto tp = std::chrono::system_clock::now();
        // auto tt = std::chrono::system_clock::to_time_t(tp);
        // std::tm tm = *std::localtime(&tt);
        // std::stringstream buffer;
        // buffer << std::put_time(&tm, "%FT%T ");
        // std::string tstr = buffer.str();
        std::string tstr = time_to_string(std::chrono::system_clock::now());
        this->inner_print(tstr); // print without extra lead in spaces
        _before += (torem = tstr.size()); // add to lead in spacing
      }
      this->inner_print(l...);
      _before -= torem;
    }
  }
  template<typename... L> void println(const std::string& fnc, L... l){
    if (!_silenced){
      this->print(fnc, l...);
      std::cout << std::endl;
    }
  }
private:
  [[nodiscard]] std::string lead_in() const {
    std::stringstream buffer;
    for (size_t i=0; i<_before; ++i) buffer << " ";
    return buffer.str();
  }
  template<typename... L>
  void inner_print(const std::string& x, L... l){
    std::cout << x;
    this->inner_print(l...);
  }
  template<typename... L>
  void inner_print(const char* x, L... l){
    std::cout << x;
    this->inner_print(l...);
  }
  template<typename T, typename... L>
  enable_if_t<!is_container<T>::value, void> inner_print(const T& x, L... l){
    std::cout << my_to_string(x);
    this->inner_print(l...);
  }
  template<typename T, typename... L>
  enable_if_t<!is_container<T>::value, void> inner_print(const std::vector<T>& x, L... args){
    size_t l = max_element_length(x);
    // size_t n = std::sqrt(x.size());
    int w = terminal_width();
    if (static_cast<int>(_before) < w) w -= static_cast<int>(_before);
    if (l) w /= static_cast<int>(l)+1;
    int count = 0;
    std::ostringstream s;
    bool o{false};
    for (const auto & y: x){
      s << (o ? ", " : (o=true, "")) << my_to_string(y, l) << (++count % w ? "" : "\n" + this->lead_in());
    }
    this->inner_print(s.str(), args...);
  }
  template<typename T, size_t N, typename... L>
  enable_if_t<!is_container<T>::value, void> inner_print(const std::array<T,N>& x, L... args){
    size_t l= max_element_length(x);
    auto w = static_cast<size_t>(terminal_width());
    if (_before < w) w -= _before;
    if (l) w /= l+1;
    size_t count = 0;
    std::ostringstream s;
    bool o{false};
    for (const auto & y: x){
      s << (o ? ", " : (o=true, "")) << my_to_string(y, l) << (++count % w ? "" : "\n" + this->lead_in());
    }
    this->inner_print(s.str(), args...);
  }
  template<typename T, typename... L>
  void inner_print(const std::vector<std::vector<T>>& vv, L... args){
    size_t l = max_element_length(vv);
    auto w = static_cast<size_t>(terminal_width());
    if (_before < w) w -= _before;
    size_t num;
    std::string s;
    if (l) w /= l+1;
    s += "[";
    ++_before;
    for (size_t i=0; i < vv.size(); ++i){
      const auto & v{vv[i]};
      num = 0;
      s += "[";
      for (size_t j=0; j < v.size(); ++j){
        s += my_to_string(v[j], l);
        if (j + 1 < v.size()) s += ',';
        if (!(++num %w)) s += "\n" + this->lead_in();
      }
      s += "]";
      s += (i + 1 < vv.size()) ? "," : "]";
      s += "\n" + this->lead_in();
    }
    --_before;
    this->inner_print(s, args...);
  }
  template<typename T, size_t N, typename... L>
  enable_if_t<!is_container<T>::value, void> inner_print(const std::vector<std::array<T,N>>& x, L... args){
    size_t l = max_element_length(x);
    auto w = static_cast<size_t>(terminal_width());
    if (_before < w) w -= _before;
    size_t num;
    std::string s;
    if (l) w /= l+1;
    for (size_t i=0; i<x.size(); num=0, ++i){
      for (auto y: x[i]){
        s += my_to_string(y, l) + ",";
        if (!(++num % w)) s += "\n" + this->lead_in();
      }
      s += "\n" + this->lead_in();
    }
    this->inner_print(s, args...);
  }
  void inner_print(){};
};

//! The single namespace wide `DebugPrinter` used with the logging macros.
// static DebugPrinter printer("");
extern DebugPrinter printer;

/*! \brief A simple timer for use in debugging and profiling

A resetable stopwatch which can be used to repeatedly time the same operation
and easily provide average per-operation time as well as uncertainty in the
computed average.
*/
template<typename TimeT = std::chrono::milliseconds>
class Stopwatch{
  typedef std::chrono::high_resolution_clock ClockT;
private:
    std::chrono::time_point<ClockT> _start, _end, _split;
    size_t presses;
public:
  /*! \brief Initialize and start the timer */
    Stopwatch(): presses(0u){
      tic();
    }
  /*! \brief Reset the timer

  Sets the start and end time of the timer to now, and resets the number of
  timed iterations to zero.
  */
    void tic(){
      presses = 0u;
      _start = _end = _split = ClockT::now();
    }
  /*! \brief Stop the timer

  Sets the end timepoint of the timer, increments the number of iterations
  timed, and returns the difference between the end and start times.
  */
    double toc(){
      _end = ClockT::now();
      ++presses;
      return elapsed();
    }
    /*! \brief Calculate the time elapsed betwen starting and stopping the timer */
    [[nodiscard]] double elapsed() const {
      auto delta = std::chrono::duration_cast<TimeT>(_end - _start);
      return static_cast<double>(delta.count());
    }
    /*! \brief Calculate the average time per iteration of the timer */
    [[nodiscard]] double average() const {
      return elapsed()/static_cast<double>(presses);
    }
    /*! \brief Calculate the uncertainty in the average time per iteration */
    [[nodiscard]] double jitter() const {
      return std::sqrt(elapsed())/static_cast<double>(presses);
    }
    /*! \brief Return the time since the timer started or last split call

    Increments the iteration counter but does not set the timer end timepoint.
    */
    double split(){
      auto new_split = ClockT::now();
      auto delta = std::chrono::duration_cast<TimeT>(new_split - _split);
      _split = new_split;
      ++presses;
      return static_cast<double>(delta.count());
    }
};


} //end namespace brille

// replace GNU __PRETTY_FUNCTION__ by __FUNCSIG__ on Windows (MSVC)
#ifdef _MSC_VER
  #define __PRETTY_FUNCTION__ __FUNCSIG__
#endif
#ifdef VERBOSE_DEBUG
  //! Prepend the containing function name to logging messages when in VERBOSE_DEBUG mode
  #define PREPEND_LOG __PRETTY_FUNCTION__
#else
  //! Prepend nothing before logging messages if not in VERBOSE_DEBUG mode
  #define PREPEND_LOG ""
#endif

//! A macro for informational logging to standard output
#define info_update(...) brille::printer.println(PREPEND_LOG, __VA_ARGS__)
//! A macro for conditional informational logging to standard output
#define info_update_if(condition, ...) if (condition) brille::printer.println(PREPEND_LOG, __VA_ARGS__)

#if defined(VERBOSE_DEBUG) || defined(DEBUG)
  //! A macro to execute instructions when in DEBUG mode
  #define debug_exec(...) __VA_ARGS__
  //! A macro for debug logging to standard output
  #define debug_update(...) brille::printer.println(PREPEND_LOG, __VA_ARGS__)
  //! A macro for conditional debug logging to standard output
  #define debug_update_if(condition, ...) if (condition) brille::printer.println(PREPEND_LOG, __VA_ARGS__)
#else
  //! A macro to elide instructions when not in DEBUG mode
  #define debug_update_if(...)
  //! A macro to elide debug logging when not in DEBUG mode
  #define debug_update(...)
  //! A macro to elide conditional debug logging when not in DEBUG mode
  #define debug_exec(...)
#endif
#ifdef VERBOSE_DEBUG
  //! A macro for verbose debug logging to standard output
  #define verbose_update(...) brille::printer.println(PREPEND_LOG, __VA_ARGS__)
  //! A macro for conditional verbose debug logging to standard output
  #define verbose_update_if(condition, ...) if (condition) brille::printer.println(PREPEND_LOG, __VA_ARGS__)
#else
  //! A macro to elide verbose debug logging when not in VERBOSE_DEBUG mode
  #define verbose_update(...)
  //! A macro to elide conditional verbose debug logging when not in VERBOSE_DEBUG mode
  #define verbose_update_if(...)
#endif

#if defined(PROFILING)
  //! A macro for printing profiling messages to standard output
  #define profile_update(...) brille::printer.println(PREPEND_LOG, __VA_ARGS__)
  //! A macro for conditionally printing profiling messages to standard output
  #define profile_update_if(condition, ...) if (condition) brille::printer.println(PREPEND_LOG, __VA_ARGS__)
#else
  //! A macro to elide profiling output when not in PROFILING mode
  #define profile_update(...)
  //! A macro to elide conditional profiling output when not in PROFILING mode
  #define profile_update_if(...)
#endif

#endif //_DEBUG_H_
