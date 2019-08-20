#ifndef _DEBUG_H_
#define _DEBUG_H_
#include <stdio.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <complex>

// #define VERBOSE_DEBUG
#define DEBUG // comment-out for no debugging output

int terminal_width(void);
int terminal_height(void);

template <typename T> struct is_AV {
  enum { value = false };
};
template <typename T> struct is_AV<std::vector<T>> {
  enum { value = true };
};
template <typename T, size_t N> struct is_AV<std::array<T, N>> {
  enum { value = true };
};
template <typename T> struct is_not_AV {
  enum { value = true };
};
template <typename T> struct is_not_AV<std::vector<T>> {
  enum { value = false };
};
template <typename T, size_t N> struct is_not_AV<std::array<T, N>> {
  enum { value = false };
};


template<bool C, typename T> using enable_if_t = typename std::enable_if<C,T>::type;
template<typename T> static enable_if_t< std::is_integral<T>::value && std::is_unsigned<T>::value, T> local_abs(T x) { return x; }
template<typename T> static enable_if_t<!std::is_integral<T>::value ||   std::is_signed<T>::value, T> local_abs(T x) { return std::abs(x); }

template<typename T, typename=typename std::enable_if<is_not_AV<T>::value>::type>
const std::string my_to_string(const T x, const size_t w=0){
  std::ostringstream streamobj;
  if (!std::is_integral<T>::value){
    streamobj << std::fixed;
    streamobj << std::setprecision(4);
  }
  // char may or may not be signed, depending on the system
  if (std::is_base_of<char,T>::value || (std::is_integral<T>::value && !std::is_signed<T>::value) ){
    if (w) streamobj << std::setw(w);
    streamobj << x;
  } else {
    if (w) streamobj << std::setw(w-1); // -1 to account for the sign
    streamobj << (x<0 ? "-" : " ") << local_abs(x);
  }
  return streamobj.str();
}
template<typename T, typename=typename std::enable_if<is_not_AV<T>::value>::type>
const std::string my_to_string(const std::complex<T> x, const size_t w=0){
  T r = std::real(x), i=std::imag(x);
  std::ostringstream streamobj;
  if (!std::is_integral<T>::value){
    streamobj << std::fixed;
    streamobj << std::setprecision(4);
  }
  if (!std::is_integral<T>::value || std::is_signed<T>::value){
    if (w>3) streamobj << std::setw(w-3); // -3 for -±i
    streamobj << (r<0 ? "-" : " ") << local_abs(r);
    streamobj << (std::signbit(i) ? "-i" : "+i") << local_abs(i);
  } else {
    if (w>2) streamobj << std::setw(w-2); // -2 for +i
    streamobj << r << "+i" << i;
  }
  return streamobj.str();
}
template<typename T, template<class> class C,
        typename=typename std::enable_if<is_not_AV<T>::value>::type,
        typename=typename std::enable_if<is_AV<C<T>>::value>::type>
const std::string my_to_string(const std::vector<C<T>>& v, const size_t w=0){
  std::string s;
  for (C<T> x: v) s += my_to_string(x) + "\n";
  return s;
}
template<typename T, typename=typename std::enable_if<is_AV<T>::value>::type>
const std::string my_to_string(const T & a, const size_t w=0){
  std::string s;
  for (auto x: a) s += my_to_string(x, w);
  return s;
}

template<typename T> static enable_if_t<is_not_AV<T>::value, size_t> max_element_length(const T& v){
  return my_to_string(v).size();
}
template<typename T> static enable_if_t<is_AV<T>::value, size_t> max_element_length(const T& v){
  size_t l=0, t=0;
  for (auto x: v){
    t = max_element_length(x);
    if (t > l) l = t;
  }
  return l;
}

class DebugPrinter{
  std::string last_function; // replace this with the stack?
public:
  DebugPrinter(const std::string& s): last_function(s) {};
  template<typename... L> void print(const std::string& fnc, L... l){
    if (last_function.compare(fnc)){
      last_function = fnc;
      std::cout << fnc << std::endl;
    }
    this->inner_print(l...);
  }
  template<typename... L> void println(const std::string& fnc, L... l){
    if (last_function.compare(fnc)){
      last_function = fnc;
      std::cout << fnc << std::endl;
    }
    this->inner_print(l...);
    std::cout << std::endl;
  }
private:
  template<typename T, typename... L>
  enable_if_t<is_not_AV<T>::value, void> inner_print(const T& x, L... l){
    // std::cout << std::to_string(x);
    std::cout << x;
    this->inner_print(l...);
  }
  // template<typename... L> void inner_print(const std::string& x, L... l){
  //   std::cout << x;
  //   this->inner_print(l...);
  // }
  // template<typename... L> void inner_print(const char* x, L... l){
  //   std::cout << x;
  //   this->inner_print(l...);
  // }
  template<typename T, typename... L>
  enable_if_t<is_not_AV<T>::value, void> inner_print(const std::vector<T>& x, L... args){
    size_t l = max_element_length(x);
    int w = terminal_width();
    if (l) w /= static_cast<int>(l)+1;
    int count = 0;
    std::string s;
    for (auto y: x){
      s += " " + my_to_string(y, l);
      if (!(++count % w)) s += "\n";
    }
    this->inner_print(s, args...);
  }
  template<typename T, typename... L>
  void inner_print(const std::vector<std::vector<T>>& vv, L... args){
    for (auto v: vv) this->inner_print(v,"\n");
    this->inner_print(args...);
  }
  template<typename T, size_t N, typename... L>
  enable_if_t<is_not_AV<T>::value, void> inner_print(const std::vector<std::array<T,N>>& x, L... args){
    size_t l = max_element_length(x);
    size_t w = static_cast<size_t>(terminal_width());
    size_t num;
    std::string s;
    if (N==9){
      num = (l) ? w/(3*l+4) : w/3;
      for (size_t i=0; i<x.size(); i+=num){
        for (int a=0; a<3; ++a){
          for (size_t j=0; j<num && (i+j)<x.size(); ++j){
            for (int b=0; b<3; ++b) s += " " + my_to_string(x[i+j][a*3+b], l);
            s += " ";
          }
          s += "\n";
        }
        s += "\n";
      }
    } else {
      if (l) w /= l+1;
      for (size_t i=0; i<x.size(); num=0, ++i){
        for (auto y: x[i]){
          s += " " + my_to_string(y, l);
          if (!(++num % w)) s += "\n";
        }
        s += "\n";
      }
    }
    this->inner_print(s, args...);
  }
  void inner_print(void){};
};

// replace GNU __PRETTY_FUNCTION__ by __FUNCSIG__ on Windows (MSVC)
#ifdef _MSC_VER
  #define __PRETTY_FUNCTION__ __FUNCSIG__
#endif

#if defined(VERBOSE_DEBUG) || defined(DEBUG)
  static DebugPrinter _debug_printer("");
  #ifdef DEBUG
    #define status_update_if(tf, ...) if (tf) _debug_printer.println("", __VA_ARGS__)
    #define status_update(...) _debug_printer.println("", __VA_ARGS__)
    #define verbose_status_update(...)
  #endif
  #ifdef VERBOSE_DEBUG
    #define status_update_if(tf, ...) if (tf) _debug_printer.println(__PRETTY_FUNCTION__, __VA_ARGS__)
    #define status_update(...) _debug_printer.println(__PRETTY_FUNCTION__, __VA_ARGS__)
    #define verbose_status_update(...) _debug_printer.println(__PRETTY_FUNCTION__, __VA_ARGS__)
  #endif
  #define debug_exec(...) __VA_ARGS__
#else
  #define status_update_if(...)
  #define status_update(...)
  #define verbose_status_update(...)
  #define debug_exec(...)
#endif

#endif //_DEBUG_H_
