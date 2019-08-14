#ifndef _DEBUG_H_
#define _DEBUG_H_
#include <iostream>
#include <string>

// #define VERBOSE_DEBUG
#define DEBUG // comment-out for no debugging output

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
  };
  template<typename... L> void println(const std::string& fnc, L... l){
    if (last_function.compare(fnc)){
      last_function = fnc;
      std::cout << fnc << std::endl;
    }
    this->inner_print(l...);
    std::cout << std::endl;
  };
private:
  template<typename T, typename... L> void inner_print(const T& x, L... l){
    std::cout << std::to_string(x);
    this->inner_print(l...);
  };
  template<typename... L> void inner_print(const std::string& x, L... l){
    std::cout << x;
    this->inner_print(l...);
  };
  template<typename... L> void inner_print(const char* x, L... l){
    std::cout << x;
    this->inner_print(l...);
  };
  void inner_print(void){};
};

// replace GNU __PRETTY_FUNCTION__ by __FUNCSIG__ on Windows (MSVC)
#ifdef _WIN32
  #define __PRETTY_FUNCTION__ __FUNCSIG__
#endif

#if defined(VERBOSE_DEBUG) || defined(DEBUG)
  static DebugPrinter _debug_printer("");
  #ifdef DEBUG
    #define status_update(...) _debug_printer.println("", __VA_ARGS__)
    #define verbose_status_update(...)
  #endif
  #ifdef VERBOSE_DEBUG
    #define status_update(...) _debug_printer.println(__PRETTY_FUNCTION__, __VA_ARGS__)
    #define verbose_status_update(...) _debug_printer.println(__PRETTY_FUNCTION__, __VA_ARGS__)
  #endif
  #define debug_exec(...) __VA_ARGS__
#else
  #define status_update(...)
  #define verbose_status_update(...)
  #define debug_exec(...)
#endif

#endif //_DEBUG_H_
