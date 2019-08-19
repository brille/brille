#ifndef _DEBUG_H_
#define _DEBUG_H_
#include <iostream>
#include <string>

// #define VERBOSE_DEBUG
#define DEBUG // comment-out for no debugging output

#ifdef _MSC_VER
  #include <windows.h>
  inline int terminal_width(void){
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    return csbi.srWindow.Right - csbi.srWindow.Left +1;
  }
  inline int terminal_height(void){
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    return csbi.srWindow.Bottom - csbi.srWindow.Top + 1;
  }
#else
  #include <sys/ioctl.h>
  #include <unistd.h>
  inline int terminal_width(void){
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return w.ws_col;
  }
  inline int terminal_height(void){
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return w.ws_row;
  }
#endif

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
  template<typename T, typename... L> void inner_print(const T& x, L... l){
    std::cout << std::to_string(x);
    this->inner_print(l...);
  }
  template<typename... L> void inner_print(const std::string& x, L... l){
    std::cout << x;
    this->inner_print(l...);
  }
  template<typename... L> void inner_print(const char* x, L... l){
    std::cout << x;
    this->inner_print(l...);
  }
  template<typename T, typename... L> void inner_print(const std::vector<T>& x, L... l){
    int w = terminal_width()/10;
    int count = 0;
    for (auto y: x){
      std::cout << " " << std::to_string(y);
      if (count++ >= w){
        count=0;
        std::cout << std::endl;
      }
    }
    this->inner_print(l...);
  }
  template<typename T, int N, typename... L> void inner_print(const std::vector<std::array<T,N>>& x, L... l){
    int w = terminal_width()/10;
    int i=0, num;
    if (N==9){
      num=w/3;
      while (i<x.size()){
        for (int a=0; a<3; ++a){
          for (int j=0; j<num; ++j) if(i+j < x.size()){
            for (int b=0; b<3; ++b) std::cout << " " << std::to_string(x[i+j][a*3+b]);
            std::cout << " ";
          }
        std::cout << std::endl;
        }
        i += num;
      }
    } else {
      for (i=0; i<x.size(); ++i){
        num = 0;
        for (auto y: x[i]){
          std::cout << " " << std::to_string(y);
          if (num++ >= w){
            num=0;
            std::cout << std::endl;
          }
        }
        std::cout << std::endl;
      }
    }
    this->inner_print(l...);
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
