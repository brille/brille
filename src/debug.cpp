#include "debug.h"
/*******************************************************************************
| These functions could be defined static and/or inline and then placed in the |
| debug header file, however doing so forces all other source files to include |
| the windows.h header file. Such behaviour is undesirable as it wreaks havok  |
| on standard functions (e.g., std::max) by replacing them with macros.        |
*******************************************************************************/
#if defined(_MSC_VER) || defined(__MINGW32__)
// #ifdef _MSC_VER
  // #define NOMINMAX
  #include <windows.h>
  int terminal_width(void){
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    return csbi.srWindow.Right - csbi.srWindow.Left;
  }
  int terminal_height(void){
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    return csbi.srWindow.Bottom - csbi.srWindow.Top;
  }
#else
  #include <sys/ioctl.h>
  #include <unistd.h>
  int terminal_width(void){
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return w.ws_col;
  }
  int terminal_height(void){
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return w.ws_row;
  }
#endif
