/* This file is part of brille.

Copyright © 2019,2020 Greg Tucker <greg.tucker@stfc.ac.uk>

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

#include "debug.hpp"
using namespace brille;
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
  int brille::terminal_width(void){
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    auto width = csbi.srWindow.Right - csbi.srWindow.Left;
    return width > 1 ? width : 1<<15;
  }
  int brille::terminal_height(void){
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    return csbi.srWindow.Bottom - csbi.srWindow.Top;
  }
#else
  #include <sys/ioctl.h>
  #include <unistd.h>
  int brille::terminal_width(void){
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return w.ws_col > 0 ? w.ws_col : std::numeric_limits<int>::max();
  }
  int brille::terminal_height(void){
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return w.ws_row > 0 ? w.ws_row : std::numeric_limits<int>::max();
  }
#endif
