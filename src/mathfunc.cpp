
#include "linear_algebra.h"

// special casting of floating point values to integers
int    cast_to_int(const int    a){ return a; }
int    cast_to_int(const double a){ return int( a<0 ? a-0.5 : a+0.5); }
double cast_to_dbl(const double a){ return a; }
double cast_to_dbl(const int    a){ return double(a); }
