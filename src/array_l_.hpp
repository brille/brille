#ifndef _BRILLE_ARRAY_L__HPP_
#define _BRILLE_ARRAY_L__HPP_
// Declare the array types without defining them
#include "array_l_pre.hpp"
// The arrays need to know about each other so their constructors can act as converts
#include "array_lvec.hpp"
// The following includes must happen *after* Array and Array2 are fully defined
// (so array.hpp and array2.hpp can not include circular references back to array_.hpp)
#include "array_l_post.hpp"

#endif