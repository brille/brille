/* This file is part of brille.

Copyright © 2020 Greg Tucker <greg.tucker@stfc.ac.uk>

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

#ifndef _BRILLE_ARRAY__HPP_
#define _BRILLE_ARRAY__HPP_
/*! \file
    \brief Forward declarations for Array* classes which can be constructed
           from each other.
*/
// Declare the array types without defining them
#include "array_pre.hpp"
// The arrays need to know about each other so their constructors can act as converts
#include "array.hpp"
#include "array2.hpp"
// The following includes must happen *after* Array and Array2 are fully defined
// (so array.hpp and array2.hpp can not include circular references back to array_.hpp)
#include "array_post.hpp"

#endif
