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

#ifndef BRILLE_ARRAY__HPP_
#define BRILLE_ARRAY__HPP_

namespace brille {
  // declare both Array and Array2 so that they can each contain conversion
  // constructors for the other.
  template<class T> class Array;
  template<class T> class Array2;

  // declare their itterators here
  template<class T> class ArrayIt;
  template<class T> class Array2It;
}
#endif
