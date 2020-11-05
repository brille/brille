/* Copyright 2019-2020 Greg Tucker
//
// This file is part of brille.
//
// brille is free software: you can redistribute it and/or modify it under the
// terms of the GNU Affero General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// brille is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with brille. If not, see <https://www.gnu.org/licenses/>.            */

#ifndef _NEIGHBOURS_H_
#define _NEIGHBOURS_H_
#include "array_latvec.hpp" // defines bArray

/*! Construct an bArray with 3 elements per array of all (2N+1)³-1
combinations of {-N,-N+1,…,N-1,N} for N=`extent`, skipping over (0,0,0).

Typically `extent`=1 and the returned bArray is
[(-1,-1,-1,),(-1,-1,0),(-1,-1,1),(-1,0,-1),…,(0,0,-1),(0,0,1),…,(1,0,1),(1,1,-1),(1,1,0),(1,1,1)]
*/
bArray<int> make_relative_neighbour_indices(const int extent=1);
/*! Construct an bArray with 4 elements per array of all (2N+1)⁴-1
combinations of {-N,-N+1,…,N-1,N} for N=`extent`, skipping over (0,0,0,0).
*/
bArray<int> make_relative_neighbour_indices4(const int extent=1);

bArray<int> make_relative_neighbour_indices_prime(const int extent=1);

#endif
