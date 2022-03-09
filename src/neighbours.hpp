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

#ifndef BRILLE_NEIGHBOURS_H_
#define BRILLE_NEIGHBOURS_H_
/*! \file
    \author Greg Tucker
    \brief Functions to build lists of neighbour subscripted indices
*/
#include "array_.hpp" // defines bArray
namespace brille {

/*! \brief Build a list of relative neighbour subscripted indices

\param extent The magnitude of the relative extent to include

Construct an Array2 with 3 elements per array of all (2N+1)³-1
combinations of {-N,-N+1,…,N-1,N} for N=`extent`, skipping over (0,0,0).

Typically `extent`=1 and the returned bArray is
[(-1,-1,-1,),(-1,-1,0),(-1,-1,1),(-1,0,-1),…,(0,0,-1),(0,0,1),…,(1,0,1),(1,1,-1),(1,1,0),(1,1,1)]
*/
bArray<int> make_relative_neighbour_indices(int extent=1);
/*! \brief Build a 4-D list of relative neighbour subscripted indices

\param extent The magnitude of the relative extent to include

Construct an Array2 with 4 elements per array of all (2N+1)⁴-1
combinations of {-N,-N+1,…,N-1,N} for N=`extent`, skipping over (0,0,0,0).
*/
bArray<int> make_relative_neighbour_indices4(int extent=1);
/*! \brief Build a special list of relative neighbour subscripted indices

\param extent The maginitude of the relative extent to include

Starting with the neighbours differing in one index only:
  [(-1,0,0),(1,0,0),(0,-1,0),(0,1,0),(0,0,-1),(0,0,1)]

then those differing in two indices:
  [(-1,-1,0),(-1,1,0),(1,-1,0),(1,1,0), (-1,0,-1),...(1,0,1), (0,-1,-1),...,(0,1,1)]

and finally the remaining neighbours differing in all three indices
  [(-1,-1,-1),(-1,-1,1),(-1,1,-1),(-1,1,1),...,(1,1,1)]

build up a list of all (2N+1)³-1
combinations of {-N,-N+1,…,N-1,N} for N=`extent`, skipping over (0,0,0).

\note This function was intended to help produce first Brillouin zones faster
      by providing the list of neighbouring reciprocal lattice points in a
      specific order. As the voro_search algorithm sorts the neighbour list this
      special ordering is no longer necessary so this function should typically
      not be used.
\see make_relative_neighbour_indices
*/
bArray<int> make_relative_neighbour_indices_prime(int extent=1);

} // namespace brille
#endif
