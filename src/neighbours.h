/*! \file */
#ifndef _NEIGHBOURS_H_
#define _NEIGHBOURS_H_
#include "arrayvector.h"

/*! Construct an ArrayVector with 3 elements per array of all (2N+1)³-1
combinations of {-N,-N+1,…,N-1,N} for N=`extent`, skipping over (0,0,0).

Typically `extent`=1 and the returned ArrayVector is
[(-1,-1,-1,),(-1,-1,0),(-1,-1,1),(-1,0,-1),…,(0,0,-1),(0,0,1),…,(1,0,1),(1,1,-1),(1,1,0),(1,1,1)]
*/
ArrayVector<int> make_relative_neighbour_indices(const int extent=1);
/*! Construct an ArrayVector with 4 elements per array of all (2N+1)⁴-1
combinations of {-N,-N+1,…,N-1,N} for N=`extent`, skipping over (0,0,0,0).
*/
ArrayVector<int> make_relative_neighbour_indices4(const int extent=1);

#endif
