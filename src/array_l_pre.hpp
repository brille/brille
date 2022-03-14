#ifndef BRILLE_ARRAY_L_PRE_HPP_
#define BRILLE_ARRAY_L_PRE_HPP_
/*! \file
    \brief Forward declarations for Array* classes which can be constructed
           from each other.
*/

namespace brille::lattice {
  // declare both Array and Array2 so that they can each contain conversion
  // constructors for the other.
  template<class T> class LVec;
}
#endif