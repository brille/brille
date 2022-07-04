#ifndef BRILLE_ARRAY_PRE_HPP_
#define BRILLE_ARRAY_PRE_HPP_
/*! \file
    \brief Forward declarations for Array* classes which can be constructed
           from each other.
*/

namespace brille {
  // declare both Array and Array2 so that they can each contain conversion
  // constructors for the other.
  template<class T> class Array;
  template<class T> class Array2;

  // declare their iterators here
  template<class T> class ArrayIt;
  template<class T> class Array2It;

//  namespace lattice {
//    template<class T> class LVec;
//  }

  template<class T> using bArray = Array2<T>;
}
//template<class T> using bArray = brille::Array2<T>;

#endif