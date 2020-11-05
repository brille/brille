#ifndef ARRAY__HPP
#define ARRAY__HPP

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
