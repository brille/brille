#ifndef BRILLE_ARRAY_L_ATTRIBUTES_HPP_
#define BRILLE_ARRAY_L_ATTRIBUTES_HPP_
#include "array_attributes.hpp"
#include "array_lvec.hpp"
namespace brille {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T>
  struct ArrayTraits<lattice::LVec<T>> {
  static constexpr bool array = true;
  static constexpr bool latvec = true;
};
#endif

}

#endif
