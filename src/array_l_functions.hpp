#ifndef BRILLE_ARRAY_L_FUNCTIONS_HPP_
#define BRILLE_ARRAY_L_FUNCTIONS_HPP_

#include "array_l_attributes.hpp"
#include "array_l_.hpp"

namespace brille {
  template <class T, template <class> class B,
            class S = std::common_type_t<T, double>>
  std::enable_if_t<isBareArray<T, B>, lattice::LVec<S>>
  from_xyz_like(const LengthUnit lu, const lattice::Lattice<T> &lat,
                const B<T> &b) {
    auto inv_xyz = lat.from_xyz(lu);
    lattice::LVec<S> coords(lu, lat, b.shape());
    auto x = b.shape();
    x.back() = 0u;
    for (auto i : b.subItr(x))
      utils::multiply_matrix_vector(coords.ptr(i), inv_xyz.data(), b.ptr(i));
    return coords;
  }
}
#endif