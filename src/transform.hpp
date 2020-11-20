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

/*! \file */

#ifndef BRILLE_TRANSFORM_HPP_
#define BRILLE_TRANSFORM_HPP_
#include "array_latvec.hpp" // defines bArray
namespace brille {

using Ptype = PrimitiveTraits::P;
using invPtype = PrimitiveTraits::invP;

/*! \brief transform_to_primitive(Reciprocal, LQVec)

  Takes a Reciprocal lattice, plus a reciprocal lattice
  vector expressed in the same standard (non-primitive) lattice and returns
  an equivalent lattice vector expressed in the primitive lattice.
  The two LQVec objects have different indices but have the same Å⁻¹
  representation.
*/
template<class T, typename S=typename std::common_type<T,Ptype>::type>
LQVec<S> transform_to_primitive(const Reciprocal& lat, const LQVec<T>& a){
  if (lat.primitive().issame(a.get_lattice())) return a;
  if (!lat.issame(a.get_lattice()))
    throw std::runtime_error("transform_to_primitive requires a common Standard lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_bravais_type());
  if (PT.does_nothing()) return LQVec<S>(a);
  assert(a.stride().back() == 1u && a.size(a.ndim()-1)==3);
  std::array<Ptype,9> Pmat = PT.get_Pt(); // the transpose of the P matrix
  auto sh = a.shape();
  LQVec<S> out(lat.primitive(), sh);
  sh.back() = 0;
  for (auto x: a.subItr(sh))
    brille::utils::multiply_matrix_vector<S,Ptype,T,3>(out.ptr(x), Pmat.data(), a.ptr(x));
  return out;
}

/*! \brief transform_from_primitive(Reciprocal, LQVec)

  Takes a Reciprocal lattice, plus a reciprocal lattice
  vector expressed in the primitive version of the lattice and returns an
  equivalent Lattice vector expressed in units of the passed lattice.
*/
template<class T, typename S=typename std::common_type<T,invPtype>::type>
LQVec<S> transform_from_primitive(const Reciprocal& lat, const LQVec<T>& a){
  if (lat.issame(a.get_lattice())) return a;
  if (!lat.primitive().issame(a.get_lattice()))
    throw std::runtime_error("transform_from_primitive requires a common primitive lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_bravais_type());
  if (PT.does_nothing()) return LQVec<S>(a);
  assert(a.stride().back() == 1u && a.size(a.ndim()-1)==3);
  std::array<invPtype,9> Pmat = PT.get_invPt(); // the inverse of the transpose of P (or the transpose of the inverse of P)
  auto sh = a.shape();
  LQVec<S> out(lat, sh);
  sh.back() = 0;
  for (auto x: a.subItr(sh))
    brille::utils::multiply_matrix_vector<S,invPtype,T,3>(out.ptr(x), Pmat.data(), a.ptr(x));
  return out;
}

/*! \brief transform_to_primitive(Direct, LDVec)

  Takes a Direct lattice, plus a direct lattice
  vector expressed in the same standard (non-primitive) lattice and returns
  an equivalent lattice vector expressed in the primitive lattice.
  The two lattice vectors have different indices but have the same Å⁻¹
  representation.
*/
template<class T, typename S=typename std::common_type<T,invPtype>::type>
LDVec<S> transform_to_primitive(const Direct& lat, const LDVec<T>& a){
  if (lat.primitive().issame(a.get_lattice())) return a;
  if (!lat.issame(a.get_lattice()))
    throw std::runtime_error("transform_to_primitive requires a common Standard lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_bravais_type());
  if (PT.does_nothing()) return LDVec<S>(a);
  assert(a.stride().back() == 1u && a.size(a.ndim()-1)==3);
  std::array<invPtype,9> Pmat = PT.get_invP(); // xₚ = P⁻¹ xₛ
  auto sh = a.shape();
  LDVec<S> out(lat.primitive(), sh);
  sh.back() = 0;
  for (auto x: a.subItr(sh))
    brille::utils::multiply_matrix_vector<S,invPtype,T,3>(out.ptr(x), Pmat.data(), a.ptr(x));
  return out;
}

/* \brief transform_from_primitive(Direct, LDVec)

  Takes a Direct lattice, plus a direct lattice
  vector expressed in the primitive version of the lattice and returns an
  equivalent Lattice vector expressed in units of the passed lattice.
*/
template<class T, typename S=typename std::common_type<T,Ptype>::type>
LDVec<S> transform_from_primitive(const Direct& lat, const LDVec<T>& a){
  if (lat.issame(a.get_lattice())) return a;
  if (!lat.primitive().issame(a.get_lattice()))
    throw std::runtime_error("transform_from_primitive requires a common primitive lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_bravais_type());
  if (PT.does_nothing()) return LDVec<S>(a);
  assert(a.stride().back() == 1u && a.size(a.ndim()-1)==3);
  std::array<Ptype,9> Pmat = PT.get_P(); // xₛ = P xₚ
  auto sh = a.shape();
  LDVec<S> out(lat, sh);
  sh.back() = 0;
  for (auto x: a.subItr(sh))
    brille::utils::multiply_matrix_vector<S,Ptype,T,3>(out.ptr(x), Pmat.data(), a.ptr(x));
  return out;
}

// utility functions for conversion of lattice vectors where only their components are stored
template<class T, typename S=typename std::common_type<T,double>::type>
bArray<S> xyz_to_hkl(const Reciprocal& lat, const bArray<T>& xyz){
  assert(xyz.stride().back() == 1u && xyz.size(xyz.ndim()-1)==3);
  auto fromxyz = lat.get_inverse_xyz_transform();
  auto sh = xyz.shape();
  bArray<S> hkl(sh);
  sh.back() = 0;
  for (auto x: xyz.subItr(sh))
    brille::utils::multiply_matrix_vector<S,double,T,3>(hkl.ptr(x), fromxyz.data(), xyz.ptr(x));
  return hkl;
}
template<class T, typename S=typename std::common_type<T,double>::type>
bArray<S> hkl_to_xyz(const Reciprocal& lat, const bArray<T>& hkl){
  assert(hkl.stride().back() == 1u && hkl.size(hkl.ndim()-1)==3);
  auto toxyz = lat.get_xyz_transform();
  auto sh = hkl.shape();
  bArray<S> xyz(sh);
  sh.back() = 0;
  for (auto x: hkl.subItr(sh))
    brille::utils::multiply_matrix_vector<S,double,T,3>(xyz.ptr(x), toxyz.data(), hkl.ptr(x));
  return xyz;
}

} // end namespace brille
#endif
