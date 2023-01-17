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
/*! \file
    \author Greg Tucker
    \brief Functions to transform Direct and Reciprocal lattice vectors between
           primitive and conventional unit cell representations
*/
#include "array_.hpp" // defines bArray
#include "primitive.hpp"
#include "array_l_.hpp"
#include "omp.h"
namespace brille::lattice {

//! The datatype of the P PrimitiveTransform matrix
using sixPtype = PrimitiveTraits::sixP;
//! The datatype of the inverse P PrimitiveTransform matrix
using invPtype = PrimitiveTraits::invP;

/*! \brief transform_to_primitive(Reciprocal, LQVec)

  Takes a Reciprocal lattice, plus a reciprocal lattice
  vector expressed in the same standard (non-primitive) lattice and returns
  an equivalent lattice vector expressed in the primitive lattice.
  The two LQVec objects have different indices but have the same Å⁻¹
  representation.
*/
template<class T, typename S=typename std::common_type_t<T,int>>
LVec<S> transform_to_primitive(const Lattice<double>& lat, const LVec<T>& a){
  if (lat.primitive().is_same(a.lattice())) return a;
  if (!lat.is_same(a.lattice()))
    throw std::runtime_error("transform_to_primitive requires a common Standard lattice");
  // different lattices can/should we check if the new lattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.bravais());
  if (PT.does_nothing()) return LVec<S>(a);
  assert(a.stride().back() == 1u && a.size(a.ndim()-1)==3);
  LengthUnit lu = a.type();
  std::array<int, 9> transform{0,0,0,0,0,0,0,0,0};
  switch (lu) {
  case LengthUnit::angstrom:
    // For a column vector matrix of basis vectors in the conventional lattice A
    // The tranformation matrix P gives the primitive basis vectors Aₚ as
    //      A P = Aₚ
    // A point in space can be represented in units of either basis, with the
    // conventional representation, x, and the primitive representation, xₚ.
    // If these points are the same, then A x == Aₚ xₚ.
    //      A x = Aₚ  xₚ
    //      A x = A P xₚ
    //        x =   P xₚ
    //    P⁻¹ x =     xₚ
    transform = PT.get_invP();
    break;
  case LengthUnit::inverse_angstrom:
    // For a row vector matrix of reciprocal basis vectors in the conventional
    // lattice, B, the same P transformation matrix relates them to the
    // primitive lattice basis vectors Bₚ
    //      P⁻¹ B = Bₚ
    // And a point in space can be represented, again, in either basis with the
    // conventional representation, q, and the primitive representation, qₚ.
    // With both representing the same point in space Bᵀ q = Bₚᵀ qₚ
    //      Bᵀ q = Bₚᵀ qₚ
    //      Bᵀ q = (P⁻¹ B)ᵀ qₚ
    //      Bᵀ q = Bᵀ(P⁻¹)ᵀ qₚ
    //         q =   (P⁻¹)ᵀ qₚ
    //      Pᵀ q = Pᵀ(P⁻¹)ᵀ qₚ
    //      Pᵀ q = (P⁻¹P)ᵀ qₚ
    //      Pᵀ q = qₚ
    transform = PT.get_6Pt();
    break;
  default:
    throw std::runtime_error("Not implemented");
  }
  auto sh = a.shape();
  LVec<S> out(lu, lat.primitive(), sh);
  sh.back() = 0;
  for (auto x: a.subItr(sh))
    brille::utils::multiply_matrix_vector(out.ptr(x), transform.data(), a.ptr(x));
  if (LengthUnit::inverse_angstrom == lu){
    out /= S(6); // correct for having used 6 * Pt
  }
  return out;
}

/*! \brief transform_from_primitive(Reciprocal, LQVec)

  Takes a Reciprocal lattice, plus a reciprocal lattice
  vector expressed in the primitive version of the lattice and returns an
  equivalent Lattice vector expressed in units of the passed lattice.
*/
template<class T, typename S=typename std::common_type_t<T,int>>
LVec<S> transform_from_primitive(const Lattice<double>& lat, const LVec<T>& a){
  if (lat.is_same(a.lattice())) return a;
  if (!lat.primitive().is_same(a.lattice()))
    throw std::runtime_error("transform_from_primitive requires a common primitive lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.bravais());
  if (PT.does_nothing()) return LVec<S>(a);
  assert(a.stride().back() == 1u && a.size(a.ndim()-1)==3);
  LengthUnit lu = a.type();
  std::array<int, 9> transform{0,0,0,0,0,0,0,0,0};
  switch (lu) {
  case LengthUnit::angstrom:
    // For a column vector matrix of basis vectors in the conventional lattice A
    // The tranformation matrix P gives the primitive basis vectors Aₚ as
    //      A P = Aₚ
    // A point in space can be represented in units of either basis, with the
    // conventional representation, x, and the primitive representation, xₚ.
    // If these points are the same, then A x == Aₚ xₚ.
    //      A x = Aₚ  xₚ
    //      A x = A P xₚ
    //        x =   P xₚ
    transform = PT.get_6P();
    break;
  case LengthUnit::inverse_angstrom:
    // For a row vector matrix of reciprocal basis vectors in the conventional
    // lattice, B, the same P transformation matrix relates them to the
    // primitive lattice basis vectors Bₚ
    //      P⁻¹ B = Bₚ
    // And a point in space can be represented, again, in either basis with the
    // conventional representation, q, and the primitive representation, qₚ.
    // With both representing the same point in space Bᵀ q = Bₚᵀ qₚ
    //      Bᵀ q = Bₚᵀ qₚ
    //      Bᵀ q = (P⁻¹ B)ᵀ qₚ
    //      Bᵀ q = Bᵀ(P⁻¹)ᵀ qₚ
    //         q =   (P⁻¹)ᵀ qₚ
    transform = PT.get_invPt();
    break;
  default:
    throw std::runtime_error("Not implemented");
  }
  auto sh = a.shape();
  LVec<S> out(lu, lat, sh);
  sh.back() = 0;
  for (auto x: a.subItr(sh))
    brille::utils::multiply_matrix_vector(out.ptr(x), transform.data(), a.ptr(x));
  if (LengthUnit::angstrom == lu){
    out /= S(6); // since we used 6P instead of P
  }
  return out;
}

template<class T, typename S=typename std::common_type_t<T,int>>
LVec<S> parallel_transform_to_primitive(const Lattice<double>& lat, const LVec<T>& a, int threads){
  if (lat.primitive().is_same(a.lattice())) return a;
  if (!lat.is_same(a.lattice()))
    throw std::runtime_error("transform_to_primitive requires a common Standard lattice");
  // different lattices can/should we check if the new lattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.bravais());
  if (PT.does_nothing()) return LVec<S>(a);
  assert(a.stride().back() == 1u && a.size(a.ndim()-1)==3);
  LengthUnit lu = a.type();
  std::array<int, 9> transform{0,0,0,0,0,0,0,0,0};
  switch (lu) {
  case LengthUnit::angstrom: transform = PT.get_invP(); break;
  case LengthUnit::inverse_angstrom: transform = PT.get_6Pt(); break;
  default: throw std::runtime_error("Not implemented");
  }
  auto sh = a.shape();
  LVec<S> out(lu, lat.primitive(), sh);
  auto na = static_cast<int64_t>(a.size(0));
  auto t_ptr = transform.data();
  if (threads < 1) threads = omp_get_max_threads();
  omp_set_num_threads(threads);
#pragma omp parallel for default(none) shared(out, t_ptr, a, na)
  for (int64_t si=0; si<na; ++si){
    auto i = static_cast<ind_t>(si);
    brille::utils::multiply_matrix_vector(out.ptr(i), t_ptr, a.ptr(i));
  }

  if (LengthUnit::inverse_angstrom == lu){
    out /= S(6); // correct for having used 6 * Pt
  }
  return out;
}

/*! \brief transform_from_primitive(Reciprocal, LQVec)

  Takes a Reciprocal lattice, plus a reciprocal lattice
  vector expressed in the primitive version of the lattice and returns an
  equivalent Lattice vector expressed in units of the passed lattice.
*/
template<class T, typename S=typename std::common_type_t<T,int>>
LVec<S> parallel_transform_from_primitive(const Lattice<double>& lat, const LVec<T>& a, int threads){
  if (lat.is_same(a.lattice())) return a;
  if (!lat.primitive().is_same(a.lattice()))
    throw std::runtime_error("transform_from_primitive requires a common primitive lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.bravais());
  if (PT.does_nothing()) return LVec<S>(a);
  assert(a.stride().back() == 1u && a.size(a.ndim()-1)==3);
  LengthUnit lu = a.type();
  std::array<int, 9> transform{0,0,0,0,0,0,0,0,0};
  switch (lu) {
  case LengthUnit::angstrom: transform = PT.get_6P(); break;
  case LengthUnit::inverse_angstrom: transform = PT.get_invPt(); break;
  default: throw std::runtime_error("Not implemented");
  }
  auto sh = a.shape();
  LVec<S> out(lu, lat, sh);
  auto na = static_cast<int64_t>(a.size(0));
  auto t_ptr = transform.data();
  if (threads < 1) threads = omp_get_max_threads();
  omp_set_num_threads(threads);
#pragma omp parallel for default(none) shared(out, t_ptr, a, na)
  for (int64_t si=0; si<na; ++si){
    auto i = static_cast<ind_t>(si);
    brille::utils::multiply_matrix_vector(out.ptr(i), t_ptr, a.ptr(i));
  }

  if (LengthUnit::angstrom == lu){
    out /= S(6); // since we used 6P instead of P
  }
  return out;
}

//! utility functions for conversion of lattice vectors where only their components are stored
template<class T, typename S=typename std::common_type_t<T,double>>
bArray<S> xyz_to_hkl(const Lattice<double>& lat, const LengthUnit lu, const bArray<T>& xyz){
  assert(xyz.stride().back() == 1u && xyz.size(xyz.ndim()-1)==3);
  typename Lattice<double>::matrix_t tmp, inv;
  switch (lu) {
  case LengthUnit::angstrom: tmp = lat.real_basis_vectors(); break;
  case LengthUnit::inverse_angstrom: tmp = lat.reciprocal_basis_vectors(); break;
  default: throw std::logic_error("Not implemented");
  }
  utils::matrix_inverse(inv.data(), tmp.data());

  auto sh = xyz.shape();
  bArray<S> hkl(sh);
  sh.back() = 0;
  for (auto x: xyz.subItr(sh))
    brille::utils::multiply_matrix_vector(hkl.ptr(x), inv.data(), xyz.ptr(x));
  return hkl;
}
//! utility functions for conversion of lattice vectors where only their components are stored
template<class T, typename S=typename std::common_type<T,double>::type>
bArray<S> hkl_to_xyz(const Lattice<double>& lat, const LengthUnit lu, const bArray<T>& hkl){
  assert(hkl.stride().back() == 1u && hkl.size(hkl.ndim()-1)==3);
  typename Lattice<double>::matrix_t tmp;
  switch (lu) {
  case LengthUnit::angstrom: tmp = lat.real_basis_vectors(); break;
  case LengthUnit::inverse_angstrom: tmp = lat.reciprocal_basis_vectors(); break;
  default: throw std::logic_error("Not implemented");
  }

  auto sh = hkl.shape();
  bArray<S> xyz(sh);
  sh.back() = 0;
  for (auto x: hkl.subItr(sh))
    brille::utils::multiply_matrix_vector(xyz.ptr(x), tmp.data(), hkl.ptr(x));
  return xyz;
}

} // end namespace brille
#endif
