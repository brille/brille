/* Copyright 2019 Greg Tucker
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

/*! \file */
#ifndef _TRANSFORM_H
#define _TRANSFORM_H

#include "latvec.hpp"

using Ptype = PrimitiveTraits::P;
using invPtype = PrimitiveTraits::invP;

/*! \brief transform_to_primitive(Reciprocal, LQVec)

  Takes a Reciprocal lattice, plus a reciprocal lattice
  vector expressed in the same standard (non-primitive) lattice and returns
  an equivalent lattice vector expressed in the primitive lattice.
  The two LQVec objects have different indices but have the same Å⁻¹
  representation.
*/
template<class T,typename S=typename std::common_type<T,Ptype>::type>
LQVec<S> transform_to_primitive(const Reciprocal& lat, const LQVec<T>& a){
  if (lat.primitive().issame(a.get_lattice())) return a;
  if (!lat.issame(a.get_lattice())) throw std::runtime_error("transform_to_primitive requires a common Standard lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_spacegroup_object().get_bravais_type());
  if (PT.does_nothing()) return a;
  std::array<Ptype,9> P = PT.get_Pt(); // the transpose of the P matrix
  LQVec<S> out(lat.primitive(), a.size());
  for (size_t i=0; i<a.size(); ++i)
    multiply_matrix_vector<S,Ptype,T,3>(out.data(i), P.data(), a.data(i) );
  return out;
}

/*! \brief transform_from_primitive(Reciprocal, LQVec)

  Takes a Reciprocal lattice, plus a reciprocal lattice
  vector expressed in the primitive version of the lattice and returns an
  equivalent Lattice vector expressed in units of the passed lattice.
*/
template<class T,typename S=typename std::common_type<T,invPtype>::type>
LQVec<S> transform_from_primitive(const Reciprocal& lat, const LQVec<T>& a){
  if (lat.issame(a.get_lattice())) return a;
  if (!lat.primitive().issame(a.get_lattice())) throw std::runtime_error("transform_from_primitive requires a common primitive lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_spacegroup_object().get_bravais_type());
  if (PT.does_nothing()) return a;
  std::array<invPtype,9> P = PT.get_invPt(); // the inverse of the transpose of P (or the transpose of the inverse of P)
  LQVec<S> out(lat, a.size());
  for (size_t i=0; i<a.size(); ++i)
    multiply_matrix_vector<S,invPtype,T,3>(out.data(i), P.data(), a.data(i) );
  return out;
}

/*! \brief transform_to_primitive(Direct, LDVec)

  Takes a Direct lattice, plus a direct lattice
  vector expressed in the same standard (non-primitive) lattice and returns
  an equivalent lattice vector expressed in the primitive lattice.
  The two lattice vectors have different indices but have the same Å⁻¹
  representation.
*/
template<class T,typename S=typename std::common_type<T,invPtype>::type>
LDVec<S> transform_to_primitive(const Direct& lat, const LDVec<T>& a){
  if (lat.primitive().issame(a.get_lattice())) return a;
  if (!lat.issame(a.get_lattice())) throw std::runtime_error("transform_to_primitive requires a common Standard lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_spacegroup_object().get_bravais_type());
  if (PT.does_nothing()) return a;
  std::array<invPtype,9> P = PT.get_invP(); // xₚ = P⁻¹ xₛ
  LDVec<S> out(lat.primitive(), a.size());
  for (size_t i=0; i<a.size(); ++i)
    multiply_matrix_vector<S,invPtype,T,3>(out.data(i), P.data(), a.data(i) );
  return out;
}

/* \brief transform_from_primitive(Direct, LDVec)

  Takes a Direct lattice, plus a direct lattice
  vector expressed in the primitive version of the lattice and returns an
  equivalent Lattice vector expressed in units of the passed lattice.
*/
template<class T,typename S=typename std::common_type<T,Ptype>::type>
LDVec<S> transform_from_primitive(const Direct& lat, const LDVec<T>& a){
  if (lat.issame(a.get_lattice())) return a;
  if (!lat.primitive().issame(a.get_lattice())) throw std::runtime_error("transform_from_primitive requires a common primitive lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_spacegroup_object().get_bravais_type());
  if (PT.does_nothing()) return a;
  std::array<Ptype,9> P = PT.get_P(); // xₛ = P xₚ
  LDVec<S> out(lat, a.size());
  for (size_t i=0; i<a.size(); ++i)
    multiply_matrix_vector<S,Ptype,T,3>(out.data(i), P.data(), a.data(i) );
  return out;
}

// utility functions for conversion of lattice vectors where only their components are stored
template<class T,typename S=typename std::common_type<T,double>::type>
ArrayVector<S> xyz_to_hkl(const Reciprocal& lat, const ArrayVector<T>& xyz){
  // this error could be relaxed by changing to hkl(xyz.numel(), xyz.size()) below
  if (xyz.numel()!=3u) throw std::runtime_error("coordinate transformation only defined for three-vectors");
  double toxyz[9], fromxyz[9];
  lat.get_xyz_transform(toxyz);
  if (!matrix_inverse(fromxyz,toxyz)) throw std::runtime_error("transform matrix toxyz has zero determinant");
  ArrayVector<S> hkl(3, xyz.size());
  for (size_t i=0; i<xyz.size(); ++i) multiply_matrix_vector<S,double,T,3>(hkl.data(i), fromxyz, xyz.data(i));
  return hkl;
}
template<class T,typename S=typename std::common_type<T,double>::type>
ArrayVector<S> hkl_to_xyz(const Reciprocal& lat, const ArrayVector<T>& hkl){
  // this error could be relaxed by changing to xyz(hkl.numel(), hkl.size()) below
  if (hkl.numel()!=3u) throw std::runtime_error("coordinate transformation only defined for three-vectors");
  double toxyz[9];
  lat.get_xyz_transform(toxyz);
  ArrayVector<S> xyz(3, hkl.size());
  for (size_t i=0; i<hkl.size(); ++i) multiply_matrix_vector<S,double,T,3>(xyz.data(i), toxyz, hkl.data(i));
  return xyz;
}


#endif
