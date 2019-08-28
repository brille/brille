/*! \file */
#ifndef _TRANSFORM_H
#define _TRANSFORM_H

#include "latvec.h"

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
LQVec<S> transform_to_primitive(const Reciprocal lat, const LQVec<T>& a){
  if (lat.primitive().issame(a.get_lattice())) return a;
  if (!lat.issame(a.get_lattice())) throw std::runtime_error("transform_to_primitive requires a common Standard lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_hall());
  if (PT.does_nothing()) return a;
  std::array<Ptype,9> P = PT.get_Pt(); // the transpose of the P matrix
  LQVec<S> out(lat.primitive(), a.size());
  for (size_t i=0; i<a.size(); ++i)
    multiply_matrix_vector<S,Ptype,T,3>(out.datapointer(i), P.data(), a.datapointer(i) );
  return out;
}

/*! \brief transform_from_primitive(Reciprocal, LQVec)

  Takes a Reciprocal lattice, plus a reciprocal lattice
  vector expressed in the primitive version of the lattice and returns an
  equivalent Lattice vector expressed in units of the passed lattice.
*/
template<class T,typename S=typename std::common_type<T,invPtype>::type>
LQVec<S> transform_from_primitive(const Reciprocal lat, const LQVec<T>& a){
  if (lat.issame(a.get_lattice())) return a;
  if (!lat.primitive().issame(a.get_lattice())) throw std::runtime_error("transform_from_primitive requires a common primitive lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_hall());
  if (PT.does_nothing()) return a;
  std::array<invPtype,9> P = PT.get_invPt(); // the inverse of the transpose of P (or the transpose of the inverse of P)
  LQVec<S> out(lat, a.size());
  for (size_t i=0; i<a.size(); ++i)
    multiply_matrix_vector<S,invPtype,T,3>(out.datapointer(i), P.data(), a.datapointer(i) );
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
LDVec<S> transform_to_primitive(const Direct lat, const LDVec<T>& a){
  if (lat.primitive().issame(a.get_lattice())) return a;
  if (!lat.issame(a.get_lattice())) throw std::runtime_error("transform_to_primitive requires a common Standard lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_hall());
  if (PT.does_nothing()) return a;
  std::array<invPtype,9> P = PT.get_invP(); // xₚ = P⁻¹ xₛ
  LDVec<S> out(lat.primitive(), a.size());
  for (size_t i=0; i<a.size(); ++i)
    multiply_matrix_vector<S,invPtype,T,3>(out.datapointer(i), P.data(), a.datapointer(i) );
  return out;
}

/* \brief transform_from_primitive(Direct, LDVec)

  Takes a Direct lattice, plus a direct lattice
  vector expressed in the primitive version of the lattice and returns an
  equivalent Lattice vector expressed in units of the passed lattice.
*/
template<class T,typename S=typename std::common_type<T,Ptype>::type>
LDVec<S> transform_from_primitive(const Direct lat, const LDVec<T>& a){
  if (lat.issame(a.get_lattice())) return a;
  if (!lat.primitive().issame(a.get_lattice())) throw std::runtime_error("transform_from_primitive requires a common primitive lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_hall());
  if (PT.does_nothing()) return a;
  std::array<Ptype,9> P = PT.get_P(); // xₛ = P xₚ
  LDVec<S> out(lat, a.size());
  for (size_t i=0; i<a.size(); ++i)
    multiply_matrix_vector<S,Ptype,T,3>(out.datapointer(i), P.data(), a.datapointer(i) );
  return out;
}


#endif
