/*! \file */
#ifndef _TRANSFORM_H
#define _TRANSFORM_H

#include "latvec.h"

using fromtype = int;
using totype = double;

/*! \brief transform_to_primitive(Reciprocal, LQVec)

  Takes a Reciprocal lattice, plus a reciprocal lattice
  vector expressed in the same standard (non-primitive) lattice and returns
  an equivalent lattice vector expressed in the primitive lattice.
  The two LQVec objects have different indices but have the same Å⁻¹
  representation.
*/
template<class T,typename S=typename std::common_type<T,typename PrimitiveTransformTraits::to>::type>
LQVec<S> transform_to_primitive(const Reciprocal lat, const LQVec<T>& a){
  if (lat.primitive().issame(a.get_lattice())) return a;
  if (!lat.issame(a.get_lattice())) throw std::runtime_error("transform_to_primitive requires a common Standard lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_hall());
  if (PT.does_nothing()) return a;
  std::array<totype,9> P = PT.get_to_primitive();
  LQVec<S> out(lat.primitive(), a.size());
  for (size_t i=0; i<a.size(); ++i)
    multiply_matrix_vector<S,totype,T,3>(out.datapointer(i), P.data(), a.datapointer(i) );
  return out;
}

/*! \brief transform_from_primitive(Reciprocal, LQVec)

  Takes a Reciprocal lattice, plus a reciprocal lattice
  vector expressed in the primitive version of the lattice and returns an
  equivalent Lattice vector expressed in units of the passed lattice.
*/
template<class T,typename S=typename std::common_type<T,typename PrimitiveTransformTraits::from>::type>
LQVec<S> transform_from_primitive(const Reciprocal lat, const LQVec<T>& a){
  if (lat.issame(a.get_lattice())) return a;
  if (!lat.primitive().issame(a.get_lattice())) throw std::runtime_error("transform_from_primitive requires a common primitive lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_hall());
  if (PT.does_nothing()) return a;
  std::array<fromtype,9> P = PT.get_from_primitive();
  LQVec<S> out(lat, a.size());
  for (size_t i=0; i<a.size(); ++i)
    multiply_matrix_vector<S,fromtype,T,3>(out.datapointer(i), P.data(), a.datapointer(i) );
  return out;
}

/*! \brief transform_to_primitive(Direct, LDVec)

  Takes a Direct lattice, plus a direct lattice
  vector expressed in the same standard (non-primitive) lattice and returns
  an equivalent lattice vector expressed in the primitive lattice.
  The two lattice vectors have different indices but have the same Å⁻¹
  representation.
*/
template<class T,typename S=typename std::common_type<T,typename PrimitiveTransformTraits::from>::type>
LDVec<S> transform_to_primitive(const Direct lat, const LDVec<T>& a){
  if (lat.primitive().issame(a.get_lattice())) return a;
  if (!lat.issame(a.get_lattice())) throw std::runtime_error("transform_to_primitive requires a common Standard lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_hall());
  if (PT.does_nothing()) return a;
  std::array<fromtype,9> P = PT.get_from_primitive();
  LDVec<S> out(lat.primitive(), a.size());
  for (size_t i=0; i<a.size(); ++i)
    multiply_matrix_vector<S,fromtype,T,3>(out.datapointer(i), P.data(), a.datapointer(i) );
  return out;
}

/* \brief transform_from_primitive(Direct, LDVec)

  Takes a Direct lattice, plus a direct lattice
  vector expressed in the primitive version of the lattice and returns an
  equivalent Lattice vector expressed in units of the passed lattice.
*/
template<class T,typename S=typename std::common_type<T,typename PrimitiveTransformTraits::to>::type>
LDVec<S> transform_from_primitive(const Direct lat, const LDVec<T>& a){
  if (lat.issame(a.get_lattice())) return a;
  if (!lat.primitive().issame(a.get_lattice())) throw std::runtime_error("transform_from_primitive requires a common primitive lattice");
  // different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
  PrimitiveTransform PT(lat.get_hall());
  if (PT.does_nothing()) return a;
  std::array<totype,9> P = PT.get_to_primitive();
  LDVec<S> out(lat, a.size());
  for (size_t i=0; i<a.size(); ++i)
    multiply_matrix_vector<S,totype,T,3>(out.datapointer(i), P.data(), a.datapointer(i) );
  return out;
}


#endif
