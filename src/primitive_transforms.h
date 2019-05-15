#ifndef _PRIMITIVE_TRANSFORMS_H
#define _PRIMITIVE_TRANSFORMS_H

#include <array>
#include "linear_algebra.h"
#include "arrayvector.h"
#include "lattice.h"
#include "latvec.h"
#include "spg_database.h"
#include "primitive.h"

template<class L, class T, template<class> class V,
         typename S=typename std::common_type<T,typename PrimitiveTransformTraits::to>::type,
         typename=typename std::enable_if<std::is_base_of<LatVec,V<T>>::value>::type,      // only Lattice Vectors
				 typename=typename std::enable_if<std::is_same<L,typename LatticeTraits<V<T>>::type>::value>::type // which maintain their Lattice type
				>
V<S> transform_to_primitive(const L lat, const V<T>& a){
	if (lat.primitive().issame(a.get_lattice())) return a;
	// different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
	SpacegroupType spgt = spgdb_get_spacegroup_type(lat.get_hall());
	PrimitiveTransform PT(spgt.centering);
	if (PT.does_nothing()) return a;
	std::array<double,9> P = PT.get_to_primitive();
	V<S> out(lat.primitive(), a.size());
	for (size_t i=0; i<a.size(); ++i)
		multiply_matrix_vector<S,double,T,3>(out.datapointer(i), P.data(), a.datapointer(i) );
	return out;
}

template<class L, class T, template<class> class V,
         typename S = typename std::common_type<T,typename PrimitiveTransformTraits::from>::type,
         typename=typename std::enable_if<std::is_base_of<LatVec,V<T>>::value>::type,      // only Lattice Vectors
				 typename=typename std::enable_if<std::is_same<L,typename LatticeTraits<V<T>>::type>::value>::type // which maintain their Lattice type
				>
V<S> transform_from_primitive(const L lat, const V<T>& a){
	if (lat.issame(a.get_lattice())) return a;
	// different lattices can/should we check if the newlattice is the primitive lattice of the input lattice?
	SpacegroupType spgt = spgdb_get_spacegroup_type(lat.get_hall());
	PrimitiveTransform PT(spgt.centering);
	if (PT.does_nothing()) return a;
	std::array<int,9> P = PT.get_from_primitive();
	V<S> out(lat, a.size());
	for (size_t i=0; i<a.size(); ++i)
		multiply_matrix_vector<S,int,T,3>(out.datapointer(i), P.data(), a.datapointer(i) );
	return out;
}



template<class L, class T, template<class> class A,
         typename S = typename std::common_type<T,typename PrimitiveTransformTraits::to>::type,
				 typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type
				>
typename LatVecTraits<L,S>::type transform_to_primitive(const L lat, const A<T>& a){
	if (a.numel()!=3u) throw std::runtime_error("An ArrayVector must have numel==3 to be transformed by a lattice");
	SpacegroupType spgt = spgdb_get_spacegroup_type(lat.get_hall());
	PrimitiveTransform PT(spgt.centering);
	typename LatVecTraits<L,S>::type out(lat.primitive(), a);
	if (PT.does_anything()){
		std::array<double,9> P = PT.get_to_primitive();
		for (size_t i=0; i<a.size(); ++i)
			multiply_matrix_vector<S,double,T,3>(out.datapointer(i), P.data(), a.datapointer(i) );
	}
	return out;
}
template<class L, class T, template<class> class A,
         typename S = typename std::common_type<T,typename PrimitiveTransformTraits::from>::type,
				 typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type
				>
typename LatVecTraits<L,S>::type transform_from_primitive(const L lat, const A<T>& a){
	if (a.numel()!=3u) throw std::runtime_error("An ArrayVector must have numel==3 to be transformed by a lattice");
	SpacegroupType spgt = spgdb_get_spacegroup_type(lat.get_hall());
	PrimitiveTransform PT(spgt.centering);
	typename LatVecTraits<L,S>::type out(lat, a);
	if (PT.does_anything()){
		std::array<int,9> P = PT.get_from_primitive();
		for (size_t i=0; i<a.size(); ++i)
			multiply_matrix_vector<S,int,T,3>(out.datapointer(i), P.data(), a.datapointer(i) );
	}
	return out;
}

#endif
