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
/*! \file
    \author Greg Tucker
    \brief Lattice vector class definitions
*/
#ifndef BRILLE_LATVEC_CLASS_H_
#define BRILLE_LATVEC_CLASS_H_

#include <typeinfo> // for std::bad_cast
#include <exception>
#include <utility>
#include "lattice_dual.hpp"
#include "tetgen.h"
#include "types.hpp"
#include "enums.hpp"
#include "array_.hpp"
// #include "array2.hpp"

//! An alias used while deciding between Array and Array2 for data storage
//template<class T>
//using bArray = brille::Array2<T>;

namespace brille::lattice {

//! Check if both vectors are the same
template<class T>
bool equal_shapes(const std::vector<T>& a, const std::vector<T>&b){
  bool ok = a.size() == b.size();
  if (ok) ok = std::equal(a.begin(), a.end(), b.begin());
  return ok;
}
//! Check if both arrays are the same
template<class T, size_t Nel>
bool equal_shapes(const std::array<T,Nel>& a, const std::array<T,Nel>&b){
  return std::equal(a.begin(), a.end(), b.begin());
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// these are all undef'd later:
#define LVEC_SCALAR_INPLACE_OP_DEF(L,X) \
L <T>&\
operator X (const T& v){\
  for (auto x: this->subItr()) this->_data[this->s2l_d(x)] X v;\
  return *this;\
}
#define LVEC_ARRAY_INPLACE_OP_DEF(L,X) \
template<template<class> class A>\
std::enable_if_t<isArray<T,A>, L <T>>&\
operator X (const A<T>& b){\
  auto itr = this->broadcastItr(b.shape());\
  if (!equal_shapes(itr.shape(), this->shape())) { \
    std::string msg="Incompatible shaped " + std::string(#L) + " arrays ";  \
    msg += list_to_string(this->shape()) + " and " + list_to_string(b.shape()); \
    msg += " for operation " + std::string(#X);\
    throw std::runtime_error(msg); \
  }\
  for (auto [ox, ax, bx]: itr) this->_data[this->s2l_d(ax)] X b[bx];\
  return *this;\
}
#define LVEC_INIT_INT(N,X) N(type_t typ, lattice_t lat, const X &n)\
: bArray<T>(static_cast<brille::ind_t>(n),3), _type(std::move(typ)), _lattice(std::move(lat))\
{}

#endif

/*! \brief 3-vector(s) expressed in units of a Direct lattice

By adding a Direct lattice to a 3-element bArray this class represents one
or more 3-vector in units of a real-space-spanning lattice.
*/
template<class T>
class LVec: public bArray<T>{
public:
  using metric_t = double; // this should probably be a template parameter
  using lattice_t = Lattice<metric_t>;
  using type_t = LengthUnit;
private:
  type_t _type = LengthUnit::inverse_angstrom;
  lattice_t _lattice;
public:
  LVec() = default;
  // default constructor for zero three-vectors:
  LVec(type_t typ, lattice_t lat): bArray<T>(0,3), _type(typ), _lattice(std::move(lat))
  {}
  // (macroed as templates can't be distinguished)
  //! integer number of three-vector constructor
  LVEC_INIT_INT(LVec, int)
  LVEC_INIT_INT(LVec, long)
  LVEC_INIT_INT(LVec, long long)
  LVEC_INIT_INT(LVec, unsigned)
  LVEC_INIT_INT(LVec, unsigned long)
  LVEC_INIT_INT(LVec, unsigned long long)
  //! Forwarding constructor to let bArray deal with everything else
  template<typename... Args> explicit LVec(type_t typ, lattice_t lat, Args... args)
  : bArray<T>(args...), _type(typ), _lattice(std::move(lat))
  {
  }
  template<class R>
  explicit LVec(const LVec<R>& other)
  : bArray<T>(other.hkl()), _type(other.type()), _lattice(other.lattice())
  {
  }

  [[nodiscard]] lattice_t lattice() const { return _lattice; }
  [[nodiscard]] type_t type() const {return _type;}
  template<class R> bool same_lattice(const LVec<R> *vec) const {return _lattice == vec->lattice() && _type == vec->type();}
  template<class R> bool star_lattice(const LVec<R> *vec) const {return _lattice == vec->lattice() && _type != vec->type();}
  template<class R> bool same_lattice(const LVec<R> &vec) const {return _lattice == vec.lattice() && _type == vec.type();}
  template<class R> bool star_lattice(const LVec<R> &vec) const {return _lattice == vec.lattice() && _type != vec.type();}

  // extract overloads preserving lattice information
  //! Return a non-copying view into the  LDVec
  template<typename... A> LVec<T> view(A... args) const;
  //! Return a copied subset of the LDVec
  template<typename... A> LVec<T> extract(A... args) const;
  //! Extract just the coordinates in units of the Direct lattice (strip off the lattice information)
  [[nodiscard]] bArray<T> hkl() const;
  //! Extract the coordinates in *an* orthonormal frame
  [[nodiscard]] bArray<double> xyz() const;
  //! Return the vector(s) expressed in units of the Reciprocal lattice
  [[nodiscard]] LVec<double> star() const;

  //! Determine the scalar product between two vectors in the object
  [[nodiscard]] double dot(ind_t i, ind_t j) const;
  //! Determine the absolute length of a vector in the object
  [[nodiscard]] double norm(ind_t i) const { return sqrt(this->dot(i,i)); }
  //! Determine the cross product of two vectors in the object
  [[nodiscard]] LVec<double> cross(ind_t i, ind_t j) const;

  LVEC_SCALAR_INPLACE_OP_DEF(LVec,+=)
  LVEC_SCALAR_INPLACE_OP_DEF(LVec,-=)
  LVEC_SCALAR_INPLACE_OP_DEF(LVec,*=)
  LVEC_SCALAR_INPLACE_OP_DEF(LVec,/=)

  LVEC_ARRAY_INPLACE_OP_DEF(LVec,+=)
  LVEC_ARRAY_INPLACE_OP_DEF(LVec,-=)
  LVEC_ARRAY_INPLACE_OP_DEF(LVec,*=)
  LVEC_ARRAY_INPLACE_OP_DEF(LVec,/=)

static LVec<T> from_std(type_t typ, lattice_t lat, const std::vector<T>& data) {
  return LVec<T>(std::move(typ), std::move(lat), bArray<T>::from_std(data));
}

template<class R, size_t N>
static LVec<T> from_std(type_t typ, lattice_t lat, const std::vector<std::array<R,N>>& data) {
  return LVec<T>(std::move(typ), std::move(lat), bArray<T>::from_std(data));
}

LVec<T> operator ^ (const T& val) {
  LVec<T> out(_type, _lattice, this->shape());
  for (auto & v: this->subItr()) out[v] = std::pow(this->val(v), val);
  return out;
}

//  template<class R>
//  void binary_operation_check(const LVec<R>& b) const{
//    assert(_lattice == b.lattice());
//  }

  //! Check whether a second LDVec is approximately the same as this object
  template<class R>
  bool is(const LVec<R>& that){
    return (same_lattice(that) && this->bArray<T>::is(that));
  }
  //! Round all elements using std::round
  [[nodiscard]] LVec<int> round() const {
    return LVec<int>(_type, _lattice, this->bArray<T>::round());
  }
  //! Find the floor of all elements using std::floor
  [[nodiscard]] LVec<int> floor() const {
    return LVec<int>(_type, _lattice, this->bArray<T>::floor());
  }
  //! Find the ceiling of all elements using std::ceil
  [[nodiscard]] LVec<int> ceil() const {
    return LVec<int>(_type, _lattice, this->bArray<T>::ceil());
  }

  LVec<T> decouple() {
    return LVec<T>(_type, _lattice, this->bArray<T>::decouple());
  }

#ifdef USE_HIGHFIVE
template<class H>
std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, bool>
 to_hdf(H& obj, const std::string& entry) const {
  auto group = overwrite_group(obj, entry);
  bool ok{true};\
  group.createAttribute("length_unit", _type);
  ok &= _lattice.to_hdf(group, "lattice");
  ok &= bArray<T>::to_hdf(group, "components");
  return ok;
}
template<class H>\
static std::enable_if_t<std::is_base_of_v<HighFive::Object,H>, LVec<T>>
from_hdf(H& obj, const std::string& entry){
  auto group = obj.getGroup(entry);
  LengthUnit lu;
  group.getAttribute("length_unit").read(lu);
  auto lat = lattice_t::from_hdf(group, "lattice");
  auto val = bArray<T>::from_hdf(group, "components");
  return LVec<T>(lu, lat, val);
}
[[nodiscard]] bool to_hdf(const std::string& filename, const std::string& dataset, unsigned perm=HighFive::File::OpenOrCreate) const {
  HighFive::File file(filename, perm);
  return this->to_hdf(file, dataset);
}
static LVec<T> from_hdf(const std::string& filename, const std::string& dataset){
  HighFive::File file(filename, HighFive::File::ReadOnly);
  return LVec<T>::from_hdf(file, dataset);
}
#endif

protected:
  void check_array();
};

#undef LVEC_SCALAR_INPLACE_OP_DEF
#undef LVEC_ARRAY_INPLACE_OP_DEF
#undef LVEC_INIT_INT


// convenience creation functions for Real and Reciprocal space vectors
template<class T, class ... P> LVec<T> LDVec(P ... args) {return LVec<T>(LengthUnit::angstrom, args ...);}
template<class T, class ... P> LVec<T> LQVec(P ... args) {return LVec<T>(LengthUnit::inverse_angstrom, args ...);}

}

#endif
