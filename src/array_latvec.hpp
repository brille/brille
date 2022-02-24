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
#include "lattice.hpp"
#include "tetgen.h"
// #include "array2.hpp"

//! An alias used while deciding between Array and Array2 for data storage
template<class T>
using bArray = brille::Array2<T>;

namespace brille {

/*! \brief Superclass to identify both LDVec and LQVec

The two Lattice vector classes, LDVec and LQVec, are subclasses of bArray.
In order to distinguish between the Lattice vector types and "bare" bArray
objects in operator- and function-overloading it is advantageous to have a
shared superclass on which to enable templates. Thus LatVec is a superclass used
for logic only which has no properties and defines no methods.
*/
class LatVec{};
template<class T> class LDVec;
template<class T> class LQVec;

/*! \brief Templated struct to help differentiate between Array2 and its derived
           classes LQVec and LDVec

Since the derived classes add Lattice information some functions behave
differently for 'raw' Array2 versus LQVec or LDVec objects.
This struct is used in template arguments to cause substitution failure-based
differentiation of otherwise identical function signatures.
*/
template<class... T> struct ArrayTraits{
  static constexpr bool array = false;  //!< is this an Array2, LQVec, or LDVec
  static constexpr bool latvec = false; //!< is this an LQVec or LDVec
};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T> struct ArrayTraits<bArray<T>>{
  static constexpr bool array = true;
  static constexpr bool latvec = false;
};
template<class T> struct ArrayTraits<LDVec<T>>{
  static constexpr bool array = true;
  static constexpr bool latvec = true;
};
template<class T> struct ArrayTraits<LQVec<T>>{
  static constexpr bool array = true;
  static constexpr bool latvec = true;
};
#endif

/*! Easier access to ArrayTraits for template substitution

true for Array2, LQVec, LDVec
*/
template<class T, template<class> class A>
inline constexpr bool isArray = ArrayTraits<A<T>>::array;
/*! Easier access to ArrayTraits for template substitution

true for LQVec, LDVec
*/
template<class T, template<class> class A>
inline constexpr bool isLatVec = ArrayTraits<A<T>>::latvec;
/*! Easier access to ArrayTraits for template substitution

true for Array2
*/
template<class T, template<class> class A>
inline constexpr bool isBareArray = isArray<T,A> && !isLatVec<T,A>;

/*! Easier access to ArrayTraits for double template substitution

true for any combination of two Array2, LQVec, LDVec objects
*/
template<class T, template<class> class A, class R, template<class> class B>
inline constexpr bool bothArrays = isArray<T,A> && isArray<R,B>;
/*! Easier access to ArrayTraits for double template substitution

true for two Array2 objects
*/
template<class T, template<class> class A, class R, template<class> class B>
inline constexpr bool bareArrays = isBareArray<T,A> && isBareArray<R,B>;
/*! Easier access to ArrayTraits for double template substitution

true for any combination of two LQVec, LDVec objects
*/
template<class T, template<class> class A, class R, template<class> class B>
inline constexpr bool bothLatVecs = isLatVec<T,A> && isLatVec<R,B>;

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
#define LVEC_INIT_INT(N,L,X) N(L lat, const X &n)\
: bArray<T>(static_cast<brille::ind_t>(n),3), lattice(std::move(lat))\
{}

#define LVEC_FROM_STD(A, L, CONTAINER)\
static A<T> from_std(L lat, const CONTAINER& data) {\
  return A<T>(std::move(lat), bArray<T>::from_std(data));\
}

#define LVEC_FROM_STD_VECTOR_ARRAY(A, L) \
template<class R, size_t N>\
static A<T> from_std(L lat, const std::vector<std::array<R,N>>& data) { \
  return A<T>(std::move(lat), bArray<T>::from_std(data));\
}

#define LVEC_SCALAR_POWER(L) \
L<T> operator ^ (const T& val) {\
  L<T> out(this->get_lattice(), this->shape());\
  for (auto & v: this->subItr()) out[v] = std::pow(this->val(v), val);\
  return out;\
}


#ifdef USE_HIGHFIVE

#define TO_HDF_OBJ(INDEX_NAMES)\
template<class H>\
std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, bool>\
 to_hdf(H& obj, const std::string& entry) const {\
  auto group = overwrite_group(obj, entry);\
  bool ok{true};\
  ok &= lattice.to_hdf(group, "lattice");\
  ok &= bArray<T>::to_hdf(group, "INDEX_NAMES");\
  return ok;\
}
#define FROM_HDF_OBJ(L, LATTICE_TYPE, INDEX_NAMES)\
template<class H>\
static std::enable_if_t<std::is_base_of_v<HighFive::Object,H>, L<T>>\
from_hdf(H& obj, const std::string& entry){\
  auto group = obj.getGroup(entry);\
  auto lat = LATTICE_TYPE::from_hdf(group, "lattice");\
  auto val = bArray<T>::from_hdf(group, "INDEX_NAMES");\
  return L<T>(lat, val);\
}

#define TO_HDF_STR \
[[nodiscard]] bool to_hdf(const std::string& filename, const std::string& dataset, unsigned perm=HighFive::File::OpenOrCreate) const {\
  HighFive::File file(filename, perm);\
  return this->to_hdf(file, dataset);\
}

#define FROM_HDF_STR(L)\
static L<T> from_hdf(const std::string& filename, const std::string& dataset){\
  HighFive::File file(filename, HighFive::File::ReadOnly);\
  return L<T>::from_hdf(file, dataset);\
}

#endif

#endif

/*! \brief 3-vector(s) expressed in units of a Direct lattice

By adding a Direct lattice to a 3-element bArray this class represents one
or more 3-vector in units of a real-space-spanning lattice.
*/
template<class T>
class LDVec: public LatVec, public bArray<T>{
  Direct lattice;
public:
  // default constructor for zero three-vectors:
  explicit LDVec(Direct lat=Direct())
  : bArray<T>(0,3), lattice(std::move(lat))
  {
  }
  // (macroed as templates can't be distinguished)
  //! integer number of three-vector constructor
  LVEC_INIT_INT(LDVec,Direct,int)
  LVEC_INIT_INT(LDVec,Direct,long)
  LVEC_INIT_INT(LDVec,Direct,long long)
  LVEC_INIT_INT(LDVec,Direct,unsigned)
  LVEC_INIT_INT(LDVec,Direct,unsigned long)
  LVEC_INIT_INT(LDVec,Direct,unsigned long long)
  //! Forwarding constructor to let bArray deal with everything else
  template<typename... Args> explicit LDVec(Direct lat, Args... args)
  : bArray<T>(args...), lattice(std::move(lat))
  {
    this->check_array();
  }
  template<class R>
  explicit LDVec(const LDVec<R>& other)
  : bArray<T>(other.get_hkl()), lattice(other.get_lattice())
  {
  }

  [[nodiscard]] Direct get_lattice() const { return lattice; }
  template<class R> bool samelattice(const LDVec<R> *vec) const { return lattice.issame(vec->get_lattice()); }
  template<class R> bool samelattice(const LQVec<R> *)    const { return false; }
  template<class R> bool starlattice(const LDVec<R> *)    const { return false; }
  template<class R> bool starlattice(const LQVec<R> *vec) const { return lattice.isstar(vec->get_lattice()); }
  template<class R> bool samelattice(const LDVec<R> &vec) const { return lattice.issame(vec.get_lattice()); }
  template<class R> bool samelattice(const LQVec<R> &)    const { return false; }
  template<class R> bool starlattice(const LDVec<R> &)    const { return false; }
  template<class R> bool starlattice(const LQVec<R> &vec) const { return lattice.isstar(vec.get_lattice()); }
  // extract overloads preserving lattice information
  //! Return a non-copying view into the  LDVec
  template<typename... A> LDVec<T> view(A... args) const;
  //! Return a copied subset of the LDVec
  template<typename... A> LDVec<T> extract(A... args) const;
  //! Extract just the coordinates in units of the Direct lattice (strip off the lattice information)
  [[nodiscard]] bArray<T> get_hkl() const;
  //! Extract the coordinates in *an* orthonormal frame
  [[nodiscard]] bArray<double> get_xyz() const;
  //! Return the vector(s) expressed in units of the Reciprocal lattice
  [[nodiscard]] LQVec<double> star() const;

  //! Determine the scalar product between two vectors in the object
  [[nodiscard]] double dot(ind_t i, ind_t j) const;
  //! Determine the absolute length of a vector in the object
  [[nodiscard]] double norm(ind_t i) const { return sqrt(this->dot(i,i)); }
  //! Determine the cross product of two vectors in the object
  [[nodiscard]] LDVec<double> cross(ind_t i, ind_t j) const;

  // LDVec<T>& operator-=(const LDVec<T>& av);
  // LDVec<T>& operator+=(const LDVec<T>& av);

  LVEC_SCALAR_INPLACE_OP_DEF(LDVec,+=)
  LVEC_SCALAR_INPLACE_OP_DEF(LDVec,-=)
  LVEC_SCALAR_INPLACE_OP_DEF(LDVec,*=)
  LVEC_SCALAR_INPLACE_OP_DEF(LDVec,/=)

  LVEC_ARRAY_INPLACE_OP_DEF(LDVec,+=)
  LVEC_ARRAY_INPLACE_OP_DEF(LDVec,-=)
  LVEC_ARRAY_INPLACE_OP_DEF(LDVec,*=)
  LVEC_ARRAY_INPLACE_OP_DEF(LDVec,/=)

  LVEC_FROM_STD(LDVec, Direct, std::vector<T>)
  LVEC_FROM_STD_VECTOR_ARRAY(LDVec, Direct)

  LVEC_SCALAR_POWER(LDVec)

  template<class R>
  void binary_operation_check(const LDVec<R>& b) const{
    assert(this->samelattice(b));
  }
  template<class R>
  void binary_operation_check(const LQVec<R>& b) const{
    assert(this->starlattice(b));
  }
  template<class R, template<class> class A,
  typename=typename std::enable_if<!std::is_base_of<LatVec,A<R>>::value>::type>
  void binary_operation_check(const A<R>&) const {}

  //! Check whether a second LDVec is approximately the same as this object
  template<class R>
  bool is(const LDVec<R>& that){
    return (this->samelattice(that) && this->bArray<T>::is(that));
  }
  //! Round all elements using std::round
  [[nodiscard]] LDVec<int> round() const {
    return LDVec<int>(this->lattice, this->bArray<T>::round());
  }
  //! Find the floor of all elements using std::floor
  [[nodiscard]] LDVec<int> floor() const {
    return LDVec<int>(this->lattice, this->bArray<T>::floor());
  }
  //! Find the ceiling of all elements using std::ceil
  [[nodiscard]] LDVec<int> ceil() const {
    return LDVec<int>(this->lattice, this->bArray<T>::ceil());
  }

  LDVec<T> decouple() {
    return LDVec<T>(lattice, this->bArray<T>::decouple());
  }

  TO_HDF_OBJ(xyz)
  TO_HDF_STR
  FROM_HDF_OBJ(LDVec, Direct, xyz)
  FROM_HDF_STR(LDVec)
protected:
  void check_array();
};


/*! \brief 3-vector(s) expressed in units of a Reciprocal lattice

By adding a Reciprocal lattice to a 3-element bArray this class represents one
or more 3-vector in units of a reciprocal-space-spanning lattice.
*/
template<class T>
class LQVec:  public LatVec, public bArray<T>{
  Reciprocal lattice;
public:
  // default constructor for zero three-vectors:
  explicit LQVec(Reciprocal  lat=Reciprocal())
  : bArray<T>(0,3), lattice(std::move(lat))
  {
  }
  //! integer number of three-vector constructor (macroed as templates can't be distinguished)
  LVEC_INIT_INT(LQVec,Reciprocal,int)
  LVEC_INIT_INT(LQVec,Reciprocal,long)
  LVEC_INIT_INT(LQVec,Reciprocal,long long)
  LVEC_INIT_INT(LQVec,Reciprocal,unsigned)
  LVEC_INIT_INT(LQVec,Reciprocal,unsigned long)
  LVEC_INIT_INT(LQVec,Reciprocal,unsigned long long)
  //! Fowarding constructor to let bArray deal with everything else
  template<typename... Args>
  explicit LQVec(Reciprocal  lat, Args... args)
  : bArray<T>(args...), lattice(std::move(lat))
  {
    this->check_array();
  }

  template<class R>
  explicit LQVec(const LQVec<R>& other)
  : bArray<T>(other.get_hkl()), lattice(other.get_lattice())
  {
  }

  [[nodiscard]] Reciprocal get_lattice() const { return lattice; }
  template<class R> bool samelattice(const LQVec<R> *vec) const { return lattice.issame(vec->get_lattice()); }
  template<class R> bool samelattice(const LDVec<R> *)    const { return false; }
  template<class R> bool starlattice(const LQVec<R> *)    const { return false; }
  template<class R> bool starlattice(const LDVec<R> *vec) const { return lattice.isstar(vec->get_lattice()); }
  template<class R> bool samelattice(const LQVec<R> &vec) const { return lattice.issame(vec.get_lattice()); }
  template<class R> bool samelattice(const LDVec<R> &)    const { return false; }
  template<class R> bool starlattice(const LQVec<R> &)    const { return false; }
  template<class R> bool starlattice(const LDVec<R> &vec) const { return lattice.isstar(vec.get_lattice()); }
  // extract overloads preserving lattice information
  //! Return a non-copying view into the  LQVec
  template<typename... A> LQVec<T> view(A... args) const;
  //! Return a copied subset of the LQVec
  template<typename... A> LQVec<T> extract(A... args) const;
  //! Extract just the coordinates in units of the Reciprocal lattice (strip off the lattice information)
  [[nodiscard]] bArray<T> get_hkl() const;
  /*! Extract the coordinates in an orthonormal frame with its first axis, x,
  along a*, its second, y, perpendicular with y⋅b*>0 , and it's third forming
  the right-handed set z=x×y.
  */
  [[nodiscard]] bArray<double> get_xyz() const;
  //! Return the vector(s) expressed in units of the Direct lattice
  [[nodiscard]] LDVec<double> star() const;

  //! Determine the scalar product between two vectors in the object.
  [[nodiscard]] double dot(ind_t i, ind_t j) const;
  //! Determine the absolute length of a vector in the object.
  [[nodiscard]] double norm(ind_t i) const { return sqrt(this->dot(i,i)); }
  //! Determine the cross product of two vectors in the object.
  [[nodiscard]] LQVec<double> cross(ind_t i, ind_t j) const;

  LVEC_SCALAR_INPLACE_OP_DEF(LQVec,+=)
  LVEC_SCALAR_INPLACE_OP_DEF(LQVec,-=)
  LVEC_SCALAR_INPLACE_OP_DEF(LQVec,*=)
  LVEC_SCALAR_INPLACE_OP_DEF(LQVec,/=)

  LVEC_ARRAY_INPLACE_OP_DEF(LQVec,+=)
  LVEC_ARRAY_INPLACE_OP_DEF(LQVec,-=)
  LVEC_ARRAY_INPLACE_OP_DEF(LQVec,*=)
  LVEC_ARRAY_INPLACE_OP_DEF(LQVec,/=)
  // LQVec<T> operator -();

  LVEC_FROM_STD(LQVec, Reciprocal, std::vector<T>)
  LVEC_FROM_STD_VECTOR_ARRAY(LQVec, Reciprocal)

  LVEC_SCALAR_POWER(LQVec)

  template<class R>
  void binary_operation_check(const LQVec<R>& b) const{
    assert(this->samelattice(b));
  }
  template<class R>
  void binary_operation_check(const LDVec<R>& b) const{
    assert(this->starlattice(b));
  }
  template<class R, template<class> class A,
  typename=typename std::enable_if<!std::is_base_of<LatVec,A<R>>::value>::type>
  void binary_operation_check(const A<R>&) const {}

  //! Check whether a second LQVec is approximately the same as this object
  template<class R>
  bool is(const LQVec<R>& that) const {
    return (this->samelattice(that) && this->bArray<T>::is(that));
  }
  //! Round all elements using std::round
  [[nodiscard]] LQVec<int> round() const {
    return LQVec<int>(this->lattice, this->bArray<T>::round());
  }
  //! Find the floor of all elements using std::floor
  [[nodiscard]] LQVec<int> floor() const {
    return LQVec<int>(this->lattice, this->bArray<T>::floor());
  }
  //! Find the ceiling of all elements using std::ceil
  [[nodiscard]] LQVec<int> ceil() const {
    return LQVec<int>(this->lattice, this->bArray<T>::ceil());
  }
  //
  static LQVec<T> from_invA(const Reciprocal& lat, const bArray<T>& vec, const int=1){
    assert(vec.is_row_ordered() && vec.is_contiguous() && vec.size(vec.ndim()-1)==3);
    // initialise an empy array for the relative lattice unit components
    bArray<T> rlu(vec.shape(), vec.stride());
    auto vsh = vec.shape();
    vsh.back() = 0;
    // grab the transformation matric from the lattice
    std::vector<double> fromxyz = lat.get_inverse_xyz_transform();
    for (auto i: vec.subItr(vsh))
      brille::utils::multiply_matrix_vector(rlu.ptr(i), fromxyz.data(), vec.ptr(i));
    // finally make the vector(s) with lattice and components:
    return LQVec<T>(lat, rlu);
  }

  LQVec<T> decouple() {
    return LQVec<T>(lattice, this->bArray<T>::decouple());
  }

  TO_HDF_OBJ(hkl)
  TO_HDF_STR
  FROM_HDF_OBJ(LQVec, Reciprocal, hkl)
  FROM_HDF_STR(LQVec)
protected:
  void check_array();
};

#undef LVEC_SCALAR_INPLACE_OP_DEF
#undef LVEC_ARRAY_INPLACE_OP_DEF
#undef LVEC_INIT_INT
#undef LVEC_FROM_STD
#undef LVEC_FROM_STD_VECTOR_ARRAY

#ifdef USE_HIGHFIVE
#undef TO_HDF_OBJ
#undef TO_HDF_STR
#undef FROM_HDF_OBJ
#undef FROM_HDF_STR
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// extend the lattice traits structs
template<class T> struct LatticeTraits<LDVec<T>>{
  using type = Direct;
  using star = Reciprocal;
};
template<class T> struct LatticeTraits<LQVec<T>>{
  using type = Reciprocal;
  using star = Direct;
};
#endif

/*! \brief Vector type information for Lattice and LatVec objects

Some templated functions require internal variables or return types which
depend on *which* subtype of Lattice of LatVec are provided. This traits struct
provides the typename of an appropriate LatVec subclass and its inverse for
those cases.

The two `using` typnames `type` and `star` are defined based on the templated
typename as

| templated typename | type | star |
| --- | --- | --- |
| Direct | LDVec<R> | LQVec<R> |
| Reciprocal | LQVec<R> | LDVec<R> |
| LDVec | LDVec<R> | LQVec<R> |
| LQVec | LQVec<R> | LDVec<R> |
*/
template<class, class> struct LatVecTraits{
  using type = void; //< LDVec<R> or LQVec<R>
  using star = void; //< LQVec<R> or LDVec<R>
};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template <class R, class S> struct LatVecTraits<LDVec<R>,S>{
  using type = LDVec<S>;
  using star = LQVec<S>;
};
template <class R, class S> struct LatVecTraits<LQVec<R>,S>{
  using type = LQVec<S>;
  using star = LDVec<S>;
};
template <class S> struct LatVecTraits<Direct,S>{
  using type = LDVec<S>;
  using star = LQVec<S>;
};
template <class S> struct LatVecTraits<Reciprocal,S>{
  using type = LQVec<S>;
  using star = LDVec<S>;
};
#endif
/*! \brief Easier access to LatVecTraits defined types

Is the Lattice vector type associated with the object
*/
template<class... T> using LatVecType = typename LatVecTraits<T...>::type;
/*! \brief Easier access to LatVecTraits defined types

Is the Lattice vector type associated with the dual/inverse of the object
*/
template<class... T> using LatVecStar = typename LatVecTraits<T...>::star;

#include "array_latvec.tpp"
#include "array_ldvec.tpp"
#include "array_lqvec.tpp"

}

#endif
