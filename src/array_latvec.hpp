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

#ifndef BRILLE_LATVEC_CLASS_H_
#define BRILLE_LATVEC_CLASS_H_

#include <typeinfo> // for std::bad_cast
#include <exception>
#include "lattice.hpp"
#include "array2.hpp"

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


template<class... T> struct ArrayTraits{
  static constexpr bool array = false;
  static constexpr bool latvec = false;
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


template<class T, template<class> class A>
inline constexpr bool isArray = ArrayTraits<A<T>>::array;
template<class T, template<class> class A>
inline constexpr bool isLatVec = ArrayTraits<A<T>>::latvec;
template<class T, template<class> class A>
inline constexpr bool isBareArray = isArray<T,A> && !isLatVec<T,A>;

template<class T, template<class> class A, class R, template<class> class B>
inline constexpr bool bothArrays = isArray<T,A> && isArray<R,B>;
template<class T, template<class> class A, class R, template<class> class B>
inline constexpr bool bareArrays = isBareArray<T,A> && isBareArray<R,B>;
template<class T, template<class> class A, class R, template<class> class B>
inline constexpr bool bothLatVecs = isLatVec<T,A> && isLatVec<R,B>;

template<class T>
bool equal_shapes(const std::vector<T>& a, const std::vector<T>&b){
  bool ok = a.size() == b.size();
  if (ok) ok = std::equal(a.begin(), a.end(), b.begin());
  return ok;
}
template<class T, size_t Nel>
bool equal_shapes(const std::array<T,Nel>& a, const std::array<T,Nel>&b){
  return std::equal(a.begin(), a.end(), b.begin());
}

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
  if (!equal_shapes(itr.shape(), this->shape())) throw std::runtime_error("Incompatible shaped Array for binary operation");\
  for (auto [ox, ax, bx]: itr) this->_data[this->s2l_d(ax)] X b[bx];\
  return *this;\
}

#define LVEC_INIT_INT(N,L,X) N(const L &lat, const X &n)\
: bArray<T>(static_cast<brille::ind_t>(n),3), lattice(lat)\
{}

/*! \brief 3-vector(s) expressed in units of a Direct lattice

By adding a Direct lattice to a 3-element bArray this class represents one
or more 3-vector in units of a real-space-spanning lattice.
*/
template<class T>
class LDVec: public LatVec, public bArray<T>{
  Direct lattice;
public:
  // default constructor for zero three-vectors:
  explicit LDVec(const Direct& lat=Direct())
  : bArray<T>(0,3), lattice(lat)
  {
  }
  //! integer number of three-vector constructor (macroed as templates can't be distinguished)
  LVEC_INIT_INT(LDVec,Direct,int)
  LVEC_INIT_INT(LDVec,Direct,long)
  LVEC_INIT_INT(LDVec,Direct,long long)
  LVEC_INIT_INT(LDVec,Direct,unsigned)
  LVEC_INIT_INT(LDVec,Direct,unsigned long)
  LVEC_INIT_INT(LDVec,Direct,unsigned long long)
  //! Fowarding constructor to let bArray deal with everything else
  template<typename... Args> LDVec(const Direct& lat, Args... args)
  : bArray<T>(args...), lattice(lat)
  {
    this->check_array();
  }
  template<class R>
  LDVec(const LDVec<R>& other)
  : bArray<T>(other.get_hkl()), lattice(other.get_lattice())
  {
  }


  Direct get_lattice() const { return lattice; }
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
  bArray<T> get_hkl() const;
  //! Extract the coordinates in *an* orthonormal frame
  bArray<double> get_xyz() const;
  //! Return the vector(s) expressed in units of the Reciprocal lattice
  LQVec<double> star() const;

  //! Determine the scalar product between two vectors in the object
  double dot(const size_t i, const size_t j) const;
  //! Determine the absolute length of a vector in the object
  double norm(const size_t i) const { return sqrt(this->dot(i,i)); }
  //! Determine the cross product of two vectors in the object
  LDVec<double> cross(const size_t i, const size_t j) const;

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
  LDVec<int> round() const {
    return LDVec<int>(this->lattice, this->bArray<T>::round());
  }
  //! Find the floor of all elements using std::floor
  LDVec<int> floor() const {
    return LDVec<int>(this->lattice, this->bArray<T>::floor());
  }
  //! Find the ceiling of all elements using std::ceil
  LDVec<int> ceil() const {
    return LDVec<int>(this->lattice, this->bArray<T>::ceil());
  }

  LDVec<T> decouple() {
    return LDVec<T>(lattice, this->bArray<T>::decouple());
  }
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
  explicit LQVec(const Reciprocal& lat=Reciprocal())
  : bArray<T>(0,3), lattice(lat)
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
  LQVec(const Reciprocal& lat, Args... args)
  : bArray<T>(args...), lattice(lat)
  {
    this->check_array();
  }

  template<class R>
  LQVec(const LQVec<R>& other)
  : bArray<T>(other.get_hkl()), lattice(other.get_lattice())
  {
  }

  //! Explicit copy constructor: required in gcc 9+ since we define our own operator= below:
  // LQVec(const LQVec<T>& other, const bool m, const bool d): bArray<T>(other.get_hkl(), m, d), lattice{other.get_lattice()} {}
  // LQVec(const LQVec<T>& other, const bool m=false, const bool d=false): bArray<T>(other.get_hkl(),m,d), lattice(other.get_lattice()) {}

  // //! Type conversion assignment:
  // template<class R> LQVec<T>& operator=(const LQVec<R>& other){
  //   this->lattice = other.get_lattice();
  //   *this = bArray<T>(other.get_hkl());
  //   return *this;
  // }
  //! Assignment operator reusing data if we can
  // LQVec<T>& operator=(const LQVec<T>& other){
  //   if (this != &other){ // do nothing if called by, e.g., a = a;
  //     this->lattice = other.get_lattice();
  //     this->bArray<T>::operator=(other);
  //   }
  //   return *this;
  // }
  // LQVec<T>& operator=(LQVec<T>&& other) noexcept {
  //
  // }

  Reciprocal get_lattice() const { return lattice; }
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
  bArray<T> get_hkl() const;
  /*! Extract the coordinates in an orthonormal frame with its first axis, x,
  along a*, its second, y, perpendicular with y⋅b*>0 , and it's third forming
  the right-handed set z=x×y.
  */
  bArray<double> get_xyz() const;
  //! Return the vector(s) expressed in units of the Direct lattice
  LDVec<double> star() const;

  //! Determine the scalar product between two vectors in the object.
  double dot(const size_t i, const size_t j) const;
  //! Determine the absolute length of a vector in the object.
  double norm(const size_t i) const { return sqrt(this->dot(i,i)); }
  //! Determine the cross product of two vectors in the object.
  LQVec<double> cross(const size_t i, const size_t j) const;

  LVEC_SCALAR_INPLACE_OP_DEF(LQVec,+=)
  LVEC_SCALAR_INPLACE_OP_DEF(LQVec,-=)
  LVEC_SCALAR_INPLACE_OP_DEF(LQVec,*=)
  LVEC_SCALAR_INPLACE_OP_DEF(LQVec,/=)

  LVEC_ARRAY_INPLACE_OP_DEF(LQVec,+=)
  LVEC_ARRAY_INPLACE_OP_DEF(LQVec,-=)
  LVEC_ARRAY_INPLACE_OP_DEF(LQVec,*=)
  LVEC_ARRAY_INPLACE_OP_DEF(LQVec,/=)
  // LQVec<T> operator -();

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
  LQVec<int> round() const {
    return LQVec<int>(this->lattice, this->bArray<T>::round());
  }
  //! Find the floor of all elements using std::floor
  LQVec<int> floor() const {
    return LQVec<int>(this->lattice, this->bArray<T>::floor());
  }
  //! Find the ceiling of all elements using std::ceil
  LQVec<int> ceil() const {
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
protected:
  void check_array();
};

#undef LVEC_SCALAR_INPLACE_OP_DEF
#undef LVEC_ARRAY_INPLACE_OP_DEF
#undef LVEC_INIT_INT

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

template<class... T> using LatVecType = typename LatVecTraits<T...>::type;
template<class... T> using LatVecStar = typename LatVecTraits<T...>::star;

#include "array_latvec.tpp"
#include "array_ldvec.tpp"
#include "array_lqvec.tpp"

}

#endif
