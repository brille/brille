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

#ifndef _LATVEC_CLASS_H_
#define _LATVEC_CLASS_H_

#include <typeinfo> // for std::bad_cast
#include <exception>
#include "lattice.hpp"
#include "array.hpp"

/*! \brief Superclass to identify both LDVec and LQVec

The two Lattice vector classes, LDVec and LQVec, are subclasses of brille::Array.
In order to distinguish between the Lattice vector types and "bare" brille::Array
objects in operator- and function-overloading it is advantageous to have a
shared superclass on which to enable templates. Thus LatVec is a superclass used
for logic only which has no properties and defines no methods.
*/
class LatVec{};
template<class T, class P> class LDVec;
template<class T, class P> class LQVec;


template<class... T> struct ArrayTraits{
  static constexpr bool array = false;
  static constexpr bool latvec = false;
};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T, class P> struct ArrayTraits<brille::Array<T,P>>{
  static constexpr bool array = true;
  static constexpr bool latvec = false;
};
template<class T, class P> struct ArrayTraits<LDVec<T,P>>{
  static constexpr bool array = true;
  static constexpr bool latvec = true;
};
template<class T, class P> struct ArrayTraits<LQVec<T,P>>{
  static constexpr bool array = true;
  static constexpr bool latvec = true;
};
#endif


template<class T, class P, template<class, class> class A>
inline constexpr bool isArray = ArrayTraits<A<T,P>>::array;
template<class T, class P, template<class, class> class A>
inline constexpr bool isLatVec = ArrayTraits<A<T,P>>::latvec;
template<class T, class P, template<class, class> class A>
inline constexpr bool isBareArray = isArray<T,P,A> && !isLatVec<T,P,A>;

template<class T, class P, template<class,class> class A, class R, class Q, template<class,class> class B>
inline constexpr bool bothArrays = isArray<T,P,A> && isArray<R,Q,B>;
template<class T, class P, template<class,class> class A, class R, class Q, template<class,class> class B>
inline constexpr bool bareArrays = isBareArray<T,P,A> && isBareArray<R,Q,B>;
template<class T, class P, template<class,class> class A, class R, class Q, template<class,class> class B>
inline constexpr bool bothLatVecs = isLatVec<T,P,A> && isLatVec<R,Q,B>;

template<class T>
bool equal_shapes(const std::vector<T>& a, const std::vector<T>&b){
  bool ok = a.size() == b.size();
  if (ok) ok = std::equal(a.begin(), a.end(), b.begin());
  return ok;
}

#define LVEC_SCALAR_INPLACE_OP_DEF(L,X) \
L <T,P>&\
operator X (const T& v){\
  for (auto x: SubIt(this->shape())) this->_data[this->s2l_d(x)] X v;\
  return *this;\
}
#define LVEC_ARRAY_INPLACE_OP_DEF(L,X) \
template<class Q, template<class,class> class A>\
std::enable_if_t<isArray<T,Q,A>, L <T,P>>&\
operator X (const A<T,Q>& b){\
  auto itr = BroadcastIt(this->shape(), b.shape());\
  if (!equal_shapes(itr.shape(), this->shape())) throw std::runtime_error("Incompatible shaped Array for binary operation");\
  for (auto [ox, ax, bx]: itr) this->_data[this->s2l_d(ax)] X b[bx];\
  return *this;\
}

#define LVEC_INIT_INT(N,L,X) N(const L &lat, const X &n)\
: brille::Array<T,P>(static_cast<brille::ind_t>(n),3), lattice(lat)\
{}

/*! \brief 3-vector(s) expressed in units of a Direct lattice

By adding a Direct lattice to a 3-element brille::Array this class represents one
or more 3-vector in units of a real-space-spanning lattice.
*/
template<class T, class P=char>
class LDVec: public LatVec, public brille::Array<T,P>{
  Direct lattice;
public:
  // default constructor for zero three-vectors:
  explicit LDVec(const Direct& lat=Direct())
  : brille::Array<T,P>(0,3), lattice(lat)
  {
  }
  //! integer number of three-vector constructor (macroed as templates can't be distinguished)
  LVEC_INIT_INT(LDVec,Direct,int)
  LVEC_INIT_INT(LDVec,Direct,long)
  LVEC_INIT_INT(LDVec,Direct,long long)
  LVEC_INIT_INT(LDVec,Direct,unsigned)
  LVEC_INIT_INT(LDVec,Direct,unsigned long)
  LVEC_INIT_INT(LDVec,Direct,unsigned long long)
  //! Fowarding constructor to let brille::Array deal with everything else
  template<typename... Args> LDVec(const Direct& lat, Args... args)
  : brille::Array<T,P>(args...), lattice(lat)
  {
    this->check_array();
  }
  template<class R, class Q>
  LDVec(const LDVec<R,Q>& other)
  : brille::Array<T,P>(other.get_hkl()), lattice(other.get_lattice())
  {
  }


  Direct get_lattice() const { return lattice; }
  template<class R, class Q> bool samelattice(const LDVec<R,Q> *vec) const { return lattice.issame(vec->get_lattice()); }
  template<class R, class Q> bool samelattice(const LQVec<R,Q> *)    const { return false; }
  template<class R, class Q> bool starlattice(const LDVec<R,Q> *)    const { return false; }
  template<class R, class Q> bool starlattice(const LQVec<R,Q> *vec) const { return lattice.isstar(vec->get_lattice()); }
  template<class R, class Q> bool samelattice(const LDVec<R,Q> &vec) const { return lattice.issame(vec.get_lattice()); }
  template<class R, class Q> bool samelattice(const LQVec<R,Q> &)    const { return false; }
  template<class R, class Q> bool starlattice(const LDVec<R,Q> &)    const { return false; }
  template<class R, class Q> bool starlattice(const LQVec<R,Q> &vec) const { return lattice.isstar(vec.get_lattice()); }
  // extract overloads preserving lattice information
  //! Return a non-copying view into the  LDVec
  template<typename... A> LDVec<T,P> view(A... args) const;
  //! Return a copied subset of the LDVec
  template<typename... A> LDVec<T,brille::ref_ptr_t> extract(A... args) const;
  //! Extract just the coordinates in units of the Direct lattice (strip off the lattice information)
  brille::Array<T,P> get_hkl() const;
  //! Extract the coordinates in *an* orthonormal frame
  brille::Array<double,brille::ref_ptr_t> get_xyz() const;
  //! Return the vector(s) expressed in units of the Reciprocal lattice
  LQVec<double,brille::ref_ptr_t> star() const;

  //! Determine the scalar product between two vectors in the object
  double dot(const size_t i, const size_t j) const;
  //! Determine the absolute length of a vector in the object
  double norm(const size_t i) const { return sqrt(this->dot(i,i)); }
  //! Determine the cross product of two vectors in the object
  LDVec<double,brille::ref_ptr_t> cross(const size_t i, const size_t j) const;

  // LDVec<T,P>& operator-=(const LDVec<T,P>& av);
  // LDVec<T,P>& operator+=(const LDVec<T,P>& av);

  LVEC_SCALAR_INPLACE_OP_DEF(LDVec,+=)
  LVEC_SCALAR_INPLACE_OP_DEF(LDVec,-=)
  LVEC_SCALAR_INPLACE_OP_DEF(LDVec,*=)
  LVEC_SCALAR_INPLACE_OP_DEF(LDVec,/=)

  LVEC_ARRAY_INPLACE_OP_DEF(LDVec,+=)
  LVEC_ARRAY_INPLACE_OP_DEF(LDVec,-=)
  LVEC_ARRAY_INPLACE_OP_DEF(LDVec,*=)
  LVEC_ARRAY_INPLACE_OP_DEF(LDVec,/=)


  template<class R, class Q>
  void binary_operation_check(const LDVec<R,Q>& b) const{
    assert(this->samelattice(b));
  }
  template<class R, class Q>
  void binary_operation_check(const LQVec<R,Q>& b) const{
    assert(this->starlattice(b));
  }
  template<class R, class Q, template<class,class> class A,
  typename=typename std::enable_if<!std::is_base_of<LatVec,A<R,Q>>::value>::type>
  void binary_operation_check(const A<R,Q>&) const {}

  //! Check whether a second LDVec is approximately the same as this object
  template<class R, class Q>
  bool is(const LDVec<R,Q>& that){
    return (this->samelattice(that) && this->brille::Array<T,P>::is(that));
  }
  //! Round all elements using std::round
  LDVec<int,brille::ref_ptr_t> round() const {
    return LDVec<int,brille::ref_ptr_t>(this->lattice, this->brille::Array<T,P>::round());
  }
  //! Find the floor of all elements using std::floor
  LDVec<int,brille::ref_ptr_t> floor() const {
    return LDVec<int,brille::ref_ptr_t>(this->lattice, this->brille::Array<T,P>::floor());
  }
  //! Find the ceiling of all elements using std::ceil
  LDVec<int,brille::ref_ptr_t> ceil() const {
    return LDVec<int,brille::ref_ptr_t>(this->lattice, this->brille::Array<T,P>::ceil());
  }

  LDVec<T,brille::ref_ptr_t> decouple() {
    return LDVec<T,brille::ref_ptr_t>(lattice, this->brille::Array<T,P>::decouple());
  }
protected:
  void check_array();
};


/*! \brief 3-vector(s) expressed in units of a Reciprocal lattice

By adding a Reciprocal lattice to a 3-element brille::Array this class represents one
or more 3-vector in units of a reciprocal-space-spanning lattice.
*/
template<class T, class P=char>
class LQVec:  public LatVec, public brille::Array<T,P>{
  Reciprocal lattice;
public:
  // default constructor for zero three-vectors:
  explicit LQVec(const Reciprocal& lat=Reciprocal())
  : brille::Array<T,P>(0,3), lattice(lat)
  {
  }
  //! integer number of three-vector constructor (macroed as templates can't be distinguished)
  LVEC_INIT_INT(LQVec,Reciprocal,int)
  LVEC_INIT_INT(LQVec,Reciprocal,long)
  LVEC_INIT_INT(LQVec,Reciprocal,long long)
  LVEC_INIT_INT(LQVec,Reciprocal,unsigned)
  LVEC_INIT_INT(LQVec,Reciprocal,unsigned long)
  LVEC_INIT_INT(LQVec,Reciprocal,unsigned long long)
  //! Fowarding constructor to let brille::Array deal with everything else
  template<typename... Args>
  LQVec(const Reciprocal& lat, Args... args)
  : brille::Array<T,P>(args...), lattice(lat)
  {
    this->check_array();
  }

  template<class R, class Q>
  LQVec(const LQVec<R,Q>& other)
  : brille::Array<T,P>(other.get_hkl()), lattice(other.get_lattice())
  {
  }

  //! Explicit copy constructor: required in gcc 9+ since we define our own operator= below:
  // LQVec(const LQVec<T,P>& other, const bool m, const bool d): brille::Array<T,P>(other.get_hkl(), m, d), lattice{other.get_lattice()} {}
  // LQVec(const LQVec<T,P>& other, const bool m=false, const bool d=false): brille::Array<T,P>(other.get_hkl(),m,d), lattice(other.get_lattice()) {}

  // //! Type conversion assignment:
  // template<class R> LQVec<T,P>& operator=(const LQVec<R,P>& other){
  //   this->lattice = other.get_lattice();
  //   *this = brille::Array<T,P>(other.get_hkl());
  //   return *this;
  // }
  //! Assignment operator reusing data if we can
  // LQVec<T,P>& operator=(const LQVec<T,P>& other){
  //   if (this != &other){ // do nothing if called by, e.g., a = a;
  //     this->lattice = other.get_lattice();
  //     this->brille::Array<T,P>::operator=(other);
  //   }
  //   return *this;
  // }
  // LQVec<T,P>& operator=(LQVec<T,P>&& other) noexcept {
  //
  // }

  Reciprocal get_lattice() const { return lattice; }
  template<class R, class Q> bool samelattice(const LQVec<R,Q> *vec) const { return lattice.issame(vec->get_lattice()); }
  template<class R, class Q> bool samelattice(const LDVec<R,Q> *)    const { return false; }
  template<class R, class Q> bool starlattice(const LQVec<R,Q> *)    const { return false; }
  template<class R, class Q> bool starlattice(const LDVec<R,Q> *vec) const { return lattice.isstar(vec->get_lattice()); }
  template<class R, class Q> bool samelattice(const LQVec<R,Q> &vec) const { return lattice.issame(vec.get_lattice()); }
  template<class R, class Q> bool samelattice(const LDVec<R,Q> &)    const { return false; }
  template<class R, class Q> bool starlattice(const LQVec<R,Q> &)    const { return false; }
  template<class R, class Q> bool starlattice(const LDVec<R,Q> &vec) const { return lattice.isstar(vec.get_lattice()); }
  // extract overloads preserving lattice information
  //! Return a non-copying view into the  LQVec
  template<typename... A> LQVec<T,P> view(A... args) const;
  //! Return a copied subset of the LQVec
  template<typename... A> LQVec<T,brille::ref_ptr_t> extract(A... args) const;
  //! Extract just the coordinates in units of the Reciprocal lattice (strip off the lattice information)
  brille::Array<T,P> get_hkl() const;
  /*! Extract the coordinates in an orthonormal frame with its first axis, x,
  along a*, its second, y, perpendicular with y⋅b*>0 , and it's third forming
  the right-handed set z=x×y.
  */
  brille::Array<double,brille::ref_ptr_t> get_xyz() const;
  //! Return the vector(s) expressed in units of the Direct lattice
  LDVec<double,brille::ref_ptr_t> star() const;

  //! Determine the scalar product between two vectors in the object.
  double dot(const size_t i, const size_t j) const;
  //! Determine the absolute length of a vector in the object.
  double norm(const size_t i) const { return sqrt(this->dot(i,i)); }
  //! Determine the cross product of two vectors in the object.
  LQVec<double,brille::ref_ptr_t> cross(const size_t i, const size_t j) const;

  LVEC_SCALAR_INPLACE_OP_DEF(LQVec,+=)
  LVEC_SCALAR_INPLACE_OP_DEF(LQVec,-=)
  LVEC_SCALAR_INPLACE_OP_DEF(LQVec,*=)
  LVEC_SCALAR_INPLACE_OP_DEF(LQVec,/=)

  LVEC_ARRAY_INPLACE_OP_DEF(LQVec,+=)
  LVEC_ARRAY_INPLACE_OP_DEF(LQVec,-=)
  LVEC_ARRAY_INPLACE_OP_DEF(LQVec,*=)
  LVEC_ARRAY_INPLACE_OP_DEF(LQVec,/=)
  // LQVec<T,P> operator -();

  template<class R, class Q>
  void binary_operation_check(const LQVec<R,Q>& b) const{
    assert(this->samelattice(b));
  }
  template<class R, class Q>
  void binary_operation_check(const LDVec<R,Q>& b) const{
    assert(this->starlattice(b));
  }
  template<class R, class Q, template<class,class> class A,
  typename=typename std::enable_if<!std::is_base_of<LatVec,A<R,Q>>::value>::type>
  void binary_operation_check(const A<R,Q>&) const {}

  //! Check whether a second LQVec is approximately the same as this object
  template<class R, class Q>
  bool is(const LQVec<R,Q>& that) const {
    return (this->samelattice(that) && this->brille::Array<T,P>::is(that));
  }
  //! Round all elements using std::round
  LQVec<int,brille::ref_ptr_t> round() const {
    return LQVec<int,brille::ref_ptr_t>(this->lattice, this->brille::Array<T,P>::round());
  }
  //! Find the floor of all elements using std::floor
  LQVec<int,brille::ref_ptr_t> floor() const {
    return LQVec<int,brille::ref_ptr_t>(this->lattice, this->brille::Array<T,P>::floor());
  }
  //! Find the ceiling of all elements using std::ceil
  LQVec<int,brille::ref_ptr_t> ceil() const {
    return LQVec<int,brille::ref_ptr_t>(this->lattice, this->brille::Array<T,P>::ceil());
  }
  //
  static LQVec<T,brille::ref_ptr_t> from_invA(const Reciprocal& lat, const brille::Array<T,P>& vec, const int=1){
    assert(vec.is_row_ordered() && vec.is_contiguous() && vec.size(vec.ndim()-1)==3);
    // initialise an empy array for the relative lattice unit components
    brille::Array<T,brille::ref_ptr_t> rlu(vec.shape(), vec.stride());
    auto vsh = vec.shape();
    vsh.back() = 0;
    // grab the transformation matric from the lattice
    std::vector<double> fromxyz = lat.get_inverse_xyz_transform();
    for (auto i: SubIt(vec.shape(), vsh))
      brille::utils::multiply_matrix_vector(rlu.ptr(i), fromxyz.data(), vec.ptr(i));
    // finally make the vector(s) with lattice and components:
    return LQVec<T,brille::ref_ptr_t>(lat, rlu);
  }

  LQVec<T,brille::ref_ptr_t> decouple() {
    return LQVec<T,brille::ref_ptr_t>(lattice, this->brille::Array<T,P>::decouple());
  }
protected:
  void check_array();
};

#undef LVEC_SCALAR_INPLACE_OP_DEF
#undef LVEC_ARRAY_INPLACE_OP_DEF
#undef LVEC_INIT_INT

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// extend the lattice traits structs
template<class T, class P> struct LatticeTraits<LDVec<T,P>>{
  using type = Direct;
  using star = Reciprocal;
};
template<class T, class P> struct LatticeTraits<LQVec<T,P>>{
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
| Direct | LDVec<R,P> | LQVec<R,P> |
| Reciprocal | LQVec<R,P> | LDVec<R,P> |
| LDVec | LDVec<R,P> | LQVec<R,P> |
| LQVec | LQVec<R,P> | LDVec<R,P> |
*/
template<class T, class R> struct LatVecTraits{
  using type = void; //< LDVec<R,P> or LQVec<R,P>
  using star = void; //< LQVec<R,P> or LDVec<R,P>
};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
// the star types *always* require creating a new Array and so will always have brille::ref_ptr_t
template <class R, class Q, class S> struct LatVecTraits<LDVec<R,Q>,S>{
  using type = LDVec<S,Q>;
  using star = LQVec<S,brille::ref_ptr_t>;
};
template <class R, class Q, class S> struct LatVecTraits<LQVec<R,Q>,S>{
  using type = LQVec<S,Q>;
  using star = LDVec<S,brille::ref_ptr_t>;
};
template <class S> struct LatVecTraits<Direct,S>{
  using type = LDVec<S,brille::ref_ptr_t>;
  using star = LQVec<S,brille::ref_ptr_t>;
};
template <class S> struct LatVecTraits<Reciprocal,S>{
  using type = LQVec<S,brille::ref_ptr_t>;
  using star = LDVec<S,brille::ref_ptr_t>;
};
#endif

template<class... T> using LatVecType = typename LatVecTraits<T...>::type;
template<class... T> using LatVecStar = typename LatVecTraits<T...>::star;

#include "array_latvec.tpp"
#include "array_ldvec.tpp"
#include "array_lqvec.tpp"

#endif
