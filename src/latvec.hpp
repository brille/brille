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
#include "arrayvector.hpp"

/*! \brief Superclass to identify both LDVec and LQVec

The two Lattice vector classes, LDVec and LQVec, are subclasses of ArrayVector.
In order to distinguish between the Lattice vector types and "bare" ArrayVector
objects in operator- and function-overloading it is advantageous to have a
shared superclass on which to enable templates. Thus LatVec is a superclass used
for logic only which has no properties and defines no methods.
*/
class LatVec{};
template<typename T> class LDVec;
template<typename T> class LQVec;

/*! \brief 3-vector(s) expressed in units of a Direct lattice

By adding a Direct lattice to a 3-element ArrayVector this class represents one
or more 3-vector in units of a real-space-spanning lattice.
*/
template<typename T> class LDVec: public LatVec, public ArrayVector<T>{
  Direct lattice;
public:
  //! Default constructor, enforces that 3-vectors are held
  LDVec(const Direct& lat=Direct(), const size_t n=0, const T *d=nullptr): ArrayVector<T>(3,n,d), lattice(lat){}
  //! Constructor taking a possibly non-contigous and/or non-row-ordered input
  template<class Integer, typename=typename std::enable_if<std::is_integral<Integer>::value>::type>
  LDVec(const Direct& lat,
        const T*d,
        const std::vector<Integer>& shape,
        const std::vector<Integer>& strides,
        const int flag=1): ArrayVector<T>(d,shape,strides), lattice(lat){
    this->check_arrayvector(flag);
  }
  //! Copy constructor, optionally verifying that only 3-element arrays are provided.
  LDVec(const Direct& lat, const ArrayVector<T>& vec, const int flag=1): ArrayVector<T>(vec), lattice(lat){ this->check_arrayvector(flag); }
  //! [Optional type conversion] copy constructor
  template<class R> LDVec(const LDVec<R>& vec): ArrayVector<T>(vec.numel(),vec.size(),vec.data()), lattice(vec.get_lattice()) {}
  //! std::vector<std::array<T,3>> copy constructor
  template<class R> LDVec(const Direct& lat, const std::vector<std::array<R,3>>& va): ArrayVector<T>(va), lattice(lat){}
  //! Explicit copy constructor
  // required in gcc 9+ since we define our own operator= below:
  LDVec(const LDVec<T>& other): ArrayVector<T>(3u,other.size(),other.data()), lattice(other.get_lattice()) {}
  //! Assignment operator reusing data if we can
  LDVec<T>& operator=(const LDVec<T>& other){
    if (this != &other){ // do nothing if called by, e.g., a = a;
      this->lattice = other.get_lattice();
      size_t m = other.numel();
      size_t n = other.size();
      // reuse the data-block if we can. otherwise, refresh/resize it
      if ( m !=this->numel() ) this->refresh(m,n);
      if ( n !=this->size()  ) this->resize(n);
      // copy-over the data (if it exists)
      if (other._data && m && n)
        for(size_t i=0; i<m*n; i++)
          this->_data[i] = other._data[i];
    }
    return *this;
  }
  //! Extract the 3-vector with index `i`
  const LDVec<T> operator[](const size_t i) const{
    bool isok = i < this->size();
    LDVec<T> out(this->lattice, isok ? 1u: 0u);
    if (isok) out = this->get(i);
    return out;
  }

  Direct get_lattice() const { return lattice; }
  template<typename R> bool samelattice(const LDVec<R> *vec) const { return lattice.issame(vec->get_lattice()); }
  template<typename R> bool samelattice(const LQVec<R> *)    const { return false; }
  template<typename R> bool starlattice(const LDVec<R> *)    const { return false; }
  template<typename R> bool starlattice(const LQVec<R> *vec) const { return lattice.isstar(vec->get_lattice()); }
  template<typename R> bool samelattice(const LDVec<R> &vec) const { return lattice.issame(vec.get_lattice()); }
  template<typename R> bool samelattice(const LQVec<R> &)    const { return false; }
  template<typename R> bool starlattice(const LDVec<R> &)    const { return false; }
  template<typename R> bool starlattice(const LQVec<R> &vec) const { return lattice.isstar(vec.get_lattice()); }
  // extract overloads preserving lattice information
  //! Return the ith single-array LDVec
  LDVec<T> extract(const size_t i=0) const ;
  /*! \brief Extract the first `num` vectors from an LDVec
    @param num the number of vectors to return
    @returns a LDVec containing the `num` vectors */
  LDVec<T> first(const size_t num) const;
  /*! Return a collection of arrays from the LDVec
    @param n the number of arrays to return
    @param i a pointer to the first index of the n arrays to return
    @returns A LDVec containing the n arrays
  */
  LDVec<T> extract(const size_t n, const size_t *i) const;
  /*! Return a collection of arrays from the LDVec
    @param idx a reference to an LDVec containing the to-be-returned indices
    @returns A LDVec containing the indicies indicated by idx
  */
  LDVec<T> extract(const ArrayVector<size_t>& idx) const;
  /*! Return a collection of arrays from the LDVec
    @param tfvec a reference to an ArrayVector<bool> with true for the to-be-returned indices
    @returns An LDVec containing the indicies indicated by tfvec
  */
  LDVec<T> extract(const ArrayVector<bool>& idx) const;
  LDVec<T> extract(const std::vector<bool>& idx) const;
  //! extract the 3-vector with index `i`
  LDVec<T> get(const size_t i) const;
  //! Extract just the coordinates in units of the Direct lattice (strip off the lattice information)
  ArrayVector<T> get_hkl() const;
  //! Extract the coordinates in *an* orthonormal frame
  ArrayVector<double> get_xyz() const;
  //! Return the vector(s) expressed in units of the Reciprocal lattice
  LQVec<double> star() const;

  //! Determine the scalar product between two vectors in the object
  double dot(const size_t i, const size_t j) const;
  //! Determine the absolute length of a vector in the object
  double norm(const size_t i) const { return sqrt(this->dot(i,i)); }
  //! Determine the cross product of two vectors in the object
  LDVec<double> cross(const size_t i, const size_t j) const;

  LDVec<T>& operator-=(const LDVec<T>& av);
  LDVec<T>& operator+=(const LDVec<T>& av);

  LDVec<T>& operator+=(const T& av);
  LDVec<T>& operator-=(const T& av);
  LDVec<T>& operator*=(const T& av);
  LDVec<T>& operator/=(const T& av);
  // LDVec<T> operator -(); // unary subtraction is not actually defined anywhere

  //! Verify two LDVec objects are compatible for binary operations
  template<class R> AVSizeInfo consistency_check(const LDVec<R>& b) const {
    if (!(this->samelattice(b))) throw std::runtime_error("arithmetic between Lattice vectors requires they have the same lattice");
    return this->ArrayVector<T>::consistency_check(b);
  }
  template<class R> AVSizeInfo consistency_check(const LQVec<R>& b) const {
    if (!(this->starlattice(b))) throw std::runtime_error("Real and reciprocal space vectors must be in dual lattices.");
    return this->ArrayVector<T>::consistency_check(b);
  }
  //! Verify that a ArrayVector object is consistent for binary operations
  template<class R, template<class> class A,
    typename=typename std::enable_if<!std::is_base_of<LatVec,A<R>>::value && std::is_base_of<ArrayVector<R>,A<R>>::value >::type
    >
  AVSizeInfo consistency_check(const A<R>& b) const {
    return this->ArrayVector<T>::consistency_check(b); // b has no lattice, so nothing to check
  }
  //! Check whether a second LDVec is approximately the same as this object
  template<typename R> bool isapprox(const LDVec<R>& that){ return (this->samelattice(that) && this->ArrayVector<T>::isapprox(that)); }
  //! Check whether two vectors in the object are approximately the same
  bool isapprox(const size_t i, const size_t j) const { return this->ArrayVector<T>::isapprox(i,j);}
  //! Round all elements using std::round
  LDVec<int> round() const {return LDVec<int>(this->lattice, this->ArrayVector<T>::round()); }
  //! Find the floor of all elements using std::floor
  LDVec<int> floor() const {return LDVec<int>(this->lattice, this->ArrayVector<T>::floor()); }
  //! Find the ceiling of all elements using std::ceil
  LDVec<int> ceil() const {return LDVec<int>(this->lattice, this->ArrayVector<T>::ceil()); }
protected:
  void check_arrayvector(const int);
};


/*! \brief 3-vector(s) expressed in units of a Reciprocal lattice

By adding a Reciprocal lattice to a 3-element ArrayVector this class represents one
or more 3-vector in units of a reciprocal-space-spanning lattice.
*/
template<typename T> class LQVec:  public LatVec, public ArrayVector<T>{
  Reciprocal lattice;
public:
  //! Default constructor, enforces that 3-vectors are held
  LQVec(const Reciprocal& lat=Reciprocal(), const size_t n=0, const T *d=nullptr): ArrayVector<T>(3,n,d), lattice(lat){}
  //! Constructor taking a possibly non-contigous and/or non-row-ordered input
  template<class Integer, typename=typename std::enable_if<std::is_integral<Integer>::value>::type>
  LQVec(const Reciprocal& lat, const T*d,
        const std::vector<Integer>& shape,
        const std::vector<Integer>& strides,
        const int flag=1): ArrayVector<T>(d,shape,strides), lattice(lat){
    this->check_arrayvector(flag);
  }
  //! Copy constructor, optionally verifying that only 3-element arrays are provided.
  LQVec(const Reciprocal& lat, const ArrayVector<T>& vec, const int flag=1): ArrayVector<T>(vec), lattice(lat){  this->check_arrayvector(flag); }
  //! [Optional type conversion] copy constructor
  template<class R> LQVec(const LQVec<R>& vec): ArrayVector<T>(vec.numel(),vec.size(),vec.data()), lattice(vec.get_lattice()) {}
  //! std::vector<std::array<T,3>> copy constructor
  template<class R> LQVec(const Reciprocal& lat, const std::vector<std::array<R,3>>& va): ArrayVector<T>(va), lattice(lat){}
  //! Explicit copy constructor
  // required in gcc 9+ since we define our own operator= below:
  LQVec(const LQVec<T>& other): ArrayVector<T>(3u,other.size(),other.data()), lattice(other.get_lattice()) {}
  //! Assignment operator reusing data if we can
  LQVec<T>& operator=(const LQVec<T>& other){
    if (this != &other){ // do nothing if called by, e.g., a = a;
      this->lattice = other.get_lattice();
      size_t m = other.numel();
      size_t n = other.size();
      // reuse the data-block if we can. otherwise, refresh/resize it
      if ( m !=this->numel() ) this->refresh(m,n);
      if ( n !=this->size()  ) this->resize(n);
      // copy-over the data (if it exists)
      if (m && n)
        for(size_t i=0; i<m*n; i++)
          this->_data[i] = other._data[i];
    }
    return *this;
  }
  //! Extract the 3-vector with index `i`
  const LQVec<T> operator[](const size_t i) const{
    bool isok = i < this->size();
    LQVec<T> out(this->lattice, isok ? 1u: 0u);
    if (isok) out = this->get(i);
    return out;
  }

  Reciprocal get_lattice() const { return lattice; }
  template<typename R> bool samelattice(const LQVec<R> *vec) const { return lattice.issame(vec->get_lattice()); }
  template<typename R> bool samelattice(const LDVec<R> *)    const { return false; }
  template<typename R> bool starlattice(const LQVec<R> *)    const { return false; }
  template<typename R> bool starlattice(const LDVec<R> *vec) const { return lattice.isstar(vec->get_lattice()); }
  template<typename R> bool samelattice(const LQVec<R> &vec) const { return lattice.issame(vec.get_lattice()); }
  template<typename R> bool samelattice(const LDVec<R> &)    const { return false; }
  template<typename R> bool starlattice(const LQVec<R> &)    const { return false; }
  template<typename R> bool starlattice(const LDVec<R> &vec) const { return lattice.isstar(vec.get_lattice()); }
  // extract overloads preserving lattice information
  //! Return the ith single-array LQVec
  LQVec<T> extract(const size_t i=0) const ;
  /*! \brief Extract the first `num` vectors from an LQVec
    @param num the number of vectors to return
    @returns a LQec containing the `num` vectors */
  LQVec<T> first(const size_t num) const;
  /*! Return a collection of arrays from the LQVec
    @param n the number of arrays to return
    @param i a pointer to the first index of the n arrays to return
    @returns A LQVec containing the n arrays
  */
  LQVec<T> extract(const size_t n, const size_t *i) const;
  /*! Return a collection of arrays from the LQVec
    @param idx a reference to an LQVec containing the to-be-returned indices
    @returns A LQVec containing the indicies indicated by idx
  */
  LQVec<T> extract(const ArrayVector<size_t>& idx) const;
  /*! Return a collection of arrays from the LQVec
    @param tfvec a reference to an ArrayVector<bool> with true for the to-be-returned indices
    @returns An LQVec containing the indicies indicated by tfvec
  */
  LQVec<T> extract(const ArrayVector<bool>& idx) const;
  LQVec<T> extract(const std::vector<bool>& idx) const;
  //! Extract the 3-vector with index `i`
  LQVec<T> get(const size_t i) const;
  //! Extract just the coordinates in units of the Reciprocal lattice (strip off the lattice information)
  ArrayVector<T> get_hkl() const;
  /*! Extract the coordinates in an orthonormal frame with its first axis, x,
  along a*, its second, y, perpendicular with y⋅b*>0 , and it's third forming
  the right-handed set z=x×y.
  */
  ArrayVector<double> get_xyz() const;
  //! Return the vector(s) expressed in units of the Direct lattice
  LDVec<double> star() const;

  //! Determine the scalar product between two vectors in the object.
  double dot(const size_t i, const size_t j) const;
  //! Determine the absolute length of a vector in the object.
  double norm(const size_t i) const { return sqrt(this->dot(i,i)); }
  //! Determine the cross product of two vectors in the object.
  LQVec<double> cross(const size_t i, const size_t j) const;

  LQVec<T>& operator+=(const LQVec<T>& av);
  LQVec<T>& operator-=(const LQVec<T>& av);
  LQVec<T>& operator+=(const T& av);
  LQVec<T>& operator-=(const T& av);
  LQVec<T>& operator*=(const T& av);
  LQVec<T>& operator/=(const T& av);
  // LQVec<T> operator -();

  //! Verify that a second LQVec object is compatible for binary operations
  template<class R> AVSizeInfo consistency_check(const LQVec<R>& b) const {
    if (!(this->samelattice(b))) throw std::runtime_error("arithmetic between Lattice vectors requires they have the same lattice");
    return this->ArrayVector<T>::consistency_check(b);
  }
  template<class R> AVSizeInfo consistency_check(const LDVec<R>& b) const {
    if (!(this->starlattice(b))) throw std::runtime_error("Reciprocal and real space vectors must be in dual lattices.");
    return this->ArrayVector<T>::consistency_check(b);
  }
  //! Verify that an ArrayVector object is compatible for binary operations
  template<class R, template<class> class A,
    typename=typename std::enable_if<!std::is_base_of<LatVec,A<R>>::value && std::is_base_of<ArrayVector<R>,A<R>>::value >::type
    >
  AVSizeInfo consistency_check(const A<R>& b) const {
    return this->ArrayVector<T>::consistency_check(b); // b has no lattice, so nothing to check
  }
  //! Check whether a second LQVec is approximately the same as this object
  template<typename R> bool isapprox(const LQVec<R>& that) const { return (this->samelattice(that) && this->ArrayVector<T>::isapprox(that)); }
  //! Check whether two vectors in the object are approximately the same.
  bool isapprox(const size_t i, const size_t j) const { return this->ArrayVector<T>::isapprox(i,j);}
  //! Round all elements using std::round
  LQVec<int> round() const {return LQVec<int>(this->lattice, this->ArrayVector<T>::round()); }
  //! Find the floor of all elements using std::floor
  LQVec<int> floor() const {return LQVec<int>(this->lattice, this->ArrayVector<T>::floor()); }
  //! Find the ceiling of all elements using std::ceil
  LQVec<int> ceil() const {return LQVec<int>(this->lattice, this->ArrayVector<T>::ceil()); }
protected:
  void check_arrayvector(const int);
};


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
template<class T, class R> struct LatVecTraits{
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


#include "latvec.tpp"
#include "ldvec.tpp"
#include "lqvec.tpp"

#endif
