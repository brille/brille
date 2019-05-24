#ifndef _LATVEC_CLASS_H_
#define _LATVEC_CLASS_H_

#include <math.h>
#include <assert.h>
#include <type_traits>
#include <typeinfo> // for std::bad_cast
#include <exception>
#include "lattice.h"
#include "arrayvector.h"
#include "linear_algebra.h"
#include "safealloc.h"


class LatVec{}; // base class to identify both LDVec and LQVec
template<typename T> class LDVec;
template<typename T> class LQVec;

template<typename T> class LDVec: public LatVec, public ArrayVector<T>{
  Direct lattice;
public:
  LDVec(const Reciprocal lat=Reciprocal(), const size_t n=0, const T *d=nullptr): ArrayVector<T>(3,n,d), lattice(lat){};
  LDVec(const Reciprocal lat, const ArrayVector<T>& vec, const int flag=1): ArrayVector<T>(vec), lattice(lat){ this->check_arrayvector(flag); };
  template<class R> LDVec(const LDVec<R>& vec): ArrayVector<T>(vec.numel(),vec.size(),vec.datapointer()), lattice(vec.get_lattice()) {};

  LDVec<T>& operator=(const LDVec<T>& other){
    if (this != &other){ // do nothing if called by, e.g., a = a;
      this->lattice = other.get_lattice();
      size_t m = other.numel();
      size_t n = other.size();
      // reuse the data-block if we can. otherwise, refresh/resize it
      if ( m !=this->numel() ) this->refresh(m,n);
      if ( n !=this->size()  ) this->resize(n);
      // copy-over the data (if it exists)
      if (other.data && m*n)
        for(size_t i=0; i<m*n; i++)
          this->data[i] = other.data[i];
    }
    return *this;
  };
  const LDVec<T> operator[](const size_t i) const{
    bool isok = i < this->size();
    LDVec<T> out(this->lattice, isok ? 1u: 0u);
    if (isok) out = this->get(i);
    return out;
  }

  Direct get_lattice() const { return lattice; };
  template<typename R> bool samelattice(const LDVec<R> *vec) const { return lattice.issame(vec->get_lattice()); };
  template<typename R> bool samelattice(const LQVec<R> *vec) const { return false; };
  template<typename R> bool starlattice(const LDVec<R> *vec) const { return false; };
  template<typename R> bool starlattice(const LQVec<R> *vec) const { return lattice.isstar(vec->get_lattice()); };
  template<typename R> bool samelattice(const LDVec<R> &vec) const { return lattice.issame(vec.get_lattice()); };
  template<typename R> bool samelattice(const LQVec<R> &vec) const { return false; };
  template<typename R> bool starlattice(const LDVec<R> &vec) const { return false; };
  template<typename R> bool starlattice(const LQVec<R> &vec) const { return lattice.isstar(vec.get_lattice()); };
  LDVec get(const size_t i) const;
  ArrayVector<T> get_hkl() const;
  ArrayVector<double> get_xyz() const;
  LQVec<double> star() const;

  double dot(const size_t i, const size_t j) const;
  double norm(const size_t i) const { return sqrt(this->dot(i,i)); };
  LDVec<double> cross(const size_t i, const size_t j) const;

  LDVec<T>& operator-=(const LDVec<T>& av);
  LDVec<T>& operator+=(const LDVec<T>& av);

  LDVec<T>& operator+=(const T& av);
  LDVec<T>& operator-=(const T& av);
  LDVec<T>& operator*=(const T& av);
  LDVec<T>& operator/=(const T& av);

  LDVec<T> operator -();
  template<class R> AVSizeInfo consistency_check(const LDVec<R>& b) const {
    if (!(this->samelattice(b))) throw std::runtime_error("arithmetic between Lattice vectors requires they have the same lattice");
    return this->ArrayVector<T>::consistency_check(b);
  };
  template<class R, template<class> class A,
    typename=typename std::enable_if<!std::is_base_of<LatVec,A<R>>::value && std::is_base_of<ArrayVector<R>,A<R>>::value >::type
    >
  AVSizeInfo consistency_check(const A<R>& b) const {
    return this->ArrayVector<T>::consistency_check(b); // b has no lattice, so nothing to check
  };
  template<typename R> bool isapprox(const LDVec<R>& that){ return (this->samelattice(that) && this->ArrayVector<T>::isapprox(that)); };
  bool isapprox(const size_t i, const size_t j) const { return this->ArrayVector<T>::isapprox(i,j);};
protected:
  void check_arrayvector(const int);
};

template<typename T> class LQVec:  public LatVec, public ArrayVector<T>{
  Reciprocal lattice;
public:
  LQVec(const Reciprocal lat=Reciprocal(), const size_t n=0, const T *d=nullptr): ArrayVector<T>(3,n,d), lattice(lat){};
  LQVec(const Reciprocal lat, const ArrayVector<T>& vec, const int flag=1): ArrayVector<T>(vec), lattice(lat){  this->check_arrayvector(flag); };
  template<class R>  LQVec(const LQVec<R>& vec): ArrayVector<T>(vec.numel(),vec.size(),vec.datapointer()), lattice(vec.get_lattice()) {};

  LQVec<T>& operator=(const LQVec<T>& other){
    if (this != &other){ // do nothing if called by, e.g., a = a;
      this->lattice = other.get_lattice();
      size_t m = other.numel();
      size_t n = other.size();
      // reuse the data-block if we can. otherwise, refresh/resize it
      if ( m !=this->numel() ) this->refresh(m,n);
      if ( n !=this->size()  ) this->resize(n);
      // copy-over the data (if it exists)
      if (other.data && m*n)
        for(size_t i=0; i<m*n; i++)
          this->data[i] = other.data[i];
    }
    return *this;
  };
  const LQVec<T> operator[](const size_t i) const{
    bool isok = i < this->size();
    LQVec<T> out(this->lattice, isok ? 1u: 0u);
    if (isok) out = this->get(i);
    return out;
  }

  Reciprocal get_lattice() const { return lattice; };
  template<typename R> bool samelattice(const LQVec<R> *vec) const { return lattice.issame(vec->get_lattice()); };
  template<typename R> bool samelattice(const LDVec<R> *vec) const { return false; };
  template<typename R> bool starlattice(const LQVec<R> *vec) const { return false; };
  template<typename R> bool starlattice(const LDVec<R> *vec) const { return lattice.isstar(vec->get_lattice()); };
  template<typename R> bool samelattice(const LQVec<R> &vec) const { return lattice.issame(vec.get_lattice()); };
  template<typename R> bool samelattice(const LDVec<R> &vec) const { return false; };
  template<typename R> bool starlattice(const LQVec<R> &vec) const { return false; };
  template<typename R> bool starlattice(const LDVec<R> &vec) const { return lattice.isstar(vec.get_lattice()); };
  LQVec get(const size_t i) const;
  ArrayVector<T> get_hkl() const;
  ArrayVector<double> get_xyz() const;
  LDVec<double> star() const;

  double dot(const size_t i, const size_t j) const;
  double norm(const size_t i) const { return sqrt(this->dot(i,i)); };
  LQVec<double> cross(const size_t i, const size_t j) const;

  LQVec<T>& operator+=(const LQVec<T>& av);
  LQVec<T>& operator-=(const LQVec<T>& av);
  LQVec<T>& operator+=(const T& av);
  LQVec<T>& operator-=(const T& av);
  LQVec<T>& operator*=(const T& av);
  LQVec<T>& operator/=(const T& av);
  LQVec<T> operator -();

  template<class R> AVSizeInfo consistency_check(const LQVec<R>& b) const {
    if (!(this->samelattice(b))) throw std::runtime_error("arithmetic between Lattice vectors requires they have the same lattice");
    return this->ArrayVector<T>::consistency_check(b);
  };
  template<class R, template<class> class A,
    typename=typename std::enable_if<!std::is_base_of<LatVec,A<R>>::value && std::is_base_of<ArrayVector<R>,A<R>>::value >::type
    >
  AVSizeInfo consistency_check(const A<R>& b) const {
    return this->ArrayVector<T>::consistency_check(b); // b has no lattice, so nothing to check
  };
  template<typename R> bool isapprox(const LQVec<R>& that) const { return (this->samelattice(that) && this->ArrayVector<T>::isapprox(that)); };
  bool isapprox(const size_t i, const size_t j) const { return this->ArrayVector<T>::isapprox(i,j);};
protected:
  void check_arrayvector(const int);
};



// extend the lattice traits structs
// template<> template<typename T> struct LatticeTraits<LDVec<T>>{
template<class T> struct LatticeTraits<LDVec<T>>{
  using type = Direct;
  using star = Reciprocal;
};
// template<> template<typename T> struct LatticeTraits<LQVec<T>>{
template<class T> struct LatticeTraits<LQVec<T>>{
  using type = Reciprocal;
  using star = Direct;
};
// create new lattice vectors traits
template<class T, class R> struct LatVecTraits{
  using type = void;
  using star = void;
};
// template<> template <typename R,typename S> struct LatVecTraits<LDVec<R>,S>{
template <class R, class S> struct LatVecTraits<LDVec<R>,S>{
  using type = LDVec<S>;
  using star = LQVec<S>;
};
// template<> template <typename R,typename S> struct LatVecTraits<LQVec<R>,S>{
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


#include "latvec.hpp"
#include "ldvec.hpp"
#include "lqvec.hpp"

#endif
