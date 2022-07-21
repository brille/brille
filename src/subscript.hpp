/* This file is part of brille.

Copyright © 2020 Greg Tucker <greg.tucker@stfc.ac.uk>

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
#ifndef BRILLE_SUBSCRIPT_HPP_
#define BRILLE_SUBSCRIPT_HPP_
/*! \file
    \author Greg Tucker
    \brief Classes for implementing automatic for loops with subscripted indices
*/
#include <algorithm>
#include <array>
#include <utility>
#include <vector>
#include <cassert>
#include <iostream>
#include <tuple>
#include <numeric>
namespace brille {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class I1, class I2>
std::enable_if_t<std::is_unsigned<I1>::value, bool>
is_fixed(I1 i, I2 n){
  return i<n;
}
template<class I1, class I2>
std::enable_if_t<std::is_signed<I1>::value, bool>
is_fixed(I1 i, I2 n){
  if (i<0) return false;
  return static_cast<I2>(i)<n;
}
#else
/*! \brief A utility function for the subscript iterators to decide if an
           input index represents a fixed value for a dimension.

\param i The index provided for a dimension
\param n The size of the array along the same dimension which will be iterated
\return `true` if `i` represents a value which should not be varied

The subscript iterators SubIt and SubIt2 can be provided with two sets of
integers when constructed. The first set is the shape of the array to iterate
over and the (optional) second set is used to indicate fixed subscripts which
are not iterated.

Any valid subscript index into an array with the specified shape indicates that
the iterator should not vary along a dimension. Therefore if a dimension
*should* be iterated over the provided value of `i` must be greater or equal
to the size of the array along that dimension.
*/
template<class Integer1, class Integer2> bool is_fixed(Integer1 i, Integer2 n);
#endif


/*! \brief An N-dimensional subscript iterator class

From an input shape and, optionally, set of fixed indices an interator
is created which is suitable to use with the shorthand for loop notation.

The subscript incides are iterated in row-order starting from (0,0,0) and
ending with (N,0,0) for an (N,M,L) shaped array.
*/
template<class T> class SubIt {
public:
  typedef std::vector<T> holder;
  holder _shape;
  holder _inpt;
  holder _sub;
  std::vector<bool> _fixed;
  size_t _first;
private:
  // void find_first() {
  //   // find the first non-fixed index (or the length of the _fixed vector)
  //   auto fitr = std::find(_fixed.begin(), _fixed.end(), false);
  //   if (fitr == _fixed.end())
  //     throw std::runtime_error("The input subscripts have fixed all dimensions!");
  //   _first = std::distance(_fixed.begin(), fitr);
  // }
  void find_first() {
    // find the first non-fixed index or the length of the _fixed vector
    _first = _fixed.size();
    for (size_t i=_first; i-- > 0;) if (!_fixed[i]) _first=i;
    if (_first == _fixed.size())
      throw std::runtime_error("The input subscripts have fixed all dimensions!");
  }
public:
  explicit SubIt()
  : _shape({0}), _inpt({0}), _sub({0}), _fixed({false}), _first(0)
  {}
  explicit SubIt(const holder& _sh)
  : _shape(_sh), _first(0)
  {
    size_t n = _shape.size();
    _inpt = holder(n, T(0));
    _sub  = holder(n, T(0));
    _fixed = std::vector<bool>(n, false);
  }
  SubIt(const holder& _sh, const holder& _in)
  : _shape(_sh), _inpt(_in), _fixed(_sh.size(), false), _first(_sh.size())
  {
    assert(_shape.size() == _inpt.size());
    size_t n = _shape.size();
    _sub = holder(n, T(0));
    for (size_t i=0; i<n; ++i){
      _fixed[i] = is_fixed(_inpt[i], _shape[i]);
      _sub[i] = _fixed[i] ? _inpt[i] : T(0);
    }
    this->find_first();
  }
  SubIt(const holder& _sh, const holder& _in, const holder& _s, std::vector<bool> _f)
  : _shape(_sh), _inpt(_in), _sub(_s), _fixed(std::move(_f)), _first(_fixed.size())
  {
    this->find_first();
  }
  SubIt(const holder& _sh, const holder& _in, holder && _s, std::vector<bool> _f)
      : _shape(_sh), _inpt(_in), _sub(_s), _fixed(std::move(_f)), _first(_fixed.size())
  {
    this->find_first();
  }


  SubIt(const SubIt<T>& o)
  : _shape(o._shape), _inpt(o._inpt), _sub(o._sub), _fixed(o._fixed), _first(o._first)
  {}
  explicit SubIt(const SubIt<T>* o)
  : _shape(o->_shape), _inpt(o->_inpt), _sub(o->_sub), _fixed(o->_fixed), _first(o->first)
  {}

  SubIt& operator=(const SubIt<T>& o){
    _shape = o._shape;
    _inpt = o._inpt;
    _sub = o._sub;
    _fixed = o._fixed;
    _first = o._first;
    return *this;
  }

  const holder& shape() const {return _shape;}
  [[nodiscard]] size_t ndim() const {return _shape.size();}
  SubIt& operator++(){
    size_t n = this->ndim();
    for (size_t dim=n; dim-->0; ) if (!_fixed[dim]) {
      if (dim > _first && _sub[dim]+1 == _shape[dim]){
        _sub[dim] = 0u;
      } else {
        ++_sub[dim];
        break;
      }
    }
    return *this;
  }
  bool operator==(const SubIt<T>& other) const {
    size_t n = this->ndim();
    if (other.ndim() != n) return false;
    bool equal{true};
    for (size_t i=0; i<n; ++i) equal &= _sub[i] == other._sub[i];
    return equal;
  }
  bool operator!=(const SubIt<T>& other) const { return !(*this == other); }
  const holder& operator*() const {return _sub;}
  const holder* operator->() const {return &_sub;}
  holder& operator*() {return _sub;}
  holder* operator->() {return &_sub;}
  SubIt<T> begin() const {
    // return SubIt<T>(_shape, _inpt, _sub, _fixed);
    size_t n = this->ndim();
    holder sub(n,T(0));
    for (size_t i=0; i<sub.size(); ++i) if (_fixed[i]) sub[i] = _inpt[i];
    return {_shape, _inpt, std::move(sub), _fixed};
  }
  SubIt<T> end() const {
    size_t n = this->ndim();
    holder val(n, T(0));
    for (size_t i=0; i<n; ++i) if (_fixed[i]) val[i] = _sub[i];
    if (_first < n) val[_first] = _shape[_first];
    return {_shape, _inpt, std::move(val), _fixed};
  }
};

/*! \brief An 2-dimensional subscript iterator class

From an input shape and, optionally, set of fixed indices an interator
is created which is suitable to use with the shorthand for loop notation.

The subscript incides are iterated in row-order starting from (0,0) and
ending with (N,0) for an (N,M) shaped array.
*/
template<class T> class SubIt2 {
public:
  typedef std::array<T,2> holder;
  holder _shape;
  holder _inpt;
  holder _sub;
  std::array<bool,2> _fixed;
  size_t _first;
private:
  void find_first() {
    // find the first non-fixed index or the length of the _fixed vector
    _first = _fixed.size();
    for (size_t i=_first; i-- > 0;) if (!_fixed[i]) _first=i;
    if (_first == _fixed.size())
      throw std::runtime_error("The input subscripts have fixed all dimensions!");
  }
public:
  explicit SubIt2()
  : _shape({0,0}), _inpt({0,0}), _sub({0,0}), _fixed({false,false}), _first(0)
  {}
  explicit SubIt2(const holder& _sh)
  : _shape(_sh), _fixed{false, false}, _first(0)
  {
    _inpt = holder({T(0), T(0)});
    _sub  = holder({T(0), T(0)});
  }
  SubIt2(const holder& _sh, const holder& _in)
  : _shape(_sh), _inpt(_in), _fixed({false, false}), _first(2u)
  {
    _sub  = holder({T(0), T(0)});
    for (size_t i=0; i<2; ++i){
      _fixed[i] = is_fixed(_inpt[i], _shape[i]);
      _sub[i] = _fixed[i] ? _inpt[i] : T(0);
    }
    this->find_first();
  }
  SubIt2(const holder& _sh, const holder& _in, const holder& _s, const std::array<bool,2>& _f)
  : _shape(_sh), _inpt(_in), _sub(_s), _fixed(_f), _first(_f.size())
  {
    this->find_first();
  }
  SubIt2(const holder & _sh, const holder & _in, holder && _s, const std::array<bool,2> & _f)
  : _shape(_sh), _inpt(_in), _sub(std::move(_s)), _fixed(_f), _first(_f.size())
  {
    this->find_first();
  }

  SubIt2(const SubIt2<T>& o)
  : _shape(o._shape), _inpt(o._inpt), _sub(o._sub), _fixed(o._fixed), _first(o._first)
  {}
  explicit SubIt2(const SubIt2<T>* o)
  : _shape(o->_shape), _inpt(o->_inpt), _sub(o->_sub), _fixed(o->_fixed), _first(o->first)
  {}

  SubIt2& operator=(const SubIt2<T>& o){
    _shape = o._shape;
    _inpt = o._inpt;
    _sub = o._sub;
    _fixed = o._fixed;
    _first = o._first;
    return *this;
  }

  const holder& shape() const {return _shape;}
  [[nodiscard]] size_t ndim() const {return 2;}
  SubIt2& operator++(){
    size_t n = 2;
    for (size_t dim=n; dim-->0; ) if (!_fixed[dim]) {
      if (dim > _first && _sub[dim]+1 == _shape[dim]){
        _sub[dim] = 0u;
      } else {
        ++_sub[dim];
        break;
      }
    }
    return *this;
  }
  bool operator==(const SubIt2<T>& other) const {;
    return _sub[0] == other._sub[0] && _sub[1] == other._sub[1];
  }
  bool operator!=(const SubIt2<T>& other) const { return !(*this == other); }
  const holder& operator*() const {return _sub;}
  const holder* operator->() const {return &_sub;}
  holder& operator*() {return _sub;}
  holder* operator->() {return &_sub;}
  // std::tuple<T,T> operator*() const {return std::make_tuple(_sub[0], _sub[1])};
  SubIt2<T> begin() const {
    holder sub({T(0), T(0)});
    if (_fixed[0]) sub[0] = _inpt[0];
    if (_fixed[1]) sub[1] = _inpt[1];
    return SubIt2<T>(_shape, _inpt, std::move(sub), _fixed);
  }
  SubIt2<T> end() const {
    holder val({T(0), T(0)});
    if (_fixed[0]) val[0] = _sub[0];
    if (_fixed[1]) val[1] = _sub[1];
    if (_first < 2) val[_first] = _shape[_first];
    return SubIt2<T>(_shape, _inpt, std::move(val), _fixed);
  }
};

/*! \brief An N-dimensional broadcasting subscript iterator class

From an two input shapes an iterator is created which is suitable to use with
the shorthand for loop notation.

```
for (auto [outer_idx, a_idx, b_idx] : BroadcastIt(a_shape, b_shape)){
  …
}
```

If the shapes of A and B are consistent then they have the same number of
dimensions and either matching sizes along all dimensions or one has a
singleton dimension at any mismatch.

For consistent a_shape and b_shape (1,3,2) and (2,3,1), respectively the
iterator iterates over the values

  1.  [(0,0,0), (0,0,0), (0,0,0)]
  2.  [(0,0,1), (0,0,1), (0,0,0)]
  3.  [(0,1,0), (0,1,0), (0,1,0)]
  4.  [(0,1,1), (0,1,1), (0,1,0)]
  5.  [(0,2,0), (0,2,0), (0,2,0)]
  6.  [(0,2,1), (0,2,1), (0,2,0)]
  7.  [(1,0,0), (0,0,0), (1,0,0)]
  8.  [(1,0,1), (0,0,1), (1,0,0)]
  9.  [(1,1,0), (0,1,0), (1,1,0)]
  10. [(1,1,1), (0,1,1), (1,1,0)]
  11. [(1,2,0), (0,2,0), (1,2,0)]
  12. [(1,2,1), (0,2,1), (1,2,0)]
*/
template<class T> class BroadcastIt{
public:
  typedef std::vector<T> holder;
private:
  holder _shape0;
  holder _shape1;
  holder _shapeO;
  holder _sub0;
  holder _sub1;
  holder _subO;
public:
  //
  BroadcastIt(const holder& a, const holder& b)
  : _shape0(a), _shape1(b), _shapeO(a.size(),0), _sub0(a.size(),0), _sub1(a.size(),0), _subO(a.size(),0)
  {
    assert(_shape0.size() == _shape1.size());
    auto nd = _shape0.size();
    for (size_t i=0; i<nd; ++i)
      if (_shape0[i]!=_shape1[i] && _shape0[i]!=1 && _shape1[i]!=1){
        std::string msg = "Can not broadcast { ";
        for (auto x: _shape0) msg += std::to_string(x) + " ";
        msg += "} and { ";
        for (auto x: _shape1) msg += std::to_string(x) + " ";
        msg += "} to a common shape";
        throw std::runtime_error(msg);
      } else {
        _shapeO[i] = _shape0[i] < _shape1[i] ? _shape1[i] : _shape0[i];
      }
  }
  BroadcastIt(const holder& s0, const holder& s1, const holder & sO, const holder& i0, const holder& i1, const holder& iO)
  : _shape0(s0), _shape1(s1), _shapeO(sO), _sub0(i0), _sub1(i1), _subO(iO)
  {
  }
  BroadcastIt(const holder & s0, const holder & s1, const holder & sO,
              holder && i0, holder && i1, holder && iO)
  : _shape0(s0), _shape1(s1), _shapeO(sO),
    _sub0(std::move(i0)), _sub1(std::move(i1)), _subO(std::move(iO))
  {
  }

  holder shape() const {return _shapeO;}
  [[nodiscard]] size_t ndim() const {return _shapeO.size();}
  T size() const {
    T prod{1};
    for (const auto & x: _shapeO) prod *= x;
    return prod;
  }

  BroadcastIt<T>& operator++(){
    size_t n = this->ndim();
    for (size_t dim=n; dim-->0; ){
      if (dim > 0 && _subO[dim]+1 == _shapeO[dim]){
        _sub1[dim] = _sub0[dim] = _subO[dim] = 0u;
      } else {
        ++_subO[dim];
        if (_shape0[dim] > 1) _sub0[dim] = _subO[dim];
        if (_shape1[dim] > 1) _sub1[dim] = _subO[dim];
        break;
      }
    }
    return *this;
  }
  bool operator==(const BroadcastIt<T>& other) const {
    size_t n = this->ndim();
    if (other.ndim() != n) return false;
    bool equal{true};
    const holder& oO{other.outer()};
    for (size_t i=0; i<n; ++i) equal &= _subO[i] == oO[i];
    return equal;
  }
  bool operator!=(const BroadcastIt<T>& other) const {return !(*this==other);}

  std::tuple<holder,holder,holder> operator*() const {return triple_subscripts();}
  BroadcastIt<T> begin() const {
    size_t n = this->ndim();
    /* Something here causes a compiler warning */
    holder s0(n,0), s1(n,0), sO(n,0);
    return {_shape0, _shape1, _shapeO, std::move(s0), std::move(s1), std::move(sO)};
  }
  BroadcastIt<T> end() const {
    size_t n = this->ndim();
    holder s0(n,0), s1(n,0), sO(n,0);
    sO[0] = _shapeO[0];
    if (_shape0[0] > 1) s0[0] = sO[0]; // not used in ==
    if (_shape1[0] > 1) s1[0] = sO[0]; // not used in ==
    return {_shape0, _shape1, _shapeO, std::move(s0), std::move(s1), std::move(sO)};
  }
private:
  std::tuple<holder,holder,holder> triple_subscripts() const {
    return std::make_tuple(_subO, _sub0, _sub1);
  }
protected:
  const holder& outer() const {return _subO;}
};

/*! \brief A 2-dimensional broadcasting subscript iterator class

From an two input shapes an iterator is created which is suitable to use with
the shorthand for loop notation.

```
for (auto [outer_idx, a_idx, b_idx] : BroadcastIt(a_shape, b_shape)){
  …
}
```

If the shapes of A and B are consistent then they have the same number of
dimensions and either matching sizes along all dimensions or one has a
singleton dimension at any mismatch.

For consistent a_shape and b_shape (1,3) and (2,3), respectively the
iterator iterates over the values

  1. [(0,0), (0,0), (0,0)]
  2. [(0,1), (0,1), (0,1)]
  3. [(0,2), (0,2), (0,2)]
  4. [(1,0), (0,0), (1,0)]
  5. [(1,1), (0,1), (1,1)]
  6. [(1,2), (0,2), (1,2)]
*/
template<class T> class BroadcastIt2{
public:
  typedef std::array<T,2> holder;
private:
  holder _shape0;
  holder _shape1;
  holder _shapeO;
  holder _sub0;
  holder _sub1;
  holder _subO;
public:
  //
  BroadcastIt2(const holder& a, const holder& b)
  : _shape0(a), _shape1(b), _shapeO({0,0}), _sub0({0,0}), _sub1({0,0}), _subO({0,0})
  {
    size_t nd = 2;
    for (size_t i=0; i<nd; ++i)
      if (_shape0[i]!=_shape1[i] && _shape0[i]!=1 && _shape1[i]!=1){
        std::string msg = "Can not broadcast { ";
        for (auto x: _shape0) msg += std::to_string(x) + " ";
        msg += "} and { ";
        for (auto x: _shape1) msg += std::to_string(x) + " ";
        msg += "} to a common shape";
        throw std::runtime_error(msg);
      } else {
        _shapeO[i] = _shape0[i] < _shape1[i] ? _shape1[i] : _shape0[i];
      }
  }
  BroadcastIt2(const holder& s0, const holder& s1, const holder & sO, const holder& i0, const holder& i1, const holder& iO)
  : _shape0(s0), _shape1(s1), _shapeO(sO), _sub0(i0), _sub1(i1), _subO(iO)
  {
  }
  BroadcastIt2(const holder & s0, const holder & s1, const holder & sO, holder && i0, holder && i1, holder && iO)
      : _shape0(s0), _shape1(s1), _shapeO(sO), _sub0(std::move(i0)), _sub1(std::move(i1)), _subO(std::move(iO))
  {
  }

  holder shape() const {return _shapeO;}
  [[nodiscard]] size_t ndim() const {return 2;}
  T size() const {return _shapeO[0] * _shapeO[1];}
  //const SubIt<T>& itr() const {return _itr;}

  BroadcastIt2<T>& operator++(){
    size_t n = 2;
    for (size_t dim=n; dim-->0; ){
      if (dim > 0 && _subO[dim]+1 == _shapeO[dim]){
        _sub1[dim] = _sub0[dim] = _subO[dim] = 0u;
      } else {
        ++_subO[dim];
        if (_shape0[dim] > 1) _sub0[dim] = _subO[dim];
        if (_shape1[dim] > 1) _sub1[dim] = _subO[dim];
        break;
      }
    }
    return *this;
  }
  bool operator==(const BroadcastIt2<T>& other) const {
    const holder& oO{other.outer()};
    return (_subO[0]==oO[0] && _subO[1]==oO[1]);
  }
  bool operator!=(const BroadcastIt2<T>& other) const {return !(*this==other);}

  std::tuple<holder,holder,holder> operator*() const {
    return std::make_tuple(_subO, _sub0, _sub1);
  }
  BroadcastIt2<T> begin() const {
    holder s0({0,0}), s1({0,0}), sO({0,0});
    return {_shape0, _shape1, _shapeO, std::move(s0), std::move(s1), std::move(sO)};
  }
  BroadcastIt2<T> end() const {
    holder s0({0,0}), s1({0,0}), sO({0,0});
    sO[0] = _shapeO[0];
    if (_shape0[0] > 1) s0[0] = sO[0]; // not used in ==
    if (_shape1[0] > 1) s1[0] = sO[0]; // not used in ==
    return {_shape0, _shape1, _shapeO, std::move(s0), std::move(s1), std::move(sO)};
  }
protected:
  const holder& outer() const {return _subO;}
};


/*! \brief Return the subscript index for an element of an Array from its linear index

\param l       the linear index to the element of the Array
\param stride  the stride detailing how far apart consecutive subcript indexed
               elements are in their linear index along all dimensions of the
               related Array.
\returns the subscript indexes
*/
template <class I>
std::vector<I> lin2sub(I l, const std::vector<I>& stride){
  std::vector<I> sub;
  size_t ndim = stride.size();
  if (1 == ndim) sub.push_back(l);
  else if (1 < ndim) {
    sub.resize(ndim);
    if (stride[ndim-1] > stride[0])
    for (I i=ndim-1; i--; ){
      sub[i] = l/stride[i];
      l -= sub[i]*stride[i];
    }
    else
    for (I i=0; i<ndim; ++i){
      sub[i] = l/stride[i];
      l -= sub[i]*stride[i];
    }
  }
  return sub;
}
/*! \brief Return the subscript index for an element of an Array2 from its linear index

\param l       the linear index to the element of the Array2
\param stride  the stride detailing how far apart consecutive subcript indexed
               elements are in their linear index along both dimensions of the
               related Array2.
\returns the subscript indexes
*/
template <class I>
std::array<I,2> lin2sub(I l, const std::array<I,2>& stride){
  std::array<I,2> sub;
  if (stride[1] > stride[0]){
    sub[1] = l/stride[1];
    sub[0] = (l-sub[1]*stride[1])/stride[0];
  } else {
    sub[0] = l/stride[0];
    sub[1] = (l-sub[0]*stride[0])/stride[1];
  }
  return sub;
}

// template <class I>
// I sub2lin(const std::vector<I>& sub, const std::vector<I>& stride){
//   assert(sub.size() == stride.size());
// #if defined(__GNUC__) && (__GNUC__ < 9 || (__GNUC__ == 9 && __GNUC_MINOR__ <= 2))
//   // serial inner_product
//   return std::inner_product(sub.begin(), sub.end(), stride.begin(), I(0));
// #else
//   // parallelized inner_product
//   return std::transform_reduce(sub.begin(), sub.end(), stride.begin(), I(0));
// #endif
// }

/*! \brief Return a linear index from its subscript and Array stride

\param sub the subscript index
\param str the stride detailing how far apart consecutive subscript indexed
           elements are in their linear index along all dimensions of the
           related Array.
\returns the linear index
*/
template <class I>
I sub2lin(const std::vector<I>& sub, const std::vector<I>& str){
  I lin{0};
  for (size_t i=0; i<sub.size(); ++i) lin += sub[i]*str[i];
  return lin;
}
/*! \brief Return a linear index from its subscript and Array2 stride

\param sub the subscript index
\param str the stride detailing how far apart consecutive subscript indexed
           elements are in their linear index along both dimensions of the
           related Array2.
\returns the linear index
*/
template <class I>
I sub2lin(const std::array<I,2>& sub, const std::array<I,2>& str){
  return sub[0]*str[0] + sub[1]*str[1];
}
/*! \brief Return a linear index from its explicit subscript and Array2 stride

\param s0 the subscript index in the first dimension
\param s1 the subscript index in the second dimension
\param str the stride detailing how far apart consecutive subscript indexed
           elements are in their linear index along both dimensions of the
           related Array2.
\returns the linear index
*/
template <class I>
I sub2lin(const I s0, const I s1, const std::array<I,2>& str){
  return s0*str[0] + s1*str[1];
}

////////////////////////////////////////////////////////////////////////////////
// These were used when the Array class contained a subscripted offset vector //
// but are no longer needed now that the offset is in the linear index.       //
////////////////////////////////////////////////////////////////////////////////
// template <class I>
// I offset_sub2lin(const std::vector<I>& off, const std::vector<I>& sub, const std::vector<I>& str){
//   I lin{0};
//   for (size_t i=0; i<sub.size(); ++i) lin += (off[i]+sub[i])*str[i];
//   return lin;
// }
//
// template <class I>
// I offset_sub2lin(const std::array<I,2>& off, const std::array<I,2>& sub, const std::array<I,2>& str){
//   return (off[0]+sub[0])*str[0] + (off[1]+sub[1])*str[1];
// }
// template <class I>
// I offset_sub2lin(const std::array<I,2>& off, const I s0, const I s1, const std::array<I,2>& str){
//   return (off[0]+s0)*str[0] + (off[1]+s1)*str[1];
// }

} // end namespace brille
#endif
