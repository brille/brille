#ifndef SUBSCRIPT_HPP
#define SUBSCRIPT_HPP

#include <algorithm>
#include <array>
#include <vector>
#include <cassert>
#include <iostream>
#include <tuple>
#include <numeric>

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
  SubIt(const holder& _sh)
  : _shape(_sh), _first(0)
  {
    size_t n = _shape.size();
    _inpt = holder(n, T(0));
    _sub  = holder(n, T(0));
    _fixed = std::vector<bool>(n, false);
  }
  SubIt(const holder& _sh, const holder& _in)
  : _shape(_sh), _inpt(_in)
  {
    assert(_shape.size() == _inpt.size());
    size_t n = _shape.size();
    _sub = holder(n, T(0));
    _fixed = std::vector<bool>(n, false);
    for (size_t i=0; i<n; ++i){
      _fixed[i] = is_fixed(_inpt[i], _shape[i]);
      _sub[i] = _fixed[i] ? _inpt[i] : T(0);
    }
    this->find_first();
  }
  SubIt(const holder& _sh, const holder& _in, const holder& _s, const std::vector<bool>& _f)
  : _shape(_sh), _inpt(_in), _sub(_s), _fixed(_f)
  {
    this->find_first();
  }

  SubIt(const SubIt<T>& o)
  : _shape(o._shape), _inpt(o._inpt), _sub(o._sub), _fixed(o._fixed), _first(o._first)
  {}
  SubIt(const SubIt<T>* o)
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
  size_t ndim() const {return _shape.size();}
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
    return SubIt<T>(_shape, _inpt, sub, _fixed);
  }
  SubIt<T> end() const {
    size_t n = this->ndim();
    holder val(n, T(0));
    for (size_t i=0; i<n; ++i) if (_fixed[i]) val[i] = _sub[i];
    if (_first < n) val[_first] = _shape[_first];
    return SubIt<T>(_shape, _inpt, val, _fixed);
  }
};

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
    size_t nd = _shape0.size();
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

  const holder& shape() const {return _shapeO;}
  const size_t ndim() const {return _shapeO.size();}
  //const SubIt<T>& itr() const {return _itr;}

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

  // const subs_t& operator*() const {return subs;}
  // const subs_t* operator->() const {return &subs;}
  // subs_t& operator*() {return subs;}
  // subs_t* operator->() {return &subs;}
  std::tuple<holder,holder,holder> operator*() const {return triple_subscripts();}
  BroadcastIt<T> begin() const {
    size_t n = this->ndim();
    holder s0(n,0), s1(n,0), sO(n,0);
    return BroadcastIt<T>(_shape0, _shape1, _shapeO, s0, s1, sO);
  }
  BroadcastIt<T> end() const {
    size_t n = this->ndim();
    holder s0(n,0), s1(n,0), sO(n,0);
    sO[0] = _shapeO[0];
    if (_shape0[0] > 1) s0[0] = sO[0]; // not used in ==
    if (_shape1[0] > 1) s1[0] = sO[0]; // not used in ==
    return BroadcastIt<T>(_shape0, _shape1, _shapeO, s0, s1, sO);
  }
private:
  // std::tuple<holder,holder,holder> triple_subscripts() const {
  //   holder o{*_itr};
  //   holder a(o), b(o);
  //   for (size_t i=0; i<o.size(); ++i) if (o[i]>0) {
  //     if (1==_shape0[i]) a[i] = 0;
  //     if (1==_shape1[i]) b[i] = 0;
  //   }
  //   return std::make_tuple(o,a,b);
  // }
  std::tuple<holder,holder,holder> triple_subscripts() const {
    return std::make_tuple(_subO, _sub0, _sub1);
  }
protected:
  const holder& outer() const {return _subO;}
};

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

template <class I>
I sub2lin(const std::vector<I>& sub, const std::vector<I>& str){
  I lin{0};
  for (size_t i=0; i<sub.size(); ++i) lin += sub[i]*str[i];
  return lin;
}

template <class I>
I offset_sub2lin(const std::vector<I>& off, const std::vector<I>& sub, const std::vector<I>& str){
  I lin{0};
  for (size_t i=0; i<sub.size(); ++i) lin += (off[i]+sub[i])*str[i];
  return lin;
}

#endif
