#ifndef SHAREDSPAN_HPP
#define SHAREDSPAN_HPP

#include <memory>
// #include <span>
#include <cassert>

#include "subscript.hpp"

namespace brille {

template<class T>
class SharedSpan {
public:
  typedef unsigned ind_t;
  typedef std::vector<ind_t> shape_t;
  typedef std::shared_ptr<char> ref_t;
protected:
  T* _data;
  ind_t _size;
  shape_t _shape;
  shape_t _stride;
  ref_t _refs;
  bool _owns;
public:
  // empty initialization -- C++ doesn't 'own' the null pointer (for operator= below)
  explicit SharedSpan(): _data(nullptr), _size(0), _owns(false) {}
  // Non-owning wrapper around raw pointer and size (e.g., C++20 std::span<T>)
  SharedSpan(T* d, ind_t sz, bool o=false):
  _data(d), _size(sz), _refs(std::make_shared<char>(0)), _owns(o) {
    _shape = shape_t(1, _size);
    _stride = shape_t(1, 1);
    this->verify();
  }
  // Non-owning wrapper around a raw pointer with shape and stride
  SharedSpan(T* d, const shape_t& sh, const shape_t& st, bool o=false):
  _data(d), _shape(sh), _stride(st), _refs(std::make_shared<char>(0)), _owns(o) {
    _size = this->size_from_shape();
    this->verify();
  }
  // Non-initializing with size only (1-D)
  SharedSpan(size_t n): _shape(1,n), _stride(1,1), _owns(true) {
    this->construct();
    this->verify();
  }
  // Initializing with size only (1-D)
  SharedSpan(size_t n, T init): _shape(1,n), _stride(1,1), _owns(true) {
    this->construct(init);
    this->verify();
  }
  // Non-initializing with shape only (N-D)
  SharedSpan(const shape_t& sh): _shape(sh), _stride(1,1), _owns(true) {
    this->reset_stride();
    this->construct();
    this->verify();
  }
  // Initializing with shape only (N-D)
  SharedSpan(const shape_t& sh, T init): _shape(sh), _stride(1,1), _owns(true) {
    this->reset_stride();
    this->construct(init);
    this->verify();
  }
  // Non-initializing with shape and stride (N-D)
  SharedSpan(const shape_t& sh, const shape_t& st): _shape(sh), _stride(st), _owns(true) {
    this->construct();
    this->verify();
  }
  // Initializing with shape and stride (N-D)
  SharedSpan(const shape_t& sh, const shape_t& st, T init): _shape(sh), _stride(st), _owns(true) {
    this->construct(init);
    this->verify();
  }
  // // Copying std::vector<T> (1-D) ***** THIS CONFLICTS WITH // Non-initializing with shape only (N-D)
  // SharedSpan(const std::vector<T>& data): _shape({data.size()}), _stride({1}), _owns(true){
  //   this->construct(data);
  //   this->verify();
  // }
  // Copying std::vector<std::array<R,N>> (2-D)
  template<class R, size_t Nel>
  SharedSpan(const std::vector<std::array<R,Nel>>& data):
  _shape({static_cast<ind_t>(data.size()), static_cast<ind_t>(Nel)}),
  _stride({static_cast<ind_t>(Nel), 1}), _owns(true){
    this->construct(data);
    this->verify();
  }
  // // [Optionally decouping] Copy constructor
  // SharedSpan(const SharedSpan<T>& other, bool decouple=false):
  // _data(other._data), _size(other.size()), _shape(other.shape()), _stride(other.stride()),
  // _refs(other.refs()), _owns(other.owns()) {
  //   if (decouple){
  //     if (_data){
  //       T* dec = new T[_size]();
  //       std::copy(_data, _data+_size, dec);
  //       if (_owns && _refs.use_count()==1) delete[] _data;
  //       _data = dec;
  //     }
  //     _refs = std::make_shared<char>(0);
  //   }
  // }
  // Copying is a shallow-copy and shares resources
  SharedSpan(const SharedSpan<T>& other):
  _data(other._data), _size(other.size()), _shape(other.shape()),
  _stride(other.stride()), _refs(other.refs()), _owns(other.owns()) {}
  // destructor
  ~SharedSpan() {
    // bool del = _owns;
    // del &= _refs.use_count() == 1;
    // del &= _data != nullptr;
    // if (del) delete[] _data;
    if (_owns && _refs.use_count()==1 && _data != nullptr) delete[] _data;
  }
  // Assignment
  SharedSpan<T>& operator=(const SharedSpan<T>& other){
    if (this != &other){
      if (_owns){
        T* old_data = _data;
        ref_t old_refs = _refs;
        _data = other._data;
        _refs = other.refs();
        _owns = other.owns();
        if (old_refs.use_count()==1) delete[] old_data;
      } else {
        _data = other._data;
        _refs = other.refs();
        _owns = other.owns();
      }
      _size = other.size();
      _shape = other.shape();
      _stride = other.stride();
    }
    return *this;
  }
  // Move constructor, also shallow
  SharedSpan(SharedSpan<T>&& other) noexcept
  : _data(std::exchange(other._data, nullptr))
  {
    std::swap(_size, other._size);
    std::swap(_shape, other._shape);
    std::swap(_shape, other._shape);
    std::swap(_stride, other._stride);
    std::swap(_refs, other._refs);
    std::swap(_owns, other._owns);
  }
  // Move assignment
  SharedSpan<T>& operator=(SharedSpan<T>&& other) noexcept
  {
    std::swap(_data, other._data);
    std::swap(_size, other._size);
    std::swap(_shape, other._shape);
    std::swap(_shape, other._shape);
    std::swap(_stride, other._stride);
    std::swap(_refs, other._refs);
    std::swap(_owns, other._owns);
    return *this;
  }

  // modifiers
  SharedSpan<T>& _decouple() { // actual decouple without checks
    T* new_data = new T[_size]();
    std::copy(_data, _data+_size, new_data);
    if (_owns){
      ref_t old_ref = _refs;
      _refs = std::make_shared<char>(0);
      if (old_ref.use_count() == 1) delete[] _data;
    } else {
      _refs = std::make_shared<char>(0);
    }
    _owns = true;
    _data = new_data;
    return *this;
  }
  SharedSpan<T>& decouple(){
    // only do anything if brille doesn't own the data or it is being shared internally
    if (!_owns || _refs.use_count() > 1) this->_decouple();
    return *this;
  }
  // accessors
  // std::span<T> span() const {return std::span(_data,_size);}
  T* data() {return _data;}
  const T* data() const {return _data;}
  size_t size() const {return _size;}
  shape_t shape() const {return _shape;}
  shape_t stride() const {return _stride;}
  shape_t cstride() const {
    shape_t cs(_stride);
    for (auto& s: cs) s *= sizeof(T);
    return cs;
  }
  bool is_column_ordered() const {
    return std::is_sorted(_stride.begin(), _stride.end(), [](ind_t a, ind_t b){return a <= b;});
  }
  bool is_row_ordered() const {
    return std::is_sorted(_stride.begin(), _stride.end(), [](ind_t a, ind_t b){return a >= b;});
  }
  bool is_contiguous() const {
    shape_t expected(1, 1);
    if (this->is_row_ordered()){
      for (auto s = _shape.rbegin(); s!=_shape.rend(); ++s)
        expected.push_back(expected.back()*(*s));
      expected.pop_back();
      return std::equal(_stride.begin(), _stride.end(), expected.rbegin());
    }
    if (this->is_column_ordered()){
      for (auto s : _shape) expected.push_back(expected.back()*s);
      // expected has one too many elements
      return std::equal(_stride.begin(), _stride.end(), expected.begin());
    }
    return false;
  }
  ref_t refs() const {return _refs;}
  bool owns() const {return _owns;}
  T& operator[](size_t i){return _data[i];}
  const T& operator[](size_t i) const {return _data[i];}
  T& operator[](const shape_t& sub){
    assert(sub.size() == _stride.size());
    return _data[sub2lin(sub, _stride)];
  }
  const T& operator[](const shape_t& sub) const {
    assert(sub.size() == _stride.size());
    return _data[sub2lin(sub, _stride)];
  }

  // modifiers
  size_t reshape(const shape_t& ns){
    size_t n = std::reduce(ns.begin(), ns.end(), 1, std::multiplies<ind_t>());
    assert(_size == n);
    _shape = ns;
    this->reset_stride();
    return n;
  }
  size_t resize(const shape_t& ns, T init=T(0)){
    size_t n = std::reduce(ns.begin(), ns.end(), 1, std::multiplies<ind_t>());
    // reuse the data container if we can
    if (_size == n) return this->reshape(ns);
    // we must preserve dimensionality if not reshaping
    assert(_shape.size() == ns.size());
    std::vector<bool> shrinking, growing;
    for (size_t i=1; i<ns.size(); ++i){
      shrinking.push_back(ns[i] <= _shape[i]);
      growing.push_back(ns[i] >= _shape[i]);
    }
    bool shrink = std::all_of(shrinking.begin(), shrinking.end(), [](bool a){return a;});
    bool grow = !shrink && std::all_of(growing.begin(), growing.end(), [](bool a){return a;});
    if (!shrink && !grow){
      std::string msg = "This method is not able to resize from { ";
      for (auto x : _shape) msg += std::to_string(x) + " ";
      msg += "} to { ";
      for (auto x : ns) msg += std::to_string(x) + " ";
      msg += "}";
      throw std::runtime_error(msg);
    }
    // make the new array
    T* nd = new T[n]();
    // fill the whole thing instead of trying to figure out difference indices
    if (grow) std::fill(nd, nd+n, init);
    // loop over the subscripted incides which fit inside both arrays
    shape_t inside = shrink ? ns : _shape;
    // store the original stride for subscript conversion
    shape_t ot = _stride;
    // update the shape and stride (maintaining row or column order)
    _shape = ns;
    this->reset_stride();
    // copy into the new array if the old array isn't the null pointer
    if (_data)
    for (auto x : SubIt(inside)) nd[sub2lin(x, _stride)] = _data[sub2lin(x, ot)];
    // delete the old array if it's owned by C++ and no other objects hold its reference
    if (_owns && _refs.use_count()==1 && _data) delete[] _data;
    // set the new data and make a new reference for it
    _data = nd;
    _size = n;
    _refs = std::make_shared<char>(0);
    return n;
  }
T* ptr(const shape_t& idx){
  return _data + sub2lin(idx, _stride);
}
const T* ptr(const shape_t& idx) const {
  return _data + sub2lin(idx, _stride);
}
private:
  void construct() {
    ind_t n = this->size_from_shape();
    if (n>0){
      _refs = std::make_shared<char>(0);
      _data = new T[n]();
      _size = n;
    }
  }
  void construct(T init) {
    this->construct();
    if (_size > 0) std::fill(_data, _data+_size, init);
  }
  void construct(const std::vector<T>& d){
    this->construct();
    if (_size > 0) for (ind_t i=0; i<_size; ++i) _data[i] = d[i];
  }
  // template<size_t Nel>
  // void construct(const std::vector<std::array<T,Nel>>& d){
  //   this->construct();
  //   if (_size > 0){
  //     ind_t cnt{0};
  //     for (auto a: d) for (auto x: a) _data[cnt++] = x;
  //     if (cnt != _size)
  //       throw std::logic_error("Something has gone wrong in construction");
  //   }
  // }
  template<class R, size_t Nel>
  void construct(const std::vector<std::array<R,Nel>>& d){
    this->construct();
    if (_size > 0){
      ind_t cnt{0};
      for (auto a: d) for (auto x: a) _data[cnt++] = static_cast<T>(x);
      if (cnt != _size)
        throw std::logic_error("Something has gone wrong in construction");
    }
  }
  void reset_stride() {
    shape_t nt(_shape.size(), 1u);
    if (_stride.front() < _stride.back()){
      for (size_t i=1; i<nt.size(); ++i) nt[i] = nt[i-1]*_shape[i-1];
    } else {
      for (size_t i=nt.size()-1; i--; ) nt[i] = nt[i+1]*_shape[i+1];
    }
    _stride = nt;
  }
  void verify() const {
    assert(_shape.size() == _stride.size());
    assert(_size == this->size_from_shape());
  }
  ind_t size_from_shape() const {
    size_t sz = std::reduce(_shape.begin(), _shape.end(), ind_t(1), std::multiplies<ind_t>());
    return static_cast<ind_t>(sz);
  }
};

} // end namespace brille

template<class T>
brille::SharedSpan<T> operator-(const brille::SharedSpan<T>& x){
  brille::SharedSpan<T> nx(x.shape(), x.stride());
  for (typename brille::SharedSpan<T>::ind_t i=0; i<x.size(); ++i)
    nx[i] = -x[i];
  return nx;
}
#endif
