#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <functional>
#include <algorithm>
#include <numeric>
#include <vector>
#include <array>
// #include <cmath>
#include <math.h>
#include <cassert>
#include <iostream>

#include "sharedspan.hpp"
#include "subscript.hpp"
#include "utilities.hpp"
#include "comparisons.hpp"
#include "approx.hpp"
#include "types.hpp"

namespace brille {

/*! \brief A multidimensional shared data array with operator overloads

The `Array<T>`` object holds a multidimensional shared data object plus
information necessary to index into that array with an offset, different shape,
and/or different stride than the underlying data contains.

The underlying `SharedSpan<T>` object is itself a multidimensional array holding
a reference to the shared heap-based raw array as well as size, shape, and
stride information allowing for multidimensional subscripted indexing.

The `Array<T>` object provides overloads for basic operations, a number of
utility methods, and methods to create offset and/or reduced-range mutable
or immutable views into the shared array.
*/
template<class T>
class Array{
public:
  using data_t = brille::SharedSpan<T>;
private:
  bool  _mutable;
  shape_t _offset; // different from {0,…,0} if a view or offset-strided access
  shape_t _shape; // different from _data.shape() if a view or reshaped
  shape_t _stride; // different from _data.stride() if reshaped or strided access
protected: // to allow inheriting classes to directly modify their data
  data_t _data;
public:
  // accessors
  bool  ismutable(void) const {return _mutable;}
  bool  isconst(void) const {return !_mutable;}
  ind_t ndim(void) const {return _shape.size();}
  ind_t numel(void) const {
    ind_t nel = static_cast<ind_t>(std::accumulate(_shape.begin(), _shape.end(), 1, std::multiplies<>()));
    return this->ndim() ? nel : 0u;
  }
  shape_t offset(void) const {return _offset;}
  // ind_t size(void) const {
  //   return std::reduce(_shape.begin(), _shape.end(), ind_t(1), std::multiplies<ind_t>());
  // }
  ind_t size(const ind_t dim) const {
    assert(dim < _shape.size());
    return _shape[dim];
  }
  shape_t shape(void) const {return _shape;}
  shape_t stride(void) const {return _stride;}
  shape_t cstrides() const {
    shape_t cs(_stride);
    for (auto& s: cs) s *= sizeof(T);
    return cs;
  }
  data_t data(void) const {return _data;}
  bool is_row_ordered() const {return _data.is_row_ordered();}
  bool is_column_ordered() const {return _data.is_column_ordered();}
  bool is_contiguous() const {return _data.is_contiguous();}
  // "normal" constructors
  explicit Array(){};
  Array(
    const shape_t& shape,
    const T init=T(0),
    const bool mut=true
  ):  _mutable(mut), _offset(shape_t(shape.size(), 0)), _shape(shape), _data(shape, init){
    this->set_stride();
    this->init_check();
  }
  Array(const data_t& data, const bool mut=true)
  :_mutable(mut), _offset(data.shape().size(), 0), _shape(data.shape()),
  _stride(data.stride()), _data(data)
  {
    this->init_check();
  }
  Array(const data_t& data, const shape_t& shape, const bool mut=true)
  :_mutable(mut), _offset(shape.size(), 0), _shape(shape), _data(data)
  {
    this->set_stride();
    this->init_check();
  }
  Array(const shape_t& shape, const shape_t& stride, const T init, bool mut=true)
  :_mutable(mut), _offset(shape.size(), 0), _shape(shape), _stride(stride),
  _data(shape, init)
  {
    this->init_check();
  }
  Array(const std::vector<T>& data)
  :_mutable(true), _offset({0}),
  _shape({static_cast<ind_t>(data.size())}), _stride({1}), _data(data)
  {
    this->init_check();
  }
  template<class R, size_t Nel>
  Array(const std::vector<std::array<R,Nel>>& data)
  :_mutable(true), _offset({0,0}),
  _shape({static_cast<ind_t>(data.size()), Nel}), _stride({Nel,1}), _data(data)
  {
    this->init_check();
  }
  // // reshape constructor?
  // Array(const shape_t& shape, const data_t& data, const bool mut=true):
  // _mutable(mut), _offset(shape_t(shape.size(), 0)), _shape(shape),
  // _stride(data.stride()), _data(data) {
  //   this->init_check();
  // }
  // view constructor
  Array(const shape_t& offset, const shape_t& shape, const data_t& data, const bool mut=true):
  _mutable(mut), _offset(offset), _shape(shape), _stride(data.stride()), _data(data) {
    this->init_check();
  }
  // strided view constructor
  Array(const shape_t& offset, const shape_t& shape, const shape_t& stride, const data_t& data, const bool mut=true):
  _mutable(mut), _offset(offset), _shape(shape), _stride(stride), _data(data) {
    this->init_check();
  }
  // // explicit copy constructing breaks the underying data connection
  // Array(const Array<T>& other, bool m=false, bool decouple=false):
  // _mutable(m), _offset(other.offset()), _shape(other.shape()),
  // _stride(other.stride()) {
  //   if (decouple && other.data().size() != this->numel()){
  //     // we want to decouple but won't use the whole other array
  //     _offset = shape_t(this->ndim(), 0u);
  //     _data = data_t(_shape,_stride);
  //     for (auto x : SubIt(_shape)) _data[x] = other[x]; //_offset==0, so _data[x] is ok.
  //   } else {
  //     // either don't decouple or decouple the whole array
  //     _data = data_t(other.data(), decouple);
  //   }
  // }
  // type casting requires a copy -- we take an extra bool to support typecasting
  // copy construction with the same syntax as non-type-casting copy construction
  template<class R>
  Array(const Array<R>& other)
  : _mutable(true), _offset(other.ndim(),0), _shape(other.shape())
  {
    this->set_stride();
    ind_t count = this->init_data();
    for (ind_t i=0; i<count; ++i){
      auto osub = this->l2s_d(i);
      _data[i] = static_cast<T>(other[osub]);
    }
  }
  Array(ind_t sh0, ind_t sh1)
  : _mutable(true), _offset({0,0}), _shape({sh0,sh1}), _data(_shape)
  {
    this->set_stride();
    this->init_check();
  }
  Array(ind_t sh0, ind_t sh1, T init)
  : _mutable(true), _offset({0,0}), _shape({sh0,sh1}), _data(_shape, init)
  {
    this->set_stride();
    this->init_check();
  }
  template<class R>
  Array<T>& operator=(const Array<R>& other){
    _mutable = true;
    _shape = other.shape();
    _offset = shape_t(_shape.size(), 0);
    this->set_stride();
    ind_t count = this->init_data();
    for(ind_t i=0; i<count; ++i){
      auto osub = this->l2s_d(i);
      _data[i] = static_cast<T>(other[osub]);
    }
    return *this;
  }

  // modifiers
  bool make_mutable() {_mutable = true; return _mutable;}
  bool make_immutable() {_mutable = false; return _mutable;}
  Array<T>& decouple(const shape_t& new_stride) {
    // create the new heap array
    data_t new_data(_shape, new_stride);
    // and copy data over, using old offset and stride to (possibly) modify subscripts
    for (auto x : SubIt(_shape)) new_data[x] = _data[this->s2s_d(x)];
    // the new offset is zero:
    _offset = shape_t(this->ndim(), 0u);
    // the shape stays the same
    // the new _stride is given by the input:
    _stride = new_stride;
    // the new _data can now be assigned
    _data = new_data; // this releases the reference to the old data and deletes it if no references remain
    // and the decoupled copy shoud (probably) be mutable
    _mutable = true;
    return *this;
  }
  Array<T>& decouple() { return this->decouple(_stride); }

  // data accessors
        T& operator[](ind_t    lin)       {return _data[this->l2s_d(lin)];}
  const T& operator[](ind_t    lin) const {return _data[this->l2s_d(lin)];}

        T& operator[](shape_t& sub)       {return _data[this->s2s_d(sub)];}
  const T& operator[](shape_t& sub) const {return _data[this->s2s_d(sub)];}
protected: // so inherited classes can calculate subscript indexes into their data
  // find the subscript index for linear index l and then add our offset
  shape_t l2s_d(ind_t l) const {
    shape_t sub = lin2sub(l, _stride);
    for (size_t i=0; i<sub.size(); ++i) sub[i] += _offset[i];
    return sub;
  }
  // add our offset to the provided subscript index
  // strided views would need more work here:
  shape_t s2s_d(const shape_t& s) const {
    assert(s.size() == _offset.size());
    shape_t sub;
    for (size_t i=0; i<s.size(); ++i) sub.push_back(s[i]+_offset[i]);
    return sub;
  }
private:
  // construct the underlying shared data object
  ind_t init_data(const T init=T(0)){
    _data = SharedSpan<T>(_shape, _stride, init);
    return _data.size();
  }
  void set_stride(void){
    _stride.clear();
    _stride.push_back(1);
    for (auto s = _shape.rbegin(); s!=_shape.rend(); ++s)
      _stride.push_back(_stride.back()*(*s));
    _stride.pop_back();
    std::reverse(_stride.begin(), _stride.end());
  }
  void init_check(void){
    if (_shape.size() != _offset.size() || _shape.size() != _stride.size())
    throw std::runtime_error("Attempting to construct Array with inconsistent offset, size, and strides");
    shape_t sh;
    for (size_t i=0; i<_shape.size(); ++i) sh.push_back(_offset[i]+_shape[i]);
    auto d_shape = _data.shape();
    if (d_shape.size() != sh.size())
      throw std::runtime_error("Unequal _data.shape() and _shape dimensions!");
    for (size_t i=0; i<sh.size(); ++i) if (d_shape[i] < sh[i]) {
      std::string msg = "The offset { ";
      for (auto x: _offset) msg += std::to_string(x) + " ";
      msg += "} and size { ";
      for (auto x: _shape) msg += std::to_string(x) + " ";
      msg += "} of an Array must match its underlying SharedSpan shape { ";
      for (auto x: d_shape) msg += std::to_string(x) + " ";
      msg += "}";
      throw std::runtime_error(msg);
    }
  }
public:
  // sub-array access
  Array<T> view() const; // whole array non-owning view
  //! View a single sub-array at index `i` or an offset-smaller array if single=false;
  Array<T> view(ind_t i) const;
  //! View the sub-arrays from `i` to `j-1`
  Array<T> view(ind_t i, ind_t j) const;
  Array<T> view(const shape_t&) const;
  Array<T> extract(ind_t i) const;

  template<typename I>
  std::enable_if_t<std::is_integral_v<I>, Array<T>>
  extract(const Array<I>& i) const;

  template<typename I>
  std::enable_if_t<std::is_integral_v<I>, Array<T>>
  extract(const std::vector<I>& i) const;

  template<typename I, size_t Nel>
  std::enable_if_t<std::is_integral_v<I>, Array<T>>
  extract(const std::array<I,Nel>& i) const;

  template<typename I, size_t Nel>
  std::enable_if_t<std::is_integral_v<I>, Array<T>>
  extract(const std::vector<std::array<I,Nel>>& i) const;

  Array<T> extract(const Array<bool>& i) const;
  Array<T> extract(const std::vector<bool>& i) const;
  bool set(const ind_t i, const Array<T>& in);
  bool set(const ind_t i, const std::vector<T>& in);
  template<size_t Nel> bool set(const ind_t i, const std::array<T,Nel>& in);
  T set(const shape_t& sub, T in);
  Array<T>& append(const ind_t, const Array<T>&);
  std::string to_string() const;
  std::string to_string(const ind_t) const;

  ind_t resize(const shape_t&, T init=T(0));
  template<class I> ind_t resize(const I, T init=T(0));
  bool all(ind_t n=0) const;
  bool any(ind_t n=0) const;
  bool all(T val, ind_t n=0) const;
  bool any(T val, ind_t n=0) const;
  ind_t count(ind_t n=0) const;
  ind_t count(T val, ind_t n=0) const;
  ind_t first(ind_t n=0) const;
  ind_t first(T val, ind_t n=0) const;
  ind_t last(ind_t n=0) const;
  ind_t last(T val, ind_t n=0) const;
  Array<int> round() const;
  Array<int> floor() const;
  Array<int> ceil() const;
  Array<T> sum(ind_t dim=0) const;
  Array<T> prod(ind_t dim=0) const;
  Array<T> min(ind_t dim=0) const;
  Array<T> max(ind_t dim=0) const;
  T sum() const;
  T prod() const;
  template<class R, size_t Nel>
  bool match(ind_t i, ind_t j, const std::array<R,Nel>& rot, int order=1) const;
  bool match(ind_t i, ind_t j, ops op=ops::plus, T val=T{0}) const;
  bool all(cmp expr, T val) const;
  bool any(cmp expr, T val) const;
  ind_t first(cmp expr, T val) const;
  ind_t last(cmp expr, T val) const;
  ind_t count(cmp expr, T val) const;
  Array<bool> is(cmp expr, T val) const;
  Array<ind_t> find(cmp expr, T val) const;
  template<class R> Array<bool> is(cmp expr, const Array<R>& that) const;
  template<class R> std::vector<bool> is(cmp expr, const std::vector<R>& val) const;
  template<class R> bool is(const Array<R>& that) const;
  std::vector<bool> is_unique() const;
  std::vector<ind_t> unique_idx() const;
  Array<T> unique() const;
  Array<T>  operator-() const;
  Array<T>& operator +=(const T&);
  Array<T>& operator -=(const T&);
  Array<T>& operator *=(const T&);
  Array<T>& operator /=(const T&);
  template<class R> Array<T>& operator +=(const Array<R>&);
  template<class R> Array<T>& operator -=(const Array<R>&);
  template<class R> Array<T>& operator *=(const Array<R>&);
  template<class R> Array<T>& operator /=(const Array<R>&);
  T dot(ind_t i, ind_t j) const;
  T norm(ind_t i) const;
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value>>
  void permute(std::vector<I>& p);
  bool swap(ind_t a, ind_t b);
  bool swap(ind_t i, ind_t a, ind_t b);
  std::vector<T> to_std() const;
  T* ptr(const ind_t i0);
  T* ptr(const shape_t& partial_subscript);
  const T* ptr(const ind_t i0) const;
  const T* ptr(const shape_t& partial_subscript) const;
  T& val(const ind_t i0);
  T& val(const shape_t& partial_subscript);
  template<typename I> T& val(std::initializer_list<I> l);
  const T& val(const ind_t i0) const;
  const T& val(const shape_t& partial_subscript) const;
  template<typename I> const T& val(std::initializer_list<I> l) const;
  // ^^^^^^^^^^ IMPLEMENTED  ^^^^^^^^^^^vvvvvvvvv TO IMPLEMENT vvvvvvvvvvvvvvvvv
};

template<class T>
class ArrayIt {
public:
  Array<T> array;
  SubIt<ind_t> subit;
public:
  // constructing with array(a) does not copy the underlying data:
  ArrayIt(const Array<T>& a, const SubIt<ind_t>& s): array(a), subit(s){}
  ArrayIt(const Array<T>& a): array(a), subit(a.shape()){}
  //
  ArrayIt<T> begin() const {
    SubIt<ind_t> s(array.shape()); // initialises to first element, e.g., {0,…,0}
    return ArrayIt(array, s);
  }
  ArrayIt<T> end() const {
    return ArrayIt(array, subit.end());
  }
  ArrayIt<T>& operator++() {
    ++subit;
    return *this;
  }
  const SubIt<ind_t>& iterator() const {return subit;}
  bool operator==(const ArrayIt<T>& other) const {
    // add checking to ensure array and other.array point to the same data?
    return subit == other.iterator();
  }
  bool operator!=(const ArrayIt<T>& other) const {
    return subit != other.iterator();
  }
  const T& operator*()  const {return array[*subit];}
  const T* operator->() const {return &(array[*subit]);}
  T& operator*()  {return array[*subit];}
  T* operator->() {return &(array[*subit]);}
};


#include "array.tpp"
} // end namespace brille
#endif // ARRAY_HPP
