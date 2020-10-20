#ifndef ARRAY2_HPP
#define ARRAY2_HPP

#include <functional>
#include <algorithm>
#include <numeric>
#include <vector>
#include <array>
#include <math.h>
#include <cassert>
#include <iostream>
#include <memory>
#include <string>

#include "subscript.hpp"
#include "utilities.hpp"
#include "comparisons.hpp"
#include "approx.hpp"
#include "types.hpp"

namespace brille {

// struct ref_ptr{
//   using type = char;
// };
// using ref_ptr_t = typename ref_ptr::type; // use char by default since it's small

using ref_ptr_t = char;
template<class T,class P> class Array2It;

/*! \brief A multidimensional shared data array with operator overloads

The `Array2<T,P>`` object holds a multidimensional shared data object plus
information necessary to index into that array with an offset, different shape,
and/or different stride than the underlying data contains.

The `Array2<T,P>` object provides overloads for basic operations, a number of
utility methods, and methods to create offset and/or reduced-range mutable
or immutable views into the shared array.
*/
template<class T, class P=ref_ptr_t>
class Array2{
public:
  using ref_t = std::shared_ptr<P>;
  using shape_t = std::array<ind_t,2>;
  using bItr = BroadcastIt2<ind_t>;
  using sItr = SubIt2<ind_t>;
  using aItr = Array2It<T,P>;
protected:
  T*    _data;     //! A (possilby shared) raw pointer
  ind_t _num;      //! The number of elements stored under the raw pointer
  bool  _own;      //! Whether we need to worry about memory management of the pointer
  ref_t _ref;      //! A shared_ptr for reference counting an owned raw pointer
  bool  _mutable;  //! Whether we are allowed to modify the memory under the pointer
  shape_t _offset; //! An offset different from {0,0} if a view
  shape_t _shape;  //! The shape of this view, prod(_shape) ≤ _num
  shape_t _stride; //! The stride of our view
public:
  // accessors
  T* data() {return _data;}
  T* data(const ind_t& idx){ return _data + this->l2l_d(idx);}
  T* data(const ind_t& i0, const ind_t& i1){ return _data + this->s2l_d(i0,i1);}
  T* data(const shape_t& idx){ return _data + this->s2l_d(idx);}
  const T* data() const {return _data;}
  const T* data(const ind_t& idx) const {return _data + this->l2l_d(idx);}
  const T* data(const ind_t& i0, const ind_t& i1) const { return _data + this->s2l_d(i0,i1);}
  const T* data(const shape_t& idx) const {return _data + this->s2l_d(idx);}
  ind_t raw_size() const {return _num;}
  bool own() const {return _own;}
  ref_t ref() const {return _ref;}
  bool  ismutable(void) const {return _mutable;}
  bool  isconst(void) const {return !_mutable;}
  ind_t ndim(void) const {return 2;}
  ind_t numel(void) const {
    return _shape[0]*_shape[1];
  }
  shape_t offset(void) const {return _offset;}
  ind_t size(const ind_t dim) const {
    assert(dim < 2);
    return _shape[dim];
  }
  shape_t shape(void) const {return _shape;}
  shape_t stride(void) const {return _stride;}
  shape_t cstride() const {
    shape_t cs(_stride);
    for (auto& s: cs) s *= sizeof(T);
    return cs;
  }
  sItr subItr() const { return sItr(_shape); }
  sItr subItr(const shape_t& fix) const { return sItr(_shape, fix); }
  aItr valItr() const { return aItr(*this); }
  bItr broadcastItr(const shape_t& other) const { return bItr(_shape, other); }
  template<class R, class Q>
  bItr broadcastItr(const Array2<R,Q>& other) const {return bItr(_shape, other.shape());}
  bool is_column_ordered() const { // stride {1,2,4,8} is column ordered
    return _stride[0] <= _stride[1];
  }
  bool is_row_ordered() const { // stride {8,4,2,1} is row ordered
    return _stride[1] <= _stride[0];
  }
  bool is_contiguous() const {
    shape_t expected({1, 1});
    if (this->is_row_ordered())    expected[0] = _shape[1];
    if (this->is_column_ordered()) expected[1] = _shape[0];
    // if a dimension is 1 (or 0?) then its stride does not impact here
    return (_shape[0]<2||expected[0]==_stride[0]) && (_shape[1]<2||expected[1]==_stride[1]);
  }
  // empty initializer
  explicit Array2()
  : _data(nullptr), _num(0), _own(false), _ref(std::make_shared<P>()),
  _mutable(false), _offset({0,0}), _shape({0,0}), _stride({0,1})
  {}
  // 1D initializer
  Array2(T* data, const ind_t num, const bool own, const bool mut=true)
  : _data(data), _num(num), _own(own), _ref(std::make_shared<P>()),
    _mutable(mut), _offset({0,0}), _shape({num,1}), _stride({1,1})
  {
    this->init_check();
  }
  // 2D initializer
  Array2(T* data, const ind_t num, const bool own,
    const shape_t& shape, const shape_t& stride, const bool mut=true)
  : _data(data), _num(num), _own(own), _ref(std::make_shared<P>()),
    _mutable(mut), _offset({0,0}), _shape(shape), _stride(stride)
  {
    this->init_check();
  }
  // 2D initializer with reference-counter specified (used in pybind11 wrapper)
  Array2(T* data, const ind_t num, const bool own, ref_t ref,
    const shape_t& shape, const shape_t& stride, const bool mut=true)
  : _data(data), _num(num), _own(own), _ref(ref),
    _mutable(mut), _offset({0,0}), _shape(shape), _stride(stride)
  {
    this->init_check();
  }
  // 2D view constructor
  Array2(T* data, const ind_t num, const bool own, const ref_t& ref,
    const shape_t& offset, const shape_t& shape, const shape_t& stride, const bool mut=true)
  : _data(data), _num(num), _own(own), _ref(ref), _mutable(mut), _offset(offset),
    _shape(shape), _stride(stride)
  {
    this->init_check();
  }
  // 2D allocate new memory constructor
  Array2(const ind_t s0, const ind_t s1)
  : _mutable(true), _offset({0,0}), _shape({s0,s1}), _stride({s1,1})
  {
    this->construct();
    this->init_check();
  }
  // 2D initialize new memory constructor
  Array2(const ind_t s0, const ind_t s1, const T init)
  : _mutable(true), _offset({0,0}), _shape({s0,s1}), _stride({s1,1})
  {
    this->construct(init);
    this->init_check();
  }
  // 2D allocate new memory constructor
  Array2(const shape_t& shape)
  : _mutable(true), _offset({0,0}), _shape(shape)
  {
    this->set_stride();
    this->construct();
    this->init_check();
  }
  // ND initialize new memory constructor
  Array2(const shape_t& shape, const T init)
  : _mutable(true), _offset({0,0}), _shape(shape)
  {
    this->set_stride();
    this->construct(init);
    this->init_check();
  }
  // ND allocate new memory with specified stride constructor
  Array2(const shape_t& shape, const shape_t& stride)
  : _mutable(true), _offset({0,0}), _shape(shape), _stride(stride)
  {
    this->construct();
    this->init_check();
  }
  // ND initialize new memory with specified stride constructor
  Array2(const shape_t& shape, const shape_t& stride, const T init)
  : _mutable(true), _offset({0,0}), _shape(shape), _stride(stride)
  {
    this->construct(init);
    this->init_check();
  }
  // otherwise possibly ambiguous constructor from std::vector
  static Array2<T,P> from_std(const std::vector<T>& data){
    ind_t num = static_cast<ind_t>(data.size());
    T* d = new T[num]();
    for (ind_t i=0; i<num; ++i) d[i] = data[i];
    // take ownership of d while creating the Array
    return Array2<T,P>(d, num, true);
  }
  template<class R, size_t Nel>
  static Array2<T,P> from_std(const std::vector<std::array<R,Nel>>& data){
    shape_t shape{{static_cast<ind_t>(data.size()), static_cast<ind_t>(Nel)}};
    shape_t stride{shape[1], 1};
    ind_t num = shape[0]*shape[1];
    T* d = new T[num]();
    ind_t x{0};
    for (ind_t i=0; i<shape[0]; ++i) for (ind_t j=0; j<shape[1]; ++j)
      d[x++] = static_cast<T>(data[i][j]);
    return Array2<T,P>(d, num, true, shape, stride);
  }

  // Rule-of-five definitions since we may not own the associated raw array:
  Array2(const Array2<T,P>& o)
  : _data(o._data), _num(o._num), _own(o._own), _ref(o._ref),
    _mutable(o._mutable), _offset(o._offset), _shape(o._shape), _stride(o._stride)
  {}
  ~Array2(){
    if (_own && _ref.use_count()==1 && _data != nullptr) delete[] _data;
  }
  Array2<T,P>& operator=(const Array2<T,P>& other){
    if (this != &other){
      if (_own){
        T* old_data = _data;
        ref_t old_ref = _ref;
        _data = other._data;
        _ref = other._ref;
        if (old_ref.use_count()==1 && old_data != nullptr) delete[] old_data;
      } else {
        _data = other._data;
        _ref = other._ref;
      }
      _own = other.own();
      _num = other.raw_size();
      _mutable = other.ismutable();
      _offset = other.offset();
      _shape = other.shape();
      _stride = other.stride();
    }
    return *this;
  }
  // Array2(Array2<T,P>&& other) noexcept
  // : _data(std::exchange(other._data, nullptr))
  // {
  //   std::swap(_num, other._num);
  //   std::swap(_own, other._own);
  //   std::swap(_ref, other._ref);
  //   std::swap(_mutable, other._mutable);
  //   std::swap(_offset, other._offset);
  //   std::swap(_shape, other._shape);
  //   std::swap(_stride, other._stride);
  // }
  // Array2<T,P>& operator=(Array2<T,P>&& other) noexcept {
  //   std::swap(_data, other._data);
  //   std::swap(_num, other._num);
  //   std::swap(_own, other._own);
  //   std::swap(_ref, other._ref);
  //   std::swap(_mutable, other._mutable);
  //   std::swap(_offset, other._offset);
  //   std::swap(_shape, other._shape);
  //   std::swap(_stride, other._stride);
  //   return *this;
  // }


  // type casting requires a copy
  // the reference pointer type does not need to be P since a new raw array is made
  // handles also Array2<T,P>(Array2<T,Q>&) reference type conversion
  template<class R, class RP>
  Array2(const Array2<R,RP>& other)
  : _mutable(true), _offset({0,0}), _shape(other.shape())
  {
    this->set_stride();
    this->construct();
    for (auto x: this->subItr()) _data[s2l_d(x)] = static_cast<T>(other[x]);
  }
  template<class R, class RP>
  Array2<T,P>& operator=(const Array2<R,RP>& other){
    _mutable = true;
    _shape = other.shape();
    _offset[0] = _offset[1] = 0;
    this->set_stride();
    this->construct();
    for (auto x: this->subItr()) _data[s2l_d(x)] = static_cast<T>(other[x]);
    return *this;
  }

  // modifiers
  bool make_mutable() {_mutable = true; return _mutable;}
  bool make_immutable() {_mutable = false; return _mutable;}
  Array2<T,ref_ptr_t> decouple() const {
    if (!_own || !_mutable || _ref.use_count() > 1) return this->_decouple();
    return *this;
  }

  // data accessors
        T& operator[](ind_t    lin)       {return _data[this->l2l_d(lin)];}
  const T& operator[](ind_t    lin) const {return _data[this->l2l_d(lin)];}

        T& operator[](shape_t& sub)       {return _data[this->s2l_d(sub)];}
  const T& operator[](shape_t& sub) const {return _data[this->s2l_d(sub)];}
protected: // so inherited classes can calculate subscript indexes into their data
  ind_t l2l_d(const ind_t l) const {
    return offset_sub2lin(_offset, lin2sub(l, _stride), _stride);
  }
  ind_t ij2l_d(const ind_t x, const ind_t y) const {
    return offset_sub2lin(_offset, x, y, _stride);
  }
  ind_t s2l_d(const shape_t& s) const {
    return offset_sub2lin(_offset, s, _stride);
  }
private:
  ind_t size_from_shape(const shape_t& s) const {
    return s[0]*s[1];
  }
  ind_t size_from_shape() const {return this->size_from_shape(_shape);}
  void construct() {
    _num = this->size_from_shape();
    if (_num > 0){
      _ref = std::make_shared<P>();
      _data = new T[_num]();
      _own = true;
    } else {
      _data = nullptr;
      _own = false;
    }
  }
  void construct(const T init){
    this->construct();
    if (_num > 0 && _data != nullptr) std::fill(_data, _data+_num, init);
  }
  void set_stride(void){
    _stride[1] = 1;
    _stride[0] = _shape[1];
  }
  void init_check(void){
    shape_t sh{{_offset[0]+_shape[0], _offset[1]+_shape[1]}};
    ind_t offset_size = this->size_from_shape(sh);
    if (_num < offset_size) {
      std::string msg = "The offset { ";
      for (auto x: _offset) msg += std::to_string(x) + " ";
      msg += "} and size { ";
      for (auto x: _shape) msg += std::to_string(x) + " ";
      msg += "} of an Array must not exceed the allocated pointer size ";
      msg += std::to_string(_num);
      throw std::runtime_error(msg);
    }
  }
  shape_t calculate_stride(const shape_t& shape) const {
    shape_t stride{{1,1}};
    if (_stride[0] < _stride[1]){
      stride[1] = shape[0];
    } else {
      stride[0] = shape[1];
    }
    return stride;
  }
  void reset_stride(){
    if (!this->is_contiguous())
      throw std::runtime_error("Re-calculating non-contiguous strides is not yet working");
    _stride = this->calculate_stride(_shape);
  }
  Array2<T,ref_ptr_t> _decouple() const {
    ind_t nnum = this->size_from_shape(_shape);
    T* new_data = new T[nnum]();
    shape_t new_st{_stride};
    if (this->is_contiguous() && nnum == _num) {
      std::copy(_data, _data+_num, new_data);
    } else {
      //subscript conversion necessary due to offset or strided array
      new_st = this->calculate_stride(_shape);
      //                                  vv(no offset)vvv         vvv(offset)vvv
      for (auto x: this->subItr()) new_data[sub2lin(x,new_st)] = _data[this->s2l_d(x)];
    }
    bool new_own = true; // always take ownership of C++ allocated memory
    auto new_ref = std::make_shared<ref_ptr_t>(); // always use the default with C++ created arrays
    bool new_mut = true; // new allocated memory should be mutable
    return Array2<T,ref_ptr_t>(new_data, nnum, new_own, new_ref, _shape, new_st, new_mut);
  }
public:
  // sub-array access
  Array2<T,P> view() const; // whole array non-owning view
  //! View a single sub-array at index `i` or an offset-smaller array if single=false;
  Array2<T,P> view(ind_t i) const;
  //! View the sub-arrays from `i` to `j-1`
  Array2<T,P> view(ind_t i, ind_t j) const;
  Array2<T,P> view(const shape_t&) const;
  // duplication of one or more sub-arrays:
  Array2<T,ref_ptr_t>
  extract(ind_t i) const;
  template<class I, class IP>
  std::enable_if_t<std::is_integral_v<I>, Array2<T,ref_ptr_t>>
  extract(const Array2<I,IP>& i) const;
  template<typename I>
  std::enable_if_t<std::is_integral_v<I>, Array2<T,ref_ptr_t>>
  extract(const std::vector<I>& i) const;
  template<typename I, size_t Nel>
  std::enable_if_t<std::is_integral_v<I>, Array2<T,ref_ptr_t>>
  extract(const std::array<I,Nel>& i) const;
  template<typename I, size_t Nel>
  std::enable_if_t<std::is_integral_v<I>, Array2<T,ref_ptr_t>>
  extract(const std::vector<std::array<I,Nel>>& i) const;
  template<class RP>
  Array2<T,ref_ptr_t>
  extract(const Array2<bool,RP>& i) const;
  Array2<T,ref_ptr_t>
  extract(const std::vector<bool>& i) const;
  template<class RP>
  bool set(const ind_t i, const Array2<T,RP>& in);
  template<class R, class RP>
  bool set(const ind_t i, const Array2<R,RP>& in);
  bool set(const ind_t i, const std::vector<T>& in);
  template<size_t Nel> bool set(const ind_t i, const std::array<T,Nel>& in);
  T set(const shape_t& sub, T in);
  template<class RP>
  Array2<T,P>& append(const ind_t, const Array2<T,RP>&);
  std::string to_string() const;
  std::string to_string(const ind_t) const;

  Array2<T,P>& reshape(const shape_t& ns);
  Array2<T,P>& resize(const shape_t&, T init=T(0));
  template<class I> Array2<T,P>& resize(const I, T init=T(0));
  bool all(ind_t n=0) const;
  bool any(ind_t n=0) const;
  ind_t count(ind_t n=0) const;
  ind_t first(ind_t n=0) const;
  ind_t last(ind_t n=0) const;
  // bool all() const;
  // bool any() const;
  // ind_t count() const;
  // ind_t first() const;
  // ind_t last() const;
  bool all(T val, ind_t n=0) const;
  bool any(T val, ind_t n=0) const;
  ind_t count(T val, ind_t n=0) const;
  ind_t first(T val, ind_t n=0) const;
  ind_t last(T val, ind_t n=0) const;
  Array2<int,ref_ptr_t> round() const;
  Array2<int,ref_ptr_t> floor() const;
  Array2<int,ref_ptr_t> ceil() const;
  Array2<T,ref_ptr_t> sum(ind_t dim=0) const;
  Array2<T,ref_ptr_t> prod(ind_t dim=0) const;
  Array2<T,ref_ptr_t> min(ind_t dim=0) const;
  Array2<T,ref_ptr_t> max(ind_t dim=0) const;
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
  Array2<bool,ref_ptr_t> is(cmp expr, T val) const;
  std::vector<ind_t> find(cmp expr, T val) const;
  template<class R, class RP>
  Array2<bool,ref_ptr_t>
  is(cmp expr, const Array2<R,RP>& that) const;
  template<class R>
  std::vector<bool>
  is(cmp expr, const std::vector<R>& val) const;
  template<class R, class RP>
  bool
  is(const Array2<R,RP>& that) const;
  std::vector<bool> is_unique() const;
  std::vector<ind_t> unique_idx() const;
  Array2<T,ref_ptr_t> unique() const;
  Array2<T,ref_ptr_t>  operator-() const;
  Array2<T,P>& operator +=(const T&);
  Array2<T,P>& operator -=(const T&);
  Array2<T,P>& operator *=(const T&);
  Array2<T,P>& operator /=(const T&);
  template<class R, class RP> Array2<T,P>& operator +=(const Array2<R,RP>&);
  template<class R, class RP> Array2<T,P>& operator -=(const Array2<R,RP>&);
  template<class R, class RP> Array2<T,P>& operator *=(const Array2<R,RP>&);
  template<class R, class RP> Array2<T,P>& operator /=(const Array2<R,RP>&);
  T dot(ind_t i, ind_t j) const;
  T norm(ind_t i) const;
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value>>
  void permute(std::vector<I>& p);
  bool swap(ind_t a, ind_t b);
  bool swap(ind_t i, ind_t a, ind_t b);
  std::vector<T> to_std() const;
  T* ptr(const ind_t i0);
  T* ptr(const ind_t i0, const ind_t j0);
  T* ptr(const shape_t& partial_subscript);
  const T* ptr(const ind_t i0) const;
  const T* ptr(const ind_t i0, const ind_t j0) const;
  const T* ptr(const shape_t& partial_subscript) const;
  T& val(const ind_t i0);
  T& val(const ind_t i0, const ind_t j0);
  T& val(const shape_t& partial_subscript);
  template<typename I> T& val(std::initializer_list<I> l);
  const T& val(const ind_t i0) const;
  const T& val(const ind_t i0, const ind_t j0) const;
  const T& val(const shape_t& partial_subscript) const;
  template<typename I> const T& val(std::initializer_list<I> l) const;

  Array2<T,P> contiguous_copy() const;
  // ^^^^^^^^^^ IMPLEMENTED  ^^^^^^^^^^^vvvvvvvvv TO IMPLEMENT vvvvvvvvvvvvvvvvv

};

template<class T,class P>
class Array2It {
public:
  Array2<T,P> array;
  SubIt2<ind_t> subit;
public:
  // constructing with array(a) does not copy the underlying data:
  explicit Array2It()
  : array(), subit()
  {}
  Array2It(const Array2<T,P>& a, const SubIt2<ind_t>& s)
  : array(a), subit(s)
  {}
  Array2It(const Array2<T,P>& a)
  : array(a), subit(a.shape()) // initialises to first element, e.g., {0,…,0}
  {}
  //
  Array2It<T,P> begin() const {
    return Array2It(array);
  }
  Array2It<T,P> end() const {
    return Array2It(array, subit.end());
  }
  Array2It<T,P>& operator++() {
    ++subit;
    return *this;
  }
  const SubIt2<ind_t>& iterator() const {return subit;}
  template<class Q>
  bool operator==(const Array2It<T,Q>& other) const {
    // add checking to ensure array and other.array point to the same data?
    return subit == other.iterator();
  }
  template<class Q>
  bool operator!=(const Array2It<T,Q>& other) const {
    return subit != other.iterator();
  }
  const T& operator*()  const {return array[*subit];}
  const T* operator->() const {return &(array[*subit]);}
  T& operator*()  {return array[*subit];}
  T* operator->() {return &(array[*subit]);}
};


#include "array2.tpp"
} // end namespace brille
#endif // ARRAY_HPP
