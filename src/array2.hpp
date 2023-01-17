/* This file is part of brille

Copyright © 2020-2022 Greg Tucker <gregory.tucker@ess.eu>

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
    \brief Two dimensional shared memory Array2 class definition.
*/
#ifndef BRILLE_ARRAY2_HPP
#define BRILLE_ARRAY2_HPP
#include <memory>
#include <atomic>
#include "subscript.hpp"
#include "utilities.hpp"
#include "comparisons.hpp"
#include "types.hpp"
#include "array_pre.hpp"
#include "hdf_interface.hpp"
namespace brille {
/*! \brief A multidimensional shared data array with operator overloads

The `Array2<T>`` object holds a multidimensional shared data object plus
information necessary to index into that array with an offset, different shape,
and/or different stride than the underlying data contains.

The `Array2<T>` object provides overloads for basic operations, a number of
utility methods, and methods to create offset and/or reduced-range mutable
or immutable views into the shared array.
*/
template<class T> class Array2{
public:
  using ref_t = std::shared_ptr<void>;
  using shape_t = std::array<ind_t,2>;
  using bItr = BroadcastIt2<ind_t>;
  using sItr = SubIt2<ind_t>;
  using aItr = Array2It<T>;
protected:
  T*    _data;     //! A (possilby shared) raw pointer
  ind_t _num;      //! The number of elements stored under the raw pointer
  ind_t _shift;    //! A linear shift to where our subscript indexing begins (_data+shift)[subscript]
  bool  _own;      //! Whether we need to worry about memory management of the pointer
  ref_t _ref;      //! A shared_ptr for reference counting an owned raw pointer
  bool  _mutable;  //! Whether we are allowed to modify the memory under the pointer
  shape_t _shape;  //! The shape of this view, prod(_shape) ≤ _num
  shape_t _stride; //! The stride of our view
public:
  bool operator!=(const Array2<T>& other) const {
      // Two arrays are not equal if their *stored values* are unequal.
      // Their raw pointers, total number of elements, shifts, references, mutabilities, and ownerships don't matter
      if (_shape != other._shape) return true;
      for (const auto & sub: subItr()) if (this->operator[](sub) != other[sub]) return true;
      return false;
  }
  // accessors
  T* data() {return _data;}
  T* data(const ind_t& idx){ return _data + this->l2l_d(idx);}
  T* data(const ind_t& i0, const ind_t& i1){ return _data + this->s2l_d(i0,i1);}
  T* data(const shape_t& idx){ return _data + this->s2l_d(idx);}
  const T* data() const {return _data;}
  const T* data(const ind_t& idx) const {return _data + this->l2l_d(idx);}
  const T* data(const ind_t& i0, const ind_t& i1) const { return _data + this->s2l_d(i0,i1);}
  const T* data(const shape_t& idx) const {return _data + this->s2l_d(idx);}
  [[nodiscard]] ind_t raw_size() const {return _num;}
  [[nodiscard]] ind_t raw_shift() const {return _shift;}
  [[nodiscard]] bool own() const {return _own;}
  [[nodiscard]] ref_t ref() const {return _ref;}
  [[nodiscard]] bool  ismutable() const {return _mutable;}
  [[nodiscard]] bool  isconst() const {return !_mutable;}
  [[nodiscard]] ind_t ndim() const {return 2;}
  [[nodiscard]] ind_t numel() const {
    return _shape[0]*_shape[1];
  }
  [[nodiscard]] ind_t size(const ind_t dim) const {
    assert(dim < 2);
    return _shape[dim];
  }
  [[nodiscard]] shape_t shape() const {return _shape;}
  [[nodiscard]] shape_t stride() const {return _stride;}
  [[nodiscard]] shape_t cstride() const {
    shape_t cs(_stride);
    for (auto& s: cs) s *= sizeof(T);
    return cs;
  }
  [[nodiscard]] sItr subItr() const { return sItr(_shape); }
  [[nodiscard]] sItr subItr(const shape_t& fix) const { return {_shape, fix}; }
  aItr valItr() const { return aItr(*this); }
  [[nodiscard]] bItr broadcastItr(const shape_t& other) const { return {_shape, other}; }
  template<class R>
  bItr broadcastItr(const Array2<R>& other) const {return bItr(_shape, other.shape());}
  [[nodiscard]] bool is_column_ordered() const { // stride {1,2,4,8} is column ordered
    return _stride[0] <= _stride[1];
  }
  [[nodiscard]] bool is_row_ordered() const { // stride {8,4,2,1} is row ordered
    return _stride[1] <= _stride[0];
  }
  [[nodiscard]] bool is_contiguous() const {
    shape_t expected({1, 1});
    if (this->is_row_ordered())    expected[0] = _shape[1];
    if (this->is_column_ordered()) expected[1] = _shape[0];
    // if a dimension is 1 (or 0?) then its stride does not impact here
    return (_shape[0]<2||expected[0]==_stride[0]) && (_shape[1]<2||expected[1]==_stride[1]);
  }
  [[nodiscard]] bool is_vector() const {
    return (_shape[0] == 1u || _shape[1] == 1u);
  }
  // empty initializer
  explicit Array2()
  : _data(nullptr), _num(0), _shift(0u), _own(false), _ref(std::make_shared<std::atomic<char>>()),
  _mutable(false), _shape({0,0}), _stride({0,1})
  {}
  // 1D initializer
  Array2(T* data, const ind_t num, const bool own, const bool mut=true)
  : _data(data), _num(num), _shift(0u), _own(own), _ref(std::make_shared<std::atomic<char>>()),
    _mutable(mut), _shape({num,1}), _stride({1,1})
  {
    this->init_check();
  }
  // 2D initializer
  Array2(T* data, const ind_t num, const bool own,
    const shape_t& shape, const shape_t& stride, const bool mut=true)
  : _data(data), _num(num), _shift(0u), _own(own), _ref(std::make_shared<std::atomic<char>>()),
    _mutable(mut), _shape(shape), _stride(stride)
  {
    this->init_check();
  }
  // 2D initializer with reference-counter specified (used in pybind11 wrapper)
  template<class P>
  Array2(T* data, const ind_t num, const bool own, std::shared_ptr<P> ref,
    const shape_t& shape, const shape_t& stride, const bool mut=true)
  : _data(data), _num(num), _shift(0u), _own(own), _ref(ref),
    _mutable(mut), _shape(shape), _stride(stride)
  {
    this->init_check();
  }
  // 2D view constructor
  template<class P>
  Array2(T* data, const ind_t num, const ind_t shift, const bool own, const std::shared_ptr<P>& ref,
    const shape_t& shape, const shape_t& stride, const bool mut=true)
  : _data(data), _num(num), _shift(shift), _own(own), _ref(ref), _mutable(mut),
    _shape(shape), _stride(stride)
  {
    this->init_check();
  }
  // 2D allocate new memory constructor
  Array2(const ind_t s0, const ind_t s1)
  : _shift(0u), _mutable(true), _shape({s0,s1}), _stride({s1,1})
  {
    this->construct();
    this->init_check();
  }
  // 2D initialize new memory constructor
  Array2(const ind_t s0, const ind_t s1, const T init)
  : _shift(0u), _mutable(true), _shape({s0,s1}), _stride({s1,1})
  {
    this->construct(init);
    this->init_check();
  }
  // 2D allocate new memory constructor
  Array2(const shape_t& shape)
  : _shift(0u), _mutable(true), _shape(shape)
  {
    this->set_stride();
    this->construct();
    this->init_check();
  }
  // ND initialize new memory constructor
  Array2(const shape_t& shape, const T init)
  : _shift(0u), _mutable(true), _shape(shape)
  {
    this->set_stride();
    this->construct(init);
    this->init_check();
  }
  // ND allocate new memory with specified stride constructor
  Array2(const shape_t& shape, const shape_t& stride)
  : _shift(0u), _mutable(true), _shape(shape), _stride(stride)
  {
    this->construct();
    this->init_check();
  }
  // ND initialize new memory with specified stride constructor
  Array2(const shape_t& shape, const shape_t& stride, const T init)
  : _shift(0u), _mutable(true), _shape(shape), _stride(stride)
  {
    this->construct(init);
    this->init_check();
  }
  // otherwise possibly ambiguous constructor from std::vector
  static Array2<T> from_std(const std::vector<T>& data){
    ind_t num = static_cast<ind_t>(data.size());
    T* d = new T[num]();
    for (ind_t i=0; i<num; ++i) d[i] = data[i];
    // take ownership of d while creating the Array
    return Array2<T>(d, num, true);
  }
  template<class R, size_t Nel>
  static Array2<T> from_std(const std::vector<std::array<R,Nel>>& data){
    shape_t shape{{static_cast<ind_t>(data.size()), static_cast<ind_t>(Nel)}};
    shape_t stride{shape[1], 1};
    ind_t num = shape[0]*shape[1];
    T* d = new T[num]();
    ind_t x{0};
    for (ind_t i=0; i<shape[0]; ++i) for (ind_t j=0; j<shape[1]; ++j)
      d[x++] = static_cast<T>(data[i][j]);
    return Array2<T>(d, num, true, shape, stride);
  }

  // Construct an Array2 from an Array, unravelling higher dimensions
  Array2(const Array<T>& nd)
  : _num(0u), _shift(0u), _shape({0,1}), _stride({1,1})
  {
    // we always want to construct contiguous row-ordered arrays
    if (nd.ndim()>0){
      _shape[0] = nd.size(0);
      if (nd.ndim()>1){
        for (ind_t i=1; i<nd.ndim(); ++i) _shape[1] *= nd.size(i);
        _stride[0] = _shape[1];
      }
      // the following is a very-light (metadata) copy unless the underlying
      // data is not row-ordered contiguous.
      auto cnd = nd.contiguous_row_ordered_copy();
      _data = cnd.data();
      _num = cnd.raw_size();
      _shift = cnd.raw_shift();
      _own = cnd.own();
      _ref = cnd.ref();
      _mutable = cnd.ismutable();
    }
  }

  // Rule-of-five definitions since we may not own the associated raw array:
  Array2(const Array2<T>& o)
  : _data(o._data), _num(o._num), _shift(o._shift), _own(o._own), _ref(o._ref),
    _mutable(o._mutable), _shape(o._shape), _stride(o._stride)
  {}
  ~Array2(){
    if (_own && _ref.use_count()==1 && _data != nullptr) delete[] _data;
  }
  Array2<T>& operator=(const Array2<T>& other){
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
      _shift = other.raw_shift();
      _mutable = other.ismutable();
      _shape = other.shape();
      _stride = other.stride();
    }
    return *this;
  }

  // type casting requires a copy
  // the reference pointer type does not need to be P since a new raw array is made
  // handles also Array2<T>(Array2<T>&) reference type conversion
  template<class R>
  explicit Array2(const Array2<R>& other)
  : _shift(0u), _mutable(true), _shape(other.shape())
  {
    this->set_stride();
    this->construct();
    for (const auto & x: this->subItr()) _data[s2l_d(x)] = static_cast<T>(other[x]);
  }
  template<class R>
  Array2<T>& operator=(const Array2<R>& other){
    _mutable = true;
    _shape = other.shape();
    _shift = 0;
    this->set_stride();
    this->construct();
    for (const auto & x: this->subItr()) _data[s2l_d(x)] = static_cast<T>(other[x]);
    return *this;
  }

  // modifiers
  bool make_mutable() {_mutable = true; return _mutable;}
  bool make_immutable() {_mutable = false; return _mutable;}
  Array2<T> decouple() const {
    if (!_own || !_mutable || _ref.use_count() > 1) return this->_decouple();
    return *this;
  }

  // data accessors
        T& operator[](ind_t    lin)       {return _data[this->l2l_d(lin)];}
  const T& operator[](ind_t    lin) const {return _data[this->l2l_d(lin)];}

        T& operator[](const shape_t& sub)       {return _data[this->s2l_d(sub)];}
  const T& operator[](const shape_t& sub) const {return _data[this->s2l_d(sub)];}
protected: // so inherited classes can calculate subscript indexes into their data
  [[nodiscard]] ind_t l2l_d(const ind_t l) const {
    return l + _shift;
  }
  [[nodiscard]] ind_t ij2l_d(const ind_t x, const ind_t y) const {
    return sub2lin(x, y, _stride) + _shift;
  }
  [[nodiscard]] ind_t s2l_d(const shape_t& s) const {
    assert(brille::utils::subscript_ok(s, _shape));
    return sub2lin(s, _stride) + _shift;
  }
private:
  [[nodiscard]] ind_t size_from_shape(const shape_t& s) const {
    return s[0]*s[1];
  }
  [[nodiscard]] ind_t size_from_shape() const {return this->size_from_shape(_shape);}
  void construct() {
    _num = this->size_from_shape();
    if (_num > 0){
      _ref = std::make_shared<std::atomic<char>>();
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
  void set_stride(){
    _stride[1] = 1;
    _stride[0] = _shape[1];
  }
  void init_check(){
    ind_t offset_size = _shift + this->size_from_shape(_shape);
    if (_num < offset_size) {
      std::string msg = "The shift { " + std::to_string(_shift) + " ";
      msg += "} and size { ";
      for (const auto & x: _shape) msg += std::to_string(x) + " ";
      msg += "} of an Array2 must not exceed the allocated pointer size ";
      msg += std::to_string(_num);
      throw std::runtime_error(msg);
    }
  }
  [[nodiscard]] shape_t calculate_stride(const shape_t& shape) const {
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
  Array2<T> _decouple() const {
    ind_t nnum = this->size_from_shape(_shape);
    T* new_data = new T[nnum]();
    shape_t new_st{_stride};
    if (this->is_contiguous() && nnum == _num) {
      std::copy(_data, _data+_num, new_data);
    } else {
      //subscript conversion necessary due to offset or strided array
      new_st = this->calculate_stride(_shape);
      //                                  vv(no offset)vvv         vvv(offset)vvv
      for (const auto & x: this->subItr()) new_data[sub2lin(x,new_st)] = _data[this->s2l_d(x)];
    }
    bool new_own = true; // always take ownership of C++ allocated memory
    auto new_ref = std::make_shared<std::atomic<char>>(); // always use the default with C++ created arrays
    bool new_mut = true; // new allocated memory should be mutable
    return Array2<T>(new_data, nnum, new_own, new_ref, _shape, new_st, new_mut);
  }
public:
  // sub-array access
  Array2<T> view() const; // whole array non-owning view
  //! View a single sub-array at index `i` or an offset-smaller array if single=false;
  Array2<T> view(ind_t i) const;
  //! View the sub-arrays from `i` to `j-1`
  Array2<T> view(ind_t i, ind_t j) const;
  Array2<T> view(const shape_t&) const;
  // duplication of one or more sub-arrays:
  Array2<T>
  extract(ind_t i) const;
  template<class I> std::enable_if_t<std::is_integral_v<I>, Array2<T>> extract(const Array2<I>& i) const;
  template<class I> std::enable_if_t<std::is_integral_v<I>, Array2<T>> extract(const std::vector<I>& i) const;
  template<class I, size_t Nel> std::enable_if_t<std::is_integral_v<I>, Array2<T>> extract(const std::array<I,Nel>& i) const;
  template<class I, size_t Nel> std::enable_if_t<std::is_integral_v<I>, Array2<T>> extract(const std::vector<std::array<I,Nel>>& i) const;
  Array2<T> extract(const Array2<bool>& i) const;
  Array2<T> extract(const std::vector<bool>& i) const;
  bool set(ind_t i, const Array2<T>& in);
  template<class R>
  bool set(ind_t i, const Array2<R>& in);
  bool set(ind_t i, const std::vector<T>& in);
  template<size_t Nel> bool set(ind_t i, const std::array<T,Nel>& in);
  T set(const shape_t& sub, T in);
  template<class R> std::enable_if_t<std::is_convertible_v<R,T>, Array2<T>&> append(ind_t, const Array2<R>&);
  [[nodiscard]] std::string to_string() const;
  [[nodiscard]] std::string to_string(ind_t) const;

  Array2<T>& reshape(const shape_t& ns);
  Array2<T>& resize(const shape_t&, T init=T(0));
  template<class I> Array2<T>& resize(I, T init=T(0));
  [[nodiscard]] bool all(ind_t n=0) const;
  [[nodiscard]] bool any(ind_t n=0) const;
  [[nodiscard]] ind_t count(ind_t n=0) const;
  [[nodiscard]] ind_t first(ind_t n=0) const;
  [[nodiscard]] ind_t last(ind_t n=0) const;
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
  [[nodiscard]] Array2<int> round() const;
  [[nodiscard]] Array2<int> floor() const;
  [[nodiscard]] Array2<int> ceil() const;
  Array2<T> sum(ind_t dim) const;
  Array2<T> prod(ind_t dim) const;
  Array2<T> min(ind_t dim=0) const;
  Array2<T> max(ind_t dim=0) const;
  T sum() const;
  T prod() const;
  template<class R, size_t Nel>
  bool match(ind_t i, ind_t j, const std::array<R,Nel>& rot, int order=1, T Ttol=T(0), int tol=1) const;
  bool match(ind_t i, ind_t j, ops op=ops::plus, T val=T{0}, T Ttol=T(0), int tol=1) const;
  bool all(cmp expr, T val, T Ttol=T(0), int tol=1) const;
  bool any(cmp expr, T val, T Ttol=T(0), int tol=1) const;
  ind_t first(cmp expr, T val, T Ttol=T(0), int tol=1) const;
  ind_t last(cmp expr, T val, T Ttol=T(0), int tol=1) const;
  ind_t count(cmp expr, T val, T Ttol=T(0), int tol=1) const;
  template<class R> ind_t first(cmp expr, const Array2<R>& val, T Ttol=T(0), R Rtol=R(0), int tol=1) const;
  template<class R> ind_t last(cmp expr, const Array2<R>& val, T Ttol=T(0), R Rtol=R(0), int tol=1) const;
  template<class R> ind_t count(cmp expr, const Array2<R>& val, T Ttol=T(0), R Rtol=R(0), int tol=1) const;
  template<class R> std::vector<ind_t> find(cmp expr, const Array2<R>& val, T Ttol=T(0), R Rtol=R(0), int tol=1) const;
  Array2<bool> is(cmp expr, T val, T Ttol=T(0), int tol=1) const;
  std::vector<ind_t> find(cmp expr, T val, T Ttol=T(0), int tol=1) const;
  /*!\brief Determine the per-element truth-value of this Array2 compared with provided data

  \param expr The comparison to be performed
  \param that The values to compare against

  \return The truth value of the provided comparison between the stored and
          provided data, will be the same shape a this Array2.

  \note The shape of the compared-against Array2 must match or be broadcastab
   le
        to the shape of this Array2.
  */
  template<class R> Array2<bool> is(cmp expr, const Array2<R>& that, T Ttol=T(0), R Rtol=R(0), int tol=1) const;
  /*!\brief Determine the per-row truth-value of this Array2 compared with provided data

  \param expr The comparison to be performed
  \param row The values to compare against, with length equal to the size of the
             second dimension of this Array2.

  \return The per-row truth value of the provided comparison between the stored
          and provided data.
  */
  template<class R> std::vector<bool> row_is(cmp expr, const std::vector<R>& row, T Ttol=T(0), R Rtol=R(0), int tol=1) const;
  template<class R> Array2<bool> row_is(cmp expr, const Array2<R> & rows, T Ttol=T(0), R Rtol=R(0), int tol=1) const;
  /*!\brief Determine the element truth-value of this Array2 compared with provided data

  \param expr The comparison to be performed
  \param vals The values to compare against, with length equal to the number of
              elements in this Array2

  \return The per-element truth value of the provided comparison between the
          stored and provided data.
  */
  template<class R> std::vector<bool> each_is(cmp expr, const std::vector<R>& vals, T Ttol=T(0), R Rtol=R(0), int tol=1) const;
  template<class R> bool is(const Array2<R>& that, T Ttol=T(0), R Rtol=R(0), int tol=1) const;
  template<class R> bool operator==(const Array2<R>& that) const {return this->is(that);}
  [[nodiscard]] std::vector<bool> is_unique(T Ttol=T(0), int tol=1) const;
  [[nodiscard]] std::vector<ind_t> unique_idx(T Ttol=T(0), int tol=1) const;
  Array2<T> unique(T Ttol=T(0), int tol=1) const;
  Array2<T>  operator-() const;
  Array2<T>& operator +=(const T&);
  Array2<T>& operator -=(const T&);
  Array2<T>& operator *=(const T&);
  Array2<T>& operator /=(const T&);
  Array2<T> operator ^(const T&);
  template<class R> Array2<T>& operator +=(const Array2<R>&);
  template<class R> Array2<T>& operator -=(const Array2<R>&);
  template<class R> Array2<T>& operator *=(const Array2<R>&);
  template<class R> Array2<T>& operator /=(const Array2<R>&);
  T dot(ind_t i, ind_t j) const;
  T norm(ind_t i) const;
  template<class I, typename=std::enable_if_t<std::is_integral<I>::value>>
  void permute(std::vector<I>& p);
  bool swap(ind_t a, ind_t b);
  bool swap(ind_t i, ind_t a, ind_t b);
  std::vector<T> to_std() const;
  T* ptr(ind_t i0);
  T* ptr(ind_t i0, ind_t j0);
  T* ptr(const shape_t& partial_subscript);
  const T* ptr(ind_t i0) const;
  const T* ptr(ind_t i0, ind_t j0) const;
  const T* ptr(const shape_t& partial_subscript) const;
  T& val(ind_t i0);
  T& val(ind_t i0, ind_t j0);
  T& val(const shape_t& partial_subscript);
  template<class I> T& val(std::initializer_list<I> l);
  const T& val(ind_t i0) const;
  const T& val(ind_t i0, ind_t j0) const;
  const T& val(const shape_t& partial_subscript) const;
  template<class I> const T& val(std::initializer_list<I> l) const;

  Array2<T> contiguous_copy() const;
  Array2<T> contiguous_row_ordered_copy() const;

  Array2<T> abs() const;
  template<class R> bool is_permutation(const Array2<R>& other, T Ttol=T(0), R Rtol=R(0), int tol=1) const;
  template<class R> std::vector<ind_t> permutation_vector(const Array2<R>& other, T Ttol=T(0), R Rtol=R(0), int tol=1) const;
  //
#if USE_HIGHFIVE
  // Read in from a file
  template<class R>
  static std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, Array2<T>>
  from_hdf(R& obj, const std::string& node){
    using namespace HighFive;
    auto type = obj.getObjectType(node);
    if (type == HighFive::ObjectType::Dataset) {
      auto dataset = obj.getDataSet(node);
      std::vector<size_t> data_shape = dataset.getDimensions();
      if (data_shape.size() != 2u)
        throw std::runtime_error("Only 2D arrays can be Array2 objects");
      auto num = dataset.getElementCount();
      T *data = new T[num]();
      dataset.read(data);
      shape_t shape{{static_cast<ind_t>(data_shape[0]),
                     static_cast<ind_t>(data_shape[1])}};
      shape_t stride{{shape[1], 1}};
      return Array2<T>(data, static_cast<ind_t>(num), true, shape, stride);
    }
    if (type == HighFive::ObjectType::Group){
      ind_t no{0};
      shape_t sh{0,0};
      auto datagroup = obj.getGroup(node);
      if (!datagroup.hasAttribute("size")){
        throw std::runtime_error("How and or why is the size attribute missing?");
      }
      datagroup.getAttribute("size").read(no);
      if (!datagroup.hasAttribute("shape")){
        throw std::runtime_error("How and or why is the shape attribute missing?");
      }
      datagroup.getAttribute("shape").read(sh);
      if (no > 0u) throw std::runtime_error("size 0 Array2 should be DataSets!");
      return Array2<T>(sh);
    }
    throw std::runtime_error("The HDF file must have a Dataset or Group at the specified Array2 node!");
  }
  static Array2<T> from_hdf(const std::string& filename, const std::string& datasetname){
    using namespace HighFive;
    File file(filename, File::ReadOnly);
    return Array2<T>::from_hdf(file, datasetname);
  }
  //
  template<class R>
  std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, bool>
  add_hdf_attributes(R& obj) const;
  template<class R>
  std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, bool>
  to_hdf(R& obj, const std::string& entry) const;
  bool to_hdf(const std::string& filename, const std::string& datasetname, const unsigned permissions=HighFive::File::OpenOrCreate) const;
#endif //USE_HIGHFIVE

  friend std::ostream & operator<<(std::ostream & os, const Array2<T>& a){
    os << a.to_string();
    return os;
  }
};

/*!\brief An iterator for the Array2 class

A simple iterator class to handle the abbreviated for loop syntax to iterate
over all elements of an Array object.
*/
template<class T>
class Array2It {
public:
  Array2<T> array;     //!< The Array2 to be iterated
  SubIt2<ind_t> subit; //!< The subscript iterator of the Array
public:
  //! Explicit empty constructor
  explicit Array2It()
  : array(), subit()
  {}
  /*! \brief Constructor taking all parameters

  \param a The Array2 to be iterated
  \param s The subscript iterator which performs the iteration
  */
  Array2It(const Array2<T>& a, const SubIt2<ind_t>& s)
  : array(a), subit(s)
  {}
  /*! \brief Construct from the Array2 only

  The subscript iterator is constructed from the shape of the passed Array2.

  \param a The Array2 to be iterated
  */
  explicit Array2It(const Array2<T>& a)
  : array(a), subit(a.shape()) // initialises to first element, e.g., {0,…,0}
  {}
  //! Return this ArrayIt with the iterator starting at index {0,0}
  Array2It<T> begin() const {
    return Array2It(array);
  }
  //! Return this ArrayIt with the iterator at its end {array.size(0),0}
  Array2It<T> end() const {
    return Array2It(array, subit.end());
  }
  //! Increment the iterator
  Array2It<T>& operator++() {
    ++subit;
    return *this;
  }
  //! Return the subscript iterator
  [[nodiscard]] const SubIt2<ind_t>& iterator() const {return subit;}
  /*!\brief Check for subscript iterator equivalency

  \note There is no check for held Array2 equivalency
  */
  bool operator==(const Array2It<T>& other) const {
    // add checking to ensure array and other.array point to the same data?
    return subit == other.iterator();
  }
  //! Check for subscript iterator inequivalency
  bool operator!=(const Array2It<T>& other) const {
    return subit != other.iterator();
  }
  //! Return a constant reference to the Array2 value at the current subscript index
  const T& operator*()  const {return array[*subit];}
  //! Return a constant pointer to the Array2 element at the current subscript index
  const T* operator->() const {return &(array[*subit]);}
  //! Return a reference to the Array2 value at the current subscript index
  T& operator*()  {return array[*subit];}
  //! Return a pointer to the Array2 element at the current subscript index
  T* operator->() {return &(array[*subit]);}
};


#include "array2.tpp"
} // end namespace brille
#endif // ARRAY_HPP
