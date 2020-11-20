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

template<class T>
Array2<T> Array2<T>::view() const{
  return Array2<T>(_data,_num,_shift,_own,_ref,_shape,_stride,false);
}
template<class T>
Array2<T> Array2<T>::view(const ind_t i) const {
  if (i<_shape[0]){
    shape_t osize{_shape};
    ind_t oshft{_shift + i*_stride[0]};
    osize[0] = 1;
    return Array2<T>(_data, _num, oshft, _own, _ref, osize, _stride, false);
  }
  throw std::runtime_error("Array2 index too large");
}

template<class T>
Array2<T> Array2<T>::view(const ind_t i, const ind_t j) const {
  if (i<j && i<_shape[0] && j<=_shape[0]){
    shape_t osize{_shape};
    ind_t oshft{_shift + i*_stride[0]};
    osize[0] = j-i;
    return Array2<T>(_data, _num, oshft, _own, _ref, osize, _stride, false);
  }
  throw std::runtime_error("Array2 view indexing error");
}

template<class T>
Array2<T> Array2<T>::view(const shape_t& v) const {
  shape_t osize{_shape};
  ind_t oshft{_shift};
  for (size_t i=0; i<_shape.size(); ++i) if (v[i] < _shape[i]) {
    osize[i] -= v[i];
    oshft += v[i]*_stride[i];
  }
  return Array2<T>(_data, _num, oshft, _own, _ref, osize, _stride, false);
}

template<class T>
Array2<T>
Array2<T>::extract(const ind_t i) const{
  return this->view(i).decouple();
}

template<class T>
template<class I>
std::enable_if_t<std::is_integral_v<I>, Array2<T>>
Array2<T>::extract(const Array2<I>& i) const {
  if (i.numel() != i.size(0))
    throw std::runtime_error("Array2 extraction requires (N,{1,...,1}) shape");
  for (auto x: i.valItr()) if (!(0 <= x && static_cast<ind_t>(x) < _shape[0]))
    throw std::runtime_error("Array2 extract index must be in range");
  shape_t osize{_shape};
  osize[0] = i.size(0);
  Array2<T> out(osize);
  ind_t j{0};
  for (auto v: i.valItr()) out.set(j++, this->view(static_cast<ind_t>(v)));
  return out;
}

template<class T>
template<class I>
std::enable_if_t<std::is_integral_v<I>, Array2<T>>
Array2<T>::extract(const std::vector<I>& i) const {
  for (auto x: i) if (!(0 <= x && static_cast<ind_t>(x) < _shape[0]))
    throw std::runtime_error("Array2 extract index must be in range");
  shape_t osize{_shape};
  osize[0] = static_cast<ind_t>(i.size());
  Array2<T> out(osize);
  for (ind_t j=0; j<osize[0]; ++j) out.set(j, this->view(static_cast<ind_t>(i[j])));
  return out;
}

template<class T>
template<class I, size_t Nel>
std::enable_if_t<std::is_integral_v<I>, Array2<T>>
Array2<T>::extract(const std::array<I,Nel>& i) const
{
  for (auto x: i) if (!(0 <= x && static_cast<ind_t>(x) < _shape[0]))
    throw std::runtime_error("Array2 extract index must be in range");
  shape_t osize{_shape};
  osize[0] = static_cast<ind_t>(Nel);
  Array2<T> out(osize);
  for (ind_t j=0; j<osize[0]; ++j)
  {
    out.set(j, this->view(static_cast<ind_t>(i[j])));
  }
  return out;
}

// template<class T>
// template<class I, size_t Nel>
// std::enable_if_t<std::is_integral_v<I>, Array<T>>
// Array2<T>::extract(const std::vector<std::array<I,Nel>>& i) const
// {
//   for (auto a: i) for (auto x: a) if (!(0 <= x && static_cast<ind_t>(x) < _shape[0]))
//     throw std::runtime_error("Array2 extract index must be in range");
//   brille::shape_t osh{static_cast<ind_t>(i.size()), static_cast<ind_t>(Nel)};
//   for (ind_t i=1; i<this->ndim(); ++i)
//     osh.push_back(_shape[i]);
//   Array<T> out(osh);
//   shape_t xi{_shape};
//   for (ind_t a=0; a<i.size(); ++a)
//   {
//     osh[0] = a;
//     for (ind_t b=0; b<Nel; ++b)
//     {
//       osh[1] = b;
//       xi[0] = static_cast<ind_t>(i[a][b]);
//       for (auto sub: SubIt2(_shape, xi))
//       {
//         for (ind_t j=0; j<_shape.size(); ++j) osh[j+1]=sub[j];
//         out[sub] = _data[this->s2l_d(xi)];
//       }
//     }
//   }
//   return out;
// }

template<class T>
Array2<T>
Array2<T>::extract(const Array2<bool>& i) const {
  shape_t isize = i.shape();
  if (isize[0] > _shape[0])
    throw std::runtime_error("Boolean Array2 extraction requires no more bools than the first Array2 dimension");
  shape_t osize{_shape};
  osize[0] = i.count();
  Array2<T> out(osize);
  ind_t n = isize[0] < _shape[0] ? isize[0] : _shape[0];
  ind_t idx{0};
  for (ind_t j=0; j<n; ++j) if (i[j]) out.set(idx++, this->view(j));
  return out;
}

template<class T>
Array2<T>
Array2<T>::extract(const std::vector<bool>& i) const {
  if (i.size() > _shape[0])
    throw std::runtime_error("Boolean Array2 extraction requires no more bools than the first Array2 dimension");
  shape_t osize{_shape};
  auto count = std::count(i.begin(), i.end(), true);
  osize[0] = static_cast<ind_t>(count);
  Array2<T> out(osize);
  ind_t n = i.size() < _shape[0] ? i.size() : _shape[0];
  ind_t idx{0};
  for (ind_t j=0; j<n; ++j) if (i[j]) out.set(idx++, this->view(j));
  return out;
}

template<class T>
bool Array2<T>::set(const ind_t i, const Array2<T>& in){
  // we might be able to do this better/faster if we both *this and in
  // have the same strides_ and we account for any offset. For now calculate
  // the 'hard' way no matter what:
  shape_t inshape = in.shape();
  if (_shape[1] != inshape[1])
    throw std::runtime_error("Set requires equal dimensions beyond the first dimension");
  inshape[0] = 1u; // this *should* be the case anyway
  for(auto sub: SubIt2<ind_t>(inshape)){
    const T& tmp = in[sub];
    sub[0] = i; // [0,0,…,0],[0,0,…,1],…,[0,1,…,0],… to [i,0,…,0],[i,0,…,1],…,[i,1,…,0],…
    _data[this->s2l_d(sub)] = tmp;
  }
  return true;
}

template<class T>
template<class R>
bool Array2<T>::set(const ind_t i, const Array2<R>& in){
  // we might be able to do this better/faster if we both *this and in
  // have the same strides_ and we account for any offset. For now calculate
  // the 'hard' way no matter what:
  shape_t inshape = in.shape();
  if (_shape[1] != inshape[1])
    throw std::runtime_error("Set requires equal dimensions beyond the first dimension");
  inshape[0] = 1u; // this *should* be the case anyway
  for(auto sub: SubIt2<ind_t>(inshape)){
    const R& tmp = in[sub];
    sub[0] = i; // [0,0,…,0],[0,0,…,1],…,[0,1,…,0],… to [i,0,…,0],[i,0,…,1],…,[i,1,…,0],…
    _data[this->s2l_d(sub)] = static_cast<T>(tmp);
  }
  return true;
}

template<class T>
bool Array2<T>::set(const ind_t i, const std::vector<T>& in){
  if (this->numel() != _shape[0]*in.size())
    throw std::runtime_error("Set requires the correct number of elements");
  shape_t tsize = this->shape();
  tsize[0] = 1u;
  // in is (hopefully) a row-ordered linear indexing
  size_t idx{0};
  for (auto sub: SubIt2<ind_t>(tsize)){
    sub[0] = i;
    _data[this->s2l_d(sub)] = in[idx++];
  }
  return true;
}

template<class T>
template<size_t Nel>
bool Array2<T>::set(const ind_t i, const std::array<T, Nel>& in){
  if (this->numel() != _shape[0]*Nel)
    throw std::runtime_error("Set requires the correct number of elements");
  shape_t tsize = this->shape();
  tsize[0] = 1u;
  // in is (hopefully) a row-ordered linear indexing
  size_t idx{0};
  for (auto sub: SubIt2<ind_t>(tsize)){
    sub[0] = i;
    _data[this->s2l_d(sub)] = in[idx++];
  }
  return true;
}

template<class T>
T Array2<T>::set(const shape_t& sub, const T in){
  auto ind = this->s2l_d(sub);
  _data[ind] = in;
  return _data[ind];
}

template<class T>
std::string Array2<T>::to_string() const{
  if (this->_num == 0) return std::string("Unallocated Array2");
  size_t width{0};
  for (ind_t i=0; i<_num; ++i){
    size_t l = my_to_string(_data[i]).size();
    if (l > width) width = l;
  }
  std::string out;
  std::vector<bool> isin(this->ndim(), false);
  bool preamble{false};
  size_t ndim = _shape.size();
  for (auto sub: this->subItr()){
    if (!preamble){
      for (ind_t i=0; i<ndim-1; ++i)
      if (sub[i] == 0 && !isin[i]){
        out += "[";
        isin[i] = true;
      }
      else
        out += " ";
      out += "["; // for i=ndim-1;
      preamble=true;
    }

    out += my_to_string(_data[this->s2l_d(sub)], width);
    // out += std::to_string(_data[this->s2l_d(sub)]);
    if (sub[ndim-1]+1 < _shape[ndim-1]){
      out += ", ";
    } else {
      for (ind_t i=0; i<ndim-1; ++i){
        bool isend{true};
        for (ind_t k=i; k<ndim; k++) isend &= sub[k]+1 == _shape[k];
        if (isend){
          out += "]";
          isin[i] = false;
        }
      }
      out += "],"; // ']' for i=ndim-1;
      size_t nret = std::count(isin.begin(), isin.end(), false);
      for (size_t n=0; n<nret; ++n) out += "\n";
      preamble=false;
    }
  }
  for (size_t i=0; i<ndim; ++i) out.pop_back(); // remove the trailing '\n's
  out.pop_back(); // and the trailing ,
  out += "\n";
  return out;
}

template<class T>
std::string Array2<T>::to_string(const ind_t i) const {
  auto out = this->view(i).to_string();
  out.pop_back(); // remove the trailing \n
  return out;
}

template<class T>
Array2<T>&
Array2<T>::reshape(const shape_t& ns){
  ind_t num = this->size_from_shape(ns);
  info_update_if(num > _num, "Array2::reshape only intended for equal-element number changes.");
  // assert( num <= _num );
  if (!this->is_contiguous())
    throw std::runtime_error("Array2::reshape does not work for strided arrays");
  _shape = ns;
  _stride = this->calculate_stride(ns);
  return *this;
}

/*! \brief Modify the number of arrays

Allocate a new T* and copy the overlapping region. The new shape
*must* have the same number of dimensions as the old shape.

Resizing an array will decouple it from any other arrays sharing their
underlying data. Such an operation should not be taken lightly.

@param ns the new shape
@param init an optional initialization value for any non-overlapping region
*/
template<class T>
Array2<T>&
Array2<T>::resize(const shape_t& ns, const T init) {
  ind_t nnum = this->size_from_shape(ns);
  // determine which axes are changing and how
  std::vector<bool> shrinking, growing;
  for (size_t i=0; i<ns.size(); ++i){
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
  T* nd = new T[nnum]();
  // fill the whole thing instead of trying to figure out difference indices
  if (grow) std::fill(nd, nd+nnum, init);
  // loop over the subscripted incides which fit inside both arrays
  shape_t inside = shrink ? ns : _shape;
  // find the new stride (maintaining row or column order)
  shape_t nt = calculate_stride(ns);
  // ensure our simple indexing will work
  if (_shift > 0)
    throw std::runtime_error("Resizing only works for zero-offset Array2s at present. Please extend.");
  // copy into the new array if the old array isn't the null pointer
  if (_data != nullptr && nd != nullptr)
    for (auto x : SubIt2(inside)) nd[sub2lin(x, nt)] = _data[sub2lin(x, _stride)];
  // delete the old _data if we hold the only reference
  if (_own && _ref.use_count()==1 && _data != nullptr) delete[] _data;
  // update data, number, ownership, reference, size, stride, and mutability
  _data = nd;
  _num = nnum;
  _own = true;
  _ref = std::make_shared<char>(); // resizing in place can not change the template parameter
  _shape = ns;
  _stride = nt;
  _mutable = true;
  // pass back a reference to this object
  return *this;
}
template<class T> template<class I>
Array2<T>&
Array2<T>::resize(const I ns, const T init) {
  shape_t nshape{_shape};
  if (nshape.size()>0)
    nshape[0] = static_cast<ind_t>(ns);
  else{
    throw std::runtime_error("Resizing a null-shaped Array2 not supported (yet)");
    //nshape.push_back(static_cast<ind_t>(ns));
  }
  return this->resize(nshape, init);
}

template<class T>
Array2<T>&
Array2<T>::append(const ind_t dim, const Array2<T>& extra) {
  assert(this != &extra);
  ind_t ndim = this->ndim();
  assert(dim<ndim);
  auto eshape = extra.shape();
  for (size_t i=0; i<ndim; ++i) if (i!=dim && eshape[i]!=_shape[i])
    throw std::runtime_error("Incompatible Array2 to append");
  // the new expanded size only changes along dim
  eshape[dim] += _shape[dim];
  // stash this original shape along dim for use in the loop:
  ind_t tshapedim = _shape[dim];
  // and resizing the Array2 to the new size copies existing entries
  this->resize(eshape);
  // so we need only iterate over the extra Array2
  for (auto x : SubIt2(extra.shape())){
    // offsetting its subscript to the new position
    auto y = x;
    y[dim] += tshapedim; // offset by the original shape along the dimension
    // and copying the contents
    (*this)[y] = extra[x];
  }
  return *this;
}

template<class T>
bool Array2<T>::all(const ind_t n) const{
  ind_t nml = this->numel();
  ind_t count = (n > 0 && n < nml) ? n : nml;
  for (ind_t i=0; i<count; ++i) if (!_data[this->l2l_d(i)]) return false;
  return true;
}
template<class T>
bool Array2<T>::any(const ind_t n) const{
  ind_t nml = this->numel();
  ind_t count = (n > 0 && n < nml) ? n : nml;
  for (ind_t i=0; i<count; ++i) if (_data[this->l2l_d(i)]) return true;
  return false;
}
template<class T>
bool Array2<T>::all(const T val, const ind_t n) const{
  ind_t nml = this->numel();
  ind_t count = (n > 0 && n < nml) ? n : nml;
  for (ind_t i=0; i<count; ++i) if (val != _data[this->l2l_d(i)]) return false;
  return true;
}
template<class T>
bool Array2<T>::any(const T val, const ind_t n) const{
  ind_t nml = this->numel();
  ind_t count = (n > 0 && n < nml) ? n : nml;
  for (ind_t i=0; i<count; ++i) if (val == _data[this->l2l_d(i)]) return true;
  return false;
}

template<class T>
ind_t
Array2<T>::count(const ind_t n) const{
  ind_t el = this->numel();
  ind_t no = (n > 0 && n < el) ? n : el;
  ind_t count{0};
  for (ind_t i=0; i<no; ++i) if (_data[this->l2l_d(i)]) ++count;
  return count;
}
template<class T>
ind_t
Array2<T>::count(const T val, const ind_t n) const{
  ind_t el = this->numel();
  ind_t no = (n > 0 && n < el) ? n : el;
  ind_t count{0};
  for (ind_t i=0; i<no; ++i) if (val == _data[this->l2l_d(i)]) ++count;
  return count;
}

template<class T>
ind_t
Array2<T>::first(const ind_t n) const{
  ind_t el = this->numel();
  ind_t no = (n > 0 && n < el) ? n : el;
  for (ind_t i=0; i<no; ++i) if (_data[this->l2l_d(i)]) return i;
  return no;
}
template<class T>
ind_t
Array2<T>::first(const T val, const ind_t n) const{
  ind_t el = this->numel();
  ind_t no = (n > 0 && n < el) ? n : el;
  for (ind_t i=0; i<no; ++i) if (val == _data[this->l2l_d(i)]) return i;
  return no;
}
template<class T>
ind_t
Array2<T>::last(const ind_t n) const{
  ind_t el = this->numel();
  ind_t no = (n > 0 && n < el) ? n : el;
  for (ind_t i=no; i--;) if (_data[this->l2l_d(i)]) return i;
  return no;
}
template<class T>
ind_t
Array2<T>::last(const T val, const ind_t n) const{
  ind_t el = this->numel();
  ind_t no = (n > 0 && n < el) ? n : el;
  for (ind_t i=no; i--;) if (val == _data[this->l2l_d(i)]) return i;
  return no;
}

#define ARRAY_ELEMENTWISE_INT_TRANSFORM(X) template<class T>\
Array2<int>\
Array2<T>:: X () const {\
  Array2<int> out(_shape, _stride);\
  for (auto x: out.subItr())\
    out[x] = static_cast<int>(std:: X (_data[this->s2l_d(x)]));\
  return out;\
}
ARRAY_ELEMENTWISE_INT_TRANSFORM(round)
ARRAY_ELEMENTWISE_INT_TRANSFORM(floor)
ARRAY_ELEMENTWISE_INT_TRANSFORM(ceil)
#undef ARRAY_ELEMENTWISE_INT_TRANSFORM

template<class T>
Array2<T>
Array2<T>::sum(const ind_t dim) const{
  assert(dim < _shape.size());
  shape_t osize = this->shape();
  osize[dim] = 1u;
  shape_t ostride = this->calculate_stride(osize); // preserve orderness
  Array2<T> out(osize, ostride);
  for (auto oidx: out.subItr()) {
    T tmp{0};
    auto idx = oidx;
    for (ind_t i=0; i<_shape[dim]; ++i) {
      idx[dim] = i;
      tmp += _data[this->s2l_d(idx)];
    }
    out[oidx] = tmp;
  }
  return out;
}

template<class T>
Array2<T>
Array2<T>::prod(const ind_t dim) const{
  assert(dim < _shape.size());
  shape_t osize = this->shape();
  osize[dim] = 1u;
  shape_t ostride = this->calculate_stride(osize); // preserve orderness
  Array2<T> out(osize, ostride);
  for (auto oidx: out.subItr()) {
    T tmp{1};
    auto idx = oidx;
    for (ind_t i=0; i<_shape[dim]; ++i) {
      idx[dim] = i;
      tmp *= _data[this->s2l_d(idx)];
    }
    out[oidx] = tmp;
  }
  return out;
}

template<class T>
Array2<T>
Array2<T>::min(const ind_t dim) const{
  assert(dim < _shape.size());
  shape_t osize = this->shape();
  osize[dim] = 1u;
  shape_t ostride = this->calculate_stride(osize); // preserve orderness
  Array2<T> out(osize, ostride);
  for (auto oidx: out.subItr()) {
    T tmp{(std::numeric_limits<T>::max)()};
    auto idx = oidx;
    for (ind_t i=0; i<_shape[dim]; ++i) {
      idx[dim] = i;
      T tst = _data[this->s2l_d(idx)];
      if (tst < tmp) tmp = tst;
    }
    out[oidx] = tmp;
  }
  return out;
}

template<class T>
Array2<T>
Array2<T>::max(const ind_t dim) const{
  size_t ndim = _shape.size();
  assert(dim < ndim);
  shape_t osize = this->shape();
  osize[dim] = 1u;
  shape_t ostride = this->calculate_stride(osize); // preserve orderness
  // the stride along any dimension can not be zero. protect against it:
  for (ind_t i=0; i<ndim; ++i) if (ostride[i] < 1u) ostride[i] = 1u;

  Array2<T> out(osize, ostride);
  for (auto oidx: out.subItr()) {
    T tmp{std::numeric_limits<T>::lowest()};
    auto idx = oidx;
    for (ind_t i=0; i<_shape[dim]; ++i) {
      idx[dim] = i;
      T tst = _data[this->s2l_d(idx)];
      if (tst > tmp) tmp = tst;
    }
    out[oidx] = tmp;
  }
  return out;
}

template<class T>
T Array2<T>::sum() const{
  T out{0};
  for (auto x : this->subItr()) out += _data[this->s2l_d(x)];
  return out;
}

template<class T>
T Array2<T>::prod() const{
  T out{1};
  for (auto x : this->subItr()) out *= _data[this->s2l_d(x)];
  return out;
}

template<class T>
template<class R, size_t Nel>
bool Array2<T>::match(const ind_t i, const ind_t j, const std::array<R,Nel>& rot, const int order) const {
  assert(this->ndim() == 2u); // only defined for 2-D Array2s
  assert(Nel == _shape[1]*_shape[1]);
  ind_t n = _shape[1];
  std::vector<T> tmp(n,T(0));
  for (ind_t k=0; k<n; ++k) tmp[k] = this->val(j,k);
  brille::Comparer<T,T> eq(brille::cmp::eq);
  if (order < 0){
    int o{0};
    do{
      //check against the current tmp vector whether this order rotation has moved j to i
      if (eq(n, this->ptr(i), _stride[1], tmp.data(), 1u)) return true;
      // otherwise apply the rotation (again)
      brille::utils::mul_mat_vec_inplace(n, rot.data(), tmp.data());
    } while (o++ < std::abs(order));
    return false;
  } else {
    // rotate exactly order times
    for (int o=0; o<order; ++o)
      brille::utils::mul_mat_vec_inplace(n, rot.data(), tmp.data());
    return eq(n, this->ptr(i), _stride[1], tmp.data(), 1u);
  }
}

template<class T>
bool Array2<T>::match(const ind_t i, const ind_t j, brille::ops op, T val) const {
  auto ai= this->view(i);
  auto aj= this->view(j);
  brille::Comparer<T,T> neq(brille::cmp::neq);
  ind_t no = this->numel()/_shape[0];
  switch (op){
    case brille::ops::plus:
    for (ind_t k=0; k<no; ++k) if (neq(ai[k], aj[k]+val)) return false;
    break;
    case brille::ops::minus:
    for (ind_t k=0; k<no; ++k) if (neq(ai[k], aj[k]-val)) return false;
    break;
    case brille::ops::times:
    for (ind_t k=0; k<no; ++k) if (neq(ai[k], aj[k]*val)) return false;
    break;
    case brille::ops::rdiv:
    for (ind_t k=0; k<no; ++k) if (neq(ai[k], aj[k]/val)) return false;
    break;
    case brille::ops::ldiv:
    for (ind_t k=0; k<no; ++k) if (neq(ai[k], val/aj[k])) return false;
    break;
    default:
    throw std::runtime_error(std::string("Unhandled operator ")+brille::to_string(op));
  }
  return true;
}

template<class T>
bool Array2<T>::all(const brille::cmp expr, const T val) const {
  if (brille::cmp::le_ge == expr)
    return this->all(brille::cmp::le, val) || this->all(brille::cmp::ge, val);
  ind_t no = this->numel();
  brille::Comparer<T,T> op(expr);
  for (ind_t k=0; k<no; ++k) if(!op(_data[this->l2l_d(k)], val)) return false;
  return true;
}
template<class T>
bool Array2<T>::any(const brille::cmp expr, const T val) const {
  return this->first(expr, val) < this->numel();
}
template<class T>
ind_t
Array2<T>::first(const brille::cmp expr, const T val) const {
  ind_t no = this->numel();
  brille::Comparer<T,T> op(expr);
  for (ind_t k=0; k<no; ++k) if(op(_data[this->l2l_d(k)], val)) return k;
  return no;
}
template<class T>
ind_t
Array2<T>::last(const brille::cmp expr, const T val) const {
  ind_t no = this->numel();
  brille::Comparer<T,T> op(expr);
  for (ind_t k=no; k--;) if(op(_data[this->l2l_d(k)], val)) return k;
  return no;
}
template<class T>
ind_t
Array2<T>::count(const brille::cmp expr, const T val) const {
  ind_t no = this->numel();
  brille::Comparer<T,T> op(expr);
  ind_t cnt{0};
  for (ind_t k=0; k<no; ++k) if(op(_data[this->l2l_d(k)], val)) ++cnt;
  return cnt;
}

template<class T>
Array2<bool>
Array2<T>::is(const brille::cmp expr, const T val) const {
  Array2<bool> out(_shape, _stride, true);
  ind_t no = this->numel();
  brille::Comparer<T,T> op(expr);
  for (ind_t k=0; k<no; ++k) out[k] = op(_data[this->l2l_d(k)], val);
  return out;
}

template<class T>
std::vector<ind_t>
Array2<T>::find(const brille::cmp expr, const T val) const {
  Array2<bool> this_is = this->is(expr, val);
  ind_t no = this->numel();
  std::vector<ind_t> out;
  for (ind_t k=0; k<no; ++k) if (this_is[k]) out.push_back(k);
  return out;
}

template<class T>
template<class R>
Array2<bool>
Array2<T>::is(const brille::cmp expr, const Array2<R>& that) const {
  // To handle singleton-dimension broadcasting, this function needs to be split
  auto tsize = that.shape();
  if (!std::equal(_shape.begin(), _shape.end(), tsize.begin(),
      [](ind_t a, ind_t b){return 1==b || a==b;})){
    std::string msg = "An Array2 with size ( ";
    for (auto x: tsize) msg += std::to_string(x) + " ";
    msg += ") can not be broadcast to match one with size ( ";
    for (auto x: _shape) msg += std::to_string(x) + " ";
    msg += ").";
    throw std::runtime_error(msg);
  }
  Array2<bool> out(_shape, _stride, true);
  if (std::equal(_shape.begin(), _shape.end(), tsize.begin())){
    // No broadcast
    brille::Comparer<T,R> op(expr);
    // no guarantees about same stride, so use subscript iterator:
    for (auto sub: this->subItr()) out[sub] = op(_data[this->s2l_d(sub)], that[sub]);
  } else {
    // Broadcast
    size_t ndim = _shape.size();
    for (ind_t i=1; i<ndim; ++i) if (_shape[i]!=tsize[i])
      throw std::runtime_error("Broadcasting beyond the first dimension requires viewing beyond the first dimension!");
    for (ind_t i=0; i<_shape[0]; ++i)
      out.set(i, this->view(i).is(expr, that));
  }
  return out;
}

template<class T>
template<class R>
std::vector<bool>
Array2<T>::is(const brille::cmp expr, const std::vector<R>& val) const{
  assert(val.size() == _shape[1]);
  std::vector<bool> out;
  out.reserve(_shape[0]);
  brille::Comparer<T,R> op(expr);
  for (ind_t i=0; i<_shape[0]; ++i)
    out.push_back(op(_shape[1], this->ptr(i), _stride[1], val.data(), 1u));
  return out;
}

template<class T>
template<class R>
bool Array2<T>::is(const Array2<R>& that) const {
  return this->is(brille::cmp::eq, that).all(true);
}

template<class T>
std::vector<bool>
Array2<T>::is_unique() const{
  if (_shape[0] < 1u) return std::vector<bool>();
  std::vector<bool> out(1u, true);
  out.reserve(_shape[0]);
  brille::Comparer<T,T> op(brille::cmp::neq);
  ind_t sz{_shape[1]}, st{_stride[1]};
  for (ind_t i=1; i<_shape[0]; ++i){
    bool isu{true};
    for (ind_t j=0; j<i; ++j){
      isu &= op(sz, this->ptr(i), st, this->ptr(j), st);
      if (!isu) break;
    }
    out.push_back(isu);
  }
  return out;
}
template<class T>
std::vector<ind_t>
Array2<T>::unique_idx() const{
  if (_shape[0] < 1u) return std::vector<ind_t>();
  std::vector<ind_t> out(1u, 0u);
  out.reserve(_shape[0]);
  brille::Comparer<T,T> op(brille::cmp::eq);
  ind_t sz{_shape[1]}, st{_stride[1]};
  for (ind_t i=1; i<_shape[0]; ++i){
    ind_t idx{i};
    for (ind_t j=0; j<i; ++j)
    if(j==out[j] && op(sz, this->ptr(i), st, this->ptr(j), st)){
      idx=j;
      break;
    }
    out.push_back(idx);
  }
  return out;
}
template<class T>
Array2<T>
Array2<T>::unique() const {
  std::vector<bool> isu = this->is_unique();
  size_t u_count = std::count(isu.begin(), isu.end(), true);
  shape_t osize{_shape};
  osize[0] = static_cast<ind_t>(u_count);
  Array2<T> out(osize, _stride);
  for (ind_t i=0,u=0; i<_shape[0]; ++i) if (isu[i]) out.set(u++, this->view(i));
  return out;
}
template<class T>
T Array2<T>::dot(ind_t i, ind_t j) const {
  assert(i<_shape[0]);
  assert(j<_shape[0]);
  ind_t sz{_shape[1]}, st{_stride[1]};
  std::vector<T> prods(sz,0);
  brille::RawBinaryOperator<T> prod(brille::ops::times);
  // perform elementwise multiplication on the views of i and j
  prod(sz, prods.data(), 1u, this->ptr(i), st, this->ptr(j), st);
  // find the sum of the products
  // std::reduce can do this in parallel but gcc<v9.3 does not implement all of C++17
  return std::accumulate(prods.begin(), prods.end(), T(0));
}
template<class T>
T Array2<T>::norm(ind_t i) const {
  return std::sqrt(this->dot(i,i));
}

template<class T>
Array2<T>
Array2<T>::operator-() const{
  Array2<T> neg(_shape, _stride);
  for (auto s: this->subItr()) neg[s] = -_data[this->s2l_d(s)];
  return neg;
}


template<class T>
template<class I, typename>
void Array2<T>::permute(std::vector<I>& p){
  std::vector<I> s=p, o(_shape[0]);
  std::iota(o.begin(), o.end(), 0u);
  std::sort(s.begin(), s.end());
  debug_update_if(!std::includes(o.begin(),o.end(),s.begin(),s.end()),"The permutation vector ",p," is invalid. Expected permutation of ",o);
  // get the inverse permutation to enable element swapping:
  for (size_t i=0; i<p.size(); ++i) s[p[i]] = i;
  // perform the swaps until we have no more to do:
  for (size_t i=0; i<_shape[0];){
    if (s[i]!=i){
      this->swap(i,s[i]);
      std::swap(s[i], s[s[i]]);
    } else{
      ++i;
    }
  }
  debug_update_if(!std::is_sorted(s.begin(),s.end()), "Undoing the permutation ",p," failed. End result is ",s);
}

template<class T>
bool Array2<T>::swap(const ind_t a, const ind_t b){
  assert(a<_shape[0] && b<_shape[0]);
  shape_t sub{_shape};
  sub[0] = a; // fix the 0th index to 'a' in the iterator
  for (auto aidx: this->subItr(sub)){
    // for each of the [a,...] subscripted indices, construct [b,...]
    shape_t bidx{aidx};
    bidx[0] = b;
    // precalculate the offset/stride/shaped subscripts:
    auto sa = this->s2l_d(aidx);
    auto sb = this->s2l_d(bidx);
    // perform the actual swap:
    T _store_ = _data[sa];
    _data[sa] = _data[sb];
    _data[sb] = _store_;
  }
  return true;
}
template<class T>
bool Array2<T>::swap(ind_t i, ind_t a, ind_t b){
  assert(i < _shape[0] && a < _shape[1] && b < _shape[1]);
  auto la{this->ij2l_d(i,a)}, lb{this->ij2l_d(i,b)};
  T _store_ = _data[la];
  _data[la] = _data[lb];
  _data[lb] = _store_;
  return true;
}

template<class T>
std::vector<T> Array2<T>::to_std() const {
  std::vector<T> out;
  for (auto x: this->valItr()) out.push_back(x);
  return out;
}

template<class T>
T* Array2<T>::ptr(const ind_t i0){
  assert(i0 < _shape[0]);
  return _data + ij2l_d(i0, 0u);
}
template<class T>
T* Array2<T>::ptr(const ind_t i0, const ind_t j0){
  assert(i0 < _shape[0] && j0 < _shape[1]);
  return _data + ij2l_d(i0, j0);
}
template<class T>
T* Array2<T>::ptr(const shape_t& p){
  assert(p[0]<_shape[0] && p[1]<_shape[1]);
  return _data + s2l_d(p);
}

template<class T>
const T* Array2<T>::ptr(const ind_t i0) const {
  assert(i0 < _shape[0]);
  return _data + ij2l_d(i0, 0u);
}
template<class T>
const T* Array2<T>::ptr(const ind_t i0, const ind_t j0) const {
  assert(i0 < _shape[0] && j0 < _shape[1]);
  return _data + ij2l_d(i0, j0);
}
template<class T>
const T* Array2<T>::ptr(const shape_t& p) const {
  assert(p[0]<_shape[0] && p[1]<_shape[1]);
  return _data + s2l_d(p);
}


template<class T>
T& Array2<T>::val(const ind_t i0){
  assert(i0 < _shape[0]);
  return _data[ij2l_d(i0, 0u)];
}
template<class T>
T& Array2<T>::val(const ind_t i0, const ind_t j0){
  assert(i0 < _shape[0] && j0 < _shape[1]);
  return _data[ij2l_d(i0, j0)];
}
template<class T>
T& Array2<T>::val(const shape_t& p){
  assert(p[0]<_shape[0] && p[1]<_shape[1]);
  return _data[s2l_d(p)];
}
template<class T>
template<class I>
T& Array2<T>::val(std::initializer_list<I> l){
  shape_t idx{l};
  return this->val(idx);
}

template<class T>
const T& Array2<T>::val(const ind_t i0) const {
  assert(i0 < _shape[0]);
  return _data[ij2l_d(i0, 0u)];
}
template<class T>
const T& Array2<T>::val(const ind_t i0, const ind_t j0) const {
  assert(i0 < _shape[0] && j0 < _shape[1]);
  return _data[ij2l_d(i0, j0)];
}
template<class T>
const T& Array2<T>::val(const shape_t& p) const {
  assert(p[0]<_shape[0] && p[1]<_shape[1]);
  return _data[s2l_d(p)];
}
template<class T>
template<class I>
const T& Array2<T>::val(std::initializer_list<I> l) const {
  shape_t idx{l};
  return this->val(idx);
}

template<class T>
Array2<T> Array2<T>::contiguous_copy() const {
  if (this->is_contiguous()) return Array2<T>(*this);
  Array2<T> out(_shape, this->calculate_stride(_shape));
  for (auto x : this->subItr()) out[x] = (*this)[x];
  return out;
}

template<class T>
Array2<T> Array2<T>::contiguous_row_ordered_copy() const {
  if (this->is_row_ordered() && this->is_contiguous()) return Array2<T>(*this);
  Array2<T> out(_shape, this->calculate_stride(_shape));
  for (auto x : this->subItr()) out[x] = (*this)[x];
  return out;
}

template<class T>
static void mutable_check(const Array2<T>& a){
  if (!a.ismutable())
  throw std::runtime_error("Immutable Array2 objects can not be modified!");
}

#define SCALAR_INPLACE_OP(X) template<class T>\
Array2<T>& Array2<T>::operator X (const T& val){\
  mutable_check(*this);\
  for (auto& v: this->valItr()) v X val;\
  return *this;\
}
SCALAR_INPLACE_OP(+=)
SCALAR_INPLACE_OP(-=)
SCALAR_INPLACE_OP(*=)
SCALAR_INPLACE_OP(/=)
#undef SCALAR_INPLACE_OP
// template<class T>
// Array2<T>&
// Array2<T>::operator+=(const T& val){
//   mutable_check(*this);
//   for (auto& v: this->valItr()) v += val;
//   return *this;
// }
// template<class T>
// Array2<T>&
// Array2<T>::operator-=(const T& val){
//   mutable_check(*this);
//   for (auto& v: this->valItr()) v -= val;
//   return *this;
// }
// template<class T>
// Array2<T>&
// Array2<T>::operator*=(const T& val){
//   mutable_check(*this);
//   auto itr = this->valItr();
//   for (auto & v: itr){
//     v *= val;
//   }
//   return *this;
// }
// template<class T>
// Array2<T>&
// Array2<T>::operator/=(const T& val){
//   mutable_check(*this);
//   for (auto& v: this->valItr()) v /= val;
//   return *this;
// }

template<class T>
bool broadcast_shape_check(const std::array<T,2>& a, const std::array<T,2>&b){
  bool ok = (a[0]==b[0] && a[1]==b[1]);
  if (!ok){
    std::string msg = "In place broadcasting is not possible for { ";
    for (auto x: a) msg += std::to_string(x) + " ";
    msg += "} shaped Array2 and { ";
    for (auto x: b) msg += std::to_string(x) + " ";
    msg += "} shaped operand";
    throw std::runtime_error(msg);
  }
  return ok;
}

#define ARRAY_INPLACE_OP(X) template<class T> template<class R>\
Array2<T>& Array2<T>::operator X (const Array2<R>& val){\
  mutable_check(*this);\
  auto itr = this->broadcastItr(val);\
  broadcast_shape_check(_shape, itr.shape());\
  for (auto [ox, tx, vx]: itr) _data[this->s2l_d(tx)] X val[vx];\
  return *this;\
}
ARRAY_INPLACE_OP(+=)
ARRAY_INPLACE_OP(-=)
ARRAY_INPLACE_OP(*=)
ARRAY_INPLACE_OP(/=)
#undef ARRAY_INPLACE_OP
