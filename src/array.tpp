template<class T, class P>
Array<T,P> Array<T,P>::view() const{
  return Array<T,P>(_data,_num,_own,_ref,_offset,_shape,_stride,false);
}
template<class T, class P>
Array<T,P> Array<T,P>::view(const ind_t i) const {
  if (i<_shape[0]){
    shape_t osize{_shape};
    shape_t oofst{_offset};
    oofst[0] += i;
    osize[0] = 1;
    // osize[0] = single ? 1 : osize[0] - 1;
    return Array<T,P>(_data, _num, _own, _ref, oofst, osize, _stride, false);
  }
  throw std::runtime_error("Array index too large");
}

template<class T, class P>
Array<T,P> Array<T,P>::view(const ind_t i, const ind_t j) const {
  if (i<j && i<_shape[0] && j<=_shape[0]){
    shape_t osize{_shape};
    shape_t oofst{_offset};
    oofst[0] += i;
    osize[0] = j-i;
    return Array<T,P>(_data, _num, _own, _ref, oofst, osize, _stride, false);
  }
  throw std::runtime_error("Array view indexing error");
}

template<class T, class P>
Array<T,P> Array<T,P>::view(const shape_t& v) const {
  shape_t osize{_shape};
  shape_t oofst{_offset};
  for (size_t i=0; i<_shape.size(); ++i) if (v[i] < _shape[i]) {
    oofst[i] += v[i];
    osize[i] -= v[i];
  }
  return Array<T,P>(_data, _num, _own, _ref,oofst, osize, _stride, false);
}

template<class T, class P>
Array<T,ref_ptr_t>
Array<T,P>::extract(const ind_t i) const{
  return this->view(i).decouple();
}

template<class T, class P>
template<class I, class IP>
std::enable_if_t<std::is_integral_v<I>, Array<T,ref_ptr_t>>
Array<T,P>::extract(const Array<I,IP>& i) const {
  if (i.numel() != i.size(0))
    throw std::runtime_error("Array extraction requires (N,{1,...,1}) shape");
  for (auto x: ArrayIt(i)) if (!(0 <= x && static_cast<ind_t>(x) < _shape[0]))
    throw std::runtime_error("Array extract index must be in range");
  shape_t osize{_shape};
  osize[0] = i.size(0);
  Array<T,P> out(osize);
  ind_t j{0};
  for (auto v: ArrayIt(i)) out.set(j++, this->view(static_cast<ind_t>(v)));
  return out;
}

template<class T, class P>
template<class I>
std::enable_if_t<std::is_integral_v<I>, Array<T,ref_ptr_t>>
Array<T,P>::extract(const std::vector<I>& i) const {
  for (auto x: i) if (!(0 <= x && static_cast<ind_t>(x) < _shape[0]))
    throw std::runtime_error("Array extract index must be in range");
  shape_t osize{_shape};
  osize[0] = static_cast<ind_t>(i.size());
  Array<T,P> out(osize);
  for (ind_t j=0; j<osize[0]; ++j) out.set(j, this->view(i[j]));
  return out;
}

template<class T, class P>
template<typename I, size_t Nel>
std::enable_if_t<std::is_integral_v<I>, Array<T,ref_ptr_t>>
Array<T,P>::extract(const std::array<I,Nel>& i) const
{
  for (auto x: i) if (!(0 <= x && static_cast<ind_t>(x) < _shape[0]))
    throw std::runtime_error("Array extract index must be in range");
  shape_t osize{_shape};
  osize[0] = static_cast<ind_t>(Nel);
  Array<T,P> out(osize);
  for (ind_t j=0; j<osize[0]; ++j)
  {
    out.set(j, this->view(static_cast<ind_t>(i[j])));
  }
  return out;
}

template<class T, class P>
template<typename I, size_t Nel>
std::enable_if_t<std::is_integral_v<I>, Array<T,ref_ptr_t>>
Array<T,P>::extract(const std::vector<std::array<I,Nel>>& i) const
{
  for (auto a: i) for (auto x: a) if (!(0 <= x && static_cast<ind_t>(x) < _shape[0]))
    throw std::runtime_error("Array extract index must be in range");
  shape_t osh{static_cast<ind_t>(i.size()), static_cast<ind_t>(Nel)};
  for (ind_t i=1; i<this->ndim(); ++i)
    osh.push_back(_shape[i]);
  Array<T,P> out(osh);
  shape_t xi{_shape};
  for (ind_t a=0; a<i.size(); ++a)
  {
    osh[0] = a;
    for (ind_t b=0; b<Nel; ++b)
    {
      osh[1] = b;
      xi[0] = static_cast<ind_t>(i[a][b]);
      for (auto sub: SubIt(_shape, xi))
      {
        for (ind_t j=0; j<_shape.size(); ++j) osh[j+1]=sub[j];
        out[sub] = _data[this->s2l_d(xi)];
      }
    }
  }
  return out;
}

template<class T, class P>
template<class RP>
Array<T,ref_ptr_t>
Array<T,P>::extract(const Array<bool,RP>& i) const {
  shape_t isize = i.shape();
  if (isize[0] > _shape[0])
    throw std::runtime_error("Boolean Array extraction requires no more bools than the first Array dimension");
  shape_t osize{_shape};
  osize[0] = i.count();
  Array<T,P> out(osize);
  ind_t n = isize[0] < _shape[0] ? isize[0] : _shape[0];
  ind_t idx{0};
  for (ind_t j=0; j<n; ++j) if (i[j]) out.set(idx++, this->view(j));
  return out;
}

template<class T, class P>
Array<T,ref_ptr_t>
Array<T,P>::extract(const std::vector<bool>& i) const {
  if (i.size() > _shape[0])
    throw std::runtime_error("Boolean Array extraction requires no more bools than the first Array dimension");
  shape_t osize{_shape};
  auto count = std::count(i.begin(), i.end(), true);
  osize[0] = static_cast<ind_t>(count);
  Array<T,P> out(osize);
  ind_t n = i.size() < _shape[0] ? i.size() : _shape[0];
  ind_t idx{0};
  for (ind_t j=0; j<n; ++j) if (i[j]) out.set(idx++, this->view(j));
  return out;
}

template<class T, class P>
template<class RP>
bool Array<T,P>::set(const ind_t i, const Array<T,RP>& in){
  // we might be able to do this better/faster if we both *this and in
  // have the same strides_ and we account for any offset. For now calculate
  // the 'hard' way no matter what:
  shape_t inshape = in.shape();
  size_t ndim = _shape.size();
  for (ind_t dim=1; dim<ndim; ++dim) if (_shape[dim] != inshape[dim])
    throw std::runtime_error("Set requires equal dimensions beyond the first dimension");
  inshape[0] = 1u; // this *should* be the case anyway
  for(auto sub: SubIt<ind_t>(inshape)){
    const T& tmp = in[sub];
    sub[0] = i; // [0,0,…,0],[0,0,…,1],…,[0,1,…,0],… to [i,0,…,0],[i,0,…,1],…,[i,1,…,0],…
    _data[this->s2l_d(sub)] = tmp;
  }
  return true;
}

template<class T, class P>
template<class R, class RP>
bool Array<T,P>::set(const ind_t i, const Array<R,RP>& in){
  // we might be able to do this better/faster if we both *this and in
  // have the same strides_ and we account for any offset. For now calculate
  // the 'hard' way no matter what:
  shape_t inshape = in.shape();
  size_t ndim = _shape.size();
  for (ind_t dim=1; dim<ndim; ++dim) if (_shape[dim] != inshape[dim])
    throw std::runtime_error("Set requires equal dimensions beyond the first dimension");
  inshape[0] = 1u; // this *should* be the case anyway
  for(auto sub: SubIt<ind_t>(inshape)){
    const R& tmp = in[sub];
    sub[0] = i; // [0,0,…,0],[0,0,…,1],…,[0,1,…,0],… to [i,0,…,0],[i,0,…,1],…,[i,1,…,0],…
    _data[this->s2l_d(sub)] = static_cast<T>(tmp);
  }
  return true;
}

template<class T, class P>
bool Array<T,P>::set(const ind_t i, const std::vector<T>& in){
  if (this->numel() != _shape[0]*in.size())
    throw std::runtime_error("Set requires the correct number of elements");
  shape_t tsize = this->shape();
  tsize[0] = 1u;
  // in is (hopefully) a row-ordered linear indexing
  size_t idx{0};
  for (auto sub: SubIt<ind_t>(tsize)){
    sub[0] = i;
    _data[this->s2l_d(sub)] = in[idx++];
  }
  return true;
}

template<class T, class P>
template<size_t Nel>
bool Array<T,P>::set(const ind_t i, const std::array<T, Nel>& in){
  if (this->numel() != _shape[0]*Nel)
    throw std::runtime_error("Set requires the correct number of elements");
  shape_t tsize = this->shape();
  tsize[0] = 1u;
  // in is (hopefully) a row-ordered linear indexing
  size_t idx{0};
  for (auto sub: SubIt<ind_t>(tsize)){
    sub[0] = i;
    _data[this->s2l_d(sub)] = in[idx++];
  }
  return true;
}

template<class T, class P>
T Array<T,P>::set(const shape_t& sub, const T in){
  auto ind = this->s2l_d(sub);
  _data[ind] = in;
  return _data[ind];
}

template<class T, class P>
std::string Array<T,P>::to_string() const{
  size_t width{0};
  for (ind_t i=0; i<_num; ++i){
    size_t l = my_to_string(_data[i]).size();
    if (l > width) width = l;
  }
  std::string out;
  std::vector<bool> isin(this->ndim(), false);
  bool preamble{false};
  size_t ndim = _shape.size();
  for (auto sub: SubIt<ind_t>(_shape)){
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

template<class T, class P>
std::string Array<T,P>::to_string(const ind_t i) const {
  return this->view(i).to_string();
}

template<class T, class P>
Array<T,P>&
Array<T,P>::reshape(const shape_t& ns){
  ind_t num = this->size_from_shape(ns);
  info_update_if(num > _num, "Array::reshape only intended for equal-element number changes.");
  // assert( num <= _num );
  if (!this->is_contiguous())
    throw std::runtime_error("Array::reshape does not work for strided arrays");
  _shape = ns;
  ind_t linoffset = sub2lin(_offset, _stride);
  _stride = this->calculate_stride(ns);
  _offset = lin2sub(linoffset, _stride);
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
template<class T, class P>
Array<T,P>&
Array<T,P>::resize(const shape_t& ns, const T init) {
  ind_t nnum = this->size_from_shape(ns);
  //if (nnum == _num) return this->reshape(ns);
  // otherwise resizing can not change dimensionality
  assert(_shape.size() == ns.size());
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
  ind_t linoffset = sub2lin(_offset, _stride);
  if (linoffset > 0)
    throw std::runtime_error("Resizing only works for zero-offset Arrays at present. Please extend.");
  // copy into the new array if the old array isn't the null pointer
  if (_data != nullptr && nd != nullptr)
    for (auto x : SubIt(inside)) nd[sub2lin(x, nt)] = _data[sub2lin(x, _stride)];
  // delete the old _data if we hold the only reference
  if (_own && _ref.use_count()==1 && _data != nullptr) delete[] _data;
  // update data, number, ownership, reference, size, stride, and mutability
  _data = nd;
  _num = nnum;
  _own = true;
  _ref = std::make_shared<P>(); // resizing in place can not change the template parameter
  _offset = shape_t(ns.size(), 0);
  _shape = ns;
  _stride = nt;
  _mutable = true;
  // pass back a reference to this object
  return *this;
}
template<class T, class P> template<class I>
Array<T,P>&
Array<T,P>::resize(const I ns, const T init) {
  shape_t nshape{_shape};
  if (nshape.size()>0)
    nshape[0] = static_cast<ind_t>(ns);
  else{
    throw std::runtime_error("Resizing a null-shaped Array not supported (yet)");
    //nshape.push_back(static_cast<ind_t>(ns));
  }
  return this->resize(nshape, init);
}

template<class T, class P>
template<class RP>
Array<T,P>&
Array<T,P>::append(const ind_t dim, const Array<T,RP>& extra) {
  assert(this != &extra);
  ind_t ndim = this->ndim();
  assert(extra.ndim() == ndim && dim<ndim);
  auto eshape = extra.shape();
  for (size_t i=0; i<ndim; ++i) if (i!=dim && eshape[i]!=_shape[i])
    throw std::runtime_error("Incompatible Array to append");
  // the new expanded size only changes along dim
  eshape[dim] += _shape[dim];
  // stash this original shape along dim for use in the loop:
  ind_t tshapedim = _shape[dim];
  // and resizing the Array to the new size copies existing entries
  this->resize(eshape);
  // so we need only iterate over the extra Array
  for (auto x : SubIt(extra.shape())){
    // offsetting its subscript to the new position
    auto y = x;
    y[dim] += tshapedim; // offset by the original shape along the dimension
    // and copying the contents
    (*this)[y] = extra[x];
  }
  return *this;
}

template<class T, class P>
bool Array<T,P>::all(const ind_t n) const{
  ind_t nml = this->numel();
  ind_t count = (n > 0 && n < nml) ? n : nml;
  for (ind_t i=0; i<count; ++i) if (!_data[this->l2l_d(i)]) return false;
  return true;
}
template<class T, class P>
bool Array<T,P>::any(const ind_t n) const{
  ind_t nml = this->numel();
  ind_t count = (n > 0 && n < nml) ? n : nml;
  for (ind_t i=0; i<count; ++i) if (_data[this->l2l_d(i)]) return true;
  return false;
}
template<class T, class P>
bool Array<T,P>::all(const T val, const ind_t n) const{
  ind_t nml = this->numel();
  ind_t count = (n > 0 && n < nml) ? n : nml;
  for (ind_t i=0; i<count; ++i) if (val != _data[this->l2l_d(i)]) return false;
  return true;
}
template<class T, class P>
bool Array<T,P>::any(const T val, const ind_t n) const{
  ind_t nml = this->numel();
  ind_t count = (n > 0 && n < nml) ? n : nml;
  for (ind_t i=0; i<count; ++i) if (val == _data[this->l2l_d(i)]) return true;
  return false;
}

template<class T, class P>
ind_t
Array<T,P>::count(const ind_t n) const{
  ind_t el = this->numel();
  ind_t no = (n > 0 && n < el) ? n : el;
  ind_t count{0};
  for (ind_t i=0; i<no; ++i) if (_data[this->l2l_d(i)]) ++count;
  return count;
}
template<class T, class P>
ind_t
Array<T,P>::count(const T val, const ind_t n) const{
  ind_t el = this->numel();
  ind_t no = (n > 0 && n < el) ? n : el;
  ind_t count{0};
  for (ind_t i=0; i<no; ++i) if (val == _data[this->l2l_d(i)]) ++count;
  return count;
}

template<class T, class P>
ind_t
Array<T,P>::first(const ind_t n) const{
  ind_t el = this->numel();
  ind_t no = (n > 0 && n < el) ? n : el;
  for (ind_t i=0; i<no; ++i) if (_data[this->l2l_d(i)]) return i;
  return no;
}
template<class T, class P>
ind_t
Array<T,P>::first(const T val, const ind_t n) const{
  ind_t el = this->numel();
  ind_t no = (n > 0 && n < el) ? n : el;
  for (ind_t i=0; i<no; ++i) if (val == _data[this->l2l_d(i)]) return i;
  return no;
}
template<class T, class P>
ind_t
Array<T,P>::last(const ind_t n) const{
  ind_t el = this->numel();
  ind_t no = (n > 0 && n < el) ? n : el;
  for (ind_t i=no; i--;) if (_data[this->l2l_d(i)]) return i;
  return no;
}
template<class T, class P>
ind_t
Array<T,P>::last(const T val, const ind_t n) const{
  ind_t el = this->numel();
  ind_t no = (n > 0 && n < el) ? n : el;
  for (ind_t i=no; i--;) if (val == _data[this->l2l_d(i)]) return i;
  return no;
}

#define ARRAY_ELEMENTWISE_INT_TRANSFORM(X) template<class T, class P>\
Array<int,ref_ptr_t>\
Array<T,P>:: X () const {\
  Array<int,ref_ptr_t> out(_shape, _stride);\
  for (auto x: SubIt(out.shape()))\
    out[x] = static_cast<int>(std:: X (_data[this->s2l_d(x)]));\
  return out;\
}
ARRAY_ELEMENTWISE_INT_TRANSFORM(round)
ARRAY_ELEMENTWISE_INT_TRANSFORM(floor)
ARRAY_ELEMENTWISE_INT_TRANSFORM(ceil)
#undef ARRAY_ELEMENTWISE_INT_TRANSFORM

template<class T, class P>
Array<T,ref_ptr_t>
Array<T,P>::sum(const ind_t dim) const{
  size_t ndim = _shape.size();
  assert(dim < ndim);
  shape_t ostride = this->stride();
  shape_t osize = this->shape();
  osize[dim] = 1u;
  if (this->is_row_ordered()){
    for (ind_t i=0; i<dim; ++i) ostride[i] /= _shape[dim];
  } else {
    for (ind_t i=dim+1; i<ndim; ++i) ostride[i] /= _shape[dim];
  }
  Array<T,ref_ptr_t> out(osize, ostride);
  for (auto oidx: SubIt<ind_t>(osize)) {
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

template<class T, class P>
Array<T,ref_ptr_t>
Array<T,P>::prod(const ind_t dim) const{
  size_t ndim = _shape.size();
  assert(dim < ndim);
  shape_t ostride = this->stride();
  shape_t osize = this->shape();
  osize[dim] = 1u;
  if (this->is_row_ordered()){
    for (ind_t i=0; i<dim; ++i) ostride[i] /= _shape[dim];
  } else {
    for (ind_t i=dim+1; i<ndim; ++i) ostride[i] /= _shape[dim];
  }
  Array<T,ref_ptr_t> out(osize, ostride);
  for (auto oidx: SubIt<ind_t>(osize)) {
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

template<class T, class P>
Array<T,ref_ptr_t>
Array<T,P>::min(const ind_t dim) const{
  size_t ndim = _shape.size();
  assert(dim < ndim);
  shape_t ostride = this->stride();
  shape_t osize = this->shape();
  osize[dim] = 1u;
  if (this->is_row_ordered()){
    for (ind_t i=0; i<dim; ++i) ostride[i] /= _shape[dim];
  } else {
    for (ind_t i=dim+1; i<ndim; ++i) ostride[i] /= _shape[dim];
  }
  Array<T,ref_ptr_t> out(osize, ostride);
  for (auto oidx: SubIt<ind_t>(osize)) {
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

template<class T, class P>
Array<T,ref_ptr_t>
Array<T,P>::max(const ind_t dim) const{
  size_t ndim = _shape.size();
  assert(dim < ndim);
  shape_t ostride = this->stride();
  shape_t osize = this->shape();
  osize[dim] = 1u;
  if (this->is_row_ordered()){
    for (ind_t i=0; i<dim; ++i) ostride[i] /= _shape[dim];
  } else {
    for (ind_t i=dim+1; i<ndim; ++i) ostride[i] /= _shape[dim];
  }
  // the stride along any dimension can not be zero. protect against it:
  for (ind_t i=0; i<ndim; ++i) if (ostride[i] < 1u) ostride[i] = 1u;

  Array<T,ref_ptr_t> out(osize, ostride);
  for (auto oidx: SubIt<ind_t>(osize)) {
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

template<class T, class P>
T Array<T,P>::sum() const{
  T out{0};
  for (auto x : _shape) out += _data[this->s2l_d(x)];
  return out;
}

template<class T, class P>
T Array<T,P>::prod() const{
  T out{1};
  for (auto x : _shape) out *= _data[this->s2l_d(x)];
  return out;
}

template<class T, class P>
template<class R, size_t Nel>
bool Array<T,P>::match(const ind_t i, const ind_t j, const std::array<R,Nel>& rot, const int order) const {
  assert(this->ndim() == 2u); // only defined for 2-D Arrays
  assert(Nel == _shape[1]*_shape[1]);
  ind_t n = _shape[1];
  auto ai=this->view(i);
  auto aj=this->view(j);
  std::vector<T> tmp(n,T(0));
  for (ind_t k=0; k<n; ++k) tmp[k] = aj[k];
  brille::Comparer<T,T> eq(brille::cmp::eq);
  if (order < 0){
    int o{0};
    do{
      //check against the current tmp vector whether this order rotation has moved j to i
      if (eq(n, ai.data(), tmp.data())) return true;
      // otherwise apply the rotation (again)
      brille::utils::mul_mat_vec_inplace(n, rot.data(), tmp.data());
    } while (o++ < std::abs(order));
    return false;
  } else {
    // rotate exactly order times
    for (int o=0; o<order; ++o)
      brille::utils::mul_mat_vec_inplace(n, rot.data(), tmp.data());
    return eq(n, ai.data(), tmp.data());
  }
}

template<class T, class P>
bool Array<T,P>::match(const ind_t i, const ind_t j, brille::ops op, T val) const {
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

template<class T, class P>
bool Array<T,P>::all(const brille::cmp expr, const T val) const {
  if (brille::cmp::le_ge == expr)
    return this->all(brille::cmp::le, val) || this->all(brille::cmp::ge, val);
  ind_t no = this->numel();
  brille::Comparer<T,T> op(expr);
  for (ind_t k=0; k<no; ++k) if(!op(_data[this->l2l_d(k)], val)) return false;
  return true;
}
template<class T, class P>
bool Array<T,P>::any(const brille::cmp expr, const T val) const {
  return this->first(expr, val) < this->numel();
}
template<class T, class P>
ind_t
Array<T,P>::first(const brille::cmp expr, const T val) const {
  ind_t no = this->numel();
  brille::Comparer<T,T> op(expr);
  for (ind_t k=0; k<no; ++k) if(op(_data[this->l2l_d(k)], val)) return k;
  return no;
}
template<class T, class P>
ind_t
Array<T,P>::last(const brille::cmp expr, const T val) const {
  ind_t no = this->numel();
  brille::Comparer<T,T> op(expr);
  for (ind_t k=no; k--;) if(op(_data[this->l2l_d(k)], val)) return k;
  return no;
}
template<class T, class P>
ind_t
Array<T,P>::count(const brille::cmp expr, const T val) const {
  ind_t no = this->numel();
  brille::Comparer<T,T> op(expr);
  ind_t cnt{0};
  for (ind_t k=0; k<no; ++k) if(op(_data[this->l2l_d(k)], val)) ++cnt;
  return cnt;
}

template<class T, class P>
Array<bool,ref_ptr_t>
Array<T,P>::is(const brille::cmp expr, const T val) const {
  Array<bool,ref_ptr_t> out(_shape, _stride, true);
  ind_t no = this->numel();
  brille::Comparer<T,T> op(expr);
  for (ind_t k=0; k<no; ++k) out[k] = op(_data[this->l2l_d(k)], val);
  return out;
}

template<class T, class P>
std::vector<ind_t>
Array<T,P>::find(const brille::cmp expr, const T val) const {
  Array<bool,P> this_is = this->is(expr, val);
  ind_t no = this->numel();
  std::vector<ind_t> out;
  for (ind_t k=0; k<no; ++k) if (this_is[k]) out.push_back(k);
  return out;
}

template<class T, class P> template<class R, class RP>
Array<bool,ref_ptr_t>
Array<T,P>::is(const brille::cmp expr, const Array<R,RP>& that) const {
  // To handle singleton-dimension broadcasting, this function needs to be split
  auto tsize = that.shape();
  if (!std::equal(_shape.begin(), _shape.end(), tsize.begin(),
      [](ind_t a, ind_t b){return 1==b || a==b;})){
    std::string msg = "An Array with size ( ";
    for (auto x: tsize) msg += std::to_string(x) + " ";
    msg += ") can not be broadcast to match one with size ( ";
    for (auto x: _shape) msg += std::to_string(x) + " ";
    msg += ").";
    throw std::runtime_error(msg);
  }
  Array<bool,ref_ptr_t> out(_shape, _stride, true);
  size_t ndim = _shape.size();
  if (std::equal(_shape.begin(), _shape.end(), tsize.begin())){
    // No broadcast
    brille::Comparer<T,R> op(expr);
    // no guarantees about same stride, so use subscript iterator:
    for (auto sub: SubIt<ind_t>(_shape))
      out[sub] = op(_data[this->s2l_d(sub)], that[sub]);
  } else {
    // Broadcast
    for (ind_t i=1; i<ndim; ++i) if (_shape[i]!=tsize[i])
      throw std::runtime_error("Broadcasting beyond the first dimension requires viewing beyond the first dimension!");
    for (ind_t i=0; i<_shape[0]; ++i)
      out.set(i, this->view(i).is(expr, that));
  }
  return out;
}

template<class T, class P> template<class R>
std::vector<bool>
Array<T,P>::is(const brille::cmp expr, const std::vector<R>& val) const{
  assert(2==this->ndim());
  assert(val.size() == _shape[1]);
  std::vector<bool> out;
  out.reserve(_shape[0]);
  brille::Comparer<T,R> op(expr);
  for (ind_t i=0; i<_shape[0]; ++i)
    out.push_back(op(_shape[1], this->ptr(i), _stride[1], val.data(), 1u));
  return out;
}

template<class T, class P> template<class R, class RP>
bool Array<T,P>::is(const Array<R,RP>& that) const {
  return this->is(brille::cmp::eq, that).all(true);
}

template<class T, class P>
std::vector<bool>
Array<T,P>::is_unique() const{
  assert(2==this->ndim());
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
template<class T, class P>
std::vector<ind_t>
Array<T,P>::unique_idx() const{
  assert(2==this->ndim());
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
template<class T, class P>
Array<T,ref_ptr_t>
Array<T,P>::unique() const {
  std::vector<bool> isu = this->is_unique();
  size_t u_count = std::count(isu.begin(), isu.end(), true);
  shape_t osize{_shape};
  osize[0] = static_cast<ind_t>(u_count);
  Array<T,ref_ptr_t> out(osize, _stride);
  for (ind_t i=0,u=0; i<_shape[0]; ++i) if (isu[i]) out.set(u++, this->view(i));
  return out;
}
template<class T, class P>
T Array<T,P>::dot(ind_t i, ind_t j) const {
  assert(2==this->ndim());
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

  // auto vi = ArrayIt(this->view(i));
  // auto vj = ArrayIt(this->view(j));
  // return std::transform(vi.begin() vi.end(), vj.begin(), vj.end(), T(0), std::multiples<T>, std::plus<T>);
}
template<class T, class P>
T Array<T,P>::norm(ind_t i) const {
  return std::sqrt(this->dot(i,i));
}

template<class T, class P>
Array<T,ref_ptr_t>
Array<T,P>::operator-() const{
  Array<T,ref_ptr_t> neg(_shape, _stride);
  for (auto s: SubIt(_shape)) neg[s] = -_data[this->s2l_d(s)];
  return neg;
}


template<class T, class P>
template<typename I, typename>
void Array<T,P>::permute(std::vector<I>& p){
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

template<class T, class P>
bool Array<T,P>::swap(const ind_t a, const ind_t b){
  assert(a<_shape[0] && b<_shape[0]);
  shape_t sub{_shape};
  sub[0] = a; // fix the 0th index to 'a' in the iterator
  auto itr = SubIt<ind_t>(_shape, sub);
  for (auto aidx: itr){
    // for each of the [a,...] subscripted indices, construct [b,...]
    shape_t bidx{aidx};
    bidx[0] = b;
    // precalculate the offset/stride/shaped subscripts:
    auto sa = this->s2l_d(aidx);
    auto sb = this->s2l_d(bidx);
    // perform the actual swap:
    T store = _data[sa];
    _data[sa] = _data[sb];
    _data[sb] = store;
  }
  return true;
}
template<class T, class P>
bool Array<T,P>::swap(ind_t i, ind_t a, ind_t b){
  assert(this->ndim() == 2 && i < _shape[0] && a < _shape[1] && b < _shape[1]);
  shape_t ia{i,a}, ib{i,b};
  T store = _data[this->s2l_d(ia)];
  _data[this->s2l_d(ia)] = _data[this->s2l_d(ib)];
  _data[this->s2l_d(ib)] = store;
  return true;
}

template<class T, class P>
std::vector<T> Array<T,P>::to_std() const {
  std::vector<T> out;
  for (auto x: ArrayIt<T,P>(*this)) out.push_back(x);
  return out;
}

template<class T, class P>
T* Array<T,P>::ptr(const ind_t i0){
  assert(i0 < _shape[0]);
  shape_t idx(_shape.size(), 0);
  idx[0] = i0;
  return _data + s2l_d(idx);
}
template<class T, class P>
T* Array<T,P>::ptr(const shape_t& p){
  assert(p.size() <= _shape.size());
  for (size_t i=0; i<p.size(); ++i) assert(p[i] < _shape[i]);
  shape_t idx(_shape.size(), 0);
  for (size_t i=0; i<p.size(); ++i) idx[i] = p[i];
  return _data + s2l_d(idx);
}

template<class T, class P>
const T* Array<T,P>::ptr(const ind_t i0) const {
  assert(i0 < _shape[0]);
  shape_t idx(_shape.size(), 0);
  idx[0] = i0;
  return _data + s2l_d(idx);
}
template<class T, class P>
const T* Array<T,P>::ptr(const shape_t& p) const {
  assert(p.size() <= _shape.size());
  for (size_t i=0; i<p.size(); ++i) assert(p[i] < _shape[i]);
  shape_t idx(_shape.size(), 0);
  for (size_t i=0; i<p.size(); ++i) idx[i] = p[i];
  return _data + s2l_d(idx);
}


template<class T, class P>
T& Array<T,P>::val(const ind_t i0){
  assert(i0 < _shape[0]);
  shape_t idx(_shape.size(), 0);
  idx[0] = i0;
  return _data[s2l_d(idx)];
}
template<class T, class P>
T& Array<T,P>::val(const shape_t& p){
  assert(p.size() <= _shape.size());
  for (size_t i=0; i<p.size(); ++i) assert(p[i] < _shape[i]);
  shape_t idx(_shape.size(), 0);
  for (size_t i=0; i<p.size(); ++i) idx[i] = p[i];
  return _data[s2l_d(idx)];
}
template<class T, class P> template<typename I>
T& Array<T,P>::val(std::initializer_list<I> l){
  shape_t idx;
  for (auto& x: l) idx.emplace_back(x);
  return this->val(idx);
}

template<class T, class P>
const T& Array<T,P>::val(const ind_t i0) const {
  assert(i0 < _shape[0]);
  shape_t idx(_shape.size(), 0);
  idx[0] = i0;
  return _data[s2l_d(idx)];
}
template<class T, class P>
const T& Array<T,P>::val(const shape_t& p) const {
  assert(p.size() <= _shape.size());
  for (size_t i=0; i<p.size(); ++i) assert(p[i] < _shape[i]);
  shape_t idx(_shape.size(), 0);
  for (size_t i=0; i<p.size(); ++i) idx[i] = p[i];
  return _data[s2l_d(idx)];
}
template<class T, class P> template<typename I>
const T& Array<T,P>::val(std::initializer_list<I> l) const {
  shape_t idx;
  for (auto& x: l) idx.emplace_back(x);
  return this->val(idx);
}

template<class T, class P>
Array<T,P> Array<T,P>::contiguous_copy() const {
  if (this->is_contiguous()) return Array<T,P>(*this);
  Array<T,P> out(_shape, this->calculate_stride(_shape));
  for (auto x : SubIt<ind_t>(_shape)) out[x] = (*this)[x];
  return out;
}


#define SCALAR_INPLACE_OP(X) template<class T, class P>\
Array<T,P>& Array<T,P>::operator X (const T& val){\
  if (!this->ismutable()) throw std::runtime_error("Immutable Arrays can not be modified");\
  ind_t nel=this->numel();\
  ind_t idx{0};\
  for (ind_t i=0; i<nel; ++i){\
    idx = this->l2l_d(i);\
    _data[idx] X val;\
  }\
  return *this;\
}
SCALAR_INPLACE_OP(+=)
SCALAR_INPLACE_OP(-=)
SCALAR_INPLACE_OP(*=)
SCALAR_INPLACE_OP(/=)
#undef SCALAR_INPLACE_OP

template<class T>
bool broadcast_shape_check(const std::vector<T>& a, const std::vector<T>&b){
  bool ok = a.size() == b.size();
  ok &= brille::approx::vector(a.size(), a.data(), b.data());
  if (!ok){
    std::string msg = "In place broadcasting is not possible for { ";
    for (auto x: a) msg += std::to_string(x) + " ";
    msg += "} shaped Array and { ";
    for (auto x: b) msg += std::to_string(x) + " ";
    msg += "} shaped operand";
    throw std::runtime_error(msg);
  }
  return ok;
}

#define ARRAY_INPLACE_OP(X) template<class T, class P> template<class R, class RP>\
Array<T,P>& Array<T,P>::operator X (const Array<R,RP>& val){\
  if (!this->ismutable()) throw std::runtime_error("Immutable Arrays can not be modified");\
  BroadcastIt<ind_t> itr(_shape, val.shape());\
  broadcast_shape_check(_shape, itr.shape());\
  for (auto [ox, tx, vx]: itr) _data[this->s2l_d(tx)] X val[vx];\
  return *this;\
}
ARRAY_INPLACE_OP(+=)
ARRAY_INPLACE_OP(-=)
ARRAY_INPLACE_OP(*=)
ARRAY_INPLACE_OP(/=)
#undef ARRAY_INPLACE_OP
