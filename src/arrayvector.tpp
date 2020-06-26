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

// Templated implementation code for arrayvector.hpp

template<typename T> T* ArrayVector<T>::checked_data(const size_t i, const size_t j) const {
  T *ptr = nullptr;
  if (i>=this->size() || j>=this->numel()){
    std::string msg = __PRETTY_FUNCTION__;
    msg += "\n Attempting to access the pointer to element " +std::to_string(j);
    msg += " of array " + std::to_string(i) + " but the ArrayVector holds ";
    msg += std::to_string(this->size()) + " arrays each with ";
    msg += std::to_string(this->numel()) + " elements";
    throw std::out_of_range(msg);
  }
  ptr = this->data(i, j);
  if (!ptr){
    throw std::runtime_error("Attempting to access uninitialized data");
  }
  return ptr;
}
template<typename T> T* ArrayVector<T>::data(const size_t i, const size_t j) const {
  return this->_data + (i*this->numel() + j);
}

template<typename T> T ArrayVector<T>::getvalue(const size_t i, const size_t j) const {
  T *ptr, out;
  ptr = this->data(i,j);
  out = *ptr;
  return out;
}
template<typename T> ArrayVector<T> ArrayVector<T>::extract(const size_t i) const {
  if (i<this->size()){
    ArrayVector<T> out(this->numel(),1u,this->data(i));
    return out;
  }
  std::string msg = "The requested element " + std::to_string(i);
  msg += " is out of bounds for an ArrayVector with size()= ";
  msg += std::to_string(this->size());
  throw std::out_of_range(msg);
}
template<typename T> ArrayVector<T> ArrayVector<T>::first(const size_t num) const {
  size_t stop = num < this->size() ? num : this->size();
  ArrayVector<T> out(this->numel(), stop);
  for (size_t j=0; j<stop; j++) out.set(j, this->data(j) );
  return out;
}
template<typename T> ArrayVector<T> ArrayVector<T>::extract(const size_t n, const size_t *i) const {
  bool allinbounds = true;
  ArrayVector<T> out(this->numel(),0u);
  for (size_t j=0; j<n; j++) if ( !(i[j]<this->size()) ){ allinbounds=false; break; }
  if (allinbounds){
    out.resize(n);
    for (size_t j=0; j<n; j++) out.set(j, this->data(i[j]) );
  }
  return out;
}
template<typename T>
template<typename I, typename>
ArrayVector<T> ArrayVector<T>::extract(const ArrayVector<I>& idx) const {
  bool allinbounds = true;
  ArrayVector<T> out(this->numel(), 0u);
  if (idx.numel() != 1u && idx.size() != 1u){
    std::string msg = "copying an ArrayVector by index ArrayVector<Integer>";
    msg += " requires numel()==1 or size()==1";
    msg += " and the provided input has (size,numel)=(";
    msg += std::to_string(idx.size())+","+std::to_string(idx.numel())+")";
    throw std::runtime_error(msg);
  }
  for (size_t i=0; i<idx.size(); ++i) for (size_t j=0; j<idx.numel(); ++j)
  if (static_cast<size_t>(idx.getvalue(i,j)) >= this->size() || idx.getvalue(i,j)<0)
    allinbounds = false;
  if (allinbounds){
    if (idx.numel()==1u){
      out.resize(idx.size());
      for (size_t j=0; j<idx.size(); ++j) out.set(j, this->data( idx.getvalue(j)) );
    } else {
      out.resize(idx.numel());
      for (size_t j=0; j<idx.numel(); ++j) out.set(j, this->data(idx.getvalue(0,j)));
    }
  }
  return out;
}
template<typename T>
template<typename I, typename> // restricted to integer I
ArrayVector<T> ArrayVector<T>::extract(const std::vector<I>& idx) const {
  ArrayVector<T> out(this->numel(),0u);
  size_t this_size = this->size();
  if (!std::all_of(idx.begin(), idx.end(), [this_size](I j){return j>=0 || static_cast<size_t>(j)<this_size;})){
    std::string msg = "Attempting to extract out of bounds ArrayVector(s): [";
    for (auto i: idx) msg += " " + std::to_string(i);
    msg += " ] but size() = " + std::to_string(this->size());
    throw std::out_of_range(msg);
  }
  out.resize(idx.size());
  for (size_t j=0; j<idx.size(); ++j) out.set(j, this->data(idx[j]) );
  return out;
}
template<typename T>
template<typename I, typename, size_t Nel> // restricted to integer I
ArrayVector<T> ArrayVector<T>::extract(const std::array<I, Nel>& idx) const {
  ArrayVector<T> out(this->numel(),0u);
  size_t this_size = this->size();
  if (!std::all_of(idx.begin(), idx.end(), [this_size](I j){return j>=0 || static_cast<size_t>(j)<this_size;})){
    std::string msg = "Attempting to extract out of bounds ArrayVector(s): [";
    for (auto i: idx) msg += " " + std::to_string(i);
    msg += " ] but size() = " + std::to_string(this->size());
    throw std::out_of_range(msg);
  }
  out.resize(Nel);
  for (size_t j=0; j<Nel; ++j) out.set(j, this->data(idx[j]) );
  return out;
}
template<typename T> ArrayVector<T> ArrayVector<T>::extract(const ArrayVector<bool>& tf) const{
  if (tf.numel() != 1u || tf.size() != this->size()){
    std::string msg = "Extracting an ArrayVector by logical indexing requires";
    msg += " an ArrayVector<bool> with numel()==1";
    msg += " and size()==ArrayVector.size().";
    throw std::runtime_error(msg);
  }
  size_t nout=0;
  for (size_t i=0; i<tf.size(); ++i) if (tf.getvalue(i,0)) ++nout;
  ArrayVector<T> out(this->numel(),nout);
  size_t idx = 0;
  for (size_t i=0; i<tf.size(); ++i)
    if (tf.getvalue(i,0)) out.set(idx++, this->data(i));
  return out;
}
template<typename T> ArrayVector<T> ArrayVector<T>::extract(const std::vector<bool>& t) const{
  if (t.size() != this->size()){
    std::string msg = "Extracting an ArrayVector by logical indexing requires";
    msg += " a std::vector<bool> with size()==ArrayVector.size().";
    msg += " Instead got " + std::to_string(t.size()) + " where " + std::to_string(this->size()) + " was expected";
    throw std::runtime_error(msg);
  }
  ArrayVector<T> o(this->numel(), std::count(t.begin(), t.end(), true));
  size_t j = 0;
  for (size_t i=0; i<t.size(); ++i) if (t[i]) o.set(j++, this->data(i));
  return o;
}
template<typename T> bool ArrayVector<T>::get(const size_t i, T* out) const {
  if (i>this->size()-1) return false;
  for (size_t j=0; j<this->numel(); j++) out[j]= this->getvalue(i,j);
  return true;
}
template<typename T> bool ArrayVector<T>::set(const size_t i, const T* in){
  if (i>this->size()-1) return false;
  for (size_t j=0; j<this->numel(); j++) this->_data[i*this->numel()+j] = in[j];
  return true;
}
template<typename T> bool ArrayVector<T>::set(const size_t i, const ArrayVector<T>* in){
  if (i>this->size()-1 || this->numel()!=in->numel() || in->size()<1u ) return false;
  for (size_t j=0; j<this->numel(); j++) this->insert( in->getvalue(0,j), i,j) ;
  return true;
}
template<typename T> bool ArrayVector<T>::set(const size_t i, const ArrayVector<T>& in){
  if (i>this->size()-1 || this->numel()!=in.numel() || in.size()<1u ) return false;
  for (size_t j=0; j<this->numel(); j++) this->insert( in.getvalue(0,j), i,j) ;
  return true;
}
template<typename T> bool ArrayVector<T>::set(const size_t i, const std::vector<T>& in){
  if (i>this->size()-1 || this->numel()!=in.size()) return false;
  for (size_t j=0; j<this->numel(); ++j) this->insert(in[j], i, j);
  return true;
}
template<typename T>template<size_t Nel> bool ArrayVector<T>::set(const size_t i, const std::array<T,Nel>& in){
  if (i>this->size()-1 || this->numel()!=Nel) return false;
  for (size_t j=0; j<this->numel(); ++j) this->insert(in[j], i, j);
  return true;
}
template<typename T> bool ArrayVector<T>::insert(const T in, const size_t i, const size_t j){
  bool inrange = i<this->size() && j<this->numel();
  if (inrange) this->_data[i*this->numel()+j] = in;
  return inrange;
}
template<typename T> void ArrayVector<T>::printformatted(const char * fmt,const size_t first, const size_t last, const char * after) const {
  size_t i,j,b=this->numel();
  for (i=first;i<last;i++){ for (j=0;j<b;j++) printf(fmt,this->getvalue(i,j)); printf(after);  }
}
template<typename T> void ArrayVector<T>::print() const {
  const char * fmt = std::is_floating_point<T>::value ? " % g " : " % d ";
  this->printformatted(fmt,0,this->size(),"\n");
}
template<typename T> void ArrayVector<T>::print(const size_t i) const {
  this->print(i,i,"\0");
}
template<typename T> void ArrayVector<T>::print(const size_t first, const size_t last, const char *after) const {
  const char * fmt = std::is_floating_point<T>::value ? " % g " : " % d ";
  if (first<this->size() && last<this->size())
    this->printformatted(fmt,first,last+1,after);
  else {
    std::string msg = "Attempted to print elements ";
    msg += std::to_string(first) + " to " + std::to_string(last);
    msg += " of size()= " + std::to_string(this->size()) + " ArrayVector";
    throw std::out_of_range(msg);
  }
}

template<typename T> void ArrayVector<T>::printheader(const char* name) const {
  std::string header = std::string(name) + " numel "
                     + std::to_string(this->numel()) + ", size "
                     + std::to_string(this->size());
  std::cout << header << std::endl;
}

template<typename T> std::string ArrayVector<T>::unsafe_to_string(const size_t first, const size_t last, const std::string &after) const {
  size_t i,j,b=this->numel();
  std::string str;
  for (i=first;i<last;i++){
    for (j=0;j<b;j++) {
      str += my_to_string( this->getvalue(i,j) );
      // if ( str.find_last_not_of('.') ){
      //   str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
      //   str.erase ( str.find_last_not_of('.') + 1, std::string::npos );
      // }
      str += " ";
    }
    str += after;
  }
  return str;
}
template<typename T> std::string ArrayVector<T>::to_string() const {
  return this->unsafe_to_string(0,this->size(),"\n");
}
template<typename T> std::string ArrayVector<T>::to_string(const size_t i) const {
  return this->unsafe_to_string(i,i+1,"");
}
template<typename T> std::string ArrayVector<T>::to_string(const std::string &after) const {
  return this->unsafe_to_string(0,this->size(), after);
}
template<typename T> std::string ArrayVector<T>::to_string(const size_t i, const std::string &after) const {
  return this->unsafe_to_string(i,i+1, after);
}
template<typename T> std::string ArrayVector<T>::to_string(const size_t first, const size_t last, const std::string &after) const {
  if (first<this->size() && last<this->size())
    return this->unsafe_to_string(first,last+1,after);
  std::string msg = "Attempted to print elements " + std::to_string(first)
                  + " to " + std::to_string(last) + " of size("
                  + std::to_string(this->size()) + ")"
                  + " ArrayVector!";
  throw std::domain_error(msg);
}
template<typename T> template<class R> std::string ArrayVector<T>::to_string(const ArrayVector<R>& other, const size_t num) const {
  if (other.size() != this->size())
    throw std::runtime_error("ArrayVector::to_string : Equal-length ArrayVectors required.");
  size_t n = (num && num<this->size()) ? num : this->size();
  std::string s;
  for (size_t i=0; i<n; ++i){
    for (size_t j=0; j<other.numel(); ++j)
      s += my_to_string(other.getvalue(i,j)) + " ";
    s += " ";
    for (size_t j=0; j<this->numel(); ++j)
      s += my_to_string(this->getvalue(i,j)) + " ";
    s += "\n";
  }
  return s;
}

template<typename T> size_t ArrayVector<T>::resize(size_t newsize){
  bool std = (newsize*this->numel())>0;
  T * newdata = nullptr;
  // allocate a new block of memory
  if (std) newdata = new T[newsize*this->numel()]();;
  if (this->size() && this->numel()) { // copy-over _data :(
    size_t smallerN = (this->size() < newsize) ? this->size() : newsize;
    for (size_t i=0; i<smallerN*this->numel(); i++) newdata[i] = this->_data[i];
    // hand-back the chunk of memory which _data points to
    delete[] this->_data;
  }
  // and set _data to the newdata pointer;
  this->N = newsize;
  if (std) this->_data = newdata;
  return newsize;
}
template<typename T> size_t ArrayVector<T>::refresh(size_t newnumel, size_t newsize){
  // first off, remove the old _data block, if it exists
  if (this->size() && this->numel())  delete[] this->_data;
  bool std = (newsize*newnumel)>0;
  T * newdata = nullptr;
  // allocate a new block of memory
  if (std) newdata = new T[newsize*newnumel]();
  // and set _data to the newdata pointer;
  this->M = newnumel;
  this->N = newsize;
  this->_data = std ? newdata : nullptr;
  return newnumel*newsize;
}


template<typename T> template<typename R> bool ArrayVector<T>::isapprox(const ArrayVector<R> &that) const {
  AVSizeInfo si = this->consistency_check(that);
  if (si.scalara^si.scalarb) return false; // only one is an "ArrayScalar"
  bool c, u;
  T Tt;
  R Rt;
  std::tie(c,u,Tt,Rt) = determine_tols<T,T>();
  for (size_t i=0; i<si.n; i++) for (size_t j=0; j<si.m; j++)
    if ( !_approx_scalar(this->getvalue(si.oneveca?0:i,j), that.getvalue(si.onevecb?0:i,si.singular?0:j), u,Tt,Rt) )
      return false;
  return true;
}
template<typename T> bool ArrayVector<T>::isapprox(const size_t i, const size_t j) const {
  return approx_vector(this->numel(), this->data(i), this->data(j));
  // for (size_t k=0; k<this->numel(); k++) if (!approx_scalar(this->getvalue(i,k),this->getvalue(j,k))) return false;
  // return true;
}

template<typename T> void ArrayVector<T>::cross(const size_t i, const size_t j, T* out) const {
  if (this->numel()!=3u){
    throw std::domain_error("cross is only defined for 3-D vectors");
  }
  if (i<this->size()&&j<this->size())  vector_cross(out,this->data(i,0),this->data(j,0));
  else {
    std::string msg = "Attempted to access elements " + std::to_string(i);
    msg += " and " + std::to_string(j) + " of a " + std::to_string(this->size());
    msg += "-element ArrayVector";
    throw std::out_of_range(msg);
  }
}
template<typename T> T ArrayVector<T>::dot(const size_t i, const size_t j) const {
  T out = 0;
  for (size_t k=0; k<this->numel(); k++) out += this->getvalue(i,k)*this->getvalue(j,k);
  return out;
}
template<typename T> T ArrayVector<T>::norm(const size_t i) const {
  return sqrt(this->dot(i,i));
}


template<typename T> bool ArrayVector<T>::all_true(const size_t n) const {
  size_t upto = (n>0 && n <= this->size()) ? n : this->size();
  for (size_t i=0; i<upto; i++)
    for (size_t j=0; j<this->numel(); j++)
      if (!this->getvalue(i,j)) return false;
  return true;
}
template<typename T> size_t ArrayVector<T>::count_true(const size_t n) const {
  size_t upto = (n>0 && n <= this->size()) ? n : this->size();
  size_t count = 0;
  for (size_t i=0; i<upto; i++)
    for (size_t j=0; j<this->numel(); j++)
      if (this->getvalue(i,j)) ++count;
  return count;
}
template<typename T> size_t ArrayVector<T>::first_true(const size_t n) const {
  size_t upto = (n>0 && n <= this->size()) ? n : this->size();
  for (size_t i=0; i<upto; i++)
    for (size_t j=0; j<this->numel(); j++)
      if (this->getvalue(i,j)) return i;
  return upto; // an invalid indexing
}
template<typename T> size_t ArrayVector<T>::last_true(const size_t n) const {
  size_t upto = (n>0 && n <= this->size()) ? n : this->size();
  for (size_t i=upto; i--;) // evaluates for i=0, then stops
    for (size_t j=0; j<this->numel(); j++)
      if (this->getvalue(i,j)) return i;
  return upto; // an invalid indexing
}
template<typename T> bool ArrayVector<T>::any_true(const size_t n) const {
  size_t upto = (n>0 && n <= this->size()) ? n : this->size();
  for (size_t i=0; i<upto; i++)
    for (size_t j=0; j<this->numel(); j++)
      if (this->getvalue(i,j)) return true;
  return false;
}
// template<typename T> bool ArrayVector<T>::all_positive(const size_t n) const {
//   size_t upto = (n>0 && n <= this->size()) ? n : this->size();
//   for (size_t i=0; i<upto; i++)
//     for (size_t j=0; j<this->numel(); j++)
//       if (this->getvalue(i,j)<0) return false;
//   return true;
// }
template<typename T> bool ArrayVector<T>::all_zero(const size_t n) const {
  size_t upto = (n>0 && n <= this->size()) ? n : this->size();
  for (size_t i=0; i<upto; i++)
    for (size_t j=0; j<this->numel(); j++)
      if (this->getvalue(i,j)) return false;
  return true;
}
template<typename T> bool ArrayVector<T>::all_approx(const T val, const size_t n) const {
  size_t upto = (n>0 && n <= this->size()) ? n : this->size();
  bool c, u;
  T Tt, Rt;
  std::tie(c,u,Tt,Rt) = determine_tols<T,T>();
  for (size_t i=0; i<upto; i++) for (size_t j=0; j<this->numel(); j++)
    if (!_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt)) return false;
  return true;
}
template<typename T> bool ArrayVector<T>::none_approx(const T val, const size_t n) const{
  size_t upto = (n>0 && n<=this->size()) ? n : this->size();
  bool c, u;
  T Tt, Rt;
  std::tie(c,u,Tt,Rt) = determine_tols<T,T>();
  for (size_t i=0; i<upto; ++i)
  for (size_t j=0; j<this->numel(); ++j)
  if (_approx_scalar(this->getvalue(i,j), val,u,Tt,Rt)) return false;
  return true;
}
template<typename T> bool ArrayVector<T>::all_approx(const Comp expr, const T val, const size_t n) const{
  size_t upto = (n>0 && n<=this->size()) ? n : this->size();
  bool c, u;
  T Tt, Rt;
  std::tie(c,u,Tt,Rt) = determine_tols<T,T>();
  switch (expr){
    case Comp::lt:{
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt) || this->getvalue(i,j) > val) return false;
    return true;}
    case Comp::gt:{
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt) || this->getvalue(i,j) < val) return false;
    return true;}
    case Comp::le:{
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (!_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt) && this->getvalue(i,j) > val) return false;
    return true;}
    case Comp::ge:{
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (!_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt) && this->getvalue(i,j) < val) return false;
    return true;}
    case Comp::nle:{
    size_t n_approx=0, n_more=0;
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt)) ++n_approx;
    else if (this->getvalue(i,j) > val)  ++n_more;
    return (n_more > 0 || n_approx==upto);}
    case Comp::nge:{
    size_t n_approx=0, n_less=0;
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt)) ++n_approx;
    else if (this->getvalue(i,j) < val)  ++n_less;
    return (n_less > 0 || n_approx==upto);}
    case Comp::eq:{
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (!_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt)) return false;
    return true;
    case Comp::le_ge:
    bool allle=true, allge=true, ijneq;
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j){
      ijneq = !_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt);
      if (allle && ijneq && this->getvalue(i,j) > val) allle = false;
      if (allge && ijneq && this->getvalue(i,j) < val) allge = false;
      if (!(allle||allge)) return false;
    }
    return true;}
    default:
    std::string msg = __PRETTY_FUNCTION__;
    msg += ": Unknown comparator";
    throw std::runtime_error(msg);
  }
}
template<typename T> bool ArrayVector<T>::any_approx(const Comp expr, const T val, const size_t n) const{
  size_t upto = (n>0 && n<=this->size()) ? n : this->size();
  bool c, u;
  T Tt, Rt;
  std::tie(c,u,Tt,Rt) = determine_tols<T,T>();
  switch(expr){
    case Comp::lt:
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (!_approx_scalar(this->getvalue(i,j), val, u, Tt, Rt) && this->getvalue(i,j) < val) return true;
    break;
    case Comp::gt:
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (!_approx_scalar(this->getvalue(i,j), val, u, Tt, Rt) && this->getvalue(i,j) > val) return true;
    break;
    case Comp::le:
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (_approx_scalar(this->getvalue(i,j), val, u, Tt, Rt) || this->getvalue(i,j) < val) return true;
    break;
    case Comp::ge:
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (_approx_scalar(this->getvalue(i,j), val, u, Tt, Rt) || this->getvalue(i,j) < val) return true;
    break;
    case Comp::eq:
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (_approx_scalar(this->getvalue(i,j), val, u, Tt, Rt)) return true;
    break;
    default:
    std::string msg = __PRETTY_FUNCTION__;
    msg += ": Unknown comparator";
    throw std::runtime_error(msg);
  }
  return false;
}
template<typename T> ArrayVector<bool> ArrayVector<T>::is_approx(const Comp expr, const T val, const size_t n) const{
  size_t upto = (n>0 && n<=this->size()) ? n : this->size();
  ArrayVector<bool> out(1u, this->size());
  for (size_t i=0; i<this->size(); ++i) out.insert(false, i);
  bool onearray, c, u;
  T Tt, Rt;
  std::tie(c,u,Tt,Rt) = determine_tols<T,T>();
  switch (expr){
    case Comp::lt: for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt) || this->getvalue(i,j) > val)
        onearray = false;
      out.insert(onearray, i);
    } break;
    case Comp::gt: for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt) || this->getvalue(i,j) < val)
        onearray = false;
      out.insert(onearray, i);
    } break;
    case Comp::le: for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (!_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt) && this->getvalue(i,j) > val)
        onearray = false;
      out.insert(onearray, i);
    } break;
    case Comp::ge: for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (!_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt) && this->getvalue(i,j) < val)
        onearray = false;
      out.insert(onearray, i);
    } break;
    case Comp::eq: for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (!_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt))
        onearray = false;
      out.insert(onearray, i);
    } break;
    case Comp::neq: for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (_approx_scalar(this->getvalue(i,j), val, u,Tt,Rt))
        onearray = false;
      out.insert(onearray, i);
    } break;
    default:
    std::string msg = __PRETTY_FUNCTION__;
    msg += ": Unknown comparator";
    throw std::runtime_error(msg);
  }
  return out;
}
template<typename T> ArrayVector<bool> ArrayVector<T>::is_approx(const Comp expr, const std::vector<T>& vals) const{
  size_t n = vals.size();
  size_t upto = (n>0 && n<=this->size()) ? n : this->size();
  ArrayVector<bool> out(1u, this->size());
  for (size_t i=0; i<this->size(); ++i) out.insert(false, i);
  bool onearray, c, u;
  T Tt, Rt;
  std::tie(c,u,Tt,Rt) = determine_tols<T,T>();
  switch(expr){
    case Comp::lt: for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (_approx_scalar(this->getvalue(i,j), vals[i], u,Tt,Rt) || this->getvalue(i,j) > vals[i])
      onearray = false;
      out.insert(onearray, i);
    } break;
    case Comp::gt: for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (_approx_scalar(this->getvalue(i,j), vals[i], u,Tt,Rt) || this->getvalue(i,j) < vals[i])
      onearray = false;
      out.insert(onearray, i);
    } break;
    case Comp::le: for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (!_approx_scalar(this->getvalue(i,j), vals[i], u,Tt,Rt) && this->getvalue(i,j) > vals[i])
      onearray = false;
      out.insert(onearray, i);
    } break;
    case Comp::ge: for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (!_approx_scalar(this->getvalue(i,j), vals[i], u,Tt,Rt) && this->getvalue(i,j) < vals[i])
      onearray = false;
      out.insert(onearray, i);
    } break;
    case Comp::eq: for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (!_approx_scalar(this->getvalue(i,j), vals[i], u,Tt,Rt))
      onearray = false;
      out.insert(onearray, i);
    } break;
    case Comp::neq: for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (_approx_scalar(this->getvalue(i,j), vals[i], u,Tt,Rt))
      onearray = false;
      out.insert(onearray, i);
    } break;
    default:
    std::string msg = __PRETTY_FUNCTION__;
    msg += ": Unknown comparator";
    throw std::runtime_error(msg);
  }
  return out;
}

template<typename T> bool ArrayVector<T>::vector_approx(const size_t i, const size_t j, const Comp op, const T val) const{
  if (i>=this->size() || j>=this->size())
    throw std::out_of_range("ArrayVector range indices out of range");
  bool ok = true, c, u;
  T Tt, Rt;
  std::tie(c,u,Tt,Rt) = determine_tols<T,T>();
  switch (op){
    case Comp::plus:
    for (size_t k=0; k<this->numel(); ++k)
    if(!_approx_scalar(this->getvalue(i,k),this->getvalue(j,k)+val,u,Tt,Rt)) ok = false;
    return ok;
    case Comp::minus:
    for (size_t k=0; k<this->numel(); ++k)
    if(!_approx_scalar(this->getvalue(i,k),this->getvalue(j,k)-val,u,Tt,Rt)) ok = false;
    return ok;
    case Comp::times:
    for (size_t k=0; k<this->numel(); ++k)
    if(!_approx_scalar(this->getvalue(i,k),this->getvalue(j,k)*val,u,Tt,Rt)) ok = false;
    return ok;
    case Comp::rdiv:
    for (size_t k=0; k<this->numel(); ++k)
    if(!_approx_scalar(this->getvalue(i,k),this->getvalue(j,k)/val,u,Tt,Rt)) ok = false;
    return ok;
    case Comp::ldiv:
    for (size_t k=0; k<this->numel(); ++k)
    if(!_approx_scalar(this->getvalue(i,k),val/this->getvalue(j,k),u,Tt,Rt)) ok = false;
    return ok;
    case Comp::eq:
    for (size_t k=0; k<this->numel(); ++k)
    if (!_approx_scalar(this->getvalue(i,k), this->getvalue(j,k),u,Tt,Rt)) ok = false;
    return ok;
    default:
    std::string msg = __PRETTY_FUNCTION__;
    msg += ": Unknown operator";
    throw std::runtime_error(msg);
  }
}
template<typename T> template<class R, size_t Nel>
bool ArrayVector<T>::rotate_approx(const size_t i, const size_t j, const std::array<R,Nel>& mat, const int order) const{
  // if (!std::is_convertible<typename std::common_type<T,R>,T>::value)
  //   throw std::runtime_error("Incompatible types.");
  size_t n = this->numel();
  if (Nel != n*n)
    throw std::runtime_error("Wrong size matrix input.");
  if (i>=this->size() || j>=this->size())
    throw std::out_of_range("ArrayVector range indices out of range");
  std::vector<T> tmpA(n), tmpB(n);
  for (size_t k=0; k<n; ++k) tmpA[k] = this->getvalue(j,k);
  bool same = true, c, u;
  T Tt, Rt;
  std::tie(c,u,Tt,Rt) = determine_tols<T,T>();
  if (order<0){
    int o=0;
    do{
      same = true;
      // check against the current tmp vector whether this order rotations has
      // moved j to i -- if it has, same stays true and we can return true
      for (size_t k=0; k<n; ++k) if (!_approx_scalar(this->getvalue(i,k), tmpA[k],u,Tt,Rt)) same = false;
      if (same) return true;
      // otherwise we need to perform the next rotation
      mul_mat_vec(tmpB.data(), n, mat.data(), tmpA.data());
      // and assign the rotated vector back to the first temporary array
      tmpA = tmpB;
    } while (o++<std::abs(order)); // check if v[j], Rv[j], R²v[j], …, Rᵒ⁻¹v[j] ≡ v[i]
    return same;
  } else {
    // rotate exactly order times
    for (int o=0; o<order; ++o){
      mul_mat_vec(tmpB.data(), n, mat.data(), tmpA.data());
      tmpA = tmpB;
    }
    for (size_t k=0; k<n; ++k) if (!_approx_scalar(this->getvalue(i,k), tmpA[k],u,Tt,Rt)) same = false;
    return same;
  }
}

template<typename T> ArrayVector<int> ArrayVector<T>::round() const{
  ArrayVector<int> out(this->numel(),this->size());
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0; j<this->numel(); j++)
      out.insert( (int)std::round(this->getvalue(i,j)), i,j);
  return out;
}
template<typename T> ArrayVector<int> ArrayVector<T>::floor() const{
  ArrayVector<int> out(this->numel(),this->size());
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0; j<this->numel(); j++)
      out.insert( std::floor(this->getvalue(i,j)), i,j);
  return out;
}
template<typename T> ArrayVector<int> ArrayVector<T>::ceil() const{
  ArrayVector<int> out(this->numel(),this->size());
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0; j<this->numel(); j++)
      out.insert( std::ceil(this->getvalue(i,j)), i,j);
  return out;
}

template<typename T> ArrayVector<T> ArrayVector<T>::sum( const int dim ) const {
  T tmp;
  ArrayVector<T> out;
  switch (dim){
    case 1:
      out.refresh(1u,this->size());
      for (size_t i=0; i<this->size(); i++){
        tmp = T(0);
        for (size_t j=0; j<this->numel(); j++) tmp += this->getvalue(i,j);
        out.insert(tmp, i,0);
      }
      break;
    default:
      out.refresh(this->numel(),1u);
      for (size_t j=0; j<this->numel(); j++){
        tmp = T(0);
        for (size_t i=0; i<this->size(); i++) tmp += this->getvalue(i,j);
        out.insert(tmp, 0,j);
      }
      break;
  }
  return out;
}
template<typename T> ArrayVector<T> ArrayVector<T>::prod( const int dim ) const {
  T tmp;
  ArrayVector<T> out;
  switch (dim){
    case 1:
      out.refresh(1u,this->size());
      for (size_t i=0; i<this->size(); i++){
        tmp = T(1);
        for (size_t j=0; j<this->numel(); j++) tmp *= this->getvalue(i,j);
        out.insert(tmp, i,0);
      }
      break;
    default:
      out.refresh(this->numel(),1u);
      for (size_t j=0; j<this->numel(); j++){
        tmp = T(1);
        for (size_t i=0; i<this->size(); i++) tmp *= this->getvalue(i,j);
        out.insert(tmp, 0,j);
      }
      break;
  }
  return out;
}
template<typename T> ArrayVector<T> ArrayVector<T>::min( const int dim ) const {
  T tmp;
  ArrayVector<T> out;
  switch (dim){
    case 1:
      out.refresh(1u,this->size());
      for (size_t i=0; i<this->size(); i++){
        tmp = (std::numeric_limits<T>::max)();
        for (size_t j=0; j<this->numel(); j++)
        if (this->getvalue(i,j) < tmp) tmp = this->getvalue(i,j);
        out.insert(tmp, i,0);
      }
      break;
    default:
      out.refresh(this->numel(),1u);
      for (size_t j=0; j<this->numel(); j++){
        tmp = (std::numeric_limits<T>::max)();
        for (size_t i=0; i<this->size(); i++)
        if (this->getvalue(i,j) < tmp) tmp = this->getvalue(i,j);
        out.insert(tmp, 0,j);
      }
      break;
  }
  return out;
}
template<typename T> ArrayVector<T> ArrayVector<T>::max( const int dim ) const {
  T tmp;
  ArrayVector<T> out;
  switch (dim){
    case 1:
      out.refresh(1u,this->size());
      for (size_t i=0; i<this->size(); i++){
        tmp = std::numeric_limits<T>::lowest();
        for (size_t j=0; j<this->numel(); j++)
        if (this->getvalue(i,j) > tmp) tmp = this->getvalue(i,j);
        out.insert(tmp, i,0);
      }
      break;
    default:
      out.refresh(this->numel(),1u);
      for (size_t j=0; j<this->numel(); j++){
        tmp = std::numeric_limits<T>::lowest();
        for (size_t i=0; i<this->size(); i++)
        if (this->getvalue(i,j) > tmp) tmp = this->getvalue(i,j);
        out.insert(tmp, 0,j);
      }
      break;
  }
  return out;
}
template<typename T> ArrayVector<bool> ArrayVector<T>::is_unique(void) const{
  // assume all are unique to start
  ArrayVector<bool> isu(1u,this->size(),true);
  bool c, u;
  T Tt, Rt;
  std::tie(c,u,Tt,Rt) = determine_tols<T,T>();
  // and only check from the second array onwards against those of lower index
  for (size_t i=1; i<this->size(); ++i) for (size_t j=0; j<i; ++j)
  if (isu.getvalue(j) && _approx_array(this->numel(), this->data(i), this->data(j), u, Tt, Rt)){
    isu.insert(false,i);
    break;
  }
  return isu;
}
template<typename T> ArrayVector<size_t> ArrayVector<T>::unique_idx(void) const{
  ArrayVector<size_t> isu(1u,this->size());
  // assume all are unique to start
  for (size_t i=0; i<this->size(); ++i) isu.insert(i,i);
  bool c, u;
  T Tt, Rt;
  std::tie(c,u,Tt,Rt) = determine_tols<T,T>();
  // and only check from the second array onwards against those of lower index
  for (size_t i=1; i<this->size(); ++i) for (size_t j=0; j<i; ++j)
  if (j==isu.getvalue(j) && _approx_array(this->numel(), this->data(i), this->data(j), u, Tt, Rt)){
    isu.insert(j,i);
    break;
  }
  return isu;
}

/*! Extract elements of an ArrayVector while preserving its (sub)class type */
template<class T, template<class> class L,typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,L<T>>::value>::type>
L<T> extract(L<T>& source, ArrayVector<bool>& idx){
  if (idx.numel()!=1u || idx.size() != source.size()){
    std::string msg = "Extracting an ArrayVector by logical indexing requires";
    msg += " an ArrayVector<bool> with numel()==1";
    msg += " and size()==ArrayVector.size().";
    throw std::runtime_error(msg);
  }
  size_t nout=0;
  for (size_t i=0; i<idx.size(); ++i) if (idx.getvalue(i,0)) ++nout;
  L<T> sink(source);
  sink.resize(nout);
  size_t at = 0;
  for (size_t i=0; i<idx.size(); ++i)
    if (idx.getvalue(i,0)) sink.set(at++, source.data(i));
  return sink;
}
template<class T, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type>
A<T> unique(A<T>& source){
  ArrayVector<bool> uniquesorce = source.is_unique();
  return extract(source, uniquesorce);
}

// cross (ArrayVector × ArrayVector)
template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value&&!std::is_base_of<LatVec,A<T>>::value>::type,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<R>,A<R>>::value&&!std::is_base_of<LatVec,A<R>>::value>::type,
         class S=typename std::common_type<T,R>::type >
A<S> cross(const A<T>& a, const A<R>& b) {
  AVSizeInfo si = a.consistency_check(b);
  if (si.m!=3u) throw std::domain_error("cross product is only defined for three vectors");
  A<S> out( 3u, si.n);
  for (size_t i=0; i<si.n; i++)
    vector_cross<S,T,R,3>(out.data(i), a.data(si.oneveca?0:i), b.data(si.onevecb?0:i));
  return out;
}

// template<class T, class R, template<class> class A,
//          typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
//          class S = typename std::common_type<T,R>::type
//          >
// A<S> dot(const A<T>& a, const A<R>& b){
//   AVSizeInfo si = a.consistency_check(b);
//   if (si.scalara^si.scalarb)
//     throw std::runtime_error("ArrayVector dot requires equal numel()");
//   A<S> out(1u,si.n);
//   S d;
//   for (size_t i=0; i<si.n; ++i){
//     d = S(0);
//     for (size_t j=0; j<si.m; ++j)
//       d+= a.getvalue((si.oneveca ? 0 : i), j) * b.getvalue((si.onevecb ? 0 : i), j);
//     out.insert(d, i, 0);
//   }
//   return out;
// }
template<class T, class R, template<class> class A,
         typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type
         >
A<S> dot(const A<T>& a, const A<R>& b){
  AVSizeInfo si = a.consistency_check(b);
  if (si.scalara^si.scalarb)
    throw std::runtime_error("ArrayVector dot requires equal numel()");
  A<S> out(1u, si.n, S(0));
  S* d = out.data();
  if (si.oneveca && si.onevecb){
    for (size_t j=0; j<si.m; ++j) d[0] += a.getvalue(0, j) * b.getvalue(0, j);
    return out;
  }
  if (si.oneveca){
    for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) d[i] += a.getvalue(0, j) * b.getvalue(i, j);
    return out;
  }
  if (si.onevecb){
    for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) d[i] += a.getvalue(i, j) * b.getvalue(0, j);
    return out;
  }
  T* ad = a.data();
  R* bd = b.data();
  for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) d[i] += ad[i*si.m+j]*bd[i*si.m+j];
  return out;
}

template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<R>,A<R>>::value && !std::is_base_of<LatVec,A<R>>::value>::type,
         class S = typename std::common_type<T,R>::type>
A<S> cat(const A<T>& a, const A<R>& b){
  if (a.numel() != b.numel())
    throw std::runtime_error("ArrayVector cat requies equal numel()");
  A<S> out(a.numel(), a.size()+b.size());
  for (size_t i=0; i<a.size(); ++i) for (size_t j=0; j<a.numel(); ++j)
    out.insert( static_cast<S>(a.getvalue(i,j)), i, j);
  for (size_t i=0; i<b.size(); ++i) for (size_t j=0; j<b.numel(); ++j)
    out.insert( static_cast<S>(b.getvalue(i,j)), a.size()+i, j);
  return out;
}

template<class T, template<class> class A, class ...Targs,
         // typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type>
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type>
A<T> cat(const A<T>& a, const A<T>& b, const A<T>& c, Targs ...Fargs){
  return cat(cat(a,b), c, Fargs...); // recursively concatenate
}

template<class T, template<class> class L,typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,L<T>>::value>::type>
L<int> round(const L<T>& a){
  L<int> out(a);
  for (size_t i=0; i<a.size(); i++) for (size_t j=0; j<a.numel(); j++) out.insert( std::round(a.getvalue(i,j)), i,j);
  return out;
}
// floor(LatVec)
template<class T, template<class> class L,typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,L<T>>::value>::type>
L<int> floor(const L<T>& a){
  L<int> out(a);
  for (size_t i=0; i<a.size(); i++) for (size_t j=0; j<a.numel(); j++) out.insert( std::floor(a.getvalue(i,j)), i,j);
  return out;
}
// ceil(LatVec)
template<class T, template<class> class L,typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,L<T>>::value>::type>
L<int> ceil(const L<T>& a){
  L<int> out(a);
  for (size_t i=0; i<a.size(); i++) for (size_t j=0; j<a.numel(); j++) out.insert( std::ceil(a.getvalue(i,j)), i,j);
  return out;
}

template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,
         typename=typename std::enable_if<!std::is_base_of<ArrayVector<R>,R>::value>::type,
         class S = typename std::common_type<T,R>::type >
A<S> operator+(const A<T>& a, const R& b){
  A<S> out(a);
  for (size_t i=0; i<out.size(); i++) for(size_t j=0; j<out.numel(); j++) out.insert( a.getvalue(i,j) + b, i,j );
  return out;
}
template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,
         typename=typename std::enable_if<!std::is_base_of<ArrayVector<R>,R>::value>::type,
         class S = typename std::common_type<T,R>::type >
A<S> operator-(const A<T>& a, const R& b){
  A<S> out(a);
  for (size_t i=0; i<out.size(); i++) for(size_t j=0; j<out.numel(); j++) out.insert( a.getvalue(i,j) - b, i,j );
  return out;
}
template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,
         typename=typename std::enable_if<!std::is_base_of<ArrayVector<R>,R>::value>::type,
         class S = typename std::common_type<T,R>::type >
A<S> operator*(const A<T>& a, const R& b){
  A<S> out(a);
  for (size_t i=0; i<out.size(); i++) for(size_t j=0; j<out.numel(); j++) out.insert( a.getvalue(i,j) * b, i,j );
  return out;
}
// template<class T, class R, template<class> class A,
//          typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,
//          typename=typename std::enable_if<!std::is_base_of<ArrayVector<R>,R>::value>::type,
//          class S = typename std::common_type<T,R>::type,
//          typename=typename std::enable_if<std::is_floating_point<S>::value>::type >
template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,
         typename=typename std::enable_if<!std::is_base_of<ArrayVector<R>,R>::value>::type,
         class S = typename std::common_type<T,R>::type> // leave off the is_floating_point restriction on S for the special case used by halfN
A<S> operator/(const A<T>& a, const R& b){
  A<S> out(a);
  for (size_t i=0; i<out.size(); i++) for(size_t j=0; j<out.numel(); j++) out.insert( a.getvalue(i,j) / b, i,j );
  return out;
}
template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,
         typename=typename std::enable_if<!std::is_base_of<ArrayVector<R>,R>::value>::type,
         class S = typename std::common_type<T,R>::type >
A<S> operator+(const R& b, const A<T>& a){
  A<S> out(a);
  for (size_t i=0; i<out.size(); i++) for(size_t j=0; j<out.numel(); j++) out.insert( b + a.getvalue(i,j), i,j );
  return out;
}
template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,
         typename=typename std::enable_if<!std::is_base_of<ArrayVector<R>,R>::value>::type,
         class S = typename std::common_type<T,R>::type >
A<S> operator-(const R& b, const A<T>& a){
  A<S> out(a);
  for (size_t i=0; i<out.size(); i++) for(size_t j=0; j<out.numel(); j++) out.insert( b - a.getvalue(i,j), i,j );
  return out;
}
template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,
         typename=typename std::enable_if<!std::is_base_of<ArrayVector<R>,R>::value>::type,
         class S = typename std::common_type<T,R>::type >
A<S> operator*(const R& b, const A<T>& a){
  A<S> out(a);
  for (size_t i=0; i<out.size(); i++) for(size_t j=0; j<out.numel(); j++) out.insert( b * a.getvalue(i,j), i,j );
  return out;
}
template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,
         typename=typename std::enable_if<!std::is_base_of<ArrayVector<R>,R>::value>::type,
         class S = typename std::common_type<T,R>::type,
         typename=typename std::enable_if<std::is_floating_point<S>::value>::type >
A<S> operator/(const R& b, const A<T>& a){
  A<S> out(a);
  for (size_t i=0; i<out.size(); i++) for(size_t j=0; j<out.numel(); j++) out.insert( b / a.getvalue(i,j), i,j );
  return out;
}


template<typename T> ArrayVector<T> ArrayVector<T>:: operator -() const {
  ArrayVector<T> out(this->numel(),this->size());
  for (size_t i=0; i<(this->numel()*this->size()); i++) out._data[i] = -(this->_data[i]);
  return out;
}


/*! Combine multiple arrays from one ArrayVector into a single-array ArrayVector
  @param av The ArrayVector from which arrays will be extracted
  @param n The number of arrays to be extraced
  @param i A pointer to the indices of the arrays to be extracted
  @returns a single-array ArrayVector with elements that are the sum of the
           extracted arrays' elements
*/
template<typename T> ArrayVector<T> accumulate(const ArrayVector<T>& av, const size_t n, const size_t *i) {
  ArrayVector<T> out(av.numel(),1u);
  // for (size_t j=0; j<av.numel(); ++j) out.insert(T(0), 0,j);
  for (size_t j=0; j<n; j++) out += av.extract(i[j]);
  return out;
}
/*! Combine multiple weighted arrays from one ArrayVector into a single-array ArrayVector
  @param av The ArrayVector from which arrays will be extracted
  @param n The number of arrays to be extraced
  @param i A pointer to the indices of the arrays to be extracted
  @param w A pointer to the weights used in combining the extracted arrays
  @returns a single-array ArrayVector with elements that are the weighted sum
           of the extracted arrays' elements
*/
template<class T, class R, template<class> class A,
         typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type
         >
A<S> accumulate(const A<T>& av, const size_t n, const size_t *i, const R *w) {
  A<S> out(av.numel(),1u);
  // for (size_t j=0; j<av.numel(); ++j) out.insert(S(0), 0,j);
  for (size_t j=0; j<n; j++) out += av.extract(i[j]) *w[j];
  return out;
}
/*! Combine multiple weighted arrays from one ArrayVector into a single-array ArrayVector,
    storing the result in the specified ArrayVector at the specified index
  @param av The ArrayVector from which arrays will be extracted
  @param n The number of arrays to be extraced
  @param i A pointer to the indices of the arrays to be extracted
  @param w A pointer to the weights used in combining the extracted arrays
  @param[out] out A reference to the ArrayVector where the result will be stored
  @param j The index into out where the array will be stored
*/
template<class T, class R, template<class> class A,
         typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type
         >
void accumulate_to(const A<T>& av, const size_t n, const size_t *i, const R *w, A<S>& out, const size_t j) {
  if (av.numel() != out.numel()){
    throw std::runtime_error("source and sink ArrayVectors must have same number of elements");
  }
  if ( j >= out.size() ){
    throw std::out_of_range("sink index out of range");
  }
  for (size_t k=0;k<n;++k) if (i[k]>=av.size()){
    throw std::out_of_range("source index out of range");
  }
  unsafe_accumulate_to(av,n,i,w,out,j);
}
/*! Combine multiple weighted arrays from one ArrayVector into a single-array ArrayVector,
    storing the result in the specified ArrayVector at the specified index
  @param av The ArrayVector from which arrays will be extracted
  @param n The number of arrays to be extraced
  @param i A pointer to the indices of the arrays to be extracted
  @param w A pointer to the weights used in combining the extracted arrays
  @param[out] out A reference to the ArrayVector where the result will be stored
  @param j The index into out where the array will be stored
  @note This function performs no bounds checking. Use accumulate_to if there is
        a need to ensure no out-of-bounds access is performed.
*/
template<class T, class R, template<class> class A,
         typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type
         >
void unsafe_accumulate_to(const A<T>& av, const size_t n, const size_t *i, const R *w, A<S>& out, const size_t j) {
  S *outdata = out.data(j);
  size_t m=av.numel();
  for (size_t x=0; x<n; ++x){
    T *avidata = av.data(i[x]);
    for (size_t y=0; y<m; ++y)
      outdata[y] += avidata[y]*w[x];
  }
}


/*! Combine multiple weighted arrays from one ArrayVector into a single-array ArrayVector,
    treating the elements of each vector as a series of scalars, eigenvectors,
    vectors, and matrices,
    storing the result in the specified ArrayVector at the specified index
  @param av The ArrayVector from which arrays will be extracted
  @param nS The number of scalar elements
  @param nE The number of eigenvector elements
  @param nV The number of vector elements
  @param nM The number of matrix elements
  @param nB The number of branches per array
  @param n The number of arrays to be extraced
  @param i A pointer to the indices of the arrays to be extracted
  @param w A pointer to the weights used in combining the extracted arrays
  @param[out] out A reference to the ArrayVector where the result will be stored
  @param j The index into out where the array will be stored
*/
template<class T, class R, template<class> class A,
         typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type
         >
void interpolate_to(const A<T>& av,
                    const std::array<unsigned,4>& Nel,
                    const size_t nB,
                    const size_t n,
                    const size_t *i,
                    const R *w,
                    A<S>& out,
                    const size_t j) {
  if (av.numel() != out.numel()){
    throw std::runtime_error("source and sink ArrayVectors must have same number of elements");
  }
  if ( j >= out.size() ){
    throw std::out_of_range("sink index out of range");
  }
  if (av.numel() != (Nel[0]+Nel[1]+Nel[2]+Nel[3]*Nel[3])*nB){
    throw std::runtime_error("Wrong number of scalar/eigenvector/vector/matrix elements or branches.");
  }
  for (size_t k=0;k<n;++k) if (i[k]>=av.size()){
    throw std::out_of_range("source index out of range");
  }
  unsafe_interpolate_to(av,Nel,nB,n,i,w,out,j);
}
template<class T, class R, class I, template<class> class A,
    typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
    class S = typename std::common_type<T,R>::type
>
void unsafe_interpolate_to(const A<T>& source,
                           const std::array<I,3>& Nel,
                           const size_t Nobj,
                           const size_t Narr,
                           const size_t *Isrc,
                           const R *weights,   /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!!!! THIS IS WHY THERE ARE TWO VERSION OF unsafe_interpolate_to!!!!!!!!!!!!!!*/
                           A<T>& sink,
                           const size_t Jsnk) {
  info_update("The function unsafe_interpolate_to should not be used any longer. Please move to the new InterpolationData class");
  T *sink_j = sink.data(Jsnk);
  T *source_0 = source.data(Isrc[0]);
  size_t offset, span = static_cast<size_t>(Nel[0])+static_cast<size_t>(Nel[1])+static_cast<size_t>(Nel[2]);
  T e_i_theta;
  for (size_t x=0; x<Narr; ++x){
    T *source_i = source.data(Isrc[x]);
    // loop over the modes. they are the first index and farthest in memory
    for (size_t Iobj=0; Iobj<Nobj; ++Iobj){
      // Scalars are first, nothing special to do:
      for (size_t Iscl=0; Iscl<Nel[0]; ++Iscl)
        sink_j[Iobj*span + Iscl] += weights[x]*source_i[Iobj*span + Iscl];
      if (Nel[1]){
        // Eigenvectors are next
        offset = Iobj*span + Nel[0];
        // find the arbitrary phase eⁱᶿ between different-object eigenvectors
        e_i_theta = antiphase(hermitian_product(Nel[1], source_0+offset, source_i+offset));
        // remove the arbitrary phase as we add the weighted value
        for(size_t Jeig=0; Jeig<Nel[1]; ++Jeig)
          sink_j[offset + Jeig] += weights[x]*(e_i_theta*source_i[offset+Jeig]);
      }
      // Vector and Matrix parts of each object are treated as scalars:
      for (size_t Ivecmat = Nel[0]+Nel[1]; Ivecmat<span; ++Ivecmat)
        sink_j[Iobj*span + Ivecmat] += weights[x]*source_i[Iobj*span + Ivecmat];
    }
  }
  // make sure each eigenvector is normalized
  // if (Nel[1]){
  //   for (size_t Iobj=0; Iobj<Nobj; ++Iobj){
  //     offset = Iobj*span + Nel[0];
  //     auto normI = std::sqrt(inner_product(Nel[1], sink_j+offset, sink_j+offset));
  //     for (size_t Jeig=0; Jeig<Nel[1]; ++Jeig) sink_j[offset+Jeig] /= normI;
  //   }
  // }
}

template<class T, class R, class I, template<class> class A,
         typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type
         >
void new_unsafe_interpolate_to(const A<T>& source,
                           const std::array<I, 4>& nEl,
                           const size_t nObj,
                           const std::vector<size_t>& iSrc,
                           const std::vector<R>& weights,
                           A<S>& sink,
                           const size_t iSnk)
{
  info_update("The function unsafe_interpolate_to should not be used any longer. Please move to the new InterpolationData class");
  S *sink_i = sink.data(iSnk);
  T *source_0 = source.data(iSrc[0]);
  size_t offset, span = static_cast<size_t>(nEl[0])+static_cast<size_t>(nEl[1])+static_cast<size_t>(nEl[2]);
  T e_i_theta;
  for (size_t x=0; x<iSrc.size(); ++x){
    T *source_i = source.data(iSrc[x]);
    // loop over the objects (modes)
    for (size_t iObj=0; iObj<nObj; ++iObj){
      // find the weighted sum of each scalar
      for (size_t iSc=0; iSc < nEl[0]; ++iSc)
        sink_i[iObj*span + iSc] += weights[x]*source_i[iObj*span + iSc];
      if (nEl[1]){
        // eigenvectors require special treatment since they can have an arbitrary phase difference
        offset = iObj*span + nEl[0];
        // find the arbitrary phase eⁱᶿ between different-object eigenvectors
        e_i_theta = antiphase(hermitian_product(nEl[1], source_0+offset, source_i+offset));
        // remove the arbitrary phase while adding the weighted value
        for (size_t iEv=0; iEv<nEl[1]; ++iEv)
          sink_i[offset+iEv] += weights[x]*(e_i_theta*source_i[offset+iEv]);
      }
      // vector and matrix elements are treated as scalars
      for (size_t iVM=nEl[0]+nEl[1]; iVM<span; ++iVM)
        sink_i[iObj*span + iVM] += weights[x]*source_i[iObj*span + iVM];
    }
  }
  // ensure the eigenvectors are still normalized
  // if (nEl[1]){
  //   for (size_t iObj=0; iObj<nObj; ++iObj){
  //     offset = iObj*span + nEl[0];
  //     auto normI = std::sqrt(inner_product(nEl[1], sink_i+offset, sink_i+offset));
  //     for (size_t iEv=0; iEv<nEl[1]; ++iEv) sink_i[offset+iEv]/=normI;
  //   }
  // }
}

template<class T> void ArrayVector<T>::permute(const std::vector<size_t>& p){
  debug_exec(std::string msg;)
  std::vector<size_t> s=p, o(this->size());
  std::iota(o.begin(), o.end(), 0u);
  std::sort(s.begin(), s.end());
  debug_exec(\
    if (!std::includes(o.begin(), o.end(), s.begin(), s.end()) || p.size()!=this->size()){\
      msg = "The provided permutation vector [";\
      for (auto x: p) msg += " " + std::to_string(x);\
      msg += " ] is invalid";\
      msg += " A permutation of [";\
      for (auto x: o) msg += " " + std::to_string(x);\
      msg += " ] was expected.";\
    }\
  )
  // get the inverse permutation so we can swap elements
  for (size_t i=0; i<p.size(); ++i) s[p[i]] = i;
  // Now perform all swapping of Arrays until everything is in order
  ArrayVector<T> store(this->numel(), 1u);
  for (size_t i=0; i<this->size();){
    if (s[i]!=i){
      store.set(0, this->extract(i));
      this->set(i, this->extract(s[i]));
      this->set(s[i], store);
      std::swap(s[i], s[s[i]]);
    } else
      ++i;
  }
  // if debugging, confirm that the permutation worked
  debug_exec(\
    if (!std::is_sorted(s.begin(), s.end())){\
      msg = "Undoing the permutation [";\
      for (auto x: p) msg += " " + std::to_string(x);\
      msg += " ] failed. End result is [";\
      for (auto x: s) msg += " " + std::to_string(x);\
      msg += " ]";\
    }\
  )
}

template<class T> bool ArrayVector<T>::swap(const size_t i, const size_t j){
  if (i<this->size() && j<this->size()){
    ArrayVector<T> store(this->numel(), 1u);
    store.set(0, this->extract(i));
    this->set(i, this->extract(j));
    this->set(j, store);
    return true;
  }
  return false;
}
template<class T> bool ArrayVector<T>::swap(const size_t i, const size_t a, const size_t b){
  if (i<this->size() && a<this->numel() && b<this->numel()){
    T tmp = this->getvalue(i,a);
    this->insert(this->getvalue(i,b), i,a);
    this->insert(tmp, i,b);
    return true;
  }
  return false;
}

// template<class T> std::vector<size_t> ArrayVector<T>::find(void) const {
//   if (this->size() != 1u || this->numel() != 1u)
//     throw std::logic_error("Find only makes sense for ArrayScalar or ScalarVector inputs");
//   std::vector<size_t> out;
//   for (size_t i=0; i<this->size(); ++i) for (size_t j=0; j<this->numel(); ++j)
//   if (this->getvalue(i,j)) out.push_back(i+j); // either i or j is zero.
//   return out;
// }
template<class T, template<class> class A,
         typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type
         >
std::vector<size_t> find(const A<T>& a){
  if (a.size() != 1u && a.numel() != 1u)
    throw std::logic_error("Find only makes sense for ArrayScalar or ScalarVector inputs");
  std::vector<size_t> out;
  for (size_t i=0; i<a.size(); ++i) for (size_t j=0; j<a.numel(); ++j)
  if (a.getvalue(i,j)) out.push_back(i+j); // either i or j is zero.
  return out;
}

template<class T> std::vector<T> ArrayVector<T>::to_std() const {
  if (this->size() != 1u && this->numel() != 1u)
    throw std::logic_error("to_std only supports ArrayScalar or ScalarVector inputs");
  std::vector<T> out;
  for (size_t i=0; i<this->size(); ++i) for (size_t j=0; j<this->numel(); ++j)
  out.push_back(this->getvalue(i,j));
  return out;
}


template<class T> template<class I>
void
ArrayVector<T>::permute_modes(const size_t idx, std::vector<I>& indices){
  assert(idx < this->size());
  size_t n_modes = indices.size();
  size_t span = this->numel()/n_modes;
  assert(this->numel() == span*n_modes);
  // we need temporary storage
  auto holder = std::unique_ptr<T[]>(new T[span]);
  // now the actual work
  for (size_t m=0; m<n_modes; ++m) if (m != static_cast<size_t>(indices[m])){
    // move-out the current mᵗʰ mode into the temporary holder
    for(size_t j=0; j<span; ++j)
      holder[j] = this->getvalue(idx, m*span + j);
    // loop through swapping as we go
    size_t current{m};
    while (m != static_cast<size_t>(indices[current])) {
      // the next index along the chain
      size_t next = static_cast<size_t>(indices[current]);
      // copy from next to current
      for(size_t j=0; j<span; ++j)
        this->insert(this->getvalue(idx, next*span+j), idx, current*span+j);
      // setting the current destination index to the current index
      indices[current] = static_cast<I>(current);
      // update the current index to the next one
      current = next;
    }
    // finished a loop, so move the original mᵗʰ into the last place
    for(size_t j=0; j<span; ++j)
      this->insert(holder[j], idx, current*span+j);
    indices[current] = static_cast<I>(current);
  }
}

template<class T> template<class I>
void
ArrayVector<T>::inverse_permute_modes(const size_t idx, std::vector<I>& indices){
  assert(idx < this->size());
  size_t n_modes = indices.size();
  size_t span = this->numel()/n_modes;
  assert(this->numel() == span*n_modes);
  // we need temporary storage
  auto holder = std::unique_ptr<T[]>(new T[span]);
  // the inverse permutation is just a bunch of swaps.
  for (size_t m=0; m<n_modes; ++m)
  while (m != static_cast<size_t>(indices[m])) {
    size_t target = static_cast<size_t>(indices[m]);
    // move-out the targeted mode into the temporary holder
    for(size_t j=0; j<span; ++j)
      holder[j] = this->getvalue(idx, target*span + j);
    // also keep track of its target's index
    I pop_idx = indices[target];
    // move the current object and index to the target
    for(size_t j=0; j<span; ++j)
      this->insert(this->getvalue(idx, m*span+j), idx, target*span+j);
    indices[target] = target;
    // and put the temporarily held data and index back at m
    for(size_t j=0; j<span; ++j)
      this->insert(holder[j], idx, m*span+j);
    indices[m] = pop_idx;
  }
}
