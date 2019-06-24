template<typename T> T* ArrayVector<T>::datapointer(size_t i, size_t j) const {
  T *ptr = nullptr;
  if (i>=this->size() || j>=this->numel()){
    std::string msg = "ArrayVector<T>::datapointer(" + std::to_string(i)
    + "," + std::to_string(j) +")" + " but size()="
    + std::to_string(this->size()) + " and numel()="
    + std::to_string(this->numel());
    throw std::domain_error(msg);
  }
  ptr = this->data + (i*this->numel() + j);
  if (!ptr) throw std::runtime_error("Attempting to access uninitialized datapointer");
  return ptr;
}

template<typename T> T ArrayVector<T>::getvalue(const size_t i, const size_t j) const {
  T *ptr, out;
  ptr = this->datapointer(i,j);
  out = *ptr;
  return out;
}
template<typename T> ArrayVector<T> ArrayVector<T>::extract(const size_t i) const {
  if (i<this->size()){
    ArrayVector<T> out(this->numel(),1u,this->datapointer(i));
    return out;
  }
  throw std::out_of_range("The requested element of the ArrayVector does not exist");
}
template<typename T> ArrayVector<T> ArrayVector<T>::extract(const size_t n, const size_t *i) const {
  bool allinbounds = true;
  ArrayVector<T> out(this->numel(),0u);
  for (size_t j=0; j<n; j++) if ( !(i[j]<this->size()) ){ allinbounds=false; break; }
  if (allinbounds){
    out.resize(n);
    for (size_t j=0; j<n; j++) out.set(j, this->datapointer(i[j]) );
  }
  return out;
}
template<typename T> ArrayVector<T> ArrayVector<T>::extract(const ArrayVector<size_t>& idx) const{
  bool allinbounds = true;
  ArrayVector<T> out(this->numel(),0u);
  if (idx.numel() != 1u) throw std::runtime_error("copying an ArrayVector by index requires ArrayVector<size_t> with numel()==1 [i.e., an ArrayScalar]");
  for (size_t j=0; j<idx.size(); ++j) if (idx.getvalue(j)>=this->size()){allinbounds=false; break;}
  if (allinbounds){
    out.resize(idx.size());
    for (size_t j=0; j<idx.size(); ++j) out.set(j, this->datapointer( idx.getvalue(j)) );
  }
  return out;
}
template<typename T> bool ArrayVector<T>::get(const size_t i, T* out) const {
  if (i>this->size()-1) return false;
  for (size_t j=0; j<this->numel(); j++) out[j]= this->getvalue(i,j);
  return true;
}
template<typename T> bool ArrayVector<T>::set(const size_t i, const T* in){
  if (i>this->size()-1) return false;
  for (size_t j=0; j<this->numel(); j++) this->data[i*this->numel()+j] = in[j];
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
template<typename T> bool ArrayVector<T>::insert(const T in, const size_t i, const size_t j){
  bool inrange = i<this->size() && j<this->numel();
  if (inrange) this->data[i*this->numel()+j] = in;
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
  else
    printf("Attempted to print elements %u to %u of size()=%u ArrayVector!\n",first,last,this->size());
}

template<typename T> void ArrayVector<T>::printheader(const char* name) const {
  printf("%s numel %u, size %u\n", name, this->numel(), this->size());
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
template<typename T> std::string ArrayVector<T>::to_string(const size_t first, const size_t last, const std::string &after) const {
  if (first<this->size() && last<this->size())
    return this->unsafe_to_string(first,last+1,after);
  std::string msg = "Attempted to print elements " + std::to_string(first)
                  + " to " + std::to_string(last) + " of size()=" + this->size()
                  + " ArrayVector!";
  throw std::domain_error(msg);
}

template<typename T> size_t ArrayVector<T>::resize(size_t newsize){
  bool std = (newsize*this->numel())>0;
  T * newdata;
  // allocate a new block of memory
  if (std) newdata = safealloc<T>(newsize*this->numel());
  if (this->size()*this->numel()) { // copy-over data :(
    size_t smallerN = (this->size() < newsize) ? this->size() : newsize;
    for (size_t i=0; i<smallerN*this->numel(); i++) newdata[i] = this->data[i];
    // hand-back the chunk of memory which data points to
    delete[] this->data;
  }
  // and set data to the newdata pointer;
  this->N = newsize;
  if (std) this->data = newdata;
  return newsize;
}
template<typename T> size_t ArrayVector<T>::refresh(size_t newnumel, size_t newsize){
  // first off, remove the old data block, if it exists
  if (this->size()*this->numel())  delete[] this->data;
  bool std = (newsize*newnumel)>0;
  T * newdata;
  // allocate a new block of memory
  if (std) newdata = safealloc<T>(newsize*newnumel);
  // and set data to the newdata pointer;
  this->M = newnumel;
  this->N = newsize;
  this->data = std ? newdata : nullptr;
  return newnumel*newsize;
}






template<typename T> template<typename R> bool ArrayVector<T>::isapprox(const ArrayVector<R> &that) const {
  AVSizeInfo si = this->consistency_check(that);
  if (si.scalara^si.scalarb) return false; // only one is an "ArrayScalar"
  for (size_t i=0; i<si.n; i++)
    for (size_t j=0; j<si.m; j++)
      if ( !approx_scalar(this->getvalue(si.oneveca?0:i,j), that.getvalue(si.onevecb?0:i,si.singular?0:j)) )
        return false;
  return true;
}
template<typename T> bool ArrayVector<T>::isapprox(const size_t i, const size_t j) const {
  for (size_t k=0; k<this->numel(); k++) if (!approx_scalar(this->getvalue(i,k),this->getvalue(j,k))) return false;
  return true;
}

template<typename T> void ArrayVector<T>::cross(const size_t i, const size_t j, T* out) const {
  if (this->numel()!=3u) throw std::domain_error("cross is only defined for 3-D vectors");
  if (i<this->size()&&j<this->size())  vector_cross(out,this->datapointer(i,0),this->datapointer(j,0));
  else throw std::domain_error("attempted to access out of bounds memory");
}
template<typename T> T ArrayVector<T>::dot(const size_t i, const size_t j) const {
  T out = 0;
  for (size_t k=0; k<this->numel(); k++) out += this->getvalue(i,k)*this->getvalue(j,k);
  return out;
}
template<typename T> T ArrayVector<T>::norm(const size_t i) const {
  return sqrt(this->dot(i,i));
}


template<typename T> bool ArrayVector<T>::arealltrue(void) const {
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0; j<this->numel(); j++)
      if (!this->getvalue(i,j)) return false;
  return true;
}
template<typename T> bool ArrayVector<T>::areanytrue(void) const {
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0; j<this->numel(); j++)
      if (this->getvalue(i,j)) return true;
  return false;
}
template<typename T> bool ArrayVector<T>::areallpositive(void) const {
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0; j<this->numel(); j++)
      if (this->getvalue(i,j)<0) return false;
  return true;
}
template<typename T> bool ArrayVector<T>::areallzero(void) const {
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0; j<this->numel(); j++)
      if (this->getvalue(i,j)) return false;
  return true;
}
template<typename T> bool ArrayVector<T>::areallapprox(const T val) const {
  T p,m, tol=2*std::numeric_limits<T>::epsilon();
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0; j<this->numel(); j++){
      m = std::abs(this->getvalue(i,j) - val);
      p = std::abs(this->getvalue(i,j) + val);
      if (m>p*tol && m>tol) return false;
    }
  return true;
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

template<class T, class R, template<class> class A,
         typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type
         >
A<S> dot(const A<T>& a, const A<R>& b){
  AVSizeInfo si = a.consistency_check(b);
  if (si.scalara^si.scalarb) throw std::runtime_error("ArrayVector dot requires equal numel()");
  A<S> out(1u,si.n);
  S tmp;
  for (size_t i=0; i<si.n; ++i){
    tmp = S(0);
    for (size_t j=0; j<si.m; ++j) tmp+= a.getvalue(si.oneveca?0:i,j) * b.getvalue(si.onevecb?0:i,j);
    out.insert(tmp,i,0);
  }
  return out;
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



// In Place arithmetic ArrayVector +-*/ ArrayVector
template<typename T> ArrayVector<T>& ArrayVector<T>:: operator +=(const ArrayVector<T> &av){
  AVSizeInfo si = this->inplace_consistency_check(av);
  for (size_t i=0; i<si.n; i++) for(size_t j=0; j<si.m; j++) this->insert( this->getvalue(i,j) + av.getvalue(si.onevecb?0:i,si.singular?0:j), i,j );
  return *this;
}
template<typename T> ArrayVector<T>& ArrayVector<T>:: operator -=(const ArrayVector<T> &av){
  AVSizeInfo si = this->inplace_consistency_check(av);
  for (size_t i=0; i<si.n; i++) for(size_t j=0; j<si.m; j++) this->insert( this->getvalue(i,j) - av.getvalue(si.onevecb?0:i,si.singular?0:j), i,j );
  return *this;
}
template<typename T> ArrayVector<T>& ArrayVector<T>:: operator *=(const ArrayVector<T> &av){
  AVSizeInfo si = this->inplace_consistency_check(av);
  for (size_t i=0; i<si.n; i++) for(size_t j=0; j<si.m; j++) this->insert( this->getvalue(i,j) * av.getvalue(si.onevecb?0:i,si.singular?0:j), i,j );
  return *this;
}
template<typename T> ArrayVector<T>& ArrayVector<T>:: operator /=(const ArrayVector<T> &av){
  AVSizeInfo si = this->inplace_consistency_check(av);
  for (size_t i=0; i<si.n; i++) for(size_t j=0; j<si.m; j++) this->insert( this->getvalue(i,j) / av.getvalue(si.onevecb?0:i,si.singular?0:j), i,j );
  return *this;
}
// In-place binary operators with scalars
template<typename T> ArrayVector<T>& ArrayVector<T>:: operator +=(const T& av){
  for (size_t i=0; i<this->size(); i++) for(size_t j=0; j<this->numel(); j++) this->insert( this->getvalue(i,j) + av, i,j );
  return *this;
}
template<typename T> ArrayVector<T>& ArrayVector<T>:: operator -=(const T& av){
  for (size_t i=0; i<this->size(); i++) for(size_t j=0; j<this->numel(); j++) this->insert( this->getvalue(i,j) - av, i,j );
  return *this;
}
template<typename T> ArrayVector<T>& ArrayVector<T>:: operator *=(const T& av){
  for (size_t i=0; i<this->size(); i++) for(size_t j=0; j<this->numel(); j++) this->insert( this->getvalue(i,j) * av, i,j );
  return *this;
}
template<typename T> ArrayVector<T>& ArrayVector<T>:: operator /=(const T& av){
  for (size_t i=0; i<this->size(); i++) for(size_t j=0; j<this->numel(); j++) this->insert( this->getvalue(i,j) / av, i,j );
  return *this;
}


template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type >
A<S> operator+(const A<T>& a, const A<R>& b){
  AVSizeInfo si = a.consistency_check(b);
  A<S> out = si.aorb ? A<S>(a) : A<S>(b);
  out.refresh(si.m,si.n); // in case a.size == b.size but one is singular, or a.numel == b.numel but one is scalar
  // if (si.oneveca || si.onevecb || si.scalara || si.scalarb){
  //   printf("=======================\n            %3s %3s %3s\n","A","B","A+B");
  //   printf("OneVector   %3d %3d\n",si.oneveca?1:0,si.onevecb?1:0);
  //   printf("ArrayScalar %3d %3d\n",si.scalara?1:0,si.scalarb?1:0);
  //   printf("-----------------------\n");
  //   printf("chosen      %3d %3d\n",si.aorb?1:0,si.aorb?0:1);
  //   printf("size()      %3u %3u %3u\n",a.size(), b.size(), out.size());
  //   printf("numel()     %3u %3u %3u\n",a.numel(), b.numel(), out.numel());
  // }
  for (size_t i=0; i<si.n; i++) for(size_t j=0; j<si.m; j++) out.insert( a.getvalue(si.oneveca?0:i,si.scalara?0:j) + b.getvalue(si.onevecb?0:i,si.scalarb?0:j), i,j );
  return out;
}
template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type >
A<S> operator-(const A<T>& a, const A<R>& b){
  AVSizeInfo si = a.consistency_check(b);
  A<S> out = si.aorb ? A<S>(a) : A<S>(b);
  out.refresh(si.m,si.n); // in case a.size == b.size but one is singular, or a.numel == b.numel but one is scalar
  for (size_t i=0; i<si.n; i++) for(size_t j=0; j<si.m; j++) out.insert( a.getvalue(si.oneveca?0:i,si.scalara?0:j) - b.getvalue(si.onevecb?0:i,si.scalarb?0:j), i,j );
  return out;
}
template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type >
A<S> operator*(const A<T>& a, const A<R>& b){
  AVSizeInfo si = a.consistency_check(b);
  A<S> out = si.aorb ? A<S>(a) : A<S>(b);
  out.refresh(si.m,si.n); // in case a.size == b.size but one is singular, or a.numel == b.numel but one is scalar
  for (size_t i=0; i<si.n; i++) for(size_t j=0; j<si.m; j++) out.insert( a.getvalue(si.oneveca?0:i,si.scalara?0:j) * b.getvalue(si.onevecb?0:i,si.scalarb?0:j), i,j );
  return out;
}
template<class T, class R, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type,
         typename=typename std::enable_if<std::is_floating_point<S>::value>::type >
A<S> operator/(const A<T>& a, const A<R>& b){
  AVSizeInfo si = a.consistency_check(b);
  A<S> out = si.aorb ? A<S>(a) : A<S>(b);
  out.refresh(si.m,si.n); // in case a.size == b.size but one is singular, or a.numel == b.numel but one is scalar
  for (size_t i=0; i<si.n; i++) for(size_t j=0; j<si.m; j++) out.insert( a.getvalue(si.oneveca?0:i,si.scalara?0:j) / b.getvalue(si.onevecb?0:i,si.scalarb?0:j), i,j );
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
  for (size_t i=0; i<(this->numel()*this->size()); i++) out.data[i] = -(this->data[i]);
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
  if (av.numel() != out.numel()) throw std::runtime_error("source and sink ArrayVectors must have same number of elements");
  if ( j >= out.size() ) throw std::out_of_range("sink index out of range");
  for (size_t k=0;k<n;++k) if (i[k]>=av.size()) throw std::out_of_range("source index out of range");
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
  S *outdata = out.datapointer(j);
  T *avidata;
  size_t m=av.numel();
  for (size_t x=0; x<n; ++x){
    avidata = av.datapointer(i[x]);
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
                    const size_t nS,
                    const size_t nE,
                    const size_t nV,
                    const size_t nM,
                    const size_t nB,
                    const size_t n,
                    const size_t *i,
                    const R *w,
                    A<S>& out,
                    const size_t j) {
  if (av.numel() != out.numel())
    throw std::runtime_error("source and sink ArrayVectors must have same number of elements");
  if ( j >= out.size() )
    throw std::out_of_range("sink index out of range");
  if (av.numel() != (nS+nE+nV+nM*nM)*nB)
    throw std::runtime_error("Wrong number of scalar/eigenvector/vector/matrix elements or branches.");
  for (size_t k=0;k<n;++k) if (i[k]>=av.size())
    throw std::out_of_range("source index out of range");
  unsafe_interpolate_to(av,nS,nE,nV,nM,nB,n,i,w,out,j);
}
/*! Combine multiple weighted arrays from one ArrayVector into a single-array ArrayVector,
    treating the elements of each vector as a series of scalars, eigenvectors,
    vectors, and matrices,
    storing the result in the specified ArrayVector at the specified index
  @param source The ArrayVector from which arrays will be extracted
  @param Nscl The number of scalar elements
  @param Neig The number of eigenvectors per array
  @param Deig The dimensionality of the eigenvectors
  @param Nvec The number of vector elements
  @param Nmat The number of matrix elements
  @param Nobj The number of branches per array
  @param Narr The number of arrays to be extraced
  @param Isrc A pointer to the indices of the arrays to be extracted
  @param weights A pointer to the weights used in combining the extracted arrays
  @param[out] sink A reference to the ArrayVector where the result will be stored
  @param Jsnk The index into sink where the array will be stored
  @note This function performs no bounds checking. Use interpolate_to if there is
        a need to ensure no sink-of-bounds access is performed.
*/

// template<class T, class R, template<class> class A,
//          typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
//          class S = typename std::common_type<T,R>::type
//          >
// void unsafe_interpolate_to(const A<T>& source,
//                            const size_t Nscl,
//                            const size_t Neig,
//                            const size_t Deig,
//                            const size_t Nvec,
//                            const size_t Nmat,
//                            const size_t Nobj,
//                            const size_t Narr,
//                            const size_t *Isrc,
//                            const R *weights,
//                            A<S>& sink,
//                            const size_t Jsnk) {
//   S *sink_j = sink.datapointer(Jsnk);
//   T *source_i, *source_0 = source.datapointer(Isrc[0]);
//   size_t offset, span = Nscl+Neig*Deig+Nvec+Nmat*Nmat;
//   T e_i_theta;
//   for (size_t x=0; x<Narr; ++x){
//     source_i = source.datapointer(Isrc[x]);
//     // loop over the modes. they are the first index and farthest in memory
//     for (size_t Iobj=0; Iobj<Nobj; ++Iobj){
//       // Scalars are first, nothing special to do:
//       for (size_t Iscl=0; Iscl<Nscl; ++Iscl)
//         sink_j[Iobj*span + Iscl] += weights[x]*source_i[Iobj*span + Iscl];
//       // Eigenvectors are next
//       for (size_t Ieig=0; Ieig<Neig; ++Ieig){
//         offset = Iobj*span + Nscl + Ieig*Deig;
//         // find the arbitrary phase eⁱᶿ between different-object eigenvectors
//         e_i_theta = antiphase(hermitian_product(Deig, source_0+offset, source_i+offset));
//         // remove the arbitrary phase as we add the weighted value
//         for(size_t Jeig=0; Jeig<Deig; ++Jeig)
//           sink_j[offset+Jeig] += weights[x]*(e_i_theta*source_i[offset+Jeig]);
//       }
//       // Vector and Matrix parts of each object are treated as scalars:
//       for (size_t Ivecmat = Nscl+Neig*Deig; Ivecmat<span; ++Ivecmat)
//         sink_j[Iobj*span + Ivecmat] += weights[x]*source_i[Iobj*span + Ivecmat];
//     }
//   }
//   // make sure each eigenvector is normalized
//   if (Neig*Deig){
//     for (size_t Iobj=0; Iobj<Nobj; ++Iobj){
//       for (size_t Ieig=0; Ieig<Neig; ++Ieig){
//         offset = Iobj*span + Nscl +Ieig*Deig;
//         auto normI = std::sqrt(inner_product(Deig, sink_j+offset, sink_j+offset));
//         for (size_t Jeig=0; Jeig<Deig; ++Jeig) sink_j[offset+Jeig] /= normI;
//       }
//     }
//   }
// }
template<class T, class R, template<class> class A,
         typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type
         >
void unsafe_interpolate_to(const A<T>& source,
                           const size_t Nscl,
                           const size_t Neig,
                           const size_t Deig,
                           const size_t Nvec,
                           const size_t Nmat,
                           const size_t Nobj,
                           const size_t Narr,
                           const size_t *Isrc,
                           const R *weights,
                           A<S>& sink,
                           const size_t Jsnk) {
  S *sink_j = sink.datapointer(Jsnk);
  T *source_i, *source_0 = source.datapointer(Isrc[0]);
  size_t offset, span = Nscl+Neig*Deig+Nvec+Nmat*Nmat;
  T e_i_theta;
  for (size_t x=0; x<Narr; ++x){
    source_i = source.datapointer(Isrc[x]);
    // loop over the modes. they are the first index and farthest in memory
    for (size_t Iobj=0; Iobj<Nobj; ++Iobj){
      // Scalars are first, nothing special to do:
      for (size_t Iscl=0; Iscl<Nscl; ++Iscl)
        sink_j[Iobj*span + Iscl] += weights[x]*source_i[Iobj*span + Iscl];
      // Eigenvectors are next
      offset = Iobj*span + Nscl;
      // find the arbitrary phase eⁱᶿ between different-object eigenvectors
      e_i_theta = antiphase(hermitian_product(Neig*Deig, source_0+offset, source_i+offset));
      // remove the arbitrary phase as we add the weighted value
      for(size_t Jeig=0; Jeig<Neig*Deig; ++Jeig)
        sink_j[offset+Jeig] += weights[x]*(e_i_theta*source_i[offset+Jeig]);
      // Vector and Matrix parts of each object are treated as scalars:
      for (size_t Ivecmat = Nscl+Neig*Deig; Ivecmat<span; ++Ivecmat)
        sink_j[Iobj*span + Ivecmat] += weights[x]*source_i[Iobj*span + Ivecmat];
    }
  }
  // make sure each eigenvector is normalized
  if (Neig*Deig){
    for (size_t Iobj=0; Iobj<Nobj; ++Iobj){
      offset = Iobj*span + Nscl;
      auto normI = std::sqrt(inner_product(Neig*Deig, sink_j+offset, sink_j+offset));
      for (size_t Jeig=0; Jeig<Neig*Deig; ++Jeig) sink_j[offset+Jeig] /= normI;
    }
  }
}
