
template<typename T> T* ArrayVector<T>::datapointer(size_t i, size_t j) const {
  T *ptr = nullptr;
  if (i>=this->size() || j>=this->numel()){
    std::string msg = __PRETTY_FUNCTION__;
    msg += "\n Attempting to access the pointer to element " +std::to_string(j);
    msg += " of array " + std::to_string(i) + " but the ArrayVector holds ";
    msg += std::to_string(this->size()) + " arrays each with ";
    msg += std::to_string(this->numel()) + " elements";
    throw std::out_of_range(msg);
  }
  ptr = this->data + (i*this->numel() + j);
  if (!ptr){
    throw std::runtime_error("Attempting to access uninitialized datapointer");
  }
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
  std::string msg = "The requested element " + std::to_string(i);
  msg += " is out of bounds for an ArrayVector with size()= ";
  msg += std::to_string(this->size());
  throw std::out_of_range(msg);
}
template<typename T> ArrayVector<T> ArrayVector<T>::first(const size_t num) const {
  size_t stop = num < this->size() ? num : this->size();
  ArrayVector<T> out(this->numel(), stop);
  for (size_t j=0; j<stop; j++) out.set(j, this->datapointer(j) );
  return out;
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
  if (idx.numel() != 1u){
    throw std::runtime_error("copying an ArrayVector by index requires ArrayVector<size_t> with numel()==1 [i.e., an ArrayScalar]");
  }
  for (size_t j=0; j<idx.size(); ++j) if (idx.getvalue(j)>=this->size()){allinbounds=false; break;}
  if (allinbounds){
    out.resize(idx.size());
    for (size_t j=0; j<idx.size(); ++j) out.set(j, this->datapointer( idx.getvalue(j)) );
  }
  return out;
}
template<typename T> ArrayVector<T> ArrayVector<T>::extract(const std::vector<size_t>& idx) const{
  bool allinbounds = true;
  ArrayVector<T> out(this->numel(),0u);
  for (auto j: idx) if (j>=this->size()) {
    std::string msg = "Attempting to extract out of bounds ArrayVector(s): [";
    for (auto i: idx) msg += " " + std::to_string(i);
    msg += " ] but size() = " + std::to_string(this->size());
    throw std::out_of_range(msg);
  }
  if (allinbounds){
    out.resize(idx.size());
    for (size_t j=0; j<idx.size(); ++j) out.set(j, this->datapointer( idx[j]) );
  }
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
    if (tf.getvalue(i,0)) out.set(idx++, this->datapointer(i));
  return out;
}
template<typename T> ArrayVector<T> ArrayVector<T>::extract(const std::vector<bool>& t) const{
  if (t.size() != this->size()){
    std::string msg = "Extracting an ArrayVector by logical indexing requires";
    msg += " a std::vector<bool> with size()==ArrayVector.size().";
    throw std::runtime_error(msg);
  }
  size_t nout=0;
  for (auto i: t) if (i) ++nout;
  ArrayVector<T> o(this->numel(),nout);
  size_t j = 0;
  for (size_t i=0; i<t.size(); ++i) if (t[i]) o.set(j++, this->datapointer(i));
  return o;
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
  else {
    std::string msg = "Attempted to print elements ";
    msg += std::to_string(first) + " to " + std::to_string(last);
    msg += " of size()= " + std::to_string(this->size()) + " ArrayVector";
    throw std::out_of_range(msg);
  }
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
  return approx_vector(this->numel(), this->datapointer(i), this->datapointer(j));
  // for (size_t k=0; k<this->numel(); k++) if (!approx_scalar(this->getvalue(i,k),this->getvalue(j,k))) return false;
  // return true;
}

template<typename T> void ArrayVector<T>::cross(const size_t i, const size_t j, T* out) const {
  if (this->numel()!=3u){
    throw std::domain_error("cross is only defined for 3-D vectors");
  }
  if (i<this->size()&&j<this->size())  vector_cross(out,this->datapointer(i,0),this->datapointer(j,0));
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
template<typename T> bool ArrayVector<T>::any_true(const size_t n) const {
  size_t upto = (n>0 && n <= this->size()) ? n : this->size();
  for (size_t i=0; i<upto; i++)
    for (size_t j=0; j<this->numel(); j++)
      if (this->getvalue(i,j)) return true;
  return false;
}
template<typename T> bool ArrayVector<T>::all_positive(const size_t n) const {
  size_t upto = (n>0 && n <= this->size()) ? n : this->size();
  for (size_t i=0; i<upto; i++)
    for (size_t j=0; j<this->numel(); j++)
      if (this->getvalue(i,j)<0) return false;
  return true;
}
template<typename T> bool ArrayVector<T>::all_zero(const size_t n) const {
  size_t upto = (n>0 && n <= this->size()) ? n : this->size();
  for (size_t i=0; i<upto; i++)
    for (size_t j=0; j<this->numel(); j++)
      if (this->getvalue(i,j)) return false;
  return true;
}
template<typename T> bool ArrayVector<T>::all_approx(const T val, const size_t n) const {
  T p,m, tol=2*std::numeric_limits<T>::epsilon();
  size_t upto = (n>0 && n <= this->size()) ? n : this->size();
  for (size_t i=0; i<upto; i++)
    for (size_t j=0; j<this->numel(); j++){
      m = std::abs(this->getvalue(i,j) - val);
      p = std::abs(this->getvalue(i,j) + val);
      if (m>p*tol && m>tol) return false;
    }
  return true;
}
template<typename T> bool ArrayVector<T>::none_approx(const T val, const size_t n) const{
  size_t upto = (n>0 && n<=this->size()) ? n : this->size();
  for (size_t i=0; i<upto; ++i)
  for (size_t j=0; j<this->numel(); ++j)
  if (approx_scalar(this->getvalue(i,j), val)) return false;
  return true;
}
template<typename T> bool ArrayVector<T>::all_approx(const std::string& expr, const T val, const size_t n) const{
  size_t upto = (n>0 && n<=this->size()) ? n : this->size();
  if (!expr.compare("lt") || !expr.compare("<")){
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (approx_scalar(this->getvalue(i,j), val) || this->getvalue(i,j) > val) return false;
    return true;
  }
  if (!expr.compare("gt") || !expr.compare(">")){
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (approx_scalar(this->getvalue(i,j), val) || this->getvalue(i,j) < val) return false;
    return true;
  }
  if (!expr.compare("le") || !expr.compare("<=") || !expr.compare("≤")){
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (!approx_scalar(this->getvalue(i,j), val) && this->getvalue(i,j) > val) return false;
    return true;
  }
  if (!expr.compare("ge") || !expr.compare(">=") || !expr.compare("≥")){
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (!approx_scalar(this->getvalue(i,j), val) && this->getvalue(i,j) < val) return false;
    return true;
  }
  if (!expr.compare("!le") || !expr.compare("!<=") || !expr.compare("!≤")){
    size_t n_approx=0, n_more=0;
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (approx_scalar(this->getvalue(i,j), val)) ++n_approx;
    else if (this->getvalue(i,j) > val)  ++n_more;
    return (n_more > 0 || n_approx==upto);
  }
  if (!expr.compare("!ge") || !expr.compare("!>=") || !expr.compare("!≥")){
    size_t n_approx=0, n_less=0;
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (approx_scalar(this->getvalue(i,j), val)) ++n_approx;
    else if (this->getvalue(i,j) < val)  ++n_less;
    return (n_less > 0 || n_approx==upto);
  }
  if (!expr.compare("eq") || !expr.compare("==")){
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (!approx_scalar(this->getvalue(i,j), val)) return false;
    return true;
  }
  if (!expr.compare("<=|>=") || !expr.compare("≤|≥") || !expr.compare(">=|<=") || !expr.compare("≥|≤")){
    bool allle=true, allge=true, ijneq;
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j){
      ijneq = !approx_scalar(this->getvalue(i,j), val);
      if (allle && ijneq && this->getvalue(i,j) > val) allle = false;
      if (allge && ijneq && this->getvalue(i,j) < val) allge = false;
      if (!(allle||allge)) return false;
    }
    return true;
  }
  std::string msg = __PRETTY_FUNCTION__;
  msg += ": Unknown comparator " + expr;
  throw std::runtime_error(msg);
}
template<typename T> bool ArrayVector<T>::any_approx(const std::string& expr, const T val, const size_t n) const{
  T p,m, tol=2*std::numeric_limits<T>::epsilon();
  size_t upto = (n>0 && n<=this->size()) ? n : this->size();
  if (!expr.compare("lt") || !expr.compare("<")){
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (!approx_scalar(this->getvalue(i,j), val) && this->getvalue(i,j) < val) return true;
    return false;
  }
  if (!expr.compare("gt") || !expr.compare(">")){
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (!approx_scalar(this->getvalue(i,j), val) && this->getvalue(i,j) > val) return true;
    return false;
  }
  if (!expr.compare("le") || !expr.compare("<=") || !expr.compare("≤")){
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (approx_scalar(this->getvalue(i,j), val) || this->getvalue(i,j) < val) return true;
    return false;
  }
  if (!expr.compare("ge") || !expr.compare(">=") || !expr.compare("≥")){
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (approx_scalar(this->getvalue(i,j), val) || this->getvalue(i,j) < val) return true;
    return false;
  }
  if (!expr.compare("eq") || !expr.compare("==")){
    for (size_t i=0; i<upto; ++i) for (size_t j=0; j<this->numel(); ++j)
    if (approx_scalar(this->getvalue(i,j), val)) return true;
    return false;
  }
  std::string msg = __PRETTY_FUNCTION__;
  msg += ": Unknown comparator " + expr;
  throw std::runtime_error(msg);
}
template<typename T> ArrayVector<bool> ArrayVector<T>::is_approx(const std::string& expr, const T val, const size_t n) const{
  T p,m, tol=2*std::numeric_limits<T>::epsilon();
  size_t upto = (n>0 && n<=this->size()) ? n : this->size();
  ArrayVector<bool> out(1u, this->size());
  for (size_t i=0; i<this->size(); ++i) out.insert(false, i);
  bool onearray;
  if (!expr.compare("lt") || !expr.compare("<")){
    for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (approx_scalar(this->getvalue(i,j), val) || this->getvalue(i,j) > val)
      onearray = false;
      out.insert(onearray, i);
    }
    return out;
  }
  if (!expr.compare("gt") || !expr.compare(">")){
    for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (approx_scalar(this->getvalue(i,j), val) || this->getvalue(i,j) < val)
      onearray = false;
      out.insert(onearray, i);
    }
    return out;
  }
  if (!expr.compare("le") || !expr.compare("<=") || !expr.compare("≤")){
    for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (!approx_scalar(this->getvalue(i,j), val) && this->getvalue(i,j) > val)
      onearray = false;
      out.insert(onearray, i);
    }
    return out;
  }
  if (!expr.compare("ge") || !expr.compare(">=") || !expr.compare("≥")){
    for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (!approx_scalar(this->getvalue(i,j), val) && this->getvalue(i,j) < val)
      onearray = false;
      out.insert(onearray, i);
    }
    return out;
  }
  if (!expr.compare("eq") || !expr.compare("==")){
    for (size_t i=0; i<upto; ++i){
      onearray = true;
      for (size_t j=0; j<this->numel(); ++j)
      if (!approx_scalar(this->getvalue(i,j), val))
      onearray = false;
      out.insert(onearray, i);
    }
    return out;
  }
  std::string msg = __PRETTY_FUNCTION__;
  msg += ": Unknown comparator " + expr;
  throw std::runtime_error(msg);
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
template<typename T> ArrayVector<bool> ArrayVector<T>::is_unique(void) const{
  ArrayVector<bool> isu(1u,this->size());
  // assume all are unique to start
  for (size_t i=0; i<this->size(); ++i) isu.insert(true,i);
  // and only check from the second array onwards against those of lower index
  for (size_t i=1; i<this->size(); ++i) for (size_t j=0; j<i; ++j)
  if (isu.getvalue(j) && approx_vector(this->numel(), this->datapointer(i), this->datapointer(j))){
    isu.insert(false,i);
    break;
  }
  return isu;
}
template<typename T> ArrayVector<size_t> ArrayVector<T>::unique_idx(void) const{
  ArrayVector<size_t> isu(1u,this->size());
  // assume all are unique to start
  for (size_t i=0; i<this->size(); ++i) isu.insert(i,i);
  // and only check from the second array onwards against those of lower index
  for (size_t i=1; i<this->size(); ++i) for (size_t j=0; j<i; ++j)
  if (j==isu.getvalue(j) && approx_vector(this->numel(), this->datapointer(i), this->datapointer(j))){
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
    if (idx.getvalue(i,0)) sink.set(at++, source.datapointer(i));
  return sink;
}
template<class T, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type>
A<T> unique(A<T>& source){
  ArrayVector<bool> uniquesorce = source.is_unique();
  return extract(source, uniquesorce);
}

template<class T, class R, template<class> class A,
         typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type
         >
A<S> dot(const A<T>& a, const A<R>& b){
  AVSizeInfo si = a.consistency_check(b);
  if (si.scalara^si.scalarb){
    throw std::runtime_error("ArrayVector dot requires equal numel()");
  }
  A<S> out(1u,si.n);
  S tmp;
  for (size_t i=0; i<si.n; ++i){
    tmp = S(0);
    for (size_t j=0; j<si.m; ++j) tmp+= a.getvalue(si.oneveca?0:i,j) * b.getvalue(si.onevecb?0:i,j);
    out.insert(tmp,i,0);
  }
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
  #ifdef SUPER_VERBOSE
  if (si.oneveca || si.onevecb || si.scalara || si.scalarb){
    printf("=======================\n            %3s %3s %3s\n","A","B","A+B");
    printf("OneVector   %3d %3d\n",si.oneveca?1:0,si.onevecb?1:0);
    printf("ArrayScalar %3d %3d\n",si.scalara?1:0,si.scalarb?1:0);
    printf("-----------------------\n");
    printf("chosen      %3d %3d\n",si.aorb?1:0,si.aorb?0:1);
    printf("size()      %3u %3u %3u\n",a.size(), b.size(), out.size());
    printf("numel()     %3u %3u %3u\n",a.numel(), b.numel(), out.numel());
  }
  #endif
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
template<class T, class R, template<class> class A,
         typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type
         >
void unsafe_interpolate_to(const A<T>& source,
                           const std::array<unsigned,4>& Nel,
                           const size_t Nobj,
                           const size_t Narr,
                           const size_t *Isrc,
                           const R *weights,   /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!!!! THIS IS WHY THERE ARE TWO VERSION OF unsafe_interpolate_to!!!!!!!!!!!!!!*/
                           A<S>& sink,
                           const size_t Jsnk) {
  S *sink_j = sink.datapointer(Jsnk);
  T *source_i, *source_0 = source.datapointer(Isrc[0]);
  size_t offset, span = static_cast<size_t>(Nel[0])+static_cast<size_t>(Nel[1])+static_cast<size_t>(Nel[2])+static_cast<size_t>(Nel[3])*static_cast<size_t>(Nel[3]);
  T e_i_theta;
  for (size_t x=0; x<Narr; ++x){
    source_i = source.datapointer(Isrc[x]);
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
          sink_j[offset+Jeig] += weights[x]*(e_i_theta*source_i[offset+Jeig]);
      }
      // Vector and Matrix parts of each object are treated as scalars:
      for (size_t Ivecmat = Nel[0]+Nel[1]; Ivecmat<span; ++Ivecmat)
        sink_j[Iobj*span + Ivecmat] += weights[x]*source_i[Iobj*span + Ivecmat];
    }
  }
  // make sure each eigenvector is normalized
  if (Nel[1]){
    for (size_t Iobj=0; Iobj<Nobj; ++Iobj){
      offset = Iobj*span + Nel[0];
      auto normI = std::sqrt(inner_product(Nel[1], sink_j+offset, sink_j+offset));
      for (size_t Jeig=0; Jeig<Nel[1]; ++Jeig) sink_j[offset+Jeig] /= normI;
    }
  }
}

template<class T, class R, template<class> class A,
         typename=typename std::enable_if< std::is_base_of<ArrayVector<T>,A<T>>::value && !std::is_base_of<LatVec,A<T>>::value>::type,
         class S = typename std::common_type<T,R>::type
         >
void new_unsafe_interpolate_to(const A<T>& source,
                           const std::array<unsigned, 4>& nEl,
                           const size_t nObj,
                           const std::vector<size_t>& iSrc,
                           const std::vector<R>& weights,
                           A<S>& sink,
                           const size_t iSnk)
{
  S *sink_i = sink.datapointer(iSnk);
  T *source_i, *source_0 = source.datapointer(iSrc[0]);
  size_t offset, span = static_cast<size_t>(nEl[0])+static_cast<size_t>(nEl[1])+static_cast<size_t>(nEl[2])+static_cast<size_t>(nEl[3])*static_cast<size_t>(nEl[3]);
  T e_i_theta;
  for (size_t x: iSrc){
    source_i = source.datapointer(iSrc[x]);
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
  if (nEl[1]){
    for (size_t iObj=0; iObj<nObj; ++iObj){
      offset = iObj*span + nEl[0];
      auto normI = std::sqrt(inner_product(nEl[1], sink_i+offset, sink_i+offset));
      for (size_t iEv=0; iEv<nEl[1]; ++iEv) sink_i[offset+iEv]/=normI;
    }
  }
}
