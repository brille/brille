// MapGrid3 Methods

template<class T> void MapGrid3<T>::print_N(const bool nl) const {
  printf("[");
  for (size_t i=0; i<3u; i++) printf(" %3u",this->N[i]);
  printf("]");
  if (nl) printf("\n");
}
template<class T> void MapGrid3<T>::print_span(const bool nl) const {
  printf("[");
  for (size_t i=0; i<3u; i++) printf(" %3u",this->span[i]);
  printf("]");
  if (nl) printf("\n");
}
template<class T> void MapGrid3<T>::print_map(void) const {
  for (size_t k=0; k<this->N[2]; ++k){
    for (size_t i=0; i<this->N[0]; ++i){
      for (size_t j=0; j<this->N[1]; ++j){
        printf(" %4d",this->map[sub2lin(i,j,k)]);
      }
      printf("\n");
    }
    printf("\n");
  }
}
//
template<class T> int MapGrid3<T>::set_map(void){
  for (size_t l=0; l<this->numel(); l++) this->map[l]= (slong)(l);
  return this->numel()-1 < this->data.size() ? 0 : 1;
}
template<class T> int MapGrid3<T>::set_map(const slong* inmap, const size_t* n, const size_t d){
  if ( d!=3u ) return -3;
  if ( n[0]*n[1]*n[2] != this->numel() ) return -2;
  if ( n[0]!=this->size(0) || n[1]!=this->size(1) || n[2]!=this->size(2)) return -1;
  for (size_t i=0; i<this->numel(); ++i) this->map[i] = inmap[i];
  return this->check_map();
}
// only call this if inmap has been allocated *for sure* with enough space to hold the map
template<class T> int MapGrid3<T>::unsafe_set_map(slong *inmap){
  for (size_t i=0; i<this->numel(); ++i) this->map[i] = inmap[i];
  return this->check_map();
}
// only call this if outmap has enough space *for sure* to hold the map
template<class T> size_t MapGrid3<T>::unsafe_get_map(slong *outmap) const {
  size_t i=0;
  for (i=0; i<this->numel(); ++i) outmap[i] = this->map[i];
  return i;
}
//
template<class T> size_t MapGrid3<T>::maximum_mapping(const slong *map2check, const size_t num2check) const {
  slong maxmap=0;
  for (size_t i=0; i<num2check; ++i) if (map2check[i]>maxmap) maxmap = map2check[i];
  return (size_t) maxmap; // maxmap (should always be) >=0
}
template<class T> size_t MapGrid3<T>::maximum_mapping(const slong *map2check) const {
   return this->maximum_mapping(map2check, this->numel());
   // size_t maxmap = 0;
   // for (size_t i=0; i<this->numel(); ++i) if (map2check[i]>maxmap) maxmap=map2check[i];
   // return maxmap;
}
template<class T> size_t MapGrid3<T>::maximum_mapping(void) const {
  return this->maximum_mapping(this->map, this->numel());
  // size_t maxmap = 0;
  // for (size_t i=0; i<this->numel(); ++i) if (this->map[i]>maxmap) maxmap=this->map[i];
  // return maxmap;
}
//
template<class T> size_t MapGrid3<T>::valid_mapping_count(void) const {
  size_t count = 0;
  for (size_t i=0; i<this->numel(); ++i) if ( this->valid_mapping(i) ) ++count;
  return count;
}
//
template<class T> int MapGrid3<T>::check_map(const ArrayVector<T>& data2check) const {
  return ( this->maximum_mapping() < data2check.size() ) ? 0 : 1;
}
template<class T> int MapGrid3<T>::check_map(void) const {
  return this->check_map(this->data);
}
//
template<class T> int MapGrid3<T>::replace_data(const ArrayVector<T>& newdata, const ArrayVector<size_t>& newshape){
  this->data =  newdata;
  this->shape = newshape;
  return this->check_map();
}
template<class T> int MapGrid3<T>::replace_data(const ArrayVector<T>& newdata){
  ArrayVector<size_t> shape(1,2);
  shape.insert(0,newdata.size());
  shape.insert(1,newdata.numel());
  return this->replace_data(newdata, shape);
}
//
template<class T> size_t MapGrid3<T>::sub2lin(const size_t i, const size_t j, const size_t k) const {
  size_t l = 0;
  if (this->is_inbounds(i,j,k)){
    l += i*this->span[0];
    l += j*this->span[1];
    l += k*this->span[2];
  } else {
    throw std::runtime_error("accessing out of bounds subindex?");
  }
  return l;
}
template<class T> int MapGrid3<T>::sub2lin(const size_t* s, size_t *l) const {
  if (this->is_inbounds(s)){
    *l = 0;
    for (size_t i=0; i<3; i++) *l += s[i]*this->span[i];
    return 0;
  }
  return 1;
}
template<class T> int MapGrid3<T>::lin2sub(const size_t l, size_t *s) const {
  if (l < this->numel() ){
    size_t lin = l;
    for (size_t i=0; i<3u; i++){
      s[i] = lin/this->span[i];
      lin -= s[i]*this->span[i];
    }
    return 0;
  }
  return 1;
}
template<class T> int MapGrid3<T>::sub2map(const size_t *s, size_t *m) const {
  size_t l;
  if (this->sub2lin(s,&l)) return 1;
  if (!this->valid_mapping(l)) return -1;
  *m = size_t(this->map[l]);
  return 0;
}
template<class T> int MapGrid3<T>::lin2map(const size_t l, size_t *m) const {
  if ( l+1 > this->numel() ) return 1;
  if (!this->valid_mapping(l)) return -1;
  *m = size_t(this->map[l]);
  return 0;
}
//
template<class T> size_t MapGrid3<T>::numel(void) const {
  return (N==nullptr) ? 0u : N[0]*N[1]*N[2];
}
template<class T> size_t MapGrid3<T>::size(const size_t i) const {
  return (i<3u && N!=nullptr) ? N[i] : 0;
}
//
template<class T> size_t MapGrid3<T>::resize(const size_t n0, const size_t n1, const size_t n2) {
  size_t old_numel = this->numel();
  this->N[0] = n0;
  this->N[1] = n1;
  this->N[2] = n2;
  this->calc_span();
  // don't bother deleting the memory if the resize is actually a reshape
  if ( old_numel != this->numel() && this->map!=nullptr) {
    delete[] map;
    this->map = nullptr;
  }
  this->instantiate_map();
  return this->numel();
}
template<class T> size_t MapGrid3<T>::resize(const size_t *n){
  size_t old_numel = this->numel();
  for (size_t i=0; i<3u; i++) this->N[i] = n[i];
  this->calc_span();
  // don't bother deleting the memory if the resize is actually a reshape
  if ( old_numel != this->numel() && this->map!=nullptr){
    delete[] map;
    this->map = nullptr;
  }
  this->instantiate_map();
  return this->numel();
}
//
template<class T> size_t MapGrid3<T>::data_ndim(void) const {
  return this->shape.size();
}
template<class T> size_t MapGrid3<T>::num_data(void) const {
  return this->data.size();
}
template<class T> ArrayVector<size_t> MapGrid3<T>::data_shape(void) const {
  return this->shape;
}
template<class T> ArrayVector<size_t> MapGrid3<T>::get_N(void) const {
  ArrayVector<size_t> out(1u,3u, this->N);
  return out;
}
// protected methods:
template<class T> void MapGrid3<T>::set_size(const size_t *n){
  for (size_t i=0; i<3u; i++) this->N[i] = n[i];
  this->calc_span();
  this->instantiate_map();
}
template<class T> void MapGrid3<T>::calc_span(){
  if (N==nullptr || this->numel()==0){
    this->span[0]=0u; this->span[1]=0u; this->span[2]=0u;
  } else {
    this->span[0] = this->N[2]*this->N[1];
    this->span[1] = this->N[2];
    this->span[2] = 1;
  }
}
template<class T> void MapGrid3<T>::instantiate_map(){
  if ( this->map == nullptr && this->numel()>0 ) this->map = new slong[this->numel()]();
}
template<class T> bool MapGrid3<T>::valid_mapping(const size_t l) const {
  return this->map[l] >= 0;
}
template<class T> bool MapGrid3<T>::valid_mapping(const size_t i, const size_t j, const size_t k) const {
  return this->valid_mapping( this->sub2lin(i,j,k) );
}
template<class T> bool MapGrid3<T>::is_inbounds(const size_t i, const size_t j, const size_t k) const {
  std::string msg = "index (" + std::to_string(i)
                        + "," + std::to_string(j)
                        + "," + std::to_string(k) + ")";
  if (i >= this->size(0) ) msg += " is out of bounds along axis 0";
  if (j >= this->size(1) ) msg += " is out of bounds along axis 1";
  if (k >= this->size(2) ) msg += " is out of bounds along axis 2";
  if (i<this->size(0) && j<this->size(1) && k<this->size(2)) return true;
  msg += ": [" + std::to_string(this->size(0))
         + " " + std::to_string(this->size(1))
         + " " + std::to_string(this->size(2)) + "]";
  msg += "/ [" + std::to_string(this->N[0])
         + " " + std::to_string(this->N[1])
         + " " + std::to_string(this->N[2]) + "]";
  throw std::runtime_error(msg);
  return (i<this->size(0) && j<this->size(1) && k<this->size(2));
}
template<class T> bool MapGrid3<T>::is_inbounds(const size_t* s) const {
  std::string msg = "index (" + std::to_string(s[0])
                        + "," + std::to_string(s[1])
                        + "," + std::to_string(s[2]) + ")";
  if (s[0] >= this->size(0) ) msg += " is out of bounds along axis 0";
  if (s[1] >= this->size(1) ) msg += " is out of bounds along axis 1";
  if (s[2] >= this->size(2) ) msg += " is out of bounds along axis 2";
  if (s[0]<this->size(0) && s[1]<this->size(1) && s[2]<this->size(2)) return true;
  msg += ": [" + std::to_string(this->size(0))
         + " " + std::to_string(this->size(1))
         + " " + std::to_string(this->size(2)) + "]";
  msg += "/ [" + std::to_string(this->N[0])
      + " " + std::to_string(this->N[1])
      + " " + std::to_string(this->N[2]) + "]";
  throw std::runtime_error(msg);
  return (s[0]<this->size(0) && s[1]<this->size(1) && s[2]<this->size(2));
}
