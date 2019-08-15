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
template<class T> void MapGrid3<T>::check_elements(void){
  size_t total_elements = 1u;
  // scalar + eigenvector + vector + matrix*matrix elements
  size_t known_elements = static_cast<size_t>(this->elements[0])
                        + static_cast<size_t>(this->elements[1])
                        + static_cast<size_t>(this->elements[2])
                        + static_cast<size_t>(this->elements[3])*static_cast<size_t>(this->elements[3]);
  // no matter what, shape[0] should be the number of gridded points
  if (shape.size()>2){
    // if the number of dimensions of the shape array is greater than two,
    // the second element is the number of modes per point                    */
    this->branches = shape.getvalue(1u);
    for (size_t i=2u; i<this->shape.size(); ++i) total_elements *= shape.getvalue(i);
  } else {
    // shape is [n_points, n_elements] or [n_points,], so there is only one mode
    this->branches = 1u;
    total_elements = shape.size() > 1 ? shape.getvalue(1u) : 1u;
  }
  if (0 == known_elements)
    this->elements[0] = total_elements;
  if (known_elements && known_elements != total_elements){
    std::string msg ="Inconsistent element counts: "
                    + std::to_string(known_elements) + " = "
                    + std::to_string(this->elements[0]) + "+"
                    + std::to_string(this->elements[1]) + "+"
                    + std::to_string(this->elements[2]) + "+"
                    + std::to_string(this->elements[3]) + "² ≠ "
                    + std::to_string(total_elements);
    throw std::runtime_error(msg);
  }
}
template<class T>
int MapGrid3<T>::replace_data(const ArrayVector<T>& newdata,
                              const ArrayVector<size_t>& newshape,
                              const std::array<unsigned, 4>& new_elements){
  this->data =  newdata;
  this->shape = newshape;
  this->elements = new_elements;
  this->check_elements();
  return this->check_map();
}
template<class T>
int MapGrid3<T>::replace_data(const ArrayVector<T>& newdata,
                              const std::array<unsigned, 4>& new_elements){
  ArrayVector<size_t> shape(1,2);
  shape.insert(newdata.size(),0);
  shape.insert(newdata.numel(),1);
  return this->replace_data(newdata, shape, new_elements);
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
template<class T> int MapGrid3<T>::sub2map(const size_t *s, size_t& m) const {
  size_t l;
  if (this->sub2lin(s,&l)) return 1;
  if (!this->valid_mapping(l)) return -1;
  m = size_t(this->map[l]);
  return 0;
}
template<class T> size_t MapGrid3<T>::sub2map(const size_t *s) const {
  size_t m=this->maximum_mapping()+1;
  this->sub2map(s,m);
  return m;
}
// template<class T> int MapGrid3<T>::lin2map(const size_t l, size_t *m) const {
//   if ( l+1 > this->numel() ) return 1;
//   if (!this->valid_mapping(l)) return -1;
//   *m = size_t(this->map[l]);
//   return 0;
// }
template<class T> int MapGrid3<T>::lin2map(const size_t l, size_t& m) const {
  if ( l+1 > this->numel() ) return 1;
  if (!this->valid_mapping(l)) return -1;
  m = size_t(this->map[l]);
  return 0;
}
template<class T> int MapGrid3<T>::map2lin(const size_t m, size_t& l) const {
  if ( m > this->maximum_mapping() ) return 1;
  for (size_t i=0; i<this->numel(); ++i)
  if (this->map[i]>-1 && static_cast<size_t>(this->map[i])==m){
    l = i;
    return 0;
  }
  return -1;
}
//
template<class T> size_t MapGrid3<T>::numel(void) const {
  return (N==nullptr) ? 0u : N[0]*N[1]*N[2];
}
template<class T> size_t MapGrid3<T>::size(const size_t i) const {
  return (i<3u && N!=nullptr) ? N[i] : 0;
}
// template<class T> size_t MapGrid3<T>::span(const size_t i) const {
//   return (i<3u && span!=nullptr) ? span[i] : 0;
// }
template<class T> size_t MapGrid3<T>::stride(const size_t i) const {
  return sizeof(T)*((i<3u && span!=nullptr) ? span[i] : 0);
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
  // std::string msg = "index (" + std::to_string(i)
  //                       + "," + std::to_string(j)
  //                       + "," + std::to_string(k) + ")";
  // if (i >= this->size(0) ) msg += " is out of bounds along axis 0";
  // if (j >= this->size(1) ) msg += " is out of bounds along axis 1";
  // if (k >= this->size(2) ) msg += " is out of bounds along axis 2";
  // if (i<this->size(0) && j<this->size(1) && k<this->size(2)) return true;
  // msg += ": [" + std::to_string(this->size(0))
  //        + " " + std::to_string(this->size(1))
  //        + " " + std::to_string(this->size(2)) + "]";
  // msg += "/ [" + std::to_string(this->N[0])
  //        + " " + std::to_string(this->N[1])
  //        + " " + std::to_string(this->N[2]) + "]";
  // throw std::runtime_error(msg);
  return (i<this->size(0) && j<this->size(1) && k<this->size(2));
}
template<class T> bool MapGrid3<T>::is_inbounds(const size_t* s) const {
  // std::string msg = "index (" + std::to_string(s[0])
  //                       + "," + std::to_string(s[1])
  //                       + "," + std::to_string(s[2]) + ")";
  // if (s[0] >= this->size(0) ) msg += " is out of bounds along axis 0";
  // if (s[1] >= this->size(1) ) msg += " is out of bounds along axis 1";
  // if (s[2] >= this->size(2) ) msg += " is out of bounds along axis 2";
  // if (s[0]<this->size(0) && s[1]<this->size(1) && s[2]<this->size(2)) return true;
  // msg += ": [" + std::to_string(this->size(0))
  //        + " " + std::to_string(this->size(1))
  //        + " " + std::to_string(this->size(2)) + "]";
  // msg += "/ [" + std::to_string(this->N[0])
  //     + " " + std::to_string(this->N[1])
  //     + " " + std::to_string(this->N[2]) + "]";
  // throw std::runtime_error(msg);
  return (s[0]<this->size(0) && s[1]<this->size(1) && s[2]<this->size(2));
}


template<class T> ArrayVector<size_t> MapGrid3<T>::get_neighbours(const size_t centre) const {
  ArrayVector<int> mzp = make_relative_neighbour_indices_prime(1); // all combinations of [-1,0,+1] for three dimensions, skipping (0,0,0)
  ArrayVector<size_t> ijk(3u,1u);
  // get the subscripted indices of the centre position
  this->lin2sub(centre, ijk.datapointer(0));
  // flag vectors holding: is the centre index 0 (isz) or the maximum (ism)
  bool isz[3], ism[3];
  for (size_t i=0; i<3u; ++i) isz[i] = 0==ijk.getvalue(0,i);
  for (size_t i=0; i<3u; ++i) ism[i] = this->size(i)-1 <= ijk.getvalue(0,i);
  // keep track of if we *can* (or should) add each mzp vector to the centre index
  ArrayVector<bool> is_valid(1u,mzp.size());
  for (size_t i=0; i<mzp.size(); ++i){
    is_valid.insert(true,i);
    for (size_t j=0; j<mzp.numel(); ++j){
      if (isz[j] && mzp.getvalue(i,j)<0 ) is_valid.insert(false,i);
      if (ism[j] && mzp.getvalue(i,j)>0 ) is_valid.insert(false,i);
    }
  }
  // Check whether adding each mzp vector yields an in-bounds neighbour
  // and whether the remaining in-bounds neighbours contain a valid mapping
  // sub2map returns 1 if sub is out of bounds, -1 if sub is not a valid mapping
  ArrayVector<size_t> tmp(3u,1u);
  size_t _i;
  for (size_t i=0; i<mzp.size(); ++i){
    if (is_valid.getvalue(i)){
      for (size_t j=0; j<3u; ++j)
        tmp.insert( ijk.getvalue(0,j) + mzp.getvalue(i,j), 0, j);
      is_valid.insert((this->sub2map(tmp.datapointer(0), _i)) ? false : true, i);
    }
  }
  // Count how many valid neighbours are left
  size_t valid_neighbours = 0;
  for (size_t i=0; i<is_valid.size(); ++i) if (is_valid.getvalue(i)) ++valid_neighbours;
  // So we can allocate our output
  ArrayVector<size_t> neighbours(1u,valid_neighbours);
  int oob = 0;
  size_t valid_neighbour=0;
  for (size_t i=0; i<mzp.size(); ++i){
    if (is_valid.getvalue(i)){
      // we can't use
      //    tmp = mzp[i] + ijk;
      // because the compiler doesn't know what to do with ArrayVector<int> + ArrayVector<size_t>
      for (size_t j=0; j<3u; ++j) tmp.insert( ijk.getvalue(0,j) + mzp.getvalue(i,j), 0, j);
      oob += this->sub2lin(tmp.datapointer(0),neighbours.datapointer(valid_neighbour++));
    }
  }
  if (oob) throw std::runtime_error("Out-of-bounds points found when there should be none.");
  return neighbours;
};
