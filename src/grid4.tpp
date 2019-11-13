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

// MapGrid4 Methods
template<class T> void MapGrid4<T>::print_N(const bool nl) const {
  std::cout << "[";
  for (size_t i=0; i<4u; ++i) std::cout << " " << this->N[i];
  std::cout << " ]";
  if (nl) std::cout << std::endl;
}
template<class T> void MapGrid4<T>::print_span(const bool nl) const {
  std::cout << "[";
  for (size_t i=0; i<4u; ++i) std::cout << " " << this->span[i];
  std::cout << " ]";
  if (nl) std::cout << std::endl;
}
template<class T> void MapGrid4<T>::print_map(void) const {
  for (size_t l=0; l<this->N[3]; ++l){
    for (size_t k=0; k<this->N[2]; ++k){
      for (size_t i=0; i<this->N[0]; ++i){
        for (size_t j=0; j<this->N[1]; ++j)
          std::cout << " " << this->map[sub2lin(i,j,k,l)];
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}
//
template<class T> int MapGrid4<T>::set_map(void){
  for (size_t l=0; l<this->numel(); l++) this->map[l]= (slong)(l);
  return this->numel()-1 < this->data.size() ? 0 : 1;
}
template<class T> int MapGrid4<T>::set_map(const slong* inmap, const size_t* n, const size_t d){
  if ( d!=4u ) return -3;
  if ( n[0]*n[1]*n[2]*n[3] != this->numel() ) return -2;
  if ( n[0]!=this->size(0) || n[1]!=this->size(1) || n[2]!=this->size(2) || n[3]!=this->size(3)) return -1;
  for (size_t i=0; i<this->numel(); ++i) this->map[i] = inmap[i];
  return this->check_map();
}
// only call this if inmap has been allocated *for sure* with enough space to hold the map
template<class T> int MapGrid4<T>::unsafe_set_map(slong *inmap){
  for (size_t i=0; i<this->numel(); ++i) this->map[i] = inmap[i];
  return this->check_map();
}
// only call this if outmap has enough space *for sure* to hold the map
template<class T> size_t MapGrid4<T>::unsafe_get_map(slong *outmap) const {
  size_t i=0;
  for (i=0; i<this->numel(); ++i) outmap[i] = this->map[i];
  return i;
}
//
template<class T> size_t MapGrid4<T>::maximum_mapping(const slong *map2check, const size_t num2check) const {
  size_t maxmap=0;
  for (size_t i=0; i<num2check; ++i)
    if (static_cast<size_t>(map2check[i])>maxmap)
      maxmap = static_cast<size_t>(map2check[i]);
  return maxmap;
}
template<class T> size_t MapGrid4<T>::maximum_mapping(const slong *map2check) const {
   return this->maximum_mapping(map2check, this->numel());
}
template<class T> size_t MapGrid4<T>::maximum_mapping(void) const {
  return this->maximum_mapping(this->map, this->numel());
}
//
template<class T> size_t MapGrid4<T>::valid_mapping_count(void) const {
  size_t count = 0;
  for (size_t i=0; i<this->numel(); ++i) if ( this->valid_mapping(i) ) ++count;
  return count;
}
//
template<class T> int MapGrid4<T>::check_map(const ArrayVector<T>& data2check) const {
  return ( this->maximum_mapping() < data2check.size() ) ? 0 : 1;
}
template<class T> int MapGrid4<T>::check_map(void) const {
  return this->check_map(this->data);
}
//
template<class T> int MapGrid4<T>::replace_data(const ArrayVector<T>& newdata, const ArrayVector<size_t>& newshape){
  this->data =  newdata;
  this->shape = newshape;
  return this->check_map();
}
template<class T> int MapGrid4<T>::replace_data(const ArrayVector<T>& newdata){
  ArrayVector<size_t> shape(1,2);
  shape.insert(0,newdata.size());
  shape.insert(1,newdata.numel());
  return this->replace_data(newdata, shape);
}
//
template<class T> size_t MapGrid4<T>::sub2lin(const size_t i0, const size_t i1, const size_t i2, const size_t i3) const {
  size_t l = 0;
  if (this->is_inbounds(i0,i1,i2,i3)){
    l += i0*this->span[0];
    l += i1*this->span[1];
    l += i2*this->span[2];
    l += i3*this->span[3];
  } else {
    throw std::runtime_error("accessing out of bounds subindex?");
  }
  return l;
}
template<class T> int MapGrid4<T>::sub2lin(const size_t* s, size_t *l) const {
  if (this->is_inbounds(s)){
    *l = 0;
    for (size_t i=0; i<4; i++) *l += s[i]*this->span[i];
    return 0;
  }
  return 1;
}
template<class T> int MapGrid4<T>::lin2sub(const size_t l, size_t *s) const {
  if (l < this->numel() ){
    size_t lin = l;
    for (size_t i=0; i<4u; i++){
      s[i] = lin/this->span[i];
      lin -= s[i]*this->span[i];
    }
    return 0;
  }
  return 1;
}
template<class T> int MapGrid4<T>::sub2map(const size_t* s, size_t& m) const {
  size_t l;
  if (this->sub2lin(s,&l)) return 1;
  if (!this->valid_mapping(l)) return -1;
  m = size_t(this->map[l]);
  return 0;
}
template<class T> int MapGrid4<T>::lin2map(const size_t l, size_t& m) const {
  if ( l+1 > this->numel() ) return 1;
  if (!this->valid_mapping(l)) return -1;
  m = size_t(this->map[l]);
  return 0;
}
//
template<class T> size_t MapGrid4<T>::numel(void) const {
  return (N==nullptr) ? 0u : this->N[0]*this->N[1]*this->N[2]*this->N[3];
}
template<class T> size_t MapGrid4<T>::size(const size_t i) const {
  return (i<4u && N!=nullptr) ? this->N[i] : 0;
}
// template<class T> size_t MapGrid4<T>::span(const size_t i) const {
//   return (i<4u && span!=nullptr) ? span[i] : 0;
// }
template<class T> size_t MapGrid4<T>::stride(const size_t i) const {
  return sizeof(T)*( (i<4u && span!=nullptr) ? span[i] : 0);
}
//
template<class T> size_t MapGrid4<T>::resize(const size_t n0, const size_t n1, const size_t n2, const size_t n3) {
  size_t old_numel = this->numel();
  this->N[0] = n0;
  this->N[1] = n1;
  this->N[2] = n2;
  this->N[3] = n3;
  this->calc_span();
  // don't bother deleting the memory if the resize is actually a reshape
  if ( old_numel != this->numel() && this->map!=nullptr) {
    delete[] map;
    this->map = nullptr;
  }
  this->instantiate_map();
  return this->numel();
}
template<class T> size_t MapGrid4<T>::resize(const size_t *n){
  size_t old_numel = this->numel();
  for (size_t i=0; i<4u; i++) this->N[i] = n[i];
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
template<class T> size_t MapGrid4<T>::data_ndim(void) const {
  return this->shape.size();
}
template<class T> size_t MapGrid4<T>::num_data(void) const {
  return this->data.size();
}
template<class T> ArrayVector<size_t> MapGrid4<T>::data_shape(void) const {
  return this->shape;
}
template<class T> ArrayVector<size_t> MapGrid4<T>::get_N(void) const {
  ArrayVector<size_t> out(1u,4u, this->N);
  return out;
}
// protected methods:
template<class T> void MapGrid4<T>::set_size(const size_t *n){
  for (size_t i=0; i<4u; i++) this->N[i] = n[i];
  this->calc_span();
  this->instantiate_map();
}
template<class T> void MapGrid4<T>::calc_span(){
  if (N==nullptr || this->numel()==0){
    this->span[0]=0u; this->span[1]=0u; this->span[2]=0u; this->span[3]=0u;
  } else {
    this->span[0] = this->N[3]*this->N[2]*this->N[1];
    this->span[1] = this->N[3]*this->N[2];
    this->span[2] = this->N[3];
    this->span[3] = 1;
  }
}
template<class T> void MapGrid4<T>::instantiate_map(){
  if ( this->map == nullptr && this->numel()>0 ) this->map = new slong[this->numel()]();
}
template<class T> bool MapGrid4<T>::valid_mapping(const size_t l) const {
  return this->map[l] >= 0;
}
template<class T> bool MapGrid4<T>::valid_mapping(const size_t i, const size_t j, const size_t k, const size_t l) const {
  return this->valid_mapping( this->sub2lin(i,j,k,l) );
}
template<class T> bool MapGrid4<T>::is_inbounds(const size_t i, const size_t j, const size_t k, const size_t l) const {
  return (i<this->size(0) && j<this->size(1) && k<this->size(2) && l<this->size(3));
}
template<class T> bool MapGrid4<T>::is_inbounds(const size_t* s) const {
  return (s[0]<this->size(0) && s[1]<this->size(1) && s[2]<this->size(2) && s[3]<this->size(3));
}
