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

template<typename T> LQVec<T> LQVec<T>::extract(const size_t i) const {
  if (i<this->size()){
    LQVec<T> out(this->get_lattice(),1u,this->data(i));
    return out;
  }
  std::string msg = "The requested element " + std::to_string(i);
  msg += " is out of bounds for a LQVec with size()= ";
  msg += std::to_string(this->size());
  throw std::out_of_range(msg);
}
template<typename T> LQVec<T> LQVec<T>::first(const size_t num) const {
  size_t stop = num < this->size() ? num : this->size();
  LQVec<T> out(this->get_lattice(), stop);
  for (size_t j=0; j<stop; j++) out.set(j, this->data(j) );
  return out;
}
template<typename T> LQVec<T> LQVec<T>::extract(const size_t n, const size_t *i) const {
  bool allinbounds = true;
  LQVec<T> out(this->get_lattice(),0u);
  for (size_t j=0; j<n; j++) if ( !(i[j]<this->size()) ){ allinbounds=false; break; }
  if (allinbounds){
    out.resize(n);
    for (size_t j=0; j<n; j++) out.set(j, this->data(i[j]) );
  }
  return out;
}
template<typename T> LQVec<T> LQVec<T>::extract(const ArrayVector<size_t>& idx) const{
  bool allinbounds = true;
  LQVec<T> out(this->get_lattice(),0u);
  if (idx.numel() != 1u) throw std::runtime_error("copying an ArrayVector by index requires ArrayVector<size_t> with numel()==1 [i.e., an ArrayScalar]");
  for (size_t j=0; j<idx.size(); ++j) if (idx.getvalue(j)>=this->size()){allinbounds=false; break;}
  if (allinbounds){
    out.resize(idx.size());
    for (size_t j=0; j<idx.size(); ++j) out.set(j, this->data( idx.getvalue(j)) );
  }
  return out;
}
template<typename T> LQVec<T> LQVec<T>::extract(const ArrayVector<bool>& tf) const{
  if (tf.numel() != 1u || tf.size() != this->size()){
    std::string msg = "Extracting an LQVec by logical indexing requires";
    msg += " an ArrayVector<bool> with numel()==1";
    msg += " and size()==LQVec.size().";
    throw std::runtime_error(msg);
  }
  size_t nout=0;
  for (size_t i=0; i<tf.size(); ++i) if (tf.getvalue(i,0)) ++nout;
  LQVec<T> out(this->get_lattice(),nout);
  size_t idx = 0;
  for (size_t i=0; i<tf.size(); ++i)
    if (tf.getvalue(i,0)) out.set(idx++, this->data(i));
  return out;
}
template<typename T> LQVec<T> LQVec<T>::extract(const std::vector<bool>& tf) const{
  if (tf.size() != this->size()){
    std::string msg = "Extracting an LQVec by logical indexing requires";
    msg += " an std::vector<bool> with size()==LQVec.size().";
    throw std::runtime_error(msg);
  }
  size_t nout=0;
  for (size_t i=0; i<tf.size(); ++i) if (tf[i]) ++nout;
  LQVec<T> out(this->get_lattice(),nout);
  size_t idx = 0;
  for (size_t i=0; i<tf.size(); ++i)
    if (tf[i]) out.set(idx++, this->data(i));
  return out;
}

template<typename T> LQVec<T> LQVec<T>::get(const size_t i) const {
  LQVec<T> out(this->get_lattice(), i<this->size() ? 1u : 0u );
  if (i<this->size()) this->ArrayVector<T>::get(i, out.data() );
  return out;
}

template<typename T> ArrayVector<T> LQVec<T>::get_hkl() const {
  return ArrayVector<T>(this->numel(),this->size(),this->_data); // strip off the Lattice information
}
template<typename T> ArrayVector<double> LQVec<T>::get_xyz() const {
  double toxyz[9];
  Reciprocal lat = this->get_lattice();
  lat.get_xyz_transform(toxyz);
  ArrayVector<double> xyz(this->numel(),this->size());
  for (size_t i=0; i<this->size(); i++) multiply_matrix_vector(xyz.data(i), toxyz, this->data(i));
  return xyz;
}
template<typename T> LDVec<double> LQVec<T>::star() const {
  double cvmt[9];
  (this->lattice).get_covariant_metric_tensor(cvmt);
  LDVec<double> slv( (this->lattice).star(), this->size() );
  for (size_t i=0; i<this->size(); i++) multiply_matrix_vector(slv.data(i), cvmt, this->data(i));
  slv /= 2.0*PI; // ai= gij/2/pi * ai_star
  return slv;
}
// template<typename T> template<typename R> LQVec<double> LQVec<T>::cross(const LQVec<R> * vec) const {
//   assert( this->samelattice(vec) );
//   AVSizeInfo si = av_consistency_check(this,vec);
//   if (si.m!=3u) throw "cross product is only defined for three vectors";
//   LDVec<double> tmp((vec->get_lattice()).star(), si.n);
//   for (size_t i=0; i<si.n; i++) vector_cross<double,T,R,3>(tmp.data(i), this->data(si.a?0:i), vec->data(si.b?0:i));
//     tmp *= (vec->get_lattice()).get_volume()/2.0/PI;
//     return tmp.star();
// }
template<typename T> LQVec<double> LQVec<T>::cross(const size_t i, const size_t j) const {
  bool bothok = (i<this->size() && j<this->size() && 3u==this->numel());
  LQVec<double> out(this->get_lattice(), bothok? 1u : 0u);
  if (bothok){
    double rlucross[3];
    vector_cross<double,T,T,3>(rlucross, this->data(i), this->data(j));

    Reciprocal rlat = this->get_lattice();
    LDVec<double> ldv( rlat.star(), 1u, rlucross);
    ldv *= rlat.get_volume()/2.0/PI;
    out =  ldv.star();
  }
  return out;
}

template<typename T> double LQVec<T>::dot(const size_t i, const size_t j) const {
  if (i>=this->size() || j>=this->size())
    throw std::out_of_range("attempted out of bounds access by dot");
  Reciprocal lat = this->get_lattice();
  double len[3] = {lat.get_a(),lat.get_b(),lat.get_c()};
  double ang[3] = {lat.get_alpha(),lat.get_beta(),lat.get_gamma()};
  return same_lattice_dot(this->data(i),this->data(j),len,ang);
}

template<typename T> LQVec<T>& LQVec<T>:: operator+=(const LQVec<T>& av){
  assert( this->samelattice(av) );
  AVSizeInfo si = this->inplace_consistency_check(av);
  if (si.m != 3u) throw std::runtime_error("LQVecs should always have numel()==3!");
  for (size_t i=0; i<si.n; i++)
    for(size_t j=0; j<si.m; j++)
      this->insert(  this->getvalue(i,j) + av.getvalue(si.onevecb?0:i,si.singular?0:j), i,j );
  return *this;
}
template<typename T> LQVec<T>& LQVec<T>:: operator-=(const LQVec<T>& av){
  assert( this->samelattice(av) );
  AVSizeInfo si = this->inplace_consistency_check(av);
  if (si.m != 3u) throw std::runtime_error("LQVecs should always have numel()==3!");
  for (size_t i=0; i<si.n; i++)
    for(size_t j=0; j<si.m; j++)
      this->insert(  this->getvalue(i,j) - av.getvalue(si.onevecb?0:i,si.singular?0:j), i,j );
  return *this;
}


template<typename T> LQVec<T>& LQVec<T>:: operator +=(const T& av){
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0;j<this->numel(); j++)
      this->insert( this->getvalue(i,j) + av, i, j);
  return *this;
}
template<typename T> LQVec<T>& LQVec<T>:: operator -=(const T& av){
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0;j<this->numel(); j++)
      this->insert( this->getvalue(i,j) - av, i, j);
  return *this;
}
template<typename T> LQVec<T>& LQVec<T>:: operator *=(const T& av){
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0;j<this->numel(); j++)
      this->insert( this->getvalue(i,j) * av, i, j);
  return *this;
}
template<typename T> LQVec<T>& LQVec<T>:: operator /=(const T& av){
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0;j<this->numel(); j++)
      this->insert( this->getvalue(i,j) / av, i, j);
  return *this;
}

template<typename T> void LQVec<T>::check_arrayvector(const int flag){
  // Lattice vectors must always have 3 elements per vector. If we create
  // a LDVec from a lattice and a Direct lattice, we need to check this:
  size_t nel = this->numel();
  if (nel > 3u){
    if (flag) throw std::runtime_error("Lattice vectors require 3 elements -- if constructing LQVec(Reciprocal,ArrayVector), set optional flag to 0 to truncate input");
    this->removeelements(4u,nel);
  }
  if (nel < 3u) {
    if (flag) throw std::runtime_error("Lattice vectors require 3 elements -- if constructing LQVec(Reciprocal,ArrayVector), set optional flag to 0 to pad input");
    this->addelements(3u-nel);
  }

}
