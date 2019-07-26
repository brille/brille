template<typename T> LDVec<T> LDVec<T>::extract(const size_t i) const {
  if (i<this->size()){
    LDVec<T> out(this->get_lattice(),1u,this->datapointer(i));
    return out;
  }
  throw std::out_of_range("The requested element of the ArrayVector does not exist");
}
template<typename T> LDVec<T> LDVec<T>::first(const size_t num) const {
  size_t stop = num < this->size() ? num : this->size();
  LDVec<T> out(this->get_lattice(), stop);
  for (size_t j=0; j<stop; j++) out.set(j, this->datapointer(j) );
  return out;
}
template<typename T> LDVec<T> LDVec<T>::extract(const size_t n, const size_t *i) const {
  bool allinbounds = true;
  LDVec<T> out(this->get_lattice(),0u);
  for (size_t j=0; j<n; j++) if ( !(i[j]<this->size()) ){ allinbounds=false; break; }
  if (allinbounds){
    out.resize(n);
    for (size_t j=0; j<n; j++) out.set(j, this->datapointer(i[j]) );
  }
  return out;
}
template<typename T> LDVec<T> LDVec<T>::extract(const ArrayVector<size_t>& idx) const{
  bool allinbounds = true;
  LDVec<T> out(this->get_lattice(),0u);
  if (idx.numel() != 1u) throw std::runtime_error("copying an ArrayVector by index requires ArrayVector<size_t> with numel()==1 [i.e., an ArrayScalar]");
  for (size_t j=0; j<idx.size(); ++j) if (idx.getvalue(j)>=this->size()){allinbounds=false; break;}
  if (allinbounds){
    out.resize(idx.size());
    for (size_t j=0; j<idx.size(); ++j) out.set(j, this->datapointer( idx.getvalue(j)) );
  }
  return out;
}
template<typename T> LDVec<T> LDVec<T>::extract(const ArrayVector<bool>& tf) const{
  if (tf.numel() != 1u || tf.size() != this->size()){
    std::string msg = "Extracting an LQVec by logical indexing requires";
    msg += " an ArrayVector<bool> with numel()==1";
    msg += " and size()==LQVec.size().";
    throw std::runtime_error(msg);
  }
  size_t nout=0;
  for (size_t i=0; i<tf.size(); ++i) if (tf.getvalue(i,0)) ++nout;
  LDVec<T> out(this->get_lattice(),nout);
  size_t idx = 0;
  for (size_t i=0; i<tf.size(); ++i)
    if (tf.getvalue(i,0)) out.set(idx++, this->datapointer(i));
  return out;
}

template<typename T> LDVec<T> LDVec<T>::get(const size_t i) const {
  LDVec<T> out(this->get_lattice(), i<this->size() ? 1u : 0u );
  if (i<this->size()) this->ArrayVector<T>::get(i, out.datapointer() );
  return out;
}

template<typename T> ArrayVector<T> LDVec<T>::get_hkl() const {
  return ArrayVector<T>(this->numel(),this->size(),this->data); // strip off the Lattice information
}
template<typename T> ArrayVector<double> LDVec<T>::get_xyz() const {
  double *toxyz = safealloc<double>(9);
  Reciprocal lat = this->get_lattice();
  lat.get_xyz_transform(toxyz);
  ArrayVector<double> xyz(this->numel(),this->size());
  for (size_t i=0; i<this->size(); i++) multiply_matrix_vector(xyz.datapointer(i), toxyz, this->datapointer(i));
  delete[] toxyz;
  return xyz;
}
template<typename T> LQVec<double> LDVec<T>::star() const {
  double *cvmt = safealloc<double>(9);
  (this->lattice).get_covariant_metric_tensor(cvmt);
  LQVec<double> slv( (this->lattice).star(), this->size() );
  for (size_t i=0; i<this->size(); i++) multiply_matrix_vector(slv.datapointer(i), cvmt, this->datapointer(i));
  slv /= 2.0*PI; // ai= gij/2/pi * ai_star
  delete[] cvmt;
  return slv;
}
template<typename T> LDVec<double> LDVec<T>::cross(const size_t i, const size_t j) const {
  bool bothok = (i<this->size() && j<this->size() && 3u==this->numel());
  LDVec<double> out(this->get_lattice(), bothok? 1u : 0u);
  if (bothok){
    double *rlucross = safealloc<double>(3);
    vector_cross<double,T,T,3>(rlucross, this->datapointer(i), this->datapointer(j));

    Direct dlat = this->get_lattice();
    LQVec<double> lqv( dlat.star(), 1u, rlucross);
    lqv *= dlat.get_volume()/2.0/PI;
    out =  lqv.star();

    delete[] rlucross;
  }
  return out;
}

template<typename T> double LDVec<T>::dot(const size_t i, const size_t j) const {
  Direct lat = this->get_lattice();
  double len[3] = {lat.get_a(),lat.get_b(),lat.get_c()};
  double ang[3] = {lat.get_alpha(),lat.get_beta(),lat.get_gamma()};
  if (i<this->size()&&j<this->size())
    return same_lattice_dot(this->datapointer(i),this->datapointer(j),len,ang);
  throw "attempted out of bounds access by dot";
}

template<typename T> LDVec<T>& LDVec<T>:: operator+=(const LDVec<T>& av){
  assert( this->samelattice(av) );
  AVSizeInfo si = this->inplace_consistency_check(av);
  if (si.m != 3u) throw "LDVecs should always have numel()==3!\n";
  for (size_t i=0; i<si.n; i++)
    for(size_t j=0; j<si.m; j++)
      this->insert(  this->getvalue(i,j) + av.getvalue(si.onevecb?0:i,si.singular?0:j), i,j );
  return *this;
}
template<typename T> LDVec<T>& LDVec<T>:: operator-=(const LDVec<T>& av){
  assert( this->samelattice(av) );
  AVSizeInfo si = this->inplace_consistency_check(av);
  if (si.m != 3u) throw "LDVecs should always have numel()==3!\n";
  for (size_t i=0; i<si.n; i++)
    for(size_t j=0; j<si.m; j++)
      this->insert(  this->getvalue(i,j) - av.getvalue(si.onevecb?0:i,si.singular?0:j), i,j );
  return *this;
}

template<typename T> LDVec<T>& LDVec<T>:: operator +=(const T& av){
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0;j<this->numel(); j++)
      this->insert( this->getvalue(i,j) + av, i, j);
  return *this;
}
template<typename T> LDVec<T>& LDVec<T>:: operator -=(const T& av){
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0;j<this->numel(); j++)
      this->insert( this->getvalue(i,j) - av, i, j);
  return *this;
}
template<typename T> LDVec<T>& LDVec<T>:: operator *=(const T& av){
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0;j<this->numel(); j++)
      this->insert( this->getvalue(i,j) * av, i, j);
  return *this;
}
template<typename T> LDVec<T>& LDVec<T>:: operator /=(const T& av){
  for (size_t i=0; i<this->size(); i++)
    for (size_t j=0;j<this->numel(); j++)
      this->insert( this->getvalue(i,j) / av, i, j);
  return *this;
}

template<typename T> void LDVec<T>::check_arrayvector(const int flag){
  // Lattice vectors must always have 3 elements per vector. If we create
  // a LDVec from a lattice and a Direct lattice, we need to check this:
  size_t nel = this->numel();
  if (nel > 3u){
    if (flag) throw std::runtime_error("Lattice vectors require 3 elements -- if constructing LDVec(Direct,ArrayVector), set optional flag to 0 to truncate input");
    this->removeelements(4u,nel);
  }
  if (nel < 3u) {
    if (flag) throw std::runtime_error("Lattice vectors require 3 elements -- if constructing LDVec(Direct,ArrayVector), set optional flag to 0 to pad input");
    this->addelements(3u-nel);
  }

}
