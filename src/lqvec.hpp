

template<typename T> LQVec<T> LQVec<T>::get(const size_t i) const {
	LQVec<T> out(this->get_lattice(), i<this->size() ? 1u : 0u );
	if (i<this->size()) this->ArrayVector<T>::get(i, out.datapointer() );
	return out;
}

template<typename T> ArrayVector<T> LQVec<T>::get_hkl() const {
	return ArrayVector<T>(this->numel(),this->size(),this->data); // strip off the Lattice information
}
template<typename T> ArrayVector<double> LQVec<T>::get_xyz() const {
	double *toxyz = safealloc<double>(9);
	Reciprocal lat = this->get_lattice();
	lat.get_xyz_transform(toxyz);
	ArrayVector<double> xyz(this->numel(),this->size());
	for (size_t i=0; i<this->size(); i++) multiply_matrix_vector(xyz.datapointer(i), toxyz, this->datapointer(i));
	delete[] toxyz;
	return xyz;
}
template<typename T> LDVec<double> LQVec<T>::star() const {
	double *cvmt = safealloc<double>(9);
    (this->lattice).get_covariant_metric_tensor(cvmt);
    LDVec<double> slv( (this->lattice).star(), this->size() );
	for (size_t i=0; i<this->size(); i++) multiply_matrix_vector(slv.datapointer(i), cvmt, this->datapointer(i));
    slv /= 2.0*PI; // ai= gij/2/pi * ai_star
	delete[] cvmt;
	return slv;
}
// template<typename T> template<typename R> LQVec<double> LQVec<T>::cross(const LQVec<R> * vec) const {
// 	assert( this->samelattice(vec) );
// 	AVSizeInfo si = av_consistency_check(this,vec);
// 	if (si.m!=3u) throw "cross product is only defined for three vectors";
// 	LDVec<double> tmp((vec->get_lattice()).star(), si.n);
// 	for (size_t i=0; i<si.n; i++) vector_cross<double,T,R,3>(tmp.datapointer(i), this->datapointer(si.a?0:i), vec->datapointer(si.b?0:i));
//     tmp *= (vec->get_lattice()).get_volume()/2.0/PI;
//     return tmp.star();
// }
template<typename T> LQVec<double> LQVec<T>::cross(const size_t i, const size_t j) const {
    bool bothok = (i<this->size() && j<this->size() && 3u==this->numel());
	LQVec<double> out(this->get_lattice(), bothok? 1u : 0u);

    if (bothok){
        double *rlucross = safealloc<double>(3);
        vector_cross<double,T,T,3>(rlucross, this->datapointer(i), this->datapointer(j));

        Reciprocal rlat = this->get_lattice();

        LDVec<double> ldv( rlat.star(), 1u, rlucross);
        ldv *= rlat.get_volume()/2.0/PI;
        out =  ldv.star();

        delete[] rlucross;
    }
	return out;
}

template<typename T> double LQVec<T>::dot(const size_t i, const size_t j) const {
	Reciprocal lat = this->get_lattice();
	double len[3] = {lat.get_a(),lat.get_b(),lat.get_c()};
	double ang[3] = {lat.get_alpha(),lat.get_beta(),lat.get_gamma()};
	if (i<this->size()&&j<this->size())
		return same_lattice_dot(this->datapointer(i),this->datapointer(j),len,ang);
	throw "attempted out of bounds access by dot";
}

template<typename T> LQVec<T>& LQVec<T>:: operator+=(const LQVec<T>& av){
	assert( this->samelattice(av) );
	AVSizeInfo si = this->inplace_consistency_check(av);
	if (si.m != 3u) throw "LQVecs should always have numel()==3!\n";
	for (size_t i=0; i<si.n; i++)
		for(size_t j=0; j<si.m; j++)
			this->insert(  this->getvalue(i,j) + av.getvalue(si.onevecb?0:i,si.singular?0:j), i,j );
	return *this;
}
template<typename T> LQVec<T>& LQVec<T>:: operator-=(const LQVec<T>& av){
	assert( this->samelattice(av) );
	AVSizeInfo si = this->inplace_consistency_check(av);
	if (si.m != 3u) throw "LQVecs should always have numel()==3!\n";
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
