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
  size_t elements = 1u;
  size_t total_elements = this->scalar_elements
                        + this->eigenvector_num * this->eigenvector_dim
                        + this->vector_elements
                        + this->matrix_elements * this->matrix_elements;
  // no matter what, shape[0] should be the number of gridded points
  if (shape.size()>2){
    // if the number of dimensions of the shape array is greater than two,
    // the last element is the number of modes per point                     */
    for (size_t i=1; i<this->shape.size()-1; ++i) elements *= shape.getvalue(i);
    this->branches = shape.getvalue(shape.size()-1);
  } else {
    // shape is [n_points, n_elements], so there is only one mode
    elements = shape.getvalue(1u);
    this->branches = 1u;
  }
  if (0 == total_elements)
    this->scalar_elements = elements;
  if (total_elements && total_elements != elements){
    std::string msg ="Inconsistent element counts: "
                    + std::to_string(this->scalar_elements) + "+"
                    + std::to_string(this->eigenvector_dim) + "×"
                    + std::to_string(this->eigenvector_num) + "+"
                    + std::to_string(this->vector_elements) + "+"
                    + std::to_string(this->matrix_elements) + "² ≠ "
                    + std::to_string(elements);
    throw std::runtime_error(msg);
  }
}
template<class T>
int MapGrid3<T>::replace_data(const ArrayVector<T>& newdata,
                              const ArrayVector<size_t>& newshape,
                              const size_t new_scalar_elements,
                              const size_t new_eigenvector_num,
                              const size_t new_eigenvector_dim,
                              const size_t new_vector_elements,
                              const size_t new_matrix_elements){
  this->data =  newdata;
  this->shape = newshape;
  this->scalar_elements = new_scalar_elements;
  this->eigenvector_num = new_eigenvector_num;
  this->eigenvector_dim = new_eigenvector_dim;
  this->vector_elements = new_vector_elements;
  this->matrix_elements = new_matrix_elements;
  this->check_elements();
  return this->check_map();
}
template<class T>
int MapGrid3<T>::replace_data(const ArrayVector<T>& newdata,
                              const size_t new_scalar_elements,
                              const size_t new_eigenvector_num,
                              const size_t new_eigenvector_dim,
                              const size_t new_vector_elements,
                              const size_t new_matrix_elements){
  ArrayVector<size_t> shape(1,2);
  shape.insert(0,newdata.size());
  shape.insert(1,newdata.numel());
  return this->replace_data(newdata,
                            shape,
                            new_scalar_elements,
                            new_eigenvector_num,
                            new_eigenvector_dim,
                            new_vector_elements,
                            new_matrix_elements);
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
  for (size_t i=0; i<this->numel(); ++i) if (this->map[i]==m){
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
  ArrayVector<int> mzp = make_relative_neighbour_indices(1); // all combinations of [-1,0,+1] for three dimensions, skipping (0,0,0)
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

template<class T> T abs_diff(const T& A, const T& B){
  return std::abs( A - B );
}
template<class T> T abs_diff(const std::complex<T>& A, const std::complex<T>& B){
  return std::abs( std::real(A)-std::real(B) )+std::abs(std::imag(A)-std::imag(B));
}

/*! \brief Determine the sorting permutation connecting neighbouring gridded
           scalars, vectors, and/or matrices.

When multiple scalars/vectors/matrices are stored at each grid point of a
MapGrid3 object, it may be useful to identify which neighbouring values can be
ascribed to one-another. For each pair of neighbouring values this function
defines a cost and then uses the Munkres' Assignment algorithm to find the
minimum-cost value-value assignment.
In the case of scalar values, the cost Cᵃᵇᵢⱼ is |Vᵇⱼ-Vᵃᵢ|.
For vector values, the cost is the distance between the vectors √∑ₖ(Vᵇⱼₖ-Vᵃᵢₖ)².
And for matrix values, the cost is

It may also be necessary to combine information from scalars, vectors, and or
matrices in order to make a unique assignment. In such a case, the information
stored at each grid point should take the form of vector values with each vector
having elements [0,1,…,S-1,S,…,S+V-1,S+V,…,S+V+M-1] where S=`n_scalar`,
V=`n_vector`, and M=`n_matrix`.

@param n_scalar Number of scalar elements in the value vectors
@param n_vector Number of vector elements in the value vectors
@param n_matrix Number of matrix elements in the value vectors
*/
template <class T> template <class R>
ArrayVector<size_t> MapGrid3<T>::sort_perm(const size_t n_scalar,
                                           const size_t n_vector,
                                           const size_t n_matrix,
                                           const R scalar_weight,
                                           const R vector_weight,
                                           const R matrix_weight) const {
  /*
    The data contained in the MapGrid3 is often multiple scalars, vectors,
    matrices, or higher-order tensors. It is often important to identify which
    elements on neighbouring modes are *most likely* related to enable accurate
    interpolation between them.
    This function produces a sorting permutation of the stored values at each
    point in the mapped grid. The last dimension of the data is interpreted as
    representing different scalars/vectors/tensors to be sorted.
  */
  ArrayVector<size_t> elshape = this->data_shape(); // (data.size(),m,...,nobj) stored at each index of the data ArrayVector
  size_t span = 1;
  for(size_t i=1; i<elshape.size()-1; ++i) span *= elshape.getvalue(i);

  size_t nS=0, nV=0, nM=0;
  if (n_scalar + n_vector + n_matrix == span){
    nS = n_scalar;
    nV = n_vector;
    // check that the number of matrix elements describe a square matrix
    double matroot = std::sqrt( (double)n_matrix );
    if (matroot*matroot == (double)n_matrix){
      nM = (size_t)matroot;
    } else {
      nM = 0;
      nV += n_matrix;
    }
  } else {
    // try to be clever about what we're spanning:
    if (span < 11u){
      switch ( (int) span){
        case 1:  nS=1; nV=0; nM=0; break;
        case 2:  nS=0; nV=2; nM=0; break;
        case 3:  nS=0; nV=3; nM=0; break;
        case 4:  nS=1; nV=3; nM=0; break;
        case 5:  nS=1; nV=0; nM=2; break;
        case 6:  nS=0; nV=6; nM=0; break;
        case 7:  nS=1; nV=6; nM=0; break;
        case 8:  nS=2; nV=6; nM=0; break;
        case 9:  nS=0; nV=0; nM=3; break;
        case 10: nS=1; nV=0; nM=3; break;
      }
    } else {
      double spanroot = std::sqrt( (double)span );
      if (spanroot*spanroot == (double)span){
        nS=0; nV=0; nM=(size_t)spanroot;
      } else {
        spanroot = std::sqrt((double)(span-1));
        if (spanroot*spanroot == (double)(span-1)){
          nS=1; nV=0; nM=(size_t)spanroot;
        } else {
          nS=1; nV=span-1; nM=0;
        }
      }
    }
  }
  if ( nS+nV+nM*nM != span) throw std::runtime_error("problem determining number of scalar, vector, matrix, elements");

  size_t nobj = elshape.getvalue(elshape.size()-1);
  // within each index of the data ArrayVector there are span*nobj entries.
  // Each object consists of span elements, and we need to skip by span to get to the next one.

  // We will return the permutation of 0:nobj-1 which sorts the objects between neighbouring mapped points.
  ArrayVector<size_t> perm( nobj, this->data.size() );
  ArrayVector<bool> sorted(   1u, this->data.size() );
  for(size_t i=0; i<this->data.size(); ++i) sorted.insert(false,i);

  bool foundneighbour, firstnotfound = true;
  size_t this_idx, that_idx;
  ArrayVector<size_t> neighbours;

  // typename GridDiffTraits<T>::type *distance = new typename GridDiffTraits<T>::type[nobj*nobj]();
  typename GridDiffTraits<T>::type gdzero = typename GridDiffTraits<T>::type(0);
  typename GridDiffTraits<T>::type tval = gdzero, vval = gdzero, mval = gdzero;
  Munkres<typename GridDiffTraits<T>::type> munkres(nobj);

  typename GridDiffTraits<T>::type sscale = typename GridDiffTraits<T>::type(scalar_weight);
  typename GridDiffTraits<T>::type vscale = typename GridDiffTraits<T>::type(vector_weight);
  typename GridDiffTraits<T>::type mscale = typename GridDiffTraits<T>::type(matrix_weight);

  bool assigned_ok;
  size_t* assignment = new size_t[nobj]();

  T* Amat = new T[nM*nM]();
  T* Bmat = new T[nM*nM]();
  T* Avec = new T[nV]();
  T* Bvec = new T[nV]();
  T A, B;
  for(size_t idx=0; idx<this->numel(); ++idx){
    if (this->valid_mapping(idx)){
      this_idx = this->map[idx];
      if (firstnotfound){
        // the first valid mapping gets an arbitrarily assigned permutation
          for(size_t j=0; j<nobj; ++j) perm.insert(j, this_idx, j);
          firstnotfound = false;
          sorted.insert(true,this_idx);
      } else { //the normal part of the loop
        // look for a neighbouring already-sorted point
        neighbours = this->get_neighbours(idx);
        // printf("linear indexed neighbours of %u are:\n",idx); neighbours.print();
        foundneighbour = false;
        for(size_t j=0; j<neighbours.size(); ++j){
          if (this->valid_mapping(neighbours.getvalue(j))){
            that_idx = this->map[neighbours.getvalue(j)];
            if (sorted.getvalue(that_idx)) foundneighbour = true;
            // printf("linear indexed neighbour %u has mapping %u %s been sorted\n",neighbours.getvalue(j),that_idx, foundneighbour ? "and has" : "but has not");
          }
          if (foundneighbour) break;
        }
        if (!foundneighbour) throw std::runtime_error("something has gone wrong in sort_perm. no already-sorted neighbours!");
        // we have the index into this->data.size() for our current point
        // *AND* an already-sorted neighbouring point.

        // calculate the cost to assign each mode from the neighbour
        /* A note about indexing the data ArrayVector:
          An earlier version of this function used indexing into the ArrayVector
          elements like an array in Fortran,
          e.g., [mode number]*[number of elements]+[mode element]
          Instead, the correct way to index the ArrayVector elements is by
          [mode_number]+[mode_element]*[number of modes]
          due to the fact that a (N_points, N_elements, N_modes) C-style array
          was used to create the ArrayVector
          {with size=N_points and numel=N_elements*N_modes}
        */
        for (size_t i=0; i<nobj; ++i)
          for (size_t j=0; j<nobj; ++j){
            tval = gdzero;
            if (nS){
              for (size_t k=0; k<nS; ++k){
                A = this->data.getvalue(this_idx, i+k*nobj);
                B = this->data.getvalue(that_idx, j+k*nobj);
                tval += std::sqrt(squared_distance(A,B));
              }
            }
            if (nV){
              for (size_t k=0; k<nV; ++k){
                Avec[k] = this->data.getvalue(this_idx, i+(nS+k)*nobj);
                Bvec[k] = this->data.getvalue(that_idx, j+(nS+k)*nobj);
              }
              /* Modify the function to take a flag which selects the vector
                 cost function? Something like:
                    0 --> vector_angle
                    1 --> vector_distance
                    2 --> 1-vector_product
              */
              // vval = vector_angle(nV,Avec,Bvec);
              // vval = vector_distance(nV,Avec,Bvec);
              vval = 1-vector_product(nV,Avec,Bvec);  // 1 - |A*⋅B|² is 0 if A==B, 1 if A⟂B
            }
            if (nM){
              for (size_t ki=0; ki<nM; ++ki)
                for (size_t kj=0; kj<nM; ++kj){
                  Amat[ki*nM+kj] = this->data.getvalue(this_idx, i+(nS+nV+ki*nM+kj)*nobj);
                  Bmat[ki*nM+kj] = this->data.getvalue(that_idx, j+(nS+nV+ki*nM+kj)*nobj);
                }
              mval = frobenius_distance(nM,Amat,Bmat);
            }
            // for each i determine the cheapest j
            munkres.get_cost()[i*nobj+j] = sscale*tval + vscale*vval + mscale*mval;
            // // for each j determine the cheapest i
            // munkres.get_cost()[j*nobj+i] = sscale*tval + vscale*vval + mscale*mval;
          }
        // and use the Munkres' algorithm to determine the optimal assignment
        munkres.run_assignment();
        assigned_ok = munkres.get_assignment(assignment);

        if (!assigned_ok) throw std::runtime_error("The Munkres' assignment algorithm failed?!");
        // We want to return a *global* sorting permutation S₀ᵢ but what we
        // have from the Munkres' algorithm is a *local* permutation mapping
        // elements at index i (this_idx) onto elements at index j (that_idx).
        // Thankfully the local and global permutations have a handy
        // relationship
        //                Sₒᵢ = S₀ⱼ[Sⱼᵢ]
        // So, along with the arbitrary choice of S₀₁ made previously, we can
        // find S₀ᵢ for all i.
        for (size_t i=0; i<nobj; ++i)
          for (size_t j=0; j<nobj; ++j)
            if (perm.getvalue(that_idx,i)==assignment[j])
              perm.insert(j,this_idx,i);
        // for (size_t i=0; i<nobj; ++i)
        //   perm.insert(assignment[perm.getvalue(that_idx, i)], this_idx, i);
        sorted.insert(true,this_idx);
      }
    }
  }
  delete[] assignment;
  delete[] Amat;
  delete[] Bmat;
  delete[] Avec;
  delete[] Bvec;
  return perm;
}
