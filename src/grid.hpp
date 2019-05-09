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
  this->lin2sub(centre, ijk.datapointer(0)); // get the subscripted indices of the centre position
  bool isz[3], ism[3]; // is the centre index 0 (isz) or the maximum (ism)
  for (size_t i=0; i<3u; ++i) isz[i] = 0==ijk.getvalue(0,i);
  for (size_t i=0; i<3u; ++i) ism[i] = this->size(i)-1 <= ijk.getvalue(0,i);
  ArrayVector<bool> is_valid(1u,mzp.size());
  for (size_t i=0; i<mzp.size(); ++i){
    // keep track of if we *can* (or should) add each mzp vector to the centre index
    is_valid.insert(true,i);
    for (size_t j=0; j<mzp.numel(); ++j){
      if (isz[j] && mzp.getvalue(i,j)<0 ) is_valid.insert(false,i);
      if (ism[j] && mzp.getvalue(i,j)>0 ) is_valid.insert(false,i);
    }
  }
  ArrayVector<size_t> tmp(3u,1u);
  for (size_t i=0; i<mzp.size(); ++i){
    if (is_valid.getvalue(i)){
      for (size_t j=0; j<3u; ++j) tmp.insert( ijk.getvalue(0,j) + mzp.getvalue(i,j), 0, j);
      is_valid.insert( this->is_inbounds(tmp.datapointer(0)) ,i); //ensure we only check in-bounds neighbours
    }
  }
  size_t valid_neighbours = 0;
  for (size_t i=0; i<is_valid.size(); ++i) if (is_valid.getvalue(i)) ++valid_neighbours;
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

template<class T> T squared_distance(const T&A, const T& B){
  return (A-B)*(A-B);
}
template<class T> T squared_distance(const std::complex<T>& A, const std::complex<T>& B){
  T r = std::real(A)-std::real(B);
  T i = std::imag(A)-std::imag(B);
  return r*r + i*i;
}


// template<class T> ArrayVector<size_t> MapGrid3<T>::sort_perm() const {
//   /*
//     The data contained in the MapGrid3 is often multiple scalars, vectors,
//     matrices, or higher-order tensors. It is often important to identify which
//     elements on neighbouring modes are *most likely* related to enable accurate
//     interpolation between them.
//     This function produces a sorting permutation of the stored values at each
//     point in the mapped grid. The last dimension of the data is interpreted as
//     representing different scalars/vectors/tensors to be sorted.
//   */
//   ArrayVector<size_t> elshape = this->data_shape(); // (data.size(),m,...,nobj) stored at each index of the data ArrayVector
//   size_t span = 1;
//
//   for(size_t i=1; i<elshape.size()-1; ++i) span *= elshape.getvalue(i);
//   size_t nobj = elshape.getvalue(elshape.size()-1);
//   // within each index of the data ArrayVector there are span*nobj entries.
//   // Each object consists of span elements, and we need to skip by span to get to the next one.
//
//   // We will return the permutation of 0:nobj-1 which sorts the objects between neighbouring mapped points.
//   ArrayVector<size_t> perm( nobj, this->data.size() );
//   ArrayVector<bool> hasbeensorted(1u, this->data.size() );
//   for(size_t i=0; i<this->data.size(); ++i) hasbeensorted.insert(false,i);
//
//   bool isvalid, foundneighbour, foundfirst = false;
//   size_t linidx;
//   size_t mapidx;
//   size_t neigmi;
//   ArrayVector<size_t> neighbours;
//
//   size_t *mink  = new size_t[nobj]();
//   size_t *neqv  = new size_t[nobj]();
//   typename GridDiffTraits<T>::type *sumabsdiff = new typename GridDiffTraits<T>::type[nobj*nobj]();
//   typename GridDiffTraits<T>::type gdzero = typename GridDiffTraits<T>::type(0);
//   typename GridDiffTraits<T>::type tmind;
//   size_t tmink;
//   bool anydegenerate = false;
//   bool *unassigned = new bool[nobj*nobj]();
//   T A, B;
//   for(size_t idx=0; idx<this->numel(); ++idx){
//     if (this->valid_mapping(idx)){
//       mapidx = this->map[idx];
//       if (foundfirst){ //the normal part of the loop
//         // look for a neighbouring already-sorted point
//         neighbours = this->get_neighbours(idx);
//         // printf("linear indexed neighbours of %u are:\n",idx); neighbours.print();
//         foundneighbour = false;
//         for(size_t j=0; j<neighbours.size(); ++j){
//           if (this->valid_mapping(neighbours.getvalue(j))){
//             neigmi = this->map[neighbours.getvalue(j)];
//             if (hasbeensorted.getvalue(neigmi)) foundneighbour = true;
//             // printf("linear indexed neighbour %u has mapping %u %s been sorted\n",neighbours.getvalue(j),neigmi, foundneighbour ? "and has" : "but has not");
//           }
//           if (foundneighbour) break;
//         }
//         if (!foundneighbour) throw std::runtime_error("something has gone wrong in sort_perm. no already-sorted neighbours!");
//         // we have the index into this->data.size() for our current point
//         // *AND* an already-sorted neighbouring point.
//
//         // now *all* we have to do is the actual sorting :/
//         for (size_t i=0; i<nobj*nobj; ++i) sumabsdiff[i]=gdzero;
//         for (size_t i=0; i<nobj; ++i)
//           for (size_t j=0; j<nobj; ++j)   // std::abs( std::complex<R> ) returns a type R value
//             for (size_t k=0; k<span; ++k) {
//               A = this->data.getvalue(mapidx, i*span+k);
//               B = this->data.getvalue(neigmi, j*span+k);
//               sumabsdiff[i*nobj+j] += abs_diff(A,B);
//             }
//         // we want to select out the minimum index, k, for each i
//         // hopefully only one k will be smallest for each i, and the mapping will be singular
//         anydegenerate = false;
//         for (size_t i=0; i<nobj*nobj; ++i) unassigned[i]=true;
//         size_t total=0, count=0;
//         long long Nunassigned = (long long)nobj;
//         const size_t max_evals{1000000};
//         do{
//           if (anydegenerate && count++ > Nunassigned) {
//             /* for there to be a degenerate set of sumabsdiff[i,:] either
//                1) there are equal-distance modes on either side
//             or 2) there are degenerate modes in j
//
//             In the first case selecting a j for the other i's should resolve
//             the degeneracy.
//             -- try not picking any at random Nunassigned times
//             In the second case, we can safely pick a j at random.
//             */
//             /* TODO far far in the future:
//                   A more sophisticated version of this algorithm would attempt
//                   to pull together information from other neighbours to resolve
//                   the degeneracy; effectively computing the partial derivatives
//             */
//             for (size_t i=0; i<nobj; ++i){
//               if (neqv[i] > 1){
//                 typename GridDiffTraits<T>::type *tmp_sad = new typename GridDiffTraits<T>::type[neqv[i]-1]();
//                 size_t tmp_neq = 0;
//                 // the first *should* be mink[i], but it's possible that mink[i]
//                 // has already been assigned to another mode; so we need to check
//                 if (!unassigned[i*nobj+mink[i]])
//                   for (size_t j=0; j<nobj; ++j)
//                     if (unassigned[i*nobj+j] && sumabsdiff[i*nobj+j] == sumabsdiff[i*nobj+mink[i]]){
//                       mink[i] = j;
//                       break;
//                     }
//                 // now mink[i] is the first unassigned mode with this sumabsdiff
//                 for (size_t j=mink[i]+1; j<nobj; ++j){
//                   if (unassigned[i*nobj+j] && sumabsdiff[i*nobj+j] == sumabsdiff[i*nobj+mink[i]]){
//                     tmp_sad[tmp_neq] = gdzero;
//                     for (size_t k=0; k<span; ++k){
//                       A = this->data.getvalue(neigmi, mink[i]*span+k);
//                       B = this->data.getvalue(neigmi, j*span+k);
//                       tmp_sad[i*nobj+j] += abs_diff(A,B);
//                     }
//                     tmp_neq++;
//                   }
//                 }
//                 delete[] tmp_sad;
//                 if (tmp_neq > 0){ // tmp_neq+1 equivalent modes at the neighbouring site.
//                   // pick the first one
//                   perm.insert( perm.getvalue(neigmi, mink[i]), mapidx, i);
//                   for (size_t j=0; j<nobj; ++j) unassigned[j*nobj+mink[i]]=false;
//                   Nunassigned -= 1;
//                   neqv[i] = 0;
//                   count = 0;
//                 }
//               }
//             }
//           }
//           anydegenerate = false;
//           if (Nunassigned > 0){
//             for (size_t i=0; i<nobj; ++i){
//               tmind=GridDiffTraits<T>::max;
//               for (size_t j=0; j<nobj; ++j){
//                 if (unassigned[i*nobj+j]){
//                   if (sumabsdiff[i*nobj+j] < tmind){
//                     tmind = sumabsdiff[i*nobj+j];
//                     tmink = j;
//                     neqv[i] = 1;
//                   } else {
//                     if (sumabsdiff[i*nobj+j] == tmind) neqv[i]+=1;
//                   }
//                   mink[i] = tmink;
//                 }
//               }
//               if (neqv[i]>1) anydegenerate=true;
//             }
//             // assign any non-degenerate cases:
//             for (size_t i=0; i<nobj; ++i){
//               if (neqv[i] == 1){
//                 perm.insert( perm.getvalue(neigmi,mink[i]), mapidx, i );
//                 for (size_t j=0; j<nobj; ++j) unassigned[j*nobj+mink[i]] = false;
//                 Nunassigned -= 1;
//                 neqv[i] = 0;
//               }
//             }
//           }
//         } while (anydegenerate && ++total < max_evals);
//         // exiting the while loop means *either* all permutations have been
//         // assigned or we have reached the maximum number of evaluations.
//         if (total >= max_evals) {
//           throw std::runtime_error("maximum loop evaluations reached");
//         }
//         hasbeensorted.insert(true,mapidx);
//         total = 0;
//       } else { // the first valid mapping gets arbitrarily assigned permutation
//         for(size_t j=0; j<nobj; ++j) perm.insert(j, mapidx, j);
//         foundfirst = true;
//         hasbeensorted.insert(true,mapidx);
//       }
//     }
//   }
//   delete[] mink;
//   delete[] neqv;
//   delete[] sumabsdiff;
//   delete[] unassigned;
//
//   return perm;
// }


template<class T> ArrayVector<size_t> MapGrid3<T>::sort_perm() const {
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
  // try to be clever about what we're spanning:
  bool matspan = false;
  size_t matsize = 0, matstart=0;
  if (span < 11u){
    switch ( (int) span){
      case 4:  matspan=true; matstart=0u; matsize=2u; break;
      case 5:  matspan=true; matstart=1u; matsize=2u; break;
      case 9:  matspan=true; matstart=0u; matsize=3u; break;
      case 10: matspan=true; matstart=1u; matsize=3u; break;
    }
  } else {
    double spanroot = std::sqrt( (double)span );
    matspan = (spanroot*spanroot) == (double)span;
    if (matspan) matsize = (size_t) spanroot;
    if (!matspan){
      spanroot = std::sqrt( (double)(span-1) );
      matspan = (spanroot*spanroot) == (double)(span-1);
      if (matspan) {
        matsize = (size_t) spanroot;
        matstart = 1u;
      }
    }
  }


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
  typename GridDiffTraits<T>::type tval = gdzero;
  Munkres<typename GridDiffTraits<T>::type> munkres(nobj);

  bool assigned_ok;
  size_t* assignment = new size_t[nobj]();

  T* Amat = new T[matsize*matsize]();
  T* Bmat = new T[matsize*matsize]();
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

        // calculate the cost (distance) to assign each mode from the neighbour
        for (size_t i=0; i<nobj; ++i)
          for (size_t j=0; j<nobj; ++j){
            tval = gdzero;
            if (matspan){
              if (matstart > 0){
                A = this->data.getvalue(this_idx, i*span);
                B = this->data.getvalue(that_idx, j*span);
                tval += squared_distance(A,B);
              }
              for (size_t ki=0; ki<matspan; ++ki)
                for (size_t kj=0; kj<matspan; ++kj){
                  Amat[ki*matspan+kj] = this->data.getvalue(this_idx, i*span+matstart+ki*matspan+kj);
                  Bmat[ki*matspan+kj] = this->data.getvalue(that_idx, j*span+matstart+ki*matspan+kj);
                }
              tval += frobenius_distance(Amat,Bmat,matspan);
            } else {
              for (size_t k=0; k<span; ++k) {
                A = this->data.getvalue(this_idx, i*span+k);
                B = this->data.getvalue(that_idx, j*span+k);
                tval += squared_distance(A,B);
              }
            }
            munkres.get_cost()[i*nobj+j] = std::sqrt(tval);
          }
        // and use the Munkres' algorithm to determine the optimal assignment
        assigned_ok = munkres.get_assignment(assignment);
        if (!assigned_ok) throw std::runtime_error("The Munkres' assignment algorithm failed?!");
        for (size_t i=0; i<nobj; ++i) perm.insert( assignment[i], this_idx, i);
        sorted.insert(true,this_idx);
      }
    }
  }
  delete[] assignment;
  delete[] Amat;
  delete[] Bmat;
  return perm;
}
