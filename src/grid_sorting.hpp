/*! Sorting values on a grid is non-trivial. This header file contains routines
    to help with this task.
*/

template<class T> std::vector<size_t> MapGrid3<T>::find_unsorted_neighbours(
  const std::vector<bool>& sorted, const size_t centre_map_idx) const {
  size_t centre_lin_idx;
  if (this->map2lin(centre_map_idx,centre_lin_idx))
    throw std::runtime_error("Mapping index has no corresponding linear index.");
  // get a listing of all in-bounds mapped neighbour linear indices
  ArrayVector<size_t> neighbours = this->get_neighbours(centre_lin_idx);
  size_t n_map_idx;
  std::vector<size_t> out;
  out.reserve(neighbours.size());
  for (size_t i=0; i<neighbours.size(); ++i){
    // sorted only contains valid-mapped point information
    // so convert from linear index to mapping index
    int ret = this->lin2map(neighbours.getvalue(i), n_map_idx);
    if (ret !=0){
      std::string msg = "Could not find map index of neighbour "
                      + std::to_string(neighbours.getvalue(i)) + " as ";
      if (ret < 0) msg += "the mapping is invalid.";
      if (ret > 0) msg += "the linear index is out of bounds.";
      throw std::runtime_error(msg);
    }
    if (!sorted[n_map_idx]) out.push_back(n_map_idx);
  }
  return out;
}

template<class T> std::vector<size_t> MapGrid3<T>::find_sorted_neighbours(
  const std::vector<bool>& sorted, const size_t clin) const {
  // get a listing of all in-bounds neighbour linear indexes
  ArrayVector<size_t> neighbours = this->get_neighbours(clin);
  size_t nidx, nnidx;
  size_t csub[3], nsub[3], nnsub[3];
  int dir[3];
  if (this->lin2sub(clin, csub))
    throw std::runtime_error("Could not find subscripted index for centre point.");

  std::vector<size_t> out;
  out.clear();
  bool possible;
  for (size_t i=0; i<neighbours.size(); ++i){
    // sorted only contains valid-mapped point information
    // so convert from linear index to mapping index
    int ret = this->lin2map(neighbours.getvalue(i), nidx);
    if (ret !=0){
      std::string msg = "Could not find map index of neighbour "
                      + std::to_string(neighbours.getvalue(i)) + " as ";
      if (ret < 0) msg += "the mapping is invalid.";
      if (ret > 0) msg += "the linear index is out of bounds.";
      throw std::runtime_error(msg);
    }
    if (sorted[nidx]){
      // make sure the output vector is empty to start:
      out.clear();
      // store this already sorted neighbour's mapped index
      out.push_back(nidx);
      // get the subscripted index for the sorted neighbour
      if (this->lin2sub(neighbours.getvalue(i), nsub))
        throw std::runtime_error("Could not find subscripted index for neighbour.");
      // determine the direction from the centre to the neighbour:
      // for (int j=0; j<3; ++j) dir[j] = nsub[j]-csub[j];
      for (int j=0; j<3; ++j) dir[j] = (nsub[j]>csub[j]) ? 1 : (nsub[j]<csub[j]) ? -1 : 0;
      // chech whether there *could* be a second nearest neighbour one more
      // step in the same direction as the first:
      possible = true;
      // first, can we represent nsub+dir as a size_t (will it stay positive)
      for (int j=0; j<3; ++j) possible &= (nsub[j]>0 || dir[j]>=0);
      if (possible){
        for (int j=0; j<3; ++j) nnsub[j] = nsub[j] + dir[j];
        // second, is nnsub a valid mapping index:
        possible = 0 == this->sub2map(nnsub, nnidx);
      }
      // Now, possible tells us if nnindx has been set to a valid mapping index
      // Check whether nnidx points to an already sorted mapped point
      if (possible && sorted[nnidx]){
        // std::cout << "Centre " << std::to_string(clin) << " (";
        // for (int j=0; j<3; ++j) std::cout << std::to_string(csub[j]) << " ";
        // std::cout << ") ";
        // std::cout << "Neighbour ( ";
        // for (int j=0; j<3; ++j) std::cout << std::to_string(nsub[j]) << " ";
        // std::cout << ") ";
        // std::cout << "Next Neighbour ( ";
        // for (int j=0; j<3; ++j) std::cout << std::to_string(nnsub[j]) << " ";
        // std::cout << ")" << std::endl;

        // if it does, add it to the output array and finish.
        out.push_back(nnidx);
        return out;
      }
    }
  }
  // We have now checked all neighbours and, at most, found sorted neighbours
  // but no sorted next neighbours. So we return the out vector which is empty
  // in the case of no sorted neighbours, or has the last-found neighbour
  // without a sorted next neighbour.

  // std::cout << "Centre " << std::to_string(clin) << " (";
  // for (int j=0; j<3; ++j) std::cout << std::to_string(csub[j]) << " ";
  // std::cout << ") ";
  // std::cout << "Neighbour ( ";
  // for (int j=0; j<3; ++j) std::cout << std::to_string(nsub[j]) << " ";
  // std::cout << ") " << std::endl;
  return out;
}



template<class T> template<class R>
bool MapGrid3<T>::sort_difference(const R scaleS,
                                  const R scaleE,
                                  const R scaleV,
                                  const R scaleM,
                                  const size_t span,
                                  const size_t nobj,
                                  ArrayVector<size_t>& perm,
                                  const size_t cidx,
                                  const size_t nidx,
                                  const int ecf,
                                  const int vcf) const {
return jv_permutation(this->data.datapointer(cidx,0),
                           this->data.datapointer(nidx,0),
                           this->scalar_elements,
                           this->eigvec_elements,
                           this->vector_elements,
                           this->matrix_elements,
                           scaleS, scaleE, scaleV, scaleM,
                           span, nobj, perm, cidx, nidx, ecf, vcf);
}

template<class T> template<class R>
bool MapGrid3<T>::sort_derivative(const R scaleS,
                                  const R scaleE,
                                  const R scaleV,
                                  const R scaleM,
                                  const size_t span,
                                  const size_t nobj,
                                  ArrayVector<size_t>& perm,
                                  const size_t cidx,
                                  const size_t nidx,
                                  const size_t nnidx,
                                  const int ecf,
                                  const int vcf) const {
ArrayVector<T> sorted(this->data.numel(),2u); // sorted neighbours
// Copy the data at each point, ensuring that the global permutation for
// nidx and nnidx are respected
size_t nn_i=0;
bool nn_i_found;
// std::cout << "Estimate the values of modes at the centre" << std::endl;
for (size_t i=0; i<nobj; ++i){
  nn_i_found = false;
  for (size_t j=0; j<nobj; ++j)
    if (perm.getvalue(nnidx,j)==perm.getvalue(nidx,i)){
      nn_i = j;
      nn_i_found = true;
    }
  if (!nn_i_found){
    std::string msg = "Next neighbour index not found.\n";
    msg += "perm[n]  = " + perm.to_string(nidx) + "\n";
    msg += "perm[nn] = " + perm.to_string(nnidx);
    throw std::runtime_error(msg);
  }
  for (size_t j=0; j<span; ++j){
    sorted.insert(this->data.getvalue(nidx,    i*span+j), 0u, i*span+j);
    sorted.insert(this->data.getvalue(nnidx,nn_i*span+j), 1u, i*span+j);
  }
}
// std::cout << sorted.to_string();
// calculate the predicted value of each element of the data.
// making the prediction is a special case of linear interpolation:
ArrayVector<T> predic(this->data.numel(),1u); // and prediction
size_t cnt = 2u;
size_t corners[2]{0u,1u};
T weights[2]{2,-1};
size_t out_i = 0u;
const size_t eigvec_dim = 3u;
const size_t eigvec_cnt = this->eigvec_elements/eigvec_dim;
if (eigvec_cnt*eigvec_dim != this->eigvec_elements)
  throw std::runtime_error("Expected 3-D eigenvectors. Please extend.");
unsafe_interpolate_to(sorted,
                      this->scalar_elements,
                      eigvec_cnt,
                      eigvec_dim,
                      this->vector_elements,
                      this->matrix_elements,
                      this->branches,
                      cnt,corners,weights,predic,out_i);
bool rslt;
// std::cout << predic.to_string();
// std::cout << this->data.to_string(cidx) << std::endl;
// find the assignment of each *predicted* object value to those at cidx:
rslt = jv_permutation(this->data.datapointer(cidx,0),
                           predic.datapointer(out_i,0),
                           this->scalar_elements,
                           this->eigvec_elements,
                           this->vector_elements,
                           this->matrix_elements,
                           scaleS, scaleE, scaleV, scaleM,
                           span, nobj, perm, cidx, nidx, ecf, vcf);
// std::cout << ((rslt)?"permutation determined":"permutation failed") << std::endl;
return rslt;
}

template<class T> template<typename R>
ArrayVector<size_t> MapGrid3<T>::new_sort_perm(const R scalar_weight,
                                               const R eigenv_weight,
                                               const R vector_weight,
                                               const R matrix_weight,
                                               const int ecf,
                                               const int vcf
                                             ) const {
//
// elshape → (data.size(),Nobj, Na, Nb, ..., Nz) stored at each grid point
// or (data.size(), Na, Nb, ..., Nz) ... but we have properties to help us out!
size_t nobj = this->branches;
size_t span = this->scalar_elements
            + this->eigvec_elements
            + this->vector_elements
            + this->matrix_elements*this->matrix_elements;
// within each index of the data ArrayVector there are span*nobj entries.

// We will return the permutations of 0:nobj-1 which sort the objects globally
ArrayVector<size_t> perm( nobj, this->data.size() );
// We need to keep track of which objects have been sorted thus far
std::vector<bool> sorted(this->data.size());
// To start with, none have been:
for (size_t i=0; i<this->data.size(); ++i) sorted[i] = false;

bool firstnotfound = true;
std::vector<size_t> nidx;

typename CostTraits<T>::type wS, wE, wV, wM;
wS = typename CostTraits<T>::type(scalar_weight);
wE = typename CostTraits<T>::type(eigenv_weight);
wV = typename CostTraits<T>::type(vector_weight);
wM = typename CostTraits<T>::type(matrix_weight);

size_t midx;
// std::cout << "Start the sorting" << std::endl;
for(size_t idx=0; idx<this->numel(); ++idx){
  if (this->valid_mapping(idx)){
    midx = this->map[idx];
    // std::cout << "linear index " << std::to_string(idx) << " mapping index "
    //           << std::to_string(midx) << std::endl;
    if (firstnotfound){
      // the first valid mapping gets an arbitrarily assigned permutation
        for(size_t j=0; j<nobj; ++j) perm.insert(j, midx, j);
        firstnotfound = false;
        sorted[midx]=true;
    } else { //the normal part of the loop
      // std::cout << "look for neighbours" << std::endl;
      // look for a neighbouring already-sorted point
      nidx = this->find_sorted_neighbours(sorted, idx);
      if (nidx.size()<1 || nidx.size()>2)
        throw std::runtime_error(
          "Logic error. Too few or too many sorted neighbours found.");
      // only one sorted neighbour, so punt
      // if (nidx.size()>0){
      if (nidx.size()==1){
        // std::cout << "one neighbour" << std::endl;
        sorted[midx] = this->sort_difference(wS,wE,wV,wM,span,nobj,perm,
                                             midx,nidx[0],ecf,vcf);
      }
      // two sorted neighbours in a line, use the drivative method
      if (nidx.size()==2){
        // std::cout << "two neighbours" << std::endl;
        sorted[midx] = this->sort_derivative(wS,wE,wV,wM,span,nobj,perm,
                                             midx,nidx[0],nidx[1],ecf,vcf);
      }
    }
  }
}
return perm;
}


template<class T> template<typename R>
ArrayVector<size_t> MapGrid3<T>::centre_sort_perm(const R scalar_weight,
                                                  const R eigenv_weight,
                                                  const R vector_weight,
                                                  const R matrix_weight,
                                                  const int ecf,
                                                  const int vcf
) const {
//
// elshape → (data.size(),Nobj, Na, Nb, ..., Nz) stored at each grid point
// or (data.size(), Na, Nb, ..., Nz) ... but we have properties to help us out!
size_t nobj = this->branches;
size_t span = this->scalar_elements
            + this->eigvec_elements
            + this->vector_elements
            + this->matrix_elements*this->matrix_elements;
// within each index of the data ArrayVector there are span*nobj entries.

// We will return the permutations of 0:nobj-1 which sort the objects globally
ArrayVector<size_t> perm( nobj, this->data.size() );
// We need to keep track of which objects have been sorted thus far
std::vector<bool> sorted(this->data.size());
for (size_t i=0; i<this->data.size(); ++i) sorted[i] = false;

typename CostTraits<T>::type wS, wE, wV, wM;
wS = typename CostTraits<T>::type(scalar_weight);
wE = typename CostTraits<T>::type(eigenv_weight);
wV = typename CostTraits<T>::type(vector_weight);
wM = typename CostTraits<T>::type(matrix_weight);

size_t mapidx;
size_t subidx[3];
// Start from the Γ point. It should always be a valid mapped point.
// The grid centre along each direction is the ((2N+1)/2 + 1)ᵗʰ entry, but we count from 0 so (2N+1)/2
for (int i=0; i<3; ++i) subidx[i] = this->size(i) >> 1;
if (this->sub2map(subidx, mapidx))
  throw std::runtime_error("Γ is not a mapped point.");
// The first point gets an arbitrary assignment
for (size_t j=0; j<nobj; ++j) perm.insert(j, mapidx, j);
// Sort ±neighbours the same way; [1 0 0] & [-1 0 0], [1 1 0] & [-1 -1 0] …
/* We could try some fancy stuff looking for q = -q pairs, but
   a) the way find_unsorted_neighbours works, they should always be equal
      distance from the start and end of the returned vector, e.g., if
      all three dimensions are present, i & 26-i are -q & +q pairs for the
      neighbours of Γ
   b) Treating them the same way only entails finding their sorting permutation
      *relative* to Γ (and not the reverse) -- so the order in which we visit
      the neighbours doesn't matter.                                          */
std::vector<size_t> unsorted_neighbours = this->find_unsorted_neighbours(sorted, mapidx);
size_t num_sorted = 0; // Γ already sorted, but this is effectively an index
for (size_t n_mapidx: unsorted_neighbours){
  sorted[n_mapidx] = this->sort_difference(wS,wE,wV,wM,span,nobj,perm,n_mapidx,mapidx,ecf,vcf);
  if (!sorted[n_mapidx]) throw std::runtime_error("Failed to sort Γ neighbour.");
  ++num_sorted;
}
// with the immediate neighbours sorted, start a recursive search (from each to
// be extra sure that we reach everywhere in the Brillouin zone.
std::vector<bool> locked(this->data.size());
for (size_t i=0; i<this->data.size(); ++i) locked[i] = false;
for (size_t n_mapidx: unsorted_neighbours){
  num_sorted += this->sort_recursion(n_mapidx,wS,wE,wV,wM,ecf,vcf,span,nobj,perm,sorted,locked);
}
if (num_sorted < this->data.size())
  throw std::runtime_error("Recursive sorting failed to find all grid points.");
if (num_sorted > this->data.size())
  std::cout << "Sorted " << std::to_string(num_sorted) << " grid points but only "
            << std::to_string(this->data.size()) << " points are mapped."
            << std::endl;
  // throw std::runtime_error("Recursive sorting visited some grid points more than once.");
return perm;
}

template<class T> template<class R>
size_t MapGrid3<T>::sort_recursion(const size_t centre,
                                   const R wS, const R wE, const R wV, const R wM,
                                   const int ecf, const int vcf,
                                   const size_t span, const size_t nobj,
                                   ArrayVector<size_t>& perm,
                                   std::vector<bool>& sorted,
                                   std::vector<bool>& locked) const {
// std::cout << "Recursive sorting" << std::endl;
std::vector<size_t> unsorted_neighbours = this->find_unsorted_neighbours(sorted, centre);
std::vector<size_t> sorted_neighbours;
size_t nlin, num_sorted=0;
bool success=false;
for (size_t nmap: unsorted_neighbours){
  // When this gets parallelised two threads might have overlapping neighbours
  // Don't bother doing anything if another thread got to this neighbour after
  // unsorted_neighbours was called.
  if (sorted[nmap]||locked[nmap]) break;
  locked[nmap] = true;
  this->map2lin(nmap, nlin);
  sorted_neighbours = this->find_sorted_neighbours(sorted, nlin);
  if (sorted_neighbours.size()<1)
    throw std::runtime_error("No sorted neighbours.");
  if (sorted_neighbours.size()>2)
    throw std::runtime_error("Too many sorted neighbours.");
  success = this->sort_difference(wS, wE, wV, wM, span, nobj, perm, nmap,
                                  sorted_neighbours[0],
                                  ecf, vcf);
  // if (sorted_neighbours.size()==1)
  //   success = this->sort_difference(wS, wE, wV, wM, span, nobj, perm, nmap,
  //                                   sorted_neighbours[0],
  //                                   ecf, vcf);
  // if (sorted_neighbours.size()==2)
  //   success = this->sort_derivative(wS, wE, wV, wM, span, nobj, perm, nmap,
  //                                   sorted_neighbours[0], sorted_neighbours[1],
  //                                   ecf, vcf);

  // if (!success) throw std::runtime_error("Failed to find permutation.");
  // Don't throw here, maybe we can sort from a different direction?
  sorted[nmap] = success;
  if (success) num_sorted++;
  locked[nmap] = false;
}
// All neighbours of map_idx sorted, so now sort all of their neighbours in turn
for (size_t un: unsorted_neighbours)
  num_sorted += this->sort_recursion(un, wS, wE, wV, wM, ecf, vcf, span, nobj, perm, sorted, locked);

return num_sorted;
}
