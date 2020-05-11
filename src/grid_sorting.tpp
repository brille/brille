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

/*! Sorting values on a grid is non-trivial. This header file contains
    templated implementation routines to help with this task. */
template<class T,class R>
std::tuple<std::vector<size_t>,std::vector<size_t>,std::vector<size_t>>
MapGrid3<T,R>::classify_neighbours(const std::vector<bool>& sorted, const size_t centre_map_idx) const {
  size_t centre_lin_idx, n_map_idx, nn_map_idx;
  size_t csub[3], nsub[3], nnsub[3];
  if (this->map2lin(centre_map_idx,centre_lin_idx))
    throw std::runtime_error("Mapping index has no corresponding linear index.");
  if (this->lin2sub(centre_lin_idx, csub))
    throw std::runtime_error("Could not find subscripted index for centre point.");

  // get a listing of all in-bounds mapped neighbour linear indices
  ArrayVector<size_t> neighbours = this->get_neighbours(centre_lin_idx);
  std::vector<size_t> sorted_neighbours;
  std::vector<size_t> sorted_next_neighbours;
  std::vector<size_t> unsorted_neighbours;
  sorted_neighbours.reserve(neighbours.size());
  sorted_next_neighbours.reserve(neighbours.size());
  unsorted_neighbours.reserve(neighbours.size());
  int dir[3];
  bool possible;
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
    if (sorted[n_map_idx]){
      // store this already sorted neighbour's mapped index
      sorted_neighbours.push_back(n_map_idx);
      // get the subscripted index for the sorted neighbour
      if (this->lin2sub(neighbours.getvalue(i), nsub))
        throw std::runtime_error("Could not find subscripted index for neighbour.");
      // determine the direction from the centre to the neighbour:
      for (int j=0; j<3; ++j) dir[j] = (nsub[j]>csub[j]) ? 1 : (nsub[j]<csub[j]) ? -1 : 0;
      // chech whether there *could* be a second nearest neighbour one more
      // step in the same direction as the first:
      possible = true;
      // first, can we represent nsub+dir as a size_t (will it stay positive)
      for (int j=0; j<3; ++j) possible &= (nsub[j]>0 || dir[j]>=0);
      if (possible){
        for (int j=0; j<3; ++j) nnsub[j] = nsub[j] + dir[j];
        // second, is nnsub a valid mapping index:
        possible = 0 == this->sub2map(nnsub, nn_map_idx);
      }
      // Now, possible tells us if nnindx has been set to a valid mapping index
      // Check whether nnidx points to an already sorted mapped point
      sorted_next_neighbours.push_back(
        (possible && sorted[nn_map_idx]) ? nn_map_idx : n_map_idx
      );
    }
    else unsorted_neighbours.push_back(n_map_idx);
  }
  /* We now have three vectors containing the mapped indices of unsorted nearest
     neighbours, sorted nearest neighbours, and sorted next-nearest neighbours
     and will return all three as a single tuple */
  return std::make_tuple(unsorted_neighbours, sorted_neighbours, sorted_next_neighbours);
}

template<class T,class R>
std::tuple<std::vector<size_t>,std::vector<size_t>>
MapGrid3<T,R>::classify_sorted_neighbours(const std::vector<bool>& sorted, const size_t centre_map_idx) const {
  size_t centre_lin_idx, n_map_idx, nn_map_idx{0};
  size_t csub[3], nsub[3], nnsub[3];
  if (this->map2lin(centre_map_idx,centre_lin_idx))
    throw std::runtime_error("Mapping index has no corresponding linear index.");
  if (this->lin2sub(centre_lin_idx, csub))
    throw std::runtime_error("Could not find subscripted index for centre point.");
  // get a listing of all in-bounds mapped neighbour linear indices
  ArrayVector<size_t> neighbours = this->get_neighbours(centre_lin_idx);
  std::vector<size_t> sorted_neighbours;
  std::vector<size_t> sorted_next_neighbours;
  sorted_neighbours.reserve(neighbours.size());
  sorted_next_neighbours.reserve(neighbours.size());
  int dir[3];
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
    if (sorted[n_map_idx]){
      // store this already sorted neighbour's mapped index
      sorted_neighbours.push_back(n_map_idx);
      // get the subscripted index for the sorted neighbour
      if (this->lin2sub(neighbours.getvalue(i), nsub))
        throw std::runtime_error("Could not find subscripted index for neighbour.");
      // determine the direction from the centre to the neighbour:
      for (int j=0; j<3; ++j) dir[j] = (nsub[j]>csub[j]) ? 1 : (nsub[j]<csub[j]) ? -1 : 0;
      // chech whether there *could* be a second nearest neighbour one more
      // step in the same direction as the first:
      bool possible = true;
      // first, can we represent nsub+dir as a size_t (will it stay positive)
      for (int j=0; j<3; ++j) possible &= (nsub[j]>0 || dir[j]>=0);
      if (possible){
        for (int j=0; j<3; ++j) nnsub[j] = nsub[j] + dir[j];
        // second, is nnsub a valid mapping index:
        possible = 0 == this->sub2map(nnsub, nn_map_idx);
      }
      // Now, possible tells us if nnindx has been set to a valid mapping index
      // Check whether nnidx points to an already sorted mapped point
      sorted_next_neighbours.push_back(
        (possible && sorted[nn_map_idx]) ? nn_map_idx : n_map_idx
      );
    }
  }
  /* We now have three vectors containing the mapped indices of sorted nearest
     neighbours, and sorted next-nearest neighbours and will return both as a
     single tuple */
  return std::make_tuple(sorted_neighbours, sorted_next_neighbours);
}

template<class T,class R>
std::vector<size_t>
MapGrid3<T,R>::find_unsorted_neighbours(
  const std::vector<bool>& sorted, const size_t centre_map_idx) const {
  size_t centre_lin_idx;
  if (this->map2lin(centre_map_idx,centre_lin_idx))
    throw std::runtime_error("Mapping index has no corresponding linear index.");
  // get a listing of all in-bounds mapped neighbour linear indices
  ArrayVector<size_t> neighbours = this->get_neighbours(centre_lin_idx);
  size_t n_map_idx{0};
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

template<class T,class R>
std::vector<size_t>
MapGrid3<T,R>::find_sorted_neighbours(
  const std::vector<bool>& sorted, const size_t clin) const {
  // get a listing of all in-bounds neighbour linear indexes
  ArrayVector<size_t> neighbours = this->get_neighbours(clin);
  size_t nidx{0}, nnidx{0};
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
  return out;
}



template<class R, class T> template<class S>
bool
MapGrid3<R,T>::sort_difference(const S scaleS,
                               const S scaleV,
                               const S scaleM,
                               const size_t valspn,
                               const size_t vecspn,
                               const size_t nobj,
                               ArrayVector<size_t>& perm,
                               const size_t cidx,
                               const size_t nidx,
                               const int vcf) const {
return jv_permutation(
  data_.values().data().data(cidx, 0), data_.vectors().data().data(cidx, 0),
  data_.values().data().data(nidx, 0), data_.vectors().data().data(nidx, 0),
  data_.values().elements(), data_.vectors().elements(), scaleS, scaleV, scaleM,
  valspn, vecspn, nobj, perm, cidx, nidx, vcf);
}

template<class T,class R> template<class S>
bool
MapGrid3<T,R>::sort_derivative(const S scaleS,
                               const S scaleV,
                               const S scaleM,
                               const size_t valspn,
                               const size_t vecspn,
                               const size_t nobj,
                               ArrayVector<size_t>& perm,
                               const size_t cidx,
                               const size_t nidx,
                               const size_t nnidx,
                               const int vcf) const {
ArrayVector<T> valsorted(data_.values().data().numel(),2u); // sorted neighbours
ArrayVector<R> vecsorted(data_.vectors().data().numel(),2u);
// Copy the data at each point, ensuring that the global permutation for
// nidx and nnidx are respected
size_t nn_i=0;
// std::cout << "Estimate the values of modes at the centre" << std::endl;
for (size_t i=0; i<nobj; ++i){
  bool nn_i_found = false;
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
  for (size_t j=0; j<valspn; ++j){
    valsorted.insert(data_.values().data().getvalue(nidx,    i*valspn+j), 0u, i*valspn+j);
    valsorted.insert(data_.values().data().getvalue(nnidx,nn_i*valspn+j), 1u, i*valspn+j);
  }
  for (size_t j=0; j<vecspn; ++j){
    vecsorted.insert(data_.vectors().data().getvalue(nidx,    i*vecspn+j), 0u, i*vecspn+j);
    vecsorted.insert(data_.vectors().data().getvalue(nnidx,nn_i*vecspn+j), 1u, i*vecspn+j);
  }
}
// std::cout << sorted.to_string();
// calculate the predicted value of each element of the data.
// making the prediction is a special case of linear interpolation:
ArrayVector<T> valpredic(valsorted.numel(),1u); // and prediction
ArrayVector<R> vecpredic(vecsorted.numel(),1u);
size_t cnt = 2u;
size_t corners[2]{0u,1u};
double weights[2]{2,-1};
size_t out_i = 0u;
unsafe_interpolate_to(valsorted,
                      data_.values().elements(),
                      data_.values().branches(),
                      cnt,corners,weights,valpredic,out_i);
unsafe_interpolate_to(vecsorted,
                      data_.vectors().elements(),
                      data_.vectors().branches(),
                      cnt,corners,weights,vecpredic,out_i);
bool rslt;
// std::cout << predic.to_string();
// std::cout << data_.data().to_string(cidx) << std::endl;
// find the assignment of each *predicted* object value to those at cidx:
rslt = jv_permutation(
  data_.values().data().data(cidx,0), data_.vectors().data().data(cidx,0),
  valpredic.data(out_i,0), vecpredic.data(out_i,0), scaleS, scaleV, scaleM,
  valspn, vecspn, nobj, perm, cidx, nidx, vcf
);
// std::cout << ((rslt)?"permutation determined":"permutation failed") << std::endl;
return rslt;
}

template<class R, class T> template<typename S>
ArrayVector<size_t>
MapGrid3<R,T>::centre_sort_perm(const S scalar_weight, const S vector_weight,
                                const S matrix_weight, const int vcf) const {
//
// elshape → (data.size(),Nobj, Na, Nb, ..., Nz) stored at each grid point
// or (data.size(), Na, Nb, ..., Nz) ... but we have properties to help us out!
size_t nobj = data_.branches();
size_t valspn = static_cast<size_t>(data_.values().branch_span());
size_t vecspn = static_cast<size_t>(data_.vectors().branch_span());
// within each index of the data ArrayVector there are spn*nobj entries.

// We will return the permutations of 0:nobj-1 which sort the objects globally
ArrayVector<size_t> perm( nobj, data_.size() );
// We need to keep track of which objects have been sorted thus far
std::vector<bool> sorted(data_.size(), false);

typename CostTraits<T>::type wS, wV, wM;
wS = typename CostTraits<T>::type(scalar_weight);
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
  sorted[n_mapidx] = this->sort_difference(wS,wV,wM, valspn, vecspn,nobj,perm,n_mapidx,mapidx,vcf);
  if (!sorted[n_mapidx]) throw std::runtime_error("Failed to sort Γ neighbour.");
  ++num_sorted;
}
// with the immediate neighbours sorted, start a recursive search (from each to
// be extra sure that we reach everywhere in the Brillouin zone.
std::vector<bool> locked(data_.size());
for (size_t i=0; i<data_.size(); ++i) locked[i] = false;
for (size_t n_mapidx: unsorted_neighbours){
  num_sorted += this->sort_recursion(n_mapidx,wS,wV,wM,vcf,valspn,vecspn,nobj,perm,sorted,locked);
}
if (num_sorted < data_.size())
  throw std::runtime_error("Recursive sorting failed to find all grid points.");
if (num_sorted > data_.size())
  std::cout << "Sorted " << std::to_string(num_sorted) << " grid points but only "
            << std::to_string(data_.size()) << " points are mapped."
            << std::endl;
  // throw std::runtime_error("Recursive sorting visited some grid points more than once.");
return perm;
}

template<class R, class T> template<class S>
size_t
MapGrid3<R,T>::sort_recursion(const size_t centre,
                              const S wS, const S wV, const S wM, const int vcf,
                              const size_t valspn, const size_t vecspn, const size_t nobj,
                              ArrayVector<size_t>& perm,
                              std::vector<bool>& sorted,
                              std::vector<bool>& locked) const {
// std::cout << "Recursive sorting" << std::endl;
std::vector<size_t> unsorted_neighbours = this->find_unsorted_neighbours(sorted, centre);
std::vector<size_t> sorted_neighbours;
size_t nlin, num_sorted=0;
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
  bool success = this->sort_difference(wS, wV, wM, valspn, vecspn, nobj, perm, nmap,
                                       sorted_neighbours[0], vcf);
  // if (sorted_neighbours.size()==1)
  //   success = this->sort_difference(wS, wV, wM, spn, nobj, perm, nmap,
  //                                   sorted_neighbours[0],  vcf);
  // if (sorted_neighbours.size()==2)
  //   success = this->sort_derivative(wS, wV, wM, spn, nobj, perm, nmap,
  //                                   sorted_neighbours[0], sorted_neighbours[1], vcf);

  // if (!success) throw std::runtime_error("Failed to find permutation.");
  // Don't throw here, maybe we can sort from a different direction?
  sorted[nmap] = success;
  if (success) num_sorted++;
  locked[nmap] = false;
}
// All neighbours of map_idx sorted, so now sort all of their neighbours in turn
for (size_t un: unsorted_neighbours)
  num_sorted += this->sort_recursion(un, wS, wV, wM, vcf, valspn, vecspn, nobj, perm, sorted, locked);

return num_sorted;
}



template<class R, class T> template<typename S>
ArrayVector<size_t>
MapGrid3<R,T>::multi_sort_perm(const S scalar_weight, const S vector_weight,
                               const S matrix_weight, const int vcf) const {
//
typename CostTraits<T>::type weights[3];
weights[0] = typename CostTraits<T>::type(scalar_weight);
weights[1] = typename CostTraits<T>::type(vector_weight);
weights[2] = typename CostTraits<T>::type(matrix_weight);
//
int func{vcf};
// elshape → (data.size(),Nobj, Na, Nb, ..., Nz) stored at each grid point
// or (data.size(), Na, Nb, ..., Nz) ... but we have properties to help us out!
size_t nobj = static_cast<size_t>(data_.branches());
size_t valspn = static_cast<size_t>(data_.values().branch_span());
size_t vecspn = static_cast<size_t>(data_.vectors().branch_span());
size_t spobj[3]{valspn, vecspn, nobj};
// within each index of the data ArrayVector there are spn*nobj entries.

// We will return the permutations of 0:nobj-1 which sort the objects globally
ArrayVector<size_t> perm( nobj, data_.size() );
// We need to keep track of which objects have been sorted thus far
std::vector<bool> sorted(data_.size(), false);
std::vector<bool> locked(data_.size(), false);
std::vector<size_t> visited(data_.size(), 0);

size_t mapidx;
size_t subidx[3];
// Start from the Γ point. It should always be a valid mapped point.
// The grid centre along each direction is the ((2N+1)/2 + 1)ᵗʰ entry, but we count from 0 so (2N+1)/2
for (int i=0; i<3; ++i) subidx[i] = this->size(i) >> 1u;
if (this->sub2map(subidx, mapidx))
  throw std::runtime_error("Γ is not a mapped point.");
// The first point gets an arbitrary assignment
for (size_t j=0; j<nobj; ++j) perm.insert(j, mapidx, j);
sorted[mapidx] = true;
++visited[mapidx];
size_t num_sorted = 1;
// kick off the recursive sorting from here:
// num_sorted += this->multi_sort_recursion(mapidx, weights, func, spobj, perm, sorted, locked, visited);
num_sorted += this->multi_sort(mapidx, weights, func, spobj, perm, sorted, locked, visited);
if (num_sorted < data_.size()){
  std::string msg = "Sorting found only "
                  + std::to_string(num_sorted) + " of "
                  + std::to_string(data_.size()) + " grid points.";
  throw std::runtime_error(msg);
}
if (num_sorted > data_.size())
  std::cout << "Sorted " << std::to_string(num_sorted) << " grid points but only "
            << std::to_string(data_.size()) << " points are mapped."
            << std::endl;
bool all_sorted = true;
size_t count_sorted=0, count_unsorted=0;
for (auto i: sorted){
  all_sorted &= i;
  i ? ++count_sorted : ++count_unsorted;
}
if (!all_sorted){
  std::cout << "Only " << count_sorted << " of " << count_sorted+count_unsorted
            << " were sorted." << std::endl;
}
// int ntimes = 1;
// while (!all_sorted && ntimes<100){
//   std::cout << "Not all grid points sorted " << std::to_string(ntimes) << "." <<std::endl;
//   for (size_t i=0; i<sorted.size(); ++i)
//     if (!sorted[i])
//       // num_sorted += this->multi_sort_recursion(i,weights,func,spobj,perm,sorted,locked, visited);
//       num_sorted += this->multi_sort(i,weights,func,spobj,perm,sorted,locked,visited);
//   //
//   all_sorted = true;
//   for (size_t i=0; i<sorted.size(); ++i) all_sorted &= sorted[i];
//   ++ntimes;
// }
return perm;
}

/*! \brief Sort objects on a grid without recursion.

The first point provided should be sorted. This function finds its unsorted
neighbours and adds them to a queue (first-in first-out). For each mapping
index in the queue it finds a consensus sorting permutation and then adds its
not-yet-queued unsorted neighbours to the end of the queue.
If there is no consensus sorting the most frequent sorting is chosen and
neighbours offering a different permutation are flagged as unsorted. This opens
up the possibility of an infinite loop of competing sorting domains so the
number of times a mapped point has been sorted is followed as well and if it
exceeds a threshold value the point is locked to kill the infinite loop.

Two neighbouring points have 4ᴰ⁻¹ shared neighbours in D dimensions (this
isn't right for 1-D, where it is 0, but is right for 2- and 3-D; and counting
overlapping unsorted neighbours in 4-D is beyond me at the moment). So the
length of the queue grows (up to) geometrically. To prevent horrendous memory
use for large grids this method adds only not-yet-queued mapping indexes to the
queue; in one example of a 101×101×3 grid this reduces the maximum queue length
from ~170k to <1k. The uniqueness check is fast but requires a vector of bools
the same size as the mapping indexes which should always be less memory than
the unchecked geometric growth.

If sorting doesn't work well, we might find that a different strategy for moving
hrough the grid might help.
For example, a first-in last-out stack would have the effect of sorting from the
centre to a boundary first rather than all of the centre then outwards in shells.
*/
template<class R, class T> template<class S>
size_t
MapGrid3<R,T>::multi_sort(
  const size_t centre, const S weights[3], const int func, const size_t spobj[3],
  ArrayVector<size_t>& perm, std::vector<bool>& sorted, std::vector<bool>& locked,
  std::vector<size_t>& visited
) const {
  std::vector<size_t> n_idx, nn_idx;
  std::vector<bool> queued;
  queued.reserve(this->maximum_mapping());
  for (size_t i=0; i<this->maximum_mapping(); ++i) queued.push_back(false);
  std::queue<size_t> to_visit;
  to_visit.push(centre);
  size_t num_sorted = 0;
  size_t count=0u, refresh=1u;
  bool more_to_do = true;
  while (more_to_do){
    size_t current = to_visit.front();
    to_visit.pop();
    queued[current] = false;
    if (!locked[current]){
      ++visited[current];
      if (visited[current]>300) locked[current]=true;
      if (!sorted[current]){
        std::tie(n_idx, nn_idx) = this->classify_sorted_neighbours(sorted, current);
        if (n_idx.size()){
          bool success;
          size_t num_derivative = 0;
          if (n_idx.size()==nn_idx.size())
            for (size_t i=0; i<n_idx.size(); ++i)
              if(nn_idx[i]!=n_idx[i]) ++num_derivative;
          if (num_derivative){ // use derivative based sorting whenever possible
            success = this->multi_sort_derivative_all(weights, func, spobj, perm, sorted, current, n_idx, nn_idx/*, num_derivative*/);
          }
          else{
            success = this->multi_sort_difference(weights, func, spobj, perm, sorted, current, n_idx);
          }
          sorted[current] = success;
          if (success) ++num_sorted;
        } else {
          // no sorted neighbours found?!
        }
      }
      if (sorted[current])
        for (size_t i: this->find_unsorted_neighbours(sorted, current))
          if (!queued[i]){
            to_visit.push(i);
            queued[i] = true;
          }
    }
    more_to_do = !to_visit.empty();
    if (++count >= refresh){
      for (int i=0; i<80; i++) std::cout << " ";
      std::cout << "\rPoints to visit: " << to_visit.size();
      count = 0u;
      refresh = to_visit.size() >> 4;
      if (refresh < 1u) refresh = 1u;
      more_to_do = to_visit.size()>0;
    }

  }
  std::cout << std::endl;
  return num_sorted;
}

template<class R, class T> template<class S>
bool
MapGrid3<R,T>::multi_sort_difference(
  const S weights[3], const int func, const size_t spobj[3],
  ArrayVector<size_t>& perm, std::vector<bool>& sorted,
  const size_t cidx, const std::vector<size_t> nidx
) const {
  if (!nidx.size()) return false;
  ArrayVector<size_t> tperm(perm.numel(), nidx.size()+1);
  for (size_t i=0; i<nidx.size(); ++i){
    for (size_t j=0; j<perm.numel(); ++j)
      tperm.insert(perm.getvalue(nidx[i],j),nidx.size(),j);
    jv_permutation(
      data_.values().data().data(cidx,0), data_.vectors().data().data(cidx,0),
      data_.values().data().data(nidx[i],0), data_.vectors().data().data(nidx[i],0),
      data_.values().elements(), data_.vectors().elements(),
      weights[0], weights[1], weights[2], spobj[0], spobj[1], spobj[2],
      tperm, i, nidx.size(), func
    );
  }
  // now tperm has the permutations calculated for the centre position based
  // on each sorted neighbour in turn stored from 0 to nidx.size()-1.
  /* In a perfect world all permutations would agree with each other, but this
  seems unlikely for all cases. So we need to check. */
  std::vector<bool> uncounted(nidx.size(), true);
  std::vector<size_t> frequency(nidx.size(), 0u), equiv_to(nidx.size(),0u);
  bool all_agree;
  for (size_t i=0; i<nidx.size()-1; ++i){
    if (uncounted[i]){
      equiv_to[i] = i;
      for (size_t j=i+1; j<nidx.size(); ++j){
        if (uncounted[j]){
          all_agree = true;
          for (size_t k=0; k<perm.numel(); ++k){
            all_agree &= tperm.getvalue(i,k) == tperm.getvalue(j,k);
            if (!all_agree) break;
          }
          if (all_agree){
            uncounted[j] = false;
            ++frequency[i];
            equiv_to[j] = i;
          }
        }
      }
    }
  }
  size_t hfidx=0, hf=0;
  for (size_t i=0; i<nidx.size(); ++i)
    if (frequency[i]>hf){
      hf = frequency[i];
      hfidx = i;
    }
  // we pick the highest-frequency permutation to be the right one
  for (size_t i=0; i<perm.numel(); ++i) perm.insert(tperm.getvalue(hfidx,i),cidx,i);
  // if the frequency is nidx.size()-1 (since we skipped j==i) then all
  // permutations agree; otherwise we need to reset the sort-flag for some
  // neighbours since we want to be able to interpolate between *any* 2+ points
  // in the grid, so all sorting permutations *must* be in agreement.
  /*
  This seems dangerous as there might be situations where two or more enclaves
  develop and a never-ending skirmish develops on their border(s).
  */
  if (hf < nidx.size()-1){
    for (size_t i=0; i<nidx.size(); ++i){
      if (equiv_to[i]!=hfidx) sorted[nidx[i]] = false;
    }
  }
  return true;
}

template<class R, class T> template<class S>
bool
MapGrid3<R,T>::multi_sort_derivative(
  const S scales[3], const int func, const size_t spobj[3],
  ArrayVector<size_t>& perm, std::vector<bool>& ,//sorted,
  const size_t cidx, const std::vector<size_t> nidx,
  const std::vector<size_t> nnidx, const size_t no_pairs
) const {
ArrayVector<R> valsdat(data_.values().data().numel(),2*no_pairs);
ArrayVector<T> vecsdat(data_.vectors().data().numel(),2*no_pairs); // sorted n and nn
// Copy the data at each point, ensuring that the global permutation for
// nidx and nnidx are respected
size_t nn_i, cnt=0;
bool nn_i_found;
size_t nobj = spobj[2], vecspn = spobj[1], valspn = spobj[0];
std::vector<size_t> n2p_idx;
n2p_idx.reserve(no_pairs);
for (size_t p=0; p<nidx.size(); ++p){
  if (nidx[p]!=nnidx[p]){ // the only indication if a next neighbour exists
    n2p_idx.push_back(p);
    if (n2p_idx.size()<no_pairs){
      for (size_t i=0; i<nobj; ++i){
        nn_i_found = false;
        for (size_t j=0; j<nobj; ++j){
          if (perm.getvalue(nnidx[p],j) == perm.getvalue(nidx[p],i)){
            nn_i = j;
            nn_i_found = true;
          }
        }
        if (!nn_i_found){
          std::string msg = "Next neighbour index not found.\n";
          msg += "perm[n]  = " + perm.to_string(nidx[p]) + "\n";
          msg += "perm[nn] = " + perm.to_string(nnidx[p]);
          throw std::runtime_error(msg);
        }
        for (size_t j=0; j<valspn; ++j){
          valsdat.insert(data_.values().data().getvalue( nidx[p],   i*valspn+j), 2*cnt+0u, i*valspn+j);
          valsdat.insert(data_.values().data().getvalue(nnidx[p],nn_i*valspn+j), 2*cnt+1u, i*valspn+j);
        }
        for (size_t j=0; j<vecspn; ++j){
          vecsdat.insert(data_.vectors().data().getvalue( nidx[p],   i*vecspn+j), 2*cnt+0u, i*vecspn+j);
          vecsdat.insert(data_.vectors().data().getvalue(nnidx[p],nn_i*vecspn+j), 2*cnt+1u, i*vecspn+j);
        }
      }
    }
  }
}
if (n2p_idx.size()!=no_pairs){
  for(auto gst: nidx){ std::cout<<" "<<std::to_string(gst);} std::cout<<std::endl;
  for(auto gst: nnidx){ std::cout<<" "<<std::to_string(gst);} std::cout<<std::endl;
  std::string msg = "Too ";
  msg += ( n2p_idx.size()>no_pairs ? "many" : "few");
  msg += " derivative pairs found!";
  throw std::runtime_error(msg);
}
ArrayVector<T> vecpredic(vecsdat.numel(),nidx.size());
ArrayVector<R> valpredic(valsdat.numel(),nidx.size());
cnt = 2u;
size_t corners[2]{0u,1u};
double weights[2]{2,-1};
for (size_t i=0; i<no_pairs; ++i){
  corners[0] = 2*i;
  corners[1] = 2*i+1;
  unsafe_interpolate_to(valsdat, data_.values().elements(), data_.values().branches(), cnt, corners, weights, valpredic, i);
  unsafe_interpolate_to(vecsdat, data_.vectors().elements(), data_.vectors().branches(), cnt, corners, weights, vecpredic, i);
}
// The derivative-based permutations have been performed, now determine the
// permutations for each neighbour/next-neighbour pair.
ArrayVector<size_t> tperm(perm.numel(), no_pairs+1);
for (size_t i=0; i<no_pairs; ++i){
  for (size_t j=0; j<perm.numel(); ++j)
    tperm.insert(perm.getvalue(nidx[n2p_idx[i]],j),no_pairs,j);
  jv_permutation(
    data_.values().data().data(cidx,0), data_.vectors().data().data(cidx,0),
    valpredic.data(i,0), vecpredic.data(i,0),
    data_.values().elements(), data_.vectors().elements(),
    scales[0], scales[1], scales[2], spobj[0], spobj[1], spobj[2],
    tperm, i, no_pairs, func
  );
}
// now tperm has the permutations calculated for the centre position based
// on each sorted neighbour in turn stored from 0 to nidx.size()-1.
/* In a perfect world all permutations would agree with each other, but this
seems unlikely for all cases. So we need to check. */
auto uncounted = std::unique_ptr<bool[]>(new bool[no_pairs]);
auto frequency = std::unique_ptr<size_t[]>(new size_t[no_pairs]);
auto equiv_to  = std::unique_ptr<size_t[]>(new size_t[no_pairs]);
bool all_agree;
for (size_t i=0; i<no_pairs; ++i){
  uncounted[i] = true;
  frequency[i] = 0u;
  equiv_to[i] = 0u;
}
for (size_t i=0; i<no_pairs-1; ++i){
  if (uncounted[i]){
    equiv_to[i] = i;
    for (size_t j=i+1; j<no_pairs; ++j){
      if (uncounted[j]){
        all_agree = true;
        for (size_t k=0; k<perm.numel(); ++k){
          all_agree &= tperm.getvalue(i,k) == tperm.getvalue(j,k);
          if (!all_agree) break;
        }
        if (all_agree){
          uncounted[j] = false;
          ++frequency[i];
          equiv_to[j] = i;
        }
      }
    }
  }
}
size_t hfidx=0, hf=0;
for (size_t i=0; i<no_pairs; ++i)
  if (frequency[i]>hf){
    hf = frequency[i];
    hfidx = i;
  }
// we pick the highest-frequency permutation to be the right one
for (size_t i=0; i<perm.numel(); ++i) perm.insert(tperm.getvalue(hfidx,i),cidx,i);
// if (hf < no_pairs-1){
//   for (size_t i=0; i<no_pairs; ++i){
//     if (equiv_to[i]!=hfidx) sorted[nidx[n2p_idx[i]]] = false;
//   }
// }
return true;
}

template<class R, class T> template<class S>
bool
MapGrid3<R,T>::multi_sort_derivative_all(
  const S scales[3], const int func, const size_t spobj[3],
  ArrayVector<size_t>& perm, std::vector<bool>& sorted,
  const size_t cidx, const std::vector<size_t> nidx,
  const std::vector<size_t> nnidx //, const size_t no_pairs
) const {
// int out_count=0;
// std::cout<< "multi_sort_derivative " << std::to_string(++out_count) << std::endl;
ArrayVector<R> valsdat(data_.values().data().numel(),2*nidx.size()); // sorted n and nn
ArrayVector<T> vecsdat(data_.vectors().data().numel(),2*nidx.size());
// Copy the data at each point, ensuring that the global permutation for
// nidx and nnidx are respected
size_t nn_i{0}, cnt=0;
bool nn_i_found;
size_t nobj = spobj[2], vecspn = spobj[1], valspn = spobj[0];
// std::cout<< "multi_sort_derivative " << std::to_string(++out_count) << std::endl;
for (size_t p=0; p<nidx.size(); ++p){
  for (size_t i=0; i<nobj; ++i){
    nn_i_found = false;
    for (size_t j=0; j<nobj; ++j){
      if (perm.getvalue(nnidx[p],j) == perm.getvalue(nidx[p],i)){
        nn_i = j;
        nn_i_found = true;
      }
    }
    if (!nn_i_found){
      std::string msg = "Next neighbour index not found.\n";
      msg += "perm[n]  = " + perm.to_string(nidx[p]) + "\n";
      msg += "perm[nn] = " + perm.to_string(nnidx[p]);
      throw std::runtime_error(msg);
    }
    for (size_t j=0; j<valspn; ++j){
      valsdat.insert(data_.values().data().getvalue( nidx[p],   i*valspn+j), 2*cnt+0u, i*valspn+j);
      valsdat.insert(data_.values().data().getvalue(nnidx[p],nn_i*valspn+j), 2*cnt+1u, i*valspn+j);
    }
    for (size_t j=0; j<vecspn; ++j){
      vecsdat.insert(data_.vectors().data().getvalue( nidx[p],   i*vecspn+j), 2*cnt+0u, i*vecspn+j);
      vecsdat.insert(data_.vectors().data().getvalue(nnidx[p],nn_i*vecspn+j), 2*cnt+1u, i*vecspn+j);
    }
  }
}
// std::cout<< "multi_sort_derivative " << std::to_string(++out_count) << std::endl;
ArrayVector<T> vecpredic(vecsdat.numel(),nidx.size());
ArrayVector<R> valpredic(valsdat.numel(),nidx.size());
cnt = 2u;
size_t corners[2]{0u,1u};
double weights[2]{2,-1};
// std::cout<< "multi_sort_derivative " << std::to_string(++out_count) << std::endl;
for (size_t i=0; i<nidx.size(); ++i){
  corners[0] = 2*i;
  corners[1] = 2*i+1;
  unsafe_interpolate_to(valsdat, data_.values().elements(), data_.values().branches(), cnt, corners, weights, valpredic, i);
  unsafe_interpolate_to(vecsdat, data_.vectors().elements(), data_.vectors().branches(), cnt, corners, weights, vecpredic, i);
}
// std::cout<< "multi_sort_derivative " << std::to_string(++out_count) << std::endl;
// The derivative-based permutations have been performed, now determine the
// permutations for each neighbour/next-neighbour pair.
ArrayVector<size_t> tperm(perm.numel(), nidx.size()+1);
for (size_t i=0; i<nidx.size(); ++i){
  for (size_t j=0; j<perm.numel(); ++j)
    tperm.insert(perm.getvalue(nidx[i],j),nidx.size(),j);
  jv_permutation(
    data_.values().data().data(cidx,0), data_.vectors().data().data(cidx,0),
    valpredic.data(i,0), vecpredic.data(i,0),
    data_.values().elements(), data_.vectors().elements(),
    scales[0], scales[1], scales[2], spobj[0], spobj[1], spobj[2],
    tperm, i, nidx.size(), func
  );
}
// std::cout<< "multi_sort_derivative " << std::to_string(++out_count) << std::endl;
// now tperm has the permutations calculated for the centre position based
// on each sorted neighbour in turn stored from 0 to nidx.size()-1.
/* In a perfect world all permutations would agree with each other, but this
seems unlikely for all cases. So we need to check. */
auto uncounted = std::unique_ptr<bool[]>(new bool[nidx.size()]);
auto frequency = std::unique_ptr<size_t[]>(new size_t[nidx.size()]);
auto equiv_to  = std::unique_ptr<size_t[]>(new size_t[nidx.size()]);
bool all_agree;
for (size_t i=0; i<nidx.size(); ++i){
  uncounted[i] = true;
  frequency[i] = 0u;
  equiv_to[i] = 0u;
}
for (size_t i=0; i<nidx.size()-1; ++i){
  if (uncounted[i]){
    equiv_to[i] = i;
    for (size_t j=i+1; j<nidx.size(); ++j){
      if (uncounted[j]){
        all_agree = true;
        for (size_t k=0; k<perm.numel(); ++k){
          all_agree &= tperm.getvalue(i,k) == tperm.getvalue(j,k);
          if (!all_agree) break;
        }
        if (all_agree){
          uncounted[j] = false;
          ++frequency[i];
          equiv_to[j] = i;
        }
      }
    }
  }
}
// std::cout<< "multi_sort_derivative " << std::to_string(++out_count) << std::endl;
size_t hfidx=0, hf=0;
for (size_t i=0; i<nidx.size(); ++i)
  if (frequency[i]>hf){
    hf = frequency[i];
    hfidx = i;
  }
// std::cout<< "multi_sort_derivative " << std::to_string(++out_count) << std::endl;
// we pick the highest-frequency permutation to be the right one
for (size_t i=0; i<perm.numel(); ++i) perm.insert(tperm.getvalue(hfidx,i),cidx,i);
// if the frequency is nidx.size()-1 (since we skipped j==i) then all
// permutations agree; otherwise we need to reset the sort-flag for some
// neighbours since we want to be able to interpolate between *any* 2+ points
// in the grid, so all sorting permutations *must* be in agreement.
/*
This seems dangerous as there might be situations where two or more enclaves
develop and a never-ending skirmish develops on their border(s).
*/
// std::cout<< "multi_sort_derivative " << std::to_string(++out_count) << std::endl;
if (hf < nidx.size()-1){
  for (size_t i=0; i<nidx.size(); ++i){
    if (equiv_to[i]!=hfidx) sorted[nidx[i]] = false;
  }
}
return true;
}
