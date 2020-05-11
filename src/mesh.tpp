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

// Perform sanity checks before attempting to interpolate
template<class T, class S> template<typename R>
unsigned int
Mesh3<T,S>::check_before_interpolating(const ArrayVector<R>& x) const{
  unsigned int mask = 0u;
  if (data_.size()==0)
    throw std::runtime_error("The mesh must be filled before interpolating!");
  if (x.numel()!=3u)
    throw std::runtime_error("Mesh3 requires x values which are three-vectors.");
  return mask;
}
//! Perform linear interpolating at the specified points in the mesh's orthonormal frame
template<class T, class S> template<typename R>
std::tuple<ArrayVector<T>,ArrayVector<S>>
Mesh3<T,S>::interpolate_at(const ArrayVector<R>& x) const{
  this->check_before_interpolating(x);
  ArrayVector<T> vals(data_.values().numel(), x.size());
  ArrayVector<S> vecs(data_.vectors().numel(), x.size());
  std::vector<size_t> vertices;
  std::vector<double> weights;
  size_t found_tet, max_valid_tet = this->mesh.number_of_tetrahedra()-1;
  for (size_t i=0; i<x.size(); ++i){
    verbose_update("Locating ",x.to_string(i));
    found_tet = this->mesh.locate(x.extract(i), vertices, weights);
    debug_update_if(found_tet > max_valid_tet,"Point ",x.to_string(i)," not found in tetrahedra!");
    if (found_tet > max_valid_tet)
      throw std::runtime_error("Point not found in tetrahedral mesh");
    verbose_update("Interpolate between vertices ", vertices," with weights ",weights);
    data_.interpolate_at(vertices, weights, vals, vecs, i);
  }
  return std::make_tuple(vals, vecs);
}
template<class T, class S> template<typename R>
std::tuple<ArrayVector<T>,ArrayVector<S>>
Mesh3<T,S>::parallel_interpolate_at(const ArrayVector<R>& x, const int threads) const{
  omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
  this->check_before_interpolating(x);
  // shared between threads
  ArrayVector<T> vals(data_.values().numel(), x.size());
  ArrayVector<S> vecs(data_.vectors().numel(), x.size());
  // private to each thread
  std::vector<size_t> indexes;
  std::vector<double> weights;
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long xsize = unsigned_to_signed<long, size_t>(x.size());
#pragma omp parallel for default(none) shared(x, vals, vecs, xsize) private(indexes, weights) schedule(dynamic)
  for (long si=0; si<xsize; ++si){
    size_t i = signed_to_unsigned<size_t, long>(si);
    this->mesh.locate(x.extract(i), indexes, weights);
    data_.interpolate_at(indexes, weights, vals, vecs, i);
  }
  return std::make_tuple(vals, vecs);
}

template<class T, class S> template<typename R>
std::vector<size_t>
Mesh3<T,S>::which_neighbours(const std::vector<R>& t, const R value, const size_t v) const{
  std::vector<size_t> out;
  for (size_t n: this->mesh.neighbours(v)) if (t[n] == value) out.push_back(n);
  return out;
}

template<class T, class S> template<typename R>
ArrayVector<size_t>
Mesh3<T,S>::multi_sort_perm(
  const R scalar_weight, const R vector_weight, const R matrix_weight, const int vf
) const {
  typename CostTraits<T>::type weights[3];
  weights[0] = typename CostTraits<T>::type(scalar_weight);
  weights[1] = typename CostTraits<T>::type(vector_weight);
  weights[2] = typename CostTraits<T>::type(matrix_weight);
  // elshape → (data.size(),Nobj, Na, Nb, ..., Nz) stored at each grid point
  // or (data.size(), Na, Nb, ..., Nz) ... but we have properties to help us out!
  size_t nobj = static_cast<size_t>(data_.branches());
  size_t valspan = static_cast<size_t>(data_.values().branch_span());
  size_t vecspan = static_cast<size_t>(data_.values().branch_span());
  size_t spobj[3]{valspan, vecspan, nobj};
  // within each index of the data ArrayVector there are span*nobj entries.

  // We will return the permutations of 0:nobj-1 which sort the objects globally
  ArrayVector<size_t> perm( nobj, data_.size() );
  // We need to keep track of which objects have been sorted thus far
  std::vector<bool> sorted(data_.size(), false);
  std::vector<bool> locked(data_.size(), false);
  std::vector<size_t> visited(data_.size(), 0);
  for (size_t i=0; i<data_.size(); ++i){
    sorted[i] = false;
    locked[i] = false;
    visited[i] = 0;
  }
  // Start from the Γ point (which we previously need to ensure is in the mesh)
  ArrayVector<double> Gamma(1u, 3u, 0.);
  std::vector<size_t> verts;
  this->mesh.locate(Gamma, verts);
  if (verts.size() != 1)
    throw std::runtime_error("The Gamma points is not a mesh point?!");
  size_t idx = verts[0];
  // assign an arbitrary sorting permulation for the Γ point
  for (size_t j=0; j<nobj; ++j) perm.insert(j, idx, j);
  sorted[idx] = true;
  ++visited[idx];
  size_t num_sorted = 1;
  // Now do all of the sorting:
  num_sorted += this->consensus_sort_from(idx, weights, vf, spobj, perm, sorted, locked, visited);
  if (num_sorted != data_.size()){
    std::string msg;
    msg = "Sorting visited " + std::to_string(num_sorted) + " of "
        + std::to_string(data_.size()) + " mesh points.";
    if (num_sorted < data_.size()) throw std::runtime_error(msg);
    std::cout << msg << std::endl;
  }
  size_t count_sorted = std::count(sorted.begin(), sorted.end(), true);
  if (count_sorted != data_.size()){
    std::cout << "Successfully sorted " << count_sorted << " of ";
    std::cout << data_.size() << " mesh points.";
  }
  return perm;
}

/*! \brief Consensus sorting of objects on a relational mesh

Starting from a provided mesh-vertex index which has an arbitrary sorting
permutation for the objects it contains, find all neighbouring vertices and add
them to a unique-element first-in first-out queue.

Dequeue a vertex and find all of its unsorted neighbouring vertices and add them
to the queue. Find all sorted vertices neighbouring the dequeued vertex and use
them to perform a consensus-based sorting of the objects at the dequeued vertex.
If no consensus exists accept the most-popular sorting and requeue all
neighbours offering an alternative permutation.
Move on to the next queued vertex.

In order to shortcut deadlocked situtations where no global consensus can be
achieved, keep track of how many times a vertex has been visited and lock its
sorting permutation after some threshhold to escape the infinite loop.

The unique-element queue has an overhead of either
  1) O(queue.size()) comparisons, or
  2) all_vertices.size() bytes of memory and an O(1) check
per to-be-added vertex. It's not clear which approach is more appropriate.
In either case, this over head is likely much smaller than the memory use of a
freely-growing queue since there is no bound on the number of shared
neighbouring vertices for two connected vertices in a relational mesh.
*/
template<class T, class S> template<typename R>
size_t
Mesh3<T,S>::consensus_sort_from(
  const size_t s, const R weights[3], const int func, const size_t spobj[3],
  ArrayVector<size_t>& perm, std::vector<bool>& done, std::vector<bool>& lock,
  std::vector<size_t>& visits
) const {
  std::vector<size_t> neighbours;
  std::vector<bool> queued(data_.size(), false);
  std::queue<size_t> queue;
  queue.push(s);
  size_t max_visits = 300;
  size_t num_sorted=0, count=0u, refresh=1u;
  bool more_to_do=true;
  while (more_to_do){
    size_t idx = queue.front();
    queue.pop();
    queued[idx] = false;
    if (!lock[idx]){
      ++visits[idx];
      if (visits[idx] > max_visits) lock[idx] = true;
      if (!done[idx]){
        neighbours = this->which_neighbours(done, true, idx); // find sorted neighbours
        done[idx] = this->consensus_sort_difference(weights, func, spobj, perm, done, idx, neighbours);
        // implement derivative-based sorting for mesh vertices :/
        if (done[idx]) ++num_sorted;
      }
      if (done[idx]) for (size_t i: this->which_neighbours(done, false, idx)) if (!queued[i]) {
        queue.push(i);
        queued[i] = true;
      }
    }
    more_to_do = !queue.empty();
    if (++count >= refresh){
      for (int i=0; i<80; ++i) std::cout << " ";
      std::cout << "\rPoints queued: " << queue.size();
      count = 0u;
      refresh = queue.size() >> 4;
      if (refresh < 1u) refresh = 1u;
      more_to_do = queue.size() > 0;
    }
    std::cout << std::endl;
  }
  return num_sorted;
}

template<class T, class S> template<typename R>
bool
Mesh3<T,S>::consensus_sort_difference(
  const R w[3], const int func, const size_t spobj[3],
  ArrayVector<size_t>& p, std::vector<bool>& done, const size_t idx,
  const std::vector<size_t> neighbours
) const {
  if (!neighbours.size()) return false;

  size_t nn = neighbours.size();
  // For each sorted neighbour, find the sorting permutation which sorts the
  // data at idx the same as the data at the neighbour.
  // Store the results into a temporary permutation array, abusing the last
  // element to hold each neighbour's permutation in turn.
  ArrayVector<size_t> t(p.numel(), nn+1);
  for (size_t i=0; i<nn; ++i){
    for (size_t j=0; j<p.numel(); ++j) t.insert(p.getvalue(neighbours[i], j), nn, j);
    jv_permutation(
      data_.values().data().data(idx, 0),           data_.vectors().data().data(idx, 0),
      data_.values().data().data(neighbours[i], 0), data_.vectors().data().data(neighbours[i], 0),
      data_.values().elements(), data_.vectors().elements(), w[0], w[1], w[2], spobj[0], spobj[1], spobj[2],
      t, i, nn, func
    );
  }
  // The first neighbours.size() elements of tperm now contain the permutation
  // for idx determined by the equal-index neighbour.
  std::vector<bool> uncounted(nn, true);
  std::vector<size_t> freq(nn, 1u), equiv_to(nn, 0u);
  for (size_t i=0; i<nn; ++i) equiv_to[i] = i;
  //bool all_agreee;
  for (size_t i=0; i<nn-1; ++i) if (uncounted[i])
  for (size_t j=i+1; j<nn; ++j) if (uncounted[j] && t.vector_approx(i, j)) {
    uncounted[j] = false;
    ++freq[i];
    freq[j] = 0;
    equiv_to[j] = i;
  }
  size_t hfidx = std::distance(freq.begin(), std::max_element(freq.begin(), freq.end()));
  // Pick the highest-frequency permutation as the right one
  for (size_t j=0; j<p.numel(); ++j) p.insert(t.getvalue(hfidx, j), idx, j);
  // If the highest frequency is neighbours.size() then all permutations agree.
  // Otherwise we need to indicate that the disagreeing neighbours are not done.
  if (freq[hfidx] < nn)
  for (size_t i=0; i<nn; ++i) if (equiv_to[i] != hfidx) done[neighbours[i]] = false;
  return true;
}
