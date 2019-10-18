
template<class T> void Mesh3<T>::check_elements(void){
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
template<class T> int Mesh3<T>::replace_data(const ArrayVector<T>& newdata,
                                             const ArrayVector<size_t>& newshape,
                                             const std::array<unsigned,4>& newelements)
{
  this->data = newdata;
  this->shape = newshape;
  this->elements = newelements;
  this->check_elements();
  return 0;
}
template<class T> int Mesh3<T>::replace_data(const ArrayVector<T>& newdata,
                                             const std::array<unsigned,4>& newelements)
{
  ArrayVector<size_t> shape(1,2);
  shape.insert(newdata.size(),0);
  shape.insert(newdata.numel(),1);
  return this->replace_data(newdata, shape, newelements);
}


// Perform sanity checks before attempting to interpolate
template<typename T> template<typename R> unsigned int Mesh3<T>::check_before_interpolating(const ArrayVector<R>& x) const{
  unsigned int mask = 0u;
  if (this->data.size()==0)
    throw std::runtime_error("The mesh must be filled before interpolating!");
  if (x.numel()!=3u)
    throw std::runtime_error("Mesh3 requires x values which are three-vectors.");
  return mask;
}
//! Perform linear interpolating at the specified points in the mesh's orthonormal frame
template<typename T> template<typename R> ArrayVector<T> Mesh3<T>::interpolate_at(const ArrayVector<R>& x) const{
  this->check_before_interpolating(x);
  ArrayVector<T> out(this->data.numel(), x.size());

  const ArrayVector<double>& positions = this->mesh.get_vertex_positions();
  std::vector<size_t> vertices; // should become corners

  ArrayVector<double> xi(3u, 1u);
  std::vector<double> weights;
  for (size_t i=0; i<x.size(); ++i){
    xi = x.extract(i);
    vertices = this->mesh.locate_for_interpolation(xi);
    weights = tetrahedron_weights(xi, positions.extract(vertices));
    info_update("Now call new_unsafe_interpolate_to");
    new_unsafe_interpolate_to(this->data, this->elements, this->branches, vertices, weights, out, i);
  }
  return out;
}
template<typename T> template<typename R> ArrayVector<T> Mesh3<T>::parallel_interpolate_at(const ArrayVector<R>& x, const int threads) const{
  omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
  this->check_before_interpolating(x);
  // shared between threads
  ArrayVector<T> out(this->data.numel(), x.size());
  const ArrayVector<double>& positions = this->mesh.get_vertex_positions();
  // private to each thread
  ArrayVector<double> xi(3u, 1u);
  std::vector<size_t> indexes;
  std::vector<double> weights;
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  slong xsize = unsigned_to_signed<slong, size_t>(x.size());
#pragma omp parallel for shared(x, out, positions) private(xi, indexes, weights)
  for (slong si=0; si<xsize; ++si){
    size_t i = signed_to_unsigned<size_t, slong>(si);
    xi = x.extract(i);
    indexes = this->mesh.locate_for_interpolation(xi);
    weights = tetrahedron_weights(xi, positions.extract(indexes));
    new_unsafe_interpolate_to(this->data, this->elements, this->branches, indexes, weights, out, i);
  }
  return out;
}


template<class T>
template<class R, class S>
ArrayVector<S> Mesh3<T>::debye_waller_sum(const LQVec<R>& Q, const R t_K) const{
  return this->debye_waller_sum(Q.get_xyz(), t_K);
}

template<class T>
template<class R, class S>
ArrayVector<S> Mesh3<T>::debye_waller_sum(const ArrayVector<R>& Q, const R t_K) const{
  const S hbar = 6.582119569E-13; // meV⋅s
  const S kB   = 8.617333252E-2; // meV⋅K⁻¹
  if (Q.numel() != 3)
    throw std::runtime_error("Debye-Waller factor requires 3-vector Q.");
  if (this->elements[0] != 1u)
    throw std::runtime_error("Debye-Waller factor requires one scalar (energy) per mode.");
  size_t nIons = this->elements[1] / 3u;
  if (0 == nIons || this->elements[2]*3u != nIons)
    throw std::runtime_error("Debye-Waller factor requires 3-vector eigenvector(s).");
  size_t nQ = Q.size();
  ArrayVector<S> WdQ(nIons,nQ); // Wᵈ(Q) has nIons entries per Q point

  S coth_en, Q_dot_e_2;
  size_t span = 1u + nIons*3u + this->elements[2] + this->elements[3]*this->elements[3];
  size_t nq = this->shape.getvalue(0u);

  const S beta = kB*t_K; // meV
  const S pref{hbar*hbar/static_cast<S>(2*nq)}; // meV²⋅s²

  S qj_sum;
  // for each input Q point
  for (size_t Qidx=0; Qidx<nQ; ++Qidx){
    // and each ion
    for (size_t d=0; d<nIons; ++d){
      qj_sum = S(0);
      // sum over all reduced q in the first Brillouin zone
      for (size_t q=0; q<nq; ++q){
        // and over all 3*nIon branches at each q
        for (size_t j=0; j<this->branches; ++j){
          // for each branch energy, find <2nₛ+1>/ħωₛ ≡ coth(2ħωₛβ)/ħωₛ
          coth_en = coth_over_en(this->data.getvalue(q,j*span), beta);
          // and find |Q⋅ϵₛ|². Note: vector_product(x,y) *is* |x⋅y|²
          Q_dot_e_2 = vector_product(3u, Q.data(Qidx), this->data.data(q,j*span+1u+3u*d));
          // adding |Q⋅ϵₛ|²coth(2ħωₛβ)/ħωₛ to the sum over s for [Qidx, d]
          qj_sum += Q_dot_e_2 * coth_en;
        }
      }
      // with the sum over s complete, normalize by ħ²/2 divided by the number
      // of points in the Brillouin zone and store the result at W[Qidx, d];
      WdQ.insert(qj_sum*pref, Qidx, d);
    }
  }
  return WdQ;
}

template<class T>
template<class R, template<class> class A, class S>
ArrayVector<S> Mesh3<T>::debye_waller(const A<R>& Q, const std::vector<R>& M, const R t_K) const{
  size_t nIons = this->elements[1] / 3u;
  if (0 == nIons || this->elements[1]*3u != nIons)
    throw std::runtime_error("Debye-Waller factor requires 3-vector eigenvector(s).");
  if (M.size() != nIons)
    throw std::runtime_error("Debye-Waller factor requires an equal number of ions and masses.");
  ArrayVector<S> WdQ = this->debye_waller_sum(Q, t_K);
  ArrayVector<S> factor(1u, Q.size());
  S d_sum;
  for (size_t Qidx=0; Qidx<Q.size(); ++Qidx){
    d_sum = S(0);
    for (size_t d=0; d<nIons; ++d){
      d_sum += std::exp(WdQ.getvalue(Qidx, d)/M[d]);
    }
    factor.insert(d_sum*d_sum, Qidx);
  }
  return factor;
}

template<class T> template<typename R>
std::vector<size_t> Mesh3<T>::which_neighbours(const std::vector<R>& t, const R value, const size_t v) const{
  std::vector<size_t> out;
  for (size_t n: this->mesh.neighbours(v)) if (t[n] == value) out.push_back(n);
  return out;
}

template<class T> template<typename R>
ArrayVector<size_t> Mesh3<T>::multi_sort_perm(
  const R scalar_weight, const R eigenv_weight, const R vector_weight,
  const R matrix_weight, const int ef, const int vf
) const {
  typename CostTraits<T>::type weights[4];
  weights[0] = typename CostTraits<T>::type(scalar_weight);
  weights[1] = typename CostTraits<T>::type(eigenv_weight);
  weights[2] = typename CostTraits<T>::type(vector_weight);
  weights[3] = typename CostTraits<T>::type(matrix_weight);
  int funcs[2]{ef, vf};
  // elshape → (data.size(),Nobj, Na, Nb, ..., Nz) stored at each grid point
  // or (data.size(), Na, Nb, ..., Nz) ... but we have properties to help us out!
  size_t nobj = this->branches;
  size_t span = static_cast<size_t>(this->elements[0])
              + static_cast<size_t>(this->elements[1])
              + static_cast<size_t>(this->elements[2])
              + static_cast<size_t>(this->elements[3])*static_cast<size_t>(this->elements[3]);
  size_t spobj[2]{span, nobj};
  // within each index of the data ArrayVector there are span*nobj entries.

  // We will return the permutations of 0:nobj-1 which sort the objects globally
  ArrayVector<size_t> perm( nobj, this->data.size() );
  // We need to keep track of which objects have been sorted thus far
  std::vector<bool> sorted(this->data.size());
  std::vector<bool> locked(this->data.size());
  std::vector<size_t> visited(this->data.size());
  for (size_t i=0; i<this->data.size(); ++i){
    sorted[i] = false;
    locked[i] = false;
    visited[i] = 0;
  }
  // Start from the Γ point (which we previously need to ensure is in the mesh)
  ArrayVector<double> Gamma(1u, 3u, 0.);
  std::vector<size_t> verts = this->mesh.locate_for_interpolation(Gamma);
  if (verts.size() != 1)
    throw std::runtime_error("The Gamma points is not a mesh point?!");
  size_t idx = verts[0];
  // assign an arbitrary sorting permulation for the Γ point
  for (size_t j=0; j<nobj; ++j) perm.insert(j, idx, j);
  sorted[idx] = true;
  ++visited[idx];
  size_t num_sorted = 1;
  // Now do all of the sorting:
  num_sorted += this->consensus_sort_from(idx, weights, funcs, spobj, perm, sorted, locked, visited);
  if (num_sorted != this->data.size()){
    std::string msg;
    msg = "Sorting visited " + std::to_string(num_sorted) + " of "
        + std::to_string(this->data.size()) + " mesh points.";
    if (num_sorted < this->data.size()) throw std::runtime_error(msg);
    std::cout << msg << std::endl;
  }
  size_t count_sorted=0;
  for (auto i: sorted) if (i) ++count_sorted;
  if (count_sorted != this->data.size()){
    std::cout << "Successfully sorted " << count_sorted << " of ";
    std::cout << this->data.size() << " mesh points.";
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
template<class T> template<typename R>
size_t Mesh3<T>::consensus_sort_from(
  const size_t s, const R weights[4], const int funcs[2], const size_t spobj[2],
  ArrayVector<size_t>& perm, std::vector<bool>& done, std::vector<bool>& lock,
  std::vector<size_t>& visits
) const {
  std::vector<size_t> neighbours;
  std::vector<bool> queued(this->data.size(), false);
  std::queue<size_t> queue;
  queue.push(s);
  size_t max_visits = 300;
  size_t num_sorted=0, idx=0u, count=0u, refresh=1u;
  bool more_to_do=true;
  while (more_to_do){
    idx = queue.front();
    queue.pop();
    queued[idx] = false;
    if (!lock[idx]){
      ++visits[idx];
      if (visits[idx] > max_visits) lock[idx] = true;
      if (!done[idx]){
        neighbours = this->which_neighbours(done, true, idx); // find sorted neighbours
        done[idx] = this->consensus_sort_difference(weights, funcs, spobj, perm, done, idx, neighbours);
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

template<class T> template<typename R>
bool Mesh3<T>::consensus_sort_difference(
  const R w[4], const int funcs[2], const size_t spobj[2],
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
    jv_permutation(this->data.data(idx, 0),
                   this->data.data(neighbours[i], 0),
                   this->elements, w[0], w[1], w[2], w[3], spobj[0], spobj[1],
                   t, i, nn, funcs[0], funcs[1]);
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
