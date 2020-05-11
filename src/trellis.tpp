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

template<class T, class R>
PolyhedronTrellis<T,R>::PolyhedronTrellis(const Polyhedron& poly, const double max_volume):
  polyhedron_(poly), vertices_({3,0})
{
  // find the extents of the polyhedron
  std::array<std::array<double,2>,3> minmax;
  for (int i=0; i<3; ++i){
    minmax[i][0] = (std::numeric_limits<double>::max)();
    minmax[i][1] = std::numeric_limits<double>::lowest();
  }
  const ArrayVector<double>& pv{poly.get_vertices()};
  for (size_t i=0; i<pv.size(); ++i) for (int j=0; j<3; ++j) {
    if (pv.getvalue(i,j) < minmax[j][0]) minmax[j][0] = pv.getvalue(i,j);
    if (pv.getvalue(i,j) > minmax[j][1]) minmax[j][1] = pv.getvalue(i,j);
  }
  // try to make an integer number of nodes fit along each dimension
  // If the Polyhedron does not have a face perpendicular to the given direction
  // this will make no difference for build time.
  double intended_length = std::cbrt(max_volume);
  std::array<double,3> node_length;
  for (int i=0; i<3; ++i){
    double len = minmax[i][1] - minmax[i][0];
    node_length[i] = len/std::ceil(len/intended_length);
  }
  // construct the trellis boundaries:
  for (int i=0; i<3; ++i){
    boundaries_[i].push_back(minmax[i][0]);
    while (boundaries_[i].back() < minmax[i][1])
      boundaries_[i].push_back(boundaries_[i].back()+node_length[i]);
    debug_update("PolyhedronTrellis has ",boundaries_[i].size()-1," bins along axis ",i,", with boundaries ",boundaries_[i]);
  }
  // find which trellis intersections are inside of the polyhedron
  std::vector<std::array<double,3>> va_int;
  for (double z: boundaries_[2]) for (double y: boundaries_[1]) for (double x: boundaries_[0])
    va_int.push_back({x,y,z});
  // the definition of this intersections span needs to match the looping order above!
  std::array<size_t,3> intersections_span{1,boundaries_[0].size(),boundaries_[0].size()*boundaries_[1].size()};
  ArrayVector<double> all_intersections(va_int);
  ArrayVector<bool> are_inside = poly.contains(all_intersections);
  // these will be retained as node vertices
  size_t kept_idx{0}, n_mapped = are_inside.count_true();
  std::vector<size_t> map_idx(all_intersections.size(), n_mapped+1);
  for (size_t i=0; i<all_intersections.size(); ++i) if (are_inside.getvalue(i))
    map_idx[i] = kept_idx++;

  // find the nodes which are fully inside the polyhedron
  /* Each node with linear index idx has a subscripted index (i,j,k)
     and is surrounded by the trellis intersections of boundaries
     (i,j,k) + { (000), (100), (110), (010), (101), (001), (011), (111)};
  */
  index_t nNodes = this->node_count();
  // the order of the cube node intersections is paramount:
  std::vector<std::array<index_t,3>> node_intersections{{{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},{{1,0,1}},{{0,0,1}},{{0,1,1}},{{1,1,1}}};
  std::vector<bool> node_is_cube(nNodes, true), node_is_outside(nNodes, false);
  for (index_t i=0; i<nNodes; ++i){
    std::array<index_t,3> node_ijk = this->idx2sub(i);
    // this node is a cube if all node intersection vertices are mapped
    for (auto ni: node_intersections) if (node_is_cube[i]) {
      size_t intersection_idx = 0;
      for (int j=0; j<3; ++j) intersection_idx += (node_ijk[j]+ni[j])*intersections_span[j];
      node_is_cube[i] = map_idx[intersection_idx] < n_mapped;
    }
    if (!node_is_cube[i]){
      bool any_mapped = false;
      for (auto ni: node_intersections) {
        size_t intersection_idx = 0;
        for (int j=0; j<3; ++j) intersection_idx += (node_ijk[j]+ni[j])*intersections_span[j];
        any_mapped |= map_idx[intersection_idx] < n_mapped;
      }
      node_is_outside[i] = !any_mapped;
    }
  }
  // Pull out the intersection points which we will keep as vertices:
  ArrayVector<double> kept_intersections = all_intersections.extract(are_inside);
  ArrayVector<double> extra_intersections(3u, 3*(all_intersections.size()-kept_intersections.size()));
  index_t nExtra=0;
  // Now actually create the node objects (which are not fully outside)
  for (index_t i=0; i<nNodes; ++i)
    if (node_is_outside[i]) {
      nodes_.push_back(NullNode());
    } else {
    std::array<index_t,3> node_ijk = this->idx2sub(i);
    std::array<index_t,8> vert_idx; // the 8 vertex indices of the cube
    std::vector<index_t> mapped_vert_idx;
    for (int k=0; k<8; ++k){
      size_t int_idx=0;
      for (int j=0; j<3; ++j) int_idx += (node_ijk[j]+node_intersections[k][j])*intersections_span[j];
      vert_idx[k] = static_cast<index_t>(map_idx[int_idx]);
      if (map_idx[int_idx] < n_mapped) mapped_vert_idx.push_back(static_cast<index_t>(map_idx[int_idx]));
    }
    bool contains_Gamma{true};
    for (int j=0; j<3; ++j){
      double tocheck = boundaries_[j][node_ijk[j]  ];
      contains_Gamma &= tocheck < 0. || approx_scalar(tocheck, 0.);
      tocheck = boundaries_[j][node_ijk[j]+1];
      contains_Gamma &= tocheck > 0. || approx_scalar(tocheck, 0.);
    }
    if (node_is_cube[i] && !contains_Gamma) {
      nodes_.push_back(CubeNode(vert_idx));
    } else {
      // This node intersects the polyhedron. First, find the interior part
      std::array<double,3> min_corner, max_corner;
      for (int j=0; j<3; ++j){
        min_corner[j] = boundaries_[j][node_ijk[j]  ];
        max_corner[j] = boundaries_[j][node_ijk[j]+1];
      }
      Polyhedron cube = polyhedron_box(min_corner, max_corner);
      debug_update("Node ",i," has bounding box\n",cube.get_vertices().to_string());
      // the cubic node extends beyond the bounding polyhedron so we must truncate it
      // Polyhedron pbc = poly.intersection(cube);
      Polyhedron cbp = cube.intersection(poly); // <- this one should (probably) always be faster since the cube is smaller
      double cut_volume = cbp.get_volume();
      if (cbp.get_vertices().size() < 4 || cut_volume < 0 || approx_scalar(cut_volume, 0.)){
        // less than four vertices can not be a polyhedron
        // negative volume means something went wrong
        // zero volume means the node is actually outside and can remain null
        nodes_.push_back(NullNode());
        continue;
      }
      if (cut_volume > cube.get_volume())
        throw std::runtime_error("Cutting the node increased its volume?!");
      // cut the larger polyhedron by the smaller one:
      // Then triangulate it into tetrahedra
      SimpleTet tri_cut(cbp, -1., contains_Gamma);
      // SimpleTet tri_cut(cbp, max_volume, contains_Gamma);
      if (tri_cut.get_vertices().size()<4){
        //something went wrong.
        /* A (somehow) likely cuplrit is that a face is missing from the cut
        cube and therefor is not a piecewise linear complex. try to re-form
        the input polyhedron and then re-triangulate.*/
        tri_cut = SimpleTet(Polyhedron(cbp.get_vertices()), -1, contains_Gamma);
        // tri_cut = SimpleTet(Polyhedron(cbp.get_vertices()), max_volume, contains_Gamma);
        if (tri_cut.get_vertices().size()<4)
          throw std::runtime_error("Error determining cut cube triangulation");
      }
      // we need to retain the *additional* vertices from the triangulation
      // and figure-out a local vertex index mapping
      std::vector<index_t> local_map;
      const ArrayVector<double>& triverts{tri_cut.get_vertices()};
      for (size_t j=0; j<triverts.size(); ++j){
        const ArrayVector<double> trij{triverts.extract(j)};
        debug_update("checking vertex ", trij.to_string(""));
        auto cube_idx = find(norm(kept_intersections.extract(mapped_vert_idx)-trij).is_approx(Comp::eq,0.));
        if (cube_idx.size()>1) throw std::logic_error("Too many matching vertices");
        if (cube_idx.size()==1){
          local_map.push_back(mapped_vert_idx[cube_idx[0]]);
        } else {
          auto extra_idx = find(norm(extra_intersections.first(nExtra)-trij).is_approx(Comp::eq,0.));
          if (extra_idx.size()>1)
            throw std::logic_error("How does one point match multiple points when all points should be unique?");
          if (extra_idx.size()>0){
            // info_update("Polyhedron vertex in extra_intersections: idx = ",cube_idx);
            local_map.push_back(static_cast<index_t>(kept_intersections.size() + extra_idx[0]));
          } else {
            // info_update("Polyhedron vertex not in kept or extra intersections, so add it.");
            // make sure we have room to store this new intersection
            // if we don't, make room; but since this involves a memory copy make lots of room
            if (extra_intersections.size() < nExtra+1) extra_intersections.resize(2*nExtra);
            // store the extra vertex
            extra_intersections.set(nExtra, trij);
            // and its mapping information
            local_map.push_back(static_cast<index_t>(kept_intersections.size()) + nExtra);
            ++nExtra;
          }
        }
      }
      std::vector<std::array<index_t,4>> idx_per_tet;
      const ArrayVector<size_t>& local_ipt{tri_cut.get_vertices_per_tetrahedron()};
      for (size_t j=0; j<local_ipt.size(); ++j){
        std::array<index_t,4> one_tet{0,0,0,0};
        for (int k=0; k<4; ++k) one_tet[k] = local_map[local_ipt.getvalue(j,k)];
        idx_per_tet.push_back(one_tet);
      }
      std::vector<std::array<double,4>> cci_per_tet;
      std::vector<double> vol_per_tet;
      for (size_t j=0; j<tri_cut.number_of_tetrahedra(); ++j){
        cci_per_tet.push_back(tri_cut.circumsphere_info(j));
        vol_per_tet.push_back(tri_cut.volume(j));
      }
      if (idx_per_tet.size()<1){
        nodes_.push_back(NullNode());
      }
      else{
        nodes_.push_back(PolyNode(idx_per_tet, cci_per_tet, vol_per_tet));
      }
    }
  }
  // Now all non-null nodes have been populated with the indices of their vertices
  // Combine the retained trellis vertices and the extra triangulated vertices
  vertices_ = cat(kept_intersections, extra_intersections.first(nExtra));
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
template<class T,class S> template<typename R>
size_t PolyhedronTrellis<T,S>::consensus_sort_nodes(
  const index_t first_idx, const R weights[3], const int func, const size_t spobj[3],
  ArrayVector<size_t>& perm,
  std::vector<bool>& node_done,
  std::vector<bool>& node_lock,
  std::vector<size_t>& node_visits,
  std::vector<bool>& vertex_done,
  std::vector<bool>& vertex_lock,
  std::vector<size_t>& vertex_visits
) const {
  std::vector<size_t> neighbours;
  std::vector<bool> queued(nodes_.size(), false);
  std::queue<index_t> queue;
  queue.push(first_idx);
  const size_t max_visits{10};
  size_t num_sorted=0, count=0u, refresh=1u;
  bool more_to_do=true;
  size_t blanks = 0;
  std::string prompt = "Nodes queued: ";
  std::cout << prompt;
  while (more_to_do){
    index_t idx = queue.front();
    queue.pop();
    queued[idx] = false;
    if (!node_lock[idx]){
      node_lock[idx] = ++node_visits[idx] > max_visits;
      if (!node_done[idx]){
        node_done[idx] = this->consensus_sort_node(weights, func, spobj, perm, vertex_done, vertex_lock, vertex_visits, idx);
        if (node_done[idx]) ++num_sorted;
      }
      // enqueue the neighbouring nodes no matter what
      for (index_t n_idx: this->which_node_neighbours(node_done, false, idx))
        if (!queued[n_idx]) {
          queue.push(n_idx);
          queued[n_idx] = true;
        }
      // and then re-add this node if it wasn't finished
      if (!node_done[idx]) {
        queue.push(idx);
        queued[idx] = true;
      }
    }
    more_to_do = !queue.empty();
    if (++count >= refresh){
      for (size_t i=0; i<blanks; ++i) std::cout << "\b";
      size_t new_blanks{1};
      if (queue.size()>0) new_blanks += static_cast<size_t>(std::floor(std::log10(queue.size())));
      if (new_blanks < blanks){
        for (size_t i=0; i<blanks; ++i) std::cout << " ";
        for (size_t i=0; i<blanks; ++i) std::cout << "\b";
      }
      std::cout << queue.size();
      blanks = new_blanks;
      count = 0u;
      refresh = queue.size() >> 4;
      if (refresh < 1u) refresh = 1u;
      more_to_do = queue.size() > 0;
    }
  }
  std::cout << "\r";
  blanks += prompt.size();
  for (size_t i=0; i<blanks; ++i) std::cout << " ";
  std::cout << "\r";
  return num_sorted;
}


template<class T,class S> template<typename R>
bool PolyhedronTrellis<T,S>::consensus_sort_node(
  const R w[3], const int func, const size_t spobj[3],
  ArrayVector<size_t>& p,
  std::vector<bool>& done,
  std::vector<bool>& lock,
  std::vector<size_t>& visits,
  const index_t node
) const {
  const size_t max_visits{10}; // this is bad.
  std::vector<index_t> verts = nodes_.vertices(node);
  for (index_t & i: verts)
  if (!lock[i]){
    if (!done[i]){
      std::vector<index_t> are_done;
      std::copy_if(verts.begin(), verts.end(), std::back_inserter(are_done), [&](index_t i){return done[i];});
      done[i] = this->consensus_sort_difference(w, func, spobj, p, done, i, are_done);
    }
    lock[i] = ++visits[i] > max_visits;
  }
  size_t count_done = std::count_if(verts.begin(), verts.end(), [&](index_t i){return done[i];});
  return count_done == verts.size();
}

template<class T,class S> template<typename R>
bool PolyhedronTrellis<T,S>::consensus_sort_difference(
  const R w[3], const int func, const size_t spobj[3],
  ArrayVector<size_t>& p, std::vector<bool>& done, const index_t idx,
  const std::vector<index_t>& presorted
) const {
  if (!presorted.size()) return false;
  // Find the *closest* point(s) which are sorted
  ArrayVector<double> vectors = vertices_.extract(presorted) - vertices_.extract(idx);
  std::vector<double> lengths;
  for (size_t i=0; i<vectors.size(); ++i) lengths.push_back(vectors.norm(i));
  double min_length = *(std::min_element(lengths.begin(), lengths.end()));
  std::vector<index_t> nearest;
  for(size_t i=0; i<lengths.size(); ++i)
    if (approx_scalar(lengths[i], min_length))
      nearest.push_back(presorted[i]);


  size_t nn = nearest.size();
  // For each sorted neighbour, find the sorting permutation which sorts the
  // data at idx the same as the data at the neighbour.
  // Store the results into a temporary permutation array, abusing the last
  // element to hold each neighbour's permutation in turn.
  ArrayVector<size_t> t(p.numel(), nn+1);
  for (size_t i=0; i<nn; ++i){
    for (size_t j=0; j<p.numel(); ++j) t.insert(p.getvalue(nearest[i], j), nn, j);
    jv_permutation(data_.values().data().data(idx, 0),        data_.vectors().data().data(idx, 0),
                   data_.values().data().data(nearest[i], 0), data_.vectors().data().data(nearest[i], 0),
                   data_.values().elements(), data_.vectors().elements(), w[0], w[1], w[2], spobj[0], spobj[1], spobj[2],
                   t, i, nn, func);
  }
  // The first presorted.size() elements of tperm now contain the permutation
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
  // If the highest frequency is nearest.size() then all permutations agree.
  // Otherwise we need to indicate that the disagreeing nearest are not done.
  if (freq[hfidx] < nn)
  for (size_t i=0; i<nn; ++i) if (equiv_to[i] != hfidx) done[nearest[i]] = false;
  return true;
}
