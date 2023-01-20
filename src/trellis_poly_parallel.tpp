/* This file is part of brille.

Copyright © 2019,2020 Greg Tucker <greg.tucker@stfc.ac.uk>

brille is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

brille is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with brille. If not, see <https://www.gnu.org/licenses/>.            */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <cassert>
#endif

#include "is_musl.h"

template<class I>
static void add_to_maps(const I n_points, I& n_kept, std::vector<I>& map_idx, std::vector<I>& map, const I i)
{
  if (map_idx[i] > n_points) {
    map_idx[i] = n_kept++;
    debug_update("vertex ", i, " added to kept list; total kept now ", n_kept);
  }
  map.push_back(map_idx[i]);
}

template<class I, class T, template<class> class A>
static bool add_vertex(const I n_points,
                       I& n_kept,
                       std::vector<I>& map_idx,
                       const A<T>& points,
                       A<T>& extra,
                       I& n_extra,
                       const T s_tol,
                       const int d_tol,
                       const A<T>& vertex,
                       std::vector<I>& map
                       )
{
  auto equals = points.row_is(brille::cmp::eq, vertex, s_tol, s_tol, d_tol);
  auto no = equals.count();
  if (no) {
    if (no > 1) throw std::runtime_error("Too many matches to vertex");
    // modifies n_kept, map_index & map
    add_to_maps(n_points, n_kept, map_idx, map, equals.first());
    return false;
  }
  equals = extra.row_is(brille::cmp::eq, vertex, s_tol, s_tol, d_tol);
  no = equals.count(n_extra); // only count set points (protect against matching an uninitialized extra entry)
  if (no) {
    if (no > 1) throw std::runtime_error("Too many extra matches to vertex");
    map.push_back(n_points + equals.first(n_extra));
    return false;
  }
  if (extra.size(0) < n_extra + 1)
    extra.resize(2 * n_extra);
  extra.set(n_extra, vertex);
  map.push_back(n_points + n_extra++);
  return true;
}


template<class T, class R, class S, template<class> class A>
void PolyTrellis<T,R,S,A>::construct(const polyhedron::Poly<S,A>& poly,
                                  const double max_volume,
                                  const bool always_triangulate,
                                  const approx_t cfg)
{
  profile_update("Start of PolyTrellis construction");
  assert(poly.face_count() > 3 && poly.volume() > 0);

  approx_ = cfg;
  S s_tol = cfg.reciprocal<S>();
  int d_tol = cfg.digit();
  // find the extents of the polyhedron
  const auto v_hkl{poly.vertices()};
  const auto pv = get_xyz(poly.vertices());
  auto min = pv.min(0).to_std();
  auto max = pv.max(0).to_std();
  auto lengths = pv.max(0) - pv.min(0);
  lengths /= (lengths / std::cbrt(max_volume)).ceil();
  auto nl = lengths.to_std();

  // build-up the trellis intersection-point knots:
  for (int i = 0; i < 3; ++i) {
    knots_[i].reserve(static_cast<size_t>(std::ceil((max[i] - min[i]) / nl[i]) + 1));
    knots_[i].push_back(min[i]);
    while (knots_[i].back() < max[i])
      knots_[i].push_back(knots_[i].back() + nl[i]);
    debug_update("PolyTrellis has ", knots_[i].size() - 1, " bins along axis ",
                 i, ", with boundaries ", knots_[i]);
  }

  auto node_centres =
      from_xyz_like(v_hkl, from_std_like(pv, this->trellis_centres()));
  double max_dist =
      this->maximum_node_circumsphere_radius() + poly.circumsphere_radius();

  std::vector<NodeType> node_type;
  auto is_null = norm(node_centres - poly.centroid()).is(brille::cmp::gt, max_dist).to_std();
  std::transform(is_null.begin(), is_null.end(), std::back_inserter(node_type),
               [](const auto & b){return b ? NodeType::assumed_null : NodeType::cube;});

  // pull together the trellis knots
  auto all_points = from_xyz_like(v_hkl, from_std_like(pv, this->trellis_intersections()));

//  VertexMapSet vertex_set(all_points, s_tol, d_tol);

  // Go through all nodes and determine if they are null, cube, or poly
  // For cube and poly nodes, map their vertices to those in the knots
  // keeping track of 'extra' non-knot polyhedron vertices, and their number
  // the number of knots used in any node vertex mappings
  auto [stash, vertex_set, node_index_map] = part_one(poly, all_points, node_type, always_triangulate, s_tol, d_tol);

  // find the bounding polyhedra vertices in the knots or extra points
  std::vector<std::pair<MapVertexType, ind_t>> boundary_map;
  boundary_map.reserve(v_hkl.size(0));
  for (ind_t i = 0; i < v_hkl.size(0); ++i) {
    boundary_map.push_back(vertex_set.add(v_hkl.view(i)));
  }

  //auto consolidated_vertex_set = vertex_set.consolidate();
  vertices_ = vertex_set.consolidate().pristine();

  auto n_kept = vertex_set.preserved_count();
//  auto n_lost = vertex_set.pristine_count() - n_kept;

  // allocates and fills the NodeContainer
  part_two(stash, node_type, n_kept, node_index_map, s_tol, d_tol);
  // Now all non-null nodes have been populated with the indices of their vertices

  // create the Faces object too:
  // update the boundary map to account for point extraction:
  std::vector<ind_t> update_boundary_map;
  for (auto x: boundary_map){
//    update_boundary_map.push_back(idx < n_kept ? idx : idx - n_lost);
    update_boundary_map.push_back(x.first == MapVertexType::Pristine ? x.second : n_kept + x.second);
  }
  auto p_faces = poly.faces(); // Faces(std::vector<std::vector>>)
  auto pf_faces = p_faces.faces(); // std::vector<std::vector>>
  for (auto & face: pf_faces) for (auto & index: face) index = update_boundary_map[index];
  boundary_ = polyhedron::Faces(pf_faces);


  profile_update("  End of PolyTrellis Construction");
  // the data_ PermutationTable should be initialised now:
  data_.initialize_permutation_table(vertices_.size(0), this->collect_keys());
}

template<class T, class R, class S, template<class> class A>
std::tuple<
    std::map<size_t, typename PolyTrellis<T,R,S,A>::poly_t>,
    VertexMapSet<S,A>, VertexIndexMap
    >
PolyTrellis<T,R,S,A>::part_one(const poly_t& poly, const A<S>& all_points, std::vector<NodeType>& node_type, const bool always_triangulate, const S s_tol, const int d_tol) {
#ifdef __MUSL__
  info_update("musl libc and OpenMP cause a segmentation violation in tests -- forcing single-threaded triangulation");
  omp_set_num_threads(1);
#endif
  profile_update("Starting PolyTrellis part_one");
  //  info_update("Using tolerance ", s_tol, " and digits ", d_tol);
  ind_t nNodes = this->node_count();
  /*
  Find intersections corresponding to nodes that intersect with the polyhedron:
  The original approach was too simplistic as it *only* considered whether an
  intersection point is *inside* of the polyhedron. Such a criteria will not
  capture nodes where there is overlap but no vertices from the node are within
  the polyhedron.
  */
  std::vector<std::pair<VertexMapSet<S,A>, VertexIndexMap>> thread_pairs;

  auto Gamma = 0 * all_points.view(0);

//  VertexMapSet vertex_set(all_points, s_tol, d_tol);
//  auto vertex_index_map = VertexIndexMap();
  std::map<size_t, poly_t> poly_stash;
  std::mutex stash_mutex;

  debug_update_if(
    std::find(node_type.begin(), node_type.end(), NodeType::assumed_null) == node_type.end(),
    "No null nodes!");

  // Collect the cube node indexes to improve parallel loop work distribution
  std::vector<ind_t> cube_indexes;
  cube_indexes.reserve(nNodes);
  for (ind_t i=0; i < nNodes; ++i) if (NodeType::cube == node_type[i]) cube_indexes.push_back(i);

//#pragma omp parallel master default(none) shared(all_points, s_tol, d_tol, thread_pairs)
#pragma omp parallel default(none) shared(all_points, s_tol, d_tol, thread_pairs)
  {
    // Perform per-thread initialization *inside* a parallel region, because
    // we don't know how many threads the pool will have outside.
    // Annoyingly this can't be done in the following parallel region
    // because the non-master threads pass by a master region as if it isn't
    // there, so they may not wait for their storage to actually be initialized!
#pragma omp master
    {
    int n_threads = omp_get_num_threads();
    for (int i = 0; i < n_threads; ++i) {
      thread_pairs.push_back(std::make_pair(VertexMapSet<S,A>(all_points, s_tol, d_tol),
                              VertexIndexMap()));
      }
    }
  }

  auto test0 = !always_triangulate;
  auto s_count = utils::u2s<long long>(cube_indexes.size());
  ThreadException thread_ex;
#pragma omp parallel default(none)\
                     shared(poly, poly_stash, stash_mutex, test0, node_type,\
                     s_count, cube_indexes, s_tol, d_tol,\
                     all_points, thread_pairs, Gamma, thread_ex)
  {
    // stash our thread number to access this thread's pair
    int thread = omp_get_thread_num();
#pragma omp for
    for (long long s_i = 0; s_i < s_count; ++s_i) {
      auto i = cube_indexes[utils::s2u<ind_t>(s_i)];
      auto this_node_faces = trellis_node_faces(i);
      // this limits all_points to have the knots *first*
      auto this_node_poly = polyhedron::Poly(all_points, this_node_faces);
      auto intersection = poly.intersection(this_node_poly, s_tol, d_tol);
      /* FIXME A less-strict comparison might be useful but, at present, causes
       *       not-in-trellis errors */
      auto test1 = this_node_poly.volume() == intersection.volume();
      auto test2 = !this_node_poly.contains(Gamma)[0];
      node_type[i] = (test0 && test1 && test2) ? NodeType::cube : NodeType::poly;
      auto indexes = this_node_faces.indexes();
      if (NodeType::cube == node_type[i]) {
        auto & node_map = thread_pairs[thread].second.get(i);
        node_map.resize(indexes.size());
        for (ind_t j = 0; j < indexes.size(); j++)
          node_map[j] = thread_pairs[thread].first.preserve(indexes[j]);
      } else {
        if (intersection.face_count() < 4 || intersection.volume() <= 0.) {
          node_type[i] = NodeType::found_null;
          continue; // we don't want to re-wrap the following code block
        }
        // find which of the vertices of the intersection are knots as well.
        const auto &iv{intersection.vertices()};

        // Determine which of the intersection vertices are from the pristine set
        std::map<ind_t, ind_t> i2p;
        for (auto idx : indexes) {
          if (auto ip = iv.row_is(brille::cmp::eq, all_points.view(idx), s_tol, s_tol, d_tol); ip.count() == 1) {
            i2p.emplace(ip.first(), idx);
          }
        }
        // *SUPER IMPORTANT* Go through the vertices of the intersection IN
        // ORDER Check for the index in the already-constructed map
        auto & node_map = thread_pairs[thread].second.get(i);
        node_map.resize(iv.size(0));
        for (ind_t iv_index = 0; iv_index < iv.size(0); ++iv_index) {
          thread_ex.run([&]{
          if (auto search = i2p.find(iv_index); search != i2p.end()) {
            node_map[iv_index] = thread_pairs[thread].first.preserve(search->second);
          } else {
					  // might-throw:
            node_map[iv_index] = thread_pairs[thread].first.add(iv.view(iv_index), AddVertexType::Crafted);
          }
					}); // end might-throw handler ThreadException runner
        }
      }
      // this was protected by a not-null check; but now if the node is null
      // there will have been a runtime error raised above, so there's no need
      // for a guard
      auto gamma_in = intersection.contains(Gamma);
      if (std::count(gamma_in.begin(), gamma_in.end(), true)) {
        auto int_vertices_are_gamma = intersection.vertices().row_is(brille::cmp::eq, Gamma, s_tol, s_tol, d_tol);
        auto no = int_vertices_are_gamma.count();
        // if the gamma point *is* an intersection point
        // it is already in the node_index map.
        if (no == 0) thread_pairs[thread].second.append(i, thread_pairs[thread].first.origin_index());
      }
      // Storing the intersection into the std::map is likely not thread safe
      stash_mutex.lock();
      poly_stash.emplace(i, intersection);
      stash_mutex.unlock();
    } // end parallel for-loop
  } // end parallel region -- back to single-thread execution
  thread_ex.rethrow(); //Now handle any encountered errors if they exist

  /* Now combine the per-thread VertexMapSets and VertexIndexMaps */
//  for (const auto & x: thread_pairs) std::cout << x.first << "\n" << x.second;

  profile_update(" Start vertex maps reduction");
  auto comb = vertex_maps::parallel_reduce(thread_pairs);
  profile_update("  End of PolyTrellis part_one");
//  std::cout << comb.second;
  return std::make_tuple(poly_stash, comb.first, comb.second);
}

template<class T, class R, class S, template<class> class A>
void
PolyTrellis<T,R,S,A>::part_two(
    const std::map<size_t, poly_t>& poly_stash,
    const std::vector<NodeType>& node_type,
    const ind_t n_kept,
    const VertexIndexMap& node_index_map,
    const S s_tol,
    const int d_tol
) {
#ifdef __MUSL__
  info_update("musl libc and OpenMP cause a segmentation violation in tests -- forcing single-threaded triangulation");
  omp_set_num_threads(1);
#endif
  profile_update("Starting PolyTrellis part_two");
  // go through all cells again and construct the actual nodes:

  /* You might be tempted to do this in parallel, but the way NodeContainer
   * is implemented makes this impossible at present.
   *
   * To setup the NodeContainer for parallelisation we would need:
   *    - the total number of nodes (nNodes) ✅
   *    - the multiplicity of each node type (count using node_type) ✅
   *    - a pre-filling constructor for NodeType ✅
   *      + this would take the counts of each type and pre-allocate their
   *        (pointer|reference) storage ✅
   *      + it could take node_type as input and even fill-in the nodes_ array!✅
   *    - a thread-safe setter ❓ Maybe the implementation is ok
   * */
  // pre-allocate the NodeContainer storage and fill-in its mapping vector
  nodes_ = NodeContainer(node_type);
  auto count = this->node_count();

  // collect the cube and polyhedron node indexes to simplify the parallel loop
  // Collect the cube node indexes to improve parallel loop work distribution
  std::vector<ind_t> cube_indexes, poly_indexes;
  cube_indexes.reserve(count);
  poly_indexes.reserve(count);
  for (ind_t i=0; i < count; ++i){
    switch (node_type[i]){
    case NodeType::cube: cube_indexes.push_back(i); break;
    case NodeType::poly: poly_indexes.push_back(i); break;
    default: continue;
    }
  }

  profile_update("Cube and Poly node indexes collected");

  // Handle all cubic nodes (this probably doesn't need a parallel loop)
  auto c_count = utils::u2s<long long>(cube_indexes.size());
#pragma omp parallel for default(none) shared(n_kept, poly_stash, cube_indexes, c_count, node_index_map)
  for (long long s_i=0; s_i<c_count; ++s_i) {
    auto i = static_cast<ind_t>(cube_indexes[utils::s2u<ind_t>(s_i)]);
    // we already ensured we don't need to worry about Γ if this is a cube
    std::array<ind_t,8> fvi;
    for (ind_t j=0; j<8u; ++j) fvi[j] = node_index_map.decode(n_kept, i, j);
    // FIXME trellis_node_faces and CubeNode do not agree on vertex ordering!
    // CubeNode requires: [0,0,0], [1,0,0], [1,1,0], [0,1,0], [1,0,1], [0,0,1], [0,1,1], [1,1,1]
    //     faces assumes: [0,0,0], [0,1,0], [0,1,1], [0,0,1], [1,0,0], [1,1,0], [1,1,1], [1,0,1]
    // The Faces object has its vertices ordered to form a valid cube polyhedron via
    //    [3,0,4,7], [3,2,1,0], [0,1,5,4], [3,7,6,2], [7,4,5,6], [2,6,5,1]
    // and then the unique indexes are extracted in order [3,0,4,7,2,1,5,6] before being mapped into node_index_map
    // In the CubeNode vertex order indexing scheme this is [5,0,1,4,6,3,2,7] which requires the permutation
    // [1,2,6,5,3,0,4,7] to put the vertices in the 'correct' order:
    std::array<ind_t,8> cube_vert_idx{{fvi[1], fvi[2], fvi[6], fvi[5], fvi[3], fvi[0], fvi[4], fvi[7]}};
    nodes_.set(i, CubeNode(cube_vert_idx));
  }

  profile_update("Cube node vertices stored");

  // Handle all polyhedron nodes (this almost certainly needs to be parallel)
  ind_t fatal_tri{0}, fatal_miss{0}, fatal_match{0}, hiccups{0};
  auto p_count = utils::u2s<long long>(poly_indexes.size());
#pragma omp parallel for default(none) shared(std::cout, n_kept, poly_stash, poly_indexes, p_count, node_index_map, s_tol, d_tol) reduction(+:fatal_tri, fatal_miss, fatal_match, hiccups)
  for (long long s_i=0; s_i<p_count; ++s_i) {
    auto i = poly_indexes[utils::s2u<ind_t>(s_i)];
    //
    auto Gamma = 0 * vertices_.view(0);
    // combine the vertex indexes for preserved (kept pristine) and appended
    auto node_verts = vertices_.extract(node_index_map.decode(n_kept, i));

    // get a reference the stashed Polyhedron from before:
    const auto & this_node{poly_stash.at(i)};
    auto gin = this_node.contains(Gamma);
    bool contains_Gamma = std::count(gin.begin(), gin.end(), true) == 1;
    // Triangulate the node polyhedron into tetrahedra (the class name LQPolyTet is misleading...)
    polyhedron::LQPolyTet<S,A> tri_cut{};
    // this uses modified TetGen, which is now thread safe :)
    tri_cut = polyhedron::LQPolyTet(this_node, contains_Gamma);
    if (tri_cut.get_vertices().size(0)<4){
      //something went wrong.
      /* A (somehow) likely culprit is that a face is missing from the cut
      cube and therefor is not a piecewise linear complex. try to re-form
      the input polyhedron and then re-triangulate.*/
      tri_cut = polyhedron::LQPolyTet(this_node.convex_hull(), contains_Gamma);
      if (tri_cut.get_vertices().size(0)<4)
        fatal_tri += 1;
      else
        hiccups += 1;
    }
    int added_in_triangulation = tri_cut.any() ? tri_cut.count() : 0;
    // make sure we can match-up the triangulated polyhedron vertices to the
    // known ones:
    std::vector<ind_t> local_map;
    const auto& triverts{tri_cut.get_vertices()};
    for (ind_t j=0; j<triverts.size(0); ++j){
      const auto trij{triverts.view(j)};
      auto known = node_verts.row_is(brille::cmp::eq, trij, s_tol, s_tol, d_tol);
      auto no = known.count();
      if (no){
        if (no > 1){
          std::cout << trij << " matches\n" << cat(1, known, node_verts) << no << " times?!" << std::endl;
//          for (const auto & x: node_index_map.get(i)) std::cout << " " << x.first << x.second;
          for (const auto & x: node_index_map.get(i)) std::cout << " " << x.second;
          std::cout << std::endl;
          fatal_match += 1;
        }
        //local_map.push_back(poly_vert_idx[known.first()]);
        local_map.push_back(node_index_map.decode(n_kept, i, known.first()));
      } else {
        auto punt = vertices_.row_is(brille::cmp::eq, trij, s_tol, s_tol, d_tol);
        if (punt.count() < 1 && added_in_triangulation > 0){
          // we've *ADDED* a new point to the triangulation, so we need to
          // add it to all points as well. hopefully this doesn't happen often
          local_map.push_back(vertices_.size(0));
          vertices_ = cat(0, vertices_, trij);
          --added_in_triangulation;
        } else {
          if (punt.count() < 1) {
            fatal_miss += 1;
          }
          local_map.push_back(punt.first());
        }
      }
    }

    // find the indices per tetrahedron, the tetrahedron circumsphere info,
    // and the tetrahedron volumes
    std::vector<std::array<ind_t,4>> idx_per_tet;
    const auto& local_ipt{tri_cut.get_vertices_per_tetrahedron()};
    for (ind_t j=0; j<local_ipt.size(0); ++j){
      std::array<ind_t,4> one_tet{0,0,0,0};
      for (ind_t k=0; k<4; ++k) one_tet[k] = local_map[local_ipt.val(j,k)];
      idx_per_tet.push_back(one_tet);
    }
    std::vector<std::array<double,4>> cci_per_tet;
    std::vector<double> vol_per_tet;
    for (ind_t j=0; j<tri_cut.number_of_tetrahedra(); ++j){
      cci_per_tet.push_back(tri_cut.circumsphere_info(j));
      vol_per_tet.push_back(tri_cut.volume(j));
    }
    if (idx_per_tet.size()<1){
      throw std::runtime_error("Triangulated node is actually Null!");
    }
    nodes_.set(i, PolyNode(idx_per_tet, cci_per_tet, vol_per_tet));
  }

  profile_update("Poly node vertexes stored");

  if (fatal_tri + fatal_miss + fatal_match){
    std::stringstream msg;
    if (fatal_tri) msg << fatal_tri << " Error(s) determining cut cube triangulation; ";
    if (fatal_match) msg << fatal_match << " Multiple matches of a triangulated vertex; ";
    if (fatal_miss) msg << fatal_miss << " Missing known vertex for triangulated vertex";
    throw std::runtime_error(msg.str());
  }
  debug_update_if(hiccups, "Bad vertex indexing occurred ", hiccups, " times");
  profile_update("  End of PolyTrellis part_two");
}


template<class T, class R, class S, template<class> class A>
std::set<size_t>
PolyTrellis<T,R,S,A>::collect_keys() {
  profile_update("Start of PolyTrellis permutation key collection");
  std::set<size_t> keys;
  long long nnodes = brille::utils::u2s<long long, size_t>(nodes_.size());
  #pragma omp parallel for default(none) shared(keys, nnodes)
  for (long long sni=0; sni<nnodes; ++sni){
    ind_t ni = brille::utils::s2u<ind_t, long long>(sni);
    std::set<size_t> t = this->collect_keys_node(ni);
    if (t.size()){
      #pragma omp critical
      {
        keys.insert(t.begin(), t.end());
      }
    }
  }
  profile_update("  End of PolyTrellis permutation key collection");
  return keys;
}

template<class T, class R, class S, template<class> class A>
std::set<size_t>
PolyTrellis<T,R,S,A>::collect_keys_node(const ind_t node){
  std::set<size_t> keys;
  if (!nodes_.is_null(node)){
    if (nodes_.is_poly(node)){
      // in the case of a polynode we must exploit the inner connectivity
      std::vector<std::array<ind_t,4>> tets = nodes_.vertices_per_tetrahedron(node);
      for (auto vt: tets){
        std::set<size_t> tmp = permutation_table_keys_from_indicies(vt.begin(), vt.end(), vertices_.size(0));
        keys.insert(tmp.begin(), tmp.end());
      }
    } else {
      std::vector<ind_t> vn = nodes_.vertices(node);
      keys = permutation_table_keys_from_indicies(vn.begin(), vn.end(), vertices_.size(0));
    }
  }
  return keys;
}
