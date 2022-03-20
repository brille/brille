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

template<class T, class R, class S, template<class> class A>
PolyTrellis<T,R,S,A>::PolyTrellis(const polyhedron::Poly<S,A>& poly, const double max_volume, const bool always_triangulate)
: boundary_(poly.faces()), vertices_(poly.vertices())
{
  profile_update("Start of PolyTrellis construction");
  assert(poly.face_count() > 3 && poly.volume() > 0);
  // find the extents of the polyhedron
  const auto v_hkl{poly.vertices()};
  const auto pv = get_xyz(poly.vertices());
  auto min = pv.min(0).to_std();
  auto max = pv.max(0).to_std();
  auto lengths = pv.max(0) - pv.min(0);
  lengths /= (lengths / std::cbrt(max_volume)).ceil();
  auto nl = lengths.to_std();

  // build-up the trellis intersection-point knots:
  for (int i=0; i<3; ++i){
    knots_[i].reserve(std::ceil((max[i] - min[i]) / nl[i]) + 1);
    knots_[i].push_back(min[i]);
    while (knots_[i].back() < max[i]) knots_[i].push_back(knots_[i].back()+nl[i]);
    debug_update("PolyTrellis has ", knots_[i].size()-1," bins along axis ",i,", with boundaries ",knots_[i]);
  }
  ind_t nNodes = this->node_count();

  auto node_centres = from_xyz_like(v_hkl, from_std_like(pv, this->trellis_centres()));
  double max_dist = this->maximum_node_circumsphere_radius() + poly.circumsphere_radius();
  std::vector<bool> node_is_null = norm(node_centres-poly.centroid()).is(brille::cmp::gt, max_dist).to_std();

  // pull together the trellis knots
  auto all_points = from_xyz_like(v_hkl, from_std_like(pv, this->trellis_intersections()));
  /*
  Find intersections corresponding to nodes that intersect with the polyhedron:
  The original approach was too simplistic as it *only* considered whether an
  intersection point is *inside* of the polyhedron. Such a criteria will not
  capture nodes where there is overlap but no vertices from the node are within
  the polyhedron.
  */
  ind_t n_kept{0}, n_extra{0}, n_points{all_points.size(0)};
  auto extra_points = 0 * all_points.view(0);
  extra_points.resize(n_points >> 1u);
  std::vector<ind_t> map_idx(n_points, n_points+1);
  std::vector<std::vector<ind_t>> node_index_map(n_points);
  std::vector<bool> node_is_cube(nNodes, false);
  auto Gamma = 0 * all_points.view(0);
  std::map<size_t, polyhedron::Poly<S,A>> poly_stash;
  ind_t gamma_at;
  bool gamma_found;
  {
    auto equals = all_points.row_is(brille::cmp::eq, Gamma);
    gamma_found = equals.count() == 1;
    gamma_at = gamma_found ? equals.first() : n_points + 1;
  }

  auto add_to_maps = [&](const ind_t i, std::vector<ind_t> & map){
    if (map_idx[i] > n_points) {
      map_idx[i] = n_kept++;
      debug_update("vertex ", i, " added to kept list; total kept now ", n_kept);
    }
    map.push_back(map_idx[i]);
  };

  auto add_vertex = [&](const auto & vertex, std::vector<ind_t> & map){
    auto equals = all_points.row_is(brille::cmp::eq, vertex);
    auto no = equals.count();
    if (no) {
      if (no>1)
        throw std::runtime_error("Too many matches to vertex");
      add_to_maps(equals.first(), map); // modifies map_index & n_kept too
      return false;
    }
    equals = extra_points.row_is(brille::cmp::eq, vertex);
    no = equals.count(n_extra); // only count set points (protect against matching an uninitialized extra_points entry)
    if (no) {
      if (no>1)
        throw std::runtime_error("Too many extra matches to vertex");
      map.push_back(n_points + equals.first(n_extra));
      return false;
    }
    if (extra_points.size(0) < n_extra + 1) extra_points.resize(2 * n_extra);
    extra_points.set(n_extra, vertex);
    map.push_back(n_points + n_extra++);
    return true;
  };

  debug_update_if(std::find(node_is_null.begin(), node_is_null.end(), true) == node_is_null.end(), "No null nodes!");

  for (ind_t i=0; i<nNodes; ++i) if (!node_is_null[i]) {
    auto this_node_faces = trellis_node_faces(i);
    // this limits all_points to have the knots *first*
    auto this_node_poly = polyhedron::Poly(all_points, this_node_faces);
//    verbose_update("Look for intersection of\n", poly.python_string(), " and node\n", this_node_poly.python_string());
    auto intersection = poly.intersection(this_node_poly);
//    if (!always_triangulate) {
//      node_is_cube[i] = this_node_poly.volume() == intersection.volume() &&
//                        !this_node_poly.contains(Gamma)[0];
//    }
    auto test0 = !always_triangulate;
    /* FIXME A less-strict comparison might be useful but, at present, causes not-in-trellis errors */
    // auto test1 = approx_float::scalar(this_node_poly.volume(), intersection.volume());
    auto test1 = this_node_poly.volume() == intersection.volume();
    auto test2 = !this_node_poly.contains(Gamma)[0];
    verbose_update_if(test1 && !test2, "Node would be cubic but contains Gamma");
    verbose_update_if(test1, "Node is cubic because ", my_to_string(this_node_poly.volume() - intersection.volume()), " is zero");
    node_is_cube[i] = test0 && test1 && test2;
    if (node_is_cube[i]) {
      verbose_update("Node ", i, " is a cube");
      auto indexes = this_node_faces.indexes();
      node_index_map[i].reserve(indexes.size());
      for (auto idx : indexes) add_to_maps(idx, node_index_map[i]); // modifies map_index too
    } else {
      node_is_null[i] = intersection.face_count() < 4 || intersection.volume() <= 0.;
      if (node_is_null[i]) continue; // we don't want to re-wrap the following code block
      debug_update("Intersection of node", i, " is \n", intersection.python_string());
      // find which of the vertices of the intersection are knots as well.
      const auto &iv{intersection.vertices()};
      node_index_map[i].reserve(iv.size(0));
      for (ind_t j=0; j < iv.size(0); ++j) add_vertex(iv.view(j), node_index_map[i]); // might modify map_index too
    }
    // this was protected by a not-null check; but now if the node is null
    // there will have been a runtime error raised above, so there's no need
    // for a guard
    auto gamma_in = intersection.contains(Gamma);
    if (std::count(gamma_in.begin(), gamma_in.end(), true)) {
        verbose_update("Intersection contains the Gamma point");
        auto int_vertices_are_gamma =
            intersection.vertices().row_is(brille::cmp::eq, Gamma);
        auto no = int_vertices_are_gamma.count();
        verbose_update(no == 0, "Gamma point not a vertex of the Intersection");
        if (!gamma_found) {
          verbose_update(
              "Maybe the Gamma point is one of the extra intersections?");
          auto extra_points_equals_gamma =
              extra_points.row_is(brille::cmp::eq, Gamma);
          if (extra_points_equals_gamma.count(n_extra /*count only the added points*/) == 1) {
            gamma_at = extra_points_equals_gamma.first(n_extra) + n_points;
            gamma_found = true;
          }
        }
        // only it's still not found do we need to do anything?
        if (!gamma_found) {
          verbose_update("Gamma not located in knots or extra intersections");
          if (extra_points.size(0) < n_extra + 1)
            extra_points.resize(2 * n_extra);
          extra_points.set(n_extra, Gamma);
          gamma_at = n_points + n_extra++;
          gamma_found = true;
        }
        // if the gamma point *is* an intersection point
        // it is already in the node_index map.
        if (no == 0 && gamma_found) node_index_map[i].push_back(gamma_at);
      }
    poly_stash.emplace(i, intersection);
  }

  // find the bounding polyhedra vertices in the knots or extra points
  std::vector<ind_t> boundary_map;
  boundary_map.reserve(v_hkl.size(0));
  for (ind_t i=0; i < v_hkl.size(0); ++i) add_vertex(v_hkl.view(i), boundary_map); // might modify map_index too

  // we now know which of all_points which are needed by the trellis
  std::vector<size_t> to_extract;
  for (size_t i=0; i<n_kept; ++i){
    size_t j = std::distance(map_idx.begin(),std::find_if(map_idx.begin(), map_idx.end(), [i](const size_t t){return t==i;}));
    debug_update_if(j >= n_points, "Could not find index with mapped value ",i,"?!");
    to_extract.push_back(j);
  }
  verbose_update("   map_idx:",map_idx,"\nto_extract:",to_extract);
  // map_idx [0, 3, 1, x, x, 2, 4, 5, …] extracts [0, 2, 5, 1, 6, 7, …] from all_points
  // combine the retained intersection vertices and the added polynode vertices
  if (n_extra > 0) /* protect against view(0,0) */{
    vertices_ = cat(0, all_points.extract(to_extract), extra_points.view(0,n_extra));
  } else {
    vertices_ = all_points.extract(to_extract);
  }

  // go through all cells again and construct the actual nodes:
  for (size_t i=0; i<nNodes; ++i){
    if (node_is_null[i]){
      nodes_.push_back(NullNode());
    } else if (node_is_cube[i]) {
      // we already ensured we don't need to worry about Γ if node_is_cube is true
      std::array<ind_t,8> cube_vert_idx;
      for (size_t j=0; j<8u; ++j) cube_vert_idx[j] = node_index_map[i][j];
      nodes_.push_back(CubeNode(cube_vert_idx));
    } else {
      std::vector<ind_t> poly_vert_idx;
      for (auto idx: node_index_map[i]) {
        debug_update_if(idx > n_kept && idx + n_kept < n_points, "Node ", i, " has an index, ", idx, ", which is out of range!");
        poly_vert_idx.push_back(idx < n_kept ? idx : idx - n_points + n_kept);
      }
      auto node_verts = vertices_.extract(poly_vert_idx);
      // get the stashed Polyhedron from before:
      auto this_node = poly_stash[i];
      auto gin = this_node.contains(Gamma);
      bool contains_Gamma = std::count(gin.begin(), gin.end(), true) == 1;
      // Triangulate the node polyhedron into tetrahedra (the class name LQPolyTet is misleading...)
      auto tri_cut = polyhedron::LQPolyTet(this_node, contains_Gamma);
      if (tri_cut.get_vertices().size(0)<4){
        //something went wrong.
        /* A (somehow) likely cuplrit is that a face is missing from the cut
        cube and therefor is not a piecewise linear complex. try to re-form
        the input polyhedron and then re-triangulate.*/
        tri_cut = polyhedron::LQPolyTet(this_node.convex_hull(), contains_Gamma);
        if (tri_cut.get_vertices().size(0)<4)
          throw std::runtime_error("Error determining cut cube triangulation");
      }

      // make sure we can match-up the triangulated polyhedron vertices to the
      // known ones:
      std::vector<ind_t> local_map;
      const auto& triverts{tri_cut.get_vertices()};
      for (ind_t j=0; j<triverts.size(0); ++j){
        const auto trij{triverts.view(j)};
        auto known = node_verts.row_is(brille::cmp::eq, trij);
        auto no = known.count();
        if (no){
          if (no > 1){
            info_update("For node ", i, ":");
            info_update("Multiple matches of ", trij.to_string(0), " to known vertices");
            info_update(node_verts.to_string());
            throw std::logic_error("Too many matching triangulated vertex");
          }
          local_map.push_back(poly_vert_idx[known.first()]);
        } else {
          auto punt = vertices_.row_is(brille::cmp::eq, trij);
          if (punt.count() < 1){
            info_update("For node ", i, ":");
            info_update("No match of ", trij.to_string(0), " to known vertices");
            info_update(node_verts.to_string(), "or extra points\n", extra_points.to_string());
            throw std::runtime_error("No match to triangulated vertex");
          }
          local_map.push_back(punt.first());
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
        info_update("Triangulated node produces NullNode for node ", i);
        nodes_.push_back(NullNode());
      } else {
        nodes_.push_back(PolyNode(idx_per_tet, cci_per_tet, vol_per_tet));
      }
    }
  }
  // create the Faces object too:
  // update the boundary map to account for point extraction:
  std::vector<ind_t> update_boundary_map;
  for (auto idx: boundary_map){
    update_boundary_map.push_back(idx < n_kept ? idx : idx - n_points + n_kept);
  }
  auto p_faces = poly.faces(); // Faces(std::vector<std::vector>>)
  auto pf_faces = p_faces.faces(); // std::vector<std::vector>>
  for (auto & face: pf_faces) for (auto & index: face) index = update_boundary_map[index];
  boundary_ = polyhedron::Faces(pf_faces);

  // Now all non-null nodes have been populated with the indices of their vertices
  profile_update("  End of PolyTrellis Construction");
  // the data_ PermutationTable should be initialised now:
  data_.initialize_permutation_table(vertices_.size(0), this->collect_keys());
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
