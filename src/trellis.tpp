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
#include <cassert>

template<class T, class R>
PolyhedronTrellis<T,R>::PolyhedronTrellis(const Polyhedron& poly, const double max_volume, const bool always_triangulate):
  polyhedron_(poly), vertices_({3,0})
{
  assert(poly.num_vertices() > 3 && poly.get_volume() > 0);
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
  index_t nNodes = this->node_count();

  ArrayVector<double> node_centres(this->trellis_centres());
  double max_dist = this->trellis_node_circumsphere_radius() + poly.get_circumsphere_radius();
  std::vector<bool> node_is_null = norm(node_centres-poly.get_centroid()).is_approx(Comp::gt, max_dist).to_std();

  ArrayVector<double> all_intersections(this->trellis_intersections());
  auto intersections_span = this->trellis_intersections_span();
  auto node_intersections = this->trellis_local_cube_indices();
  /*
  Find intersections corresponding to nodes that intersect with the polyhedron:
  The original approach was too simplistic as it *only* considered whether an
  intersection point is *inside* of the polyhedron. Such a criteria will not
  capture nodes where there is overlap but no vertices from the node are within
  the polyhedron.
  */
  size_t n_kept{0}, n_extra{0}, n_intersections{all_intersections.size()};
  ArrayVector<double> extra_intersections(3u, n_intersections >> 1u);
  std::vector<size_t> map_idx(n_intersections, n_intersections+1);
  std::vector<std::vector<index_t>> node_index_map(n_intersections);
  std::vector<bool> node_is_cube(nNodes, false);
  ArrayVector<double> Gamma(3u, 1u, 0.);
  std::map<size_t, Polyhedron> poly_stash;
  Polyhedron node_zero = this->trellis_local_cube();
  for (index_t i=0; i<nNodes; ++i) if (!node_is_null[i]) {
    std::array<index_t,3> node_ijk = this->idx2sub(i);
    std::vector<index_t> this_node_idx;
    for (auto ni: node_intersections){
      size_t intersection_idx = 0;
      for (int j=0; j<3; ++j)
        intersection_idx += (node_ijk[j]+ni[j])*intersections_span[j];
      this_node_idx.push_back(intersection_idx);
    }
    ArrayVector<double> this_node_int = all_intersections.extract(this_node_idx);

    if (!always_triangulate && poly.contains(this_node_int).all_true()
      && !(this_node_int.min(0).all_approx(Comp::le, 0.) && this_node_int.max(0).all_approx(Comp::ge, 0.)) )
      node_is_cube[i] = true;

    if (node_is_cube[i]){
      debug_update("Node ",i," is a cube");
      for (auto idx: this_node_idx) if (map_idx[idx] > n_intersections)
        map_idx[idx] = n_kept++;
      for (auto idx: this_node_idx) node_index_map[i].push_back(map_idx[idx]);
    } else if (!node_is_null[i]) {
      // Find the intersection of the Node Polyhedron and the input Polyhedron
      Polyhedron this_node = node_zero.translate(node_centres.extract(i)).intersection(poly);
      node_is_null[i] = this_node.get_vertices().size() < 4u;
      if (!node_is_null[i]){
        double this_node_volume = this_node.get_volume();
        node_is_null[i] = this_node_volume < 0. || approx_scalar(this_node_volume, 0.);
      }
      // find which of the vertices of the this_node polyhedron are intersection
      // vertices as well.
      const ArrayVector<double>& this_node_verts{this_node.get_vertices()};
      for (size_t j=0; j<this_node_verts.size(); ++j){
        const ArrayVector<double> tnvj{this_node_verts.extract(j)};
        debug_update("checking vertex ", tnvj.to_string(""));
        auto cube_idx = find(norm(this_node_int-tnvj).is_approx(Comp::eq,0.));
        if (cube_idx.size()>1) throw std::logic_error("Too many matching vertices");
        if (cube_idx.size()==1){
          index_t local_idx = this_node_idx[cube_idx[0]];
          if (map_idx[local_idx] > n_intersections) map_idx[local_idx] = n_kept++;
          node_index_map[i].push_back(map_idx[local_idx]);
        } else {
          auto extra_idx = find(norm(extra_intersections.first(n_extra)-tnvj).is_approx(Comp::eq,0.));
          if (extra_idx.size()>1)
            throw std::logic_error("How does one point match multiple points when all points should be unique?");
          if (extra_idx.size()>0){
            // info_update("Polyhedron vertex in extra_intersections: idx = ",cube_idx);
            node_index_map[i].push_back(static_cast<index_t>(n_intersections + extra_idx[0]));
          } else {
            // info_update("Polyhedron vertex not in kept or extra intersections, so add it.");
            // make sure we have room to store this new intersection
            // if we don't, make room; but since this involves a memory copy make lots of room
            if (extra_intersections.size() < n_extra+1) extra_intersections.resize(2*n_extra);
            // store the extra vertex
            extra_intersections.set(n_extra, tnvj);
            // and its mapping information
            node_index_map[i].push_back(static_cast<index_t>(n_intersections + n_extra));
            n_extra++;
          }
        }
      }
      if (!node_is_null[i] && this_node.contains(Gamma).count_true() > 0){
        auto gamma_idx = find(norm(this_node_int-Gamma).is_approx(Comp::eq, 0.));
        if (gamma_idx.size() == 0) {
          gamma_idx = find(norm(extra_intersections.first(n_extra)-Gamma).is_approx(Comp::eq, 0.));
          if (gamma_idx.size() == 0) {
            if (extra_intersections.size() < n_extra + 1) extra_intersections.resize(2*n_extra);
            extra_intersections.set(n_extra, Gamma);
            node_index_map[i].push_back(static_cast<index_t>(n_intersections + n_extra++));
          }
        }
      }
      if (!node_is_null[i]) poly_stash.emplace(i, this_node);
    }
  }

  // we now know which of all_intersections which are needed by the trellis
  std::vector<size_t> to_extract;
  for (size_t i=0; i<n_kept; ++i){
    size_t j = std::distance(map_idx.begin(),std::find_if(map_idx.begin(), map_idx.end(), [i](const size_t t){return t==i;}));
    debug_update_if(j >= n_intersections, "Could not find index with mapped value ",i,"?!");
    to_extract.push_back(j);
  }
  debug_update("   map_idx:",map_idx,"\nto_extract:",to_extract);
  // map_idx [0, 3, 1, x, x, 2, 4, 5, …] extracts [0, 2, 5, 1, 6, 7, …] from all_intersections
  // combine the retained intersection vertices and the added polynode vertices
  vertices_ = cat(all_intersections.extract(to_extract), extra_intersections.first(n_extra));

  // go through all cells again and construct the actual nodes:
  for (size_t i=0; i<nNodes; ++i){
    if (node_is_null[i]){
      nodes_.push_back(NullNode());
    } else if (node_is_cube[i]) {
      // we already ensured we don't need to worry about Γ if node_is_cube is true
      std::array<index_t,8> cube_vert_idx;
      for (size_t j=0; j<8u; ++j) cube_vert_idx[j] = node_index_map[i][j];
      nodes_.push_back(CubeNode(cube_vert_idx));
    } else {
      std::vector<index_t> poly_vert_idx;
      for (auto idx: node_index_map[i])
        poly_vert_idx.push_back( idx < n_kept ? idx : idx - n_intersections + n_kept);
      auto node_verts = vertices_.extract(poly_vert_idx);
      // get the stashed Polyhedron from before:
      auto this_node = poly_stash[i];
      bool contains_Gamma = this_node.contains(Gamma).count_true() == 1u;
      // Triangulate the node polyhedron into tetrahedra
      SimpleTet tri_cut(this_node, -1., contains_Gamma);
      if (tri_cut.get_vertices().size()<4){
        //something went wrong.
        /* A (somehow) likely cuplrit is that a face is missing from the cut
        cube and therefor is not a piecewise linear complex. try to re-form
        the input polyhedron and then re-triangulate.*/
        tri_cut = SimpleTet(Polyhedron(this_node.get_vertices()), -1, contains_Gamma);
        // tri_cut = SimpleTet(Polyhedron(cbp.get_vertices()), max_volume, contains_Gamma);
        if (tri_cut.get_vertices().size()<4)
          throw std::runtime_error("Error determining cut cube triangulation");
      }
      // make sure we can match-up the triangualted polyhedron vertices to the
      // known ones:
      std::vector<index_t> local_map;
      const ArrayVector<double>& triverts{tri_cut.get_vertices()};
      for (size_t j=0; j<triverts.size(); ++j){
        const ArrayVector<double> trij{triverts.extract(j)};
        auto known_idx = find(norm(node_verts - trij).is_approx(Comp::eq, 0.));
        if (known_idx.size() > 1){
          info_update("For node ",i,":\n","Multiple matches of ",trij.to_string(0)," to known vertices\n",node_verts.to_string());
          throw std::logic_error("Too many matching to triangulated vertex");
        }
        if (known_idx.size() == 1u){
          local_map.push_back(poly_vert_idx[known_idx[0]]);
        } else {
          // check against all vertices. maybe this is the added Gamma point.
          auto punt_idx = find(norm(vertices_ - trij).is_approx(Comp::eq, 0.));
          if (punt_idx.size() < 1){
            info_update("For node ",i,":\n","No match of ",trij.to_string(0),
                        " to known vertices\n",node_verts.to_string(),
                        " or extra_intersections\n",extra_intersections.to_string());
            throw std::runtime_error("No match to triangulated vertex");
          }
          local_map.push_back(punt_idx[0]);
        }
      }
      // find the indices per tetrahedron, the tetrahedron circumsphere info,
      // and the tetrahedron volumes
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
      } else {
        nodes_.push_back(PolyNode(idx_per_tet, cci_per_tet, vol_per_tet));
      }
    }
  }
  // Now all non-null nodes have been populated with the indices of their vertices

  // the InterpolationData PermutationTable should be initialised now:
  data_.initialize_permutation_table(vertices_.size(), this->collect_keys());
}

template<class T,class S>
std::set<size_t> PolyhedronTrellis<T,S>::collect_keys() {
  std::set<size_t> keys;
  long long nnodes = unsigned_to_signed<long long, size_t>(nodes_.size());
  #pragma omp parallel for default(none) shared(keys, nnodes)
  for (long long sni=0; sni<nnodes; ++sni){
    size_t ni = signed_to_unsigned<size_t, long long>(sni);
    std::set<size_t> t = this->collect_keys_node(ni);
    if (t.size()){
      #pragma omp critical
      {
        keys.insert(t.begin(), t.end());
      }
    }
  }
  return keys;
}

template<class T, class S>
std::set<size_t>
PolyhedronTrellis<T,S>::collect_keys_node(const index_t node){
  std::set<size_t> keys;
  if (!nodes_.is_null(node)){
    if (nodes_.is_poly(node)){
      // in the case of a polynode we must exploit the inner connectivity
      std::vector<std::array<index_t,4>> tets = nodes_.vertices_per_tetrahedron(node);
      for (auto vt: tets){
        std::set<size_t> tmp = permutation_table_keys_from_indicies(vt.begin(), vt.end(), vertices_.size());
        keys.insert(tmp.begin(), tmp.end());
      }
    } else {
      std::vector<index_t> vn = nodes_.vertices(node);
      keys = permutation_table_keys_from_indicies(vn.begin(), vn.end(), vertices_.size());
    }
  }
  return keys;
}
