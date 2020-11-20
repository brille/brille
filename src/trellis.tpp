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
#include <cassert>

template<class T, class R>
PolyhedronTrellis<T,R>::PolyhedronTrellis(const Polyhedron& poly, const double max_volume, const bool always_triangulate)
: polyhedron_(poly), vertices_(0,3)
{
  profile_update("Start of PolyhedronTrellis construction");
  assert(poly.num_vertices() > 3 && poly.get_volume() > 0);
  // find the extents of the polyhedron
  std::array<std::array<double,2>,3> minmax;
  for (int i=0; i<3; ++i){
    minmax[i][0] = (std::numeric_limits<double>::max)();
    minmax[i][1] = std::numeric_limits<double>::lowest();
  }
  const auto& pv{poly.get_vertices()};
  // auto pvsh = pv.shape();
  // pvsh.back() = 0;
  // for (auto x: SubIt(pv.shape(), pvsh)) for (unsigned j=0; j<3u; ++j){
  //   x.back() = j;
  //   if (pv[x] < minmax[j][0]) minmax[j][0] = pv[x];
  //   if (pv[x] > minmax[j][1]) minmax[j][1] = pv[x];
  // }
  for (ind_t i=0; i<pv.size(0); ++i) for (ind_t j=0; j<3; ++j) {
    if (pv.val(i,j) < minmax[j][0]) minmax[j][0] = pv.val(i,j);
    if (pv.val(i,j) > minmax[j][1]) minmax[j][1] = pv.val(i,j);
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
  ind_t nNodes = this->node_count();

  auto node_centres = bArray<double>::from_std(this->trellis_centres());
  double max_dist = this->trellis_node_circumsphere_radius() + poly.get_circumsphere_radius();
  std::vector<bool> node_is_null = norm(node_centres-poly.get_centroid()).is(brille::cmp::gt, max_dist).to_std();

  auto all_intersections = bArray<double>::from_std(this->trellis_intersections());
  auto intersections_span = this->trellis_intersections_span();
  auto node_intersections = this->trellis_local_cube_indices();
  /*
  Find intersections corresponding to nodes that intersect with the polyhedron:
  The original approach was too simplistic as it *only* considered whether an
  intersection point is *inside* of the polyhedron. Such a criteria will not
  capture nodes where there is overlap but no vertices from the node are within
  the polyhedron.
  */
  ind_t n_kept{0}, n_extra{0}, n_intersections{all_intersections.size(0)};
  bArray<double> extra_intersections(n_intersections>>1u, 3u);
  std::vector<ind_t> map_idx(n_intersections, n_intersections+1);
  std::vector<std::vector<ind_t>> node_index_map(n_intersections);
  std::vector<bool> node_is_cube(nNodes, false);
  bArray<double> Gamma(1u, 3u, 0.);
  std::map<size_t, Polyhedron> poly_stash;
  Polyhedron node_zero = this->trellis_local_cube();
  for (ind_t i=0; i<nNodes; ++i) if (!node_is_null[i]) {
    std::array<ind_t,3> node_ijk = this->idx2sub(i);
    std::vector<ind_t> this_node_idx;
    for (auto ni: node_intersections){
      size_t intersection_idx = 0;
      for (int j=0; j<3; ++j)
        intersection_idx += (node_ijk[j]+ni[j])*intersections_span[j];
      this_node_idx.push_back(intersection_idx);
    }
    auto this_node_int = all_intersections.extract(this_node_idx);

    if (!always_triangulate){
      auto in = poly.contains(this_node_int);
      brille::cmp l{brille::cmp::le}, g{brille::cmp::ge};
      if (std::count(in.begin(), in.end(), false) == 0 && !(this_node_int.min(0).all(l,0.) && this_node_int.max(0).all(g,0.)) )
        node_is_cube[i] = true; // node_is_cube used again, so can't be combined
    }

    if (node_is_cube[i]){
      verbose_update("Node ",i," is a cube");
      for (auto idx: this_node_idx) if (map_idx[idx] > n_intersections)
        map_idx[idx] = n_kept++;
      for (auto idx: this_node_idx) node_index_map[i].push_back(map_idx[idx]);
    } else if (!node_is_null[i]) {
      // Find the intersection of the Node Polyhedron and the input Polyhedron
      Polyhedron this_node = node_zero.translate(node_centres.view(i)).intersection(poly);
      node_is_null[i] = this_node.get_vertices().size(0) < 4u;
      if (!node_is_null[i]){
        double this_node_volume = this_node.get_volume();
        node_is_null[i] = this_node_volume < 0. || brille::approx::scalar(this_node_volume, 0.);
      }
      // find which of the vertices of the this_node polyhedron are intersection
      // vertices as well.
      const auto& this_node_verts{this_node.get_vertices()};
      for (size_t j=0; j<this_node_verts.size(0); ++j){
        const auto tnvj{this_node_verts.view(j)};
        verbose_update("checking vertex ", tnvj.to_string());
        auto cube_idx = norm(this_node_int-tnvj).find(brille::cmp::eq,0.);
        if (cube_idx.size()>1) throw std::logic_error("Too many matching vertices");
        if (cube_idx.size()==1){
          verbose_update("Polyhedron vertex in node cube vertices: idx = ",cube_idx);
          ind_t local_idx = this_node_idx[cube_idx[0]];
          if (map_idx[local_idx] > n_intersections) map_idx[local_idx] = n_kept++;
          node_index_map[i].push_back(map_idx[local_idx]);
        } else {
          bool notlocated{true};
          if (n_extra > 0) /*protect against view(0,0)*/{
            auto extra_idx = norm(extra_intersections.view(0,n_extra)-tnvj).find(brille::cmp::eq,0.);
            if (extra_idx.size()>1)
              throw std::logic_error("How does one point match multiple points when all points should be unique?");
            if (extra_idx.size()>0){
              verbose_update("Polyhedron vertex in extra_intersections: idx = ",extra_idx);
              node_index_map[i].push_back(static_cast<ind_t>(n_intersections + extra_idx[0]));
              notlocated = false;
            }
          }
          if (notlocated){
            verbose_update("Polyhedron vertex not in node or extra intersections, so add it.");
            // make sure we have room to store this new intersection
            // if we don't, make room; but since this involves a memory copy make lots of room
            if (extra_intersections.size(0) < n_extra+1) extra_intersections.resize(2*n_extra);
            // store the extra vertex
            extra_intersections.set(n_extra, tnvj);
            // and its mapping information
            node_index_map[i].push_back(static_cast<ind_t>(n_intersections + n_extra));
            n_extra++;
          }
        }
      }
      if (!node_is_null[i]){
        auto gin = this_node.contains(Gamma);
        if (std::count(gin.begin(), gin.end(), true)>0){
          verbose_update("Non-null node contains Gamma");
          auto gamma_idx = norm(this_node_int).find(brille::cmp::eq, 0.);
          if (gamma_idx.size() == 0) {
            verbose_update("Gamma not located in node cube vertices");
            bool gamma_not_found{true};
            if (n_extra > 0) /*protect against view(0,0)*/ {
              gamma_idx = norm(extra_intersections.view(0,n_extra)).find(brille::cmp::eq, 0.);
              gamma_not_found = gamma_idx.size() == 0;
            }
            if (gamma_not_found) {
              verbose_update("Gamma not located in extra intersections");
              if (extra_intersections.size(0) < n_extra + 1) extra_intersections.resize(2*n_extra);
              extra_intersections.set(n_extra, Gamma);
              node_index_map[i].push_back(static_cast<ind_t>(n_intersections + n_extra++));
            }
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
  verbose_update("   map_idx:",map_idx,"\nto_extract:",to_extract);
  // map_idx [0, 3, 1, x, x, 2, 4, 5, …] extracts [0, 2, 5, 1, 6, 7, …] from all_intersections
  // combine the retained intersection vertices and the added polynode vertices
  if (n_extra > 0) /* protect against view(0,0) */{
    vertices_ = cat(0,all_intersections.extract(to_extract), extra_intersections.view(0,n_extra));
  } else {
    vertices_ = all_intersections.extract(to_extract);
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
      for (auto idx: node_index_map[i])
        poly_vert_idx.push_back( idx < n_kept ? idx : idx - n_intersections + n_kept);
      auto node_verts = vertices_.extract(poly_vert_idx);
      // get the stashed Polyhedron from before:
      auto this_node = poly_stash[i];
      auto gin = this_node.contains(Gamma);
      bool contains_Gamma = std::count(gin.begin(), gin.end(), true) == 1;
      // Triangulate the node polyhedron into tetrahedra
      SimpleTet tri_cut(this_node, -1., contains_Gamma);
      if (tri_cut.get_vertices().size(0)<4){
        //something went wrong.
        /* A (somehow) likely cuplrit is that a face is missing from the cut
        cube and therefor is not a piecewise linear complex. try to re-form
        the input polyhedron and then re-triangulate.*/
        tri_cut = SimpleTet(Polyhedron(this_node.get_vertices()), -1, contains_Gamma);
        // tri_cut = SimpleTet(Polyhedron(cbp.get_vertices()), max_volume, contains_Gamma);
        if (tri_cut.get_vertices().size(0)<4)
          throw std::runtime_error("Error determining cut cube triangulation");
      }
      // make sure we can match-up the triangualted polyhedron vertices to the
      // known ones:
      std::vector<ind_t> local_map;
      const auto& triverts{tri_cut.get_vertices()};
      for (ind_t j=0; j<triverts.size(0); ++j){
        const auto trij{triverts.view(j)};
        auto known_idx = norm(node_verts - trij).find(brille::cmp::eq, 0.);
        if (known_idx.size() > 1){
          info_update("For node ",i,":\n","Multiple matches of ",trij.to_string(0)," to known vertices\n",node_verts.to_string());
          throw std::logic_error("Too many matching to triangulated vertex");
        }
        if (known_idx.size() == 1u){
          local_map.push_back(poly_vert_idx[known_idx[0]]);
        } else {
          // check against all vertices. maybe this is the added Gamma point.
          auto punt_idx = norm(vertices_ - trij).find(brille::cmp::eq, 0.);
          if (punt_idx.size() < 1){
            info_update("For node ",i,":\n","No match of ",trij.to_string(),
                        " to known vertices\n",node_verts.to_string(),
                        " or extra_intersections\n",extra_intersections.to_string());
            throw std::runtime_error("No match to triangulated vertex");
          }
          local_map.push_back(punt_idx[0]);
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
  profile_update("  End of PolyhedronTrellis Construction");
  // the data_ PermutationTable should be initialised now:
  data_.initialize_permutation_table(vertices_.size(0), this->collect_keys());
}

template<class T, class S>
std::set<size_t>
PolyhedronTrellis<T,S>::collect_keys() {
  profile_update("Start of PolyhedronTrellis permutation key collection");
  std::set<size_t> keys;
  long long nnodes = brille::utils::u2s<long long, size_t>(nodes_.size());
  #pragma omp parallel for default(none) shared(keys, nnodes)
  for (long long sni=0; sni<nnodes; ++sni){
    size_t ni = brille::utils::s2u<size_t, long long>(sni);
    std::set<size_t> t = this->collect_keys_node(ni);
    if (t.size()){
      #pragma omp critical
      {
        keys.insert(t.begin(), t.end());
      }
    }
  }
  profile_update("  End of PolyhedronTrellis permutation key collection");
  return keys;
}

template<class T, class S>
std::set<size_t>
PolyhedronTrellis<T,S>::collect_keys_node(const ind_t node){
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
