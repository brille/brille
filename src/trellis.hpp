template<typename T>
PolyhedronTrellis<T>::PolyhedronTrellis(const Polyhedron& poly, const double fraction):
  polyhedron_(poly), vertices_({3,0})
{
  double node_volume = ((fraction> 0. && fraction<1.) ? fraction : 1.)*poly.get_volume();
  double node_length = std::cbrt(node_volume);
  // find the extents of the polyhedron
  std::array<std::array<double,2>,3> minmax;
  for (size_t i=0; i<3u; ++i){
    minmax[i][0] = (std::numeric_limits<double>::max)();
    minmax[i][1] = std::numeric_limits<double>::lowest();
  }
  const ArrayVector<double>& pv{poly.get_vertices()};
  for (size_t i=0; i<pv.size(); ++i){
    for (size_t j=0; i<3u; ++j){
      if (pv.getvalue(i,j) < minmax[j][0]) minmax[j][0] = pv.getvalue(i,j);
      if (pv.getvalue(i,j) > minmax[j][1]) minmax[j][1] = pv.getvalue(i,j);
    }
  }
  // construct the trellis boundaries:
  // boundaries_ = std::array<std::vector<double>,3>({std::vector<double>(),std::vector<double>(),std::vector<double>()});
  for (size_t i=0; i<3u; ++i){
    boundaries_[i].push_back(minmax[i][0]);
    while (boundaries_[i].back() < minmax[i][1])
      boundaries_[i].push_back(boundaries_[i].back()+node_length);
    debug_update("PolyhedronTrellis has ",boundaries_[i].size()-1," bins along axis ",i,", with boundaries ",boundaries_[i]);
  }
  // find which trellis intersections are inside of the polyhedron
  std::vector<std::array<double,3>> va_int;
  for (double x: boundaries_[0]) for (double y: boundaries_[1]) for (double z: boundaries_[2])
    va_int.push_back({x,y,z});
  ArrayVector<double> all_intersections(va_int);
  ArrayVector<bool> are_inside = poly.contains(all_intersections);
  // these will be retained as node vertices
  std::vector<size_t> map_idx, keep_idx = find(are_inside);
  for (size_t i=0; i<all_intersections.size(); ++i)
    map_idx.push_back(std::distance(keep_idx.begin(), std::find(keep_idx.begin(), keep_idx.end(), i)));

  // find the nodes which are fully inside the polyhedron
  /* Each node with linear index idx has a subscripted index (i,j,k)
     and is surrounded by the trellis intersections of boundaries
     (i,j,k) + { (000), (100), (110), (010), (101), (001), (011), (111)};
  */
  size_t nNodes = this->node_count();
  // the order of the cube node intersections is paramount:
  std::vector<std::array<size_t,3>> node_intersections{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{1,0,1},{0,0,1},{0,1,1},{1,1,1}};
  std::array<size_t,3> intersections_span{1,boundaries_[0].size(),boundaries_[0].size()*boundaries_[1].size()};
  std::vector<bool> node_is_cube(nNodes, true), node_is_outside(nNodes, false);
  for (size_t i=0; i<nNodes; ++i){
    std::array<size_t,3> node_ijk = this->idx2sub(i);
    // this node is a cube if all node intersection vertices are mapped
    for (auto ni: node_intersections) if (node_is_cube[i]) {
      size_t intersection_idx = 0;
      for (size_t j=0; j<3; ++j) intersection_idx += (node_ijk[j]+ni[j])*intersections_span[j];
      node_is_cube[i] = map_idx[intersection_idx] < keep_idx.size();
    }
    if (!node_is_cube[i]){
      bool any_mapped = false;
      for (auto ni: node_intersections) {
        size_t intersection_idx = 0;
        for (size_t j=0; j<3; ++j) intersection_idx += (node_ijk[j]+ni[j])*intersections_span[j];
        any_mapped |= map_idx[intersection_idx] < keep_idx.size();
      }
      node_is_outside[i] = !any_mapped;
    }
  }
  // Pull out the intersection points which we will keep as vertices:
  ArrayVector<double> kept_intersections = all_intersections.extract(are_inside);
  ArrayVector<double> extra_intersections(3u, 3*(all_intersections.size()-kept_intersections.size()));
  size_t nExtra=0;
  std::vector<size_t> extra_map;
  // Now actually create the node objects (which are not fully outside)
  for (size_t i=0; i<nNodes; ++i) if (!node_is_outside[i]) {
    std::array<size_t,3> node_ijk = this->idx2sub(i);
    std::array<size_t,8> vert_idx; // the 8 vertex indices of the cube
    for (size_t k=0; k<8u; ++k){
      size_t int_idx=0;
      for (size_t j=0; j<3; ++j) int_idx += (node_ijk[j]+node_intersections[k][j])*intersections_span[j];
      vert_idx[k] = map_idx[int_idx];
    }
    if (node_is_cube[i]) {
      nodes_[i] = CubeNode(vert_idx);
    } else {
      // This node intersects the polyhedron. First, find the interior part
      std::array<double,3> min_corner, max_corner;
      for (size_t j=0; j<3u; ++j){
        min_corner[j] = boundaries_[j][node_ijk[j]  ];
        max_corner[j] = boundaries_[j][node_ijk[j]+1];
      }
      Polyhedron cube = polyhedron_box(min_corner, max_corner);
      Polyhedron cut_cube = Polyhedron::bisect(cube, poly.get_normals(), poly.get_points());
      // Then triangulate it into tetrahedra
      SimpleTet tri_cut(cut_cube);
      // we need to retain the *additional* vertices from the triangulation
      // and figure-out a local vertex index mapping
      std::vector<size_t> local_map;
      const ArrayVector<double>& triverts{tri_cut.get_vertices()};
      for (size_t j=0; j<triverts.size(); ++j){
        std::vector<size_t> cube_idx = find(norm(cube.get_vertices()-triverts.extract(j)).is_approx("=",0.));
        if (cube_idx.size()>1) throw std::logic_error("Too many matching vertices");
        if (cube_idx.size()==1){
          // we could try to use the cube polyhedron vertices to map-back to the full trellis intersections
          // but that seems hard, so as a first attempt check this vertex against *all* kept intersections
          cube_idx = find(norm(kept_intersections-triverts.extract(j)).is_approx("=",0.));
          if (cube_idx.size()==1) local_map.push_back(cube_idx[0]);
          else throw std::logic_error("Problem finding index for triangulated vertex");
        } else {
          cube_idx = find(norm(extra_intersections-triverts.extract(j)).is_approx("=",0.));
          if (cube_idx.size()>0) local_map.push_back(extra_map[cube_idx[0]]);
          else {
            // make sure we have room to store this new intersection
            // if we don't, make room; but since this involves a memory copy make lots of room
            if (extra_intersections.size() < nExtra) extra_intersections.resize(2*nExtra);
            // store the extra intersection
            extra_intersections.set(nExtra, triverts.extract(j));
            // and its mapping information
            extra_map.push_back(kept_intersections.size()+nExtra);
            local_map.push_back(kept_intersections.size()+nExtra);
            ++nExtra;
          }
        }
      }
      std::vector<std::array<size_t,4>> idx_per_tet;
      const ArrayVector<size_t>& local_ipt{tri_cut.get_vertices_per_tetrahedron()};
      for (size_t j=0; j<local_ipt.size(); ++j){
        std::array<size_t,4> one_tet{0,0,0,0};
        for (size_t k=0; k<4u; ++k) one_tet[k] = local_map[local_ipt.getvalue(j,k)];
        idx_per_tet.push_back(one_tet);
      }
      nodes_[i] = PolyNode(idx_per_tet);
    }
  }
  // Now all non-empty nodes have been populated with the indices of their vertices
  // Combine the retained trellis vertices and the extra triangulated vertices
  vertices_ = cat(kept_intersections, extra_intersections);
}


template<class T>
ArrayVector<double> TrellisData<T>::debye_waller_sum(const LQVec<double>& Q, const double t_K) const{
  return this->debye_waller_sum(Q.get_xyz(), t_K);
}

template<class T>
ArrayVector<double> TrellisData<T>::debye_waller_sum(const ArrayVector<double>& Q, const double t_K) const{
  const double hbar = 6.582119569E-13; // meV⋅s
  const double kB   = 8.617333252E-2; // meV⋅K⁻¹
  size_t nQ = Q.size();
  size_t nIons = elements_[1] / 3u; // already checked to be correct
  ArrayVector<double> WdQ(nIons,nQ); // Wᵈ(Q) has nIons entries per Q point
  double coth_en, Q_dot_e_2;
  size_t span = 1u + nIons*3u + elements_[2] + elements_[3]*elements_[3];
  size_t nq = shape_[0];

  const double beta = kB*t_K; // meV
  const double pref{hbar*hbar/static_cast<double>(2*nq)}; // meV²⋅s²

  double qj_sum;
  // for each input Q point
  for (size_t Qidx=0; Qidx<nQ; ++Qidx){
    // and each ion
    for (size_t d=0; d<nIons; ++d){
      qj_sum = double(0);
      // sum over all reduced q in the first Brillouin zone
      for (size_t q=0; q<nq; ++q){
        // and over all 3*nIon branches at each q
        for (size_t j=0; j<branches_; ++j){
          // for each branch energy, find <2nₛ+1>/ħωₛ ≡ coth(2ħωₛβ)/ħωₛ
          coth_en = coth_over_en(data_.getvalue(q,j*span), beta);
          // and find |Q⋅ϵₛ|². Note: vector_product(x,y) *is* |x⋅y|²
          Q_dot_e_2 = vector_product(3u, Q.data(Qidx), data_.data(q,j*span+1u+3u*d));
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
template<template<class> class A>
ArrayVector<double> TrellisData<T>::debye_waller(const A<double>& Q, const std::vector<double>& M, const double t_K) const{
  size_t nIons = elements_[1] / 3u;
  if (0 == nIons || elements_[1]*3u != nIons)
    throw std::runtime_error("Debye-Waller factor requires 3-vector eigenvector(s).");
  if (M.size() != nIons)
    throw std::runtime_error("Debye-Waller factor requires an equal number of ions and masses.");
  ArrayVector<double> WdQ = this->debye_waller_sum(Q, t_K);
  ArrayVector<double> factor(1u, Q.size());
  double d_sum;
  for (size_t Qidx=0; Qidx<Q.size(); ++Qidx){
    d_sum = double(0);
    for (size_t d=0; d<nIons; ++d){
      d_sum += std::exp(WdQ.getvalue(Qidx, d)/M[d]);
    }
    factor.insert(d_sum*d_sum, Qidx);
  }
  return factor;
}
