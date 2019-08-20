void BrillouinZone::wedge_search(const bool pbv, const bool pok){
  debug_exec(std::string update_msg;)
  // Get the full pointgroup symmetry information
  PointSymmetry fullps = this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
  // And use it to find only the highest-order rotation operation along each
  // unique stationary axis.
  // PointSymmetry rotps = fullps.nfolds(1); // 1 to request only orders>1
  PointSymmetry rotps = fullps.higher(1); // 1 to request only orders>1
  // Get the vectors pointing to each full Brillouin zone facet cetre
  LQVec<double> xyz = cat(cat(this->get_points(), this->get_vertices()), this->get_half_edges());

  // status_update("xyz=\n", xyz.to_string());

  // if rotps is empty, there are no 2+-fold rotations, act like we have ̄1:
  if (rotps.size()==0){
    status_update("No 2+-fold rotation operations");
    this->ir_wedge_is_ok(xyz.first(1u)); // for ̄1, assigns this->ir_polyhedron
    return;
  }

  int max_order=rotps.order(rotps.size()-1); // since the rotations are sorted

  // rotps is sorted in increasing rotation-order order, but we want to sort
  // by increasing stationary-axis length (so that, e.g, [100] is dealt with
  // before [111]). We could ammend the PointSymmetry.sort method but since the
  // rotations are defined in a lattice that would only work for cubic systems;
  // e.g., in a hexagonal lattice [100], [010], and [1-10] are all the same
  // length. So we need to create a sorting permutation here for rotps to use.
  std::vector<size_t> perm(rotps.size());
  std::iota(perm.begin(), perm.end(), 0u); // 0u, 1u, 2u,...
  std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b){
    LQVec<int> lq(this->outerlattice, 2u);
    lq.set(0u, rotps.axis(a));
    lq.set(1u, rotps.axis(b));
    return lq.dot(0u,0u) < lq.dot(1u,1u);
  });
  // with the permutation found, permute rotps:
  rotps.permute(perm);

  status_update("Rotations:\n", rotps.getall());

  // make lattice vectors from the stationary axes
  LQVec<int> z(this->outerlattice, rotps.axes());
  // Find ̂zᵢ⋅ẑⱼ
  std::vector<std::vector<int>> dotij(z.size());
  double dottmp;
  for (size_t i=0; i<z.size(); ++i) for (size_t j=0; j<z.size(); ++j){
    dottmp = std::abs(z.dot(i,j)/z.norm(i)/z.norm(j));
    dotij[i].push_back( approx_scalar(dottmp,0.) ? -1: approx_scalar(dottmp, 1.) ? 1 : 0);
  }
  status_update("dot(i,j) flag\n", dotij);

  // Find a suitable in-rotation-plane vector for each stationary axis
  LQVec<double> x(this->outerlattice, rotps.perpendicular_axes());
  // or pick one of the BZ facet-points if the basis vector preffered flag
  if (pbv) for (size_t i=0; i<x.size(); ++i) for (size_t j=0; j<xyz.size(); ++j)
  if (dot(xyz.get(j), z.get(i)).all_approx(0.)){
    x.set(i, xyz.get(j));
    break;
  }

  /* Two rotations with ẑᵢ⋅ẑⱼ≠0 should have (Rⁿx̂ᵢ)⋅(Rᵐx̂ⱼ)=1 for some n,m.
     To start, assume that n and m are 0 and find the x̂ᵢ such that x̂ᵢ⋅(ẑᵢ×ẑⱼ)=1
  */
  std::vector<bool> handled(z.size(), false);
  bool flag;
  for (size_t i=0; i<z.size()-1; ++i) for (size_t j=i+1; j<z.size(); ++j)
  // if zᵢ and zⱼ are neither parallel or perpendicular
  if (0 == dotij[i][j]){
    if (!handled[i]){
      if (!pbv /*basis vectors not preferred*/) x.set(i, z.cross(i, j));
      if (handled[j] && !approx_scalar(x.dot(i,j), 0.) && x.dot(i,j)<0) x.set(i, -x.get(i));
      handled[i] = true;
    }
    if (!handled[j]){
      if (!pbv /*basis vectors not preferred*/){
        flag = norm(cross(x.get(i), z.get(j))).all_approx(">",1e-10);
        // if both or neither parallel is ok (pok) and zⱼ∥xᵢ, xⱼ=zᵢ×zⱼ; otherwise xⱼ=xᵢ
        // x.set(j, (pok^u_parallel_v) ? z.cross(i,j) : x.get(i));
        x.set(j, (pok||flag) ? x.get(i) : z.cross(i,j));
      }
      if (!approx_scalar(x.dot(i,j), 0.) && x.dot(i,j)<0) x.set(j, -x.get(j));
      handled[j] = true;
    }
  }

  LQVec<double> y = cross(z, x); // complete a right-handed coordinate system
  x= x/norm(x);
  y= y/norm(y);
  status_update("    z                x                        y           ");
  debug_exec(update_msg = "--------- -----------------------  -----------------------\n";)
  for (size_t i=0; i<z.size(); ++i){
    debug_exec(update_msg += z.to_string(i) +" "+  x.to_string(i) +" "+ y.to_string(i) +"\n";)
  }
  status_update(update_msg);

  /* Each symmetry operation in rotps is guaranteed to be a proper rotation,
     but there is no guarantee that it represents a *right handed* rotation.
     For the normal vectors to point the right way, we need to know if it is. */
  std::vector<bool> is_right_handed(z.size(), true);
  LQVec<double> Rv(this->outerlattice, 1u);
  for (size_t i=0; i<z.size(); ++i) if (rotps.order(i)>2){
    multiply_matrix_vector(Rv.datapointer(0), rotps.data(i), x.datapointer(i));
    if (dot(y.get(i), Rv.get(0)).getvalue(0) < 0)
      is_right_handed[i] = false;
  }
  status_update("Right handed:", is_right_handed);

  /* We now have for every symmetry operation a consistent vector in the
     rotation plane. We can now use z, x, and rotps to find and add the wedge
     normals. We should find one normal for each 2-fold axis and two normals
     for each 3-, 4-, and 6-fold axis. Plus we need space for one extra normal
     if there is only one rotation axis and the pointgroup has inversion.     */
  size_t found=0, expected = (rotps.size()==1) ? 1 : 0;
  for (size_t i=0; i<rotps.size(); ++i) expected += (rotps.order(i)>2) ? 2 : 1;
  LQVec<double> normals(this->outerlattice, expected);
  if (rotps.size() == 1){
    // We are assuming the system has inversion symmetry. We handle the case
    // where it doesn't elsewhere.
    this->wedge_normal_check(z.get(0), normals, found);
    this->ir_wedge_is_ok(normals.first(found)); // updates this->wedge_normals
  }
  bool accepted;
  int order;
  LQVec<double> vi(this->outerlattice, max_order), zi(this->outerlattice, 1u);
  LQVec<double> vxz, zxv;
  for (size_t i=0; i<rotps.size(); ++i){
    zi = is_right_handed[i] ? z.get(i) : -z.get(i);
    order = rotps.order(i);
    status_update("\nOrder ", order, ", z=", zi.to_string(0));
    accepted = false;
    vi.set(0, x.get(i)); // do something better here?
    for (int j=1; j<order; ++j) multiply_matrix_vector(vi.datapointer(j), rotps.data(i), vi.datapointer(j-1));
    zxv = cross(zi, vi.first(order));
    zxv = zxv/norm(zxv);
    debug_exec(update_msg ="          R^n v                 z x (R^n v)      \n";)
    debug_exec(update_msg+="------------------------ ------------------------\n";)
    for (int j=0; j<order; ++j){
      debug_exec(update_msg += vi.to_string(j) +" "+ zxv.to_string(j)+"\n");
    }
    status_update(update_msg);
    if (2==order){
      // one-normal version of wedge_normal_check allows for either ±n
      // → no need to check z×Rv = z×(-x) = -(z×x)
      accepted = this->wedge_normal_check(zxv.get(0), normals, found);
      /* if we couldn't add a 2-fold normal, we have more work to do. */
      if (!accepted){
        // for now, hope that 90° away is good enough
        this->wedge_normal_check(vi.get(0), normals, found);
      }
    } else {
      /* Consecutive acceptable normals *must* point into the irreducible wedge
         otherwise they destroy the polyhedron.
         Check between each pair of in-plane vectors (Rⁿx, Rⁿ⁺¹x),
         for n=[0,order) with the last check (Rᵒ⁻¹x,R⁰x)≡(Rᵒ⁻¹x,x)            */
      // if (this->isinside_wedge(zi, /*constructing=*/false).getvalue(0)){
        for (int j=0; j<order; ++j){
          if (accepted) break;
          accepted = this->wedge_normal_check(zxv.get(j), -zxv.get((j+1)%order), normals, found);
        }
      // } // isinside_wedge
    } // order>2
  }
  this->ir_wedge_is_ok(normals.first(found));
  status_update("wedge_search finished");
}
