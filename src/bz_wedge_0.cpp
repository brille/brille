#ifdef _MSC_VER
void BrillouinZone::wedge_search(const bool prefer_basis_vectors, const bool parallel_ok){
  debug_exec(std::string update_msg;)
  /*
  The Pointgroup symmetry information comes from, effectively, spglib which
  has all rotation matrices defined in the conventional unit cell -- which is
  our `outerlattice`. Consequently we must work in the outerlattice here.
  */
  PointSymmetry psym = this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
  // extract the symmetry operations with rotation axes: (not 1 or ̄1)
  std::vector<std::array<int,9>> rotations;
  for (size_t i=0; i<psym.size(); ++i)
    if (rotation_order(psym.data(i)) > 1)
      rotations.push_back(psym.get(i));
  // we can stop now if there are no 2+-fold rotations:
  if (rotations.size()<1){
    debug_update("No 2+-fold operations");
    return;
  }
  // sort the rotations by their orders
  std::sort(rotations.begin(), rotations.end(), [](std::array<int,9> a, std::array<int,9> b){
    return rotation_order(a.data()) < rotation_order(b.data()); // > for high to low, < for low to high
  });
  std::vector<int> orders;
  for (std::array<int,9> R: rotations) orders.push_back(rotation_order(R.data()));
  int max_order=0;
  for (int o: orders) if (o>max_order) max_order = o;

  debug_exec(update_msg = "Rotation orders:"; for (int o: orders) update_msg += " " + std::to_string(o);)
  debug_update(update_msg);

  //
  std::array<std::array<int,3>,3> axis_plane_vecs;
  LQVec<double> u(this->outerlattice,rotations.size());
  LQVec<double> x(this->outerlattice,rotations.size());
  LQVec<double> y(this->outerlattice,rotations.size());
  //
  for (size_t i=0; i<rotations.size(); ++i){
    axis_plane_vecs = rotation_axis_and_perpendicular_vectors(rotations[i].data());
    for (size_t j=0; j<3; ++j){
      u.insert(static_cast<double>(axis_plane_vecs[0][j]),i,j);
      x.insert(static_cast<double>(axis_plane_vecs[1][j]),i,j);
      // y.insert(static_cast<double>(axis_plane_vecs[2][j]),i,j);
    }
  }
  // check for zero-length rotation axes *before* trying to normalize them
  ArrayVector<bool> nonzero = norm(u).is_approx(">",0.);
  u = u.extract(nonzero);
  x = x.extract(nonzero);
  y = cross(u,x); // in case the basis is not orthogonal

  // find the unique rotation axes: (here we want to equate u and -u later)
  ArrayVector<size_t> u_equiv_idx = u.unique_idx();
  std::vector<size_t> unique_idx;
  unique_idx.push_back(u_equiv_idx.getvalue(0)); // the first vector is always unique
  bool is_in_unique_idx=false;
  for (size_t i=1; i<u.size(); ++i){
    is_in_unique_idx=false;
    for (size_t idx: unique_idx) if (u_equiv_idx.getvalue(i)==idx) is_in_unique_idx=true;
    if (!is_in_unique_idx) unique_idx.push_back(i);
  }
  debug_exec(\
  if (unique_idx.size() != u.is_unique().count_true())\
    throw std::runtime_error("Unique count is off.");\
  )

  // Verify that we have a right-handed rotation axis for each matrix
  LQVec<double> Rx(this->outerlattice, 1u);
  for (size_t i=0; i<rotations.size(); ++i){
    // Calculate Rx
    multiply_matrix_vector(Rx.data(0), rotations[i].data(), x.data(i));
    // y ≡ u × x, so if rotated x points away from y we need to flip u.
    if (dot(Rx, y.get(i)).getvalue(0) < 0) u.set(i, -u.get(i));
  }

  // Find ûᵢ⋅ûⱼ
  double* ui_dot_uj = new double[u.size()*u.size()]();
  for (size_t i=0; i<u.size(); ++i) for (size_t j=0; j<u.size(); ++j) ui_dot_uj[i*u.size()+j] = u.dot(i,j)/10;

  debug_exec(\
  update_msg = "unique(u):\n";\
  for (size_t i=0; i<u.size(); ++i)\
  if (u_equiv_idx.getvalue(i)==i) update_msg += u.to_string(i) + "\n";\
  update_msg += "dot(ui, uj)\n";\
  for (size_t i=0; i<u.size(); ++i)\
    if (u_equiv_idx.getvalue(i)==i){\
      for (size_t j=0; j<u.size(); ++j)\
      if (u_equiv_idx.getvalue(j)==j)\
      update_msg += approx_scalar(ui_dot_uj[i*u.size()+j],0.) ? " 0" : " x";\
      update_msg += "\n";\
    }\
  )
  debug_update(update_msg);

  // Two rotations with ûᵢ⋅ûⱼ ≠ 0 should have (Rⁿv̂ᵢ)⋅(Rᵐv̂ⱼ) = 1 for some n,m.
  // To start with, assume that n and m are 0 and find the v̂ᵢ such
  // that v̂ᵢ⋅(ûᵢ×ûⱼ)=1
  std::vector<bool> handled;
  for (size_t i=0; i<u.size(); ++i) handled.push_back(false);
  // we need a place to stash the cross product of uᵢ and uⱼ
  LQVec<double> primitive_basis(this->lattice, 3u);
  for (size_t i=0; i<3u; ++i) for (size_t j=0; j<3u; ++j) primitive_basis.insert(i==j?1:0,i,j);
  LQVec<double> xyz = transform_from_primitive(this->outerlattice, primitive_basis);
  debug_update("basis vectors\n", xyz.to_string());

  // why not just stash away the v̂ᵢ now?
  LQVec<double> v(x); // copy x
  bool u_parallel_v;
  for (size_t i=  0; i<u.size()-1; ++i) if(u_equiv_idx.getvalue(i) == i)
  for (size_t j=i+1; j<u.size()  ; ++j) if(u_equiv_idx.getvalue(j) == j)
  if (!approx_scalar(ui_dot_uj[i*u.size()+j], 0.0)){
    if (!handled[i]){
      if (prefer_basis_vectors){
      for (int k=3; k--;)
        if (!handled[i] && dot(u.extract(i), xyz.extract(k)).all_approx(0.)){
          v.set(i, xyz.extract(k));
          handled[i] = true;
        }
      }
      if (!handled[i]){
        v.set(i, u.cross(i,j));
        handled[i] = true;
      }
      if (handled[j] && !approx_scalar(v.dot(i, j), 0.) && v.dot(i, j)<0)
        v.set(i, -v.get(i));
    }
    if (!handled[j]) {
      if (prefer_basis_vectors){
        for (int k=3; k--;)
        if (!handled[j] && dot(u.extract(j), xyz.extract(k)).all_approx(0.)){
          v.set(j, xyz.extract(k));
          handled[j] = true;
        }
      }
      if (!handled[j]){
        /* If j has been handled already, we don't want to just set  vᵢ to vⱼ
        in case of, e.g., [100],[010],[111] where [100]×[111] is [0̄11]
        and [010]×[111] is [10̄1]                                           */
        // protect against v ∥ u
        u_parallel_v = approx_scalar(norm(cross(v.get(i), u.get(j))).getvalue(0)/10, 0.0);
        if (!(parallel_ok ^ u_parallel_v)){
          v.set(j, v.get(i));
        } else {
          v.set(j, u.cross(i,j));
        }
        handled[j] = true;
      }
      if (!approx_scalar(v.dot(i, j), 0.) && v.dot(i, j) < 0)
        v.set(j, -v.get(j));
    }
  }
  delete[] ui_dot_uj;
  /* for the case where three rotations have ûᵢ⋅ûⱼ=0, ûᵢ⋅ûₖ≠0, ûⱼ⋅ûₖ≠0
    the above will set vᵢ=v̂ₖ=ûᵢ×ûₖ and v̂ⱼ=ûⱼ×ûₖ and we assume (hope?) that
    there exists some Rₙv̂ₖ=ûⱼ×ûₖ.     */
  // Now go back through and copy the unique (a,b) values to the non-uniques.
  for (size_t i=0; i<u.size(); ++i)
  if (u_equiv_idx.getvalue(i)!=i) v.set(i,v.get(u_equiv_idx.getvalue(i)));

  debug_update("u_equiv_idx = " + u_equiv_idx.to_string(" "));
  debug_update("u =\n" +  u.to_string());
  debug_update("v =\n" + v.to_string());

  // We now have for every rotation a consistent vector in the rotation plane.
  // It's time to use ̂u, ̂v, R, and order(R) to find/add the wedge normals
  LQVec<double> normals(this->outerlattice, 2*rotations.size());
  LQVec<double> vj(this->outerlattice, max_order);
  size_t total_found=0;
  /* If there is only one unique rotation axis then we are guaranteed to get
     the *wrong* irreducible wedge by this method if the pointgroup has ̄1 as
     a symmetry elelement. Since we're implicitly assuming this we will always
     end up with double the real irreducible zone whether ̄1 is present or not.*/
  size_t unq_count=0;
  for (size_t i=0; i<u.size(); ++i) if (u_equiv_idx.getvalue(i)==i) ++unq_count;
  if (unq_count == 1){
    // find the unique one, and insert it as a wedge normal:
    for (size_t i=0; i<u.size(); ++i) if (u_equiv_idx.getvalue(i)==i)
    this->wedge_normal_check(u.get(i), normals, total_found);
  }
  bool accepted=false;
  int order;
  for (size_t j=0; j<rotations.size(); ++j){
    order = orders[j];
    debug_update("r, u: " + std::to_string(order) + " (" + u.to_string(j,")"));
    accepted=false;
    if (2==order){
      // the single normal version of add_wedge_normal_check allows for the
      // possibility of either n or -n being added, so we don't need to check
      // both u×v and u×Rv (as Rv = -v for a 2-fold axis).
      accepted = this->wedge_normal_check(cross(u.get(j), v.get(j)), normals, total_found); // increments total_found if adding was a success
      //// but a 2-fold axis could always(?) be folded 90 degrees away too(?)
      if (!accepted) this->wedge_normal_check(v.get(j), normals, total_found);
    } else {
      vj.set(0, v.get(j));
      for (int k=1; k<order; ++k)
        multiply_matrix_vector(vj.data(k), rotations[j].data(), vj.data(k-1));
      debug_update("R^n v\n", vj.to_string());
      // consecutive acceptable normals *must* point into the irreducible wedge
      // and we need to check between Rⁿ⁻¹v and Rⁿv=Iv=v, so k and (k+1)%n
      for (int k=0; k<order; ++k){
        if (accepted) break;
        accepted = this->wedge_normal_check(cross(u.get(j), vj.get(k)),
                                            cross(vj.get((k+1)%order), u.get(j)),
                                            normals, total_found);
      }
    }
  }
  this->ir_wedge_is_ok(normals.first(total_found)); // assigns the ir_polyhedron
  // otherwise ir_wedge_is_ok is called by wedge_normal_check
  debug_update("wedge_search finished");
}

#else // if not _WIN32

void BrillouinZone::wedge_search(const bool prefer_basis_vectors, const bool parallel_ok){
  debug_exec(std::string update_msg;)
  /*
  The Pointgroup symmetry information comes from, effectively, spglib which
  has all rotation matrices defined in the conventional unit cell -- which is
  our `outerlattice`. Consequently we must work in the outerlattice here.
  */
  PointSymmetry psym = this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
  // extract the symmetry operations with rotation axes: (not 1 or ̄1)
  std::vector<std::array<int,9>> rotations;
  for (size_t i=0; i<psym.size(); ++i)
    if (rotation_order(psym.data(i)) > 1)
      rotations.push_back(psym.get(i));
  //
  LQVec<double> primitive_basis(this->lattice, 3u);
  for (size_t i=0; i<3u; ++i) for (size_t j=0; j<3u; ++j) primitive_basis.insert(i==j?1:0,i,j);
  LQVec<double> xyz = transform_from_primitive(this->outerlattice, primitive_basis);
  debug_update("basis vectors\n", xyz.to_string());
  // we can stop now if there are no 2+-fold rotations:
  if (rotations.size()<1){
    debug_update("No 2+-fold operations");
    this->ir_wedge_is_ok(xyz.first(1u)); // assigns the ir_polyhedron
    return;
  }
  // temporary storage for the stationary vector and two perpendicular vectors
  // for each symmetry operations
  std::array<std::array<int,3>,3> axis_plane_vecs;
  // sort the rotations by their orders
  std::sort(rotations.begin(), rotations.end(), [&](std::array<int,9> a, std::array<int,9> b){
    int ra,rb;
    ra = rotation_order(a.data());
    rb = rotation_order(b.data());
    if (ra == rb){
      LQVec<int> lu(this->outerlattice, 2u);
      axis_plane_vecs = rotation_axis_and_perpendicular_vectors(a.data());
      for (size_t i=0; i<3; ++i) lu.insert(axis_plane_vecs[0][i], 0u, i);
      axis_plane_vecs = rotation_axis_and_perpendicular_vectors(b.data());
      for (size_t i=0; i<3; ++i) lu.insert(axis_plane_vecs[0][i], 1u, i);
      return lu.norm(0) < lu.norm(1); // always sort shorter rotation axes first
    }
    return ra < rb; // > for high to low, < for low to high
  });
  std::vector<int> orders;
  for (std::array<int,9> R: rotations) orders.push_back(rotation_order(R.data()));

  // storage for the stationary vector (u) and two perpendicular vectors (x & y)
  LQVec<double> u(this->outerlattice,rotations.size());
  LQVec<double> x(this->outerlattice,rotations.size());
  LQVec<double> y(this->outerlattice,rotations.size());
  // get the temporary array, then fill u and x
  for (size_t i=0; i<rotations.size(); ++i){
    axis_plane_vecs = rotation_axis_and_perpendicular_vectors(rotations[i].data());
    for (size_t j=0; j<3; ++j){
      u.insert(static_cast<double>(axis_plane_vecs[0][j]),i,j);
      x.insert(static_cast<double>(axis_plane_vecs[1][j]),i,j);
    }
  }
  // double check that we haven't let through any vector-less operations
  ArrayVector<bool> nonzero = norm(u).is_approx(">",0.);
  u = u.extract(nonzero);
  x = x.extract(nonzero);
  // and complete a right-handed orthogonal coordinate system for each
  y = cross(u,x);

  // find the unique rotation axes: (here we want to equate u and -u later)
  ArrayVector<size_t> u_equiv_idx = u.unique_idx();
  std::vector<size_t> unique_idx;
  unique_idx.push_back(u_equiv_idx.getvalue(0)); // the first vector is always unique
  bool is_in_unique_idx=false;
  for (size_t i=1; i<u.size(); ++i){
    is_in_unique_idx=false;
    for (size_t idx: unique_idx) if (u_equiv_idx.getvalue(i)==idx) is_in_unique_idx=true;
    if (!is_in_unique_idx) unique_idx.push_back(i);
  }
  // make sure we have the right-highest-order operation in unique_idx:
  size_t idx;
  LQVec<double> Rx(this->outerlattice, 1u);
  for (size_t i=0; i<unique_idx.size(); ++i){
    idx = unique_idx[i];
    for (size_t j=0; j<u.size(); ++j) if (idx == u_equiv_idx.getvalue(j)){
      // check if R*x points in the same direction as u×x≡y
      multiply_matrix_vector(Rx.data(0), rotations[j].data(), x.data(j));
      if (dot(Rx, y.get(j)).getvalue(0)>0) unique_idx[i] = j;
    }
  }

  // debug_exec(update_msg = "All rotation orders:"; for (int o:orders) update_msg += " "+std::to_string(o);)
  // debug_update(update_msg);
  // debug_exec(update_msg = "All equivalent idx:";)
  // for (size_t i=0;i<u.size(); ++i) debug_exec(update_msg += " "+std::to_string(u_equiv_idx.getvalue(i));)
  // debug_update(update_msg);
  // debug_update("all u\n", u.to_string());

  LQVec<double> newu(this->outerlattice, unique_idx.size());
  LQVec<double> newx(this->outerlattice, unique_idx.size());
  LQVec<double> newy(this->outerlattice, unique_idx.size());
  std::vector<std::array<int, 9>> unique_rotations(unique_idx.size());
  std::vector<int> unique_orders(unique_idx.size());
  for (size_t i=0; i<unique_idx.size(); ++i){
    newu.set(i, u.get(unique_idx[i]));
    newx.set(i, x.get(unique_idx[i]));
    newy.set(i, y.get(unique_idx[i]));
    unique_rotations[i]=rotations[unique_idx[i]];
    unique_orders[i]=orders[unique_idx[i]];
  }
  u = newu;
  x = newx;
  y = newy;
  rotations = unique_rotations;
  orders = unique_orders;

  int max_order=0;
  for (int o: orders) if (o>max_order) max_order = o;

  debug_exec(update_msg = "Rotation orders:"; for (int o: orders) update_msg += " " + std::to_string(o);)
  debug_update(update_msg);

  // Find unique ûᵢ⋅ûⱼ
  size_t count = u.size();
  double* ui_dot_uj = new double[count*count]();
  for (size_t i=0; i<count; ++i) for (size_t j=0; j<count; ++j)
  ui_dot_uj[count*i+j] = u.dot(i,j)/10;

  debug_update("dot(ui, uj)");
  debug_exec(update_msg = "");
  for(size_t i=0; i<count; ++i){
    for (size_t j=0; j<count; ++j)
      debug_exec(update_msg+=approx_scalar(ui_dot_uj[i*count+j],0.)?" 0":" x");
    debug_exec(update_msg+="\n");
  }
  debug_update(update_msg);

  // Two rotations with ûᵢ⋅ûⱼ ≠ 0 should have (Rⁿv̂ᵢ)⋅(Rᵐv̂ⱼ) = 1 for some n,m.
  // To start with, assume that n and m are 0 and find the v̂ᵢ such
  // that v̂ᵢ⋅(ûᵢ×ûⱼ)=1
  std::vector<bool> handled;
  for (size_t i=0; i<count; ++i) handled.push_back(false);

  LQVec<double> v(x); // copy x in case any of the rotations are fully independent
  bool u_parallel_v;
  for (size_t i=0; i<count-1; ++i)
  for (size_t j=i+1; j<count; ++j)
  if (!approx_scalar(ui_dot_uj[i*count+j], 0.0)){
    if (!handled[i]){
      if (prefer_basis_vectors){
      for (int k=3; k--;)
        if (!handled[i] && dot(u.extract(i), xyz.extract(k)).all_approx(0.)){
          v.set(i, xyz.extract(k));
          handled[i] = true;
        }
      }
      if (!handled[i]){
        v.set(i, u.cross(i, j));
        handled[i] = true;
      }
      if (handled[j] && !approx_scalar(v.dot(i, j), 0.) && v.dot(i, j)<0)
        v.set(i, -v.get(i));
    }
    if (!handled[j]) {
      if (prefer_basis_vectors){
        for (int k=3; k--;)
        if (!handled[j] && dot(u.extract(j), xyz.extract(k)).all_approx(0.)){
          v.set(j, xyz.extract(k));
          handled[j] = true;
        }
      }
      if (!handled[j]){
        /* If j has been handled already, we don't want to just set  vᵢ to vⱼ
        in case of, e.g., [100],[010],[111] where [100]×[111] is [0̄11]
        and [010]×[111] is [10̄1]                                           */
        // protect against v ∥ u
        u_parallel_v = approx_scalar(norm(cross(v.get(i), u.get(j))).getvalue(0)/10, 0.0);
        if (!(parallel_ok ^ u_parallel_v)){
          v.set(j, v.get(i));
        } else {
          v.set(j, u.cross(i, j));
        }
        handled[j] = true;
      }
      if (!approx_scalar(v.dot(i, j), 0.) && v.dot(i, j) < 0)
        v.set(j, -v.get(j));
    }
  }
  delete[] ui_dot_uj;

  debug_update("u =\n" +  u.to_string());
  debug_update("v =\n" + v.to_string());

  // We now have for every rotation a consistent vector in the rotation plane.
  // It's time to use ̂u, ̂v, R, and order(R) to find/add the wedge normals
  /*Allocate space for an extra normal in case there is one unique axis and
    the system has inversion symmetry and 3-fold-or-higher rotational symmetry
    in which case we need 2+1 normals to describe the irreducible wedge */
  LQVec<double> normals(this->outerlattice, 2*count+1);
  LQVec<double> vj(this->outerlattice, max_order);
  LQVec<double> uj(this->outerlattice, 1u);
  size_t total_found=0;
  /* If there is only one unique rotation axis then we are guaranteed to get
     the *wrong* irreducible wedge by this method if the pointgroup has ̄1 as
     a symmetry elelement. Since we're implicitly assuming this we will always
     end up with double the real irreducible zone whether ̄1 is present or not.*/
  if (count == 1){
    // insert the unique axis it as a wedge normal:
    this->wedge_normal_check(u.get(0), normals, total_found);
    this->ir_wedge_is_ok(normals.first(total_found)); // to ensure the wedge normals are up to date
    debug_update("_____ ", total_found, " ____ normals found");
  }
  bool accepted=false;
  int order;
  for (size_t j=0; j<count; ++j){
    uj = u.get(j);
    order = orders[j];
    debug_update("r, u: " + std::to_string(order) + " (" + uj.to_string(0,")"));
    accepted=false;
    vj.set(0, v.get(j));
    for (int k=1; k<order; ++k)
    multiply_matrix_vector(vj.data(k), rotations[j].data(), vj.data(k-1));
    debug_update("R^n v\n", vj.first(order).to_string());
    if (2==order){
      // the single normal version of add_wedge_normal_check allows for the
      // possibility of either n or -n being added, so we don't need to check
      // both u×v and u×Rv (as Rv = -v for a 2-fold axis).
      accepted = this->wedge_normal_check(cross(uj, vj.get(0)), normals, total_found); // increments total_found if adding was a success
      //// but a 2-fold axis could always(?) be folded 90 degrees away too(?)
      if (!accepted) this->wedge_normal_check(vj.get(0), normals, total_found);
    } else {
      // consecutive acceptable normals *must* point into the irreducible wedge
      // and we need to check between Rⁿ⁻¹v and Rⁿv=Iv=v, so k and (k+1)%n
      if (this->isinside_wedge(uj,/*constructing=*/true).getvalue(0))
      for (int k=0; k<order; ++k){
        if (accepted) break;
        accepted = this->wedge_normal_check(cross(uj, vj.get(k)),
                                            cross(vj.get((k+1)%order), uj),
                                            normals, total_found);
      }
      if (this->isinside_wedge(-uj,/*constructing=*/true).getvalue(0))
      for (int k=0; k<order; ++k){
        if (accepted) break;
        accepted = this->wedge_normal_check(cross(-uj, vj.get(k)),
                                            cross(vj.get((k+1)%order), -uj),
                                            normals, total_found);
      }
    }
    debug_update("_____ ", total_found, " ____ normals found");
  }
  this->ir_wedge_is_ok(normals.first(total_found)); // assigns the ir_polyhedron
  // otherwise ir_wedge_is_ok is called by wedge_normal_check
  debug_update("wedge_search finished");
}

#endif
