// template<class T>
// static bool next_set(const std::vector<std::vector<T>>& equiv, std::vector<size_t>& idx, std::vector<T>& set){
//   bool goon = true;
//   size_t n = idx.size();
//   while (goon && n-- > 0u){
//     if (++idx[n] < equiv[n].size())
//       goon = false;
//     else
//       idx[n] = 0u;
//   }
//   for (size_t i=0; i<idx.size(); ++i) set[i] = equiv[i][idx[i]];
//   return goon; // overflow if goon is still true
// }

template <class T> ArrayVector<bool> keep_if(LQVec<T>& normals, LQVec<T>& points){
  // determine whether we should keep points based on the provided normals.
  // We want to only keep those points in the positive half-space for any one
  // normal, but if all normals contribute to reduce the remaining points to
  // lie in a single plane we instead want to keep all points.
  ArrayVector<bool> keep(1u, points.size(), true);
  std::vector<size_t> nop(normals.size(), 0); // number of not-on-plane points
  for (size_t i=0; i<normals.size(); ++i)
    nop[i] = dot(normals.extract(i), points).is_approx(">",0.).count_true();
  // If there are no planes with 0 off-plane points, divide the space
  if (std::find(nop.begin(),nop.end(),0u)==nop.end())
    for (size_t i=0; i<points.size(); ++i)
      keep.insert(dot(normals, points.extract(i)).all_approx(">=",0.), i);
  return keep;
}

void BrillouinZone::wedge_brute_force(void){
  debug_exec(std::string msg;)
  Spacegroup sg = this->outerlattice.get_spacegroup_object();
  // For P1 the irreducible Brillouin zone *is* the first Brillouin zone
  // which has already been set
  if (/*P 1*/sg.get_international_table_number() < 2) return;
  // The other triclinic spacegroup, P -1, is too complicated for this algorithm
  if (/*P 1 or P-1*/sg.get_international_table_number() < 3){
    this->wedge_triclinic();
    return;
  }
  // get the characteristic points of the first Brillouin zone
  // the face centers and corners:
  LQVec<double> special = cat(this->get_points(), this->get_vertices());
  // add the mid-facet-edge points:
  special = cat(special, this->get_half_edges());

  // // Since we're assuming inversion symmetry (rectified later if necessary),
  // // for every Friedel-pair of points we can immediately discard one
  // // (half of all points for a crystallographic spacegroup should be discarded)
  // std::vector<bool> isok(special.size(), true);
  // for (size_t i=0; i<special.size()-1; ++i) if (isok[i])
  // for (size_t j=i+1; j<special.size(); ++j) if (isok[j])
  // if (special.vector_approx(i,j,"*",-1.0)) isok[j]=false;
  // // isok[i] == true → point i is the first half of a Friedel-pair or a loner
  // special = special.extract(isok);

  // Grab the pointgroup symmetry operations
  PointSymmetry fullps = this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
  // Now restrict the symmetry operations to those with order > 1.
  PointSymmetry ps = fullps.higher(1);
  // PointSymmetry ps = fullps.nfolds(1);
  // ps is sorted in increasing rotation-order order, but we want to sort
  // by increasing stationary-axis length (so that, e.g, [100] is dealt with
  // before [111]).
  std::vector<size_t> perm(ps.size());
  std::iota(perm.begin(), perm.end(), 0u); // 0u, 1u, 2u,...
  std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b){
    LQVec<int> lq(this->outerlattice, 2u);
    lq.set(0u, ps.axis(a));
    lq.set(1u, ps.axis(b));
    return lq.dot(0u,0u) < lq.dot(1u,1u);
  });
  // std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b){
  //   return ps.order(a) > ps.order(b);
  // });
  ps.permute(perm); // ps now sorted with shortest stationary axis first

  ArrayVector<bool> keep(1u, 0u);
  // if (fullps.has_space_inversion()){
  //   // ps does not include 1 or -1, so all operations have stationary axes
  //   // the stationary axis of each rotation is, of course, express in real space!
  //   LDVec<int> all_axes(this->outerlattice.star(), ps.axes());
  //   if (all_axes.is_unique().count_true() == 1u){
  //     status_update("Deal with -1 since there is only one stationary axis");
  //     // use the stationary axis as the normal of a mirror plane
  //     keep = dot(all_axes.extract(0).star(), special).is_approx(">=", 0.);
  //     special = special.extract(keep);
  //   }
  // }
  status_update("Deal with 2-fold axes along highest-symmetry directions first");

  std::vector<bool> sym_unused(ps.size(), true);
  // deal with any two-fold axes along êᵢ first:
  // The stationary vector of each rotation is a real space vector!
  LDVec<int> vec(this->outerlattice.star(), 1u);// must be int since ps.axis returns array<int,3>
  LQVec<double> nrm(this->outerlattice, 1u);
  ArrayVector<bool> is_ei(1u, 9u);
  std::vector<std::array<double,3>> eiv{{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,-1,0},{1,0,1},{0,1,1},{1,0,-1},{0,1,-1}};
  LQVec<double> eis(this->outerlattice, 9u);
  for (size_t i=0; i<eis.size(); ++i) eis.set(i, eiv[i]);
  LDVec<double> reis(this->outerlattice.star(), 2u);
  for (size_t i=0; i<reis.size(); ++i) reis.set(i, eiv[i]);

  for (size_t i=0; i<ps.size(); ++i) if (ps.order(i)==2){
    vec.set(0, ps.axis(i));
    // Stationary axis along reciprocal space basis vector
    is_ei = norm(cross(eis, vec.star())).is_approx("==", 0.);
    if (is_ei.count_true() == 1u){
      status_update("Is 2-fold axis ",i," ei*?: ",is_ei.to_string(" "));
      if (is_ei.getvalue(0)) nrm.set(0, eis.cross(2,0)); /* (100)* → n=(001)*×(100)* */
      if (is_ei.getvalue(1)) nrm.set(0, eis.cross(0,1)); /* (010)* → n=(100)*×(010)* */
      if (is_ei.getvalue(2)) nrm.set(0, eis.cross(1,2)); /* (001)* → n=(010)*×(001)* */
      if (is_ei.getvalue(3)) nrm.set(0, eis.cross(3,2)); /* (110)* → n=(110)*×(001)* */
      if (is_ei.getvalue(4)) nrm.set(0, eis.cross(2,4)); /* (1̄10)* → n=(001)*×(1̄10)* */
      if (is_ei.getvalue(5)) nrm.set(0, eis.cross(1,5)); /* (101)* → n=(010)*×(101)* */
      if (is_ei.getvalue(6)) nrm.set(0, eis.cross(6,0)); /* (011)* → n=(011)*×(100)* */
      if (is_ei.getvalue(7)) nrm.set(0, eis.cross(7,1)); /* (1̄0̄1)* → n=(10̄1)*×(010)* */
      if (is_ei.getvalue(8)) nrm.set(0, eis.cross(0,8)); /* (01̄1)* → n=(100)*×(01̄1)* */
      // keep any special points beyond the bounding plane
      keep = dot(nrm, special).is_approx(">=", 0.);
      status_update("Keeping special points with\n",nrm.to_string(0)," dot p >= 0:\n", special.to_string(keep));
      special = special.extract(keep);
      sym_unused[i] = false;
    }
    // Stationary axis along real space basis vector
    is_ei = norm(cross(reis, vec)).is_approx("==", 0.);
    if (sym_unused[i] && is_ei.count_true() == 1u){
      status_update("Is 2-fold axis ",i," ei?: ",is_ei.to_string(" "));
      if (is_ei.getvalue(0)) nrm.set(0, eiv[1]); /* (100) → n = (010)* */
      if (is_ei.getvalue(1)) nrm.set(0, eiv[0]); /* (010) → n = (100)* */
      // keep any special points beyond the bounding plane
      keep = dot(nrm, special).is_approx(">=", 0.);
      status_update("Keeping special points with\n",nrm.to_string(0)," dot p >= 0:\n", special.to_string(keep));
      special = special.extract(keep);
      sym_unused[i] = false;
    }
  }

  status_update("Now figure out how all special points are related for each symmetry operation");
  // Find which points are mapped onto equivalent points by each symmetry operation
  std::vector<std::vector<size_t>> one_sym;
  std::vector<size_t> one_type, type_order;
  std::vector<bool> unfound(special.size(), true), type_unfound;
  for (size_t i=0; i<ps.size(); ++i) if (sym_unused[i]){
    status_update("Unused symmetry ",i," with order ", ps.order(i));
    one_sym.clear();
    for (auto b: unfound) b = true;
    for (size_t j=0; j<special.size(); ++j) if (unfound[j]){
      one_type.clear();
      one_type.push_back(j);
      unfound[j] = false;
      for (size_t k=j+1; k<special.size(); ++k)
      // if (unfound[k] && special.rotate_approx(k, j, ps.get(i), -ps.order(i))){ // -order checks all possible rotations
      if (unfound[k] && special.rotate_approx(k, j, transpose(ps.get_proper(i)), -ps.order(i))){ // -order checks all possible rotations
        one_type.push_back(k);
        unfound[k] = false;
      }
      status_update("Point equivalent to ",j," for symmetry ",i,":",one_type);
      // sort the equivalent points by their relative order for this operation
      // such that Rⁱj ≡ type_order[i]
      type_order.clear();
      type_order.insert(type_order.begin(), ps.order(i), special.size());
      type_unfound.clear(); type_unfound.insert(type_unfound.begin(), one_type.size(), true);
      for (int o=0; o<ps.order(i); ++o)
      for (size_t k=0; k<one_type.size(); ++k) if (type_unfound[k]){
        // if (special.rotate_approx(one_type[k], j, ps.get(i), o)){
        if (special.rotate_approx(one_type[k], j, transpose(ps.get_proper(i)), o)){
          type_order[o] = one_type[k];
          type_unfound[k] = false;
        }
      }
      // and store the sorted equivalent indices
      status_update("Which are sorted by their rotation order:",type_order);
      one_sym.push_back(type_order);
    }
    std::sort(one_sym.begin(), one_sym.end(),
      [&](std::vector<size_t>& a, std::vector<size_t>& b){
        size_t as = a.size() - std::count(a.begin(), a.end(), special.size());
        size_t bs = b.size() - std::count(b.begin(), b.end(), special.size());
        return as > bs;
      }
    );
    size_t keep_count;
    LQVec<double> pt0(this->outerlattice, 1u), pt1(this->outerlattice, 1u);
    for (size_t s=0; s<one_sym.size(); ++s) if (sym_unused[i]){
      // we need at least two equivalent points, ideally there will be the same number as the order
      if (one_sym[s].size() > static_cast<size_t>(std::count(one_sym[s].begin(), one_sym[s].end(), special.size())+1)){
        type_order = one_sym[s];
        status_update("Highest-multiplicity distinct point type:",type_order);
      // auto loc = std::find_if(one_sym.begin(), one_sym.end(), [&](std::vector<size_t>x){return x.size()==static_cast<size_t>(ps.order(i));});
      // if (loc != one_sym.end()){
      //   size_t idx = loc - one_sym.begin();
        // type_order = one_sym[idx];
        // grab the rotation stationary axis
        vec.set(0, ps.axis(i));
        for (size_t j=0, k=1; j<type_order.size(); ++j, k=(j+1)%type_order.size())
        if (sym_unused[i] && type_order[j]<special.size() && type_order[k]<special.size()){
          if (ps.order(i)>2){
            // hold the two special points in their own LQVec<double> (!not <int>!)
            pt0 = special.extract(type_order[j]);
            pt1 = special.extract(type_order[k]);
            // we have two plane normals to worry about:
            status_update("Stationary vector",vec.to_string(0),"and special points",pt0.to_string(0),"and",pt1.to_string(0));
            // find both cross products, remembering that we want normals
            // pointing *into* the wedge.
            nrm.resize(2);
            if ( dot(pt1, cross(vec.star(), pt0)).all_approx("<",0.) ){
              // the rotation is left handed, so swap the special points
              nrm.set(0, cross(vec.star(), pt1));
              nrm.set(1, cross(pt0, vec.star()));
            } else {
              nrm.set(0, cross(vec.star(), pt0));
              nrm.set(1, cross(pt1, vec.star()));
            }
            status_update("give normals:", nrm.to_string(";"));
            // now check that all special points are inside of the wedge defined by the normals
          } else {
            // order == 2, so only one normal to worry about:
            nrm = cross(vec.star(), special.extract(type_order[j]));
            // make sure we don't remove all points out of the plane containing
            // the rotation axis and the two special points
            if (dot(nrm, special).is_approx(">",0.).count_true() == 0)
              nrm *= -1; // switch the cross product
          }
          // check to make sure that using these normals do not remove all but a plane of points
          keep = keep_if(nrm, special);
          keep_count = keep.count_true();
          // We need at least three points (plus Γ) to define a polyhedron.
          // Also skip the extraction if we are keeping all points
          if (keep_count > 2 && keep_count < keep.size()){
            status_update("Keeping special points:\n",special.to_string(keep));
            special = special.extract(keep);
            sym_unused[i]=false;
          }
        }
      }
    }
  }

  // status_update("Remaining special points\n", special.to_string());

  // append the Γ point onto the list of special points:
  special.resize(special.size()+1u);
  for (size_t i=0; i<3u; ++i) special.insert(0, special.size()-1, i);
  // and then form the irreducible polyhedron by finding the convex hull
  // of these points: (which sets the wedge normals too)
  this->set_ir_vertices(special);

  /*If the full pointsymmetry has space inversion *and* there is only a single
    stationary rotation/rotoinversion axis the thus-far-found polyhedron has
    twice the volume of the true irreducible polyhedron.
    We need to split the polyhedron in half by a plane perpendicular to *one of*
    the full-Brillouin zone normal directions.
    This is easier for spacegroups with orthogonal basis vectors -- we can
    usually pick the stationary axis -- but is trickier if the basis vectors
    are not orthogonal.
  */
  if (fullps.has_space_inversion()){
    // ps only contains operations *with* a stationary axis (described in real space)
    LDVec<int> all_axes(this->outerlattice.star(), ps.axes());
    double ir_volume = this->ir_polyhedron.get_volume();
    double goal_volume = this->polyhedron.get_volume()/static_cast<double>(fullps.size());
    if (all_axes.is_unique().count_true() == 1u || approx_scalar(ir_volume, 2.0*goal_volume) ){
      status_update("Deal with -1 since there is only one stationary axis (or doubled volume for some other reason)");
      Polyhedron div;
      LQVec<double> bz_n = this->get_normals(), wg_n = this->get_ir_wedge_normals();
      ArrayVector<double> gamma(3u, 1u);
      for (size_t i=0; i<3u; ++i) gamma.insert(0, 0, i);
      for (size_t i=0; i<bz_n.size(); ++i){
        div = this->ir_polyhedron.divide(bz_n.extract(i).get_xyz(), gamma);
        if (approx_scalar(div.get_volume(), goal_volume)){
          // set div to be the ir_polyhedron
          this->ir_polyhedron = div;
          // add the new normal to wedge normals list
          wg_n.resize(wg_n.size()+1);
          wg_n.set(wg_n.size()-1, -bz_n.extract(i));
          this->set_ir_wedge_normals(wg_n);
          break;
        }
      }
      if (approx_scalar(this->ir_polyhedron.get_volume(), ir_volume)){
        status_update("Polyhedron volume still double expected.")
        bool proceed=true;
        // check for other dividing planes :/
        LQVec<double> cij(this->outerlattice, 1u);
        for (size_t i=0; proceed && i<bz_n.size()-1; ++i)
        for (size_t j=i+1; proceed && j<bz_n.size(); ++j){
          div = this->ir_polyhedron.divide(bz_n.cross(i,j).get_xyz(), gamma);
          if (approx_scalar(div.get_volume(), goal_volume)){
            this->ir_polyhedron = div;
            wg_n.resize(wg_n.size()+1);
            wg_n.set(wg_n.size()-1, -bz_n.cross(i,j));
            this->set_ir_wedge_normals(wg_n);
            proceed = false;
          }
        }
        // ArrayVector<double> xyz(3u, 3u);
        // for (size_t i=0; i<3u; ++i) for (size_t j=0; j<3u; ++j) xyz.insert(i==j?1.:0., i, j);
        // for (size_t i=0; i<xyz.size(); ++i){
        //   div = this->ir_polyhedron.divide(xyz.extract(i), gamma);
        //   if (approx_scalar(div.get_volume(), goal_volume)){
        //     // set the polyhedron
        //     this->ir_polyhedron = div;
        //     // set the wedge normals
        //     LQVec<double> ir_n = this->get_ir_normals();
        //     LQVec<double> ir_p = this->get_ir_points();
        //     LQVec<double> bz_p = this->get_points();
        //     ArrayVector<bool> not_bz(1u, 0u);
        //     for (size_t i=0; i<bz_n.size(); ++i){
        //       // check if the irBZ face point is on a first BZ zone face too
        //       not_bz = dot(bz_n.extract(i), ir_p - bz_p.extract(i)).is_approx("!=", 0.);
        //       ir_n = ir_n.extract(not_bz);
        //       ir_p = ir_p.extract(not_bz);
        //     }
        //     // the remaining irBZ faces are the irreducible reciprocal space wedge
        //     // which we store with the opposite sign in this->ir_wedge_normals
        //     this->set_ir_wedge_normals(-ir_n);
        //     break;
        //   }
        // }
      }
    }
  }
}


void BrillouinZone::wedge_triclinic(void){
  /* Assuming that this is spacegroup P -1 we have inversion symmetry only.
     We always want to find a set of planes that divides symmetry equivalent
     regions of space. With only -1 symmetry we have to make a choice as to
     *which* space inversion we care about.
     We could restrict ourselves to {x≥0, y, z}, {x, y≥0, z}, {x, y, z≥0} or
     their opposites, or any other subspace for which xyz ≥ 0.
     Equivalently then, we can resstrict ourselves to the subspace where
     x̂⋅(111) ≥ 0.
  */
  LQVec<double> nrm(this->outerlattice, 1u);
  nrm.set(0, std::array<double,3>({1,1,1}));
  this->set_ir_wedge_normals(nrm);
  this->irreducible_vertex_search();
}
