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
  LQVec<double> centres = this->get_points();
  LQVec<double> corners = this->get_vertices();
  LQVec<double> special = cat(centres, corners);
  // if this spacegroup is Rhombohedral, Monoclinic, or Triclinic
  // add the mid-facet-edge points:
  Pointgroup pg = this->outerlattice.get_pointgroup_object();
  // Holohedry h = pg.get_holohedry();
  // if (Holohedry::triclinic == h || Holohedry::monoclinic == h ||
  //     // some ofthe trigonal pointgroups are associated with rhombohedral spacegroups
  //     (Holohedry::trigonal == h && !sg.get_choice().compare("R")))
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
  // The first Brillouin zone volume
  double volume_goal = this->get_polyhedron().get_volume();
  // Divided by the number of pointgroup symmetry operations is the Goal volume
  volume_goal /= fullps.size();
  // Now restrict the symmetry operations to those with order > 1.
  PointSymmetry ps = fullps.higher(1);
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
  ps.permute(perm); // ps now sorted with shortest stationary axis first

  ArrayVector<bool> keep(1u, 0u);
  if (fullps.has_space_inversion()){
    // ps does not include 1 or -1, so all operations have stationary axes
    LQVec<int> all_axes(this->outerlattice, ps.axes());
    if (all_axes.is_unique().count_true() == 1u){
      status_update("Deal with -1 since there is only one stationary axis");
      // use the stationary axis as the normal of a mirror plane
      keep = dot(all_axes.extract(0), special).is_approx(">=", 0.);
      special = special.extract(keep);
    }
  }

  status_update("Deal with 2-fold axes along highest-symmetry directions first");

  std::vector<bool> sym_unused(ps.size(), true);
  // deal with any two-fold axes along êᵢ first:
  LQVec<int> eis(this->outerlattice, 9u);
  LQVec<int> vec(this->outerlattice, 1u);// must be int since ps.axis returns array<int,3>
  LQVec<double> nrm(this->outerlattice, 1u), pts(this->outerlattice, 2u);
  ArrayVector<bool> is_ei(1u, 9u);
  std::vector<std::array<int,3>> eiv{{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,-1,0},{1,0,1},{0,1,1},{1,0,-1},{0,1,-1}};
  for (size_t i=0; i<9u; ++i) eis.set(i, eiv[i]);

  for (size_t i=0; i<ps.size(); ++i) if (ps.order(i)==2){
    vec.set(0, ps.axis(i));
    is_ei = norm(cross(eis, vec)).is_approx("==", 0.);
    status_update("Is 2-fold axis ",i," ei?: ",is_ei.to_string(" "));
    if (is_ei.count_true() == 1u){
      if (is_ei.getvalue(0)) // (100) → n = (001)×(100) { =(010) if cubic}
        nrm.set(0, eis.cross(2,0));
      if (is_ei.getvalue(1)) // (010) → n = (100)×(010) { =(001) if cubic}
        nrm.set(0, eis.cross(0,1));
      if (is_ei.getvalue(2)) // (001) → n = (010)×(001) { =(100) if cubic}
        nrm.set(0, eis.cross(1,2));
      if (is_ei.getvalue(3)) // (110) → n = (110)×(001) { =(1̄10) if cubic}
        nrm.set(0, eis.cross(3,2));
      if (is_ei.getvalue(4)) // (1̄10) → n = (001)×(1̄10) { =(110) if cubic}
        nrm.set(0, eis.cross(2,4));
      if (is_ei.getvalue(5)) // (101) → n = (010)×(101) { =(10̄1) if cubic}
        nrm.set(0, eis.cross(1,5));
      if (is_ei.getvalue(6)) // (011) → n = (011)×(100) { =(01̄1) if cubic}
        nrm.set(0, eis.cross(6,0));
      if (is_ei.getvalue(7)) // (1̄0̄1) → n = (10̄1)×(010) { =(101) if cubic}
        nrm.set(0, eis.cross(7,1));
      if (is_ei.getvalue(8)) // (01̄1) → n = (100)×(01̄1) { =(011) if cubic}
        nrm.set(0, eis.cross(0,8));
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
      if (unfound[k] && special.rotate_approx(k, j, ps.get_proper(i), -ps.order(i))){ // -order checks all possible rotations
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
        if (special.rotate_approx(one_type[k], j, ps.get_proper(i), o)){
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
            pts.set(0, special.extract(type_order[j]));
            pts.set(1, special.extract(type_order[k]));
            if ( dot(eis.extract(1), cross(vec, eis.extract(0))).all_approx("<",0.) ){
              // the rotation is left handed, so swap the special points
              pts.set(1, special.extract(type_order[k]));
              pts.set(0, special.extract(type_order[j]));
            }
            // we have two plane normals to worry about:
            status_update("Stationary vector\n",vec.to_string(),"and special points\n",pts.to_string());
            // find both cross products, remembering that we want normals
            // pointing *into* the wedge.
            nrm.resize(2);
            nrm.set(0, cross(1.0*vec, pts.extract(0)));
            nrm.set(1, cross(pts.extract(1), 1.0*vec));
            status_update("give normals:\n", nrm.to_string());
            // now check that all special points are inside of the wedge defined by the normals
          } else {
            // order == 2, so only one normal to worry about:
            nrm = cross(vec, special.extract(type_order[j]));
            // make sure we don't remove all points out of the plane containing
            // the rotation axis and the two special points
            if (dot(nrm, special).is_approx(">",0.).count_true() == 0)
              nrm *= -1; // switch the cross product
          }
          // check to make sure that using these normals do not remove all but a plane of points
          keep = keep_if(nrm, special);
          if (keep.count_true() < keep.size()){
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
