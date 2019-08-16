static size_t number_of_polyhedra_to_check(size_t n){
  size_t total=0u, partial=1u;
  for (size_t i=0; i<n; ++i){
    total += partial;
    partial *= n-i;
    partial /= i+1;
  }
  return total;
}
template<class T>
static bool next_set(const std::vector<std::vector<T>>& equiv, std::vector<size_t>& idx, std::vector<T>& set){
  bool goon = true;
  size_t n = idx.size();
  while (goon && n-- > 0u){
    if (++idx[n] < equiv[n].size())
      goon = false;
    else
      idx[n] = 0u;
  }
  for (size_t i=0; i<idx.size(); ++i) set[i] = equiv[i][idx[i]];
  return goon; // overflow if goon is still true
}

void BrillouinZone::wedge_brute_force(void){
  debug_exec(std::string msg;)
  // get the characteristic points of the first Brillouin zone
  LQVec<double> centres = this->get_points();
  LQVec<double> corners = this->get_vertices();
  LQVec<double> special = cat(centres, corners);
  // if this spacegroup is Rhombohedral, Monoclinic, or Triclinic
  // add the mid-facet-edge points:
  Spacegroup sg = this->outerlattice.get_spacegroup_object();
  Pointgroup pg = this->outerlattice.get_pointgroup_object();
  Holohedry h = pg.get_holohedry();
  if (Holohedry::triclinic == h || Holohedry::monoclinic == h ||
      // some ofthe trigonal pointgroups are associated with rhombohedral spacegroups
      (Holohedry::trigonal == h && !sg.get_choice().compare("R")))
    special = cat(special, this->get_half_edges());

  // Since we're assuming inversion symmetry (rectified later if necessary),
  // for every Friedel-pair of points we can immediately discard one
  // (half of all points for a crystallographic spacegroup should be discarded)
  std::vector<bool> isok(special.size(), true);
  for (size_t i=0; i<special.size()-1; ++i) if (isok[i])
  for (size_t j=i+1; j<special.size(); ++j) if (isok[j])
  if (special.vector_approx(i,j,"*",-1.0)) isok[j]=false;
  // isok[i] == true → point i is the first half of a Friedel-pair or a loner
  special = special.extract(isok);

  // Grab the pointgroup symmetry operations
  PointSymmetry ps = this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
  // The first Brillouin zone volume
  double volume_goal = this->get_polyhedron().get_volume();
  // Divided by the number of pointgroup symmetry operations is the Goal volume
  volume_goal /= this->outerlattice.get_pointgroup_symmetry(this->time_reversal).size();
  // Now restrict the symmetry operations to those with order > 1.
  ps = ps.higher(1);
  // And use them to find which special points are symmetry equivalent
  std::vector<std::vector<size_t>> sym_equiv_idx;
  std::vector<size_t> one_type;
  std::vector<bool> unfound(special.size(), true);
  for (size_t i=0; i<special.size(); ++i) if (unfound[i]) {
    one_type.clear(); // empty the contents, reserving the memory
    one_type.push_back(i);
    unfound[i] = false; // not necessary, likely
    for (size_t k=0; k<ps.size(); ++k)
    for (size_t j=i+1; j<special.size(); ++j)
    if (unfound[j] && special.rotate_approx(i,j,ps.get(k))){ // special[i] ≡ R[k]*special[j]
      one_type.push_back(j);
      unfound[j] = false;
    }
    sym_equiv_idx.push_back(one_type);
  }
  std::sort(sym_equiv_idx.begin(), sym_equiv_idx.end(),
    [](std::vector<size_t>& a, std::vector<size_t>& b){
      return a.size() > b.size(); // descending sort the vectors by length
    }
  );
  std::cout << "sym_equiv_idx:" << std::endl;
  for (auto x: sym_equiv_idx){
    for (auto y: x) std::cout << " " << y;
    std::cout << std::endl;
  }

  /* We want to look for polyhedra that are defined by planes with normals that
  are the cross products between vectors from the Γ point to any pair of high
  symmetry points in the first Brillouin zone. We are trying to limit how many
  such planes we might need to consider by only considering cross products
  between symmetry distinct points, but this omits the possibility of normals
  defined from symmetry-equivalent points. At least in the case of positive
  isometry symmetry opperations (rotations only, no rotoinversions) the cross
  product between any two symmetry-equivalent points *is* the rotation axis.
  So if we also consider polyhedra defined by some number of rotation planes we
  might find a valid solution sooner. */

  // The nfolds method finds unique rotation-axis symmetry operations
  // of a minimum order (default 0) and returns the highest order operation
  // with each axies. We've already restricted ps to >1 order operations only,
  // so here we're abusing nfolds to find just the unique rotation axes.
  LQVec<double> rot_normals(this->outerlattice, ps.nfolds().axes() );

  debug_exec(size_t total_checks = 0;)
  double check_volume=-1;
  size_t Ndistinct = sym_equiv_idx.size(), count = 0;
  size_t expected = (Ndistinct*(Ndistinct-1))>>1;
  LQVec<double> symplanes(this->outerlattice, expected);
  LQVec<double> normals(this->outerlattice, 0u);
  // indexing for the Ndistinct symmetry-distinct points
  std::vector<size_t> distinct_idx(Ndistinct,0u);
  std::vector<size_t> special_idx;
  for (auto x: sym_equiv_idx) special_idx.push_back(x[0]);

  for (size_t i=0; i<rot_normals.size(); ++i){
    while (!approx_scalar(check_volume, volume_goal)){
      // Find this set of symmetry plane normals:
      count = 0u;
      for (size_t i=0; i<Ndistinct-1; ++i)
      for (size_t j=i+1; j<Ndistinct; ++j)
      symplanes.set(count++, special.cross(special_idx[i], special_idx[j]));
      debug_exec( if (count!=expected) status_update("Found ",count," normals but expected ",expected); )

      // copy the symmetry plane normals and add one rotation-axis
      normals = cat(symplanes, rot_normals.extract(i));
      // keep only unqiue plane normals -- maybe we should skip any rot_normals[i] which is in symplanes instead?
      normals = normals.extract(normals.is_unique());
      // From these construct a Brillouin-zone-limited polyhedron, and get its volume
      check_volume = this->ir_wedge_is_ok(normals) ? this->get_ir_polyhedron().get_volume() : 0.0;
      debug_exec(++total_checks;)
      // grab the next set of special indices, if distinct_idx overflows, break out
      if (next_set(sym_equiv_idx, distinct_idx, special_idx)) break;
    }
    // if we have found the right volue, break out of the for loop
    if (approx_scalar(check_volume, volume_goal)) break;
  }

  
  // reset special_idx (this shouldn't be necessary, since we overflow back to this condition)
  special_idx.clear();
  for (auto x: sym_equiv_idx) special_idx.push_back(x[0]);
  // to handle sign permutations of the cross products
  std::vector<std::vector<bool>> all_signs(expected, {false, true});
  all_signs[0] = {false}; // we can skip [-1,*] as it is equivalent to [1,-*]
  std::vector<size_t> signs_idx(expected, 0u);
  std::vector<bool> signs(expected);
  // run through all possible combinations until we find one with the right volume
  while (!approx_scalar(check_volume, volume_goal)){
    // Find this set of symmetry plane normals:
    count = 0u;
    for (size_t i=0; i<Ndistinct-1; ++i)
    for (size_t j=i+1; j<Ndistinct; ++j)
    symplanes.set(count++, special.cross(special_idx[i], special_idx[j]));
    debug_exec( if (count!=expected) status_update("Found ",count," normals but expected ",expected); )
    /* Each cross product could be done with the opposite handedness. So we
       need to try the 2ᴺ combinations of ±1 for each :( */
    for (size_t i=0; i<expected; ++i) signs[i] = all_signs[i][0];
    do {
      // copy the symmetry plane normals
      normals = symplanes;
      // so that we can flip the signs of some
      for (size_t i=0; i<expected; ++i) if (signs[i]) normals.set(i, -1*normals.extract(i));
      // Make sure we only keep unique plane normals:
      normals = normals.extract(normals.is_unique());
      // From these construct a Brillouin-zone-limited polyhedron, and get its volume
      check_volume = this->ir_wedge_is_ok(normals) ? this->get_ir_polyhedron().get_volume() : 0.0;
      debug_exec(++total_checks;)
      // if we have found the right volue, break out
      if (approx_scalar(check_volume, volume_goal)) break;
    } while (!next_set(all_signs, signs_idx, signs));
    // grab the next set of special indices, if distinct_idx overflows, break out
    if (next_set(sym_equiv_idx, distinct_idx, special_idx)) break;
  }


  if (!approx_scalar(check_volume, volume_goal)){
    std::string err = "Failed to find a polyhedron with the correct volume for";
    err += " the irreducible Brillouin zone";
    debug_exec(err += " in " + std::to_string(total_checks) + " attempts";)
    err += ".\nLast attempt yielded " + std::to_string(this->get_ir_polyhedron().get_volume());
    err += " when we wanted " + std::to_string(volume_goal);
    // throw std::runtime_error(err);
    std::cout << err << std::endl;
  } else {
    debug_exec(\
      msg= "A polyhedron with the same volume as the irreducible ";\
      msg += "Brillouin zone has been found after " + std::to_string(total_checks);\
      msg += " attempts.";\
    )
    status_update(msg);
  }
}
