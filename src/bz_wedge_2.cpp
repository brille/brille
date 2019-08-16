static size_t number_of_polyhedra_to_check(size_t n){
  size_t total=0u, partial=1u;
  for (size_t i=0; i<n; ++i){
    total += partial;
    partial *= n-i;
    partial /= i+1;
  }
  return total;
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

  std::cout << std::to_string(special.size()) << " total special points" << std::endl;

  // Since we're assuming inversion symmetry (rectified later if necessary),
  // for every Friedel-pair of points we can immediately discard one
  // (half of all points for a crystallographic spacegroup should be discarded)
  std::vector<bool> isok(special.size(), true);
  for (size_t i=0; i<special.size()-1; ++i) if (isok[i])
  for (size_t j=i+1; j<special.size(); ++j) if (isok[j])
  if (special.vector_approx(i,j,"*",-1.0)) isok[j]=false;
  // isok[i] == true → point i is the first half of a Friedel-pair or a loner
  special = special.extract(isok);
  std::cout << "Reduced to " << std::to_string(special.size()) << " after Friedel-pair reduction" << std::endl;

  // Grab the pointgroup symmetry operations with order>1
  PointSymmetry ps = this->outerlattice.get_pointgroup_symmetry(this->time_reversal).higher(1);
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

  std::cout << "Of which " << sym_equiv_idx.size() << " points are unique by symmetry" << std::endl;
  /* We want all combinations of two special points, which when combined with
  the Γ point, represent all possible symmetry planes in the Brillouin zone  */
  size_t count=0, expected = (special.size()*(special.size()-1))>>1;
  LQVec<double> symplanes(this->outerlattice, expected);
  for (size_t i=0; i<special.size()-1; ++i)
  for (size_t j=i+1; j<special.size(); ++j)
  symplanes.set(count++, special.cross(i,j));
  if ( count != expected ){
    msg = "Found " + std::to_string(count) + " normals but expected " + std::to_string(expected);
    throw std::runtime_error(msg);
  }
  // keep just the unique symmetry plane normals
  symplanes = symplanes.extract(symplanes.is_unique());
  size_t to_check = number_of_polyhedra_to_check(symplanes.size());
  std::cout << "With " << std::to_string(symplanes.size()) << " unique plane ";
  std::cout << "normals, we might need to check " << std::to_string(to_check);
  std::cout << std::endl;


}
