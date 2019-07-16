#include "bz.h"


void BrillouinZone::set_vertices(const ArrayVector<double> newverts){
  if (newverts.numel()!=3u) throw std::runtime_error("BrillouinZone objects only take 3-D vectors for vertices");
  this->vertices = newverts;
}
void BrillouinZone::set_faces(const ArrayVector<int> newfaces){
  if (newfaces.numel()!=3u) throw std::runtime_error("BrillouinZone objects only take 3-D vectors for faces");
  this->faces = newfaces;
}
void BrillouinZone::set_faces_per_vertex(const ArrayVector<int> newfpv){
  if (newfpv.numel()!=3u) throw std::runtime_error("BrillouinZone vertices are the intersections of three faces each");
  this->faces_per_vertex = newfpv;
}
void BrillouinZone::set_ir_wedge_normals(const LQVec<double>& x){
  bool already_same = this->outerlattice.issame(x.get_lattice());
  LQVec<double> xp(this->outerlattice);
  PrimitiveTransform PT(this->outerlattice.get_hall());
  bool transform_needed = ( PT.does_anything() && this->lattice.issame(x.get_lattice()) );
  if (!(already_same || transform_needed))
    throw std::runtime_error("ir_wedge_normals must be in the standard or primitive lattice used to define the BrillouinZone object");
  if (transform_needed)  xp = transform_from_primitive(this->outerlattice,x);
  const LQVec<double> & xsl = transform_needed ? xp : x;
  this->ir_wedge_normals = x.get_hkl();
}
void BrillouinZone::set_ir_vertices(const LQVec<double>& x){
  bool already_same = this->outerlattice.issame(x.get_lattice());
  LQVec<double> xp(this->outerlattice);
  PrimitiveTransform PT(this->outerlattice.get_hall());
  bool transform_needed = ( PT.does_anything() && this->lattice.issame(x.get_lattice()) );
  if (!(already_same || transform_needed))
    throw std::runtime_error("ir_face_normals must be in the standard or primitive lattice used to define the BrillouinZone object");
  if (transform_needed)  xp = transform_from_primitive(this->outerlattice,x);
  const LQVec<double> & xsl = transform_needed ? xp : x;
  this->ir_vertices = x.get_hkl();
}
void BrillouinZone::set_ir_face_normals(const LQVec<double>& x){
  bool already_same = this->outerlattice.issame(x.get_lattice());
  LQVec<double> xp(this->outerlattice);
  PrimitiveTransform PT(this->outerlattice.get_hall());
  bool transform_needed = ( PT.does_anything() && this->lattice.issame(x.get_lattice()) );
  if (!(already_same || transform_needed))
    throw std::runtime_error("ir_face_normals must be in the standard or primitive lattice used to define the BrillouinZone object");
  if (transform_needed)  xp = transform_from_primitive(this->outerlattice,x);
  const LQVec<double> & xsl = transform_needed ? xp : x;
  this->ir_face_normals = x.get_hkl();
}
void BrillouinZone::set_ir_face_points(const LQVec<double>& x){
  bool already_same = this->outerlattice.issame(x.get_lattice());
  LQVec<double> xp(this->outerlattice);
  PrimitiveTransform PT(this->outerlattice.get_hall());
  bool transform_needed = ( PT.does_anything() && this->lattice.issame(x.get_lattice()) );
  if (!(already_same || transform_needed))
    throw std::runtime_error("ir_face_normals must be in the standard or primitive lattice used to define the BrillouinZone object");
  if (transform_needed)  xp = transform_from_primitive(this->outerlattice,x);
  const LQVec<double> & xsl = transform_needed ? xp : x;
  this->ir_face_points = x.get_hkl();
}
// void BrillouinZone::set_ir_wedge_normals(const ArrayVector<double> newplanes){
  //   if (newplanes.numel()!=3u) throw std::runtime_error("BrillouinZone objects only take 3-D vectors for plane normals.");
  //   this->ir_wedge_normals = newplanes;
  // }
// void BrillouinZone::set_ir_vertices(const ArrayVector<double>& x){
//   if (x.numel()!=3u) throw std::runtime_error("BrillouinZone object only take 3-D vectors for vertices.");
//   this->ir_vertices = x;
// }
// void BrillouinZone::set_ir_face_normals(const ArrayVector<double>& x){
//   if (x.numel()!=3u) throw std::runtime_error("BrillouinZone objects only take 3-D vectors for face normals");
//   this->ir_face_normals = x;
// }
// void BrillouinZone::set_ir_face_points(const ArrayVector<double>& x){
//   if (x.numel()!=3u) throw std::runtime_error("BrillouinZone objects only take 3-D vectors for face points");
//   this->ir_face_points = x;
// }
void BrillouinZone::set_ir_faces_per_vertex(const ArrayVector<int>& x){
  if (x.numel()!=3u) throw std::runtime_error("BrillouinZone vertices are the intersections of three faces");
  this->ir_faces_per_vertex = x;
}

LQVec<double> BrillouinZone::get_vertices(void) const {
  LQVec<double> lqverts(this->lattice,this->vertices);
  if (this->isprimitive())
    lqverts = transform_from_primitive(this->outerlattice,lqverts);
  return lqverts;
}
LQVec<double> BrillouinZone::get_primitive_vertices(void) const {
  return LQVec<double>(this->lattice, this->vertices);
}
LQVec<int>    BrillouinZone::get_faces   (void) const {
  LQVec<int> lqfaces(this->lattice,this->faces);
  if (this->isprimitive())
    lqfaces = transform_from_primitive(this->outerlattice,lqfaces);
  return lqfaces;
}
LQVec<int> BrillouinZone::get_primitive_faces(void) const {
  return LQVec<int>(this->lattice, this->faces   );
}
ArrayVector<int> BrillouinZone::get_faces_per_vertex(void) const {
  ArrayVector<int> out = this->faces_per_vertex; // make sure we return a copy, not the internal object
  return out;
}
//
LQVec<double> BrillouinZone::get_ir_wedge_normals(void) const {
  return LQVec<double>(this->outerlattice, this->ir_wedge_normals);
}
LQVec<double> BrillouinZone::get_ir_vertices() const {
  return LQVec<double>(this->outerlattice, this->ir_vertices);
}
LQVec<double> BrillouinZone::get_ir_face_normals() const {
  return LQVec<double>(this->outerlattice, this->ir_face_normals);
}
LQVec<double> BrillouinZone::get_ir_face_points() const {
  return LQVec<double>(this->outerlattice, this->ir_face_points);
}
ArrayVector<int> BrillouinZone::get_ir_faces_per_vertex() const {
  return ArrayVector<int>(this->ir_faces_per_vertex);
}
//
LQVec<double> BrillouinZone::get_primitive_ir_wedge_normals(void) const {
  LQVec<double> lqwn(this->outerlattice, this->ir_wedge_normals);
  if (this->isprimitive())
    lqwn = transform_to_primitive(this->outerlattice, lqwn);
  return lqwn;
}
LQVec<double> BrillouinZone::get_primitive_ir_vertices() const {
  LQVec<double> x(this->outerlattice, this->ir_vertices);
  if (this->isprimitive())
    x = transform_to_primitive(this->outerlattice, x);
  return x;
}
LQVec<double> BrillouinZone::get_primitive_ir_face_normals() const {
  LQVec<double> x(this->outerlattice, this->ir_face_normals);
  if (this->isprimitive())
    x = transform_to_primitive(this->outerlattice, x);
  return x;
}
LQVec<double> BrillouinZone::get_primitive_ir_face_points() const {
  LQVec<double> x(this->outerlattice, this->ir_face_points);
  if (this->isprimitive())
    x = transform_to_primitive(this->outerlattice, x);
  return x;
}
//
void BrillouinZone::print() const {
  printf("BrillouinZone with %u vertices and %u faces\n",this->vertices_count(),this->faces_count());
}

void BrillouinZone::wedge_search(const int time_reversal){
  /*
  The Pointgroup symmetry information comes from, effectively, spglib which
  has all rotation matrices defined in the conventional unit cell -- which is
  our `outerlattice`. Consequently we must work in the outerlattice here.
  */
  PointSymmetry psym = get_pointgroup_symmetry(this->outerlattice.get_hall(), time_reversal);
  // extract the symmetry operations with rotation axes: (not 1 or ̄1)
  std::vector<std::array<int,9>> rotations;
  for (size_t i=0; i<psym.size(); ++i)
    if (rotation_order(psym.get(i)) > 1)
      rotations.push_back(psym.getarray(i));
  // sort the rotations by their orders
  std::sort(rotations.begin(), rotations.end(), [](std::array<int,9> a, std::array<int,9> b){
    return rotation_order(a.data()) < rotation_order(b.data());
  });
  std::vector<int> orders;
  for (std::array<int,9> R: rotations) orders.push_back(rotation_order(R.data()));
  int max_order=0;
  for (int o: orders) if (o>max_order) max_order = o;
  std::cout << "Rotation orders:";
  for (int o: orders) std::cout << " " << o;
  std::cout << std::endl;
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
      y.insert(static_cast<double>(axis_plane_vecs[2][j]),i,j);
    }
  }
  // check for zero-length rotation axes *before* trying to normalize them
  ArrayVector<bool> nonzero = norm(u).is_approx(">",0.);
  u = u.extract(nonzero);
  x = x.extract(nonzero);
  // y = y.extract(nonzero);
  y = cross(u,x); // in case the basis is not orthogonal
  // make sure we're dealing with unit vectors
  // u = u/norm(u);
  // x = x/norm(x);
  // y = y/norm(y);

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
  if (unique_idx.size() != u.is_unique().count_true()) //this should only be in for debugging
  throw std::runtime_error("Unique count is off.");

  // Verify that we have a right-handed rotation axis for each matrix
  LQVec<double> Rx(this->outerlattice, 1u);
  for (size_t i=0; i<rotations.size(); ++i){
    // Calculate Rx
    multiply_matrix_vector(Rx.datapointer(0), rotations[i].data(), x.datapointer(i));
    // y ≡ u × x, so if rotated x points away from y we need to flip u.
    if (dot(Rx, y.get(i)).getvalue(0) < 0) u.set(i, -u.get(i));
  }

  // Find ûᵢ⋅ûⱼ
  double* ui_dot_uj = new double[u.size()*u.size()]();
  for (size_t i=0; i<u.size(); ++i) for (size_t j=0; j<u.size(); ++j) ui_dot_uj[i*u.size()+j] = u.dot(i,j)/10;

  std::cout << "unique(u):\n";
  for (size_t i=0; i<u.size(); ++i)
  if (u_equiv_idx.getvalue(i)==i)
  std::cout << u.to_string(i) << std::endl;

  std::cout << "dot(ui, uj)\n";
  for (size_t i=0; i<u.size(); ++i){
    if (u_equiv_idx.getvalue(i)==i){
      for (size_t j=0; j<u.size(); ++j)
      if (u_equiv_idx.getvalue(j)==j)
      std::cout << " " << (approx_scalar(ui_dot_uj[i*u.size()+j],0.)?"0":"x");
      std::cout << std::endl;
    }
  }
  // Two rotations with ûᵢ⋅ûⱼ ≠ 0 should have (Rⁿv̂ᵢ)⋅(Rᵐv̂ⱼ) = 1 for some n,m.
  // To start with, assume that n and m are 0 and find the v̂ᵢ such
  // that v̂ᵢ⋅(ûᵢ×ûⱼ)=1
  std::vector<bool> handled;
  for (size_t i=0; i<u.size(); ++i) handled.push_back(false);
  // we need a place to stash the cross product of uᵢ and uⱼ
  LQVec<double> n(this->outerlattice, 1u);
  // why not just stash away the v̂ᵢ now?
  LQVec<double> v(x); // copy x
  for (size_t i=0; i<u.size()-1; ++i)
  if(u_equiv_idx.getvalue(i)==i)
  for (size_t j=i+1; j<u.size(); ++j)
  if(u_equiv_idx.getvalue(j)==j){
    if (!approx_scalar(ui_dot_uj[i*u.size()+j], 0.0)){
      if (!handled[i] && handled[j]){
        n = u.cross(i,j); // this is why we need to skip non-unique axes.
        // v.set(i, n/norm(n));
        v.set(i, n);
        // we don't want to just set vᵢ to vⱼ in case of, e.g., [100],[010],[111]
        // where [100]×[111] is [0̄11] and [010]×[111] is [10̄1]
        if (v.dot(i,j) < 0) v.set(i, -v.get(i));
        handled[i] = true;
      }
      if (!handled[i] && !handled[j]){
        n = u.cross(i,j); // this is why we need to skip non-unique axes.
        // v.set(i, n/norm(n));
        v.set(i, n);
        handled[i] = true;
      }
      if (handled[i] && !handled[j]){
        // protect against v ∥ u
        if (!approx_scalar(norm(cross(v.get(i), u.get(j))).getvalue(0)/10, 0.0))
          v.set(j, v.get(i));
        else{
          n = u.cross(i,j);
          // v.set(i, n/norm(n));
          v.set(i, n);
          if (v.dot(i,j) < 0) v.set(i, -v.get(i));
        }
        handled[j] = true;
      }
    }
  }
  delete[] ui_dot_uj;
  /* for the case where three rotations have ûᵢ⋅ûⱼ=0, ûᵢ⋅ûₖ≠0, ûⱼ⋅ûₖ≠0
    the above will set vᵢ=v̂ₖ=ûᵢ×ûₖ and v̂ⱼ=ûⱼ×ûₖ and we assume (hope?) that
    there exists some Rₖⁿv̂ₖ=ûⱼ×ûₖ.     */
  // Now go back through and copy the unique (a,b) values to the non-uniques.
  for (size_t i=0; i<u.size(); ++i)
  if (u_equiv_idx.getvalue(i)!=i) v.set(i,v.get(u_equiv_idx.getvalue(i)));

  std::cout << "u_equiv_idx = " << u_equiv_idx.to_string(" ") << std::endl;
  std::cout << "u =\n" << u.to_string();
  std::cout << "v =\n" << v.to_string();

  // We now have for every rotation a consistent vector in the rotation plane.
  // It's time to use ̂u, ̂v, R, and order(R) to find/add the wedge normals
  LQVec<double> normals(this->outerlattice, 2*rotations.size());
  LQVec<double> vj(this->outerlattice, max_order);
  size_t total_found=0;
  bool accepted=false;
  int order;
  for (size_t j=0; j<rotations.size(); ++j){
    order = orders[j];
    std::cout << "r, u: " << order << " (" << u.to_string(j,") \n");

    accepted=false;
    if (2==order){
      // the single normal version of add_wedge_normal_check allows for the
      // possibility of either n or -n being added, so we don't need to check
      // both u×v and u×Rv (as Rv = -v for a 2-fold axis).
      accepted = this->wedge_normal_check(cross(u.get(j), v.get(j)), normals, total_found); // increments total_found if adding was a success
      //// but a 2-fold axis could always(?) be folded 90 degrees away too(?)
      if (!accepted)
        accepted = this->wedge_normal_check(v.get(j), normals, total_found);
    } else {
      vj.set(0, v.get(j));
      for (int k=1; k<order; ++k)
        multiply_matrix_vector(vj.datapointer(k), rotations[j].data(), vj.datapointer(k-1));
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
  std::cout << "wedge_search finished" << std::endl;
  if (total_found == 1){
    this->set_ir_wedge_normals(normals.get(0));
    this->irreducible_vertex_search();
  }
  // otherwise set_ir_wedge_normals and irreducible_vertex_search are called by add_wedge_normal_check
}
bool BrillouinZone::wedge_normal_check(const LQVec<double>& n, LQVec<double>& normals, size_t& num){
  for(int i=0; i<30; ++i) std::cout << " ";
  std::cout << "Considering " << n.to_string(0,"... ");
  if (norm(n).all_approx(0.0)){
    std::cout << "rejected; zero-length." << std::endl;
    return false;
  }
  if (num==0){
    std::cout << "accepted; first normal." << std::endl;
    normals.set(num, n.get(0));
    num=num+1;
    return true;
  }
  if (norm(cross(normals.first(num), n)).any_approx("==",0.)){
    std::cout << "rejected; already present" << std::endl;
    return false;
  }
  normals.set(num,  n.get(0));
  if (this->ir_wedge_is_ok(normals.first(num+1))){
    std::cout << "accepted" << std::endl;
    num=num+1;
    return true;
  }
  normals.set(num, -n.get(0));
  if (this->ir_wedge_is_ok(normals.first(num+1))){
    std::cout << "accepted (*-1)" << std::endl;
    num=num+1;
    return true;
  }
  std::cout << "rejected; would cause null wedge" << std::endl;
  return false;
}
bool BrillouinZone::wedge_normal_check(const LQVec<double>& n0, const LQVec<double>& n1, LQVec<double>& normals, size_t& num){
  for(int i=0; i<30; ++i) std::cout << " ";
  std::cout << "Considering " << n0.to_string(0," and ") << n1.to_string(0,"... ");
  bool p0=false, p1=false;
  if (num>0){
    p0 = norm(cross(normals.first(num), n0)).any_approx("==",0.);
    p1 = norm(cross(normals.first(num), n1)).any_approx("==",0.);
  }
  if (p0 && p1){
    std::cout << "rejected; already present" << std::endl;
    return false;
  }
  if (num>0 && dot(normals.first(num)/norm(n0),n0/norm(n0)).any_approx("==",1.)){
    normals.set(num,  n1.get(0));
    if (this->ir_wedge_is_ok(normals.first(num+1))){
      std::cout << "accepted n1 (n0 already present)" << std::endl;
      num=num+1;
      return true;
    }
    std::cout << "rejected; n0 already present, n1 causes null wedge" << std::endl;
    return false;
  }
  if (num>0 && dot(normals.first(num)/norm(n1),n1/norm(n1)).any_approx("==",1.)){
    normals.set(num,  n0.get(0));
    if (this->ir_wedge_is_ok(normals.first(num+1))){
      std::cout << "accepted n0 (n1 already present)" << std::endl;
      num=num+1;
      return true;
    }
    std::cout << "rejected; n1 already present, n0 causes null wedge" << std::endl;
    return false;
  }
  if (num>0 && (p0 || p1)){
    std::cout << "rejected; inverse of one already present" << std::endl;
    return false;
  }
  normals.set(num,   n0.get(0));
  normals.set(num+1, n1.get(0));
  if (this->ir_wedge_is_ok(normals.first(num+2))){
    std::cout << "accepted +n0 & +n1" << std::endl;
    num=num+2;
    return true;
  }
  std::cout << "rejected; would cause null wedge" << std::endl;
  return false;
}
bool BrillouinZone::ir_wedge_is_ok(const LQVec<double>& normals){
  this->set_ir_wedge_normals(normals);
  this->irreducible_vertex_search();
  LQVec<double> v = this->get_ir_vertices();
  if (0==v.size()) return false; // with no vertices there is also no volume
  // std::cout << std::to_string(v.size());
  // std::cout << " IR vertices:\n" << v.to_string();
  ArrayVector<bool> non_zero = norm(v).is_approx(">",0.0); // is_approx(">",0) ≡ !≈ 0

  std::vector<size_t> idx;
  for (size_t i=0; i<v.size(); ++i) if (non_zero.getvalue(i)) idx.push_back(i);
  size_t nzNo = idx.size();

  // the number of unique combinations of 3-vertices is the number that we need to check
  size_t ntocheck=0;
  // 0.5*∑ᵢ₌₁ᴺ⁻¹ i²-i  -- we can skip 1 and the last value of our for loop is N-1
  for (int i=2; i<nzNo; ++i) ntocheck += (i*(i-1))>>1;

  // for all combinations of three vertices (and ⃗0) find the volume of the
  // tetrahedron they describe. The volume of a tetrahedron with one vertex at
  // ⃗0 and the other three at ⃗a, ⃗b, and ⃗c is ( ⃗a ⋅ ⃗b × ⃗c )/6.
  // Since we're just checking if there are any non-zero-volume tetrahedra
  // just checking the dot and cross products is probably enough.
  std::vector<double> volumes;
  volumes.reserve(ntocheck);
  for (int i=0; i<nzNo-2; ++i)
  for (int j=i+1; j<nzNo-1; ++j)
  for (int k=j+1; k<nzNo; ++k)
  volumes.push_back(dot(v.get(idx[i]), v.cross(idx[j], idx[k])).getvalue(0));

  double volume = 0;
  for (auto abc: volumes) if(!approx_scalar(abc, 0.)) volume += std::abs(abc);

  return !approx_scalar(volume, 0.0);
}

std::vector<std::vector<int>> fpv_to_vpf(LQVec<double>& verts, ArrayVector<int>& fpv){
  // we need to go from faces_per_vertex to vertices_per_face
  size_t max_face=0;
  for (size_t i=0; i<fpv.size(); ++i) for (size_t j=0; j<fpv.numel(); ++j) if (fpv.getvalue(i,j)>max_face) max_face = fpv.getvalue(i,j);
  std::vector<bool> present(max_face+1);
  for (size_t i=0; i<present.size(); ++i) present[i]=false;
  for (size_t i=0; i<fpv.size(); ++i) for (size_t j=0; j<fpv.numel(); ++j) present[fpv.getvalue(i,j)] = true;

  std::vector<std::vector<int>> vpf(max_face+1);
  bool tmp;
  for (size_t i=0; i<present.size(); ++i) if (present[i])
  for (size_t j=0; j<fpv.size(); ++j){
    tmp = false;
    for (size_t k=0; k<fpv.numel(); ++k) if (fpv.getvalue(j,k)==i) tmp=true;
    if (tmp) vpf[i].push_back(j);
  }
  // std::vector<ArrayVector<int>> av_vpf;
  // for (auto i: vpf){
  //   if (i.size()) av_vpf.push_back(ArrayVector<int>(1u, i.size(), i.data()));
  //   else av_vpf.push_back(ArrayVector<int>(1u,0u));
  // }
  return vpf;
}

ArrayVector<bool> keep_vertex(LQVec<double>& verts, ArrayVector<int>& fpv){
  std::vector<std::vector<int>> vpf = fpv_to_vpf(verts, fpv);
  ArrayVector<bool> vert_ok(1u, verts.size());
  for (size_t i=0; i<verts.size(); ++i) vert_ok.insert(false, i);
  LQVec<double> face_verts(verts.get_lattice(),0u);
  ArrayVector<bool> isunique;

  std::vector<std::vector<size_t>> equal(verts.size());
  for (size_t i=0; i<verts.size(); ++i) for (size_t j=0; j<verts.size(); ++j)
  if (i!=j && verts.isapprox(i,j)) equal[i].push_back(j);

  for (auto face: vpf){
    if (face.size()==0) continue;
    face_verts.resize(face.size());
    for (size_t i=0; i<face.size(); ++i) face_verts.set(i, verts.extract(face[i]));
    isunique = face_verts.is_unique();
    if (isunique.count_true() > 2){ //more checking to do
      for (size_t i=0; i<face.size(); ++i) if (isunique.getvalue(i)) vert_ok.insert(true, face[i]);
    }
  }
  // // now actually purge the equal vertices
  // for (size_t i=0; i<verts.size(); ++i) if (vert_ok.getvalue(i))
  // for (size_t j: equal[i]) vert_ok.insert(false, j);

  return vert_ok;
}

void BrillouinZone::irreducible_vertex_search(void){
  /* We need to check for three-plane intersections for all combinations of two
     1st Brillouin zone planes and one irreducible reciprocal space normal and
     two irreducible reciprocal space normals and one 1st Brillouin zone plane.
  */
  size_t Nbz = this->faces.size();
  size_t Nir = this->ir_wedge_normals.size();

  if (0==Nir){
    ArrayVector<double> no_verts(3,0), no_norms(3,0), no_point(3,0);
    ArrayVector<int> no_ijk(3,0);
    this->ir_vertices=        no_verts;
    this->ir_face_normals=    no_norms;
    this->ir_face_points=     no_point;
    this->ir_faces_per_vertex=no_ijk;
    return;
  }

  // for which there are M*(N*(N-1))/2 + N*(M*(M-1))/2 total possible combinations
  size_t n21 = ((Nbz*(Nbz-1))>>1)*Nir;
  size_t n12 = ((Nir*(Nir-1))>>1)*Nbz;
  size_t n03 = 0;
  for (int i=2; i<Nir; ++i) n03 += (i*(i-1))>>1;


  LQVec<int> bzfaces = this->get_faces();
  LQVec<double> bznormals = bzfaces / norm(bzfaces);
  LQVec<double> bzpoints = bzfaces / (double)2.0;
  LQVec<double> irnormals = this->get_ir_wedge_normals();
  LQVec<double> vertices30 = this->get_vertices();
  ArrayVector<int> i30 = this->get_faces_per_vertex();

  LQVec<double> vertices21(bzfaces.get_lattice(), n21);
  ArrayVector<int> i21(3, n21);
  LQVec<double> vertices12(bzfaces.get_lattice(), n12);
  ArrayVector<int> i12(3, n12);
  LQVec<double> vertices03(bzfaces.get_lattice(), n03);
  ArrayVector<int> i03(3, n03);

  size_t c21=0, c12=0, c03=0;
  for (size_t i=0  ; i<(Nbz-1); ++i)
  for (size_t j=i+1; j<Nbz    ; ++j)
  for (size_t k=0  ; k<Nir     ; ++k)
  if (intersect_at(bznormals.get(i), bzpoints.get(i),
                   bznormals.get(j), bzpoints.get(j),
                   irnormals.get(k),
                   vertices21, c21)){
    i21.insert(i, c21, 0); i21.insert(j, c21, 1); i21.insert(k, c21++, 2);
  }
  for (size_t i=0; i<Nbz; ++i)
  for (size_t j=0; j<(Nir-1); ++j)
  for (size_t k=j+1; k<Nir; ++k)
  if (intersect_at(bznormals.get(i), bzpoints.get(i),
                   irnormals.get(j),
                   irnormals.get(k),
                   vertices12, c12)){
    i12.insert(i, c12, 0); i12.insert(j, c12, 1); i12.insert(k, c12++, 2);
  }
  for (size_t i=0; i<(Nir-2); ++i)
  for (size_t j=i+1; j<(Nir-1); ++j)
  for (size_t k=j+1; k<Nir; ++k)
  if (intersect_at(irnormals.get(i),
                   irnormals.get(j),
                   irnormals.get(k),
                   vertices03, c03)){
    i03.insert(i, c03, 0); i03.insert(j, c03, 1); i03.insert(k, c03++, 2);
  }
  // make sure we shrink all sets of vertices to just those found!
  vertices21.resize(c21); i21.resize(c21);
  vertices12.resize(c12); i12.resize(c12);
  vertices03.resize(c03); i03.resize(c03);
  // we want to find which of the new intersection points are within the first Brillouin zone
  // but since the first Brillouin zone is already defined, we can use isinside!
  ArrayVector<bool> isin12 = this->isinside(vertices12);
  ArrayVector<bool> isin21 = this->isinside(vertices21);
  ArrayVector<bool> isin03 = this->isinside(vertices03); // should all be true, since the wedge only intersects at 0.
  vertices12 = extract(vertices12, isin12);
  vertices21 = extract(vertices21, isin21);
  vertices03 = extract(vertices03, isin03);
  i12        = extract(i12, isin12);
  i21        = extract(i21, isin21);
  i03        = extract(i03, isin03);
  // Now we have four lists of vertices, plus lists of the normal vectors
  // and on-plane points, which define the three planes that intersect at each
  // vertex.
  // We want to combine these lists:

  std::vector<bool> bz_face_present, ir_face_present;
  size_t max_bz_idx=0, max_ir_idx=0;
  for (size_t i=0; i<i30.size(); ++i){
    for (size_t j=0; j<3u; ++j)
    if (i30.getvalue(i,j)>max_bz_idx) max_bz_idx = i30.getvalue(i,j);
  }
  for (size_t i=0; i<i21.size(); ++i){
    if (i21.getvalue(i,0)>max_bz_idx) max_bz_idx = i21.getvalue(i,0);
    if (i21.getvalue(i,1)>max_bz_idx) max_bz_idx = i21.getvalue(i,1);
    if (i21.getvalue(i,2)>max_ir_idx) max_ir_idx = i21.getvalue(i,2);
  }
  for (size_t i=0; i<i12.size(); ++i){
    if (i12.getvalue(i,0)>max_bz_idx) max_bz_idx = i12.getvalue(i,0);
    if (i12.getvalue(i,1)>max_ir_idx) max_ir_idx = i12.getvalue(i,1);
    if (i12.getvalue(i,2)>max_ir_idx) max_ir_idx = i12.getvalue(i,2);
  }
  for (size_t i=0; i<i03.size(); ++i){
    for (size_t j=0; j<3u; ++j)
    if (i03.getvalue(i,j)>max_ir_idx) max_ir_idx = i03.getvalue(i,j);
  }
  size_t num_bz_idx = max_bz_idx+1, num_ir_idx = max_ir_idx+1; // since we count from 0.
  for (size_t i=0; i<num_bz_idx; ++i) bz_face_present.push_back(false);
  for (size_t i=0; i<num_ir_idx; ++i) ir_face_present.push_back(false);
  for (size_t i=0; i<i30.size(); ++i){
    for (size_t j=0; j<3u; ++j)
    bz_face_present[i30.getvalue(i,j)] = true;
  }
  for (size_t i=0; i<i21.size(); ++i){
    bz_face_present[i21.getvalue(i,0)] = true;
    bz_face_present[i21.getvalue(i,1)] = true;
    ir_face_present[i21.getvalue(i,2)] = true;
  }
  for (size_t i=0; i<i12.size(); ++i){
    bz_face_present[i12.getvalue(i,0)] = true;
    ir_face_present[i12.getvalue(i,1)] = true;
    ir_face_present[i12.getvalue(i,2)] = true;
  }
  for (size_t i=0; i<i03.size(); ++i){
    for (size_t j=0; j<3u; ++j)
    ir_face_present[i03.getvalue(i,j)] = true;
  }
  size_t bz_faces=0, ir_faces=0;
  for (bool tf: bz_face_present) if (tf) ++bz_faces;
  for (bool tf: ir_face_present) if (tf) ++ir_faces;

  size_t total_verts;
  total_verts  = vertices30.size() + vertices12.size();
  total_verts += vertices21.size() + vertices03.size();

  LQVec<double> all_verts(bzfaces.get_lattice(), total_verts);
  LQVec<double> all_norms(bzfaces.get_lattice(), bz_faces+ir_faces);
  LQVec<double> all_point(bzfaces.get_lattice(), bz_faces+ir_faces);
  ArrayVector<int> all_ijk(3u, total_verts);

  std::vector<size_t> bz_face_mapped, ir_face_mapped;
  for (size_t i=0; i<num_bz_idx; ++i) bz_face_mapped.push_back(0);
  for (size_t i=0; i<num_ir_idx; ++i) ir_face_mapped.push_back(0);

  size_t face_idx=0;
  for (size_t i=0; i<i30.size(); ++i) for (size_t j=0; j<3u; ++j){
    if (0==bz_face_mapped[i30.getvalue(i,j)]){
      all_norms.set(face_idx, bznormals.extract(i30.getvalue(i,j)));
      all_point.set(face_idx,  bzpoints.extract(i30.getvalue(i,j)));
      bz_face_mapped[i30.getvalue(i,j)] = ++face_idx; // so that bz_face_mapped is the index+1
    }
  }
  for (size_t i=0; i<i21.size(); ++i){
    if (0==bz_face_mapped[i21.getvalue(i,0)]){
      all_norms.set(face_idx, bznormals.extract(i21.getvalue(i,0)));
      all_point.set(face_idx,  bzpoints.extract(i21.getvalue(i,0)));
      bz_face_mapped[i21.getvalue(i,0)] = ++face_idx;
    }
    if (0==bz_face_mapped[i21.getvalue(i,1)]){
      all_norms.set(face_idx, bznormals.extract(i21.getvalue(i,1)));
      all_point.set(face_idx,  bzpoints.extract(i21.getvalue(i,1)));
      bz_face_mapped[i21.getvalue(i,1)] = ++face_idx;
    }
    if (0==ir_face_mapped[i21.getvalue(i,2)]){
      all_norms.set(face_idx, irnormals.extract(i21.getvalue(i,2)));
      all_point.set(face_idx, 0.0*all_point.extract(i21.getvalue(i,2)));
      ir_face_mapped[i21.getvalue(i,2)] = ++face_idx;
    }
  }
  for (size_t i=0; i<i12.size(); ++i){
    if (0==bz_face_mapped[i12.getvalue(i,0)]){
      all_norms.set(face_idx, bznormals.extract(i12.getvalue(i,0)));
      all_point.set(face_idx,  bzpoints.extract(i12.getvalue(i,0)));
      bz_face_mapped[i12.getvalue(i,0)] = ++face_idx;
    }
    if (0==ir_face_mapped[i12.getvalue(i,1)]){
      all_norms.set(face_idx, irnormals.extract(i12.getvalue(i,1)));
      all_point.set(face_idx, 0.0*all_point.extract(i12.getvalue(i,1)));
      ir_face_mapped[i12.getvalue(i,1)] = ++face_idx;
    }
    if (0==ir_face_mapped[i12.getvalue(i,2)]){
      all_norms.set(face_idx, irnormals.extract(i12.getvalue(i,2)));
      all_point.set(face_idx, 0.0*all_point.extract(i12.getvalue(i,2)));
      ir_face_mapped[i12.getvalue(i,2)] = ++face_idx;
    }
  }
  for (size_t i=0; i<i03.size(); ++i) for (size_t j=0; j<3u; ++j){
    if (0==ir_face_mapped[i03.getvalue(i,j)]){
      all_norms.set(face_idx, irnormals.extract(i03.getvalue(i,j)));
      all_point.set(face_idx, 0.0*all_point.extract(i30.getvalue(i,j)));
      ir_face_mapped[i03.getvalue(i,j)] = ++face_idx;
    }
  }

  size_t vert_idx=0;
  for (size_t i=0; i<i30.size(); ++i){
    all_ijk.insert(bz_face_mapped[i30.getvalue(i,0)]-1, vert_idx, 0);
    all_ijk.insert(bz_face_mapped[i30.getvalue(i,1)]-1, vert_idx, 1);
    all_ijk.insert(bz_face_mapped[i30.getvalue(i,2)]-1, vert_idx, 2);
    all_verts.set(vert_idx++, vertices30.extract(i));
  }
  for (size_t i=0; i<i21.size(); ++i){
    all_ijk.insert(bz_face_mapped[i21.getvalue(i,0)]-1, vert_idx, 0);
    all_ijk.insert(bz_face_mapped[i21.getvalue(i,1)]-1, vert_idx, 1);
    all_ijk.insert(ir_face_mapped[i21.getvalue(i,2)]-1, vert_idx, 2);
    all_verts.set(vert_idx++, vertices21.extract(i));
  }
  for (size_t i=0; i<i12.size(); ++i){
    all_ijk.insert(bz_face_mapped[i12.getvalue(i,0)]-1, vert_idx, 0);
    all_ijk.insert(ir_face_mapped[i12.getvalue(i,1)]-1, vert_idx, 1);
    all_ijk.insert(ir_face_mapped[i12.getvalue(i,2)]-1, vert_idx, 2);
    all_verts.set(vert_idx++, vertices12.extract(i));
  }
  for (size_t i=0; i<i03.size(); ++i){
    all_ijk.insert(ir_face_mapped[i03.getvalue(i,0)]-1, vert_idx, 0);
    all_ijk.insert(ir_face_mapped[i03.getvalue(i,1)]-1, vert_idx, 1);
    all_ijk.insert(ir_face_mapped[i03.getvalue(i,2)]-1, vert_idx, 2);
    all_verts.set(vert_idx++, vertices03.extract(i));
  }

  // Find which vertices are inside the irreducible wedge
  ArrayVector<bool> in_wedge = this->isinside_wedge(all_verts);
  all_verts = extract(all_verts, in_wedge);
  all_ijk   = extract(all_ijk,   in_wedge);
  // and which vertices only contribute to zero-area faces
  ArrayVector<bool> keep = keep_vertex(all_verts, all_ijk);
  all_verts = extract(all_verts, keep);
  all_ijk   = extract(all_ijk,   keep);
  // and remove unused normals and in-plane points:
  std::vector<bool> still_contributes;
  for (size_t i=0; i<all_norms.size(); ++i) still_contributes.push_back(false);
  for (size_t i=0; i<all_ijk.size(); ++i) for (size_t j=0; j<3u; ++j) still_contributes[all_ijk.getvalue(i,j)] = true;
  size_t n_contribute = 0;
  for (bool tf: still_contributes) if (tf) ++n_contribute;
  LQVec<double> in_norms(bzfaces.get_lattice(), n_contribute);
  LQVec<double> in_point(bzfaces.get_lattice(), n_contribute);

  std::vector<size_t> face_mapped;
  for (size_t i=0; i<all_norms.size(); ++i) face_mapped.push_back(0);
  face_idx = 0;
  size_t idx;
  for (size_t i=0; i<all_ijk.size(); ++i) for (size_t j=0; j<3u; ++j){
    idx = all_ijk.getvalue(i,j);
    if (0==face_mapped[idx]){
      in_norms.set(face_idx, all_norms.extract(idx));
      in_point.set(face_idx, all_norms.extract(idx));
      face_mapped[idx] = ++face_idx;
    }
    all_ijk.insert(face_mapped[idx]-1, i,j);
  }

  this->set_ir_vertices(all_verts);
  this->set_ir_face_normals(in_norms);
  this->set_ir_face_points(in_point);
  this->set_ir_faces_per_vertex(all_ijk);
}

void BrillouinZone::vertex_search(const int extent){
  // LQVec<int> tau(this->lattice);
  // int ntau = make_all_indices(&tau,extent);
  LQVec<int> tau(this->lattice, make_relative_neighbour_indices(extent) );
  int ntau = (int)tau.size();
  // the number of unique combinations of 3-taus is the number that we need to check
  size_t ntocheck=0;
  // // there is probably a better way to do this, but brute force never hurt anyone
  // for (int i=0; i<(ntau-2); i++) for (int j=i+1; j<(ntau-1); j++) for (int k=j+1; k<ntau; k++) ntocheck++;
  // an improved method: 0.5*∑ᵢ₌₁ᴺ⁻¹ i²-i  -- we can skip 1 and the last value of our for loop is N-1
  for (int i=2; i<ntau; ++i) ntocheck += (i*(i-1))>>1;
  // Is there a closed-form equivalent to this expression?

  LQVec<double> all_vertices(this->lattice,ntocheck);
  ArrayVector<int> all_ijk(3,ntocheck);

  // LQVec<double> tauhat(this->lattice), halftau(this->lattice);
  // ArrayVector<double> lentau = tau.norm();
  ArrayVector<double> lentau = norm(tau);
  LQVec<double> tauhat = tau/lentau;
  double two = 2;
  LQVec<double> halftau(tau/two);

  ArrayVector<double> tauhat_xyz;
  tauhat_xyz = tauhat.get_xyz();

  int count=0;
  for (int i=0; i<(ntau-2); i++){
    for (int j=i+1; j<(ntau-1); j++){
      for (int k=j+1; k< ntau   ; k++){
        if ( three_plane_intersection(tauhat, halftau, tauhat_xyz, i,j,k, all_vertices, count) ){
          all_ijk.insert(i, count, 0); //insert value i at position (count,0)
          all_ijk.insert(j, count, 1);
          all_ijk.insert(k, count, 2);
          count++;
        }
      }
    }
  }
  // there are count intersections of three planes (strictly count<=ntocheck, but probably count < ntocheck/2)

  // next we need to check for the kernel of intersection points which are closer to the origin than any (non-intersection-defining) planes
  LQVec<double> in_verts(this->lattice,count);
  ArrayVector<int> in_ijk(3,count);
  int in_cnt = 0;
  for (int i=0; i<count; i++){
    // this between_origin_and_plane expects all vectors in an orthonormal frame
    if ( between_origin_and_plane( &halftau, &all_vertices, &all_ijk, i, &in_verts, in_cnt, 1e-10 ) ){
      in_ijk.set(in_cnt++, all_ijk.datapointer(i));
    }
  }
  // if ( in_cnt > 0) in_verts.print(0,in_cnt-1);

  // it's possible that multiple three-plane intersections have given the same
  // intersection point -- e.g., for a cubic system the intersection points of
  // (100),(010),(001); (110),(010),(001); (100),(110),(001); (100),(010),(011)
  // (101),(010),(001); ... are all the same point, (111).
  // The true vertex of the Brillouin Zone is the intersection of the three
  // planes with smallest norm

  // First, find the first vertex which is unique of equivalent vertices
  bool *vertisunique = new bool[in_cnt]();
  for (int i=0; i<in_cnt; i++) vertisunique[i] = true;
  for (int i=0; i<in_cnt-1; i++){
    if (vertisunique[i]){
      for (int j=i+1;j<in_cnt; j++){
        if (vertisunique[j] && in_verts.isapprox(i,j)){
          // printf("vert %d == %d\n",i,j);
           vertisunique[j] = false;
         }
      }
    }
  }
  // count up the unique vertices and keep track of their indices
  int unqcnt=0, *unqidx = new int[in_cnt]();
  for (int i=0; i<in_cnt; i++) if (vertisunique[i]) unqidx[unqcnt++]=i;

  //printf("%d unique vertices\n",unqcnt);

  // create a mapping which holds the indices of all equivalent vertices
  // so unqidxmap[1,:] are all of the indices which equal the first unique vertex
  int *unqidxmap = new int[in_cnt*unqcnt]();
  int *numidxmap = new int[unqcnt]();
  for (int i=0; i<unqcnt; i++){
    numidxmap[i]=1;
    unqidxmap[i*in_cnt] = unqidx[i];
    for (int j=0; j<in_cnt; j++){
      if (unqidx[i]!=j && !vertisunique[j] && in_verts.isapprox(unqidx[i],j))
        unqidxmap[ i*in_cnt + numidxmap[i]++ ] = j;
    }
  }
  delete[] vertisunique;
  // and determine the "length" of the planes which define each vertex, maintaining
  // the equivalence relationship already established, but allocating only the
  // memory actually needed by first finding the maximum number intersections
  // which gave the same vertex
  int maxequiv = 0;
  for (int i=0; i<unqcnt; i++) if ( numidxmap[i]>maxequiv) maxequiv=numidxmap[i];
  double *unqlenmap = new double[maxequiv*unqcnt](); // no need to allocate in_cnt*unqcnt memory when maxequiv is known
  for (int i=0; i<unqcnt; i++){
    for (int j=0; j<numidxmap[i]; j++){
      unqlenmap[i*maxequiv + j] = 0;
      for (int k=0; k<3; k++){
        unqlenmap[i*maxequiv + j] += halftau.norm( in_ijk.getvalue(unqidxmap[i*in_cnt+j],k));
      }
    }
  }
  // use the "length" information to select which equivalent vertex we should keep
  int *minequividx = new int[unqcnt]();
  double *minequivlen = new double[unqcnt]();
  for (int i=0; i<unqcnt; i++){
    minequivlen[i] = std::numeric_limits<double>::max(); // better than 1./0.
    for (int j=0; j<numidxmap[i]; j++){
      if ( unqlenmap[i*maxequiv +j] < minequivlen[i]){
        minequividx[i] = unqidxmap[i*in_cnt+j];
        minequivlen[i] = unqlenmap[i*maxequiv +j];
      }
    }
  }
  delete[] unqidx;
  delete[] numidxmap;
  delete[] unqidxmap;
  delete[] unqlenmap;
  delete[] minequivlen;

  if (unqcnt == 0)   throw std::runtime_error("No unique vertices found?!");

  LQVec<double> unq_vrt(this->lattice, unqcnt);
  ArrayVector<int> unq_ijk(3,unqcnt);
  for (int i=0; i<unqcnt; i++){
      unq_vrt.set( i, in_verts.datapointer(minequividx[i]) );
      unq_ijk.set( i, in_ijk.datapointer(minequividx[i]) );
  }
  delete[] minequividx;

  // store the reciprocal space positions of the vertices of the first Brillouin Zone
  this->set_vertices(unq_vrt); // does this work with LQVec smlst_vrt?

  // determine which of the taus actually contribute to at least one vertex
  int ncontrib=0, *contrib = new int[ntau]();
  for (int i=0; i<ntau; i++)
  for (int j=0; j<unqcnt; j++)
  if ( unq_ijk.getvalue(j,0) == i || unq_ijk.getvalue(j,1) == i || unq_ijk.getvalue(j,2) == i ){
    contrib[ncontrib++]=i;
    break;
  }
  ArrayVector<int> faces(3,ncontrib);
  for (int i=0; i<ncontrib; i++) faces.set(i, tau.datapointer(contrib[i]));
  this->set_faces(faces);

  // Each vertex is the intersection of three faces -- smlst_ijk contains the indexes into tau
  // but since tau contains planes which do not contribute to the first Brillouin Zone
  // we still have work to do. Replace the indices in smlst_ijk with their equivalent
  // indices into this->faces (using contrib as the map)
  ArrayVector<int> fpv(3,unqcnt);
  for (int i=0; i<unqcnt; i++)
  for (int j=0; j<3; j++)
  for (int k=0; k<ncontrib; k++){
    if ( unq_ijk.getvalue(i,j) == contrib[k] ){
      fpv.insert(k,i,j);
      break;
    }
  }
  this->set_faces_per_vertex(fpv);

  delete[] contrib;
}

template<typename T> ArrayVector<bool> BrillouinZone::isinside(const LQVec<T>& p) const {
  bool isouter = this->outerlattice.issame(p.get_lattice());
  bool isinner = this->lattice.issame(p.get_lattice());
  if (!(isouter||isinner))
    throw std::runtime_error("Q points provided to BrillouinZone::isinside must be in the standard or primitive lattice used to define the BrillouinZone object");
  ArrayVector<bool> out(1u, p.size());
  LQVec<int> faces;
  if (isouter)
    faces = this->get_faces();
  else
    faces = this->get_primitive_faces();
  LQVec<double> points = faces/2.0;
  LQVec<double> normals = points/norm(points);
  for (size_t i=0; i<p.size(); ++i)
    out.insert( dot(normals, p.get(i)-points).all_approx("<=",0.), i );
  return out;
}

template<typename T> ArrayVector<bool> BrillouinZone::isinside_wedge(const LQVec<T> &p) const {
  bool isouter = this->outerlattice.issame(p.get_lattice());
  bool isinner = this->lattice.issame(p.get_lattice());
  if (!(isouter||isinner))
    throw std::runtime_error("Q points provided to BrillouinZone::isinside_wedge must be in the standard or primitive lattice used to define the BrillouinZone object");
  ArrayVector<bool> out(1u, p.size());
  LQVec<double> normals;
  if (isouter)
    normals = this->get_ir_wedge_normals();
  else
    normals = this->get_primitive_ir_wedge_normals();
  if (normals.size())
    for (size_t i=0; i<p.size(); ++i)
      out.insert(dot(normals, p.get(i)).all_approx(">=",0.), i);
  else
    for (size_t i=0; i<p.size(); ++i)
      out.insert(true, i); // with no normals *all* points are "inside" the wedge
  return out;
}

bool BrillouinZone::ir_moveinto(const LQVec<double>& Q, LQVec<double>& q, LQVec<int>& tau, std::vector<std::array<int,9>> R, const int time_reversal) const {
  bool isouter = this->outerlattice.issame(Q.get_lattice());
  bool isinner = this->lattice.issame(Q.get_lattice());
  if (!(isouter||isinner))
    throw std::runtime_error("Q points provided to BrillouinZone::ir_moveinto must be in the standard or primitive lattice used to define the BrillouinZone object");
  /*
  The Pointgroup symmetry information comes from, effectively, spglib which
  has all rotation matrices defined in the conventional unit cell -- which is
  our `outerlattice`. Consequently we must work in the outerlattice here.
  */
  // get_pointgroup_symmetry constructs the PointSymmetry object from a vector<array<int,9>>
  // the PointSymmetry constructor automatically sorts the operations in this case.
  PointSymmetry psym = get_pointgroup_symmetry(this->outerlattice.get_hall(), time_reversal);
  int max_order= rotation_order(psym.get(psym.size()-1)); // since they're sorted

  LQVec<int> bz_tau = isouter ? this->get_faces() : this->get_primitive_faces();
  ArrayVector<double> bz_tau_len = norm(bz_tau);
  LQVec<double> bz_points = bz_tau/2.0;
  LQVec<double> bz_normals = bz_points/norm(bz_points);
  LQVec<double> ir_normals = isouter ? this->get_ir_wedge_normals() : this->get_primitive_ir_wedge_normals();

  // by chance some Q points might already be in the IR-Bz:
  ArrayVector<bool> in_bz = this->isinside(Q), in_ir = this->isinside_wedge(Q);
  size_t nQ = Q.size();
  // ensure q, tau, and R can hold one for each Q.
  q.resize(nQ);
  tau.resize(nQ);
  R.resize(nQ);
  // We want to find R(q)+τ = Q with q ∈ IR-Bz and R ∈ Pointgroup operations.
  LQVec<double> qi(Q.get_lattice(), 1u), qj(Q.get_lattice(), 1u);
  LQVec<int> taui;
  std::array<int,9> Ri, Rj, RE={1,0,0, 0,1,0, 0,0,1};
  ArrayVector<double> q_dot_bz_normals;
  ArrayVector<int> Nhkl;

  size_t count, count1, count2, maxat=0, nBZtau=bz_tau.size(), nR=psym.size();
  int maxnm = 0;
  for (size_t i=0; i<nQ; ++i){
    count = 0;
    count1= 0;
    count2= 0;
    qi = Q.get(i);
    taui = 0*tau.get(i);
    Ri = RE;
    while (!(in_bz.getvalue(i) && in_ir.getvalue(i)) && count++ < 50*max_order*nR*nBZtau){
      while (!in_bz.getvalue(i) && count1++ < 50*nBZtau){
        q_dot_bz_normals = dot(qi, bz_normals);
        Nhkl = (q_dot_bz_normals/bz_tau_len).round();
        if (!Nhkl.all_zero()){ // qi is outside of the 1st Bz
          maxnm = 0;
          maxat = 0;
          for (size_t j=0; j<Nhkl.size(); ++j){
            if (Nhkl.getvalue(j)>=maxnm && (maxnm==0 || q_dot_bz_normals.getvalue(j)>q_dot_bz_normals.getvalue(maxat))){
              maxnm = Nhkl.getvalue(j);
              maxat = j;
            }
          }
          qi -= bz_tau[maxat] * static_cast<double>(maxnm);
          taui += bz_tau[maxat];
          in_bz.insert(this->isinside(qi).getvalue(0), i);
        } else {
          in_bz.insert(true, i);
        }
      }
      // moving Q into the first Bz could move q in or out of the IR wedge
      in_ir.insert(this->isinside_wedge(qi).getvalue(0), i);
      // If the rotation matrices, R, form a complete pointgroup
      // then for every Rᵢ and Rⱼ ∈ R, RᵢRⱼ≡Rₖ ∈ R. This means that we do not
      // need to check combinations of rotation matrices (or powers) for
      // complete pointgroups.
      if (!in_ir.getvalue(i)){
        for (size_t j=0; j<nR; ++j){
          // we could check if the rotation axis u∥qi, but that is probably more work than just "rotating" and getting the same point
          multiply_matrix_vector(qj.datapointer(0), psym.get(j), qi.datapointer(0));
          if(this->isinside_wedge(qj).getvalue(0)){
            qi = qj;
            Ri = psym.getarray(j);
            in_ir.insert(true, i);
            break;
          }
        }
      }
      if (!in_ir.getvalue(i)){
        std::string msg = "Either the provided point symmetry operations do not";
        msg += " form a complete group or there is a bug in this code.";
        msg += " Either way, please fix me.";
        throw std::runtime_error(msg);
      }
      // moving qi into the irreducible wedge could have moved it out of the
      // first brillouin zone.
      in_bz.insert(this->isinside(qi).getvalue(0), i);
    }
    // now we have found *a* set of q, tau, R for which R⁻¹q + tau = Q
    q.set(i, qi);
    tau.set(i, taui);
    // We *could* calculate R⁻¹ but might run into rounding/type issues.
    // We know that there must be an element of R for which RᵢRⱼ=E,
    // the identity element. Let's look for Rⱼ instead.
    for (size_t j=0; j<nR; ++j){
      multiply_matrix_matrix(Rj.data(), psym.get(j), Ri.data());
      if (approx_matrix(3, RE.data(), Rj.data())){
        R[i] = psym.getarray(j);
        break;
      }
    }
  }
  if (!in_bz.all_true() || !in_ir.all_true()){
    std::string msg;
    for (size_t i=0; i<nQ; ++i)
      if (!in_bz.getvalue(i) || !in_ir.getvalue(i)){
        msg += "Q=" + Q.to_string(i) + " is outside of the BrillouinZone "
            + " : tau = " + tau.to_string(i) + " , q = " + q.to_string(i)
            + " R = [";
        for (int r: R[i]) msg += " " + std::to_string(r);
        msg += " ]\n" ;
        throw std::runtime_error(msg);
      }
  }
  return (in_bz.all_true() && in_ir.all_true());
}

bool BrillouinZone::moveinto(const LQVec<double>& Q, LQVec<double>& q, LQVec<int>& tau){
  bool already_same = this->lattice.issame(Q.get_lattice());
  LQVec<double> Qprim(this->lattice), qprim(this->lattice);
  LQVec<int> tauprim(this->lattice);
  PrimitiveTransform PT(this->outerlattice.get_hall());
  bool transform_needed = ( PT.does_anything() && this->outerlattice.issame(Q.get_lattice()) );
  if (!(already_same || transform_needed))
    throw std::runtime_error("Q points provided to BrillouinZone::isinside must be in the standard or primitive lattice used to define the BrillouinZone object");

  if (transform_needed)  Qprim = transform_to_primitive(this->outerlattice,Q);
  const LQVec<double> & Qsl = transform_needed ? Qprim : Q;
  LQVec<double> & qsl = transform_needed ? qprim : q;
  LQVec<int> & tausl = transform_needed? tauprim : tau;

  // Determine which points in Q are already inside the first BZ
  ArrayVector<bool> allinside = this->isinside(Qsl);
  // ensure that qsl and tausl can hold each qi and taui
  qsl.resize(Qsl.size());
  tausl.resize(Qsl.size());

  LQVec<int> facehkl(this->lattice,this->faces);
  ArrayVector<double> facelen = norm(facehkl);
  LQVec<double> facenrm = facehkl/facelen;
  LQVec<double> qi;
  LQVec<int> taui;
  ArrayVector<double> q_dot_facenrm;
  ArrayVector<int> Nhkl;
  size_t maxat = 0;
  int maxnm = 0;
  size_t count =0;
  for (size_t i=0; i<Qsl.size(); i++){
    count = 0;
    qi = Qsl.get(i);
    taui = 0*tausl.get(i);
    while (!allinside.getvalue(i) && count++ < 50*facelen.size()){
      // std::cout << "Moving q = " << qi.to_string() << std::endl;
      q_dot_facenrm = dot( qi , facenrm );
      Nhkl = (q_dot_facenrm/facelen).round();
      // std::cout << "Nhkl = " << Nhkl.to_string() << std::endl;
      if ( Nhkl.all_zero() ) {allinside.insert(true,i); break;} // qi is *on* the Brilluoin Zone surface (or inside) so break.
      maxnm = 0;
      maxat = 0;
      for (size_t j=0; j<Nhkl.size(); ++j){
        if (Nhkl.getvalue(j)>=maxnm && (maxnm==0 || q_dot_facenrm.getvalue(j)>q_dot_facenrm.getvalue(maxat)) ){
          maxnm = Nhkl.getvalue(j);
          maxat = j;
        }
      }
      // std::cout << "Of which, the maximum is vector " << std::to_string(maxat);
      // std::cout << " with value " << facehkl.to_string(maxat) << " " << std::to_string(maxnm);
      // std::cout << " of which will be removed." << std::endl;

      qi -= facehkl[maxat] * (double)(maxnm); // ensure we subtract LQVec<double>
      taui += facehkl[maxat] * maxnm; // but add LQVec<int>

      allinside.insert(this->isinside(qi).getvalue(0), i);
    }
    qsl.set(i, &qi);
    tausl.set(i, &taui);
  }
  if (!allinside.all_true()){
    std::string msg;
    for (size_t i=0; i<Qsl.size(); ++i)
      if (!allinside.getvalue(i))
        msg += "Q=" + Qsl.to_string(i) + " is outside of the BrillouinZone "
            + " : tau = " + tausl.to_string(i) + " , q = " + qsl.to_string(i) + "\n";
    throw std::runtime_error(msg);
  }
  if (transform_needed){ // then we need to transform back q and tau
    q   = transform_from_primitive(this->outerlattice,qsl);
    tau = transform_from_primitive(this->outerlattice,tausl);
  }
  return allinside.all_true(); // return false if any points are still outside of the first Brilluoin Zone
}

bool three_plane_intersection(const LQVec<double>& n,                // plane normals
                              const LQVec<double>& p,                // a point on each plane
                              const ArrayVector<double>& xyz,        // the plane normals in a inverse Angstrom orthonormal coordinate system
                              const int i, const int j, const int k, // indices into n, p, xyz
                              LQVec<double>& iat,                    // output storage array of intersection points
                              const int idx)                         // the index where the intersection is inserted if found
                              {
  // we need to check whether the matrix formed by the orthonormal-frame components of the three planes is nearly-singular
  double *M = new double[9];
  xyz.get(i, M);
  xyz.get(j, M+3);
  xyz.get(k, M+6);
  double detM;
  detM = matrix_determinant(M);
  delete[] M;
  if ( std::abs(detM) > 1e-10 ){ // this 1e-10 provides a cutoff for far-from-origin intersections
    LQVec<double> ni,nj,nk, pi,pj,pk, cij,cjk,cki, tmp;
    ni=n.get(i);  nj=n.get(j);  nk=n.get(k);
    pi=p.get(i);  pj=p.get(j);  pk=p.get(k);
    // cij=ni.cross(&nj);  cjk=nj.cross(&nk);  cki=nk.cross(&ni);
    cij=cross(ni,nj);  cjk=cross(nj,nk);  cki=cross(nk,ni);

    // tmp = cjk*(pi.dot(&ni)) + cki*(pj.dot(&nj)) + cij*(pk.dot(&nk));
    tmp = cjk*dot(pi,ni) + cki*dot(pj,nj) + cij*dot(pk,nk);
    tmp /= detM;

    iat.set(idx, tmp.datapointer() );
    return true;
  }
  return false;
}

bool intersect_at(const LQVec<double>& ni, const LQVec<double>& pi,
                  const LQVec<double>& nj, const LQVec<double>& pj,
                  const LQVec<double>& nk, const LQVec<double>& pk,
                  LQVec<double>& intersect, const int idx
                 ){
  double detM, M[9];
  ni.get_xyz().get(0, M  );
  nj.get_xyz().get(0, M+3);
  nk.get_xyz().get(0, M+6);
  detM = matrix_determinant(M);
  if (std::abs(detM) > 1e-10){
    LQVec<double> tmp = cross(nj,nk)*dot(pi,ni)
                      + cross(nk,ni)*dot(pj,nj)
                      + cross(ni,nj)*dot(pk,nk);
    tmp /= detM;
    intersect.set(idx, tmp.datapointer());
    return true;
  }
  return false;
}
bool intersect_at(const LQVec<double>& ni, const LQVec<double>& pi,
                  const LQVec<double>& nj, const LQVec<double>& pj,
                  const LQVec<double>& nk,
                  LQVec<double>& intersect, const int idx
                 ){
  double detM, M[9];
  ni.get_xyz().get(0, M  );
  nj.get_xyz().get(0, M+3);
  nk.get_xyz().get(0, M+6);
  detM = matrix_determinant(M);
  if (std::abs(detM) > 1e-10){
    LQVec<double> tmp = cross(nj,nk)*dot(pi,ni) + cross(nk,ni)*dot(pj,nj);
    tmp /= detM;
    intersect.set(idx, tmp.datapointer());
    return true;
  }
  return false;
}
bool intersect_at(const LQVec<double>& ni, const LQVec<double>& pi,
                  const LQVec<double>& nj,
                  const LQVec<double>& nk,
                  LQVec<double>& intersect, const int idx
                 ){
  double detM, M[9];
  ni.get_xyz().get(0, M  );
  nj.get_xyz().get(0, M+3);
  nk.get_xyz().get(0, M+6);
  detM = matrix_determinant(M);
  if (std::abs(detM) > 1e-10){
    LQVec<double> tmp = cross(nj,nk)*dot(pi,ni);
    tmp /= detM;
    intersect.set(idx, tmp.datapointer());
    return true;
  }
  return false;
}
bool intersect_at(const LQVec<double>& ni,
                  const LQVec<double>& nj,
                  const LQVec<double>& nk,
                  LQVec<double>& intersect, const int idx
                 ){
  intersect.set(idx, 0*intersect.get(idx));
  return true;
}

bool between_origin_and_plane(const LQVec<double> *p,
                              const LQVec<double> *v,
                              const ArrayVector<int> *ijk,
                              const int idx,
                              LQVec<double> *inv,
                              const int store_at,
                              const double tol){
  // p and v should be the points defining each plane and the vertices of the intersections of three planes
  ArrayVector<double> v_p = dot(v->get(idx),*p), p_p = dot(*p,*p);
  // we want to skip over the planes which gave us this intersection point
  size_t i, skip1, skip2, skip3;
  skip1 = (size_t) ijk->getvalue(idx,0);
  skip2 = (size_t) ijk->getvalue(idx,1);
  skip3 = (size_t) ijk->getvalue(idx,2);

  double vpi, ppi, eps = std::numeric_limits<double>::epsilon();
  for (i=0; i < p->size(); i++){
    if ( !(i==skip1||i==skip2||i==skip3) ){
      vpi = v_p.getvalue(i);
      ppi = p_p.getvalue(i);
      if ((vpi-ppi) > (vpi+ppi)*eps && (vpi-ppi) > eps && (vpi-ppi) > tol)
        return false;
    }
  }
  // none of p are closer to the origin than v(i)
  inv->set(store_at, v->datapointer(idx));
  return true;
}
