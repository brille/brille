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

#include "bz.hpp"
#include <omp.h>

template<typename... A>
void BrillouinZone::set_polyhedron(const LQVec<double>& v, const LQVec<double>& p, A... args){
  bool both_same = v.get_lattice().issame(p.get_lattice());
  if (!both_same)
    throw std::runtime_error("The vertices, and points of a polyhedron must all be in the same cooridinate system");
  bool is_outer = this->outerlattice.issame(v.get_lattice());
  bool is_inner = this->lattice.issame(v.get_lattice());
  LQVec<double> vp(this->outerlattice), pp(this->outerlattice);
  PrimitiveTransform PT(this->outerlattice.get_spacegroup_object().get_bravais_type());
  is_inner &= PT.does_anything();
  if (!(is_outer || is_inner))
    throw std::runtime_error("The polyhedron must be described in the conventional or primitive lattice used to define the BrillouinZone object");
  if (is_inner){
    vp = transform_from_primitive(this->outerlattice, v);
    pp = transform_from_primitive(this->outerlattice, p);
  }
  const LQVec<double> & vref = is_inner ? vp : v;
  const LQVec<double> & pref = is_inner ? pp : p;
  this->polyhedron = Polyhedron(vref.get_xyz(), pref.get_xyz(), args...);
}
template<typename... A>
void BrillouinZone::set_ir_polyhedron(const LQVec<double>& v, const LQVec<double>& p, const LQVec<double>& n, A... args){
  bool all_same = v.get_lattice().issame(p.get_lattice()) && p.get_lattice().issame(n.get_lattice());
  if (!all_same)
    throw std::runtime_error("The vertices, points, and normals of a polyhedron must all be in the same cooridinate system");
  bool is_outer = this->outerlattice.issame(v.get_lattice());
  bool is_inner = this->lattice.issame(v.get_lattice());
  LQVec<double> vp(this->outerlattice), pp(this->outerlattice), np(this->outerlattice);
  PrimitiveTransform PT(this->outerlattice.get_spacegroup_object().get_bravais_type());
  is_inner &= PT.does_anything();
  if (!(is_outer || is_inner))
    throw std::runtime_error("The polyhedron must be described in the conventional or primitive lattice used to define the BrillouinZone object");
  if (is_inner){
    vp = transform_from_primitive(this->outerlattice, v);
    pp = transform_from_primitive(this->outerlattice, p);
    np = transform_from_primitive(this->outerlattice, n);
  }
  const LQVec<double> & vref = is_inner ? vp : v;
  const LQVec<double> & pref = is_inner ? pp : p;
  const LQVec<double> & nref = is_inner ? np : n;
  this->ir_polyhedron = Polyhedron(vref.get_xyz(), pref.get_xyz(), nref.get_xyz(), args...);
}
LQVec<double> BrillouinZone::get_ir_polyhedron_wedge_normals(void) const {
  LQVec<double> ir_n = this->get_ir_normals();
  LQVec<double> ir_p = this->get_ir_points();
  LQVec<double> bz_n = this->get_normals();
  LQVec<double> bz_p = this->get_points();
  ArrayVector<bool> not_bz(1u, 0u);
  for (size_t i=0; i<bz_n.size(); ++i){
    // check if the irBZ face point is on a first BZ zone face too
    not_bz = dot(bz_n.extract(i), ir_p - bz_p.extract(i)).is_approx(Comp::neq, 0.);
    ir_n = ir_n.extract(not_bz);
    ir_p = ir_p.extract(not_bz);
  }
  // It is possible that, lacking inversion symmetry, we have found an
  // irreducible polyhedron which is comprised of two convex polyhedra that
  // are mutualy inverse. In such a case for every normal in ir_n there is also
  // its opposite also in ir_n, and so ir_n defines a wedge and an anti-wedge in
  // which no finite point can be inside
  if (ir_n.size()%2 == 0 /* [lacks inversion] or this->no_ir_mirroring? */){
    std::vector<bool> no_inverse(ir_n.size(), true);
    ArrayVector<bool> keep(1u, ir_n.size(), true);
    for (size_t i=0; i<ir_n.size()-1; ++i) if (no_inverse[i])
    for (size_t j=i+1; j<ir_n.size(); ++j) if (no_inverse[j])
    if ((ir_n.extract(i)+ir_n.extract(j)).all_approx(Comp::eq, 0.)) {
      no_inverse[i] = no_inverse[j] = false;
      keep.insert(false, j);
      break;
    }
    if (std::count_if(no_inverse.begin(), no_inverse.end(), [](bool t){return t;})==0)
      ir_n = ir_n.extract(keep);
  }
  // the remaining irBZ faces are the irreducible reciprocal space wedge
  // which we store with the opposite sign in this->ir_wedge_normals
  return -ir_n;
}
bool BrillouinZone::set_ir_vertices(const LQVec<double>& v){
  bool is_outer = this->outerlattice.issame(v.get_lattice());
  bool is_inner = this->lattice.issame(v.get_lattice());
  LQVec<double> vp(this->outerlattice);
  PrimitiveTransform PT(this->outerlattice.get_spacegroup_object().get_bravais_type());
  is_inner &= PT.does_anything();
  if (!(is_outer || is_inner))
    throw std::runtime_error("The polyhedron must be described in the conventional or primitive lattice used to define the BrillouinZone object");
  if (is_inner)
    vp = transform_from_primitive(this->outerlattice, v);
  const LQVec<double> & vref = is_inner ? vp : v;
  debug_update("Generate a convex polyhedron from\n", vref.to_string());
  this->ir_polyhedron = Polyhedron(vref.get_xyz());
  debug_update("Generated a ", this->ir_polyhedron.string_repr()," from the specified vertices");
  // check that the polyhedron we've specified has the desired properties
  if (this->check_ir_polyhedron()){
    // set the wedge normals as well
    this->set_ir_wedge_normals(this->get_ir_polyhedron_wedge_normals());
    return true;
  }
  return false;
}
Polyhedron BrillouinZone::get_polyhedron(void) const {return this->polyhedron;}
Polyhedron BrillouinZone::get_ir_polyhedron(const bool true_ir) const {
  // If the ir polyhedron fills the first BZ using the symmetry operations
  // of the pointgroup (plus time inversion, if included), then no ir_mirroring
  // is to be performed. In this case or if the 'true_ir' was requested, return
  // the computed ir_polyhedron
  if (this->no_ir_mirroring || !true_ir) return this->ir_polyhedron;
  // Otherwise, the ir_polyedron is only half of the true (non-convex)
  // irreducible polyhedron, so add the lack of inversion symmetry explicitly.
  return this->ir_polyhedron + this->ir_polyhedron.mirror();
}
bool BrillouinZone::check_ir_polyhedron(void){
  this->check_if_mirroring_needed(); // move this to end of wedge_brute_force?
  PointSymmetry fullps = this->outerlattice.get_pointgroup_symmetry(this->time_reversal?1:0);
  double volume_goal = this->polyhedron.get_volume() / static_cast<double>(fullps.size());
  Polyhedron irbz = this->get_ir_polyhedron(), rotated;
  if (!approx_scalar(irbz.get_volume(), volume_goal)){
    debug_update("The current 'irreducible' polyhedron has the wrong volume");
    debug_update("Since ",irbz.get_volume()," != ",volume_goal);
    return false;
  }
  /* TODO FIXME -- make a LatticePolyhedron class to handle this? */
  /* Our Polyhedron class holds the vertices, and plane information of a
     convex polyhedron expressed in an orthonormal frame tied to the
     conventional reciprocal space lattice. The conversion from lattice
     vectors to absolute vectors is done with a transformation matrix, B.
     The rotation and rotoinversions, R, stored in the PointSymmetry object
     are expressed in units of the conventional real space lattice.

     Reciprocal space vectors rotate via the transpose of R, Rᵀ; and since the
     polyhedron contains vertices, points, and normals of the form Bx and we
     need to apply Rᵀ to x directly, we want to "rotate" each part of the
     polyhedron by B Rᵀ B⁻¹, e.g., B Rᵀ B⁻¹ B x = B Rᵀ x.
  */
  // double B[9], invB[9], RtinvB[8];
  std::array<double,9> B, invB, RtinvB, BRtinvB;
  this->outerlattice.get_B_matrix(B.data());
  verbose_update("B\n",B);
  matrix_inverse(invB.data(), B.data());
  verbose_update("inverse(B)\n",invB);
  std::array<int,9> Rt;
  // get the operations of the pointgroup which are not 1 or -1
  // keeping -1 would probably be ok, but hopefully it won't hurt to remove it now
  PointSymmetry ps = fullps.higher(1);
  for (size_t i=0; i<ps.size(); ++i){
    Rt = transpose(ps.get(i)); // Rᵀ
    debug_update("\n\n\n\n\ntranspose(R)\n",Rt);
    multiply_matrix_matrix(RtinvB.data(), Rt.data(), invB.data()); // Rᵀ*B⁻¹
    debug_update("transpose(R) inverse(B)\n",RtinvB);
    multiply_matrix_matrix(BRtinvB.data(), B.data(), RtinvB.data()); // B*Rᵀ*B⁻¹
    debug_update("Rotate the polyhedron by\n", BRtinvB);
    rotated = irbz.rotate(BRtinvB);
    debug_update("Checking for intersection ",i);
    if (irbz.intersects(rotated)){
      debug_update("The current 'irreducible' polyhedron intersects with itself and therefore is not correct.");
      return false;
    }
  }
  // volume is right and no intersections
  return true;
}

// first Brillouin zone
LQVec<double> BrillouinZone::get_vertices(void) const {
  ArrayVector<double> v = this->polyhedron.get_vertices(); // in the outerlattice xyz coordinate system
  LQVec<double> lv(this->outerlattice, v.size());
  double fromxyz[9];
  this->outerlattice.get_inverse_xyz_transform(fromxyz);
  for (size_t i=0; i<v.size(); i++)
    multiply_matrix_vector<double,double,double,3>(lv.data(i), fromxyz, v.data(i));
  return lv;
}
LQVec<double> BrillouinZone::get_primitive_vertices(void) const {
  LQVec<double> v = this->get_vertices();
  if (this->isprimitive()) v = transform_to_primitive(this->outerlattice, v);
  return v;
}
LQVec<double> BrillouinZone::get_points(void) const {
  ArrayVector<double> p = this->polyhedron.get_points();
  LQVec<double> lp(this->outerlattice, p.size());
  double fromxyz[9];
  this->outerlattice.get_inverse_xyz_transform(fromxyz);
  for (size_t i=0; i<p.size(); i++)
    multiply_matrix_vector<double,double,double,3>(lp.data(i), fromxyz, p.data(i));
  return lp;
}
LQVec<double> BrillouinZone::get_primitive_points(void) const {
  LQVec<double> p = this->get_points();
  if (this->isprimitive()) p = transform_to_primitive(this->outerlattice, p);
  return p;
}
LQVec<double> BrillouinZone::get_normals(void) const {
  ArrayVector<double> n = this->polyhedron.get_normals();
  LQVec<double> ln(this->outerlattice, n.size());
  double fromxyz[9];
  this->outerlattice.get_inverse_xyz_transform(fromxyz);
  for (size_t i=0; i<n.size(); i++)
    multiply_matrix_vector<double,double,double,3>(ln.data(i), fromxyz, n.data(i));
  return ln;
}
LQVec<double> BrillouinZone::get_primitive_normals(void) const {
  LQVec<double> n = this->get_normals();
  if (this->isprimitive()) n = transform_to_primitive(this->outerlattice, n);
  return n;
}
LQVec<double> BrillouinZone::get_half_edges(void) const{
  ArrayVector<double> h = this->polyhedron.get_half_edges();
  LQVec<double> lh(this->outerlattice, h.size());
  double fromxyz[9];
  this->outerlattice.get_inverse_xyz_transform(fromxyz);
  for (size_t i=0; i<h.size(); ++i)
  multiply_matrix_vector(lh.data(i), fromxyz, h.data(i));
  return lh;
}
std::vector<std::vector<int>> BrillouinZone::get_faces_per_vertex(void) const {
  return this->polyhedron.get_faces_per_vertex();
}
std::vector<std::vector<int>> BrillouinZone::get_vertices_per_face(void) const {
  return this->polyhedron.get_vertices_per_face();
}
// irreducible first Brillouin zone
LQVec<double> BrillouinZone::get_ir_vertices(void) const {
  Polyhedron irp = this->get_ir_polyhedron();
  ArrayVector<double> v = irp.get_vertices();
  LQVec<double> lv(this->outerlattice, v.size());
  double fromxyz[9];
  this->outerlattice.get_inverse_xyz_transform(fromxyz);
  for (size_t i=0; i<v.size(); i++)
    multiply_matrix_vector<double,double,double,3>(lv.data(i), fromxyz, v.data(i));
  return lv;
}
LQVec<double> BrillouinZone::get_ir_primitive_vertices(void) const {
  LQVec<double> v = this->get_ir_vertices();
  if (this->isprimitive()) v = transform_to_primitive(this->outerlattice, v);
  return v;
}
LQVec<double> BrillouinZone::get_ir_points(void) const {
  Polyhedron irp = this->get_ir_polyhedron();
  ArrayVector<double> p = irp.get_points();
  LQVec<double> lp(this->outerlattice, p.size());
  double fromxyz[9];
  this->outerlattice.get_inverse_xyz_transform(fromxyz);
  for (size_t i=0; i<p.size(); i++)
    multiply_matrix_vector<double,double,double,3>(lp.data(i), fromxyz, p.data(i));
  return lp;
}
LQVec<double> BrillouinZone::get_ir_primitive_points(void) const {
  LQVec<double> p = this->get_ir_points();
  if (this->isprimitive()) p = transform_to_primitive(this->outerlattice, p);
  return p;
}
LQVec<double> BrillouinZone::get_ir_normals(void) const {
  Polyhedron irp = this->get_ir_polyhedron();
  ArrayVector<double> n = irp.get_normals();
  LQVec<double> ln(this->outerlattice, n.size());
  double fromxyz[9];
  this->outerlattice.get_inverse_xyz_transform(fromxyz);
  for (size_t i=0; i<n.size(); i++)
    multiply_matrix_vector<double,double,double,3>(ln.data(i), fromxyz, n.data(i));
  return ln;
}
LQVec<double> BrillouinZone::get_ir_primitive_normals(void) const {
  LQVec<double> n = this->get_ir_normals();
  if (this->isprimitive()) n = transform_to_primitive(this->outerlattice, n);
  return n;
}
std::vector<std::vector<int>> BrillouinZone::get_ir_faces_per_vertex(void) const {
  return this->get_ir_polyhedron().get_faces_per_vertex();
}
std::vector<std::vector<int>> BrillouinZone::get_ir_vertices_per_face(void) const {
  return this->get_ir_polyhedron().get_vertices_per_face();
}

// irreducible reciprocal space
LQVec<double> BrillouinZone::get_ir_wedge_normals(void) const {
  LQVec<double> out(this->outerlattice, 0u);
  if (this->ir_wedge_normals.size())
    out = LQVec<double>(this->outerlattice, this->ir_wedge_normals);
  return out;
}
//
LQVec<double> BrillouinZone::get_primitive_ir_wedge_normals(void) const {
  LQVec<double> lqwn(this->outerlattice, 0u);
  if (this->ir_wedge_normals.size()){
    lqwn = LQVec<double>(this->outerlattice, this->ir_wedge_normals);
    if (this->isprimitive())
      lqwn = transform_to_primitive(this->outerlattice, lqwn);
  }
  return lqwn;
}
void BrillouinZone::set_ir_wedge_normals(const LQVec<double>& x){
  bool already_same = this->outerlattice.issame(x.get_lattice());
  LQVec<double> xp(this->outerlattice);
  PrimitiveTransform PT(this->outerlattice.get_spacegroup_object().get_bravais_type());
  bool transform_needed = ( PT.does_anything() && this->lattice.issame(x.get_lattice()) );
  if (!(already_same || transform_needed))
    throw std::runtime_error("ir_wedge_normals must be in the standard or primitive lattice used to define the BrillouinZone object");
  if (transform_needed)  xp = transform_from_primitive(this->outerlattice,x);
  const LQVec<double> & xref = transform_needed ? xp : x;
  this->ir_wedge_normals = xref.get_hkl();
}

void BrillouinZone::print() const {
  std::string msg = "BrillouinZone with ";
  msg += std::to_string(this->vertices_count()) + " vertices and ";
  msg += std::to_string(this->faces_count()) + " faces";
  std::cout << msg << std::endl;
}

bool BrillouinZone::wedge_normal_check(const LQVec<double>& n, LQVec<double>& normals, size_t& num){
  std::string msg = "Considering " + n.to_string(0,"... ");
  if (norm(n).all_approx(0.0)){
    debug_update(msg, "rejected; zero-length");
    return false;
  }
  if (num==0){
    debug_update(msg, "accepted; first normal");
    normals.set(num, n.get(0));
    num=num+1;
    return true;
  }
  if (norm(cross(normals.first(num), n)).any_approx(Comp::eq,0.)){
    debug_update(msg, "rejected; already present");
    return false;
  }
  normals.set(num,  n.get(0));
  if (this->ir_wedge_is_ok(normals.first(num+1))){
    debug_update(msg, "accepted");
    num=num+1;
    return true;
  }
  normals.set(num, -n.get(0));
  if (this->ir_wedge_is_ok(normals.first(num+1))){
    debug_update(msg, "accepted (*-1)");
    num=num+1;
    return true;
  }
  debug_update(msg, "rejected; addition causes null wedge");
  return false;
}
bool BrillouinZone::wedge_normal_check(const LQVec<double>& n0, const LQVec<double>& n1, LQVec<double>& normals, size_t& num){
  std::string msg = "Considering " + n0.to_string(0," and ") + n1.to_string(0,"... ");
  bool p0=false, p1=false;
  if (num>0){
    p0 = norm(cross(normals.first(num), n0)).any_approx(Comp::eq,0.);
    p1 = norm(cross(normals.first(num), n1)).any_approx(Comp::eq,0.);
  }
  if (p0 && p1){
    debug_update(msg, "rejected; already present");
    return false;
  }
  if (num>0 && dot(normals.first(num)/norm(n0),n0/norm(n0)).any_approx(Comp::eq,1.)){
    normals.set(num,  n1.get(0));
    if (this->ir_wedge_is_ok(normals.first(num+1))){
      debug_update(msg, "n1 accepted (n0 present)");
      num=num+1;
      return true;
    }
    debug_update(msg, "n1 rejected (n0 present); addition causes null wedge");
    return false;
  }
  if (num>0 && dot(normals.first(num)/norm(n1),n1/norm(n1)).any_approx(Comp::eq,1.)){
    normals.set(num,  n0.get(0));
    if (this->ir_wedge_is_ok(normals.first(num+1))){
      debug_update(msg, "n0 accepted (n1 present)");
      num=num+1;
      return true;
    }
    debug_update(msg, "n0 rejected (n1 present); addition causes null wedge");
    return false;
  }
  if (num>0 && (p0 || p1)){
    debug_update_if(p0, msg, "-n0 present; addition causes null wedge)");
    debug_update_if(p1, msg, "-n1 present; addition causes null wedge)");
    return false;
  }
  normals.set(num,   n0.get(0));
  normals.set(num+1, n1.get(0));
  if (this->ir_wedge_is_ok(normals.first(num+2))){
    debug_update(msg, "n0 & n1 accepted");
    num=num+2;
    return true;
  }
  debug_update(msg, "n0 & n1 rejected; adding both would cause null wedge");
  return false;
}
bool BrillouinZone::ir_wedge_is_ok(const LQVec<double>& normals){
  this->set_ir_wedge_normals(normals); // assigns this->ir_wedge_normals
  this->irreducible_vertex_search(); // assigns this->ir_polyhedron
  return !approx_scalar(this->ir_polyhedron.get_volume(), 0.0);
}

void BrillouinZone::shrink_and_prune_outside(const size_t cnt, LQVec<double>& vrt, ArrayVector<int>& ijk) const {
  verbose_update("shrinking to ",cnt);
  if(vrt.size() && ijk.size()){
    vrt.resize(cnt);
    ijk.resize(cnt);
    if (cnt){ // isinside has a problem with vrt.size()==0
      ArrayVector<bool> isin = this->isinside(vrt);
      verbose_update("and retaining ", isin.count_true(), " inside vertices");
      vrt = extract(vrt, isin);
      ijk = extract(ijk, isin);
    }
  }
}

void BrillouinZone::irreducible_vertex_search(){
  /* We need to check for three-plane intersections for all combinations of two
     1st Brillouin zone planes and one irreducible reciprocal space normal and
     two irreducible reciprocal space normals and one 1st Brillouin zone plane.
  */
  size_t Nbz = this->get_normals().size();
  size_t Nir = this->ir_wedge_normals.size();

  if (0==Nir){
    this->ir_polyhedron = this->polyhedron;
    return;
  }

  // for which there are M*(N*(N-1))/2 + N*(M*(M-1))/2 total possible combinations
  size_t n21 = ((Nbz*(Nbz-1))>>1)*Nir;
  size_t n12 = ((Nir*(Nir-1))>>1)*Nbz;
  size_t n03 = 0;
  for (size_t i=2; i<Nir; ++i) n03 += (i*(i-1))>>1;
  verbose_update("Checking {",n21,", ",n12,", ",n03,"} {2:1, 1:2, 0:3} zone:wedge 3-plane intersection points");

  LQVec<double> bznormals = this->get_normals();
  LQVec<double> bzpoints = this->get_points();

  // We will create a polyhedron using (some of) these normals. It is imperitive
  // that the polyhedron normals point *outwards* from the centre of the polyhedron
  // while we've thus far defined the wedge normals to point *into* the irreducible
  // reciprocal space wedge.
  LQVec<double> irnormals = -1.0*this->get_ir_wedge_normals();
  LQVec<double> vertices30 = this->get_vertices();
  std::vector<std::vector<int>> i30 = this->get_faces_per_vertex();

  LQVec<double> vertices21(bznormals.get_lattice(), n21);
  ArrayVector<int> i21(3, n21);
  LQVec<double> vertices12(bznormals.get_lattice(), n12);
  ArrayVector<int> i12(3, n12);
  LQVec<double> vertices03(bznormals.get_lattice(), n03);
  ArrayVector<int> i03(3, n03);

  int c21=0, c12=0, c03=0;
  if (n21){ // protect against Nbz=0, since size_t(0)-1 = 4294967294 or 18446744073709551615 if its 32- or 64-bit
    for (size_t i=0  ; i<(Nbz-1); ++i)
    for (size_t j=i+1; j< Nbz   ; ++j)
    for (size_t k=0  ; k< Nir   ; ++k)
    if (intersect_at(bznormals.get(i), bzpoints.get(i),
                     bznormals.get(j), bzpoints.get(j),
                     irnormals.get(k),
                     vertices21, c21)){
      i21.insert(i, c21, 0); i21.insert(j, c21, 1); i21.insert(k, c21++, 2);
    }
  }
  if (n12){ // protect against Nir=0, since size_t(0)-1 = 4294967294 or 18446744073709551615 if its 32- or 64-bit
    for (size_t i=0  ; i< Nbz   ; ++i)
    for (size_t j=0  ; j<(Nir-1); ++j)
    for (size_t k=j+1; k< Nir   ; ++k)
    if (intersect_at(bznormals.get(i), bzpoints.get(i),
                     irnormals.get(j),
                     irnormals.get(k),
                     vertices12, c12)){
      i12.insert(i, c12, 0); i12.insert(j, c12, 1); i12.insert(k, c12++, 2);
    }
  }
  if (n03){
    for (size_t i=0  ; i<(Nir-2); ++i)
    for (size_t j=i+1; j<(Nir-1); ++j)
    for (size_t k=j+1; k< Nir   ; ++k)
    if (intersect_at(irnormals.get(i),
                     irnormals.get(j),
                     irnormals.get(k),
                     vertices03, c03)){
      i03.insert(i, c03, 0); i03.insert(j, c03, 1); i03.insert(k, c03++, 2);
    }
  }
  verbose_update("Intersections found");
  // make sure we shrink all sets of vertices to just those found!
  // plus remove any intersection points outside of the first Brillouin zone
  this->shrink_and_prune_outside(static_cast<size_t>(c21), vertices21, i21);
  this->shrink_and_prune_outside(static_cast<size_t>(c12), vertices12, i12);
  this->shrink_and_prune_outside(static_cast<size_t>(c03), vertices03, i03);
  verbose_update("Intersections pruned");
  // Now we have four lists of vertices, plus lists of the normal vectors
  // and on-plane points, which define the three planes that intersect at each
  // vertex.
  // We want to combine these lists:

  int max_bz_idx=0, max_ir_idx=0;
  for (auto i: i30) for (int j: i) if (j > max_bz_idx) max_bz_idx = j;
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
  max_bz_idx++; max_ir_idx++; // since we count from 0, and need one more element than the highest index.
  std::vector<bool> bz_face_present(max_bz_idx, false), ir_face_present(max_ir_idx, false);
  for (auto i: i30) for (int j: i) bz_face_present[j] = true;
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

  LQVec<double> all_verts(bznormals.get_lattice(), total_verts);
  LQVec<double> all_norms(bznormals.get_lattice(), bz_faces+ir_faces);
  LQVec<double> all_point(bznormals.get_lattice(), bz_faces+ir_faces);
  ArrayVector<int> all_ijk(3u, total_verts);

  std::vector<size_t> bz_face_mapped(max_bz_idx, 0u), ir_face_mapped(max_ir_idx, 0u);

  size_t face_idx=0;
  verbose_update("Combine ", i30.size(), " 3:0 normals and plane-points");
  for (auto i: i30) for (int j: i)
    if (0==bz_face_mapped[j]){
      all_norms.set(face_idx, bznormals.extract(j));
      all_point.set(face_idx,  bzpoints.extract(j));
      bz_face_mapped[j] = ++face_idx; // so that bz_face_mapped is the index+1
    }
  verbose_update("Combine ", i21.size(), " 2:1 normals and plane-points");
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
  verbose_update("Combine ", i12.size(), " 1:2 normals and plane-points");
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
  verbose_update("Combine ", i03.size(), " 0:3 normals and plane-points");
  for (size_t i=0; i<i03.size(); ++i) for (size_t j=0; j<3u; ++j){
    if (0==ir_face_mapped[i03.getvalue(i,j)]){
      all_norms.set(face_idx, irnormals.extract(i03.getvalue(i,j)));
      all_point.set(face_idx, 0.0*all_point.extract(i03.getvalue(i,j)));
      ir_face_mapped[i03.getvalue(i,j)] = ++face_idx;
    }
  }

  verbose_update("Normals and plane-points combined");

  size_t vert_idx=0;
  verbose_update("Combine ", i30.size(), " 3:0 vertices and planes-per-vertex");
  for (size_t i=0; i<i30.size(); ++i){
    all_ijk.insert(bz_face_mapped[i30[i][0]]-1, vert_idx, 0);
    all_ijk.insert(bz_face_mapped[i30[i][1]]-1, vert_idx, 1);
    all_ijk.insert(bz_face_mapped[i30[i][2]]-1, vert_idx, 2);
    all_verts.set(vert_idx++, vertices30.extract(i));
  }
  verbose_update("Combine ", i21.size(), " 2:1 vertices and planes-per-vertex");
  for (size_t i=0; i<i21.size(); ++i){
    all_ijk.insert(bz_face_mapped[i21.getvalue(i,0)]-1, vert_idx, 0);
    all_ijk.insert(bz_face_mapped[i21.getvalue(i,1)]-1, vert_idx, 1);
    all_ijk.insert(ir_face_mapped[i21.getvalue(i,2)]-1, vert_idx, 2);
    all_verts.set(vert_idx++, vertices21.extract(i));
  }
  verbose_update("Combine ", i12.size(), " 1:2 vertices and planes-per-vertex");
  for (size_t i=0; i<i12.size(); ++i){
    all_ijk.insert(bz_face_mapped[i12.getvalue(i,0)]-1, vert_idx, 0);
    all_ijk.insert(ir_face_mapped[i12.getvalue(i,1)]-1, vert_idx, 1);
    all_ijk.insert(ir_face_mapped[i12.getvalue(i,2)]-1, vert_idx, 2);
    all_verts.set(vert_idx++, vertices12.extract(i));
  }
  verbose_update("Combine ", i03.size(), " 0:3 vertices and planes-per-vertex");
  for (size_t i=0; i<i03.size(); ++i){
    all_ijk.insert(ir_face_mapped[i03.getvalue(i,0)]-1, vert_idx, 0);
    all_ijk.insert(ir_face_mapped[i03.getvalue(i,1)]-1, vert_idx, 1);
    all_ijk.insert(ir_face_mapped[i03.getvalue(i,2)]-1, vert_idx, 2);
    all_verts.set(vert_idx++, vertices03.extract(i));
  }
  verbose_update("Vertices and planes-per-vertex combined");
  // four lists now combined into one.
  //      all_ijk   -- (N,3) array of which three planes intersected at a vertex
  //      all_verts -- (N,) vector of intersection vertex locations
  //      all_norms -- (N,) vector of plane normals (indexed by all_ijk)
  //      all_point -- (N,) vector of on-plane points (indexed by all_ijk)

  ArrayVector<bool> keep;
  // Find which vertices are inside the irreducible wedge
  const bool constructing{true}; // here we *are* building-up the irreducible Brillouin zone
  keep = this->isinside_wedge(all_verts, constructing);
  // and pull out those vertices and their intersecting plane indices
  all_verts = extract(all_verts, keep);
  all_ijk   = extract(all_ijk,   keep);

  // it is imperitive that the xyz coordinate system of the irreducible
  // polyhedron is the same as that used by the Brillouin zone polyhedron.
  // Deal with this in a special function elsewhere.
  //
  // Creating a Polyhedron object automatically keeps only unique vertices
  // and facet planes which are polygons, plus finds the vertex to facet
  // and facet to vertex indexing required for, e.g., plotting
  this->set_ir_polyhedron(all_verts, all_point, all_norms);
  // this->set_ir_polyhedron(all_verts, all_point, all_norms, all_ijk);
  verbose_update("Found a ",this->ir_polyhedron.string_repr());
}

void BrillouinZone::voro_search(const int extent){
  std::array<double, 3> bbmin{1e3,1e3,1e3}, bbmax{-1e3,-1e3,-1e3};
  LQVec<int> primtau(this->lattice, make_relative_neighbour_indices(extent));
  size_t ntau = primtau.size();
  std::vector<size_t> perm(ntau);
  std::iota(perm.begin(), perm.end(), 0u); // {0u, 1u, 2u, ..., ntau-1}
  std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b){
    return primtau.norm(a) < primtau.norm(b);
  });
  verbose_update("unsorted primtau\n",primtau.to_string(norm(primtau)));
  verbose_update("permutation:",perm);
  primtau.permute(perm);
  verbose_update("sorted primtau\n",primtau.to_string(norm(primtau)));
  LQVec<double> convtau = transform_from_primitive(this->outerlattice, primtau);
  verbose_update("convtau\n",convtau.to_string(norm(convtau)));
  // the first Brillouin zone polyhedron will be expressed in absolute units
  // in the xyz frame of the conventional reciprocal lattice
  ArrayVector<double> tau = convtau.get_xyz();
  double tij;
  for (size_t i=0; i<ntau; ++i)
    for (size_t j=0; j<3u; ++j){
      tij = tau.getvalue(i, j);
      if (tij < bbmin[j]) bbmin[j] = tij;
      if (tij > bbmax[j]) bbmax[j] = tij;
  }
  // create an initialize the voro++ voronoicell object
  Polyhedron voronoi = polyhedron_box(bbmin, bbmax);
  // and then use the reciprocal lattice points to subdivide the cell until
  // only the first Brillouin zone is left:
  Polyhedron firstbz = Polyhedron::bisect(voronoi, tau/norm(tau), tau/2.0);
  this->polyhedron = firstbz;
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
  // Yes. This is, of course, the Binomal coefficient (ntau, 3)
  debug_update_if(ntau<3, "vertex_search:: Not enough taus for binomial coefficient");
  unsigned long long bc = binomial_coefficient(static_cast<unsigned>(ntau), 3u);
  if (bc != static_cast<unsigned long long>(ntocheck))
    throw std::runtime_error("Mistake calculating the binomial coefficient?!");

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
      in_ijk.set(in_cnt++, all_ijk.data(i));
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
  bool *vertisunique = nullptr;
  vertisunique = new bool[in_cnt]();
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
  int unqcnt=0, *unqidx = nullptr;
  unqidx = new int[in_cnt]();
  for (int i=0; i<in_cnt; i++) if (vertisunique[i]) unqidx[unqcnt++]=i;

  //printf("%d unique vertices\n",unqcnt);

  // create a mapping which holds the indices of all equivalent vertices
  // so unqidxmap[1,:] are all of the indices which equal the first unique vertex
  int *unqidxmap = nullptr;
  unqidxmap = new int[in_cnt*unqcnt]();
  int *numidxmap = nullptr;
  numidxmap = new int[unqcnt]();
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
  double *unqlenmap = nullptr;
  unqlenmap = new double[maxequiv*unqcnt](); // no need to allocate in_cnt*unqcnt memory when maxequiv is known
  for (int i=0; i<unqcnt; i++){
    for (int j=0; j<numidxmap[i]; j++){
      unqlenmap[i*maxequiv + j] = 0;
      for (int k=0; k<3; k++){
        unqlenmap[i*maxequiv + j] += halftau.norm( in_ijk.getvalue(unqidxmap[i*in_cnt+j],k));
      }
    }
  }
  // use the "length" information to select which equivalent vertex we should keep
  int *minequividx = nullptr;
  minequividx = new int[unqcnt]();
  double *minequivlen = nullptr;
  minequivlen = new double[unqcnt]();
  for (int i=0; i<unqcnt; i++){
    minequivlen[i] = (std::numeric_limits<double>::max)(); // better than 1./0.
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
      unq_vrt.set( i, in_verts.data(minequividx[i]) );
      unq_ijk.set( i, in_ijk.data(minequividx[i]) );
  }
  delete[] minequividx;

  // determine which of the taus actually contribute to at least one vertex
  int ncontrib=0, *contrib = nullptr;
  contrib = new int[ntau]();
  for (int i=0; i<ntau; i++)
  for (int j=0; j<unqcnt; j++)
  if ( unq_ijk.getvalue(j,0) == i || unq_ijk.getvalue(j,1) == i || unq_ijk.getvalue(j,2) == i ){
    contrib[ncontrib++]=i;
    break;
  }
  ArrayVector<int> faces(3,ncontrib);
  for (int i=0; i<ncontrib; i++) faces.set(i, tau.data(contrib[i]));

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

  // store the unique vertices, on-face points (half of 'faces'), and faces_per_vertex
  LQVec<int> half_tau(this->lattice, faces);
  // to ensure that a common coordinate system is used by all polyhedron,
  // use an assignment function to set the bz polyhedron
  this->set_polyhedron(unq_vrt, half_tau/2.0);
  // this->set_polyhedron(unq_vrt, half_tau/2.0, fpv);

  delete[] contrib;
}

template<typename T> ArrayVector<bool> BrillouinZone::isinside(const LQVec<T>& p) const {
  bool isouter = this->outerlattice.issame(p.get_lattice());
  bool isinner = this->lattice.issame(p.get_lattice());
  if (!(isouter||isinner))
    throw std::runtime_error("Q points provided to BrillouinZone::isinside must be in the standard or primitive lattice used to define the BrillouinZone object");
  ArrayVector<bool> out(1u, p.size());
  LQVec<double> points, normals;
  if (isouter){
    points = this->get_points();
    normals = this->get_normals();
  } else {
    points = this->get_primitive_points();
    normals = this->get_primitive_normals();
  }
  for (size_t i=0; i<p.size(); ++i)
    out.insert( dot(normals, p.get(i)-points).all_approx(Comp::le,0.), i );
  return out;
}
template<typename T> std::vector<bool> BrillouinZone::isinside_std(const LQVec<T>& p) const {
  bool isouter = this->outerlattice.issame(p.get_lattice());
  bool isinner = this->lattice.issame(p.get_lattice());
  if (!(isouter||isinner))
    throw std::runtime_error("Q points must be in the standard or primitive lattice");
  std::vector<bool> out(p.size(), true);
  LQVec<double> points, normals;
  if (isouter){
    points = this->get_points();
    normals = this->get_normals();
  } else {
    points = this->get_primitive_points();
    normals = this->get_primitive_normals();
  }
  for (size_t i=0; i<p.size(); ++i)
    out[i] = dot(normals, p.get(i)-points).all_approx(Comp::le, 0.);
  return out;
}

template<typename T> ArrayVector<bool> BrillouinZone::isinside_wedge(const LQVec<T> &p, const bool constructing) const {
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
  if (normals.size()){
    // If a pointgroup has inversion symmetry then for every point, p, there is
    // an equivalent point, -p. This indicates that a point p is already in
    // the irreducible wedge if it has n̂ᵢ⋅p ≥ 0 for all irredudible-bounding-
    // plane normals, ̂nᵢ, *or* the opposite -- n̂ᵢ⋅ p ≤ 0 for all ̂nᵢ.
    // Since we are interested in enforcing a smallest-possible irreducible
    // Brillouin zone, we want to exclude the all n̂ᵢ⋅ p ≤ 0 half of the
    // reciprocal wedge precisely because they are equivalent.
    // It is only in the case where a pointgroup *does not have* space-inversion
    // symmetry and time-reversal symmetry is to be excluded that we can allow
    // n̂ᵢ⋅ p ≤ 0 to be a valid in-wedge solution.
    // The ArrayVector all_approx method has a special switched execution path
    // for checking whether all values are ≤ or ≥ a value simultaneously

    // when constructing the irreducible Brillouin zone we need to only consider
    // the ≥0 case so that we end up with a convex polyhedron. The ir_polyhedron
    // accessor method mirrors the half-polyhedron in this case, so when
    // identifying whether a point is inside of the irreducible Brillouin zone
    // we must allow for the ≤0 case as well.
    Comp cmp = constructing||this->no_ir_mirroring ? Comp::ge : Comp::le_ge;
    for (size_t i=0; i<p.size(); ++i)
      out.insert(dot(normals, p.get(i)).all_approx(cmp,0.), i);
  } else {
    for (size_t i=0; i<p.size(); ++i)
      out.insert(true, i); // with no normals *all* points are "inside" the wedge
  }
  return out;
}

template<typename T> std::vector<bool> BrillouinZone::isinside_wedge_std(const LQVec<T> &p, const bool constructing) const {
  bool isouter = this->outerlattice.issame(p.get_lattice());
  bool isinner = this->lattice.issame(p.get_lattice());
  if (!(isouter||isinner))
    throw std::runtime_error("Q points provided to BrillouinZone::isinside_wedge must be in the standard or primitive lattice used to define the BrillouinZone object");
  std::vector<bool> out(p.size(), true); // default to true in case of no normals when all points are 'inside'
  LQVec<double> normals;
  if (isouter){
    verbose_update("Check points against standard wedge normals");
    normals = this->get_ir_wedge_normals();
  } else {
    verbose_update("Check points against primitive wedge normals");
    normals = this->get_primitive_ir_wedge_normals();
  }
  if (normals.size()){
    verbose_update("Checking if points are within wedge\n",normals.to_string());
    Comp cmp = constructing||this->no_ir_mirroring ? Comp::ge : Comp::le_ge;
    for (size_t i=0; i<p.size(); ++i)
      out[i] = dot(normals, p.get(i)).all_approx(cmp,0.);
    verbose_update("Comparison to normals\n",normals.to_string()," via ",cmp," yields\n",out);
  }
  return out;
}

bool BrillouinZone::moveinto(const LQVec<double>& Q, LQVec<double>& q, LQVec<int>& tau, const int threads) const {
  verbose_update("BrillouinZone::moveinto called with ",threads," threads");
  omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
  bool already_same = this->lattice.issame(Q.get_lattice());
  LQVec<double> Qprim(this->lattice), qprim(this->lattice);
  LQVec<int> tauprim(this->lattice);
  PrimitiveTransform PT(this->outerlattice.get_spacegroup_object().get_bravais_type());
  bool transform_needed = ( PT.does_anything() && this->outerlattice.issame(Q.get_lattice()) );
  if (!(already_same || transform_needed))
    throw std::runtime_error("Q points provided to BrillouinZone::moveinto must be in the standard or primitive lattice used to define the BrillouinZone object");

  if (transform_needed)  Qprim = transform_to_primitive(this->outerlattice,Q);
  const LQVec<double> & Qsl = transform_needed ? Qprim : Q;
  LQVec<double> & qsl = transform_needed ? qprim : q;
  LQVec<int> & tausl = transform_needed? tauprim : tau;

  // the face centre points and normals in the primitive lattice:
  auto points = this->get_primitive_points();
  auto normals = this->get_primitive_normals();
  normals = normals/norm(normals); // ensure they're normalised
  auto taus = (2.0*points).round();
  auto taulen = norm(taus);
  size_t max_count = taus.size();
  // ensure that qsl and tausl can hold each qi and taui
  qsl.resize(Qsl.size());
  tausl.resize(Qsl.size());
  long long snQ = unsigned_to_signed<long long, size_t>(Qsl.size());
#pragma omp parallel for default(none)\
shared(Qsl, tausl, qsl, points, normals, taus, taulen, snQ, max_count)\
schedule(dynamic)
  for (long long si=0; si<snQ; si++){
    size_t i = signed_to_unsigned<size_t, long long>(si);
    LQVec<int> taui = Qsl.get(i).round();
    LQVec<double> qi = Qsl.get(i) - taui;
    LQVec<int> last_shift = taui;
    size_t count{0};
    while (count++ < max_count && dot(normals, qi-points).any_approx(Comp::gt,0.)){
      auto qi_dot_normals = dot(qi , normals);
      auto Nhkl = (qi_dot_normals/taulen).round().to_std();
      auto qidn = qi_dot_normals.to_std();
      if (std::any_of(Nhkl.begin(), Nhkl.end(), [](int a){return a > 0;})){
        int maxnm{0};
        size_t maxat{0};
        for (size_t j=0; j<Nhkl.size(); ++j)
                                                         // protect against oscillating by ±τ
        if (Nhkl[j]>0 && Nhkl[j]>=maxnm && (0==maxnm || (norm(taus[j]+last_shift).all_approx(Comp::gt, 0.) && qidn[j]>qidn[maxat]))){
          maxnm = Nhkl[maxat=j];
        }
        qi -= taus[maxat] * (double)(maxnm); // ensure we subtract LQVec<double>
        taui += taus[maxat] * maxnm; // but add LQVec<int>
        last_shift = taus[maxat] * maxnm;
      }
    }
    qsl.set(i, qi);
    tausl.set(i, taui);
  }
  if (transform_needed){ // then we need to transform back q and tau
    q   = transform_from_primitive(this->outerlattice,qsl);
    tau = transform_from_primitive(this->outerlattice,tausl);
  }
  ArrayVector<bool> allinside = this->isinside(q);
  if (!allinside.all_true()){
    for (size_t i=0; i<Q.size(); ++i) if (!allinside.getvalue(i)){
      info_update("Q  =",Q.to_string(i)  ," tau  =",tau.to_string(i)  ," q  =",q.to_string(i));
      info_update("Qsl=",Qsl.to_string(i)," tausl=",tausl.to_string(i)," qsl=",qsl.to_string(i),"\n");
    }
    throw std::runtime_error("Not all points inside Brillouin zone");
    // return false;
  }
  return true; // otherwise an error has been thrown
}
bool BrillouinZone::ir_moveinto(
  const LQVec<double>& Q, LQVec<double>& q, LQVec<int>& tau,
  std::vector<size_t>& Ridx, std::vector<size_t>& invRidx, const int threads
) const {
  verbose_update("BrillouinZone::ir_moveinto called with ",threads," threads");
  omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
  /* The Pointgroup symmetry information comes from, effectively, spglib which
  has all rotation matrices defined in the conventional unit cell -- which is
  our `outerlattice`. Consequently we must work in the outerlattice here.  */
  if (!this->outerlattice.issame(Q.get_lattice()))
    throw std::runtime_error("Q points provided to ir_moveinto must be in the standard lattice used to define the BrillouinZone object");
  // get the PointSymmetry object, containing all operations
  PointSymmetry psym = this->get_pointgroup_symmetry();
  // ensure q, tau, and Rm can hold one for each Q.
  size_t nQ = Q.size();
  q.resize(nQ);
  tau.resize(nQ);
  Ridx.resize(nQ);
  invRidx.resize(nQ);
  // find q₁ₛₜ in the first Brillouin zone and τ ∈ [reciprocal lattice vectors]
  // such that Q = q₁ₛₜ + τ
  this->moveinto(Q, q, tau, threads);
  // by chance some first Bz points are likely already in the IR-Bz:
  std::vector<bool> in_ir = this->isinside_wedge_std(q);
  //LQVec<double> qj(Q.get_lattice(), 1u);
  auto lat = Q.get_lattice();
  // OpenMP 2 (VS) doesn't like unsigned loop counters
  size_t n_outside{0};
  long long snQ = unsigned_to_signed<long long, size_t>(nQ);
  #pragma omp parallel for default(none) shared(psym, Ridx, invRidx, q, in_ir, lat, snQ) reduction(+:n_outside) schedule(dynamic)
  for (long long si=0; si<snQ; ++si){
    size_t i = signed_to_unsigned<size_t, long long>(si);
    // any q already in the irreducible zone need no rotation → identity, but we need to find the index of E
    bool outside=!in_ir[i];
    if (outside){
      // for others find the jᵗʰ operation which moves qᵢ into the irreducible zone
      LQVec<double> qj(lat, 1u); // a place to hold the multiplication result
      for (size_t j=0; j<psym.size(); ++j) if (outside) {
        // The point symmetry matrices relate *real space* vectors! We must use
        // their transposes' to rotate reciprocal space vectors.
        multiply_matrix_vector(qj.data(0), transpose(psym.get(j)).data(), q.data(i));
        if ( this->isinside_wedge_std(qj)[0] ){ /* store the result */
          // and (Rⱼᵀ)⁻¹ ∈ G, such that Qᵢ = (Rⱼᵀ)⁻¹⋅qᵢᵣ + τᵢ.
          q.set(i, qj); // keep Rⱼᵀ⋅qᵢ as qᵢᵣ
          invRidx[i] = j; // Rⱼ *is* the inverse of what we want for ouput
          Ridx[i] = psym.get_inverse_index(j); // find the index of Rⱼ⁻¹
          outside = false;
        }
      }
    } else {
      invRidx[i] = Ridx[i] = psym.find_index({1,0,0, 0,1,0, 0,0,1});
    }
    if (outside) {
      ++n_outside;
      in_ir[i] = false;
    }
  }
  if (n_outside > 0) for (size_t i=0; i<nQ; ++i) if (!in_ir[i]){
    std::string msg = "Q = " + Q.to_string(i);
    msg += " is outside of the irreducible BrillouinZone ";
    msg += " : tau = " + tau.to_string(i) + " , q = " + q.to_string(i);
    throw std::runtime_error(msg);
    return false;
  }
  return true; // otherwise we hit the runtime error above
}
bool BrillouinZone::ir_moveinto_wedge(const LQVec<double>& Q, LQVec<double>& q, std::vector<std::array<int,9>>& R, const int threads) const {
  omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
  /* The Pointgroup symmetry information comes from, effectively, spglib which
  has all rotation matrices defined in the conventional unit cell -- which is
  our `outerlattice`. Consequently we must work in the outerlattice here.  */
  if (!this->outerlattice.issame(Q.get_lattice()))
    throw std::runtime_error("Q points provided to ir_moveinto must be in the standard lattice used to define the BrillouinZone object");
  // get the PointSymmetry object, containing all operations
  PointSymmetry psym = this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
  // ensure q and R can hold one for each Q.
  size_t nQ = Q.size();
  q.resize(nQ);
  R.resize(nQ);
  // by chance some first Bz points are likely already in the IR-wedge:
  std::vector<bool> in_ir = this->isinside_wedge_std(Q);
  //LQVec<double> qj(Q.get_lattice(), 1u);
  auto lat = Q.get_lattice();
  // OpenMP 2 (VS) doesn't like unsigned loop counters
  size_t n_outside{0};
  long long snQ = unsigned_to_signed<long long, size_t>(nQ);
  #pragma omp parallel for default(none) shared(psym, R, q, Q, in_ir, lat, snQ) reduction(+:n_outside) schedule(dynamic)
  for (long long si=0; si<snQ; ++si){
    size_t i = signed_to_unsigned<size_t, long long>(si);
    // any q already in the irreducible zone need no rotation → identity
    bool outside=!in_ir[i];
    if (outside){
      // for others find the jᵗʰ operation which moves qᵢ into the irreducible zone
      LQVec<double> qj(lat, 1u); // a place to hold the multiplication result
      for (size_t j=0; j<psym.size(); ++j) if (outside) {
        // The point symmetry matrices relate *real space* vectors! We must use
        // their transposes' to rotate reciprocal space vectors.
        multiply_matrix_vector(qj.data(0), transpose(psym.get(j)).data(), Q.data(i));
        if ( this->isinside_wedge_std(qj)[0] ){ /* store the result */
          q.set(i, qj); // keep Rⱼᵀ⋅Qᵢ as qᵢᵣ
          R[i] = transpose(psym.get_inverse(j)); // and (Rⱼᵀ)⁻¹ ∈ G, such that Q = (Rⱼᵀ)⁻¹⋅qᵢᵣ
          outside = false;
        }
      }
    } else {
      q.set(i, Q.extract(i));
      R[i] = {1,0,0, 0,1,0, 0,0,1};
    }
    if (outside) {
      ++n_outside;
      in_ir[i] = false;
    }
  }
  if (n_outside > 0) for (size_t i=0; i<nQ; ++i) if (!in_ir[i]){
    std::string msg = "Q = " + Q.to_string(i);
    msg += " is outside of the irreducible reciprocal space wedge ";
    msg += " , irQ = " + q.to_string(i);
    throw std::runtime_error(msg);
    return false;
  }
  return true; // otherwise we hit the runtime error above
}


bool three_plane_intersection(const LQVec<double>& n,                // plane normals
                              const LQVec<double>& p,                // a point on each plane
                              const ArrayVector<double>& xyz,        // the plane normals in a inverse Angstrom orthonormal coordinate system
                              const int i, const int j, const int k, // indices into n, p, xyz
                              LQVec<double>& iat,                    // output storage array of intersection points
                              const int idx)                         // the index where the intersection is inserted if found
                              {
  // we need to check whether the matrix formed by the orthonormal-frame components of the three planes is nearly-singular
  double *M = nullptr;
  M = new double[9];
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

    iat.set(idx, tmp.data() );
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
    intersect.set(idx, tmp.data());
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
    intersect.set(idx, tmp.data());
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
    intersect.set(idx, tmp.data());
    return true;
  }
  return false;
}
bool intersect_at(const LQVec<double>& ni,
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
    intersect.set(idx, 0*intersect.get(idx));
    return true;
  }
  return false;
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
  inv->set(store_at, v->data(idx));
  return true;
}
