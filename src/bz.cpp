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

#include "bz.hpp"
using namespace brille;

LQVec<double> BrillouinZone::get_ir_polyhedron_wedge_normals() const {
  auto ir_n = this->get_ir_normals();
  auto ir_p = this->get_ir_points();
  auto bz_n = this->get_normals();
  auto bz_p = this->get_points();
  for (ind_t i=0; i<bz_n.size(0); ++i){
    // check if the irBZ face point is on a first BZ zone face too
    // FIXME Add tolerance here
    auto not_bz = dot(bz_n.view(i), ir_p - bz_p.view(i)).is(brille::cmp::neq, 0.);
    if (!not_bz.any())
      throw std::runtime_error("Why are all of the points on the first BZ surface?!");
    ir_n = ir_n.extract(not_bz);
    ir_p = ir_p.extract(not_bz);
  }
  // It is possible that, lacking inversion symmetry, we have found an
  // irreducible polyhedron which comprises two convex polyhedra that
  // are mutually inverse. In such a case for every normal in ir_n there is also
  // the opposite also in ir_n, and so ir_n defines a wedge and an anti-wedge in
  // which no finite point can be inside
  if (ir_n.size(0)%2 == 0 /* [lacks inversion] or this->no_ir_mirroring? */){
    std::vector<bool> no_inverse(ir_n.size(0), true), keep(ir_n.size(0), true);
    for (ind_t i=0; i<ir_n.size(0)-1; ++i) if (no_inverse[i])
    for (ind_t j=i+1; j<ir_n.size(0); ++j) if (no_inverse[j])
    if ((ir_n.view(i)+ir_n.view(j)).all(brille::cmp::eq, 0.)) {
      // FIXME Consider adding tolerance                             ^^^
      no_inverse[i] = no_inverse[j] = keep[j] = false;
      break;
    }
    if (0==std::count(no_inverse.begin(), no_inverse.end(), true))
      ir_n = ir_n.extract(keep);
  }
  // the remaining irBZ faces are the irreducible reciprocal space wedge
  // which we store with the opposite sign in this->ir_wedge_normals
  return -ir_n;
}

BrillouinZone::poly_t BrillouinZone::get_polyhedron() const {return first_poly;}
BrillouinZone::poly_t BrillouinZone::get_ir_polyhedron(const bool true_ir) const {
  // If the ir polyhedron fills the first BZ using the symmetry operations
  // of the pointgroup (plus time inversion, if included), then no ir_mirroring
  // is to be performed. In this case or if the 'true_ir' was requested, return
  // the computed ir_polyhedron
  if (no_ir_mirroring || !true_ir) return ir_poly;
  // Otherwise, it is only half of the true (non-convex)
  // irreducible polyhedron, so add the lack of inversion symmetry explicitly.
  return ir_poly + ir_poly.mirror();
}
bool BrillouinZone::check_ir_polyhedron(){
  profile_update("Start BrillouinZone::check_ir_polyhedron");
  this->check_if_mirroring_needed(); // move this to end of wedge_brute_force?
  auto ps = outerlattice.get_pointgroup_symmetry(time_reversal?1:0);
  auto volume_goal = first_poly.volume() / static_cast<double>(ps.size());
  if (!brille::approx::scalar(ir_poly.volume(), volume_goal, 1000)){ // setting tol=100 moves the relative difference allowed to ~2e-9
    info_update("The current 'irreducible' polyhedron has the wrong volume since ",ir_poly.volume()," != ",volume_goal);
    info_update(std::abs(ir_poly.volume() - volume_goal) / (ir_poly.volume() + volume_goal));
    return false;
  }
  ps = ps.higher(1);
  for (size_t i=0; i < ps.size(); ++i) {
    auto rotated = ir_poly.apply(ps, i);
    if (ir_poly.intersects(rotated)){
      debug_update("The trial irreducible polyhedron intersects itself.");
      return false;
    }
  }
  // volume is right and no intersections
  profile_update("  End BrillouinZone::check_ir_polyhedron");
  return true;
}

// first Brillouin zone
BrillouinZone::poly_t::vertex_t BrillouinZone::get_vertices() const {
  return first_poly.vertices();
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_primitive_vertices() const {
  auto v = this->get_vertices();
  if (this->isprimitive()) v = transform_to_primitive(outerlattice, v);
  return v;
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_points() const {
  return first_poly.face_points();
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_primitive_points() const {
  auto p = this->get_points();
  if (this->isprimitive()) p = transform_to_primitive(this->outerlattice, p);
  return p;
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_normals() const {
  return first_poly.face_normals();
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_primitive_normals() const {
  auto n = this->get_normals();
  if (this->isprimitive()) n = transform_to_primitive(this->outerlattice, n);
  return n;
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_half_edges() const{
  return first_poly.half_edges();
}
BrillouinZone::poly_t::faces_t BrillouinZone::get_faces_per_vertex() const {
  return first_poly.faces_per_vertex();
}
BrillouinZone::poly_t::faces_t BrillouinZone::get_vertices_per_face() const {
  return first_poly.faces();
}
// irreducible first Brillouin zone
BrillouinZone::poly_t::vertex_t BrillouinZone::get_ir_vertices() const {
  return ir_poly.vertices();
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_ir_primitive_vertices() const {
  auto v = this->get_ir_vertices();
  if (this->isprimitive()) v = transform_to_primitive(this->outerlattice, v);
  return v;
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_ir_points() const {
  return ir_poly.face_points();
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_ir_primitive_points() const {
  auto p = this->get_ir_points();
  if (this->isprimitive()) p = transform_to_primitive(this->outerlattice, p);
  return p;
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_ir_normals() const {
  return ir_poly.face_normals();
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_ir_primitive_normals() const {
  auto n = this->get_ir_normals();
  if (this->isprimitive()) n = transform_to_primitive(this->outerlattice, n);
  return n;
}
BrillouinZone::poly_t::faces_t BrillouinZone::get_ir_faces_per_vertex() const {
  return ir_poly.faces_per_vertex();
}
BrillouinZone::poly_t::faces_t BrillouinZone::get_ir_vertices_per_face() const {
  return ir_poly.faces();
}

// irreducible reciprocal space
LQVec<double> BrillouinZone::get_ir_wedge_normals() const {
  LQVec<double> out(outerlattice, 0u);
  if (ir_wedge_normals.size(0))
    out = LQVec<double>(outerlattice, ir_wedge_normals);
  return out;
}
//
LQVec<double> BrillouinZone::get_primitive_ir_wedge_normals(void) const {
  LQVec<double> lqwn(outerlattice, 0u);
  if (ir_wedge_normals.size(0)){
    lqwn = LQVec<double>(outerlattice, ir_wedge_normals);
    if (this->isprimitive())
      lqwn = transform_to_primitive(outerlattice, lqwn);
  }
  return lqwn;
}


void BrillouinZone::print() const {
  std::string msg = "BrillouinZone with ";
  msg += std::to_string(this->vertices_count()) + " vertices and ";
  msg += std::to_string(this->faces_count()) + " faces";
  std::cout << msg << std::endl;
}

//
//void BrillouinZone::irreducible_vertex_search(){
//  using namespace brille;
//  /* We need to check for three-plane intersections for all combinations of two
//     1st Brillouin zone planes and one irreducible reciprocal space normal and
//     two irreducible reciprocal space normals and one 1st Brillouin zone plane.
//  */
//  // TODO Why can't this use the polyhedron cutting algorithm?
//  ind_t Nbz = this->get_normals().size(0);
//  ind_t Nir = this->ir_wedge_normals.size(0);
//
//  if (0==Nir){
//    this->ir_polyhedron = this->polyhedron;
//    return;
//  }
//
//  // for which there are M*(N*(N-1))/2 + N*(M*(M-1))/2 total possible combinations
//  ind_t n21 = ((Nbz*(Nbz-1))>>1)*Nir;
//  ind_t n12 = ((Nir*(Nir-1))>>1)*Nbz;
//  ind_t n03 = 0;
//  for (ind_t i=2; i<Nir; ++i) n03 += (i*(i-1))>>1;
//  verbose_update("Checking {",n21,", ",n12,", ",n03,"} {2:1, 1:2, 0:3} zone:wedge 3-plane intersection points");
//
//  auto bznormals = this->get_normals();
//  auto bzpoints = this->get_points();
//
//  // We will create a polyhedron using (some of) these normals. It is imperitive
//  // that the polyhedron normals point *outwards* from the centre of the polyhedron
//  // while we've thus far defined the wedge normals to point *into* the irreducible
//  // reciprocal space wedge.
//  auto irnormals = -1.0*this->get_ir_wedge_normals();
//  auto vertices30 = this->get_vertices();
//  std::vector<std::vector<int>> i30 = this->get_faces_per_vertex();
//
//  LQVec<double> vertices21(bznormals.get_lattice(), n21);
//  LQVec<double> vertices12(bznormals.get_lattice(), n12);
//  LQVec<double> vertices03(bznormals.get_lattice(), n03);
//  bArray<int> i21(n21,3);
//  bArray<int> i12(n12,3);
//  bArray<int> i03(n03,3);
//
//  int c21=0, c12=0, c03=0;
//  if (n21){ // protect against Nbz=0, since size_t(0)-1 = 4294967294 or 18446744073709551615 if its 32- or 64-bit
//    for (ind_t i=0  ; i<(Nbz-1); ++i)
//    for (ind_t j=i+1; j< Nbz   ; ++j)
//    for (ind_t k=0  ; k< Nir   ; ++k)
//    if (brille::intersect_at(
//      bznormals.view(i), bzpoints.view(i),
//      bznormals.view(j), bzpoints.view(j),
//      irnormals.view(k),                   vertices21, c21)
//    ){
//      i21.val(c21,0)=i; i21.val(c21,1)=j; i21.val(c21++,2)=k;
//    }
//  }
//  if (n12){ // protect against Nir=0, since size_t(0)-1 = 4294967294 or 18446744073709551615 if its 32- or 64-bit
//    for (ind_t i=0  ; i< Nbz   ; ++i)
//    for (ind_t j=0  ; j<(Nir-1); ++j)
//    for (ind_t k=j+1; k< Nir   ; ++k)
//    if (brille::intersect_at(
//      bznormals.view(i), bzpoints.view(i),
//      irnormals.view(j),
//      irnormals.view(k),                   vertices12, c12)
//    ){
//      i12.val(c12,0)=i, i12.val(c12,1)=j, i12.val(c12++,2)=k;
//    }
//  }
//  if (n03){
//    for (ind_t i=0  ; i<(Nir-2); ++i)
//    for (ind_t j=i+1; j<(Nir-1); ++j)
//    for (ind_t k=j+1; k< Nir   ; ++k)
//    if (brille::intersect_at(
//      irnormals.view(i),
//      irnormals.view(j),
//      irnormals.view(k),                   vertices03, c03)
//    ){
//      i03.val(c03,0)=i; i03.val(c03,1)=j; i03.val(c03++,2)=k;
//    }
//  }
//  verbose_update("Intersections found");
//  // make sure we shrink all sets of vertices to just those found!
//  // plus remove any intersection points outside of the first Brillouin zone
//  this->shrink_and_prune_outside(static_cast<size_t>(c21), vertices21, i21);
//  this->shrink_and_prune_outside(static_cast<size_t>(c12), vertices12, i12);
//  this->shrink_and_prune_outside(static_cast<size_t>(c03), vertices03, i03);
//  verbose_update("Intersections pruned");
//  // Now we have four lists of vertices, plus lists of the normal vectors
//  // and on-plane points, which define the three planes that intersect at each
//  // vertex.
//  // We want to combine these lists:
//
//  int max_bz_idx=0, max_ir_idx=0;
//  for (auto i: i30) for (int j: i) if (j > max_bz_idx) max_bz_idx = j;
//  for (ind_t i=0; i<i21.size(0); ++i){
//    if (i21.val(i,0)>max_bz_idx) max_bz_idx=i21.val(i,0);
//    if (i21.val(i,1)>max_bz_idx) max_bz_idx=i21.val(i,1);
//    if (i21.val(i,2)>max_ir_idx) max_ir_idx=i21.val(i,2);
//  }
//  for (ind_t i=0; i<i12.size(0); ++i){
//    if (i12.val(i,0)>max_bz_idx) max_bz_idx=i12.val(i,0);
//    if (i12.val(i,1)>max_ir_idx) max_ir_idx=i12.val(i,1);
//    if (i12.val(i,2)>max_ir_idx) max_ir_idx=i12.val(i,2);
//  }
//  for (ind_t i=0; i<i03.size(0); ++i) for (ind_t j=0; j<3u; ++j){
//    if (i03.val(i,j)>max_ir_idx) max_ir_idx=i03.val(i,j);
//  }
//  max_bz_idx++; max_ir_idx++; // since we count from 0, and need one more element than the highest index.
//  std::vector<bool> bz_face_present(max_bz_idx, false), ir_face_present(max_ir_idx, false);
//  for (auto i: i30) for (int j: i) bz_face_present[j] = true;
//  for (ind_t i=0; i<i21.size(0); ++i){
//    bz_face_present[i21.val(i,0)] = true;
//    bz_face_present[i21.val(i,1)] = true;
//    ir_face_present[i21.val(i,2)] = true;
//  }
//  for (ind_t i=0; i<i12.size(0); ++i){
//    bz_face_present[i12.val(i,0)] = true;
//    ir_face_present[i12.val(i,1)] = true;
//    ir_face_present[i12.val(i,2)] = true;
//  }
//  for (ind_t i=0; i<i03.size(0); ++i) for (ind_t j=0; j<3u; ++j)
//    ir_face_present[i03.val(i,j)] = true;
//  ind_t bz_faces = static_cast<ind_t>(std::count(bz_face_present.begin(), bz_face_present.end(), true));
//  ind_t ir_faces = static_cast<ind_t>(std::count(ir_face_present.begin(), ir_face_present.end(), true));
//
//  ind_t total_verts;
//  total_verts  = vertices30.size(0) + vertices12.size(0);
//  total_verts += vertices21.size(0) + vertices03.size(0);
//
//  LQVec<double> all_verts(bznormals.get_lattice(), total_verts);
//  LQVec<double> all_norms(bznormals.get_lattice(), bz_faces+ir_faces);
//  LQVec<double> all_point(bznormals.get_lattice(), bz_faces+ir_faces);
//  // bArray<int> all_ijk(total_verts,3);
//
//  std::vector<ind_t> bz_face_mapped(max_bz_idx, 0u), ir_face_mapped(max_ir_idx, 0u);
//
//  ind_t face_idx=0;
//  verbose_update("Combine ", i30.size(), " 3:0 normals and plane-points");
//  for (auto i: i30) for (int j: i)
//    if (0==bz_face_mapped[j]){
//      all_norms.set(face_idx, bznormals.view(j));
//      all_point.set(face_idx,  bzpoints.view(j));
//      bz_face_mapped[j] = ++face_idx; // so that bz_face_mapped is the index+1
//    }
//  verbose_update("Combine ", i21.size(0), " 2:1 normals and plane-points");
//  ind_t jdx;
//  for (ind_t i=0; i<i21.size(0); ++i){
//    jdx = i21.val(i,0);
//    if (0==bz_face_mapped[jdx]){
//      all_norms.set(face_idx, bznormals.view(jdx));
//      all_point.set(face_idx,  bzpoints.view(jdx));
//      bz_face_mapped[jdx] = ++face_idx;
//    }
//    jdx = i21.val(i,1);
//    if (0==bz_face_mapped[jdx]){
//      all_norms.set(face_idx, bznormals.view(jdx));
//      all_point.set(face_idx,  bzpoints.view(jdx));
//      bz_face_mapped[jdx] = ++face_idx;
//    }
//    jdx = i21.val(i,2);
//    if (0==ir_face_mapped[jdx]){
//      all_norms.set(face_idx, irnormals.view(jdx));
//      all_point.view(face_idx) *= 0.0;
//      ir_face_mapped[jdx] = ++face_idx;
//    }
//  }
//  verbose_update("Combine ", i12.size(0), " 1:2 normals and plane-points");
//  for (ind_t i=0; i<i12.size(0); ++i){
//    jdx = i12.val(i,0);
//    if (0==bz_face_mapped[jdx]){
//      all_norms.set(face_idx, bznormals.view(jdx));
//      all_point.set(face_idx,  bzpoints.view(jdx));
//      bz_face_mapped[jdx] = ++face_idx;
//    }
//    jdx = i12.val(i,1);
//    if (0==ir_face_mapped[jdx]){
//      all_norms.set(face_idx, irnormals.view(jdx));
//      all_point.view(face_idx) *= 0.0;
//      ir_face_mapped[jdx] = ++face_idx;
//    }
//    jdx = i12.val(i,2);
//    if (0==ir_face_mapped[jdx]){
//      all_norms.set(face_idx, irnormals.view(jdx));
//      all_point.view(face_idx) *= 0.0;
//      ir_face_mapped[jdx] = ++face_idx;
//    }
//  }
//  verbose_update("Combine ", i03.size(0), " 0:3 normals and plane-points");
//  for (ind_t i=0; i<i03.size(0); ++i) for (ind_t j=0; j<3u; ++j){
//    jdx = i03.val(i,j);
//    if (0==ir_face_mapped[jdx]){
//      all_norms.set(face_idx, irnormals.view(jdx));
//      all_point.view(face_idx) *= 0.0;
//      ir_face_mapped[jdx] = ++face_idx;
//    }
//  }
//  verbose_update("Normals and plane-points combined");
//
//  ind_t vert_idx=0;
//  verbose_update("Combine ", i30.size(), " 3:0 vertices and planes-per-vertex");
//  for (ind_t i=0; i<i30.size(); ++i){
//    // all_ijk.val(vert_idx,0) = bz_face_mapped[i30[i][0]]-1;
//    // all_ijk.val(vert_idx,1) = bz_face_mapped[i30[i][1]]-1;
//    // all_ijk.val(vert_idx,2) = bz_face_mapped[i30[i][2]]-1;
//    all_verts.set(vert_idx++, vertices30.view(i));
//  }
//  verbose_update("Combine ", i21.size(0), " 2:1 vertices and planes-per-vertex");
//  for (ind_t i=0; i<i21.size(0); ++i){
//    // all_ijk.val(vert_idx,0) = bz_face_mapped[i21.val(i,0)]-1;
//    // all_ijk.val(vert_idx,1) = bz_face_mapped[i21.val(i,1)]-1;
//    // all_ijk.val(vert_idx,2) = ir_face_mapped[i21.val(i,2)]-1;
//    all_verts.set(vert_idx++, vertices21.view(i));
//  }
//  verbose_update("Combine ", i12.size(0), " 1:2 vertices and planes-per-vertex");
//  for (ind_t i=0; i<i12.size(0); ++i){
//    // all_ijk.val(vert_idx,0) = bz_face_mapped[i12.val(i,0)]-1;
//    // all_ijk.val(vert_idx,1) = ir_face_mapped[i12.val(i,1)]-1;
//    // all_ijk.val(vert_idx,2) = ir_face_mapped[i12.val(i,2)]-1;
//    all_verts.set(vert_idx++, vertices12.view(i));
//  }
//  verbose_update("Combine ", i03.size(0), " 0:3 vertices and planes-per-vertex");
//  for (ind_t i=0; i<i03.size(0); ++i){
//    // all_ijk.val(vert_idx,0) = ir_face_mapped[i03.val(i,0)]-1;
//    // all_ijk.val(vert_idx,1) = ir_face_mapped[i03.val(i,1)]-1;
//    // all_ijk.val(vert_idx,2) = ir_face_mapped[i03.val(i,2)]-1;
//    all_verts.set(vert_idx++, vertices03.view(i));
//  }
//  verbose_update("Vertices and planes-per-vertex combined");
//  // four lists now combined into one.
//  //      all_ijk   -- (N,3) array of which three planes intersected at a vertex
//  //      all_verts -- (N,) vector of intersection vertex locations
//  //      all_norms -- (N,) vector of plane normals (indexed by all_ijk)
//  //      all_point -- (N,) vector of on-plane points (indexed by all_ijk)
//
//  // Find which vertices are inside the irreducible wedge
//  const bool constructing{true}; // here we *are* building-up the irreducible Brillouin zone
//  auto keep = this->isinside_wedge(all_verts, constructing);
//  // and pull out those vertices and their intersecting plane indices
//  all_verts = all_verts.extract(keep);
//  // all_ijk   = all_ijk.extract(keep);
//
//  // it is imperitive that the xyz coordinate system of the irreducible
//  // polyhedron is the same as that used by the Brillouin zone polyhedron.
//  // Deal with this in a special function elsewhere.
//  //
//  // Creating a Polyhedron object automatically keeps only unique vertices
//  // and facet planes which are polygons, plus finds the vertex to facet
//  // and facet to vertex indexing required for, e.g., plotting
//  this->set_ir_polyhedron(all_verts, all_point, all_norms);
//  // this->set_ir_polyhedron(all_verts, all_point, all_norms, all_ijk);
//  verbose_update("Found a ",this->ir_polyhedron.string_repr());
//}

void BrillouinZone::voro_search(const int extent){
  profile_update("Start BrillouinZone::voro_search with ",extent," extent");
  using namespace brille;
  LQVec<int> primtau(lattice, make_relative_neighbour_indices(extent));
  std::vector<ind_t> perm(primtau.size(0));
  std::iota(perm.begin(), perm.end(), 0u); // {0u, 1u, 2u, ..., ntau-1}
  std::sort(perm.begin(), perm.end(), [&](ind_t a, ind_t b){
    return primtau.norm(a) < primtau.norm(b);
  });
  primtau.permute(perm);
//  verbose_update("sorted primtau\n",cat(1,primtau,norm(primtau)).to_string());
  // the first Brillouin zone polyhedron will be expressed in absolute units
  // in the xyz frame of the conventional reciprocal lattice
  auto tau = transform_from_primitive(outerlattice, primtau);
  auto half_tau = tau / 2.0;
  auto [plane_b, plane_c] = plane_points_from_normal(tau / norm(tau), half_tau);
  // and then use the reciprocal lattice points to subdivide the cell until
  // only the first Brillouin zone is left:
  auto box = bounding_box<double>(1.0 * tau);
  first_poly = box.cut(half_tau, plane_b, plane_c);
  profile_update("  End BrillouinZone::voro_search with ",extent," extent");
}

//
//// moved back from the header file now that their template parameters are lost again
//double
//brille::normals_matrix_determinant(
//  const LQVec<double>& a,
//  const LQVec<double>& b,
//  const LQVec<double>& c
//){
//  std::vector<double> metric(9);
//  std::vector<double> ax{a.get_xyz().to_std()}, bx{b.get_xyz().to_std()}, cx{c.get_xyz().to_std()};
//  for (size_t i=0; i<3; ++i){
//    metric[0+i] = ax[i];
//    metric[3+i] = bx[i];
//    metric[6+i] = cx[i];
//  }
//  return brille::utils::matrix_determinant(metric.data());
//}
//
//bool
//brille::intersect_at(
//  const LQVec<double>& n_i, const LQVec<double>& p_i,
//  const LQVec<double>& n_j, const LQVec<double>& p_j,
//  const LQVec<double>& n_k, const LQVec<double>& p_k,
//  LQVec<double>& intersect, const int idx
//){
//  double detM = brille::normals_matrix_determinant(n_i,n_j,n_k);
//  if (brille::approx::scalar(detM, 0.)) return false;
//  auto tmp = cross(n_j,n_k)*dot(p_i,n_i) + cross(n_k,n_i)*dot(p_j,n_j) + cross(n_i,n_j)*dot(p_k,n_k);
//  tmp /= detM;
//  intersect.set(idx, tmp);
//  return true;
//}
//
//bool
//brille::intersect_at(
//  const LQVec<double>& n_i, const LQVec<double>& p_i,
//  const LQVec<double>& n_j, const LQVec<double>& p_j,
//  const LQVec<double>& n_k,
//  LQVec<double>& intersect, const int idx
//){
//  double detM = brille::normals_matrix_determinant(n_i,n_j,n_k);
//  if (brille::approx::scalar(detM, 0.)) return false;
//  auto tmp = cross(n_j,n_k)*dot(p_i,n_i) + cross(n_k,n_i)*dot(p_j,n_j);
//  tmp /= detM;
//  intersect.set(idx, tmp);
//  return true;
//}
//
//bool
//brille::intersect_at(
//  const LQVec<double>& n_i, const LQVec<double>& p_i,
//  const LQVec<double>& n_j,
//  const LQVec<double>& n_k,
//  LQVec<double>& intersect, const int idx
//){
//  double detM = brille::normals_matrix_determinant(n_i,n_j,n_k);
//  if (brille::approx::scalar(detM, 0.)) return false;
//  auto tmp = cross(n_j,n_k)*dot(p_i,n_i);
//  tmp /= detM;
//  intersect.set(idx, tmp);
//  return true;
//}
//
//bool
//brille::intersect_at(
//  const LQVec<double>& n_i,
//  const LQVec<double>& n_j,
//  const LQVec<double>& n_k,
//  LQVec<double>& intersect, const int idx
//){
//  double detM = brille::normals_matrix_determinant(n_i,n_j,n_k);
//  if (brille::approx::scalar(detM, 0.)) return false;
//  intersect.set(idx, 0*intersect.view(idx));
//  return true;
//}
