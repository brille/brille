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
using namespace brille::lattice;

LVec<double> BrillouinZone::get_ir_polyhedron_wedge_normals() const {
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

BrillouinZone::poly_t BrillouinZone::get_polyhedron() const {return _first;}
BrillouinZone::poly_t BrillouinZone::get_ir_polyhedron(const bool true_ir) const {
  // If the ir polyhedron fills the first BZ using the symmetry operations
  // of the pointgroup (plus time inversion, if included), then no ir_mirroring
  // is to be performed. In this case or if the 'true_ir' was requested, return
  // the computed ir_polyhedron
  if (no_ir_mirroring || !true_ir) return _irreducible;
  // Otherwise, it is only half of the true (non-convex)
  // irreducible polyhedron, so add the lack of inversion symmetry explicitly.
  return _irreducible + _irreducible.mirror();
}
bool BrillouinZone::check_ir_polyhedron(){
  profile_update("Start BrillouinZone::check_ir_polyhedron");
  this->check_if_mirroring_needed(); // move this to end of wedge_brute_force?
  auto ps = _outer.pointgroup_symmetry();
  if (time_reversal) ps = ps.add_space_inversion();

  auto volume_goal = _first.volume() / static_cast<double>(ps.size());
  if (!brille::approx_float::scalar(_irreducible.volume(), volume_goal, float_tolerance, float_tolerance, approx_tolerance)){ // setting tol=100 moves the relative difference allowed to ~2e-9
    debug_update("The current 'irreducible' polyhedron,\n",_irreducible.python_string(),"\nhas the wrong volume since ",
        _irreducible.volume()," != ",volume_goal);
    debug_update(std::abs(_irreducible.volume() - volume_goal) / (_irreducible.volume() + volume_goal));
    return false;
  }
  ps = ps.higher(1);
  for (size_t i=0; i < ps.size(); ++i) {
    auto rotated = _irreducible.apply(ps, static_cast<ind_t>(i));
    if (_irreducible.intersects(rotated, float_tolerance, approx_tolerance)){
      debug_update("_irreducible\n", _irreducible.python_string(),"\nintersects\n",rotated.python_string());
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
  return _first.vertices();
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_primitive_vertices() const {
  auto v = this->get_vertices();
  if (this->isprimitive()) v = transform_to_primitive(_outer, v);
  return v;
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_points() const {
  return _first.face_points();
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_primitive_points() const {
  auto p = this->get_points();
  if (this->isprimitive()) p = transform_to_primitive(this->_outer, p);
  return p;
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_normals() const {
  return _first.face_normals();
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_primitive_normals() const {
  auto n = this->get_normals();
  if (this->isprimitive()) n = transform_to_primitive(this->_outer, n);
  return n;
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_half_edges() const{
  return _first.half_edges();
}
BrillouinZone::poly_t::faces_t::faces_t BrillouinZone::get_faces_per_vertex() const {
  return _first.faces_per_vertex();
}
BrillouinZone::poly_t::faces_t::faces_t BrillouinZone::get_vertices_per_face() const {
  return _first.faces().faces();
}
// irreducible first Brillouin zone
BrillouinZone::poly_t::vertex_t BrillouinZone::get_ir_vertices() const {
  return _irreducible.vertices();
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_ir_primitive_vertices() const {
  auto v = this->get_ir_vertices();
  if (this->isprimitive()) v = transform_to_primitive(this->_outer, v);
  return v;
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_ir_points() const {
  return _irreducible.face_points();
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_ir_primitive_points() const {
  auto p = this->get_ir_points();
  if (this->isprimitive()) p = transform_to_primitive(this->_outer, p);
  return p;
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_ir_normals() const {
  return _irreducible.face_normals();
}
BrillouinZone::poly_t::vertex_t BrillouinZone::get_ir_primitive_normals() const {
  auto n = this->get_ir_normals();
  if (this->isprimitive()) n = transform_to_primitive(this->_outer, n);
  return n;
}
BrillouinZone::poly_t::faces_t::faces_t BrillouinZone::get_ir_faces_per_vertex() const {
  return _irreducible.faces_per_vertex();
}
BrillouinZone::poly_t::faces_t::faces_t BrillouinZone::get_ir_vertices_per_face() const {
  return _irreducible.faces().faces();
}

// irreducible reciprocal space
LVec<double> BrillouinZone::get_ir_wedge_normals() const {
  if (ir_wedge_normals.size(0))
    return LVec<double>(LengthUnit::inverse_angstrom, _outer, ir_wedge_normals);
  return {LengthUnit::inverse_angstrom, _outer, 0u};
}
//
LVec<double> BrillouinZone::get_primitive_ir_wedge_normals() const {
  if (ir_wedge_normals.size(0)){
    auto x = LVec<double>(LengthUnit::inverse_angstrom, _outer, ir_wedge_normals);
    if (this->isprimitive())
      x = transform_to_primitive(_outer, x);
    return x;
  }
  return {LengthUnit::inverse_angstrom, _outer, 0u};
}


void BrillouinZone::print() const {
  std::string msg = "BrillouinZone with ";
  msg += std::to_string(this->vertices_count()) + " vertices and ";
  msg += std::to_string(this->faces_count()) + " faces";
  std::cout << msg << std::endl;
}

void BrillouinZone::voro_search(const int extent, const bool divide_primitive){
  profile_update("Start BrillouinZone::voro_search with ",extent," extent");
  using namespace brille;
  LVec<int> primitive_tau(LengthUnit::inverse_angstrom, _inner, make_relative_neighbour_indices(extent));
  std::vector<ind_t> perm(primitive_tau.size(0));
  std::iota(perm.begin(), perm.end(), 0u); // {0u, 1u, 2u, ..., ntau-1}
  std::sort(perm.begin(), perm.end(), [&](ind_t a, ind_t b){
    return primitive_tau.norm(a) < primitive_tau.norm(b);
  });
  primitive_tau.permute(perm);
//  verbose_update("sorted primitive_tau\n",cat(1,primitive_tau,norm(primitive_tau)).to_string());

  auto tau = primitive_tau;
  if (!divide_primitive){
    tau = transform_from_primitive(_outer, primitive_tau);
  }
//  // the first Brillouin zone polyhedron will be expressed in absolute units
//  // in the xyz frame of the conventional reciprocal lattice
//  auto tau = transform_from_primitive(_outer, primitive_tau);

//  verbose_update("sorted tau\n",cat(1,tau,norm(tau)).to_string());
  auto [plane_a, plane_b, plane_c] = plane_points_from_normal(tau / norm(tau), tau / 2.0);
  // and then use the reciprocal lattice points to subdivide the cell until
  // only the first Brillouin zone is left:
  auto box = polyhedron::bounding_box(1.0 * tau);
//  verbose_update("bounding box\n", box.python_string());
//  throw std::runtime_error("time to halt is now!");
  _first = box.cut(plane_a, plane_b, plane_c, float_tolerance, approx_tolerance);
  if (divide_primitive){
    auto vertices = transform_from_primitive(_outer, _first.vertices());
    _first = polyhedron::Poly(vertices, _first.faces());
  }
//  verbose_update("cut box\n", _first.python_string());
  profile_update("  End BrillouinZone::voro_search with ",extent," extent");
}
