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
#ifndef BRILLE_BZ_CLASS_H_
#define BRILLE_BZ_CLASS_H_
/*! \file
    \author Greg Tucker
    \brief Defines a Brillouin zone class
*/
#include <omp.h>

#include <utility>
#include "neighbours.hpp"
#include "transform.hpp"
#include "polyhedron_flex.hpp"
#include "phonon.hpp"
#include "hdf_interface.hpp"
#include "array_l_.hpp"
#include "bz_config.hpp"
#include "approx_float.hpp"
#include "approx_config.hpp"

namespace brille {
/*! \brief An object to hold information about the first Brillouin zone of a
  Reciprocal lattice

  The BrillouinZone object is created from a Reciprocal lattice and if that
  lattice is not a primitive type the object contains the primitive lattice as
  well.
  By looking for intersections between three (primitive) reciprocal planes, the
  class object identifies the minimal subset of vertices that define the first
  Brillouin zone as well as keeping track of the planes which contributed to
  each vertex and the plane reciprocal vectors themselves.
*/
class BrillouinZone {
public:
  using lattice_t = lattice::Lattice<double>;
  using poly_t = polyhedron::Poly<double, lattice::LVec>;

protected:
  lattice_t _inner;    //!< The primitive reciprocal lattice
  lattice_t _outer;    //!< The lattice passed in at construction
  poly_t _first;       //!< The first Brillouin zone polyhedron
  poly_t _irreducible; //!< The irreducible Brillouin zone polyhedron
  bArray<double> ir_wedge_normals; //!< The normals of the irreducible
                                   //!< reciprocal space wedge planes.
  bool time_reversal; //!< A flag to indicate if time reversal symmetry should
                      //!< be included with pointgroup operations
  bool has_inversion; //!< A computed flag indicating if the pointgroup has
                      //!< space inversion symmetry or if time reversal symmetry
                      //!< has been requested
  bool is_primitive; //!< A computed flag indicating if the primitive version of
                     //!< a conventional lattice is in use
  bool no_ir_mirroring;
  int approx_tolerance;
  double float_tolerance;

public:
  /*! \brief Construct a Brillouin zone object from all properties
   *
   * @param lat The primitive reciprocal lattice
   * @param olat The lattice (normally) passed at construction
   * @param poly The first Brillouin zone polyhedron
   * @param irp An irreducible Brillouin zone polyhedron
   * @param irn The irreducible Brillouin zone wedge normals
   * @param tr Whether time reversal symmetry has been added
   * @param hi Whether the point group symmetry has space inversion (or time
   * reversal was added)
   * @param ip Whether the primitive version of a conventional lattice is used
   * @param nim Whether the irreducible Brillouin zone requires mirroring to get
   * the correct result
   */
  BrillouinZone(lattice_t in_lat, lattice_t out_lat, poly_t f_poly,
                poly_t i_poly, const bArray<double> &irn, bool tr, bool hi,
                bool ip, bool nim, int tol, double f_tol)
      : _inner(std::move(in_lat)), _outer(std::move(out_lat)),
        _first(std::move(f_poly)), _irreducible(std::move(i_poly)),
        ir_wedge_normals(irn), time_reversal(tr), has_inversion(hi),
        is_primitive(ip), no_ir_mirroring(nim), approx_tolerance(tol),
        float_tolerance(f_tol) {}
  /*! \brief Construct a Brillouin zone object

  \param lat A `Reciprocal` lattice
  \param to_prim A flag to indicate if the primtive lattice should be used
  \param extent An integer to control how-far the vertex-finding algorithm
  should search in τ-index. The default indicates that (̄1̄1̄1), (̄1̄10), (̄1̄11),
  (̄10̄1), ..., (111) are included. \param tr controls whether time reversal
  symmetry should be added if absent \param wedge_search controls whether and
  attempt is made to find the irreducible reciprocal space wedge and, therefore,
  the irreducible Brillouin zone.
  */
  explicit BrillouinZone(const lattice_t &lat,
                         const BrillouinZoneConfig cfg = BrillouinZoneConfig())
      : _inner(cfg.primitive() ? lat.primitive() : lat), _outer(lat),
        time_reversal(cfg.time_reversal()) {
    is_primitive = !(_inner.is_same(_outer));
    has_inversion = time_reversal || lat.has_space_inversion();
    no_ir_mirroring = true;
    // the approximate configuration can not be a default function parameter
    // if we want to pick-up runtime changes in the namespace object
    approx_tolerance = approx_float::config.digit();
    float_tolerance = approx_float::config.reciprocal<double>();
    this->construct(cfg.divide_extent() - 1, cfg.divide_primitive(),
                    cfg.wedge_search());
  }
  BrillouinZone(const lattice_t &lat, const BrillouinZoneConfig cfg,
                const approx_float::Config ac)
      : _inner(cfg.primitive() ? lat.primitive() : lat), _outer(lat),
        time_reversal(cfg.time_reversal()), approx_tolerance(ac.digit()),
        float_tolerance(ac.reciprocal<double>()) {
    is_primitive = !(_inner.is_same(_outer));
    has_inversion = time_reversal || lat.has_space_inversion();
    no_ir_mirroring = true;
    this->construct(cfg.divide_extent() - 1, cfg.divide_primitive(),
                    cfg.wedge_search());
  }
  BrillouinZone(const lattice_t &lat, const approx_float::Config ac,
                const BrillouinZoneConfig cfg = BrillouinZoneConfig())
      : _inner(cfg.primitive() ? lat.primitive() : lat), _outer(lat),
        time_reversal(cfg.time_reversal()), approx_tolerance(ac.digit()),
        float_tolerance(ac.reciprocal<double>()) {
    is_primitive = !(_inner.is_same(_outer));
    has_inversion = time_reversal || lat.has_space_inversion();
    no_ir_mirroring = true;
    this->construct(cfg.divide_extent() - 1, cfg.divide_primitive(),
                    cfg.wedge_search());
  }

  void construct(int test_extent, bool divide_primitive, bool wedge_search) {
    profile_update("Start of BrillouinZone construction");
    double old_volume = -1.0, new_volume = 0.0;
    // initial test_extent based on spacegroup or pointgroup?
    do {
      old_volume = new_volume;
      this->voro_search(++test_extent, divide_primitive);
      new_volume = _first.volume();
    } while (!approx_float::scalar(new_volume, old_volume, float_tolerance,
                                   float_tolerance, approx_tolerance));

    if (brille::approx_float::scalar(new_volume, 0., float_tolerance,
                                     float_tolerance, approx_tolerance)) {
      throw std::runtime_error(
          "failed to produce a non-null first Brillouin zone.");
    }
    // in case we've been asked to perform a wedge search for, e.g., P1 or P-1,
    // set the irreducible wedge now as the search will do nothing.
    _irreducible = _first;
    if (wedge_search) {
      debug_update("Starting wedge search from ", _first.python_string());
      bool success{false};
      if (_outer.is_triclinic()) {
        success = !has_inversion || this->wedge_triclinic();
      } else {
        success =
            this->wedge_brute_force()
            // no special 2-fold or mirror handling
            || this->wedge_brute_force(false, false)
            // no special 2-fold handling (but special mirror handling)
            || this->wedge_brute_force(false, true)
            // no special mirror handling (maybe not useful)
            || this->wedge_brute_force(true, false)
            // handle non order(2) operations in decreasing order
            || this->wedge_brute_force(true, true, false)
            || this->wedge_brute_force(true, false, true, false)
            // do not handle 2-folds or mirrors first; start with highest-order (non-inverting) operation
            || this->wedge_brute_force(false, false, false)
            // also do not sort list of equivalent remaining points when identifying the order>2 wedges
            || this->wedge_brute_force(false, false, false, false);
        // other combinations of special_2_folds, special_mirrors,
        // and sort_by_length are possible but not necessarily useful.
      }
      if (!success) {
//        info_update("First Brillouin zone\n", _first.python_string(),
//                    "\nand 'irreducible' Brillouin zone\n",
//                    _irreducible.python_string());
        std::string msg = "Failed to find an irreducible Brillouin zone.\n";
        msg += "First Brillouin zone\n" + _first.python_string();
        msg += "Last tried candidate\n" + _irreducible.python_string();
        msg += " Consider increasing tolerances from ";
        msg += my_to_string(float_tolerance);
        msg += " and " + std::to_string(approx_tolerance);
        throw std::runtime_error(msg);
      }
    }
    profile_update("  End of BrillouinZone construction");
  }
  /*! \brief Decide if the held 'irreducible' Brillouin zone polyhedron must be
  mirrored

  The wedge searching algorithm may result in a polyhedron which is only half
  the volume of the true irreducible Brillouin zone. This is most likely to
  occur if the spacegroup has inversion symmetry and/or the Brillouin zone is to
  include time reversal symmetry.
  In the case where this is true we can regain *an* irreducible Brillouin zone
  by finding the union of the polyhedron and its inverse, however, this union is
  unlikely to be convex which will cause problems for other algorthims.

  The solution employed here finds the union of the polyhedron and each of the
  pointgroup operated mirrored polyhedra then selects the union with the fewest
  unique vertices as most-likely-to-be convex. If the selected union polyhedron
  has the same volume as the irreducible Brillouin zone it replaces the
  polyhedron stored in this object.
  */
  void check_if_mirroring_needed() {
    no_ir_mirroring = true;
    if (!has_inversion) {
      auto ps = _outer.pointgroup_symmetry();
      if (time_reversal)
        ps = ps.add_space_inversion();
      auto goal = _first.volume() / static_cast<double>(ps.size());
      auto found = _irreducible.volume();
      if (approx_float::scalar(goal, 2.0 * found, float_tolerance,
                               float_tolerance, approx_tolerance)) {
        debug_update("Apply mirroring strategy since ", goal, " ~= 2* ", found);
        /*The current Polyhedron at this->ir_polyhedron has half the anticipated
        volume. We can 'fix' this the easy way by mirroring the Polyhedron and
        gluing it onto the current version; the resultant polyhedron will come
        to a point at (0,0,0) and will not be convex. Instead, try to be clever
        and look for a convex combination of the found polyhedron and its
        mirror with one of the pointgroup operations applied:*/
        std::vector<poly_t> valid;
        auto mirrored = _irreducible.mirror();
        for (ind_t i = 0; i < ps.size(); ++i) {
          auto u = _irreducible.combine(mirrored.apply(ps, i), float_tolerance, approx_tolerance);
          verbose_update("Mirror union with volume ",u.volume(),":\n", u.python_string());
          // Polyhedron::combine needs to be fixed. As of now it does not look
          // for extraneous faces (like internal faces which are coplanar and
          // pointing in opposite directions) or coplanar external faces which
          // share an edge. Until this is fixed cross our fingers and hope that we
          // created a convex polyhedron such that the convex hull of its points
          // gives the same polyhedron back:
          auto c = poly_t(u.vertices(), float_tolerance, approx_tolerance);
          verbose_update("Produces a convex hull with volume ", c.volume(), ":\n", c.python_string());
          if (approx_float::scalar(c.volume(), goal, float_tolerance, float_tolerance, approx_tolerance)){
            valid.push_back(c);
          }
        }
        // any polyhedra in `valid` should be equally-good, so take the first one
        if (!valid.empty()) {
          _irreducible = valid.front();
          found = goal;
        }
      }
      // if found == goal no mirroring is required.
      // if found == goal/2, mirroring is required.
      no_ir_mirroring = approx_float::scalar(goal, found, float_tolerance,
                                             float_tolerance, approx_tolerance);
    }
  }
  /*! \brief Determine the validity of the stored irreducible Brillouin zone
  polyhedron

  An irreducible Brillouin zone should have volume equal to that of the first
  Brillouin zone divided by the number of pointgroup symmetry operations.
  Furthermore, each pointgroup operation acting on an irreducible Brillouin zone
  polyhedron should produce a new polyhedron will null intersection.

  \return Whether both requirements are fulfilled.
  */
  bool check_ir_polyhedron();
  //  bool wedge_explicit();
  //! Returns the lattice passed in at construction
  [[nodiscard]] lattice_t get_lattice() const { return _outer; };
  //! Returns the lattice actually used to find the Brillouin zone vertices,
  //! which may be a primitive lattice depending on the flag at creation
  [[nodiscard]] lattice_t get_primitive_lattice() const { return _inner; };
  //! Returns the number of vertices defining the first Brillouin zone
  [[nodiscard]] ind_t vertices_count() const {
    return this->get_vertices().size(0);
  };
  //! Returns the number of reciprocal lattice points defining the first
  //! Brillouin zone
  [[nodiscard]] ind_t faces_count() const {
    return this->get_points().size(0);
  };
  // irreducible reciprocal space wedge
  void set_ir_wedge_normals(const lattice::LVec<double> &x) {
    bool already_same = _outer.is_same(x.lattice());
    lattice::LVec<double> xp(x.type(), _outer);
    PrimitiveTransform PT(_outer.bravais());
    bool transform_needed = (PT.does_anything() && _inner.is_same(x.lattice()));
    if (!(already_same || transform_needed))
      throw std::runtime_error(
          "ir_wedge_normals must be in the standard or primitive lattice used "
          "to define the BrillouinZone object");
    if (transform_needed)
      xp = transform_from_primitive(_outer, x);
    const lattice::LVec<double> &xref = transform_needed ? xp : x;
    ir_wedge_normals = xref.hkl();
  }
  [[nodiscard]] lattice::LVec<double> get_ir_wedge_normals() const;
  [[nodiscard]] lattice::LVec<double> get_primitive_ir_wedge_normals() const;
  //  /*! \brief Set the first Brillouin zone polyhedron
  //
  //  Set the first Brillouin zone polyhedron from its vertices, central facet
  //  plane points, and either the three intersecting planes which gave each
  //  vertex or all intersecting planes which give each vertex as well as all
  //  vertices which form a corner of each facet plane polygon. \param v All
  //  vertices in the polyhedron \param p All (τ/2) plane points which define
  //  the facets of the polyhedron \param args A variable argument input passed
  //  directly to the various
  //              `Polyhedron` class constructors.
  //
  //  The (ordered) optional inputs in `args` are:
  //  1. `fpv` A `std::vector<std::vector<int>>` with *all* facet-plane indices
  //  for each vertex
  //  2. `vpf` A `std::vector<std::vector<int>>` with *all* vertex indices for
  //  each facet-plane
  //  */
  //  template<typename... A>
  //  void
  //  set_polyhedron(const lattice::LVec<double>& v, const
  //  lattice::LVec<double>& p, A... args){
  //    bool both_same = v.get_lattice().is_same(p.get_lattice());
  //    if (!both_same)
  //      throw std::runtime_error("The vertices, and points of a polyhedron
  //      must all be in the same cooridinate system");
  //    bool is_outer = this->_outer.is_same(v.get_lattice());
  //    bool is_inner = this->lattice.is_same(v.get_lattice());
  //    lattice::LVec<double> vp(this->_outer), pp(this->_outer);
  //    PrimitiveTransform PT(this->_outer.get_bravais_type());
  //    is_inner &= PT.does_anything();
  //    if (!(is_outer || is_inner))
  //      throw std::runtime_error("The polyhedron must be described in the
  //      conventional or primitive lattice used to define the BrillouinZone
  //      object");
  //    if (is_inner){
  //      vp = transform_from_primitive(this->_outer, v);
  //      pp = transform_from_primitive(this->_outer, p);
  //    }
  //    const lattice::LVec<double> & vref = is_inner ? vp : v;
  //    const lattice::LVec<double> & pref = is_inner ? pp : p;
  //    this->polyhedron = Polyhedron(vref.get_xyz(), pref.get_xyz(), args...);
  //  }
  /*! \brief Set the irreducible Brillouin zone polyhedron

  Set the irreducible first Brillouin zone polyhedron from its vertices, facet
  plane points, facet plane normals, and, optionally, all intersecting planes
  which give each vertex as well as all vertices which form a corner of each
  facet plane polygon.
  \param v All vertices in the polyhedron
  \param p All (τ/2) plane points which define the facets of the polyhedron
  \param n All facet normals of the polyhedron
  \param args A variable argument input passed directly to the various
              `Polyhedron` class constructors.

  The (ordered) optional inputs in `args` are:
  1. `fpv` A `std::vector<std::vector<int>>` with *all* facet-plane indices for
  each vertex
  2. `vpf` A `std::vector<std::vector<int>>` with *all* vertex indices for each
  facet-plane
  */
  template <typename... A>
  void set_ir_polyhedron(const lattice::LVec<double> &v,
                         const lattice::LVec<double> &p,
                         const lattice::LVec<double> &n, A... args) {
    bool all_same =
        v.lattice().is_same(p.lattice()) && p.lattice().is_same(n.lattice());
    if (!all_same)
      throw std::runtime_error(
          "The vertices, points, and normals of a polyhedron must all be in "
          "the same cooridinate system");
    bool is_outer = _outer.is_same(v.lattice());
    bool is_inner = _inner.is_same(v.lattice());
    lattice::LVec<double> vp(v.type(), _outer), pp(p.type(), _outer),
        np(n.type(), _outer);
    PrimitiveTransform PT(_outer.bravais());
    is_inner &= PT.does_anything();
    if (!(is_outer || is_inner))
      throw std::runtime_error(
          "The polyhedron must be described in the conventional or primitive "
          "lattice used to define the BrillouinZone object");
    if (is_inner) {
      vp = transform_from_primitive(_outer, v);
      pp = transform_from_primitive(_outer, p);
      np = transform_from_primitive(_outer, n);
    }
    const auto &vr{is_inner ? vp : v};
    const auto &pr{is_inner ? pp : p};
    const auto &nr{is_inner ? np : n};
    _irreducible = poly_t(vr, pr, nr, args...);
  }
  //  /*! \brief Set the irreducible first Brillouin zone polyhedron from its
  //  vertices
  //
  //  Since no facet information is provided via this method, a convex hull of
  //  the provided points is calculated and used as the irreducible polyhedron.
  //  Additionally the volume and symmetry of the calculated polyhedron is
  //  checked against the first Brillouin zone using the pointgroup symmetry
  //  operations.
  //
  //  \param v All vertices in the irreducible polyhedron
  //  \return A bool indicating if the found convex hull polyhedron is the
  //          expected volume and has the expected symmetry.
  //  */
  //  bool set_ir_vertices(const lattice::LVec<double>& v){
  //    bool is_outer = this->_outer.is_same(v.lattice());
  //    bool is_inner = this->_inner.is_same(v.lattice());
  //    lattice::LVec<double> vp(this->_outer);
  //    PrimitiveTransform PT(this->_outer.get_bravais_type());
  //    is_inner &= PT.does_anything();
  //    if (!(is_outer || is_inner))
  //      throw std::runtime_error("The polyhedron must be described in the
  //      conventional or primitive lattice used to define the BrillouinZone
  //      object");
  //    if (is_inner)
  //      vp = transform_from_primitive(this->_outer, v);
  //    const lattice::LVec<double> & vref = is_inner ? vp : v;
  //    debug_update("Generate a convex polyhedron from\n", vref.to_string());
  //    this->ir_polyhedron = Polyhedron(vref.get_xyz());
  //    debug_update("Generated a ", this->ir_polyhedron.string_repr()," from
  //    the specified vertices");
  //    // check that the polyhedron we've specified has the desired properties
  //    if (this->check_ir_polyhedron()){
  //      // set the wedge normals as well
  //      this->set_ir_wedge_normals(this->get_ir_polyhedron_wedge_normals());
  //      return true;
  //    }
  //    return false;
  //  }

  // TODO Most of these accessors can be removed and replaced with the first or
  // irreducible polyhedra properties
  //! Returns the first Brillouin zone polyhedron
  [[nodiscard]] poly_t get_polyhedron() const;
  //! Returns the vertices of the first Brillouin zone polyhedron expressed as
  //! conventional unit cell vectors
  [[nodiscard]] poly_t::vertex_t get_vertices() const;
  //! Returns the on-plane points (τ/2) of the first Brillouin zone polyhedron
  //! expressed as conventional unit cell vectors
  [[nodiscard]] poly_t::vertex_t get_points() const;
  //! Returns the normals ( ̂τ ) of the first Brillouin zone polyhedron expressed
  //! as conventional unit cell vectors
  [[nodiscard]] poly_t::vertex_t get_normals() const;
  //! Returns the half-edge points of the first Brillouin zone polyhedron
  //! expressed as conventional unit cell vectors
  [[nodiscard]] poly_t::vertex_t get_half_edges() const;
  //! Returns the vertices of the first Brillouin zone polyhedron expressed as
  //! primitive unit cell vectors
  [[nodiscard]] poly_t::vertex_t get_primitive_vertices() const;
  //! Returns the on-plane points (τ/2) of the first Brillouin zone polyhedron
  //! expressed as primitive unit cell vectors
  [[nodiscard]] poly_t::vertex_t get_primitive_points() const;
  //! Returns the normals ( ̂τ ) of the first Brillouin zone polyhedron expressed
  //! as primitive unit cell vectors
  [[nodiscard]] poly_t::vertex_t get_primitive_normals() const;
  //! Returns the `points` and `normals` indices for each vertex of the first
  //! Brillouin zone polyhedron
  [[nodiscard]] poly_t::faces_t::faces_t get_faces_per_vertex() const;
  //! Returns the vertex indices for each facet of the first Brillouin zone
  //! polyhedron
  [[nodiscard]] poly_t::faces_t::faces_t get_vertices_per_face() const;
  //! Returns the irreducible first Brillouin zone polyhedron
  [[nodiscard]] poly_t get_ir_polyhedron(bool true_ir = true) const;
  //! Returns the vertices of the irreducible first Brillouin zone polyhedron
  //! expressed as conventional unit cell vectors
  [[nodiscard]] poly_t::vertex_t get_ir_vertices() const;
  //! Returns the on-plane points of the irreducible first Brillouin zone
  //! polyhedron expressed as conventional unit cell vectors
  [[nodiscard]] poly_t::vertex_t get_ir_points() const;
  //! Returns the normals of the irreducible first Brillouin zone polyhedron
  //! expressed as conventional unit cell vectors
  [[nodiscard]] poly_t::vertex_t get_ir_normals() const;
  //! Returns the vertices of the irreducible first Brillouin zone polyhedron
  //! expressed as primitive unit cell vectors
  [[nodiscard]] poly_t::vertex_t get_ir_primitive_vertices() const;
  //! Returns the on-plane points of the irreducible first Brillouin zone
  //! polyhedron expressed as primitive unit cell vectors
  [[nodiscard]] poly_t::vertex_t get_ir_primitive_points() const;
  //! Returns the normals of the irreducible first Brillouin zone polyhedron
  //! expressed as primitive unit cell vectors
  [[nodiscard]] poly_t::vertex_t get_ir_primitive_normals() const;
  //! Returns the `points` and `normals` indices for each vertex of the
  //! irreducible first Brillouin zone polyhedron
  [[nodiscard]] poly_t::faces_t::faces_t get_ir_faces_per_vertex() const;
  //! Returns the vertex indices for each facet of the irreducible first
  //! Brillouin zone polyhedron
  [[nodiscard]] poly_t::faces_t::faces_t get_ir_vertices_per_face() const;
  //! Print a representation of the object to std::cout
  void print() const;
  //  /*! \brief Attempt to find an irreducible section of reciprocal space
  //
  //  Using the pointgroup symmetry opperations of the conventional unit cell,
  //  this method defines *an* irreducible section of reciprocal space.
  //  The irreducible section bounds a solid angle N and consequently the full
  //  reciprocal space is 4π/N larger than the irreducible wedge.
  //  As the selected irreducible section is related by symmetry to the rest of
  //  reciprocal space there are (at least) 4π/N - 1 equivalent other choices
  //  for the irreducible wedge.
  //
  //  This method finds the stationary axes of the 'proper' rotation matrix of
  //  each pointgroup operation (if an operation, R, is a rotoinversion, the
  //  'proper' rotation is ̄1R), collectively \f$\hat{z}_i\f$, and two
  //  right-handed perpendicular axes \f$\hat{x}_i\f$ and \f$\hat{y}_i\f$. It
  //  then attempts to find a combination of subsets of \f$\hat{x}_i\f$ and
  //  \f$\hat{y}_i\f$ which are the normals of a convex 'wedge' (a flat-sided
  //  cone) with the expected subtended solid angle.
  //
  //  \param prefer_basis_vectors Whether the reciprocal lattice basis vectors
  //                              should be used instead of the pointgroup
  //                              operator perpendicular axes when determining
  //                              \f$\hat{x}\f$
  //  \param parallel_ok Whether a newly found \f$\hat{x}_j\f$ is allowed to be
  //                     parallel to the \f$\hat{z}_i\f$ used in its
  //                     determination
  //
  //  \deprecated This method is less successful than wedge_brute_force
  //  \see wedge_brute_force
  //  */
  //  void wedge_search(const bool prefer_basis_vectors=true, const bool
  //  parallel_ok=true);
  /*! \brief Attempt to find an irreducible section of reciprocal space

  This method collects the vertices, facet centres, and mid-edge points of the
  first Brillouin zone polyhedron, which are all points of high symmetry, and
  attempts to find a subset of these which are the vertices of an irreducible
  section of the first Brillouin zone.

  For each of the pointgroup operations in turn other than 1 and ̄1, the
  high-symmetry points are collected into subsets of points which successive
  applications of the operation map onto each other. The subsets are ordered
  with an arbitrary first point and each successive element of the set being the
  next application of the operation, including identification of missing subset
  members. With all subsets found for a pointgroup operation the most-complete
  highest-member-count subgroup is chosen to cull the list of high-symmetry
  points. If the pointgroup operation is 2 or ̄2 then two points in the chosen
  subset (and the Γ point) define a plane -- this plane divides space and the
  dot product with its normal is used to determine which high-symmetry points
  are removed (if their dot product is negative).
  If the pointgroup operation is higher-order then two successive points in the
  chosen subset, each combined with a point along the operations characteristic
  axis and Γ, define two planes and high-symmetry points 'behind' either are
  removed.

  The high-symmetry points remaining after all pointgroup operations have been
  used in this way *might* be the vertices of an irreducible section of the
  first Brillouin zone. The method checks this by forming a Polyhedron from
  their convex hull, verifying that it has the expected volume relative to the
  first Brillouin zone, and then ensuring that every pointgroup operation acting
  on the test polyhedron produces a new polyhedron which does not intersect with
  the test polyhedron.

  \param special_2_folds If true any 2 pointgroup operations are used first
                         to cull the high-symmetry points along with preselected
                         dividing planes which ensure the an irreducible
                         Brillouin zone is found for high-symmetry cubic
                         lattices
  \param special_mirrors If true any ̄2 pointgroup operations are used directly
                         after the special 2-fold axes with their characteristic
                         axis choses as the dividing plane normal
  \param sort_by_length  If true the characteristic axis length is used to sort
                         the pointgroup operations so that, e.g., [100] is used
                         before [111] in a cubic system; otherwise the
                         operations are sorted by decreasing isometry:
                         6, 4, 3, 2, -2, -3, -4, -6
  \param sort_one_sym    If true the equivalent subsets of high-symmetry points
                         for a given pointgroup operation are sorted by how many
                         of the expected number were found to still remain in
                         the list of culled high-symmetry points; otherwise the
                         order of subsets depends on the order of the remaining
                         high-symmetry points.
  \return true if the found convex polyhedron has all of the expected properties
          of an irreducible Brillouin zone polyhedron.
  */
  bool wedge_brute_force(bool special_2_folds = true,
                         bool special_mirrors = true,
                         bool sort_by_length = true, bool sort_one_sym = true);
  bool wedge_triclinic();
  //  /*!
  //  With the first Brillouin zone and *an* irreducible section of reciprocal
  //  space already identified, this method finds all intersections of
  //  combinations of three planes (nBz,nIr) = {(0,3), (1,2), (2,1), and (3,0)}
  //  where nBz are the number of intersecting first Brillouin zone facet planes
  //  and nIr are the number of intersecting irreducible reciprocal space wedge
  //  planes. Once all combinations of intersection points are found, those
  //  within the first Brillouin zone *and* the irreducible reciprocal space
  //  wedge define the irreducible first Brillouin zone.
  //  @note the (3,0) intersection points *are* the first Brillouin zone
  //  vertices and the (0,3) intersection points are all identically ⃗0.
  //  */
  //  void irreducible_vertex_search();
  /*! \brief Construct the first Brillouin zone Polyhedron

  \param scale A scale factor for the starting box

  One common construction of the first Brillouin zone is the region of
  reciprocal space closer to one primitive lattice point than any other. In this
  construction all planes dividing space between any two primitive lattice
  points form part of the surface of *a* Brillouin zone, but not necessarily
  the *first* Brillouin zone. This method uses successive division of a
  polyhedron to identify the first Brillouin zone for a primitive lattice.

  The starting polyhedron is a box bounding all integer primitive lattice points
  (`i`,`j`,`k`) where each `i`, `j`, and `k` are ∈ `-scale`,…,-1,0,1,…,`scale`.
  Then, the (`i`,`j`,`k`) are sorted in increasing distance from (`0`,`0`,`0`)
  and used successively to divide the polyhedron by the plane passing through
  (`i`,`j`,`k`)/2 with normal (`i`,`j`,`k`), keeping the subpolyhedron which
  contains (`0`,`0`,`0`). This results in the first Brillouin zone centred on
  the origin of reciprocal space.

  The algorithm used by this method has been inspired by the Voronoi cell
  software library [Voro++]([http://math.lbl.gov/voro++/), especially
  [voro::voronoicell_base](http://math.lbl.gov/voro++/doc/refman/classvoro_1_1voronoicell__base.html).
  */
  void voro_search(int scale = 1, bool divide_primitive = false);
  //! Return true if the lattice used to find the vertices is not the same as
  //! the one passed at construction, and is therefore primitive
  [[nodiscard]] bool isprimitive() const { return is_primitive; };

  /*! \brief Determine whether points are inside of the first Brillouin zone
    @param p A reference to a lattice::LVec list of Q points to be checked
    @returns A `std::vector<bool>`` with each element indicating if the
             associated Q point is inside of the Brillouin zone.
  */
  template <class T>
  [[nodiscard]] std::vector<bool> isinside(const lattice::LVec<T> &p) const {
    bool inner = _inner.is_same(p.lattice());
    if (!(inner || _outer.is_same(p.lattice())))
      throw std::runtime_error(
          "Q points must be in the standard or primitive lattice");
    lattice::LVec<T> pp(p.type(), _outer);
    PrimitiveTransform PT(_outer.bravais());
    inner &= PT.does_anything();
    if (inner)
      pp = transform_from_primitive(_outer, p);
    return _first.contains(inner ? pp : p, float_tolerance, approx_tolerance);
  }
  /*! \brief Determine whither points are inside the irreducible reciprocal
  space wedge

  A point is within the irreducible wedge if its dot product with all wedge
  normals is greater or equal to zero *or* if all dot products are negative or
  equal to zero since the irreducible wedge is a flat-sided cone.
  This method performs the dot product for the vector from the origin to each
  provided point with all wedge normals to determine if each point is inside the
  irreducible wedge.

  @param p A reference to a lattice::LVec list of Q points to be checked
  @param pos Restrict the comparison to `dot(n,p)≥0`; useful when determining
             the normals which define the irreducible wedge.
  @returns An `std::vector<bool>`` with each element indicating if the
           associated Q point is inside our irreducible reciprocal space.
  */
  template <class T>
  std::vector<bool> isinside_wedge(const lattice::LVec<T> &p,
                                   const bool pos = false) const {
    bool isouter = this->_outer.is_same(p.lattice());
    bool isinner = this->_inner.is_same(p.lattice());
    if (!(isouter || isinner)) {
      std::string msg = "Q points provided to BrillouinZone::isinside_wedge ";
      msg += "must be in the standard or primitive lattice ";
      msg += "used to define the BrillouinZone object";
      throw std::runtime_error(msg);
    }
    std::vector<bool> out(p.size(0), true);
    auto normals =
        isouter ? get_ir_wedge_normals() : get_primitive_ir_wedge_normals();
    if (normals.size(
            0)) { // with no normals *all* points are "inside" the wedge
      // If a pointgroup has inversion symmetry then for every point, p, there
      // is an equivalent point, -p. This indicates that a point p is already in
      // the irreducible wedge if it has n̂ᵢ⋅p ≥ 0 for all irredudible-bounding-
      // plane normals, ̂nᵢ, *or* the opposite -- n̂ᵢ⋅ p ≤ 0 for all ̂nᵢ.
      // Since we are interested in enforcing a smallest-possible irreducible
      // Brillouin zone, we want to exclude the all n̂ᵢ⋅ p ≤ 0 half of the
      // reciprocal wedge precisely because they are equivalent.
      // It is only in the case where a pointgroup *does not have*
      // space-inversion symmetry and time-reversal symmetry is to be excluded
      // that we can allow n̂ᵢ⋅ p ≤ 0 to be a valid in-wedge solution. The Array
      // all method has a special switched execution path for checking whether
      // all values are ≤ or ≥ a value simultaneously

      // when constructing the irreducible Brillouin zone we need to only
      // consider the ≥0 case so that we end up with a convex polyhedron. The
      // ir_polyhedron accessor method mirrors the half-polyhedron in this case,
      // so when identifying whether a point is inside of the irreducible
      // Brillouin zone we must allow for the ≤0 case as well.
      cmp c = pos || no_ir_mirroring ? cmp::ge : cmp::le_ge;
#pragma omp parallel for default(none) shared(out, normals, p, c)              \
    schedule(dynamic)
      for (long long i = 0; i < p.size(0); ++i) // separately changed to ind_t
        //TODO Add tolerance here?
        out[i] = dot(normals, p.view(i)).all(c, 0.);
    }
    return out;
  }
  /*! \brief Find q and τ such that Q=q+τ and τ is a reciprocal lattice vector
    \param[in] Q A reference to lattice::LVec list of Q points
    \param[out] q The reduced reciprocal lattice vectors
    \param[out] tau The reciprocal lattice zone centres
    \param threads The number of OpenMP threads to use, if less than one the
                   number returned by `omp_get_max_threads()` is used instead.
    \return true
  */
  bool moveinto(const lattice::LVec<double> &Q, lattice::LVec<double> &q,
                lattice::LVec<int> &tau, int threads=0) const;
  /*! \brief Find q, τ, and R∈G such that Q = Rᵀq + τ, where τ is a reciprocal
             lattice vector and R is a pointgroup symmetry operation of the
             conventional unit cell pointgroup, G.
    @param [in] Q a refernce to a lattice::LVec list of Q points
    @param [out] q The irreducible reduced reciprocal lattice vectors
    @param [out] tau The conventional reciprocal lattice zone centres
    @param [out] Ridx The pointgroup operation index for R
    @param [out] invRidx The pointgroup operation index for R⁻¹
    @param [in] threads An optional number of OpenMP threads to use
    @return the success status
    @note `Ridx` and `invRidx` index the `PointSymmetry` object accessible via
    `BrillouinZone::get_pointgroup_symmetry()`;
  */
  bool ir_moveinto(const lattice::LVec<double> &Q, lattice::LVec<double> &q,
                   lattice::LVec<int> &tau, std::vector<size_t> &Ridx,
                   std::vector<size_t> &invRidx, int threads = 0) const;
  /*! \brief Find q and R∈G such that Q = Rᵀq such that q ∈ the irreducible
  wedge

  \param [in] Q a refernce to a lattice::LVec list of Q points
  \param [out] q The irreducible reduced reciprocal lattice vectors
  \param [out] R The pointgroup operation matrix for R
  \param [in] threads An optional number of OpenMP threads to use
  \return the success status
  */
  bool ir_moveinto_wedge(const lattice::LVec<double> &Q,
                         lattice::LVec<double> &q, std::vector<size_t> &R,
                         int threads = 0) const;

  //! \brief Get the PointSymmetry object used by this BrillouinZone object
  //! internally
  [[nodiscard]] PointSymmetry get_pointgroup_symmetry() const {
    auto pgs = _outer.pointgroup_symmetry();
    return time_reversal ? pgs.add_space_inversion() : pgs;
  }
  //! \brief Accessor for whether the BrillouinZone was constructed with
  //! additional time reversal symmetry
  [[nodiscard]] int add_time_reversal() const { return time_reversal ? 1 : 0; }

private:
//  void _moveinto_prim(const lattice::LVec<double> &Q, lattice::LVec<double> &q,
//                      lattice::LVec<int> &tau, const lattice::LVec<double> &a,
//                      const lattice::LVec<double> &b,
//                      const lattice::LVec<double> &c, int threads) const;
  template <class T>
  [[nodiscard]] bool _inside_wedge_outer(const lattice::LVec<T> &p,
                                         const bool pos = false) const {
    brille::cmp c =
        pos || no_ir_mirroring ? brille::cmp::ge : brille::cmp::le_ge;
    auto normals = get_ir_wedge_normals();
    return normals.size(0) == 0 || dot(normals, p).all(c, T(0), float_tolerance, approx_tolerance);
  }

  void shrink_and_prune_outside(const size_t cnt, lattice::LVec<double> &vrt,
                                bArray<int> &ijk) const {
    verbose_update("shrinking to ", cnt);
    if (vrt.size(0) && ijk.size(0)) {
      vrt.resize(cnt);
      ijk.resize(cnt);
      if (cnt) { // isinside has a problem with vrt.size()==0
        auto isin = this->isinside(vrt);
        verbose_update("and retaining ",
                       std::count(isin.begin(), isin.end(), true),
                       " inside vertices");
        vrt = vrt.extract(isin);
        ijk = ijk.extract(isin);
      }
    }
  }

  [[nodiscard]] lattice::LVec<double> get_ir_polyhedron_wedge_normals() const;

public:
#ifdef USE_HIGHFIVE
  // FIXME Update these to reflect class contents
  template <class HF>
  std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, bool>
  to_hdf(HF &obj, const std::string &entry) const {
    auto group = overwrite_group(obj, entry);
    bool ok{true};
    ok &= _inner.to_hdf(group, "inner");
    ok &= _outer.to_hdf(group, "outer");
    ok &= _first.to_hdf(group, "first");
    ok &= _irreducible.to_hdf(group, "irreducible");
    ok &= ir_wedge_normals.to_hdf(group, "ir_wedge_normals");
    group.createAttribute("time_reversal", time_reversal ? 1 : 0);
    group.createAttribute("has_inversion", has_inversion ? 1 : 0);
    group.createAttribute("is_primitive", is_primitive ? 1 : 0);
    group.createAttribute("no_ir_mirroring", no_ir_mirroring ? 1 : 0);
    group.createAttribute("approx_tolerance", approx_tolerance);
    group.createAttribute("float_tolerance", float_tolerance);
    return ok;
  }
  // Input from HDF5 file/object
  template <class HF>
  static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>,
                          BrillouinZone>
  from_hdf(HF &obj, const std::string &entry) {
    auto group = obj.getGroup(entry);
    auto lat = lattice_t::from_hdf(group, "inner");
    auto olat = lattice_t::from_hdf(group, "outer");
    auto poly = BrillouinZone::poly_t::from_hdf(group, "first");
    auto ir_p = BrillouinZone::poly_t::from_hdf(group, "irreducible");
    auto ir_w = bArray<double>::from_hdf(group, "ir_wedge_normals");
    int tr, hi, ip, nim, tol;
    double f_tol;
    group.getAttribute("time_reversal").read(tr);
    group.getAttribute("has_inversion").read(hi);
    group.getAttribute("is_primitive").read(ip);
    group.getAttribute("no_ir_mirroring").read(nim);
    group.getAttribute("approx_tolerance").read(tol);
    group.getAttribute("float_tolerance").read(f_tol);
    return {lat, olat, poly, ir_p, ir_w, tr > 0, hi > 0, ip > 0, nim > 0, tol, f_tol};
  }
  [[nodiscard]] bool
  to_hdf(const std::string &filename, const std::string &entry,
         const unsigned perm = HighFive::File::OpenOrCreate) const {
    HighFive::File file(filename, perm);
    return this->to_hdf(file, entry);
  }
  static BrillouinZone from_hdf(const std::string &filename,
                                const std::string &entry) {
    HighFive::File file(filename, HighFive::File::ReadOnly);
    return BrillouinZone::from_hdf(file, entry);
  }
#endif // USE_HIGHFIVE
  bool operator!=(const BrillouinZone &other) const {
    if (time_reversal != other.time_reversal)
      return true;
    if (has_inversion != other.has_inversion)
      return true;
    if (is_primitive != other.is_primitive)
      return true;
    if (no_ir_mirroring != other.no_ir_mirroring)
      return true;
    if (_inner != other._inner)
      return true;
    if (_outer != other._outer)
      return true;
    if (_first != other._first)
      return true;
    if (_irreducible != other._irreducible)
      return true;
    if (ir_wedge_normals != other.ir_wedge_normals)
      return true;
    if (approx_tolerance != other.approx_tolerance)
      return true;
    if (float_tolerance != other.float_tolerance)
      return true;
    return false;
  }
  bool operator==(const BrillouinZone &other) const {
    return !this->operator!=(other);
  }
};

} // end namespace brille
#endif
