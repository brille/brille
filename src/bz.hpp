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
#include <omp.h>
#include "neighbours.hpp"
#include "transform.hpp"
#include "polyhedron.hpp"
#include "phonon.hpp"
#include "approx.hpp"
namespace brille {

/*! \brief An object to hold information about the first Brillouin zone of a Reciprocal lattice

  The BrillouinZone object is created from a Reciprocal lattice and if that
  lattice is not a primitive type the object contains the primitive lattice as
  well.
  By looking for intersections between three (primitive) reciprocal planes, the
  class object identifies the minimal subset of vertices that define the first
  Brillouin zone as well as keeping track of the planes which contributed to
  each vertex and the plane reciprocal vectors themselves.
*/
class BrillouinZone {
  Reciprocal lattice;               //!< The primitive reciprocal lattice
  Reciprocal outerlattice;          //!< The lattice passed in at construction
  Polyhedron polyhedron; //!< The vertices, facet normals, and relation information defining the first Brillouin zone polyhedron
  Polyhedron ir_polyhedron; //!< The vertices, facet normals, facet points, and relation information defining the irreducible first Bz polyhedron
  bArray<double> ir_wedge_normals; //!< The normals of the irreducible reciprocal space wedge planes.
  bool time_reversal; //!< A flag to indicate if time reversal symmetry should be included with pointgroup operations
  bool has_inversion; //!< A computed flag indicating if the pointgroup has space inversion symmetry or if time reversal symmetry has been requested
  bool is_primitive; //!< A computed flag indicating if the primitive version of a conventional lattice is in use
  bool no_ir_mirroring;
public:
  /*!
  @param lat A Reciprocal lattice
  @param toprim A flag to indicate if the primtive lattice should be used
  @param extent An integer to control how-far the vertex-finding algorithm should
                search in τ-index. The default indicates that (̄1̄1̄1),
                (̄1̄10), (̄1̄11), (̄10̄1), ..., (111) are included.
  */
  BrillouinZone(const Reciprocal& lat,
                const bool toprim=true,
                const int extent=1,
                const bool tr=false,
                const bool wedge_search=true
               ):
  lattice(toprim ? lat.primitive() : lat), outerlattice(lat), time_reversal(tr)
  {
    profile_update("Start of BrillouinZone construction");
    this->is_primitive = !(this->lattice.issame(this->outerlattice));
    this->has_inversion = this->time_reversal || lat.has_space_inversion();
    this->no_ir_mirroring = true;
    double old_volume = -1.0, new_volume=0.0;
    int test_extent = extent-1;
    // initial test_extent based on spacegroup or pointgroup?
    while (!brille::approx::scalar(new_volume, old_volume)){
      old_volume = new_volume;
      this->voro_search(++test_extent);
      new_volume = this->polyhedron.get_volume();
    }
    // fallback in case voro_search fails for some reason?!?
    if (brille::approx::scalar(new_volume, 0.)){
      throw std::runtime_error("voro_search failed to produce a non-null first Brillouin zone.");
    } else {
      verbose_update("New polyhedron volume ", this->polyhedron.get_volume());
      verbose_update("First Brillouin zone found using extent ",test_extent,", a ",this->polyhedron.string_repr());
    }
    // in case we've been asked to perform a wedge search for, e.g., P1 or P-1,
    // set the irreducible wedge now as the search will do nothing.
    this->ir_polyhedron = this->polyhedron;
    if (wedge_search){
      bool success = this->wedge_brute_force()
                  || this->wedge_brute_force(false,false) // no special 2-fold or mirror handling
                  || this->wedge_brute_force(false,true) // no special 2-fold handling (but special mirror handling)
                  || this->wedge_brute_force(true, false) // no special mirror handling (maybe not useful)
                  || this->wedge_brute_force(true, true, false) // last ditch effort, handle non order(2) operations in decreasing order
                  || this->wedge_brute_force(true, false, true, false)
                  ;
      // other combinations of special_2_folds, special_mirrors,
      // and sort_by_length are possible but not necessarily useful.
      if (!success)
        throw std::runtime_error("Failed to find an irreducible Brillouin zone.");
    }
    profile_update("  End of BrillouinZone construction");
  }
  void check_if_mirroring_needed(void){
    this->no_ir_mirroring = true;
    if (!this->has_inversion){
      PointSymmetry ps = this->outerlattice.get_pointgroup_symmetry(this->time_reversal?1:0);
      double goal = this->polyhedron.get_volume() / static_cast<double>(ps.size());
      double found = this->ir_polyhedron.get_volume();
      if (brille::approx::scalar(goal, 2.0*found)){
        /*The current Polyhedron at this->ir_polyhedron has half the anticipated
        volume. We can 'fix' this the easy way by mirroring the Polyhedron and
        gluing it onto the current version; the resultant polyhedron will come
        to a point at (0,0,0) and will not be convex. Instead, try to be clever
        and look for a convex combination of the found polyhedron and its
        mirror with one of the pointgroup operations applied:*/
        std::vector<Polyhedron> all_unions;
        Polyhedron mirrored = this->ir_polyhedron.mirror();
        for (auto & r: ps.getall())
          all_unions.push_back(this->ir_polyhedron + mirrored.rotate(r));
        // The combination of a polyhedron and its rotated inverse which has
        // the least number of vertices is (hopefully) convex.
        auto min_vert_union = std::min_element(
          all_unions.begin(),all_unions.end(),
          [](const Polyhedron & a, const Polyhedron & b){
            return a.num_vertices() < b.num_vertices();
          }
        );
        // Polyhedron::operator+() needs to be fixed. As of now it does not look
        // for extraneous faces (like internal faces which are coplanar and
        // pointing in opposite directions) or coplanar external faces which
        // share an edge. Until this is fixed cross our fingers and hope that we
        // created a convex polyhedron such that the convex hull of its points
        // gives the same polyhedron back:
        Polyhedron mvu_convex_hull(min_vert_union->get_vertices());
        if (brille::approx::scalar(goal, mvu_convex_hull.get_volume())){
          // we found a polyhedron with the right volume which is convex!
          // so we can keep this as *the* ir_polyhedron
          this->ir_polyhedron = mvu_convex_hull;
          // and we need to update the found volume for the check below
          found = goal;
        }
      }
      // if found == goal no mirroring is required.
      // if found == goal/2, mirroring is required.
      this->no_ir_mirroring = brille::approx::scalar(goal, found);
    }
  }

  bool check_ir_polyhedron(void);
  bool wedge_explicit(void);
  //! Returns the lattice passed in at construction
  const Reciprocal get_lattice() const { return this->outerlattice;};
  //! Returns the lattice actually used to find the Brillouin zone vertices,
  //! which may be a primitive lattice depending on the flag at creation
  const Reciprocal get_primitive_lattice() const { return this->lattice;};
  //! Returns the number of vertices defining the first Brillouin zone
  size_t vertices_count() const { return this->get_vertices().size(0);};
  //! Returns the number of reciprocal lattice points defining the first Brillouin zone
  size_t faces_count() const { return this->get_points().size(0);};
  // irreducible reciprocal space wedge
  void
  set_ir_wedge_normals(const LQVec<double>& x){
    bool already_same = this->outerlattice.issame(x.get_lattice());
    LQVec<double> xp(this->outerlattice);
    PrimitiveTransform PT(this->outerlattice.get_bravais_type());
    bool transform_needed = ( PT.does_anything() && this->lattice.issame(x.get_lattice()) );
    if (!(already_same || transform_needed))
      throw std::runtime_error("ir_wedge_normals must be in the standard or primitive lattice used to define the BrillouinZone object");
    if (transform_needed)  xp = transform_from_primitive(this->outerlattice,x);
    const LQVec<double> & xref = transform_needed ? xp : x;
    this->ir_wedge_normals = xref.get_hkl();
  }
  LQVec<double> get_ir_wedge_normals() const;
  LQVec<double> get_primitive_ir_wedge_normals() const;
  /*!
  Set the first Brillouin zone polyhedron from its vertices, central facet plane
  points, and either the three intersecting planes which gave each vertex or
  all intersecting planes which give each vertex as well as all vertices which
  form a corner of each facet plane polygon.
  @param vertices All vertices in the polyhedron
  @param points All (τ/2) plane points which define the facets of the polyhedron
  @param [optional] fpv A std::vector<std::vector<int>> with *all* facet-plane indices for each vertex
  @param [optional] vpf A std::vector<std::vector<int>> with *all* vertex indices for each facet-plane
  */
  template<typename... A>
  void
  set_polyhedron(const LQVec<double>& v, const LQVec<double>& p, A... args){
    bool both_same = v.get_lattice().issame(p.get_lattice());
    if (!both_same)
      throw std::runtime_error("The vertices, and points of a polyhedron must all be in the same cooridinate system");
    bool is_outer = this->outerlattice.issame(v.get_lattice());
    bool is_inner = this->lattice.issame(v.get_lattice());
    LQVec<double> vp(this->outerlattice), pp(this->outerlattice);
    PrimitiveTransform PT(this->outerlattice.get_bravais_type());
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
  /*!
  Set the irreducible first Brillouin zone polyhedron from its vertices, facet
  plane points, facet plane normals, and, optionally, all intersecting planes
  which give each vertex as well as all vertices which form a corner of each
  facet plane polygon.
  @param vertices All vertices in the polyhedron
  @param points A point on each of the planes which define the facets of the polyhedron
  @param normals The normal direction which, along with the on-plane-point, defines each facet plane
  @param [optional] fpv A std::vector<std::vector<int>> with *all* facet-plane indices for each vertex
  @param [optional] vpf A std::vector<std::vector<int>> with *all* vertex indices for each facet-plane
  */
  template<typename... A>
  void
  set_ir_polyhedron(const LQVec<double>& v, const LQVec<double>& p, const LQVec<double>& n, A... args){
    bool all_same = v.get_lattice().issame(p.get_lattice()) && p.get_lattice().issame(n.get_lattice());
    if (!all_same)
      throw std::runtime_error("The vertices, points, and normals of a polyhedron must all be in the same cooridinate system");
    bool is_outer = this->outerlattice.issame(v.get_lattice());
    bool is_inner = this->lattice.issame(v.get_lattice());
    LQVec<double> vp(this->outerlattice), pp(this->outerlattice), np(this->outerlattice);
    PrimitiveTransform PT(this->outerlattice.get_bravais_type());
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
  /*!
  Set the irreducible first Brillouin zone polyhedron from its vertices. Since
  no facet information is provided via this method, a convex hull of the
  provided points is calculated and used as the irreducible polyhedron.
  Additionally the volume and symmetry of the calculated polyhedron is checked
  against the first Brillouin zone using the pointgroup symmetry operations.
  @param vertices All vertices in the irreducible polyhedron
  @returns A bool indicating if the found convex hull polyhedron is the
           expected volume and has the expected symmetry.                     */
  bool set_ir_vertices(const LQVec<double>& v){
    bool is_outer = this->outerlattice.issame(v.get_lattice());
    bool is_inner = this->lattice.issame(v.get_lattice());
    LQVec<double> vp(this->outerlattice);
    PrimitiveTransform PT(this->outerlattice.get_bravais_type());
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
  //! Returns the first Brillouin zone polyhedron
  Polyhedron get_polyhedron(void) const;
  //! Returns the vertices of the first Brillouin zone polyhedron expressed as conventional unit cell vectors
  LQVec<double> get_vertices(void) const;
  //! Returns the on-plane points (τ/2) of the first Brillouin zone polyhedron expressed as conventional unit cell vectors
  LQVec<double> get_points(void) const;
  //! Returns the normals ( ̂τ ) of the first Brillouin zone polyhedron expressed as conventional unit cell vectors
  LQVec<double> get_normals(void) const;
  //! Returns the half-edge points of the first Brillouin zone polyhedron expressed as conventional unit cell vectors
  LQVec<double> get_half_edges(void) const;
  //! Returns the vertices of the first Brillouin zone polyhedron expressed as primitive unit cell vectors
  LQVec<double> get_primitive_vertices(void) const;
  //! Returns the on-plane points (τ/2) of the first Brillouin zone polyhedron expressed as primitive unit cell vectors
  LQVec<double> get_primitive_points(void) const;
  //! Returns the normals ( ̂τ ) of the first Brillouin zone polyhedron expressed as primitive unit cell vectors
  LQVec<double> get_primitive_normals(void) const;
  //! Returns the `points` and `normals` indices for each vertex of the first Brillouin zone polyhedron
  std::vector<std::vector<int>> get_faces_per_vertex(void) const;
  //! Returns the vertex indices for each facet of the first Brillouin zone polyhedron
  std::vector<std::vector<int>> get_vertices_per_face(void) const;
  //! Returns the irreducible first Brillouin zone polyhedron
  Polyhedron get_ir_polyhedron(const bool true_ir=true) const;
  //! Returns the vertices of the irreducible first Brillouin zone polyhedron expressed as conventional unit cell vectors
  LQVec<double> get_ir_vertices(void) const;
  //! Returns the on-plane points of the irreducible first Brillouin zone polyhedron expressed as conventional unit cell vectors
  LQVec<double> get_ir_points(void) const;
  //! Returns the normals of the irreducible first Brillouin zone polyhedron expressed as conventional unit cell vectors
  LQVec<double> get_ir_normals(void) const;
  //! Returns the vertices of the irreducible first Brillouin zone polyhedron expressed as primitive unit cell vectors
  LQVec<double> get_ir_primitive_vertices(void) const;
  //! Returns the on-plane points of the irreducible first Brillouin zone polyhedron expressed as primitive unit cell vectors
  LQVec<double> get_ir_primitive_points(void) const;
  //! Returns the normals of the irreducible first Brillouin zone polyhedron expressed as primitive unit cell vectors
  LQVec<double> get_ir_primitive_normals(void) const;
  //! Returns the `points` and `normals` indices for each vertex of the irreducible first Brillouin zone polyhedron
  std::vector<std::vector<int>> get_ir_faces_per_vertex(void) const;
  //! Returns the vertex indices for each facet of the irreducible first Brillouin zone polyhedron
  std::vector<std::vector<int>> get_ir_vertices_per_face(void) const;
  //! Print a representation of the object to std::cout
  void print() const;
  /*!
  Using the pointgroup symmetry opperations of the conventional unit cell,
  this method defines *an* irreducible section of reciprocal space.
  The irreducible section bounds a solid angle N and consequently the full
  reciprocal space is 4π/N larger than the irreducible wedge.
  As the selected irreducible section is related by symmetry to the rest of
  reciprocal space there are (at least) 4π/N - 1 equivalent other choices for
  the irreducible wedge.
  */
  void wedge_search(const bool prefer_basis_vectors=true, const bool parallel_ok=true);
  bool wedge_brute_force(bool special_2_folds = true, bool special_mirrors = true, bool sort_by_length=true, bool sort_one_sym=true);
  void wedge_triclinic(void);
  /*!
  With the first Brillouin zone and *an* irreducible section of reciprocal space
  already identified, this method finds all intersections of combinations of
  three planes (nBz,nIr) = {(0,3), (1,2), (2,1), and (3,0)} where nBz are the
  number of intersecting first Brillouin zone facet planes and nIr are the
  number of intersecting irreducible reciprocal space wedge planes.
  Once all combinations of intersection points are found, those within the
  first Brillouin zone *and* the irreducible reciprocal space wedge define the
  irreducible first Brillouin zone.
  @note the (3,0) intersection points *are* the first Brillouin zone vertices
  and the (0,3) intersection points are all identically ⃗0.
  */
  void irreducible_vertex_search();
  /*! \brief Search for the vertices defining the first Brillouin zone

    Search through all combinations of three planes from `extent`*(̄1̄1̄1) to
    `extent`*(111) for three-plane intersection points, these are all the
    vertices of *a* Brillouin zone. Sort through the intersection points to
    find the subset that are closer to the origin than any of the
    non-intersection-defining planes, these are vertices of the *first*
    Brillouin zone but may not be unique. So determine a unique subset and
    record their locations, contributing intersecting-planes, and all planes
    which contribute to one or more vertex in the object.
  */
  void voro_search(const int extent=1);
  //! Return true if the lattice used to find the vertices is not the same as the one passed at construction, and is therefore primitive
  bool isprimitive(void) const {return this->is_primitive;};

  /*! \brief Determine whether points are inside of the first Brillouin zone
    @param p A reference to a LQVec list of Q points to be checked
    @returns A std::vector<bool> with each element indicating if the
             associated Q point is inside of the Brillouin zone.
  */
  template<class T>
  std::vector<bool>
  isinside(const LQVec<T>& p) const {
    bool isouter = this->outerlattice.issame(p.get_lattice());
    bool isinner = this->lattice.issame(p.get_lattice());
    if (!(isouter||isinner))
      throw std::runtime_error("Q points must be in the standard or primitive lattice");
    std::vector<bool> out(p.size(0), true);
    LQVec<double> points, normals;
    if (isouter){
      points = this->get_points();
      normals = this->get_normals();
    } else {
      points = this->get_primitive_points();
      normals = this->get_primitive_normals();
    }
    #pragma omp parallel for default(none) shared(out, normals, p, points) schedule(dynamic)
    for (long long i=0; i<p.size(0); ++i)
      out[i] = dot(normals, p.view(i)-points).all(brille::cmp::le, 0.);
    return out;
  }
  /*! \brief Determine whither points are inside the irreducible reciprocal space wedge
  @param p A reference to a LQVec list of Q points to be checked
  @returns An std::vector<bool> with each element indicating if the
           associated Q point is inside our irreducible reciprocal space.
  */
  template<class T>
  std::vector<bool>
  isinside_wedge(const LQVec<T> &p, const bool constructing=false) const {
    bool isouter = this->outerlattice.issame(p.get_lattice());
    bool isinner = this->lattice.issame(p.get_lattice());
    if (!(isouter||isinner)){
      std::string msg = "Q points provided to BrillouinZone::isinside_wedge ";
      msg += "must be in the standard or primitive lattice ";
      msg += "used to define the BrillouinZone object";
      throw std::runtime_error(msg);
    }
    std::vector<bool> out(p.size(0), true);
    LQVec<double> normals;
    if (isouter)
      normals = this->get_ir_wedge_normals();
    else
      normals = this->get_primitive_ir_wedge_normals();
    if (normals.size(0)){ // with no normals *all* points are "inside" the wedge
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
      // The Array all method has a special switched execution path
      // for checking whether all values are ≤ or ≥ a value simultaneously

      // when constructing the irreducible Brillouin zone we need to only consider
      // the ≥0 case so that we end up with a convex polyhedron. The ir_polyhedron
      // accessor method mirrors the half-polyhedron in this case, so when
      // identifying whether a point is inside of the irreducible Brillouin zone
      // we must allow for the ≤0 case as well.
      brille::cmp c = constructing||this->no_ir_mirroring ? brille::cmp::ge : brille::cmp::le_ge;
      #pragma omp parallel for default(none) shared(out, normals, p, c) schedule(dynamic)
      for (long long i=0; i<p.size(0); ++i)
        out[i] = dot(normals, p.view(i)).all(c, 0.);
    }
    return out;
  }
  /*! \brief Find q and τ such that Q=q+τ and τ is a reciprocal lattice vector
    @param[in] Q A reference to LQVec list of Q points
    @param[out] q The reduced reciprocal lattice vectors
    @param[out] tau The reciprocal lattice zone centres
  */
  bool
  moveinto(const LQVec<double>& Q, LQVec<double>& q, LQVec<int>& tau, const int threads=0) const
  {
    profile_update("BrillouinZone::moveinto called with ",threads," threads");
    omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
    bool already_same = this->lattice.issame(Q.get_lattice());
    LQVec<double> Qprim(this->lattice);
    LQVec<double> qprim(this->lattice);
    LQVec<int> tauprim(this->lattice);
    PrimitiveTransform PT(this->outerlattice.get_bravais_type());
    bool transform_needed = ( PT.does_anything() && this->outerlattice.issame(Q.get_lattice()) );
    if (!(already_same || transform_needed)){
      std::string msg = "Q points provided to BrillouinZone::moveinto must be ";
      msg += "in the standard or primitive lattice used to define ";
      msg += "the BrillouinZone object";
      throw std::runtime_error(msg);
    }
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
    size_t max_count = taus.size(0);
    // ensure that qsl and tausl can hold each qi and taui
    qsl.resize(Qsl.size(0));
    tausl.resize(Qsl.size(0));
    long long snQ = brille::utils::u2s<long long, size_t>(Qsl.size(0));
  #pragma omp parallel for default(none)\
  shared(Qsl, tausl, qsl, points, normals, taus, taulen, snQ, max_count)\
  schedule(dynamic)
    for (long long si=0; si<snQ; si++){
      size_t i = brille::utils::s2u<size_t, long long>(si);
      auto taui = Qsl.view(i).round();
      auto qi = Qsl.view(i) - taui;
      auto last_shift = taui;
      size_t count{0};
      while (count++ < max_count && dot(normals, qi-points).any(brille::cmp::gt,0.)){
        auto qi_dot_normals = dot(qi , normals);
        auto Nhkl = (qi_dot_normals/taulen).round().to_std();
        auto qidn = qi_dot_normals.to_std();
        if (std::any_of(Nhkl.begin(), Nhkl.end(), [](int a){return a > 0;})){
          int maxnm{0};
          size_t maxat{0};
          for (size_t j=0; j<Nhkl.size(); ++j)
                                                           // protect against oscillating by ±τ
          if (Nhkl[j]>0 && Nhkl[j]>=maxnm && (0==maxnm || (norm(taus.view(j)+last_shift).all(brille::cmp::gt, 0.) && qidn[j]>qidn[maxat]))){
            maxnm = Nhkl[maxat=j];
          }
          qi -= taus.view(maxat) * static_cast<double>(maxnm); // ensure we subtract LQVec<double>
          taui += taus.view(maxat) * maxnm; // but add LQVec<int>
          last_shift = taus.view(maxat) * maxnm;
        }
      }
      qsl.set(i, qi);
      tausl.set(i, taui);
    }
    if (transform_needed){ // then we need to transform back q and tau
      q   = transform_from_primitive(this->outerlattice,qsl);
      tau = transform_from_primitive(this->outerlattice,tausl);
    }
    auto allinside = this->isinside(q);
    if (std::count(allinside.begin(), allinside.end(), false) > 0){
      for (size_t i=0; i<Q.size(0); ++i) if (!allinside[i]){
        info_update("Q  =",Q.to_string(i)  ," tau  =",tau.to_string(i)  ," q  =",q.to_string(i));
        info_update("Qsl=",Qsl.to_string(i)," tausl=",tausl.to_string(i)," qsl=",qsl.to_string(i),"\n");
      }
      throw std::runtime_error("Not all points inside Brillouin zone");
      // return false;
    }
    return true; // otherwise an error has been thrown
  }
  /*! \brief Find q, τ, and R∈G such that Q = Rᵀq + τ, where τ is a reciprocal
             lattice vector and R is a pointgroup symmetry operation of the
             conventional unit cell pointgroup, G.
    @param [in] Q a refernce to a LQVec list of Q points
    @param [out] q The irreducible reduced reciprocal lattice vectors
    @param [out] τ The conventional reciprocal lattice zone centres
    @param [out] R The pointgroup operation index for R
    @param [out] invR The pointgroup operation index for R⁻¹
    @param [in] nthreads An optional number of OpenMP threads to use
    @note R and invR index the PointSymmetry object accessible via BrillouinZone::get_pointgroup_symmetry();
  */
  bool
  ir_moveinto(
    const LQVec<double>& Q, LQVec<double>& q, LQVec<int>& tau,
    std::vector<size_t>& Ridx, std::vector<size_t>& invRidx, const int threads=0)
  const {
    profile_update("BrillouinZone::ir_moveinto called with ",threads," threads");
    omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
    /* The Pointgroup symmetry information comes from, effectively, spglib which
    has all rotation matrices defined in the conventional unit cell -- which is
    our `outerlattice`. Consequently we must work in the outerlattice here.  */
    if (!this->outerlattice.issame(Q.get_lattice()))
      throw std::runtime_error("Q points provided to ir_moveinto must be in the standard lattice used to define the BrillouinZone object");
    // get the PointSymmetry object, containing all operations
    PointSymmetry psym = this->get_pointgroup_symmetry();
    // ensure q, tau, and Rm can hold one for each Q.
    size_t nQ = Q.size(0);
    auto Qshape = Q.shape();
    q.resize(Qshape);
    tau.resize(Qshape);
    Ridx.resize(nQ);
    invRidx.resize(nQ);
    // find q₁ₛₜ in the first Brillouin zone and τ ∈ [reciprocal lattice vectors]
    // such that Q = q₁ₛₜ + τ
    this->moveinto(Q, q, tau, threads);
    // by chance some first Bz points are likely already in the IR-Bz:
    std::vector<bool> in_ir = this->isinside_wedge(q);
    //LQVec<double> qj(Q.get_lattice(), 1u);
    auto lat = Q.get_lattice();
    // OpenMP 2 (VS) doesn't like unsigned loop counters
    size_t n_outside{0};
    long long snQ = brille::utils::u2s<long long, size_t>(nQ);
    #pragma omp parallel for default(none) shared(psym, Ridx, invRidx, q, in_ir, lat, snQ) reduction(+:n_outside) schedule(dynamic)
    for (long long si=0; si<snQ; ++si){
      size_t i = brille::utils::s2u<size_t, long long>(si);
      // any q already in the irreducible zone need no rotation → identity, but we need to find the index of E
      bool outside=!in_ir[i];
      if (outside){
        // for others find the jᵗʰ operation which moves qᵢ into the irreducible zone
        LQVec<double> qj(lat, 1u); // a place to hold the multiplication result
        for (size_t j=0; j<psym.size(); ++j) if (outside) {
          // The point symmetry matrices relate *real space* vectors! We must use
          // their transposes' to rotate reciprocal space vectors.
          brille::utils::multiply_matrix_vector(qj.ptr(0), transpose(psym.get(j)).data(), q.ptr(i));
          if ( this->isinside_wedge(qj)[0] ){ /* store the result */
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
  /*! \brief Find q and R∈G such that Q = Rᵀq such that q ∈ the irreducible wedge
  */
  bool
  ir_moveinto_wedge(const LQVec<double>& Q, LQVec<double>& q, std::vector<std::array<int,9>>& R, const int threads=0) const {
    omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
    /* The Pointgroup symmetry information comes from, effectively, spglib which
    has all rotation matrices defined in the conventional unit cell -- which is
    our `outerlattice`. Consequently we must work in the outerlattice here.  */
    if (!this->outerlattice.issame(Q.get_lattice()))
      throw std::runtime_error("Q points provided to ir_moveinto must be in the standard lattice used to define the BrillouinZone object");
    // get the PointSymmetry object, containing all operations
    PointSymmetry psym = this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
    // ensure q and R can hold one for each Q.
    size_t nQ = Q.size(0);
    auto Qshape = Q.shape();
    q.resize(Qshape);
    R.resize(nQ);
    // by chance some first Bz points are likely already in the IR-wedge:
    std::vector<bool> in_ir = this->isinside_wedge(Q);
    //LQVec<double> qj(Q.get_lattice(), 1u);
    auto lat = Q.get_lattice();
    // OpenMP 2 (VS) doesn't like unsigned loop counters
    size_t n_outside{0};
    long long snQ = brille::utils::u2s<long long, size_t>(nQ);
    #pragma omp parallel for default(none) shared(psym, R, q, Q, in_ir, lat, snQ) reduction(+:n_outside) schedule(dynamic)
    for (long long si=0; si<snQ; ++si){
      size_t i = brille::utils::s2u<size_t, long long>(si);
      // any q already in the irreducible zone need no rotation → identity
      bool outside=!in_ir[i];
      if (outside){
        // for others find the jᵗʰ operation which moves qᵢ into the irreducible zone
        LQVec<double> qj(lat, 1u); // a place to hold the multiplication result
        for (size_t j=0; j<psym.size(); ++j) if (outside) {
          // The point symmetry matrices relate *real space* vectors! We must use
          // their transposes' to rotate reciprocal space vectors.
          brille::utils::multiply_matrix_vector(qj.ptr(0), transpose(psym.get(j)).data(), Q.ptr(i));
          if ( this->isinside_wedge(qj)[0] ){ /* store the result */
            q.set(i, qj); // keep Rⱼᵀ⋅Qᵢ as qᵢᵣ
            R[i] = transpose(psym.get_inverse(j)); // and (Rⱼᵀ)⁻¹ ∈ G, such that Q = (Rⱼᵀ)⁻¹⋅qᵢᵣ
            outside = false;
          }
        }
      } else {
        q.set(i, Q.view(i));
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

  //! \brief Get the PointSymmetry object used by this BrillouinZone object internally
  PointSymmetry get_pointgroup_symmetry() const{
    return this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
  }
  //! \brief Accessor for whether the BrillouinZone was constructed with additional time reversal symmetry
  int add_time_reversal() const {
    return this->time_reversal ? 1 : 0;
  }
private:
  template<class T, class R>
  bool
  wedge_normal_check(const LQVec<T>& n, LQVec<R>& normals, size_t& num){
    std::string msg = "Considering " + n.to_string(0) + "... ";
    if (norm(n).all(brille::cmp::eq, 0.0)){
      debug_update(msg, "rejected; zero-length");
      return false;
    }
    if (num==0){
      debug_update(msg, "accepted; first normal");
      normals.set(num, n.view(0));
      num=num+1;
      return true;
    }
    // num > 0 for all following views
    if (norm(cross(normals.view(0,num), n)).any(brille::cmp::eq, 0.)){
      debug_update(msg, "rejected; already present");
      return false;
    }
    normals.set(num,  n.view(0));
    if (this->ir_wedge_is_ok(normals.view(0,num+1))){
      debug_update(msg, "accepted");
      num=num+1;
      return true;
    }
    normals.set(num, -n.extract(0));
    if (this->ir_wedge_is_ok(normals.view(0,num+1))){
      debug_update(msg, "accepted (*-1)");
      num=num+1;
      return true;
    }
    debug_update(msg, "rejected; addition causes null wedge");
    return false;
  }
  template<class T, class R, class U>
  bool
  wedge_normal_check(const LQVec<T>& n0, const LQVec<R>& n1, LQVec<U>& normals, size_t& num){
    std::string msg = "Considering " + n0.to_string(0)+ " and " + n1.to_string(0) + "... ";
    bool p0=false, p1=false;
    if (num>0){
      p0 = norm(cross(normals.view(0,num), n0)).any(brille::cmp::eq,0.);
      p1 = norm(cross(normals.view(0,num), n1)).any(brille::cmp::eq,0.);
    }
    if (p0 && p1){
      debug_update(msg, "rejected; already present");
      return false;
    }
    if (num>0 && dot(normals.view(0,num)/norm(n0),n0/norm(n0)).any(brille::cmp::eq,1.)){
      normals.set(num,  n1.view(0));
      if (this->ir_wedge_is_ok(normals.view(0,num+1))){
        debug_update(msg, "n1 accepted (n0 present)");
        num=num+1;
        return true;
      }
      debug_update(msg, "n1 rejected (n0 present); addition causes null wedge");
      return false;
    }
    if (num>0 && dot(normals.view(0,num)/norm(n1),n1/norm(n1)).any(brille::cmp::eq,1.)){
      normals.set(num,  n0.view(0));
      if (this->ir_wedge_is_ok(normals.view(0,num+1))){
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
    normals.set(num,   n0.view(0));
    normals.set(num+1, n1.view(0));
    if (this->ir_wedge_is_ok(normals.view(0,num+2))){
      debug_update(msg, "n0 & n1 accepted");
      num=num+2;
      return true;
    }
    debug_update(msg, "n0 & n1 rejected; adding both would cause null wedge");
    return false;
  }
  void
  shrink_and_prune_outside(const size_t cnt, LQVec<double>& vrt, bArray<int>& ijk) const {
    verbose_update("shrinking to ",cnt);
    if(vrt.size(0) && ijk.size(0)){
      vrt.resize(cnt);
      ijk.resize(cnt);
      if (cnt){ // isinside has a problem with vrt.size()==0
        auto isin = this->isinside(vrt);
        verbose_update("and retaining ", std::count(isin.begin(), isin.end(), true), " inside vertices");
        vrt = vrt.extract(isin);
        ijk = ijk.extract(isin);
      }
    }
  }
  template<class T>
  bool
  ir_wedge_is_ok(const LQVec<T>& normals){
    this->set_ir_wedge_normals(normals); // assigns this->ir_wedge_normals
    this->irreducible_vertex_search(); // assigns this->ir_polyhedron
    return !brille::approx::scalar(this->ir_polyhedron.get_volume(), 0.0);
  }
  LQVec<double> get_ir_polyhedron_wedge_normals(void) const;
};

/*! \brief Find the intersection point of three planes, if it exists.

From the normal vector and on-plane point describing three planes, determine
if an intersection point exists (and is not too far from the origin), if so
calculate the intersection pointd and store it in `intersect` at `idx`.
@param ni A LQVec<double> reference to the normal for plane `i`
@param pi A LQVec<double> reference to the on-plane point for plane `i`
@param nj A LQVec<double> reference to the normal for plane `j`
@param pj A LQVec<double> reference to the on-plane point for plane `j`
@param nk A LQVec<double> reference to the normal for plane `k`
@param pk A LQVec<double> reference to the on-plane point for plane `k`
@param[out] intersect A LQVec<double> reference to a list of intersection points
@param idx The location in `intersect` where the found intersection point should
           be stored, if it exists and is not too far from the origin.
*/
bool
intersect_at(const LQVec<double>& ni, const LQVec<double>& pi,
             const LQVec<double>& nj, const LQVec<double>& pj,
             const LQVec<double>& nk, const LQVec<double>& pk,
             LQVec<double>& intersect, const int idx);
/*! \brief A specialization where one plane passes through the origin

@param ni A LQVec<double> reference to the normal for plane `i`
@param pi A LQVec<double> reference to the on-plane point for plane `i`
@param nj A LQVec<double> reference to the normal for plane `j`
@param pj A LQVec<double> reference to the on-plane point for plane `j`
@param nk A LQVec<double> reference to the normal for plane `k`, which passes through the origin
@param[out] intersect A LQVec<double> reference to a list of intersection points
@param idx The location in `intersect` where the found intersection point should
           be stored, if it exists and is not too far from the origin.
*/
bool
intersect_at(const LQVec<double>& ni, const LQVec<double>& pi,
             const LQVec<double>& nj, const LQVec<double>& pj,
             const LQVec<double>& nk,
             LQVec<double>& intersect, const int idx);
/*! \brief A specialization where two planes pass through the origin

@param ni A LQVec<double> reference to the normal for plane `i`
@param pi A LQVec<double> reference to the on-plane point for plane `i`
@param nj A LQVec<double> reference to the normal for plane `j`, which passes through the origin
@param nk A LQVec<double> reference to the normal for plane `k`, which passes through the origin
@param[out] intersect A LQVec<double> reference to a list of intersection points
@param idx The location in `intersect` where the found intersection point should
           be stored, if it exists and is not too far from the origin.
*/
bool
intersect_at(const LQVec<double>& ni, const LQVec<double>& pi,
             const LQVec<double>& nj,
             const LQVec<double>& nk,
             LQVec<double>& intersect, const int idx);
/*! \brief A (trivial) specialization where three planes pass through the origin

@param ni A LQVec<double> reference to the normal for plane `i`, which passes through the origin
@param nj A LQVec<double> reference to the normal for plane `j`, which passes through the origin
@param nk A LQVec<double> reference to the normal for plane `k`, which passes through the origin
@param[out] intersect A LQVec<double> reference to a list of intersection points
@param idx The location in `intersect` where the found intersection point should
           be stored, if it exists and is not too far from the origin.
@note The intersection point of three non-degenerate planes passing through the
      origin *is* the origin. Degenerate planes do not intersect at a point, so
      this function will return false for degenerate inputs.
*/
bool
intersect_at(const LQVec<double>& ni,
             const LQVec<double>& nj,
             const LQVec<double>& nk,
             LQVec<double>& intersect, const int idx);
/*! \brief Construct and return the determinant of a 3D normals matrix

Used by the intersect_at overloaded functions to find the determinant of the
normals matrix which determines *if* three planes intersect.

@param a A LQVec<double> reference to the first plane normal
@param b A LQVec<double> reference to the second plane normal
@param c A LQVec<double> reference to the third plane normal
*/
double
normals_matrix_determinant(const LQVec<double>& a, const LQVec<double>&b, const LQVec<double>& c);

} // end namespace brille
#endif
