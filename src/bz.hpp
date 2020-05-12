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

#ifndef _BZ_CLASS_H_
#define _BZ_CLASS_H_

#include "neighbours.hpp"
#include "transform.hpp"
// #include "pointgroup.hpp"
// #include <iostream>
// #include <algorithm>
// #include <vector>
#include "polyhedron.hpp"
// #include "debug.hpp"
#include "phonon.hpp"

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
  ArrayVector<double> ir_wedge_normals; //!< The normals of the irreducible reciprocal space wedge planes.
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
    this->is_primitive = !(this->lattice.issame(this->outerlattice));
    this->has_inversion = this->time_reversal || lat.has_space_inversion();
    this->no_ir_mirroring = true;
    double old_volume = -1.0, new_volume=0.0;
    int test_extent = extent-1;
    // initial test_extent based on spacegroup or pointgroup?
    while (!approx_scalar(new_volume, old_volume)){
      old_volume = new_volume;
      this->voro_search(++test_extent);
      new_volume = this->polyhedron.get_volume();
    }
    // fallback in case voro_search fails for some reason?!?
    if (approx_scalar(new_volume, 0.)){
      info_update("voro_search failed to produce a non-null first Brillouin zone.");
      this->vertex_search(extent);
    } else {
      verbose_update("New polyhedron volume ", this->polyhedron.get_volume());
      verbose_update("First Brillouin zone found using extent ",test_extent,", a ",this->polyhedron.string_repr());
    }
    // in case we've been asked to perform a wedge search for, e.g., P1 or P-1,
    // set the irreducible wedge now as the search will do nothing.
    this->ir_polyhedron = this->polyhedron;
    if (wedge_search){
      this->wedge_brute_force();
      if (!this->check_ir_polyhedron()) this->wedge_brute_force(false,false); // no special 2-fold or mirror handling
      if (!this->check_ir_polyhedron()) this->wedge_brute_force(false,true); // no special 2-fold handling (but special mirror handling)
      if (!this->check_ir_polyhedron()) this->wedge_brute_force(true, false); // no special mirror handling (maybe not useful)
      if (!this->check_ir_polyhedron()) this->wedge_brute_force(true, true, false); // last ditch effort, handle non order(2) operations in decreasing order
      // other combinations of special_2_folds, special_mirrors,
      // and sort_by_length are possible but not necessarily useful.
      if (!this->check_ir_polyhedron())
        info_update("Failed to find an irreducible Brillouin zone.");
    }
  }
  void check_if_mirroring_needed(void){
    this->no_ir_mirroring = true;
    if (!this->has_inversion){
      PointSymmetry ps = this->outerlattice.get_pointgroup_symmetry(this->time_reversal?1:0);
      double goal = this->polyhedron.get_volume() / static_cast<double>(ps.size());
      double found = this->ir_polyhedron.get_volume();
      if (approx_scalar(goal, 2.0*found)){
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
        if (approx_scalar(goal, mvu_convex_hull.get_volume())){
          // we found a polyhedron with the right volume which is convex!
          // so we can keep this as *the* ir_polyhedron
          this->ir_polyhedron = mvu_convex_hull;
          // and we need to update the found volume for the check below
          found = goal;
        }
      }
      // if found == goal no mirroring is required.
      // if found == goal/2, mirroring is required.
      this->no_ir_mirroring = approx_scalar(goal, found);
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
  size_t vertices_count() const { return this->get_vertices().size();};
  //! Returns the number of reciprocal lattice points defining the first Brillouin zone
  size_t faces_count() const { return this->get_points().size();};
  // irreducible reciprocal space wedge
  void set_ir_wedge_normals(const LQVec<double>&);
  LQVec<double> get_ir_wedge_normals() const;
  LQVec<double> get_primitive_ir_wedge_normals() const;
  /*!
  Set the first Brillouin zone polyhedron from its vertices, central facet plane
  points, and either the three intersecting planes which gave each vertex or
  all intersecting planes which give each vertex as well as all vertices which
  form a corner of each facet plane polygon.
  @param vertices All vertices in the polyhedron
  @param points All (τ/2) plane points which define the facets of the polyhedron
  @param fpv An ArrayVector<int> with three facet-plane indices for each vertex
  @param fpv A std::vector<std::vector<int>> with *all* facet-plane indices for each vertex
  @param vpf A std::vector<std::vector<int>> with *all* vertex indices for each facet-plane
  */
  template<typename... A> void set_polyhedron(const LQVec<double>&, const LQVec<double>&, A...);
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
  template<typename... A> void set_ir_polyhedron(const LQVec<double>&, const LQVec<double>&, const LQVec<double>&, A...);
  /*!
  Set the irreducible first Brillouin zone polyhedron from its vertices. Since
  no facet information is provided via this method, a convex hull of the
  provided points is calculated and used as the irreducible polyhedron.
  Additionally the volume and symmetry of the calculated polyhedron is checked
  against the first Brillouin zone using the pointgroup symmetry operations.
  @param vertices All vertices in the irreducible polyhedron
  @returns A bool indicating if the found convex hull polyhedron is the
           expected volume and has the expected symmetry.                     */
  bool set_ir_vertices(const LQVec<double>&);
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
  void wedge_brute_force(bool special_2_folds = true, bool special_mirrors = true, bool sort_by_length=true);
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
  void vertex_search(const int extent=1);
  void voro_search(const int extent=1);
  // void vertex_search_xyz(const int extent=1);
  //! Return true if the lattice used to find the vertices is not the same as the one passed at construction, and is therefore primitive
  bool isprimitive(void) const {return this->is_primitive;};

  /*! \brief Determine whether points are inside of the first Brillouin zone
    @param p A reference to a LQVec list of Q points to be checked
    @returns An ArrayVector<bool> with each 1-element array indicating if the
             associated Q point is inside of the Brillouin zone.
  */
  template<typename T> ArrayVector<bool> isinside(const LQVec<T>& p) const ;
  template<typename T> std::vector<bool> isinside_std(const LQVec<T>& p) const ;
  /*! \brief Determine whither points are inside the irreducible reciprocal space wedge
  @param p A reference to a LQVec list of Q points to be checked
  @returns An ArrayVector<bool> with each 1-element array indicating if the
           associated Q point is inside our irreducible reciprocal space.
  */
  template<typename T> ArrayVector<bool> isinside_wedge(const LQVec<T> &p, const bool constructing=false) const;
  template<typename T> std::vector<bool> isinside_wedge_std(const LQVec<T> &p, const bool constructing=false) const;
  /*! \brief Find q and τ such that Q=q+τ and τ is a reciprocal lattice vector
    @param[in] Q A reference to LQVec list of Q points
    @param[out] q The reduced reciprocal lattice vectors
    @param[out] tau The reciprocal lattice zone centres
  */
  bool moveinto(const LQVec<double>& Q, LQVec<double>& q, LQVec<int>& tau, int nthreads=0) const;
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
  bool ir_moveinto(const LQVec<double>& Q, LQVec<double>& q, LQVec<int>& tau, std::vector<size_t>& Rm, std::vector<size_t>& invRm, int nthreads=0) const ;
  bool ir_moveinto_wedge(const LQVec<double>& Q, LQVec<double>& q, std::vector<std::array<int,9>>& R, int threads=0) const;
  //! \brief Get the PointSymmetry object used by this BrillouinZone object internally
  PointSymmetry get_pointgroup_symmetry() const{
    return this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
  }
  //! \brief Accessor for whether the BrillouinZone was constructed with additional time reversal symmetry
  int add_time_reversal() const {
    return this->time_reversal ? 1 : 0;
  }
private:
  void shrink_and_prune_outside(const size_t cnt, LQVec<double>& vrt, ArrayVector<int>& ijk) const;
  bool wedge_normal_check(const LQVec<double>& n, LQVec<double>& normals, size_t& num);
  bool wedge_normal_check(const LQVec<double>& n0, const LQVec<double>& n1, LQVec<double>& normals, size_t& num);
  bool ir_wedge_is_ok(const LQVec<double>& normals);
  LQVec<double> get_ir_polyhedron_wedge_normals(void) const;
};

/*! \brief Determine whether a given point is between a plane and the origin

For a given list of plane-defining points, which are also their normal vectors,
find whether an indicated vertex is closer to the origin than all planes which
did not define it.
The additional input complexity is to avoid memory copies in the calling code.
@param p A pointer to the LQVec array of all plane-defining centre points
@param v A pointer to a LQVec array of vertices, of which one is checked
@param ijk A pointer to an ArrayVector of the plane-indices for each vertex in v
@param idx The one index of v (and ijk) which will be checked
@param[out] inv A pointer to an LQVec where v[idx] should be stored, if it is found closer to the origin
@param store_at The index into inv where v[idx] should be stored
@param tol An additional absolute tolerance to determine closer-to-originness
@note The function compares `dot(v[idx],p)-dot(p,p)` with zero for all plane
      centres in `p`. If the difference is greater than tol then the plane is
      closer to the origin than the vertex. Rather than relying on solely the
      passed value of `tol` the function also compares the difference to machine
      precision and errs on the side of claiming a vertex is closer than a plane
      when the difference is smaller than the sum of the dot products times epsilon
      or just epsilon, whichever is greater.
*/
bool between_origin_and_plane(const LQVec<double> *p,
                              const LQVec<double> *v,
                              const ArrayVector<int> *ijk,
                              const int idx,
                                    LQVec<double> *inv,
                              const int store_at=0,
                              const double tol=0*1e-15);
/*! \brief Find the intersection point of three planes if it exists

  From vectors to points on three planes and their normals
  determine if an intersection point exists and is not too large.
  @param n A pointer to the plane normal vectors
  @param p A pointer to the points on each plane
  @param xyz A pointer to the plane normals in orthonormal units
  @param i index into `n`, `p`, and `xyz` for the first plane
  @param j index into `n`, `p`, and `xyz` for the second plane
  @param k index into `n`, `p`, and `xyz` for the third plane
  @param[out] iat A LQVec to hold the intersection point, if it exists
  @param idx Where the intersection point should be stored in `iat`
*/
bool three_plane_intersection(const LQVec<double>& n,
                              const LQVec<double>& p,
                              const ArrayVector<double>& xyz,
                              const int i, const int j, const int k,
                              LQVec<double>& iat, const int idx=0);
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
bool intersect_at(const LQVec<double>& ni, const LQVec<double>& pi,
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
bool intersect_at(const LQVec<double>& ni, const LQVec<double>& pi,
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
bool intersect_at(const LQVec<double>& ni, const LQVec<double>& pi,
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
bool intersect_at(const LQVec<double>& ni,
                  const LQVec<double>& nj,
                  const LQVec<double>& nk,
                  LQVec<double>& intersect, const int idx);
#endif
