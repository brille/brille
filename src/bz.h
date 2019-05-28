/*! \file */
#ifndef _BZ_CLASS_H_
#define _BZ_CLASS_H_

#include "neighbours.h"
#include "transform.h"

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
//  ArrayVector<double> position;
//  ArrayVector<int> type;
  ArrayVector<double> vertices;     //!< The coordinates of the first Brillouin zone vertices in an orthonormal frame
  ArrayVector<int> faces;           //!< The reciprocal lattice points defining the first Brillouin zone faces -- twice the actual face vectors
  ArrayVector<int> faces_per_vertex;//!< The indexes of the three faces providing the intersection for each vertex
public:
  /*!
  @param lat A Reciprocal lattice
  @param toprim A flag to indicate if the primtive lattice should be used
  @param extent An integer to control how-far the vertex-finding algorithm should
                search in τ-index. The default indicates that (̄1̄1̄1),
                (̄1̄10), (̄1̄11), (̄10̄1), ..., (111) are included.
  */
  BrillouinZone(Reciprocal lat, bool toprim=true, int extent=1): outerlattice(lat) {
    lattice = toprim ? lat.primitive() : lat;
    this->vertex_search(extent);
  }
  void set_vertices(ArrayVector<double> newverts);
  void set_faces(ArrayVector<int> newfaces);
  void set_faces_per_vertex(ArrayVector<int> newfpv);
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
  // void vertex_search_xyz(const int extent=1);
  //! Returns the number of vertices defining the first Brillouin zone
  size_t vertices_count() const { return vertices.size();};
  //! Returns the number of reciprocal lattice points defining the first Brillouin zone
  size_t faces_count() const { return faces.size();};
  //! Returns the lattice passed in at construction
  const Reciprocal get_lattice() const { return this->outerlattice;};
  //! Returns the lattice actually used to find the Brillouin zone vertices,
  //! which may be a primitive lattice depending on the flag at creation
  const Reciprocal get_primitive_lattice() const { return this->lattice;};
  //! Returns the vertices as LQVec objects expressed in the construction-input lattice
  LQVec<double>    get_vertices() const;
  //! Returns the defining reciprocal lattice points as LQVec objects expressed in the construction-input lattice
  LQVec<int>       get_faces() const ;
  //! Returns the vertices as LQVec objects expressed in the zone-defining lattice
  LQVec<double>    get_primitive_vertices() const;
  //! Returns the defining reciprocal lattice points as LQVec objects expressed in the zone-defining lattice
  LQVec<int>       get_primitive_faces() const ;
  /*! \brief Returns the indexes into the reciprocal lattice points array for
       the three which produced each intersection vertex

  @note In systems with more than three faces adjacent to a vertex the
        BrillouinZone object will find the correct vertices but will not return
        all neighbouring-face indexes for each vertex.
        This is problematic when attempting to invert the information returned here,
        e.g., when constructing 3-D patches for visualisation of the Brillouin zone
        facets. Thankfully, one can overcome this issue by recalling that the dot
        product of all vertices on a facet minus a point on the facet with the
        facet normal is zero -- `dot(vertices-faces/2,faces)==0`.
  */
  ArrayVector<int> get_faces_per_vertex() const;
  //! Print a representation of the object to the console
  void print() const;
  //! Return true if the lattice used to find the vertices is not the same as the one passed at construction, and is therefore primitive
  bool isprimitive(void) const {return !(lattice.issame(outerlattice));};

  /*! \brief Determine whether points are inside of the first Brillouin zone
    @param p A reference to a LQVec list of Q points to be checked
    @returns An ArrayVector<bool> with each 1-element array indicating if the
             associated Q point is inside of the Brillouin zone.
  */
  template<typename T> ArrayVector<bool> isinside(const LQVec<T>& p);
  /*! \brief Find q and τ such that Q=q+τ and τ is a reciprocal lattice vector
    @param[in] Q A reference to LQVec list of Q points
    @param[out] q The reduced reciprocal lattice vectors
    @param[out] tau The reciprocal lattice zone centres
  */
  bool moveinto(const LQVec<double>& Q, LQVec<double>& q, LQVec<int>& tau);
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
bool three_plane_intersection(const LQVec<double> *n,
                              const LQVec<double> *p,
                              const ArrayVector<double> *xyz,
                              const int i, const int j, const int k,
                              LQVec<double> *iat, const int idx=0);

#endif
