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

#ifndef BRILLE_POLYHEDRON_H_
#define BRILLE_POLYHEDRON_H_
/*! \file
    \author Greg Tucker
    \brief A class representing convex polyhedra in three dimensions
*/
#include <random>
#include <utility>
// #include <chrono>/
// #include <vector>
#include "array_latvec.hpp" // defines bArray
// #include "debug.hpp"
// #include "utilities.hpp"
// #include "approx.hpp"
namespace brille {
//! Find the unique elements elements of a vector
template<typename T> static std::vector<T> unique(const std::vector<T>& x){
    std::vector<T> out;
    out.push_back(x[0]);
    for (auto & v: x) if (std::find(out.begin(), out.end(), v) == out.end())
        out.push_back(v);
    return out;
}
//! Determine if a facet of a Polyhedron is part of a convex volume
template<typename T> static bool is_not_dangling(const std::vector<size_t>& adjacents, const std::vector<T>& face){
  bool not_dangling{true};
  for (auto & x: face) not_dangling &= adjacents[x] > 2u;
  return not_dangling;
}
//! Determine the winding angles for the vertices of a Polyhedron facet
template<class T>
std::vector<T> bare_winding_angles(const bArray<T>& vecs, const ind_t i, const bArray<T>& n){
  if (vecs.ndim() !=2 || vecs.size(1) != 3)
    throw std::runtime_error("Finding a winding angle requires the cross product, which is only defined in 3D");
  // vecs should be normalized already
  std::vector<T> angles(vecs.size(0), 0.);
  T dotij, y_len, angij;
  bArray<T> x(1u,3u), y(1u,3u); // ensure all have memory allocated
  T crsij[3]; // to hold the cross product
  for (ind_t j=0; j<vecs.size(0); ++j) if (i!=j) {
    dotij = vecs.dot(i,j);
    brille::utils::vector_cross(crsij, vecs.ptr(i), vecs.ptr(j));
    x = dotij * vecs.view(i);
    y = vecs.view(j) - x;
    y_len = y.norm(0) * (std::signbit(brille::utils::vector_dot(crsij, n.ptr(0))) ? -1 : 1);
    angij = std::atan2(y_len, dotij);
    angles[j] = angij < 0 ? angij+2*brille::pi : angij;
  }
  return angles;
}
// template<class T>
// std::vector<T> winding_angles(const LQVec<T>& vecs, const size_t i, const LQVec<T>& n){
//   // vecs should be normalized already
//   std::vector<T> angles(vecs.size());
//   T dotij, y_len, angij;
//   LQVec<T> crsij, x, y;
//   for (size_t j=0; j<vecs.size(); ++j){
//     if (j == i){
//       angles[j] = 0.0;
//       continue;
//     }
//     dotij = vecs.dot(i,j);
//     crsij = vecs.cross(i,j);
//     x = dotij * vecs.get(i);
//     y = vecs.get(j) - x;
//     y_len = y.norm(0) * (std::signbit(dot(crsij, n).getvalue(0)) ? -1 : 1);
//     angij = std::atan2(y_len, dotij);
//     angles[j] = angij < 0 ? angij+2*brille::pi : angij;
//   }
//   return angles;
// }
//! Determine the squared area of a triangle from the lengths of its sides
template<class T> static T triangle_area_squared(const T a, const T b, const T c){
  T s = (a+b+c)/2;
  return s*(s-a)*(s-b)*(s-c);
}
//! Check whether the points defining a Polyhedron facet encompass a finite area
template<class T>
int face_has_area(const bArray<T>& points){
  // first verify that all points are coplanar
  // pick the first three points to define a plane, then ensure all points are in it
  if (points.size(0)<3) return -2; // can't be a face
  auto p0 = points.view(0);
  T s=0;
  bArray<T> a,b;
  for (ind_t i=1; i<points.size(0)-1; ++i){
    a = points.view(i)-p0;
    for (ind_t j=i+1; j<points.size(0); ++j){
      b = points.view(j)-p0;
      s += triangle_area_squared(a.norm(0), b.norm(0), (a-b).norm(0));
    }
  }
  return brille::approx::scalar(s, T(0)) ? 0 : 1;
}
//! Check whether two sizes are compatible for broadcasting
static inline bool ok_size(const size_t& a, const size_t& b){return (a==1u)||(a==b);}

/*! \brief Check if a point is 'behind' a plane

\param n the plane normal
\param p a point in the plane
\param x the point to check
\returns whether `x` is on or behind the plane

\note A vector from any point on a plane to a test point has a positive dot
      product with the plane normal if the test point is 'outside' of the plane,
      zero if the test point is on the plane, and negative if it is 'behind'
      the plane.
*/
template<class T>
bArray<bool>
point_inside_plane(
  const bArray<T>& n, const bArray<T>& p, const bArray<T>& x
){
  return dot(n, x-p).is(brille::cmp::le, 0.); // true if x is closer to the origin
}
/*! \brief Find the intersection point of three planes, broadcasting over inputs

\param      n_i the normal to the first plane
\param      p_i a point on the first plane
\param      n_j the normal to the second plane
\param      p_j a point on the second plane
\param      n_k the normal to the third plane
\param      p_k a point on the third plane
\param[out] at  the intersection point of the three planes is set within, if it exists
\returns whether each intersection exists

\note If more than one plane is provided in any of the n_x, p_x combinations
      the intersection is found for all broadcast combinations of three planes.
      To successfully broadcast all i, j, and k planes must be singular or the
      same number. The output array will contain the intersection point for all
      broadcast combinations.
*/
template<class T>
std::vector<bool> intersection(
  const bArray<T>& n_i, const bArray<T>& p_i,
  const bArray<T>& n_j, const bArray<T>& p_j,
  const bArray<T>& n_k, const bArray<T>& p_k,
  bArray<T>& at
){
  using namespace brille;
  ind_t num[6]{n_i.size(0), p_i.size(0), n_j.size(0), p_j.size(0), n_k.size(0), p_k.size(0)};
  ind_t npt=0;
  for (ind_t i : num) if (i>npt) npt=i;
  for (ind_t i : num) if (!ok_size(i,npt))
    throw std::runtime_error("All normals and points must be singular or equal sized");
  // output storage to indicate if an intersection exists
  std::vector<bool> out(npt, false);
  // find the scaled intersection points
  // cross, multiply, and dot scale-up singleton inputs
  auto tmp = cross(n_j,n_k)*dot(p_i,n_i) + cross(n_k,n_i)*dot(p_j,n_j) + cross(n_i,n_j)*dot(p_k,n_k);
  T detM, mat[9];
  const T *ptr_i=n_i.ptr(0,0), *ptr_j=n_j.ptr(0,0), *ptr_k=n_k.ptr(0,0);
  bool u_i{n_i.size(0) == npt}, u_j{n_j.size(0) == npt}, u_k{n_k.size(0) == npt};
  for (ind_t a=0; a<npt; ++a){
    if (u_i) ptr_i = n_i.ptr(a,0);
    if (u_j) ptr_j = n_j.ptr(a,0);
    if (u_k) ptr_k = n_k.ptr(a,0);
    for (ind_t b=0; b<3; ++b){
      mat[b  ] = ptr_i[b];
      mat[b+3] = ptr_j[b];
      mat[b+6] = ptr_k[b];
    }
    detM = utils::matrix_determinant(mat);
    if (std::abs(detM) > 1e-10){
      at.set(a, tmp.view(a)/detM);
      out[a] = true;
    }
  }
  return out;
}
/*! \brief Find the intersection point of three planes within a set of planes, broadcasting over inputs

\param      n   the normals to the bounding plane(s)
\param      p   a point on each of the bounding plane(s)
\param      n_i the normal to the first plane
\param      p_i a point on the first plane
\param      n_j the normal to the second plane
\param      p_j a point on the second plane
\param      n_k the normal to the third plane
\param      p_k a point on the third plane
\param[out] at  the intersection point of the three planes is set within, if it exists
\returns whether each intersection point exists and is inside of the bounding
         plane(s).

\note If more than one plane is provided in any of the n_x, p_x combinations
      the intersection is found for all broadcast combinations of three planes.
      To successfully broadcast all i, j, and k planes must be singular or the
      same number. The output array will contain the intersection point for all
      broadcast combinations.
*/
template<class T>
std::vector<bool> intersection(
  const bArray<T>& n,  const bArray<T>& p,
  const bArray<T>& n_i, const bArray<T>& p_i,
  const bArray<T>& n_j, const bArray<T>& p_j,
  const bArray<T>& n_k, const bArray<T>& p_k,
  bArray<T>& at){
  std::vector<bool> ok = intersection(n_i, p_i, n_j, p_j, n_k, p_k, at);
  for (ind_t i=0; i<ok.size(); ++i) if (ok[i]){
    auto pip = point_inside_plane(n, p, at.view(i));
    ok[i] = pip.all();
  }
  return ok;
}
/*! \brief Find the intersection point of three planes within a set of planes

\param      n   the normals to the bounding plane(s)
\param      p   a point on each of the bounding plane(s)
\param      n_i the normal to the first plane
\param      p_i a point on the first plane
\param      n_j the normal to the second plane
\param      p_j a point on the second plane
\param      n_k the normal to the third plane
\param      p_k a point on the third plane
\returns the intersection point of the three planes if it exists and is inside
         of the bounding plane(s), otherwise an empty Array2.
*/
template<class T>
bArray<T>
one_intersection(
  const bArray<T>& n,  const bArray<T>& p,
  const bArray<T>& n_i, const bArray<T>& p_i,
  const bArray<T>& n_j, const bArray<T>& p_j,
  const bArray<T>& n_k, const bArray<T>& p_k)
{
  bArray<T> at(1u,3u);
  if (intersection(n,p,n_i,p_i,n_j,p_j,n_k,p_k,at)[0])
    return at;
  else
    return bArray<T>();
}
//! Return a new vector containing the reversed elements of a vector
template<class T> std::vector<T> reverse(const std::vector<T>& x){
  std::vector<T> r;
  for (size_t i = x.size(); i--;) r.push_back(x[i]);
  return r;
}
//! Return a new vector with each element-vector's elements reversed
template<class T> std::vector<std::vector<T>> reverse_each(const std::vector<std::vector<T>>& x){
  std::vector<std::vector<T>> r;
  std::transform(x.begin(), x.end(), std::back_inserter(r), [](const std::vector<T>& y){return reverse(y);});
  // for (auto i: x) r.push_back(reverse(i));
  return r;
}
/*! \brief A three-dimensional convex solid with polygonal facets

There is no mathematical restriction that a polygon be convex and concave
polygons are handled by other computational geometry libraries. These libraries
are much more complicated than is required here (as is supporting concave
polygons), so this class is a minimal implementation of the general polygon
tailored for `brille`.
*/
class Polyhedron{
public:
  using shape_t = typename bArray<double>::shape_t; //! The subscript index container of the vertices, points, and normals
protected:
  bArray<double> vertices;  //!< The vertices of the polyhedron
  bArray<double> points;    //!< A point on each facet of the polyhedron (extraneous, should be removed)
  bArray<double> normals;   //!< The normal vector of each facet (extraneous, but handy)
  std::vector<std::vector<int>> faces_per_vertex; //!< A listing for each vertex of the facets it contributes to (extraneous, should be removed)
  std::vector<std::vector<int>> vertices_per_face;//!< A listing of the vertices bounding each facet polygon
public:
    bool operator!=(const Polyhedron& other) const {
        if (vertices != other.vertices) return true;
        if (points != other.points) return true;
        if (normals != other.normals) return true;
        if (faces_per_vertex != other.faces_per_vertex) return true;
        if (vertices_per_face != other.vertices_per_face) return true;
        return false;
    }
    bool operator==(const Polyhedron& other) const { return !(this->operator!=(other)); }
  //! empty initializer
  explicit Polyhedron():
    vertices(bArray<double>()),
    points(bArray<double>()),
    normals(bArray<double>()),
    faces_per_vertex(std::vector<std::vector<int>>()),
    vertices_per_face(std::vector<std::vector<int>>())
  {}
  //! Create a convex-hull Polyhedron from a set of points
  explicit Polyhedron(const bArray<double>& v): vertices(v){
    this->keep_unique_vertices();
    if (vertices.size(0) > 3){
      this->find_convex_hull();
      this->find_all_faces_per_vertex();
      this->polygon_vertices_per_face();
      this->purge_central_polygon_vertices();
      this->sort_polygons();
      this->purge_extra_vertices();
    }
  }
  //! Build a Polyhedron from vertices and vectors pointing to face centres
  Polyhedron(const bArray<double>& v, const bArray<double>& p):
  vertices(v), points(p), normals(p/norm(p)) {
    this->keep_unique_vertices();
    this->find_all_faces_per_vertex();
    this->polygon_vertices_per_face();
    // this->sort_polygons();
    this->purge_central_polygon_vertices();
    this->sort_polygons();
    this->purge_extra_vertices();
  }
  //! initialize from vertices, points, and all relational information
  Polyhedron(const bArray<double>& v,
             const bArray<double>& p,
             std::vector<std::vector<int>> fpv,
             std::vector<std::vector<int>> vpf):
    vertices(v), points(p), normals(p/norm(p)), faces_per_vertex(std::move(fpv)), vertices_per_face(std::move(vpf)){
    this->sort_polygons();
  }
  // initalize from vertices, points, normals, and three-plane intersection information
  //! Build a Polyhedron from vertices, on-face points, and face normals
  Polyhedron(const bArray<double>& v,
             const bArray<double>& p,
             const bArray<double>& n):
  vertices(v), points(p), normals(n) {
    verbose_update("Construct a polyhedron from vertices:\n",vertices.to_string());
    verbose_update("and planes (points, normals):\n", normals.to_string(), points.to_string());
    this->keep_unique_vertices();
    this->find_all_faces_per_vertex();
    this->polygon_vertices_per_face();
    // this->sort_polygons();
    this->purge_central_polygon_vertices();
    this->sort_polygons();
    this->purge_extra_vertices();
  }
  //! initialize from vertices, points, normals, and all relational information
  Polyhedron(const bArray<double>& v,
             const bArray<double>& p,
             const bArray<double>& n,
             std::vector<std::vector<int>> fpv,
             std::vector<std::vector<int>> vpf):
    vertices(v), points(p), normals(n), faces_per_vertex(std::move(fpv)), vertices_per_face(std::move(vpf)){
      this->special_keep_unique_vertices(); // for + below, which doesn't check for vertex uniqueness
  }
  //! initialize from vertices, and vertices_per_face
  Polyhedron(const bArray<double>& v,
             std::vector<std::vector<int>> vpf):
  vertices(v), points(0,3), normals(0,3), vertices_per_face(std::move(vpf)) {
    this->find_face_points_and_normals();
    this->sort_polygons();
    this->find_all_faces_per_vertex(); // do we really need this?
  }
  //! initialize from vertices, point, normals, and vertices_per_face (which needs sorting)
  Polyhedron(
          const bArray<double>& v,
          const bArray<double>& p,
          const bArray<double>& n,
          std::vector<std::vector<int>> vpf
          ):
          vertices(v), points(p), normals(n), vertices_per_face(std::move(vpf)){
      this->sort_polygons();
      this->find_all_faces_per_vertex();
      verbose_update("Finished constructing Polyhedron with vertices\n",vertices.to_string(),"points\n",points.to_string(),"normals\n",normals.to_string());
  }
  //! copy constructor
  Polyhedron(const Polyhedron& other):
    vertices(other.get_vertices()),
    points(other.get_points()),
    normals(other.get_normals()),
    faces_per_vertex(other.get_faces_per_vertex()),
    vertices_per_face(other.get_vertices_per_face()) {}
  //! assignment from another CentredPolyhedron
  Polyhedron& operator=(const Polyhedron& other){
    this->vertices = other.get_vertices();
    this->points = other.get_points();
    this->normals = other.get_normals();
    this->faces_per_vertex = other.get_faces_per_vertex();
    this->vertices_per_face = other.get_vertices_per_face();
    return *this;
  }
  //! Return a space-inverted Polyhedron
  [[nodiscard]] Polyhedron mirror() const {
    return {-1*this->vertices, -1*this->points, -1*this->normals, this->faces_per_vertex, reverse_each(this->vertices_per_face)};

  }
  //! Apply the generalised rotation matrix to the returned Polyhedron
  template<class T> Polyhedron rotate(const std::array<T,9> rot) const {
    bArray<double> newv(vertices.size(0),3u);
    bArray<double> newp(points.size(0),  3u);
    bArray<double> newn(normals.size(0), 3u);
    for (ind_t i=0; i<vertices.size(0); ++i)
      brille::utils::multiply_matrix_vector<double,T,double>(newv.ptr(i), rot.data(), vertices.ptr(i));
    for (ind_t i=0; i<points.size(0); ++i)
      brille::utils::multiply_matrix_vector<double,T,double>(newp.ptr(i), rot.data(), points.ptr(i));
    for (ind_t i=0; i<normals.size(0); ++i)
      brille::utils::multiply_matrix_vector<double,T,double>(newn.ptr(i), rot.data(), normals.ptr(i));
    return {newv, newp, newn, this->faces_per_vertex, this->vertices_per_face};
  }
  //! Move the origin of the returned Polyhedron
  template<class T> [[nodiscard]] Polyhedron translate(const bArray<T>& vec) const {
    if (vec.numel() != 3)
      throw std::runtime_error("Translating a Polyhedron requires a single three-vector");
    // protect against + with an empty Array:
    if (0==vertices.size(0)||0==points.size(0)) return Polyhedron();
    return {vertices+vec, points+vec, normals, this->faces_per_vertex, this->vertices_per_face};
  }
  //! Combine two Polyhedron objects (not a true union)
  Polyhedron operator+(const Polyhedron& other) const {
    auto ndim = this->vertices.ndim();
    ind_t d = this->vertices.size(ndim-1);
    if (other.vertices.ndim() != ndim || other.vertices.size(ndim-1) != d)
      throw std::runtime_error("Only equal dimensionality polyhedra can be combined.");
    // combine all vertices, points, and normals; adjusting the fpv and vpf indexing
    ind_t tvn = this->vertices.size(0);
    ind_t tfn = this->points.size(0);
    ind_t ovn = other.vertices.size(0);
    ind_t ofn = other.points.size(0);
    bArray<double> v(tvn+ovn, d);
    bArray<double> p(tfn+ofn, d), n(tfn+ofn, d);
    for (ind_t i=0; i<tvn; ++i) v.set(i,    this->vertices.view(i));
    for (ind_t i=0; i<ovn; ++i) v.set(tvn+i,other.vertices.view(i));
    for (ind_t i=0; i<tfn; ++i) p.set(i,      this->points.view(i));
    for (ind_t i=0; i<ofn; ++i) p.set(tfn+i,  other.points.view(i));
    for (ind_t i=0; i<tfn; ++i) n.set(i,     this->normals.view(i));
    for (ind_t i=0; i<ofn; ++i) n.set(tfn+i, other.normals.view(i));
    std::vector<std::vector<int>> fpv(this->faces_per_vertex), vpf(this->vertices_per_face);
    fpv.resize(tvn+ovn); vpf.resize(tfn+ofn);
    for (ind_t i=0; i<ovn; ++i) for (auto j: other.faces_per_vertex[i])  fpv[tvn+i].push_back(static_cast<int>(tfn+j));
    for (ind_t i=0; i<ofn; ++i) for (auto j: other.vertices_per_face[i]) vpf[tfn+i].push_back(static_cast<int>(tvn+j));
    return {v,p,n, fpv, vpf};
  }
  //! Return the number of vertices in the Polyhedron
  [[nodiscard]] size_t num_vertices() const { return vertices.size(0); }
  //! Return the number of facets in the Polyhedron
  [[nodiscard]] size_t num_faces() const { return normals.size(0); }
  //! Return the vertex positions of the Polyhedron
  [[nodiscard]] bArray<double> get_vertices() const { return vertices; }
  //! Return the points on each facet of the Polyhedron
  [[nodiscard]] bArray<double> get_points() const { return points; }
  //! Return the normals of each facet of the Polyhedron
  [[nodiscard]] bArray<double> get_normals() const { return normals; }
  //! Return the facets per vertex
  [[nodiscard]] std::vector<std::vector<int>> get_faces_per_vertex() const { return faces_per_vertex; }
  //! Return the vertices per facet
  [[nodiscard]] std::vector<std::vector<int>> get_vertices_per_face() const {return vertices_per_face; }
  //! Return the facet polygons middle-edge points
  [[nodiscard]] bArray<double> get_half_edges() const{
    // for each face find the point halfway between each set of neighbouring vertices
    // Convex polyhedra always have edges neighbouring two faces, so we will
    // only ever find (∑pᵢ)>>1 half-edge points, where pᵢ are the number of
    // vertices (or edges) on the iᵗʰ face.
    ind_t nfv = 0;
    for (auto f: vertices_per_face) nfv += static_cast<ind_t>(f.size());
    // we don't care about point order, but we only want to have one copy
    // of each half-edge point. At the expense of memory, we can keep track
    // of which pairs we've already visited:
    std::vector<bool> unvisited(nfv*nfv, true);
    bArray<double> hep(nfv>>1, 3u);
    ind_t found=0, a,b;
    for (auto f: vertices_per_face) for (size_t i=0; i<f.size(); ++i){
      a = static_cast<ind_t>(f[i]);
      b = static_cast<ind_t>(f[(i+1)%f.size()]); // cyclic boundary condition on vertices
      // we visit the facet vertices in an arbitrary order, so we need to keep
      // track of our progress in [a,b] and [b,a]
      if (unvisited[a*nfv+b] && unvisited[b*nfv+a]){
        unvisited[a*nfv+b] = unvisited[b*nfv+a] = false;
        hep.set(found++, (vertices.view(a)+vertices.view(b))/2.0);
      }
    }
    if (found != nfv>>1){
      std::string msg = "Found " + std::to_string(found) + " half edge";
      msg += " points but expected to find " + std::to_string(nfv>>1);
      if (found < nfv>>1) msg += ". Is the polyhedron open? ";
      throw std::runtime_error(msg);
    }
    return hep;
  }
  //! Return a string representation of the Polyhedron
  [[nodiscard]] std::string string_repr() const {
    size_t nv = vertices.size(0), nf=points.size(0);
    std::string repr = "Polyhedron with ";
    repr += std::to_string(nv) + " " + (1==nv?"vertex":"vertices") + " and ";
    repr += std::to_string(nf) + " " + (1==nf?"facet":"facets");
    debug_exec(repr += "; volume " +std::to_string(this->get_volume());)
    return repr;
  }
  //! Return the volume of the Polyhedron
  [[nodiscard]] double get_volume() const {
    /* per, e.g., http://wwwf.imperial.ac.uk/~rn/centroid.pdf

    For a polyhedron with N triangular faces, each with ordered vertices
    (aᵢ, bᵢ, cᵢ), one can define nᵢ = (bᵢ-aᵢ)×(cᵢ-aᵢ) for each face and then
    find that the volume of the polyhedron is V = 1/6 ∑ᵢ₌₁ᴺ aᵢ⋅ nᵢ

    In our case here the polyhedron faces are likely not triangular, but we can
    subdivide each n-polygon face into n-2 triangles relatively easily.
    Furthermore, we can ensure that the vertex order is correct by comparing
    the triangle-normal to our already-stored facet normals.
    */
    double volume{0.}, subvol;
    double n[3];
    for (size_t f=0; f<normals.size(0); ++f){
      auto a = this->vertices.view(vertices_per_face[f][0]);
      for (size_t i=1; i<vertices_per_face[f].size()-1; ++i){ // loop over triangles
        auto ba = this->vertices.view(vertices_per_face[f][ i ]) - a;
        auto ca = this->vertices.view(vertices_per_face[f][i+1]) - a;
        brille::utils::vector_cross(n, ba.ptr(0), ca.ptr(0));
        subvol = brille::utils::vector_dot(a.ptr(0), n);
        // if (brille::utils::vector_dot(n, normals.data(f)) < 0) subvol *= -1.0;
        volume += subvol;
      }
    }
    return volume/6.0; // not-forgetting the factor of 1/6
  }
  //! Return the centre of mass of the Polyhedron (assuming uniform density)
  [[nodiscard]] bArray<double> get_centroid() const {
    // also following http://wwwf.imperial.ac.uk/~rn/centroid.pdf
    bArray<double> centroid({1u,3u}, 0.);
    double n[3];
    for (auto verts: vertices_per_face){
      auto a = vertices.view(verts[0]);
      for (size_t i=1; i<verts.size()-1; ++i){ // loop over triangles
        auto b = vertices.view(verts[ i ]);
        auto c = vertices.view(verts[i+1]);
        auto ba = b-a;
        auto ca = c-a;
        brille::utils::vector_cross(n, ba.ptr(0), ca.ptr(0));
        auto bc = b+c;
        ba = a+b;
        ca = c+a;
        for (ind_t j=0; j<3; ++j){
          centroid[j] += n[j]*(ba.val(0,j)*ba.val(0,j) + bc.val(0,j)*bc.val(0,j) + ca.val(0,j)*ca.val(0,j));
        }
      }
    }
    return centroid/(48 * this->get_volume());
  }
  //! Return the radius of the sphere encompassing the Polyhedron, centred on is centroid
  [[nodiscard]] double get_circumsphere_radius() const {
    auto centroid2vertices = vertices - this->get_centroid();
    auto c2v_lengths = norm(centroid2vertices);
    auto longest_c2v = c2v_lengths.max(0);
    double csphrad = longest_c2v[0];
    return csphrad;
    // return norm(vertices - this->get_centroid()).max(0)[0];
  }
protected:
  void keep_unique_vertices(){
    std::vector<bool> flg(vertices.size(0), true);
    for (ind_t i=1; i<vertices.size(0); ++i) for (ind_t j=0; j<i; ++j)
      if (flg[i] && flg[j]) flg[i] = !vertices.match(i,j);
    this->vertices = this->vertices.extract(flg);
  }
  void special_keep_unique_vertices(){
    std::vector<bool> flg(vertices.size(0), true);
    for (ind_t i=1; i<vertices.size(0); ++i) for (ind_t j=0; j<i; ++j)
      if (flg[i] && flg[j]) flg[i] = !vertices.match(i,j);
    if (std::count(flg.begin(), flg.end(), false)){
      // we need to correct vertices_per_face
      std::vector<int> map;
      int count{0};
      for (auto f: flg) map.push_back( f ? count++ : -1);
      // find the equivalent vertices for the non-unique ones:
      for (ind_t i=0; i<vertices.size(0); ++i) if (map[i]<0)
      for (ind_t j=0; j<vertices.size(0); ++j) if (i!=j && map[j]>=0 && vertices.match(i,j)) map[i] = map[j];
      for (auto & face: vertices_per_face){
          std::vector<int> one;
          for (auto & f: face) one.push_back(map[f]);
          face = one;
      }
      // keep only the unique vertices
      this->vertices = this->vertices.extract(flg);
    }
  }
  void find_convex_hull(){
    /* Find the set of planes which contain all vertices.
       The cross product between two vectors connecting three points defines a
       plane normal. If the plane passing through the three points partitions
       space into point-free and not-point-free then it is one of the planes of
       the convex-hull.                                                       */
    debug_update_if(vertices.size(0)<3, "find_convex_hull:: Not enough vertices for brille::utils::binomial_coefficient");
    unsigned long long bc = brille::utils::binomial_coefficient(vertices.size(0), 3u);
    if (bc > static_cast<unsigned long long>(std::numeric_limits<brille::ind_t>::max()))
      throw std::runtime_error("Too many vertices to count all possible normals with a `size_t` integer");
    bArray<double> n(static_cast<brille::ind_t>(bc), 3u);
    bArray<double> p(static_cast<brille::ind_t>(bc), 3u);
    bArray<double> ab(2u, 3u);
    ind_t count = 0;
    // The same algorithm without temporary nijk, vi, vj, vk arrays (utilizing
    // vertices.extract, p.extract, and n.extract instead) does not work.
    // Either the compiler optimizes something important away or there is a bug
    // in repeatedly extracting the same array from an Array. In either
    // case, vectors which should not have been zero ended up as such.
    for (ind_t i=0; i<vertices.size(0)-2; ++i){
      auto vi = vertices.view(i);
      for (ind_t j=i+1; j<vertices.size(0)-1; ++j){
        auto vji = vertices.view(j) - vi;
        for (ind_t k=j+1; k<vertices.size(0); ++k){
          auto vki = vertices.view(k) - vi;
          auto nijk = cross(vji, vki);
          // increment the counter only if the new normal is not ⃗0 and partitions space
          if (!brille::approx::scalar(nijk.norm(0),0.) && dot(nijk, vertices-vi).all(brille::cmp::le_ge, 0.)){
            // verify that the normal points the right way:
            if (dot(nijk, vertices-vi).all(brille::cmp::ge,0.)) nijk = -1*nijk;
            // normalize the cross product to ensure we can determine uniqueness later
            n.set(count, nijk/nijk.norm(0));
            p.set(count, vi+(vji+vki)/3.0); // the average of (vi+vj+vk)
            ++count;
          }
        }
      }
    }
    if (n.size(0) < count)
      throw std::logic_error("Too many normal vectors found");
    // check that we only keep one copy of each unique normal vector
    n.resize(count); p.resize(count);
    auto nok = n.is_unique();
    normals = n.extract(nok);
    points = p.extract(nok);
  }
  // This is the function which is bad for non-convex polyhedra, since it
  // assigns all faces with n⋅(v-p)=0 to a vertex irrespective of whether two
  // faces have opposite direction normals.
  void find_all_faces_per_vertex(){
    // std::vector<std::vector<int>> fpv(vertices.size(0));
    // for (size_t i=0; i<vertices.size(0); ++i){
    //   auto isonplane = dot(normals, vertices.view(i)-points).is(brille::cmp::eq,0.);
    //   for (size_t j=0; j<points.size(0); ++j) if (isonplane[j]) fpv[i].push_back(static_cast<int>(j));
    // }
    std::vector<std::vector<int>> fpv;
    for (ind_t i=0; i<vertices.size(0); ++i){
      std::vector<ind_t> onplane = dot(normals, vertices.view(i)-points).find(brille::cmp::eq,0.0);
      std::vector<int> vind;
      for (auto j: onplane) vind.push_back(static_cast<int>(j));
      fpv.push_back(vind);
    }
    verbose_update("Found faces per vertex array\n",fpv);
    this->faces_per_vertex = fpv;
  }

  void polygon_vertices_per_face() {
    bool flag;
    // We have 3+ faces per vertex, so we can now find the vertices per face
    std::vector<std::vector<int>> vpf(this->points.size(0));
    for (size_t i=0; i<vertices.size(0); ++i){
      for (auto facet: faces_per_vertex[i]){
        // if (std::find(vpf[facet].begin(), vpf[facet].end(), static_cast<int>(i))==vpf[facet].end()) vpf[facet].push_back(i);
        flag = true;
        for (auto vertex: vpf[facet]) if (static_cast<size_t>(vertex)==i) flag = false;
        if (flag) vpf[facet].push_back(static_cast<int>(i));
      }
    }
    verbose_update("Found vertices per face array\n", vpf);
    // additionally, we only want to keep faces which describe polygons
    std::vector<bool> is_polygon(vpf.size(), true);
    for (size_t i=0; i<vpf.size(); ++i)
      is_polygon[i] = face_has_area(vertices.extract(vpf[i])) > 0;
    verbose_update("Face is polygon:",is_polygon);

    this->points = this->points.extract(is_polygon);
    this->normals = this->normals.extract(is_polygon);

    // we should modify faces_per_vertex here, to ensure its indexing is correct
    size_t count = 0, max=vpf.size(); std::vector<size_t> map;
    for (size_t i=0; i<max; ++i) map.push_back(is_polygon[i] ? count++ : max);
    std::vector<std::vector<int>> reduced_fpv(faces_per_vertex.size());
    for (size_t i=0; i<faces_per_vertex.size(); ++i)
        for (auto facet: faces_per_vertex[i]) if (is_polygon[facet])
            reduced_fpv[i].push_back(static_cast<int>(map[facet]));
    this->faces_per_vertex = reduced_fpv;
    verbose_update("Faces per vertex reduced to\n", faces_per_vertex);

    // plus cut-down the vertices_per_face vector
    std::vector<std::vector<int>> polygon_vpf;
    for (size_t i=0; i<vpf.size(); ++i) if (is_polygon[i]) polygon_vpf.push_back(vpf[i]);
    // for (auto i: vpf) if (i.size()>2) polygon_vpf.push_back(i);
    this->vertices_per_face = polygon_vpf;
    verbose_update("Vertices per (polygon) face\n", vertices_per_face);
  }

  void purge_central_polygon_vertices(){
    /* We often build polyhedra from convex hulls of points and it is not
       uncommon for such point sets to include central face polygon points.
       Such points are not needed to describe the polyhedra and they inhibit
       sort_polygons from working, as the normalisation of face vertices causes
       a division by zero.
    */
    // Go through all faces an remove central vertices
    verbose_update("Starting vertices_per_face\n",vertices_per_face);
    for (size_t j=0; j<normals.size(0); ++j){
      //facet_normal = normals.extract(j);
      auto facet_verts = vertices.extract(vertices_per_face[j]);
      auto facet_centre = facet_verts.sum(0)/static_cast<double>(facet_verts.size(0));
      facet_verts -= facet_centre;
      auto is_centre = norm(facet_verts).is(brille::cmp::eq,0.).to_std(); // std::vector needed for erase
      verbose_update("Face ",j," has central vertices", is_centre);
      for (size_t i=0; i<vertices_per_face[j].size();){
        if (is_centre[i]){
          is_centre.erase(is_centre.begin()+i);
          vertices_per_face[j].erase(vertices_per_face[j].begin()+i);
        } else {
          ++i;
        }
      }
    }
    this->actual_vertex_purge();
  }
  void actual_vertex_purge(){
    // go through all faces again and determine whether a vertex is present
    std::vector<bool> keep(vertices.size(0), false);
    for (size_t i=0; i<keep.size(); ++i){
      for (auto v: vertices_per_face)
      if (std::find(v.begin(), v.end(), static_cast<int>(i)) != v.end()){
        keep[i] = true;
        break;
      }
    }
    size_t total = std::count(keep.begin(), keep.end(), true);
    if (total < vertices.size(0)){
      verbose_update("Keeping ", total, " of ", vertices.size(0), " vertices");
      // Remap the vertices_per_face array
      size_t count{0};
      std::vector<size_t> map;
      for (auto tf : keep) map.push_back(tf ? count++: total);
      for (auto& fv: vertices_per_face) for (auto& v: fv) if (map[v]<total) v=static_cast<int>(map[v]);
      // Remove elements of faces_per_vertex and vertices[.extract(keep)]
      for (size_t i=keep.size(); i-- > 0; ) if (!keep[i]) faces_per_vertex.erase(faces_per_vertex.begin()+i);
      vertices = vertices.extract(keep);
    }
  }
  void sort_polygons(){
    std::vector<std::vector<int>> sorted_vpp;
    std::vector<int> facet;
    std::vector<ind_t> perm, used;
    std::vector<double> angles;
    double min_angle;
    ind_t min_idx;
    auto all_normals = this->get_normals();
    for (ind_t j=0; j<this->points.size(0); ++j){
      facet = this->vertices_per_face[j];
      verbose_update("Sorting face ",j," which has vertices",facet);
      auto facet_normal = all_normals.view(j);
      auto facet_verts = vertices.extract(facet);
      auto facet_centre = facet_verts.sum(0)/static_cast<double>(facet.size());
      facet_verts -= facet_centre; // these are now on-face vectors to each vertex
      // if a point happens to be at the face centre dividing by the norm is a problem.
      facet_verts = facet_verts/norm(facet_verts); // and now just their directions;
      verbose_update("With on-plane vectors\n",facet_verts.to_string());
      perm.resize(facet.size());
      perm[0] = 0; // always start with whichever vertex is first
      used.clear();
      used.push_back(0);
      for (size_t i=1; i<facet.size(); ++i){
        angles = bare_winding_angles(facet_verts, perm[i-1], facet_normal);
        min_angle = 1e3;
        min_idx=static_cast<ind_t>(facet.size())+1;
        for (ind_t k=0; k<facet.size(); ++k)
          if ( std::find(used.begin(),used.end(),k)==used.end() // ensure the point hasn't already been picked
               && !brille::approx::scalar(angles[k], 0.0) // that its not along the same line
               && angles[k] < min_angle // and that it has a smaller winding angle
             ){
            min_idx=k;
            min_angle = angles[k];
          }
        if (min_idx >= facet.size()){
          std::string msg = "Error finding minimum winding angle polygon vertex\n";
          for (size_t d=0; d<facet.size(); ++d)
            msg += "Facet vertex " + std::to_string(d) + " " + this->vertices.view(facet[d]).to_string();
          msg += "Facet centre   " + facet_centre.to_string();
          msg += "Facet face vertices\n" + facet_verts.to_string();
          msg += "Winding angles [ ";
          for (auto ang: angles) msg += std::to_string(ang) + " ";
          msg += "]";
          throw std::runtime_error(msg);
        }
        perm[i] = min_idx;
        used.push_back(min_idx);
      }
      verbose_update("Producing sorting permutation",perm);
      std::vector<int> sorted_face_v;
      for (size_t i=0; i<facet.size(); ++i) sorted_face_v.push_back(facet[perm[i]]); // this could be part of the preceeding loop.
      // and add this face to the output collection
      sorted_vpp.push_back(sorted_face_v);
    }
    this->vertices_per_face = sorted_vpp;
  }
  void purge_extra_vertices(){
    /* If we used our convex hull algorithm to determine our polygon, it might
    have extraneous vertices within its convex polygonal faces */
    /* This method should be used only after find_all_faces_per_vertex,
      polygon_vertices_per_face, and sort_polygons have all been run as it
      assumes that the vertices per face are in increasing-winding-order.
    */
    for (ind_t n=0; n<normals.size(0); ++n){
      verbose_update("A face with vertices ", vertices_per_face[n]);
      for (size_t i=0, j; vertices_per_face[n].size()>3 && i<vertices_per_face[n].size();){
        // pull out the vector from the previous point to this one
        j = (vertices_per_face[n].size()+i-1)%vertices_per_face[n].size();
        auto prev = vertices.view(vertices_per_face[n][i]) - vertices.view(vertices_per_face[n][j]);
        // pull out the vector from this point to the next one
        j = (vertices_per_face[n].size()+i+1)%vertices_per_face[n].size();
        auto next = vertices.view(vertices_per_face[n][j]) - vertices.view(vertices_per_face[n][i]);

        // check whether the cross product points the same way as the face normal
        if (dot(normals.view(n), cross(prev, next)).all(brille::cmp::gt, 0.)) {
          // left-turn, keep this point
          ++i;
        } else {
          // right-turn, remove this point.
          vertices_per_face[n].erase(vertices_per_face[n].begin()+i);
        }
      }
      verbose_update("Is a convex polygonal face with vertices ", vertices_per_face[n]);
    }
    this->actual_vertex_purge();
  }
  void find_face_points_and_normals(){
    // if the Polyhedron is defined from its vertices and vertices_per_face,
    // then we require the input to be correct, and calculate points and normals
    this->points.resize(vertices_per_face.size());
    this->normals.resize(vertices_per_face.size());
    ind_t count = 0;
    for (const auto& face: vertices_per_face){
      // auto centre = vertices.extract(face).sum(0);
      // this->points.set(count, centre/static_cast<double>(face.size()));
      // this->normals.set(count, cross(vertices.view(face[1])-vertices.view(face[0]), vertices.view(face[2])-vertices.view(face[1])));
      auto fv = vertices.extract(face);
      auto centre = fv.sum(0)/static_cast<double>(face.size());
      this->points.set(count, centre);
      auto v01 = fv.view(1) - fv.view(0);
      auto v12 = fv.view(2) - fv.view(1);
      auto nrm = cross(v01, v12)/norm(v01)/norm(v12);
      this->normals.set(count, nrm);
      ++count;
    }
  }
public:
  //! Return a copy of this Polyhedron with its centre of mass at the origin
  [[nodiscard]] Polyhedron centre() const {
    auto centroid = this->get_centroid();
    return {vertices - centroid, points - centroid, normals, faces_per_vertex, vertices_per_face};
  }
  //! Determine if each of a set of points lies within the Polyhedron
  [[nodiscard]] std::vector<bool> contains(const std::vector<std::array<double,3>>& x) const {
    return this->contains(bArray<double>::from_std(x));
  }
  //! Determine if each of a set of points lies within the Polyhedron
  [[nodiscard]] std::vector<bool> contains(const bArray<double>& x) const {
    if (x.ndim()!=2 || x.size(x.ndim()-1)!=3) throw std::runtime_error("x must contain 3-vectors");
    std::vector<bool> out;
    for (ind_t i=0; i<x.size(0); ++i)
      out.push_back(dot(normals, x.view(i)-points).all(brille::cmp::le, 0.));
    return out;
  }
  /* Since we have the machinery to bisect a Polyhedron by a series of planes,
     the simplest way of checking for the intersection between two polyhedra is
     to bisect one of them by all of the planes of the other and if the
     resulting polyhedron has non-zero volume then the two polyhedra intersect.
     But this is almost certainly a far-from-optimal solution; especially since
     the bisect algorithm assumes that the input polyhedron contains (and will
     always contain) the origin.
  */
  //! Determine if the intersection of another Polyhedron with this one is non-null
  [[nodiscard]] bool intersects(const Polyhedron& other) const {
    verbose_update("Checking intersection of ",this->string_repr()," with ",other.string_repr());
    auto itrsct = Polyhedron::bisect(*this, other.normals, other.points);
    // If two polyhedra intersect one another, their intersection is not null.
    return !brille::approx::scalar(itrsct.get_volume(), 0.);
  }
  //! Return the intersection of another Polyhedron and this one
  [[nodiscard]] Polyhedron intersection(const Polyhedron& other) const {
    // auto centroid = this->get_centroid();
    // Polyhedron centred(vertices - centroid, points - centroid, normals, faces_per_vertex, vertices_per_face);
    // Polyhedron ipoly = Polyhedron::bisect(centred, other.normals, other.points-centroid);
    // return Polyhedron(ipoly.vertices + centroid, ipoly.points + centroid, ipoly.normals, ipoly.faces_per_vertex, ipoly.vertices_per_face);
    return Polyhedron::bisect(*this, other.normals, other.points);
  }
  //! Partition this Polyhedron by a plane and return the part 'behind' the plane
  template<class T>
  Polyhedron
  divide(const bArray<T>&n, const bArray<T>& p){
    bArray<double> centroid = this->get_centroid();
    Polyhedron centred(vertices-centroid, points-centroid, normals, faces_per_vertex, vertices_per_face);
    Polyhedron divided = Polyhedron::bisect(centred, n, p-centroid);
    return divided.translate(centroid);
  }

  /*! Find the polyhedron which results from slicing an existant polyhedron by
  one or more plane passing through its volume. The part closest to the origin
  is retained.*/
  template<class T>
  static Polyhedron
  bisect(const Polyhedron& pin, const bArray<T>& n, const bArray<T>& p) {
    assert(n.ndim()==2 && p.ndim()==2 && n.size(1)==3 && p.size(1)==(3) && n.size(0)==p.size(0));
    Polyhedron pout(pin);
    std::vector<int> vertex_map;
    // copy the current vertices, normals, and relational information
    auto pv  = pout.get_vertices();
    auto pn  = pout.get_normals();
    auto pp  = pout.get_points();
    auto fpv = pout.get_faces_per_vertex();
    auto vpf = pout.get_vertices_per_face();
    // move the output polyhedron to be centred on the origin -- which requires we move all cutting planes as well
    // auto origin = sum(pv)/static_cast<double>(pv.size(0));
    verbose_update("Cut a ",pout.string_repr()," by ",n.size(0)," planes");
    for (ind_t i=0; i<n.size(0); ++i){
      if (brille::approx::scalar(pout.get_volume(), 0.)) break; // we can't do anything with an empty polyhedron
      verbose_update("Cut with a plane passing through ",p.to_string(i)," normal to ",n.to_string(i));
      auto n_i = n.view(i);
      auto p_i = p.view(i);
      // check whether there's anything to do
      auto keep = point_inside_plane(n_i, p_i, pv).to_std(); // new keep each loop
      if (std::find(keep.begin(), keep.end(), false)!=keep.end()){
        verbose_update("Pre-cut ",i," polyhedron vertices:\n",pv.to_string());
        verbose_update("Pre-cut ",i," polyhedron planes (p,n):\n",cat(1,pn,pp).to_string());
        verbose_update("Pre-cut ",i," polyhedron faces_per_vertex:\n", fpv);
        verbose_update("Pre-cut ",i," polyhedron vertices_per_face:\n",vpf);
        // compile the list of to-be-deleted vertices
        std::vector<int> del_vertices;
        for (size_t j=0; j<keep.size(); ++j) if(!keep[j])
          del_vertices.push_back(static_cast<int>(j));
        verbose_update("Vertices beyond the cut ", del_vertices);
        // and the list of facets which need to be cut or removed
        std::vector<int> cut;
        for (size_t f=0; f<vpf.size(); ++f){
            bool tocut = false;
            for (auto x: del_vertices)
                if (std::find(vpf[f].begin(), vpf[f].end(), x) != vpf[f].end()) tocut = true;
            if (tocut) cut.push_back(static_cast<int>(f));
        }
        verbose_update("Facets ",cut," are cut by the plane");
        // find the new intersection points of two neighbouring facets and the cutting plane
        std::vector<int> new_vector;
        int last_face = static_cast<int>(pn.size(0)); // the to-be-index of the new facet (normal)
        unsigned new_face_vertex_count{0}, new_vertex_count{0};
        for (size_t j=0; j<cut.size()-1; ++j){
          auto n_j = pn.view(cut[j]);
          auto p_j = pp.view(cut[j]);
          for (size_t k=j+1; k<cut.size(); ++k){
            auto n_k = pn.view(cut[k]);
            auto p_k = pp.view(cut[k]);
            // check if the three planes intersect, and give a point not-outside of the polyhedron
            auto at = one_intersection(pn,pp, n_i,p_i, n_j,p_j, n_k,p_k);
            if ( at.size(0) == 1 ){ // there is a single intersection
              int lv;
              if (norm(pv-at).all(brille::cmp::neq, 0.)){  /* only add a point if its new */
                // grab the index of the next-added vertex
                lv = static_cast<int>(pv.size(0));
                verbose_update("Add intersection point ",at.to_string(0));//," to existing points\n",pv.to_string());
                pv.append(0, at);
                // plus add the face index to the faces_per_vertex list for this new vertex
                fpv.push_back({cut[j], cut[k], last_face});
                // track how many new vertices we add
                ++new_vertex_count;
              } else {
                // find the matching index that already is in the list:
                lv = static_cast<int>(norm(pv-at).first(brille::cmp::eq,0.));
                verbose_update("Reusing existing intersection point ",pv.to_string(lv)," for found intersection ",at.to_string(0));
              }
              // add the new vertex to the list for each existing facet -- if its not already present
              if (std::find(vpf[cut[j]].begin(), vpf[cut[j]].end(), lv)==vpf[cut[j]].end()) vpf[cut[j]].push_back(lv);
              if (std::find(vpf[cut[k]].begin(), vpf[cut[k]].end(), lv)==vpf[cut[k]].end()) vpf[cut[k]].push_back(lv);
              // and the yet-to-be-created facet
              new_vector.push_back(lv);
              // keeping track of how many points
              ++new_face_vertex_count;
            }
          }
        }
        debug_update_if(new_vector.size()!=new_face_vertex_count, "New vertices is wrong size!");
        if (!new_vector.empty()){
          // add the normal and point for the new face
          pn.append(0,n_i);
          pp.append(0,p_i);
          // extend the vertices per face vector
          vpf.push_back(new_vector);
        }
        // extend the keep std::vector<bool> to cover the new points
        for (size_t z=0; z<new_vertex_count; ++z) keep.push_back(true);
        verbose_update("keep:",keep,"vertices:\n", pv.to_string());
        // find the new indices for all vertices, and extract the kept vertices
        vertex_map.resize(pv.size(0));
        int count{0};
        for (size_t j=0; j<pv.size(0); ++j)
          vertex_map[j]= (keep[j] ? count++ : -1);
        if (count == 0){
            // Without any kept vertices we do not have a Polyhedron.
            return Polyhedron();
        }
        pv = pv.extract(keep);
        verbose_update("vertex mapping:", vertex_map);
        // go through the vertices_per_face array, replacing vertex indicies
        // and making a face map to skip now-empty faces
        for (auto & face : vpf){
          std::vector<int> nnv;
          for (auto x: face) if (keep[x]) nnv.push_back(vertex_map[x]);
          if (nnv.empty()) nnv.push_back(0); // avoid having an unallocated face. Removed later
          face = unique(nnv); // remove duplicates (where do they come from?)
        }

        // we need to calculate the face centres, but only for true faces (more than 2 vertices)
        for (ind_t j=0; j<vpf.size(); ++j) if(vpf[j].size()>2) {
          auto cen = pv.extract(vpf[j]);
          // store the centroid as the on-plane point
          // We need to be extra careful whenever we use Array.set since we may be sharing data with other Array(s)!!!
          pp = pp.decouple(); // decouple only copies memory if it is shared
          // if (cen.size(0)) // this check isn't necessary since vpf[j].size()!=0
          pp.set(j, cen.sum(0)/static_cast<double>(cen.size(0)) );
        }
        // remove any faces without three vertices
        std::vector<bool> has3v;
        for (size_t j=0; j<pp.size(0); ++j) has3v.push_back(unique(vpf[j]).size()>2);
        pp = pp.extract(has3v);
        pn = pn.extract(has3v);
        // and remove their empty vertex lists
        vpf.erase(std::remove_if(vpf.begin(), vpf.end(), [](std::vector<int> i){return unique(i).size()<3;}), vpf.end());
        // remove any faces that are not connected on all sides
        bool check_again{true};
        while (check_again){
            // find the number of adjacent faces for each vertes
            std::vector<size_t> adjacent_face_no(pv.size(0), 0u);
            for (size_t j=0; j<pv.size(0); ++j){
                for (auto & face: vpf)
                    if (std::find(face.begin(), face.end(), static_cast<int>(j))!=face.end())
                        adjacent_face_no[j]++;
            }
            // for each face now check if all vertices are adjacent to 3+ faces
            // storing the result for reducing pp & pn
            std::vector<bool> ok;
            for (auto & face: vpf)
              ok.push_back(is_not_dangling(adjacent_face_no, face));
            // and going through a second time to erase extraneous faces:
            vpf.erase(
              std::remove_if(vpf.begin(), vpf.end(),
                [adjacent_face_no](std::vector<int>& face){
                  return !is_not_dangling(adjacent_face_no, face);
                }
              ), vpf.end());
            check_again = std::find(ok.begin(), ok.end(), false) != ok.end();
            if (check_again){
                pp = pp.extract(ok);
                pn = pn.extract(ok);
            }
        }

        // remove any vertices not on a face
        std::vector<bool> kv(pv.size(0), false);
        for (auto & face: vpf) for (auto & x: face) kv[x] = true;
        if (std::count(kv.begin(), kv.end(), false)){
          std::vector<int> kv_map;
          int kv_count{0};
          for (auto && j : kv) kv_map.push_back(j ? kv_count++ : -1);
          for (auto & face: vpf){
            std::vector<int> nnv;
            for (auto & x: face) nnv.push_back(kv_map[x]);
            face = nnv;
          }
          pv = pv.extract(kv);
        }
        // with less than four points, normals, or vertices it can't be a Polyhedron
	      if (pv.size(0)<4 || pp.size(0)<4 || pn.size(0)<4){
          verbose_update("Null intersection");
          return Polyhedron();
        }
        // use the Polyhedron intializer to sort out fpv -- vpf should be correct
        pout = Polyhedron(pv, pp, pn, vpf);
        verbose_update("New ",pout.string_repr());
        // copy the updated vertices, normals, and relational information
        pv=pout.get_vertices();
        pn=pout.get_normals();
        pp=pout.get_points();
        fpv = pout.get_faces_per_vertex();
        vpf = pout.get_vertices_per_face();
      }
    }
    return pout;
  }
  //! Construst a list of random points within the volume of the Polyhedron via rejection
  [[nodiscard]] bArray<double> rand_rejection(const ind_t n, const unsigned int seed=0) const {
    // initialize the random number generator with an optional non-zero seed:
    auto tics = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed>0 ? seed : static_cast<unsigned int>(tics));
    // construct the uniform distribution spanning [0,1]
    std::uniform_real_distribution<double> distribution(0.,1.);
    // find the minimum bounding corner and the vector to the maximum bounding corner
    bArray<double> vmin = vertices.min(0),  vdif = vertices.max(0) - vmin;
    // initialize the output points array
    bArray<double> p(n,3);
    // generate random points until we have `n` which are inside the polyhedron
    for (ind_t i=0; i<n; ){
      // generate a in the box between vmin and vmin+vdif
      p.set(i, vmin + vdif * distribution(generator));
      // check if we can accept this
      if ( this->contains(p.view(i))[0] ) ++i;
    }
    return p;
  }
#ifdef USE_HIGHFIVE
public:
  template<class HF> std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, bool>
  to_hdf(HF& obj, const std::string& entry) const {
    auto group = overwrite_group(obj, entry);
    bool ok{true};
    ok &= vertices.to_hdf(group, "vertices");
    ok &= points.to_hdf(group, "points");
    ok &= normals.to_hdf(group, "normals");
    ok &= lists_to_hdf(faces_per_vertex, group, "faces_per_vertex");
    ok &= lists_to_hdf(vertices_per_face, group, "vertices_per_face");
    return ok;
  }
  [[nodiscard]] bool to_hdf(const std::string& filename, const std::string& dataset, const unsigned perm=HighFive::File::OpenOrCreate) const{
    HighFive::File file(filename, perm);
    return this->to_hdf(file, dataset);
  }
  template<class HF> static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, Polyhedron>
  from_hdf(HF& obj, const std::string& entry){
    auto group = obj.getGroup(entry);
    auto v = bArray<double>::from_hdf(group, "vertices");
    auto p = bArray<double>::from_hdf(group, "points");
    auto n = bArray<double>::from_hdf(group, "normals");
    auto fpv = lists_from_hdf<int>(group, "faces_per_vertex");
    auto vpf = lists_from_hdf<int>(group, "vertices_per_face");
    return {v, p, n, fpv, vpf};
  }
  static Polyhedron from_hdf(const std::string& filename, const std::string& dataset){
    HighFive::File file(filename, HighFive::File::ReadOnly);
    return Polyhedron::from_hdf(file, dataset);
  }
#endif //USE_HIGHFIVE
};

//! Construct a Polyhedron box from its minimal and maximal vertices
template<class T>
Polyhedron polyhedron_box(std::array<T,3>& xmin, std::array<T,3>& xmax){
  std::vector<std::array<T,3>> v{
    {xmin[0], xmin[1], xmin[2]}, // 000 0
    {xmin[0], xmax[1], xmin[2]}, // 010 1
    {xmin[0], xmax[1], xmax[2]}, // 011 2
    {xmin[0], xmin[1], xmax[2]}, // 001 3
    {xmax[0], xmin[1], xmin[2]}, // 100 4
    {xmax[0], xmax[1], xmin[2]}, // 110 5
    {xmax[0], xmax[1], xmax[2]}, // 111 6
    {xmax[0], xmin[1], xmax[2]}  // 101 7
  };
  std::vector<std::vector<int>> vpf{{3,0,4,7},{3,2,1,0},{0,1,5,4},{3,7,6,2},{7,4,5,6},{2,6,5,1}};
  return Polyhedron(bArray<double>::from_std(v), vpf);
}

} // end namespace brille
#endif // BRILLE_POLYHEDRON_H_
