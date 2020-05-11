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

#ifndef _POLYHEDRON_H_
#define _POLYHEDRON_H_

#include <random>
#include <chrono>
#include <vector>
#include "arrayvector.hpp"
#include "latvec.hpp"
#include "debug.hpp"
#include "utilities.hpp"

template<typename T> static std::vector<T> unique(const std::vector<T>& x){
    std::vector<T> out;
    out.push_back(x[0]);
    for (auto & v: x) if (std::find(out.begin(), out.end(), v) == out.end())
        out.push_back(v);
    return out;
}
template<typename T> static bool is_dangling(const std::vector<size_t>& adjacents, const std::vector<T>& face){
  bool not_dangling{true};
  for (auto & x: face) not_dangling &= adjacents[x] > 2u;
  return !not_dangling;
}

template <typename T>
std::vector<T> bare_winding_angles(const ArrayVector<T>& vecs, const size_t i, const ArrayVector<T>& n){
  if (vecs.numel()!=3u)
    throw std::runtime_error("Finding a winding angle requires the cross product, which is only defined in 3D");
  // vecs should be normalized already
  std::vector<T> angles(vecs.size(), 0.);
  T dotij, y_len, angij;
  ArrayVector<T> x(3u,1u), y(3u,1u); // ensure all have memory allocated
  T crsij[3]; // to hold the cross product
  for (size_t j=0; j<vecs.size(); ++j) if (i!=j) {
    // if (j == i){
    //   angles[j] = 0.0;
    //   continue;
    // }
    dotij = vecs.dot(i,j);
    vector_cross(crsij, vecs.data(i), vecs.data(j));
    // crsij = vecs.cross(i,j);
    x = dotij * vecs.extract(i);
    y = vecs.extract(j) - x;
    y_len = y.norm(0) * (std::signbit(vector_dot(crsij, n.data(0))) ? -1 : 1);
    angij = std::atan2(y_len, dotij);
    angles[j] = angij < 0 ? angij+2*PI : angij;
  }
  return angles;
}
// template <typename T>
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
//     angles[j] = angij < 0 ? angij+2*PI : angij;
//   }
//   return angles;
// }
template<class T> static T triangle_area_squared(const T a, const T b, const T c){
  T s = (a+b+c)/2;
  return s*(s-a)*(s-b)*(s-c);
}
template<class T> int face_has_area(const ArrayVector<T>& points){
  // first verify that all points are coplanar
  // pick the first three points to define a plane, then ensure all points are in it
  if (points.size()<3) return -2; // can't be a face
  ArrayVector<T> p0 = points.extract(0);
  // ArrayVector<double> n = cross(points.extract(1)-p0, points.extract(2)-p0);
  // if (!dot(n, points-p0).all_approx(Comp::eq,0.)) return -1; // non-coplanar
  T s=0;
  ArrayVector<T> a,b;
  for (size_t i=1; i<points.size()-1; ++i){
    a = points.extract(i)-p0;
    for (size_t j=i+1; j<points.size(); ++j){
      b = points.extract(j)-p0;
      s += triangle_area_squared(a.norm(0), b.norm(0), (a-b).norm(0));
    }
  }
  return approx_scalar(s, T(0)) ? 0 : 1;
}

static inline bool ok_size(const size_t& a, const size_t& b){return (a==1u)||(a==b);}

template<class T> ArrayVector<bool> point_inside_plane(
  const ArrayVector<T>& n, const ArrayVector<T>& p, const ArrayVector<T>& x
){
  return dot(n, x-p).is_approx(Comp::le, 0.); // true if x is closer to the origin
}
template<class T> ArrayVector<bool> intersection(
  const ArrayVector<T>& ni, const ArrayVector<T>& pi,
  const ArrayVector<T>& nj, const ArrayVector<T>& pj,
  const ArrayVector<T>& nk, const ArrayVector<T>& pk,
  ArrayVector<T>& at
){
  size_t num[6]{ni.size(), pi.size(), nj.size(), pj.size(), nk.size(), pk.size()};
  size_t npt=0;
  for (int i=0; i<6; ++i) if (num[i]>npt) npt=num[i];
  for (int i=0; i<6; ++i) if (!ok_size(num[i],npt))
    throw std::runtime_error("All normals and points must be singular or equal sized");
  // output storage to indicate if an intersection exists
  ArrayVector<bool> out(1u, npt);
  // find the scaled intersection points
  // cross, multiply, and dot scale-up singleton inputs
  ArrayVector<T> tmp = cross(nj,nk)*dot(pi,ni) + cross(nk,ni)*dot(pj,nj) + cross(ni,nj)*dot(pk,nk);
  T detM, M[9];
  for (size_t i=0; i<npt; ++i){
    ni.get(i, M  );
    nj.get(i, M+3);
    nk.get(i, M+6);
    detM = matrix_determinant(M);
    if (std::abs(detM) > 1e-10){
      at.set(i, tmp.extract(i)/detM);
      out.insert(true, i);
    } else {
      out.insert(false, i);
    }
  }
  return out;
}
template<class T>
ArrayVector<bool> intersection(
  const ArrayVector<T>& n,  const ArrayVector<T>& p,
  const ArrayVector<T>& ni, const ArrayVector<T>& pi,
  const ArrayVector<T>& nj, const ArrayVector<T>& pj,
  const ArrayVector<T>& nk, const ArrayVector<T>& pk,
  ArrayVector<T>& at){
  ArrayVector<bool> ok = intersection(ni, pi, nj, pj, nk, pk, at);
  for (size_t i=0; i<ok.size(); ++i) if (ok.getvalue(i)) {
    ok.insert( point_inside_plane(n, p, at.extract(i)).all_true() , i);
  }
  return ok;
}
// template<typename... L> bool one_intersection(L... args){
//   return intersection(args...).getvalue(0);
// }
template<class T>
bool one_intersection(const ArrayVector<T>& n,  const ArrayVector<T>& p,
                      const ArrayVector<T>& ni, const ArrayVector<T>& pi,
                      const ArrayVector<T>& nj, const ArrayVector<T>& pj,
                      const ArrayVector<T>& nk, const ArrayVector<T>& pk,
                      ArrayVector<T>& at){
  return intersection(n,p,ni,pi,nj,pj,nk,pk,at).getvalue(0);
}

template<class T> std::vector<T> reverse(const std::vector<T>& x){
  std::vector<T> r;
  for (size_t i = x.size(); i--;) r.push_back(x[i]);
  return r;
}
template<class T> std::vector<std::vector<T>> reverse_each(const std::vector<std::vector<T>>& x){
  std::vector<std::vector<T>> r;
  std::transform(x.begin(), x.end(), std::back_inserter(r), [](const std::vector<T>& y){return reverse(y);});
  // for (auto i: x) r.push_back(reverse(i));
  return r;
}

class Polyhedron{
protected:
  ArrayVector<double> vertices;
  ArrayVector<double> points;
  ArrayVector<double> normals;
  std::vector<std::vector<int>> faces_per_vertex;
  std::vector<std::vector<int>> vertices_per_face;
public:
  // empty initializer
  Polyhedron(): vertices(ArrayVector<double>(3u, 0u)),
                points(ArrayVector<double>(3u, 0u)),
                normals(ArrayVector<double>(3u, 0u)),
                faces_per_vertex(std::vector<std::vector<int>>()),
                vertices_per_face(std::vector<std::vector<int>>()) {}
  //! Create a convex-hull Polyhedron from a set of points
  Polyhedron(const ArrayVector<double>& v): vertices(v){
    this->keep_unique_vertices();
    if (vertices.size() > 3){
      this->find_convex_hull();
      this->find_all_faces_per_vertex();
      this->polygon_vertices_per_face();
      this->purge_central_polygon_vertices();
      this->sort_polygons();
      this->purge_extra_vertices();
    }
  }
  //! Build a Polyhedron from vertices and vectors pointing to face centres
  Polyhedron(const ArrayVector<double>& v, const ArrayVector<double>& p):
  vertices(v), points(p), normals(p/norm(p)) {
    this->keep_unique_vertices();
    this->find_all_faces_per_vertex();
    this->polygon_vertices_per_face();
    // this->sort_polygons();
    this->purge_central_polygon_vertices();
    this->sort_polygons();
    this->purge_extra_vertices();
  }
  // initialize from vertices, points, and all relational information
  Polyhedron(const ArrayVector<double>& v,
             const ArrayVector<double>& p,
             const std::vector<std::vector<int>>& fpv,
             const std::vector<std::vector<int>>& vpf):
    vertices(v), points(p), normals(p/norm(p)), faces_per_vertex(fpv), vertices_per_face(vpf){
    this->sort_polygons();
  }
  // initalize from vertices, points, normals, and three-plane intersection information
  //! Build a Polyhedron from vertices, on-face points, and face normals
  Polyhedron(const ArrayVector<double>& v,
             const ArrayVector<double>& p,
             const ArrayVector<double>& n):
  vertices(v), points(p), normals(n) {
    verbose_update("Construct a polyhedron from vertices:\n",vertices.to_string());
    verbose_update("and planes (points, normals):\n", normals.to_string(points));
    this->keep_unique_vertices();
    this->find_all_faces_per_vertex();
    this->polygon_vertices_per_face();
    // this->sort_polygons();
    this->purge_central_polygon_vertices();
    this->sort_polygons();
    this->purge_extra_vertices();
  }
  // initialize from vertices, points, normals, and all relational information
  Polyhedron(const ArrayVector<double>& v,
             const ArrayVector<double>& p,
             const ArrayVector<double>& n,
             const std::vector<std::vector<int>>& fpv,
             const std::vector<std::vector<int>>& vpf):
    vertices(v), points(p), normals(n), faces_per_vertex(fpv), vertices_per_face(vpf){
      this->special_keep_unique_vertices(); // for + below, which doesn't check for vertex uniqueness
  }
  // initialize from vertices, and vertices_per_face
  Polyhedron(const ArrayVector<double>& v,
             const std::vector<std::vector<int>>& vpf):
  vertices(v), points({3u,0u}), normals({3u,0u}), vertices_per_face(vpf) {
    this->find_face_points_and_normals();
    this->sort_polygons();
    this->find_all_faces_per_vertex(); // do we really need this?
  }
  // initialize from vertices, point, normals, and vertices_per_face (which needs sorting)
  Polyhedron(
          const ArrayVector<double>& v,
          const ArrayVector<double>& p,
          const ArrayVector<double>& n,
          const std::vector<std::vector<int>>& vpf
          ):
          vertices(v), points(p), normals(n), vertices_per_face(vpf){
      this->sort_polygons();
      this->find_all_faces_per_vertex();
  }
  // copy constructor
  Polyhedron(const Polyhedron& other):
    vertices(other.get_vertices()),
    points(other.get_points()),
    normals(other.get_normals()),
    faces_per_vertex(other.get_faces_per_vertex()),
    vertices_per_face(other.get_vertices_per_face()) {}
  // assignment from another CentredPolyhedron
  Polyhedron& operator=(const Polyhedron& other){
    this->vertices = other.get_vertices();
    this->points = other.get_points();
    this->normals = other.get_normals();
    this->faces_per_vertex = other.get_faces_per_vertex();
    this->vertices_per_face = other.get_vertices_per_face();
    return *this;
  }
  Polyhedron mirror(void) const {
    return Polyhedron(-1*this->vertices, -1*this->points, -1*this->normals, this->faces_per_vertex, reverse_each(this->vertices_per_face));
    // return Polyhedron(-1*this->vertices, -1*this->points, -1*this->normals, this->faces_per_vertex, this->vertices_per_face);
  }
  // template<class T> Polyhedron rotate(const std::array<T,9> rot) const {
  //   ArrayVector<double> newv(3u, this->vertices.size()), newp(3u, this->points.size()), newn(3u, this->normals.size());
  //   for (size_t i=0; i<this->vertices.size(); ++i)
  //     multiply_matrix_vector<double,T,double>(newv.data(i), rot.data(), this->vertices.data(i));
  //   for (size_t i=0; i<this->points.size(); ++i)
  //     multiply_matrix_vector<double,T,double>(newp.data(i), rot.data(), this->points.data(i));
  //   for (size_t i=0; i<this->normals.size(); ++i)
  //     multiply_matrix_vector<double,T,double>(newn.data(i), rot.data(), this->normals.data(i));
  //   return Polyhedron(newv, newp, newn, this->faces_per_vertex, this->vertices_per_face);
  //   // return Polyhedron(newv, newp, newn);
  // }
  template<class T> Polyhedron rotate(const std::array<T,9> rot) const {
    ArrayVector<double> newv(3u, vertices.size());
    ArrayVector<double> newp(3u, points.size());
    ArrayVector<double> newn(3u, normals.size());
    for (size_t i=0; i<vertices.size(); ++i)
      multiply_matrix_vector<double,T,double>(newv.data(i), rot.data(), vertices.data(i));
    for (size_t i=0; i<points.size(); ++i)
      multiply_matrix_vector<double,T,double>(newp.data(i), rot.data(), points.data(i));
    for (size_t i=0; i<normals.size(); ++i)
      multiply_matrix_vector<double,T,double>(newn.data(i), rot.data(), normals.data(i));
    return Polyhedron(newv, newp, newn, this->faces_per_vertex, this->vertices_per_face);
  }
  Polyhedron operator+(const Polyhedron& other) const {
    size_t d = this->vertices.numel();
    if (other.vertices.numel() != d) throw std::runtime_error("Only equal dimensionality polyhedra can be combined.");
    // combine all vertices, points, and normals; adjusting the fpv and vpf indexing
    int tvn = static_cast<int>(this->vertices.size());
    int tfn = static_cast<int>(this->points.size());
    size_t ovn = other.vertices.size();
    size_t ofn = other.points.size();
    ArrayVector<double> v(d, tvn+ovn);
    ArrayVector<double> p(d, tfn+ofn), n(d, tfn+ofn);
    for (size_t i=0; i<signed_to_unsigned<size_t>(tvn); ++i) v.set(i,     this->vertices.extract(i));
    for (size_t i=0; i<ovn; ++i)                             v.set(tvn+i, other.vertices.extract(i));
    for (size_t i=0; i<signed_to_unsigned<size_t>(tfn); ++i) p.set(i,     this->points.extract(i));
    for (size_t i=0; i<ofn; ++i)                             p.set(tfn+i, other.points.extract(i));
    for (size_t i=0; i<signed_to_unsigned<size_t>(tfn); ++i) n.set(i,     this->normals.extract(i));
    for (size_t i=0; i<ofn; ++i)                             n.set(tfn+i, other.normals.extract(i));
    std::vector<std::vector<int>> fpv(this->faces_per_vertex), vpf(this->vertices_per_face);
    fpv.resize(tvn+ovn); vpf.resize(tfn+ofn);
    for (size_t i=0; i<ovn; ++i) for (auto j: other.faces_per_vertex[i]) fpv[tvn+i].push_back(tfn+j);
    for (size_t i=0; i<ofn; ++i) for (auto j: other.vertices_per_face[i]) vpf[tfn+i].push_back(tvn+j);
    return Polyhedron(v,p,n, fpv, vpf);
  }
  size_t num_vertices() const { return vertices.size(); }
  size_t num_faces() const { return normals.size(); }
  ArrayVector<double> get_vertices(void) const { return vertices; }
  ArrayVector<double> get_points(void) const { return points; }
  ArrayVector<double> get_normals(void) const { return normals; }
  std::vector<std::vector<int>> get_faces_per_vertex(void) const { return faces_per_vertex; }
  std::vector<std::vector<int>> get_vertices_per_face(void) const {return vertices_per_face; }
  //
  // const ArrayVector<double>& get_vertices(void) const { return vertices; }
  // const ArrayVector<double>& get_points(void) const { return points; }
  // const ArrayVector<double>& get_normals(void) const { return normals; }
  // const std::vector<std::vector<int>>& get_faces_per_vertex(void) const { return faces_per_vertex; }
  // const std::vector<std::vector<int>>& get_vertices_per_face(void) const {return vertices_per_face; }
  ArrayVector<double> get_half_edges(void) const{
    // for each face find the point halfway between each set of neighbouring vertices
    // Convex polyhedra always have edges neighbouring two faces, so we will
    // only ever find (∑pᵢ)>>1 half-edge points, where pᵢ are the number of
    // vertices (or edges) on the iᵗʰ face.
    size_t nfv = 0;
    for (auto f: vertices_per_face) nfv += f.size();
    // we don't care about point order, but we only want to have one copy
    // of each half-edge point. At the expense of memory, we can keep track
    // of which pairs we've already visited:
    std::vector<bool> unvisited(nfv*nfv, true);
    ArrayVector<double> hep(3u, nfv>>1);
    size_t found=0, a,b;
    for (auto f: vertices_per_face) for (size_t i=0; i<f.size(); ++i){
      a = f[i];
      b = f[(i+1)%f.size()]; // cyclic boundary condition on vertices
      // we visit the facet vertices in an arbitrary order, so we need to keep
      // track of our progress in [a,b] and [b,a]
      if (unvisited[a*nfv+b] && unvisited[b*nfv+a]){
        unvisited[a*nfv+b] = unvisited[b*nfv+a] = false;
        hep.set(found++, (vertices.extract(a)+vertices.extract(b))/2.0);
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
  std::string string_repr(void) const {
    size_t nv = vertices.size(), nf=points.size();
    std::string repr = "Polyhedron with ";
    repr += std::to_string(nv) + " " + (1==nv?"vertex":"vertices") + " and ";
    repr += std::to_string(nf) + " " + (1==nf?"facet":"facets");
    debug_exec(repr += "; volume " +std::to_string(this->get_volume());)
    return repr;
  }
  double get_volume(void) const {
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
    ArrayVector<double> a(3u, 1u), ba(3u, 1u), ca(3u, 1u);
    for (size_t f=0; f<normals.size(); ++f){
      a = this->vertices.extract(vertices_per_face[f][0]);
      for (size_t i=1; i<vertices_per_face[f].size()-1; ++i){ // loop over triangles
        ba = this->vertices.extract(vertices_per_face[f][ i ]) - a;
        ca = this->vertices.extract(vertices_per_face[f][i+1]) - a;
        vector_cross(n, ba.data(0), ca.data(0));
        subvol = vector_dot(a.data(0), n);
        // if (vector_dot(n, normals.data(f)) < 0) subvol *= -1.0;
        volume += subvol;
      }
    }
    return volume/6.0; // not-forgetting the factor of 1/6
  }
  ArrayVector<double> get_centroid(void) const {
    // also following http://wwwf.imperial.ac.uk/~rn/centroid.pdf
    double centroid[3]{0.,0.,0.};
    double n[3];
    ArrayVector<double> a(3u, 1u), b(3u, 1u), c(3u, 1u), ba(3u, 1u), ca(3u, 1u), bc(3u, 1u);
    for (auto verts: vertices_per_face){
      a = vertices.extract(verts[0]);
      for (size_t i=1; i<verts.size()-1; ++i){ // loop over triangles
        b = vertices.extract(verts[ i ]);
        c = vertices.extract(verts[i+1]);
        ba = b-a;
        ca = c-a;
        vector_cross(n, ba.data(0), ca.data(0));
        ba = a+b;
        bc = b+c;
        ca = c+a;
        for (int j=0; j<3; ++j){
          centroid[j] += n[j]*(ba.getvalue(0,j)*ba.getvalue(0,j) + bc.getvalue(0,j)*bc.getvalue(0,j) + ca.getvalue(0,j)*ca.getvalue(0,j));
        }
      }
    }
    double normalisation = 1.0/(48 * this->get_volume());
    for (int j=0; j<3; ++j) centroid[j] *= normalisation;
    return ArrayVector<double>(3u, 1u, centroid);
  }
protected:
  void keep_unique_vertices(void){
    std::vector<bool> flg;
    for (size_t i=0; i<vertices.size(); ++i) flg.push_back(true);
    int t = 3; // a tolerance multiplier tuning parameter, 3 seems to work OK.
    int n = static_cast<int>(vertices.numel());
    for (size_t i=1; i<vertices.size(); ++i) for (size_t j=0; j<i; ++j)
      if (flg[i]&&flg[j]) flg[i]=!approx_vector(n, vertices.data(i), vertices.data(j), t);
    this->vertices = this->vertices.extract(flg);
  }
    void special_keep_unique_vertices(){
        std::vector<bool> flg(vertices.size(), true);
        int t = 3; // a tolerance multiplier tuning parameter, 3 seems to work OK.
        int n = static_cast<int>(vertices.numel());
        for (size_t i=1; i<vertices.size(); ++i)
            for (size_t j=0; j<i; ++j)
                if (flg[i] && flg[j])
                    flg[i] = !approx_vector(n, vertices.data(i), vertices.data(j), t);
        if (std::count(flg.begin(), flg.end(), false)){
            // we need to correct vertices_per_face
            std::vector<int> map;
            int count{0};
            for (auto f: flg) map.push_back( f ? count++ : -1);
            // find the equivalent vertices for the non-unique ones:
            for (size_t i=0; i<vertices.size(); ++i)
                if (map[i]<0)
                    for (size_t j=0; j<vertices.size(); ++j)
                        if (i != j && map[j]>=0 && approx_vector(n, vertices.data(i), vertices.data(j), t))
                            map[i] = map[j];
            for (auto & face: vertices_per_face){
                std::vector<int> one;
                for (auto & f: face) one.push_back(map[f]);
                face = one;
            }
            // keep only the unique vertices
            this->vertices = this->vertices.extract(flg);
        }
    }
  void find_convex_hull(void){
    /* Find the set of planes which contain all vertices.
       The cross product between two vectors connecting three points defines a
       plane normal. If the plane passing through the three points partitions
       space into point-free and not-point-free then it is one of the planes of
       the convex-hull.                                                       */
    debug_update_if(vertices.size()<3, "find_convex_hull:: Not enough vertices for binomial_coefficient");
    unsigned long long bc = binomial_coefficient(vertices.size(), 3u);
    if (bc > static_cast<unsigned long long>(std::numeric_limits<size_t>::max()))
      throw std::runtime_error("Too many vertices to count all possible normals with a `size_t` integer");
    ArrayVector<double> n(3u, static_cast<size_t>(bc)), p(3u, static_cast<size_t>(bc));
    ArrayVector<double> ab(3u, 2u);
    size_t count = 0;
    // The same algorithm without temporary nijk, vi, vj, vk arrays (utilizing
    // vertices.extract, p.extract, and n.extract instead) does not work.
    // Either the compiler optimizes something important away or there is a bug
    // in repeatedly extracting the same array from an ArrayVector. In either
    // case, vectors which should not have been zero ended up as such.
    ArrayVector<double> nijk(3u,1u), vi(3u,1u), vj(3u,1u), vk(3u,1u);
    for (size_t i=0; i<vertices.size()-2; ++i){
      vi = vertices.extract(i);
      for (size_t j=i+1; j<vertices.size()-1; ++j){
        vj = vertices.extract(j);
        ab.set(0, vj-vi);
        for (size_t k=j+1; k<vertices.size(); ++k){
          vk = vertices.extract(k);
          ab.set(1, vk-vi);
          ab.cross(0, 1, nijk.data());
          // debug_update(i," ",j," ",k," ", ab.to_string(" x ")," = ",nijk.to_string(0));
          // increment the counter only if the new normal is not ⃗0 and partitions space
          if (!approx_scalar(nijk.norm(0),0.) && dot(nijk, vertices - vi).all_approx(Comp::le_ge, 0.)){
            // verify that the normal points the right way:
            if (dot(nijk, vertices - vi).all_approx(Comp::ge,0.))
              nijk = -1*nijk;
            // normalize the cross product to ensure we can determine uniqueness later
            n.set(count, nijk/nijk.norm(0));
            p.set(count, (vi+vj+vk)/3.0);
            ++count;
          }
        }
      }
    }
    if (n.size() < count)
      throw std::logic_error("Too many normal vectors found");
    // check that we only keep one copy of each unique normal vector
    n.resize(count); p.resize(count);
    ArrayVector<bool> nok = n.is_unique();
    normals = n.extract(nok);
    points = p.extract(nok);
  }
  // This is the function which is bad for non-convex polyhedra, since it
  // assigns all faces with n⋅(v-p)=0 to a vertex irrespective of whether two
  // faces have opposite direction normals.
  void find_all_faces_per_vertex(void){
    ArrayVector<double> vmp;
    std::vector<std::vector<int>> fpv(vertices.size());
    ArrayVector<bool> isonplane(1u, points.size());
    for (size_t i=0; i<vertices.size(); ++i){
      isonplane = dot(normals, vertices.extract(i) - points).is_approx(Comp::eq,0.);
      for (size_t j=0; j<points.size(); ++j) if (isonplane.getvalue(j)) fpv[i].push_back(static_cast<int>(j));
    }
    verbose_update("Found faces per vertex array\n",fpv);
    this->faces_per_vertex = fpv;
  }
  void polygon_vertices_per_face(void) {
    bool flag;
    // We have 3+ faces per vertex, so we can now find the vertices per face
    std::vector<std::vector<int>> vpf(this->points.size());
    for (size_t i=0; i<vertices.size(); ++i){
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
    for (auto facet: faces_per_vertex[i]) if (is_polygon[facet]) reduced_fpv[i].push_back(static_cast<int>(map[facet]));
    this->faces_per_vertex = reduced_fpv;
    verbose_update("Faces per vertex reduced to\n", faces_per_vertex);

    // plus cut-down the vertices_per_face vector
    std::vector<std::vector<int>> polygon_vpf;
    for (size_t i=0; i<vpf.size(); ++i) if (is_polygon[i]) polygon_vpf.push_back(vpf[i]);
    // for (auto i: vpf) if (i.size()>2) polygon_vpf.push_back(i);
    this->vertices_per_face = polygon_vpf;
    verbose_update("Vertices per (polygon) face\n", vertices_per_face);
  }
  void purge_central_polygon_vertices(void){
    /* We often build polyhedra from convex hulls of points and it is not
       uncommon for such point sets to include central face polygon points.
       Such points are not needed to describe the polyhedra and they inhibit
       sort_polygons from working, as the normalisation of face vertices causes
       a division by zero.
    */
    // Go through all faces an remove central vertices
    ArrayVector<double> facet_verts(3u, 0u), /*facet_normal,*/ facet_centre;
    ArrayVector<bool> is_centre_av(1u, 0u);
    std::vector<bool> is_centre;
    verbose_update("Starting vertices_per_face\n",vertices_per_face);
    for (size_t j=0; j<normals.size(); ++j){
      size_t facet_size = vertices_per_face[j].size();
      //facet_normal = normals.extract(j);
      facet_verts.resize(facet_size);
      for (size_t i=0; i<facet_size; ++i) facet_verts.set(i, vertices.extract(vertices_per_face[j][i]));
      facet_centre = sum(facet_verts)/static_cast<double>(facet_size);
      facet_verts -= facet_centre;
      is_centre_av = norm(facet_verts).is_approx(Comp::eq,0.);
      is_centre.resize(is_centre_av.size());
      for (size_t i=0; i<is_centre_av.size(); ++i) is_centre[i] = is_centre_av.getvalue(i);
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
    // debug_update("Pre-purge vertices_per_face\n",vertices_per_face);
    this->actual_vertex_purge();
    // debug_update("Purged vertices_per_face\n",vertices_per_face);
  }
  void actual_vertex_purge(void){
    // go through all faces again and determine whether a vertex is present
    ArrayVector<bool> keep(1u, vertices.size());
    for (size_t i=0; i<vertices.size(); ++i){
      bool flag = false;
      for (auto verts: vertices_per_face)
        if (!flag && std::find(verts.begin(), verts.end(), static_cast<int>(i))!=verts.end())
          flag = true;
      keep.insert(flag, i);
    }
    if (keep.count_true() < vertices.size()){
      verbose_update("Keeping ", keep.count_true(), " of ", vertices.size(), " vertices");
      // Remap the vertices_per_face array
      size_t count = 0, max=keep.count_true();
      std::vector<size_t> map;
      for (size_t i=0; i<keep.size(); ++i) map.push_back(keep.getvalue(i) ? count++ : max+1);
      //for (auto facet: vertices_per_face)
      //for (size_t i=0; i<facet.size(); ++i) if (map[facet[i]]<max) facet[i]=map[facet[i]];
      for (size_t j=0; j<vertices_per_face.size(); ++j)
        for (size_t i=0; i<vertices_per_face[j].size(); ++i)
          if (map[vertices_per_face[j][i]] < max)
            vertices_per_face[j][i] = map[vertices_per_face[j][i]];
      // Remove elements of faces_per_vertex and vertices[.extract(keep)]
      for (size_t i=keep.size(); i-- > 0; ) if (!keep.getvalue(i)) faces_per_vertex.erase(faces_per_vertex.begin()+i);
      vertices = vertices.extract(keep);
    }
  }
  void sort_polygons(void){
    std::vector<std::vector<int>> sorted_vpp(this->points.size());
    ArrayVector<double> facet_verts(3u, 0u), facet_centre, facet_normal;
    std::vector<int> facet, perm, used;
    std::vector<double> angles;
    double min_angle;
    size_t min_idx;
    ArrayVector<double> all_normals = this->get_normals();
    for (size_t j=0; j<this->points.size(); ++j){
      facet = this->vertices_per_face[j];
      verbose_update("Sorting face ",j," which has vertices",facet);
      facet_normal = all_normals.extract(j);
      facet_verts.resize(vertices_per_face[j].size());
      for (size_t i=0; i<vertices_per_face[j].size(); ++i) facet_verts.set(i, this->vertices.extract(facet[i]));
      facet_centre = sum(facet_verts)/static_cast<double>(facet.size());
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
        min_idx=facet.size()+1;
        for (size_t k=0; k<facet.size(); ++k)
          if ( std::find(used.begin(),used.end(),k)==used.end() // ensure the point hasn't already been picked
               && !approx_scalar(angles[k], 0.0) // that its not along the same line
               && angles[k] < min_angle // and that it has a smaller winding angle
             ){
            min_idx=k;
            min_angle = angles[k];
          }
        if (min_idx >= facet.size()){
          for (size_t d=0; d<facet.size(); ++d) debug_update("Facet vertex ",d," ", this->vertices.extract(facet[d]).to_string(0));
          debug_update("Facet centre", facet_centre.to_string(0));
          debug_update("Facet face vertices\n", facet_verts.to_string());
          debug_update("Winding angles ", angles);
          throw std::runtime_error("Error finding minimum winding angle polygon vertex");
        }
        perm[i] = min_idx;
        used.push_back(static_cast<int>(min_idx));
      }
      verbose_update("Producing sorting permutation",perm);
      for (size_t i=0; i<facet.size(); ++i) sorted_vpp[j].push_back(facet[perm[i]]); // this could be part of the preceeding loop.
    }
    this->vertices_per_face = sorted_vpp;
  }
  void purge_extra_vertices(void){
    /* If we used our convex hull algorithm to determine our polygon, it might
    have extraneous vertices within its convex polygonal faces */
    /* This method should be used only after find_all_faces_per_vertex,
      polygon_vertices_per_face, and sort_polygons have all been run as it
      assumes that the vertices per face are in increasing-winding-order.
    */
    ArrayVector<double> abc(3u, 3u);
    for (size_t n=0; n<normals.size(); ++n){
      verbose_update("A face with vertices ", vertices_per_face[n]);
      for (size_t i=0, j; vertices_per_face[n].size()>3 && i<vertices_per_face[n].size();){
        // pull out the vector from the previous point to this one
        j = (vertices_per_face[n].size()+i-1)%vertices_per_face[n].size();
        abc.set(0, vertices.extract(vertices_per_face[n][i]) - vertices.extract(vertices_per_face[n][j]));
        // pull out the vector from this point to the next one
        j = (vertices_per_face[n].size()+i+1)%vertices_per_face[n].size();
        abc.set(1, vertices.extract(vertices_per_face[n][j]) - vertices.extract(vertices_per_face[n][i]));
        // and find their dot product, storing the result in the same buffer
        abc.cross(0, 1, abc.data(2));
        // check whether the cross product points the same way as the face normal
        if (dot(normals.extract(n), abc.extract(2)).all_approx(Comp::gt, 0.)) {
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
  void find_face_points_and_normals(void){
    // if the Polyhedron is defined from its vertices and vertices_per_face,
    // then we require the input to be correct, and calculate points and normals
    this->points.resize(vertices_per_face.size());
    this->normals.resize(vertices_per_face.size());
    size_t count = 0;
    for (auto face: vertices_per_face){
      ArrayVector<double> centre(3u, 1u, 0.);
      for (size_t v: face) centre += vertices.extract(v);
      this->points.set(count, centre/static_cast<double>(face.size()));
      this->normals.set(count, cross(vertices.extract(face[1])-vertices.extract(face[0]), vertices.extract(face[2])-vertices.extract(face[1])));
      ++count;
    }
  }
public:
  Polyhedron centre(void) const {
    ArrayVector<double> centroid = this->get_centroid();
    return Polyhedron(vertices - centroid, points - centroid, normals, faces_per_vertex, vertices_per_face);
  }
  ArrayVector<bool> contains(const std::vector<std::array<double,3>>& x) const {
    return this->contains(ArrayVector<double>(x));
  }
  ArrayVector<bool> contains(const ArrayVector<double>& x) const {
    if (x.numel() != 3u) throw std::runtime_error("x must contain 3-vectors");
    ArrayVector<bool> out(1u, x.size(), false);
    for (size_t i=0; i<x.size(); ++i)
      out.insert(dot(this->normals, x.extract(i)-this->points).all_approx(Comp::le,0.), i);
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
  bool intersects(const Polyhedron& other) const {
    verbose_update("Checking intersection of ",this->string_repr()," with ",other.string_repr());
    ArrayVector<double> centroid = this->get_centroid();
    Polyhedron centred(vertices - centroid, points - centroid, normals, faces_per_vertex, vertices_per_face);
    Polyhedron ipoly = Polyhedron::bisect(centred, other.normals, other.points-centroid);
    double iv = ipoly.get_volume();
    // If two polyhedra intersect one another, their intersection is not null.
    return !approx_scalar(iv, 0.);
  }
  bool fuzzy_intersects(const Polyhedron& other) const {
    verbose_update("Checking intersection of ",this->string_repr()," with ",other.string_repr());
    const ArrayVector<double>& n = other.normals;
    const ArrayVector<double>& p = other.points;
    ArrayVector<double> at(3u,1u,0.);
    for (size_t i=0; i<n.size(); ++i){
      ArrayVector<double> crit = dot(n.extract(i), vertices-p.extract(i));
      if (crit.any_approx(Comp::gt,0.)){
        std::vector<size_t> del = find(crit.is_approx(Comp::gt,0.));
        std::vector<size_t> cut;
        for (auto x: del) for (auto f: faces_per_vertex[x])
        if (std::find(cut.begin(), cut.end(), f)==cut.end()) cut.push_back(f);
        verbose_update("Facets ",cut," are cut by the plane");
        // find the new intersection points of two neighbouring facets and the cutting plane
        for (size_t j=0; j<cut.size()-1; ++j) for (size_t k=j+1; k<cut.size(); ++k)
        // check if the three planes intersect, and give a point not-outside of the polyhedron
        if (one_intersection(normals, points, /* all polyhedron planes, to check for outsideness*/
                             n.extract(i), p.extract(i), /* the cutting plane */
                             normals.extract(cut[j]), points.extract(cut[j]), /* first cut plane */
                             normals.extract(cut[k]), points.extract(cut[k]), /* second cut plane */
                             at /* the intersection point, if it exists */ )
              && norm(vertices-at.extract(0)).none_approx(0.) /* we would only add a point if its new */
           ){
               return true;
        }
      }
    }
    return false;
  }
  Polyhedron intersection(const Polyhedron& other) const {
    // ArrayVector<double> centroid = this->get_centroid();
    // Polyhedron centred(vertices - centroid, points - centroid, normals, faces_per_vertex, vertices_per_face);
    // Polyhedron ipoly = Polyhedron::bisect(centred, other.normals, other.points-centroid);
    // return Polyhedron(ipoly.vertices + centroid, ipoly.points + centroid, ipoly.normals, ipoly.faces_per_vertex, ipoly.vertices_per_face);
    return Polyhedron::bisect(*this, other.normals, other.points);
  }
  bool intersects_fast(const Polyhedron& other) const {
    // check if any of our vertices are inside of the other polyhedron
    /* if the dot product is zero it means that a point is on the surface of the
       other polyhedron, which is fine. So we're using strictly less than zero. */
    for (size_t i=0; i<vertices.size(); ++i)
      if (dot(other.normals, vertices.extract(i)-other.points).any_approx(Comp::lt, 0.)){
        // for those of our vertices *in* the other polyhedron
        // ensure that they are not actually a shared vertex
        if (norm(other.vertices - vertices.extract(i)).none_approx(0.)) return true;
      }
    // check if any of the other vertices are inside of our polyhedron
    for (size_t i=0; i<other.vertices.size(); ++i)
      if (dot(normals, other.vertices.extract(i)-points).any_approx(Comp::lt, 0.))
        if (norm(this->vertices - other.vertices.extract(i)).none_approx(0.)) return true;
    // check for intersecting planes :(
    return this->intersects(other);
  }
  template<class T> Polyhedron divide(const ArrayVector<T>&n, const ArrayVector<T>& p){
    ArrayVector<double> centroid = this->get_centroid();
    Polyhedron centred(vertices-centroid, points-centroid, normals, faces_per_vertex, vertices_per_face);
    Polyhedron divided = Polyhedron::bisect(centred, n, p-centroid);
    return Polyhedron(divided.vertices+centroid, divided.points+centroid, divided.normals, divided.faces_per_vertex, divided.vertices_per_face);
  }

  /*! Find the polyhedron which results from slicing an existant polyhedron by
  one or more plane passing through its volume. The part closest to the origin
  is retained.*/
  template<class T> static Polyhedron bisect(const Polyhedron& pin, const ArrayVector<T>& n, const ArrayVector<T>& p) {
    if (n.numel()!=3u || p.numel()!=3u)
      throw std::runtime_error("Wrong number of vector elements");
    if (n.size()!=p.size())
      throw std::runtime_error("Equal number of vectors and points required");
    ArrayVector<bool> keep, valid;
    Polyhedron pout(pin);
    ArrayVector<double> pn, pp, pv=pout.get_vertices();
    std::vector<std::vector<int>> fpv, vpf;
    std::vector<int> del_vertices, cut, face_vertices, new_vector;
    std::vector<int> vertex_map;
    ArrayVector<double> at(3u, 1u), ni(3u,1u), pi(3u,1u), cen(3u,0u);
    // copy the current vertices, normals, and relational information
    pv=pout.get_vertices();
    pn=pout.get_normals();
    pp=pout.get_points();
    fpv = pout.get_faces_per_vertex();
    vpf = pout.get_vertices_per_face();
    // move the output polyhedron to be centred on the origin -- which requires we move all cutting planes as well
    // ArrayVector<double> origin = sum(pv)/static_cast<double>(pv.size());
    verbose_update("Cut a ",pout.string_repr()," by ",n.size()," planes");
    for (size_t i=0; i<n.size(); ++i){
      if (approx_scalar(pout.get_volume(), 0.)) break; // we can't do anything with an empty polyhedron
      verbose_update("Checking a plane"," passing through ",p.to_string(i)," normal to ",n.to_string(i));
      // check whether there's anything to do
      ni = n.extract(i);
      pi = p.extract(i);
      keep = point_inside_plane(ni, pi, pv);
      if (!keep.all_true()){
        verbose_update("Pre-cut ",i," polyhedron vertices:\n",pv.to_string());
        verbose_update("Pre-cut ",i," polyhedron planes (p,n):\n",pn.to_string(pp));
        verbose_update("Pre-cut ",i," polyhedron faces_per_vertex:\n", fpv);
        verbose_update("Pre-cut ",i," polyhedron vertices_per_face:\n",vpf);
        // compile the list of to-be-deleted vertices
        del_vertices.clear();
        for (size_t j=0; j<keep.size(); ++j) if(!keep.getvalue(j))
          del_vertices.push_back(static_cast<int>(j));
        verbose_update("Vertices beyond the cut ", del_vertices);
        // and the list of facets which need to be cut or removed
        cut.clear();
        for (size_t f=0; f<vpf.size(); ++f){
            bool tocut = false;
            for (auto x: del_vertices)
                if (std::find(vpf[f].begin(), vpf[f].end(), x) != vpf[f].end()) tocut = true;
            if (tocut) cut.push_back(static_cast<int>(f));
        }
        verbose_update("Facets ",cut," are cut by the plane");
        // find the new intersection points of two neighbouring facets and the cutting plane
        new_vector.clear();
        int last_face = static_cast<int>(pn.size()); // the to-be-index of the new facet (normal)
        unsigned new_face_vertex_count{0}, new_vertex_count{0};
        for (size_t j=0; j<cut.size()-1; ++j)
        for (size_t k=j+1; k<cut.size(); ++k)
        // check if the three planes intersect, and give a point not-outside of the polyhedron
        if ( at.size() > 0 &&
            one_intersection(pn, pp, /* all polyhedron planes, to check for outsideness*/
                             ni, pi, /* the cutting plane */
                             pn.extract(cut[j]), pp.extract(cut[j]), /* first cut plane */
                             pn.extract(cut[k]), pp.extract(cut[k]), /* second cut plane */
                             at /* the intersection point, if it exists */ )
           ){
          int lv;
          if (norm(pv-at.extract(0)).none_approx(0.)){  /* only add a point if its new */
            // grab the index of the next-added vertex
            lv = static_cast<int>(pv.size());
            // add the intersection point to all vertices
            verbose_update("adding the intersection point ",at.to_string(0)," to existing points\n",pv.to_string());
            pv = cat(pv, at.extract(0));
            // plus add the face index to the faces_per_vertex list for this new vertex
            fpv.push_back({cut[j], cut[k], last_face});
            // track how many new vertices we add
            ++new_vertex_count;
          } else {
            // find the matching index that already is in the list:
            lv = static_cast<int>(find(norm(pv-at.extract(0)).is_approx(Comp::eq,0.))[0]);
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
        debug_update_if(new_vector.size()!=new_face_vertex_count, "New vertices is wrong size!");
        if (new_vector.size()){
            // add the normal and point for the new face
            pn = cat(pn, ni);
            pp = cat(pp, pi);
            // extend the vertices per face vector
            vpf.push_back(new_vector);
        }
        // extend the keep ArrayVector<bool> to cover the new points
        keep = cat(keep, ArrayVector<bool>(1u, new_vertex_count, true));
        verbose_update("keep? vertices:\n",pv.to_string(keep));
        // find the new indices for all vertices, and extract the kept vertices
        vertex_map.resize(pv.size());
        int count{0};
        for (size_t j=0; j<pv.size(); ++j)
          vertex_map[j]= (keep.getvalue(j) ? count++ : -1);
        if (count == 0){
            // Without any kept vertices we do not have a Polyhedron.
            return Polyhedron();
        }
        pv = pv.extract(keep);
        verbose_update("vertex mapping:", vertex_map);
        // go through the vertices_per_face array, replacing vertex indicies
        // and making a face map to skip now-empty faces
        for (auto & face : vpf){
          new_vector.clear();
          for (auto x: face) if (keep.getvalue(x)) new_vector.push_back(vertex_map[x]);
          if (new_vector.empty()) new_vector.push_back(0); // avoid having an unallocated face. Removed later
          face = unique(new_vector); // remove duplicates (where do they come from?)
        }

        // we need to calculate the face centres
        for (size_t j=0; j<vpf.size(); ++j){
          cen.resize(vpf[j].size());
          for (size_t k=0; k<vpf[j].size(); ++k) cen.set(k, pv.extract(vpf[j][k]));
          // store the centroid as the on-plane point, but don't divide by zero
          if (cen.size()) pp.set(j, sum(cen)/static_cast<double>(cen.size()) );
        }
        // remove any faces without three vertices
        keep.resize(pp.size());\
        for (size_t j=0; j<pp.size(); ++j) keep.insert(unique(vpf[j]).size()>2, j);
        pp = pp.extract(keep);
        pn = pn.extract(keep);
        // and remove their empty vertex lists
        vpf.erase(std::remove_if(vpf.begin(), vpf.end(), [](std::vector<int> i){return unique(i).size()<3;}), vpf.end());
        // remove any faces that are not connected on all sides
        bool check_again{true};
        while (check_again){
            // find the number of adjacent faces for each vertes
            std::vector<size_t> adjacent_face_no(pv.size(), 0u);
            for (size_t j=0; j<pv.size(); ++j){
                for (auto & face: vpf)
                    if (std::find(face.begin(), face.end(), static_cast<int>(j))!=face.end())
                        adjacent_face_no[j]++;
            }
            // for each face now check if all vertices are adjacent to 3+ faces
            // storing the result for reducing pp & pn
            keep.resize(pp.size());
            size_t face_idx{0};
            for (auto & face: vpf)
              keep.insert(!is_dangling(adjacent_face_no, face), face_idx++);
            // and going through a second time to erase extraneous faces:
            vpf.erase(
              std::remove_if(vpf.begin(), vpf.end(),
                [adjacent_face_no](std::vector<int>& face){
                  return is_dangling(adjacent_face_no, face);
                }
              ), vpf.end());
            check_again = keep.count_true() < pp.size();
            if (check_again){
                pp = pp.extract(keep);
                pn = pn.extract(keep);
            }
        }

        // remove any vertices not on a face
        std::vector<bool> kv(pv.size(), false);
        for (auto & face: vpf) for (auto & x: face) kv[x] = true;
//        info_update("vertices per face\n",vpf,"\nkeeping vertices ",kv);
        if (std::count(kv.begin(), kv.end(), false)){
//            info_update("cut vertices\n",pv.to_string());
            std::vector<int> kv_map;
            int kv_count{0};
            for (auto && j : kv) kv_map.push_back(j ? kv_count++ : -1);
            for (auto & face: vpf){
                new_vector.clear();
                for (auto & x: face) new_vector.push_back(kv_map[x]);
                face = new_vector;
            }
            pv = pv.extract(kv);
//            info_update("become\n",pv.to_string());
        }

        verbose_update("to\n",vpf);
        verbose_update("keep? plane normals:\n",pn.to_string(keep));
        verbose_update("keep? plane points:\n",pp.to_string(keep));
        // use the Polyhedron intializer to sort out fpv and vpf -- really just fpv, vpf should be correct
        pout = Polyhedron(pv, pp, pn, vpf);
        // pout = Polyhedron(pv);
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
  ArrayVector<double> rand_rejection(const size_t n, const unsigned int seed=0) const {
    // initialize the random number generator with an optional non-zero seed:
    std::default_random_engine generator(seed>0 ? seed : std::chrono::system_clock::now().time_since_epoch().count());
    // construct the uniform distribution spanning [0,1]
    std::uniform_real_distribution<double> distribution(0.,1.);
    // find the minimum bounding corner and the vector to the maximum bounding corner
    ArrayVector<double> vmin = vertices.min(),  vdif = vertices.max() - vmin;
    // initialize the output points array
    ArrayVector<double> p(3, n);
    ArrayVector<double> r(3,1);
    // generate random points until we have `n` which are inside the polyhedron
    for (size_t i=0; i<n; ){
      // generate a vector between [0,0,0] and [1,1,1]
      for (size_t j=0; j<3u; ++j) r.insert(distribution(generator), 0,j);
      p.set(i, vmin + r*vdif);
      if (this->contains(p.extract(i)).getvalue(0)) ++i;
    }
    return p;
  }
};


template<class T>
Polyhedron polyhedron_box(std::array<T,3>& xmin, std::array<T,3>& xmax){
  ArrayVector<double> v(3u, 8u);
  v.insert(xmin[0], 0, 0); v.insert(xmin[1], 0, 1); v.insert(xmin[2], 0, 2); // 000 0
  v.insert(xmin[0], 1, 0); v.insert(xmax[1], 1, 1); v.insert(xmin[2], 1, 2); // 010 1
  v.insert(xmin[0], 2, 0); v.insert(xmax[1], 2, 1); v.insert(xmax[2], 2, 2); // 011 2
  v.insert(xmin[0], 3, 0); v.insert(xmin[1], 3, 1); v.insert(xmax[2], 3, 2); // 001 3
  v.insert(xmax[0], 4, 0); v.insert(xmin[1], 4, 1); v.insert(xmin[2], 4, 2); // 100 4
  v.insert(xmax[0], 5, 0); v.insert(xmax[1], 5, 1); v.insert(xmin[2], 5, 2); // 110 5
  v.insert(xmax[0], 6, 0); v.insert(xmax[1], 6, 1); v.insert(xmax[2], 6, 2); // 111 6
  v.insert(xmax[0], 7, 0); v.insert(xmin[1], 7, 1); v.insert(xmax[2], 7, 2); // 101 7
  std::vector<std::vector<int>> vpf{{3,0,4,7},{3,2,1,0},{0,1,5,4},{3,7,6,2},{7,4,5,6},{2,6,5,1}};
  return Polyhedron(v, vpf);
}

#endif // _POLYHEDRON_H_
