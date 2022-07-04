/* This file is part of brille

Copyright © 2019-2022 Greg Tucker <gregory.tucker@ess.eu>

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
#include "tetgen.h"
#include "geometry.hpp"
namespace brille {
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
  std::vector<std::vector<int>> vertices_per_face;//!< A listing of the vertices bounding each facet polygon
public:
    bool operator!=(const Polyhedron& other) const {
      bool vertices_permuted{false};
      if (vertices != other.vertices) {
        vertices_permuted = vertices.is_permutation(other.vertices);
        if (!vertices_permuted) return true;
      }

      const auto & ovpf{other.vertices_per_face};

      if (vertices_per_face.size() != ovpf.size()) return true;

      auto add_to_no = [](const size_t & x, const std::vector<int>& v){return x + v.size();};
      auto vpf_no = std::accumulate(vertices_per_face.begin(), vertices_per_face.end(), 0u, add_to_no);
      auto ovpf_no = std::accumulate(ovpf.begin(), ovpf.end(), 0u, add_to_no);
      if (vpf_no != ovpf_no) return true;

      // now the hard work if we haven't failed yet
      std::vector<ind_t> permutation(vertices.size(0));
      if (!vertices_permuted){
        std::iota(permutation.begin(), permutation.end(), 0u);
      } else {
        // or do we need the opposite permutation vector? The first test used a permutation which was its own inverse
        permutation = vertices.permutation_vector(other.vertices);
      }

      std::vector<std::vector<int>> faces;
      faces.reserve(ovpf.size());
      auto permute = [p=permutation](const auto & i){return p[i];};
      for (const auto & face: ovpf){
        std::vector<int> one;
        one.reserve(face.size());
        std::transform(face.begin(), face.end(), std::back_inserter(one), permute);
        faces.push_back(one);
      }

      for (const auto & face: vertices_per_face){
        auto at = faces.end();
        for (auto itr=faces.begin(); itr != at; ++itr){
          if ((*itr).size() != face.size()) continue;
          for (size_t roll=0; roll < (*itr).size(); ++roll){
            std::vector<int> one((*itr).size());
            for (size_t i=0; i<one.size(); ++i) one[i] = (*itr)[(i + roll) % one.size()];
            if (std::equal(face.begin(), face.end(), one.begin())){
              at = itr;
              break;
            }
          }
          if (at == itr) break;
        }
        if (at != faces.end()){
          faces.erase(at); // remove it to avoid re-matching
        } else {
          // no equal face found, so these can't be the same
          return true;
        }
      }
      return false;
    }
    bool operator==(const Polyhedron& other) const { return !(this->operator!=(other)); }
  //! empty initializer
  explicit Polyhedron():
    vertices(bArray<double>()),
    vertices_per_face(std::vector<std::vector<int>>())
  {}
  //! Create a convex-hull Polyhedron from a set of points
  explicit Polyhedron(const bArray<double>& v): vertices(v){
    verbose_update("Construct convex hull Polyhedron from vertices:\n", vertices.to_string());
    debug_update("Convex Hull construction: \n", vertices.to_string());
    this->keep_unique_vertices();
    debug_update("Unique vertices \n", vertices.to_string());
    auto [a, b, c] = find_convex_hull_planes(vertices);
    debug_update("Bounding planes (a, b, c) \n", cat(1, a, b, c).to_string());
    // here's where we would double-check that the vertices are unique... but nothing has modified them

    auto fpv = find_planes_containing_point(a, b, c, vertices);
    debug_update("Faces per vertex\n", fpv);
    // find the vertices for each face, remove a central vertices if present, sort remaining face by to winding angle
    auto vpf = polygon_faces(a, b, c, fpv, vertices);
    debug_update("Gives vertices per face\n", vpf);
    // remove any vertices which are no longer part of a face -- takes the place of actual_vertex_purge
    if (polygon_face_vertex_purge(vertices, vpf)){
      debug_update("Some vertices purged. Now, vertices\n", vertices.to_string(), "vpf ", vpf);
    }
    // convert vpf from ind_t to int, storing it in vertices_per_face
    utils::convert_lists_type(vpf, vertices_per_face);
    verbose_update("Finished constructing convex hull Polyhedron with", this->python_string());
  }
  //! Build a Polyhedron from vertices and vectors pointing to face centres
  Polyhedron(const bArray<double>& v, const bArray<double>& p):
  vertices(v) {
//        FIXME
    this->keep_unique_vertices();
    auto [a, b, c] = plane_points_from_normal(p, p);
    auto fpv = find_planes_containing_point(a, b, c, vertices);
    auto vpf = polygon_faces(a, b, c, fpv, vertices);
    if (polygon_face_vertex_purge(vertices, vpf)){
      debug_update("Some vertices purged");
    }
    // convert vpf from ind_t to int, storing it in vertices_per_face
    utils::convert_lists_type(vpf, vertices_per_face);
    verbose_update("Finished constructing vertex, point, Polyhedron with", this->python_string());
  }

  //! Build a Polyhedron from vertices, on-face points, and face normals
  Polyhedron(const bArray<double>& v, const bArray<double>& p, const bArray<double>& n): vertices(v){
    verbose_update("Construct a polyhedron from vertices:\n",vertices.to_string());
    verbose_update("and planes (points, normals):\n", cat(1,p,n).to_string());
    this->keep_unique_vertices();
    auto [a, b, c] = plane_points_from_normal(n, p);
    auto fpv = find_planes_containing_point(a, b, c, vertices);
    auto vpf = polygon_faces(a, b, c, fpv, vertices);
    if (polygon_face_vertex_purge(vertices, vpf)){
      debug_update("Some vertices purged");
    }
    // convert vpf from ind_t to int, storing it in vertices_per_face
    utils::convert_lists_type(vpf, vertices_per_face);
    verbose_update("Finished constructing vertex, point, normal Polyhedron with", this->python_string());
  }

  //! initialize from vertices, and vertices_per_face (which must be correct)
  Polyhedron(const bArray<double>& v, std::vector<std::vector<int>> vpf): vertices(v), vertices_per_face(std::move(vpf)) {
    this->special_keep_unique_vertices(); // updates vertices_per_face too
    auto vpf_ind_t = polygon_faces(utils::convert_lists_type<ind_t>(vertices_per_face), vertices);
    if (polygon_face_vertex_purge(vertices, vpf_ind_t)){
      info_update("This shouldn't happen. Why have you provided extraneous vertices?!");
    }
    utils::convert_lists_type(vpf_ind_t, vertices_per_face);
    verbose_update("Finished constructing Polyhedron with vertices\n", this->python_string());
  }
  //! initialize from vertices, normals, and vertices_per_face (which needs sorting)
  Polyhedron(const bArray<double>& v,const bArray<double>& n,std::vector<std::vector<int>> vpf):
          vertices(v), vertices_per_face(std::move(vpf)){
      this->verify_vertices_per_face(n); // verifies the vertices_per_face is correct
      this->sort_polygons();
      verbose_update("Finished constructing Polyhedron with vertices\n", this->python_string());
  }
  //! copy constructor
  Polyhedron(const Polyhedron& other):
    vertices(other.get_vertices()),
    vertices_per_face(other.get_vertices_per_face()) {}
  //! assignment from another CentredPolyhedron
  Polyhedron& operator=(const Polyhedron& other){
    this->vertices = other.get_vertices();
    this->vertices_per_face = other.get_vertices_per_face();
    return *this;
  }
  //! Return a space-inverted Polyhedron
  [[nodiscard]] Polyhedron mirror() const {
    return {-1*this->vertices, reverse_each(this->vertices_per_face)};
  }
  //! Apply the generalised rotation matrix to the returned Polyhedron
  template<class T> Polyhedron rotate(const std::array<T,9> rot) const {
    bArray<double> newv(vertices.size(0),3u);
    for (ind_t i=0; i<vertices.size(0); ++i)
      brille::utils::multiply_matrix_vector<double,T,double>(newv.ptr(i), rot.data(), vertices.ptr(i));
    return {newv, this->vertices_per_face};
  }
  //! Move the origin of the returned Polyhedron
  template<class T> [[nodiscard]] Polyhedron translate(const bArray<T>& vec) const {
    if (vec.numel() != 3)
      throw std::runtime_error("Translating a Polyhedron requires a single three-vector");
    // protect against + with an empty Array:
    if (0==vertices.size(0)) return Polyhedron();
    return {vertices+vec, this->vertices_per_face};
  }
  //! Combine two Polyhedron objects (not a true union)
  Polyhedron operator+(const Polyhedron& other) const {
    auto ndim = this->vertices.ndim();
    ind_t d = this->vertices.size(ndim-1);
    if (other.vertices.ndim() != ndim || other.vertices.size(ndim-1) != d)
      throw std::runtime_error("Only equal dimensionality polyhedra can be combined.");
    // combine all vertices, points, and normals; adjusting the fpv and vpf indexing
    ind_t tvn = this->vertices.size(0);
    ind_t ovn = other.vertices.size(0);
    auto tfn = static_cast<ind_t>(this->num_faces());
    auto ofn = static_cast<ind_t>(other.num_faces());
    bArray<double> v(tvn+ovn, d);
    for (ind_t i=0; i<tvn; ++i) v.set(i,    this->vertices.view(i));
    for (ind_t i=0; i<ovn; ++i) v.set(tvn+i,other.vertices.view(i));
    std::vector<std::vector<int>> vpf(this->vertices_per_face);
    vpf.resize(this->num_faces() + other.num_faces());
    for (ind_t i=0; i<ofn; ++i) for (auto j: other.vertices_per_face[i]) vpf[tfn+i].push_back(static_cast<int>(tvn+j));

    // v *might* contain duplicate vertices at this point, which can not be allowed
    std::vector<bool> flg(v.size(0), true);
    for (ind_t i=1; i<v.size(0); ++i) for (ind_t j=0; j<i; ++j)
        if (flg[i] && flg[j]) flg[i] = !v.match(i,j);
    if (std::count(flg.begin(), flg.end(), false)){
      // we need to correct vertices_per_face
      std::vector<int> map;
      int count{0};
      map.reserve(flg.size());
      for (const auto & f: flg) map.push_back( f ? count++ : -1);
      // find the equivalent vertices for the non-unique ones:
      for (ind_t i=0; i<v.size(0); ++i) if (map[i]<0)
          for (ind_t j=0; j<v.size(0); ++j) if (i!=j && map[j]>=0 && v.match(i,j)) map[i] = map[j];
      // overwrite the vertices per face indexes with the (possibly) new values
      for (auto & f: vpf)
        std::transform(f.begin(), f.end(), f.begin(), [map](const int i){return map[i];});
      // keep only the unique vertices
      v = v.extract(flg);
    }
    return {v, vpf};
  }
  //! Return the number of vertices in the Polyhedron
  [[nodiscard]] size_t num_vertices() const { return vertices.size(0); }
  //! Return the number of facets in the Polyhedron
  [[nodiscard]] size_t num_faces() const { return this->vertices_per_face.size(); }
  //! Return the vertex positions of the Polyhedron
  [[nodiscard]] bArray<double> get_vertices() const { return vertices; }
  //! Return the points on each facet of the Polyhedron
  [[nodiscard]] bArray<double> get_points() const {
    bArray<double> p(this->num_faces(), 3u);
    size_t idx{0};
    for (const auto & face: this->vertices_per_face){
      auto fv = vertices.extract(face);
      p.set(idx++, fv.sum(0) / static_cast<double>(face.size()));
    }
    return p;
  }
  //! Return the normals of each facet of the Polyhedron
  [[nodiscard]] bArray<double> get_normals() const {
    bArray<double> n(this->num_faces(), 3u);
    size_t idx{0};
    for (const auto & face: vertices_per_face){
      n.set(idx++, three_point_normal(vertices, face));
    }
    return n;
  }
  //! Return the facets per vertex
  [[nodiscard]] std::vector<std::vector<int>> get_faces_per_vertex() const {
    return brille::utils::invert_lists(vertices_per_face);
  }
  //! Return the vertices per facet
  [[nodiscard]] std::vector<std::vector<int>> get_vertices_per_face() const {return vertices_per_face; }
  //! Return the facet indices required to define each plane
  [[nodiscard]] bArray<ind_t> get_facet_triplets() const {
    bArray<ind_t> triplets(vertices_per_face.size(), 3u);
    for (ind_t i=0; i<triplets.size(0); ++i){
      for (ind_t j=0; j<3u; ++j){
        triplets[{i, j}] = vertices_per_face[i][j];
      }
    }
    return triplets;
  }
  //! Return the facet polygons middle-edge points
  [[nodiscard]] bArray<double> get_half_edges() const{
    // for each face find the point halfway between each set of neighbouring vertices
    // Convex polyhedra always have edges neighbouring two faces, so we will
    // only ever find (∑pᵢ)>>1 half-edge points, where pᵢ are the number of
    // vertices (or edges) on the iᵗʰ face.
    ind_t nfv = 0;
    for (const auto & f: vertices_per_face) nfv += static_cast<ind_t>(f.size());
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
    size_t nv = vertices.size(0), nf=this->num_faces();
    std::string repr = "Polyhedron with ";
    repr += std::to_string(nv) + " " + (1==nv?"vertex":"vertices") + " and ";
    repr += std::to_string(nf) + " " + (1==nf?"facet":"facets");
    debug_exec(repr += "; volume " +std::to_string(this->get_volume());)
    return repr;
  }
  [[nodiscard]] std::string python_string() const {
      std::string s{"(np.array("};
      s += vertices.to_string() + "),\n" + lists_to_string(vertices_per_face) + ")";
      return s;
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
    for (size_t f=0; f<this->num_faces(); ++f){
      auto a = this->vertices.view(vertices_per_face[f][0]);
      for (size_t i=1; i<vertices_per_face[f].size()-1; ++i){ // loop over triangles
        auto ba = this->vertices.view(vertices_per_face[f][ i ]) - a;
        auto ca = this->vertices.view(vertices_per_face[f][i+1]) - a;
        brille::utils::vector_cross(n, ba.ptr(0), ca.ptr(0));
        subvol = brille::utils::vector_dot(a.ptr(0), n);
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

  void verify_vertices_per_face(const bArray<double>& normals){
    // get_normals calculates normal vectors from the vertex index order of each face
    auto n_dot_n = dot(normals, get_normals());
    // if the dot product is negative the vertex index order is wrong
    for (size_t i=0; i < vertices_per_face.size(); ++i){
      if (n_dot_n.val(i) < 0.){
        std::swap(vertices_per_face[i][0], vertices_per_face[i][1]);
      }
    }
    n_dot_n = dot(normals, get_normals());
    for (size_t i=0; i < vertices_per_face.size(); ++i){
      if (n_dot_n.val(i) < 0.){
        std::string msg = "Wrong normal direction for face " + std::to_string(i) + " with vertex indexes [ ";
        for (const auto & idx: vertices_per_face[i]) msg += std::to_string(idx) + " ";
        msg += "]";
        throw std::runtime_error(msg);
      }
    }
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
      // Remove elements of vertices[.extract(keep)]
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
    for (ind_t j=0; j<this->vertices_per_face.size(); ++j){
      facet = this->vertices_per_face[j];
//      verbose_update("Sorting face ",j," which has vertices",facet);
      auto facet_normal = all_normals.view(j);
      auto facet_verts = vertices.extract(facet);
      auto facet_centre = facet_verts.sum(0)/static_cast<double>(facet.size());
      facet_verts -= facet_centre; // these are now on-face vectors to each vertex
      // if a point happens to be at the face centre dividing by the norm is a problem.
      facet_verts = facet_verts/norm(facet_verts); // and now just their directions;
//      verbose_update("With on-plane vectors\n",facet_verts.to_string());
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
               && !brille::approx_float::scalar(angles[k], 0.0) // that its not along the same line
               && angles[k] < min_angle // and that it has a smaller winding angle
             ){
            min_idx=k;
            min_angle = angles[k];
          }
        if (min_idx >= facet.size()){
          std::string msg = "Error finding minimum winding angle polygon vertex\n";
          for (ind_t qq=0; qq < facet_verts.size(0); ++qq){
            for (ind_t ww=qq+1; ww < facet_verts.size(0); ++ww){
              auto val = (facet_verts.view(qq) - facet_verts.view(ww)).abs().sum();
              if (val < 1e-8) msg +=  "(" + std::to_string(qq) + "," + std::to_string(ww) + ") -> " + std::to_string(val * 1e15) + "x10^-15\n";
            }
          }
          for (size_t d=0; d<facet.size(); ++d)
            msg += "Facet vertex " + std::to_string(d) + " " + this->vertices.view(facet[d]).to_string();
          msg += "Facet centre   " + facet_centre.to_string();
          msg += "Facet face vertices\n" + facet_verts.to_string();
          msg += "Winding angles [ ";
          for (auto ang: angles) msg += std::to_string(ang) + " ";
          msg += "]";
          msg += "\n" + this->python_string();
          throw std::runtime_error(msg);
        }
        perm[i] = min_idx;
        used.push_back(min_idx);
      }
//      verbose_update("Producing sorting permutation",perm);
      std::vector<int> sorted_face_v;
      for (size_t i=0; i<facet.size(); ++i) sorted_face_v.push_back(facet[perm[i]]); // this could be part of the preceeding loop.
      // and add this face to the output collection
      sorted_vpp.push_back(sorted_face_v);
    }
    this->vertices_per_face = sorted_vpp;
  }
public:
  //! Return a copy of this Polyhedron with its centre of mass at the origin
  [[nodiscard]] Polyhedron centre() const {
    auto centroid = this->get_centroid();
    return {vertices - centroid, vertices_per_face};
  }
  //! Determine if each of a set of points lies within the Polyhedron
  [[nodiscard]] std::vector<bool> contains(const std::vector<std::array<double,3>>& x) const {
    return this->contains(bArray<double>::from_std(x));
  }
  //! Determine if each of a set of points lies within the Polyhedron
  [[nodiscard]] std::vector<bool> contains(const bArray<double>& x) const {
    if (x.ndim()!=2 || x.size(x.ndim()-1)!=3) throw std::runtime_error("x must contain 3-vectors");
    std::vector<bool> out;
    auto normals = this->get_normals();
    auto points = this->get_points();
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
    // If two polyhedra intersect one another, their intersection is not null.
    return !brille::approx_float::scalar(this->intersection(other).get_volume(), 0.);
  }
  //! Return the intersection of another Polyhedron and this one
  [[nodiscard]] Polyhedron intersection(const Polyhedron& other) const {
    auto v = other.get_vertices();
    std::vector<int> a, b, c;
    for (const auto & face: other.get_vertices_per_face()){
      a.push_back(face[0]);
      b.push_back(face[1]);
      c.push_back(face[2]);
    }
    return Polyhedron::bisect(*this, v.extract(a), v.extract(b), v.extract(c));
  }
  //! Partition this Polyhedron by a plane and return the part 'behind' the plane
  template<class T>
  Polyhedron
  divide(const bArray<T>& a, const bArray<T>& b, const bArray<T>& c){//(const bArray<T>&n, const bArray<T>& p){
    bArray<double> centroid = this->get_centroid();
    Polyhedron centred(vertices-centroid, vertices_per_face);
    Polyhedron divided = Polyhedron::bisect(centred, a - centroid, b - centroid, c - centroid);
    return divided.translate(centroid);
  }

  template<class T>
  [[nodiscard]] size_t matching_face_index(const bArray<T>& a, const bArray<T>& b, const bArray<T>& c) const {
    auto match{vertices_per_face.size()};
    for (size_t i=0; i < vertices_per_face.size(); ++i){
      const auto & face{vertices_per_face[i]};
      std::vector<double> o3d;
      o3d.reserve(face.size());
      for (const auto & vertex: face){
        o3d.push_back(orient3d(a.ptr(0), b.ptr(0), c.ptr(0), vertices.ptr(vertex)));
      }
      auto all_zero = std::all_of(o3d.begin(), o3d.end(), [](const auto & x){return approx_float::scalar(x, 0.);});
      if (all_zero){
        all_zero = dot(three_point_normal(vertices, face[0], face[1], face[2]), three_point_normal(a, b, c)).all(cmp::gt, 0.);
      }
      if (all_zero) match = i;
    }
    return match;
  }

  template<class T>
  [[nodiscard]] bool has_plane(const bArray<T>& a, const bArray<T>& b, const bArray<T>& c) const {
    return this->matching_face_index(a, b, c) < vertices_per_face.size();
  }

  template<class T>
  [[nodiscard]] bool none_beyond(const bArray<T>& a, const bArray<T>& b, const bArray<T>& c) const {
    auto match = matching_face_index(a, b, c);
    std::vector<int> v(vertices.size(0));
    std::iota(v.begin(), v.end(), 0);
    if (match < vertices_per_face.size()){
      auto itr = std::remove_if(v.begin(), v.end(), [f=vertices_per_face[match]](const auto & x){return std::find(f.begin(), f.end(), x) != f.end();});
      v.erase(itr, v.end());
    }
    auto is_negative = [&](const auto & i){
      return orient3d(a.ptr(0), b.ptr(0), c.ptr(0), vertices.ptr(i)) < 0;
    };
    return std::all_of(v.begin(), v.end(), is_negative);
  }

  /*! Find the polyhedron which results from slicing an existant polyhedron by
  one or more plane passing through its volume. The part closest to the origin
  is retained.*/
  template<class T>
  static Polyhedron
  bisect(const Polyhedron& pin, const bArray<T>& p_a, const bArray<T>& p_b, const bArray<T>& p_c) {
    assert(p_a.ndim()==2 && p_b.ndim()==2 && p_c.ndim() == 2 && p_a.size(1)==3 && p_b.size(1)==(3) && p_c.size(1) == 3);
    assert(p_a.size(0) == p_b.size(0) && p_b.size(0) == p_c.size(0));
    Polyhedron pout(pin);
    std::vector<int> vertex_map;
    // copy the current vertices, normals, and relational information
    auto pv  = pout.get_vertices();
    auto pn  = pout.get_normals();
    auto vpf = pout.get_vertices_per_face();
    // move the output polyhedron to be centred on the origin -- which requires we move all cutting planes as well
    // auto origin = sum(pv)/static_cast<double>(pv.size(0));
    debug_update("Cut a ",pout.string_repr()," by ",p_a.size(0)," planes");
    for (ind_t i=0; i<p_a.size(0); ++i){
      if (brille::approx_float::scalar(pout.get_volume(), 0.)) break; // we can't do anything with an empty polyhedron
      verbose_update("Cut with a plane passing through ",p_a.to_string(i),", ",p_b.to_string(i),", and ",p_c.to_string(i), " with normal ", three_point_normal(p_a, p_b, p_c).to_string(i));
      auto a_i = p_a.view(i);
      auto b_i = p_b.view(i);
      auto c_i = p_c.view(i);
      // check if there are any points beyond this plane, discounting any vertices which are part of a co-planar face
      if (pout.none_beyond(a_i, b_i, c_i)) continue;
      // check whether there's anything to do
      auto keep = point_inside_plane(a_i, b_i, c_i, pv);
      if (std::find(keep.begin(), keep.end(), false)!=keep.end()){
        verbose_update("Pre-cut ",i," polyhedron vertices:\n",pv.to_string());
        verbose_update("Pre-cut ",i," polyhedron plane normals:\n",pn.to_string());
        verbose_update("Pre-cut ",i," polyhedron vertices_per_face:\n",vpf);
        // compile the list of facets which need to be cut or removed
        auto cut = keep_to_cut_list(keep, vpf);
        verbose_update("Facets ",cut," are cut by the plane");
        // find the new intersection points of two neighbouring facets and the cutting plane
        std::vector<int> new_vector;
//        int last_face = static_cast<int>(pn.size(0)); // the to-be-index of the new facet (normal)
        unsigned new_face_vertex_count{0}, new_vertex_count{0};
        for (size_t j=0; j<cut.size()-1; ++j){ for (size_t k=j+1; k<cut.size(); ++k){
          // check if the facets-edge and plane intersect
          auto at = edge_plane_intersection(pv, vpf[cut[j]], vpf[cut[k]], a_i, b_i, c_i);
          if ( at.size(0) == 1 ){ // there is a single intersection
            int lv;
//            verbose_update("Cut facets ", cut[j]," and ", cut[k]);
            auto index = pv.first(brille::cmp::eq, at); // == pv.size(0) if at is not in pv
            if (index < pv.size(0)){
              verbose_update("Cut facets ", cut[j]," and ", cut[k], " use vertex at ",pv.view(index).to_string(),
                             "             for intersection at ",at.to_string(0));
              if (index >= pv.size(0) - new_vertex_count) {
                // this is a vertex added during this cut -- we want to keep its average?
                info_update("Figure out this whole averaging thing");
              }
              keep[index] = true; // ensure that we don't accidentally remove a needed vertex

            } else {
              verbose_update("Cut facets ", cut[j]," and ", cut[k], " add intersection at ",at.to_string(0));
              pv.append(0, at);
              ++new_vertex_count;
            }
            lv = static_cast<int>(index); // why is lv an int again?

            // add the new vertex to the list for each existing facet -- if it's not already present
            brille::utils::add_if_missing(vpf[cut[j]], lv);
            brille::utils::add_if_missing(vpf[cut[k]], lv);
            // and the yet-to-be-created facet
            new_vector.push_back(lv);
            // keeping track of how many points
            ++new_face_vertex_count;
          }
        }} // end of loops j and k
        debug_update_if(new_vector.size()!=new_face_vertex_count, "New vertices is wrong size!");
        if (!new_vector.empty() && !utils::unordered_list_in_lists(new_vector, vpf)){
          // add the normal for the new face
          pn.append(0,three_point_normal(a_i, b_i, c_i));
          // extend the vertices per face vector
          vpf.push_back(new_vector);
        }
        // extend the keep std::vector<bool> to cover the new points
        for (size_t z=0; z<new_vertex_count; ++z) keep.push_back(true);
//        verbose_update("keep:",keep,"vertices:\n", pv.to_string());
        // find the new indices for all vertices, and extract the kept vertices
        vertex_map.resize(pv.size(0));
        int count{0};
        for (size_t j=0; j<pv.size(0); ++j)
          vertex_map[j]= (keep[j] ? count++ : -1);
        if (count == 0){
            // Without any kept vertices we do not have a Polyhedron.
            return Polyhedron();
        }
        auto new_pv = pv.extract(keep);
//        verbose_update("vertex mapping:", vertex_map);
        // go through the vertices_per_face array, replacing vertex indicies
        // and making a face map to skip now-empty faces
//        verbose_update("vertices per face:\n", vpf);
        for (auto & face : vpf){
          std::vector<int> nnv;
          for (auto x: face) if (keep[x]) nnv.push_back(vertex_map[x]);
          if (nnv.empty()) nnv.push_back(0); // avoid having an unallocated face. Removed later
          nnv = unique(nnv); // remove duplicates (where do they come from?)
          if (nnv.size() > 2 && !(keep[face[0]] && keep[face[1]] && keep[face[2]])){
            if (dot(three_point_normal(pv, face), three_point_normal(new_pv, nnv)).all(brille::cmp::le, 0.))
              std::swap(nnv[1], nnv[2]);
          }
          face = nnv;
        }
        pv = new_pv; // done with normals, so we can switch to the reduced vertices
//        verbose_update("gets reduced to:\n", vpf);

        // remove any faces without three vertices
        std::vector<bool> has3v;
        has3v.reserve(vpf.size());
        for (const auto & j : vpf) has3v.push_back(unique(j).size()>2);
//        verbose_update("Keep faces with 3+ vertices: ", has3v);
        pn = pn.extract(has3v);
        // and remove their empty vertex lists
        vpf.erase(std::remove_if(vpf.begin(), vpf.end(), [](const std::vector<int> & i){return unique(i).size()<3;}), vpf.end());
//        verbose_update("So that vertices per face becomes\n", vpf);

        // remove any faces that are not connected on all sides
        bool check_again{true};
        while (check_again){
          // find the number of adjacent faces for each vertex
          std::vector<size_t> adjacent_face_no(pv.size(0), 0u);
          for (size_t j=0; j<pv.size(0); ++j){
              for (auto & face: vpf)
                  if (std::find(face.begin(), face.end(), static_cast<int>(j))!=face.end())
                      adjacent_face_no[j]++;
          }
          // for each face now check if all vertices are adjacent to 3+ faces
          // storing the result for reducing pn
          auto check = [adjacent_face_no](const std::vector<int> & face) -> bool {return !is_not_dangling(adjacent_face_no, face);};
          std::vector<bool> not_ok;
          std::transform(vpf.begin(), vpf.end(), std::back_inserter(not_ok), check);
          check_again = std::find(not_ok.begin(), not_ok.end(), true) != not_ok.end();
          // and going through a second time to erase extraneous faces:
          verbose_update_if(check_again, "Removing faces which are dangling");
          vpf.erase(std::remove_if(vpf.begin(), vpf.end(), check), vpf.end());
          verbose_update_if(check_again, "Vertices per face now\n", vpf);
          if (check_again){
            std::transform(not_ok.begin(), not_ok.end(), not_ok.begin(), [](const auto & b){return !b;});
            pn = pn.extract(not_ok);
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
          verbose_update("Keep face-full vertices ", kv, " of \n",pv.to_string());
          pv = pv.extract(kv);
          verbose_update("Resulting in \n", pv.to_string());
        }
        // with less than four normals, or vertices it can't be a Polyhedron
        if (pv.size(0)<4 || pn.size(0)<4){
          verbose_update("Null intersection since we have ",pv.size(0)," points ", pn.size(0), " planes");
          return Polyhedron();
        }
        // THIS USES STD::MOVE ON VPF!
        pout = Polyhedron(pv, pn, vpf);
        verbose_update("New ",pout.string_repr());
        // copy the updated vertices, normals, and relational information
        pv=pout.get_vertices();
        pn=pout.get_normals();
        vpf = pout.get_vertices_per_face();
      }
    }
    debug_update("Cut down to a ", pout.string_repr(), "\n", pout.python_string());
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
    auto vpf = lists_from_hdf<int>(group, "vertices_per_face");
    return {v, vpf};
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

  template<class T>
  Polyhedron point_bounding_box(const bArray<T>& points){
    auto xmin = points.min(0);
    auto xmax = points.max(0);
    std::vector<std::array<T,3>> v{
      {xmin[{0, 0}], xmin[{0,1}], xmin[{0,2}]}, // 000 0
      {xmin[{0, 0}], xmax[{0,1}], xmin[{0,2}]}, // 010 1
      {xmin[{0, 0}], xmax[{0,1}], xmax[{0,2}]}, // 011 2
      {xmin[{0, 0}], xmin[{0,1}], xmax[{0,2}]}, // 001 3
      {xmax[{0, 0}], xmin[{0,1}], xmin[{0,2}]}, // 100 4
      {xmax[{0, 0}], xmax[{0,1}], xmin[{0,2}]}, // 110 5
      {xmax[{0, 0}], xmax[{0,1}], xmax[{0,2}]}, // 111 6
      {xmax[{0, 0}], xmin[{0,1}], xmax[{0,2}]}  // 101 7
    };
    std::vector<std::vector<int>> vpf{{3,0,4,7},{3,2,1,0},{0,1,5,4},{3,7,6,2},{7,4,5,6},{2,6,5,1}};
    return Polyhedron(bArray<double>::from_std(v), vpf);
  }


} // end namespace brille
#endif // BRILLE_POLYHEDRON_H_
