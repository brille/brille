#ifndef _POLYHEDRON_H_
#define _POLYHEDRON_H_

#include <vector>
#include "arrayvector.h"

template <typename T>
std::vector<T> bare_winding_angles(const ArrayVector<T>& vecs, const size_t i, const ArrayVector<T>& n){
  if (vecs.numel()!=3u)
    throw std::runtime_error("Finding a winding angle requires the cross product, which is only defined in 3D");
  // vecs should be normalized already
  std::vector<T> angles(vecs.size());
  T dotij, y_len, angij;
  ArrayVector<T> x(3u,1u), y(3u,1u); // ensure all have memory allocated
  T crsij[3]; // to hold the cross product
  for (size_t j=0; j<vecs.size(); ++j){
    if (j == i){
      angles[j] = 0.0;
      continue;
    }
    dotij = vecs.dot(i,j);
    vector_cross(crsij, vecs.datapointer(i), vecs.datapointer(j));
    // crsij = vecs.cross(i,j);
    x = dotij * vecs.extract(i);
    y = vecs.extract(j) - x;
    y_len = y.norm(0) * (std::signbit(vector_dot(crsij, n.datapointer(0))) ? -1 : 1);
    angij = std::atan2(y_len, dotij);
    angles[j] = angij < 0 ? angij+2*PI : angij;
  }
  return angles;
};
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
                vertices_per_face(std::vector<std::vector<int>>()) {};
  // initalize from vertices, points, and three-plane intersection information
  Polyhedron(const ArrayVector<double>& v,
             const ArrayVector<double>& p,
             const ArrayVector<int>& fpv):
  vertices(v), points(p), normals(p/norm(p)) {
    this->extend_unique_faces_per_vertex(fpv);
    this->polygon_vertices_per_face();
    this->sort_polygons();
  };
  // initialize from vertices, points, and all relational information
  Polyhedron(const ArrayVector<double>& v,
             const ArrayVector<double>& p,
             const std::vector<std::vector<int>>& fpv,
             const std::vector<std::vector<int>>& vpf):
    vertices(v), points(p), normals(p/norm(p)), faces_per_vertex(fpv), vertices_per_face(vpf){};
  // initalize from vertices, points, normals, and three-plane intersection information
  Polyhedron(const ArrayVector<double>& v,
             const ArrayVector<double>& p,
             const ArrayVector<double>& n,
             const ArrayVector<int>& fpv):
  vertices(v), points(p), normals(n) {
    this->extend_unique_faces_per_vertex(fpv);
    this->polygon_vertices_per_face();
    this->sort_polygons();
  };
  // initialize from vertices, points, normals, and all relational information
  Polyhedron(const ArrayVector<double>& v,
             const ArrayVector<double>& p,
             const ArrayVector<double>& n,
             const std::vector<std::vector<int>>& fpv,
             const std::vector<std::vector<int>>& vpf):
    vertices(v), points(p), normals(n), faces_per_vertex(fpv), vertices_per_face(vpf){};
  // copy constructor
  Polyhedron(const Polyhedron& other):
    vertices(other.get_vertices()),
    points(other.get_points()),
    normals(other.get_normals()),
    faces_per_vertex(other.get_faces_per_vertex()),
    vertices_per_face(other.get_vertices_per_face()) {};
  // assignment from another CentredPolyhedron
  Polyhedron& operator=(const Polyhedron& other){
    this->vertices = other.get_vertices();
    this->points = other.get_points();
    this->normals = other.get_normals();
    this->faces_per_vertex = other.get_faces_per_vertex();
    this->vertices_per_face = other.get_vertices_per_face();
    return *this;
  };
  ArrayVector<double> get_vertices(void) const { return vertices; };
  ArrayVector<double> get_points(void) const { return points; };
  ArrayVector<double> get_normals(void) const { return normals; };
  std::vector<std::vector<int>> get_faces_per_vertex(void) const { return faces_per_vertex; };
  std::vector<std::vector<int>> get_vertices_per_face(void) const {return vertices_per_face; };
protected:
  void extend_unique_faces_per_vertex(const ArrayVector<int>& fpv){
    // first look for vertices that are equivalent -- that is points where more than three planes interesect
    std::vector<bool> uniqueflg;
    std::vector<size_t> uniqueidx;
    for (size_t i=0; i<this->vertices.size(); ++i){
      uniqueflg.push_back(true);
      uniqueidx.push_back(i);
    }
    int tol_multiplier = 3; // a tuning parameter, 3 seems to work OK.
    for (size_t i=1; i<this->vertices.size(); ++i)
      for (size_t j=0; j<i; ++j){
        if (uniqueflg[j] && approx_vector(this->vertices.numel(), this->vertices.datapointer(i), this->vertices.datapointer(j), tol_multiplier)){
          uniqueflg[i]=false;
          uniqueidx[i]=j;
        }
        if (!uniqueflg[i]) break;
      }
    std::vector<size_t> mapidx(this->vertices.size());
    size_t no_unique=0;
    for (size_t i=0; i<this->vertices.size(); ++i)
      mapidx[i] = uniqueflg[i] ? no_unique++ : mapidx[uniqueidx[i]];
    std::vector<std::vector<int>> fpv_ext(no_unique);
    bool already_present;
    int tmp;
    for (size_t i=0; i<this->vertices.size(); ++i)
      for (size_t j=0; j<3; ++j){
        already_present=false;
        tmp = fpv.getvalue(i,j);
        for (auto k: fpv_ext[mapidx[i]]) if (k == tmp) already_present=true;
        if (!already_present) fpv_ext[mapidx[i]].push_back(tmp);
      }

    // and extract the unique vertices, and unique planes_per_vertex:
    this->vertices = this->vertices.extract(uniqueflg);
    this->faces_per_vertex = fpv_ext;
  };
  void polygon_vertices_per_face(void) {
    bool already_present;
    // We have 3+ faces per vertex, so we can now find the vertices per face
    std::vector<std::vector<int>> vpf(this->points.size());
    for (size_t i=0; i<this->faces_per_vertex.size(); ++i)
      for (auto face: this->faces_per_vertex[i]){
        // ensure that we don't somehow list a vertex multiple times for one face
        already_present=false;
        for (auto vert: vpf[face]) if (vert == i) already_present = false;
        if (!already_present) vpf[face].push_back(i);
      }
      // additionally, we only want to keep faces which describe polygons
    ArrayVector<bool> is_polygon(1u, vpf.size());
    for (size_t i=0; i<vpf.size(); ++i)
      is_polygon.insert(vpf[i].size()>2 ? true : false, i);
    this->prune_faces(is_polygon);
    // plus cut-down the vertices_per_face vector
    std::vector<std::vector<int>> polygon_vpf;
    for (auto i: vpf) if (i.size()>2) polygon_vpf.push_back(i);
    this->vertices_per_face = polygon_vpf;
  };
  void prune_faces(ArrayVector<bool>& keep){
    this->points = this->points.extract(keep);
    this->normals = this->normals.extract(keep);
  };
  void sort_polygons(void){
    std::vector<std::vector<int>> sorted_vpp(this->points.size());
    ArrayVector<double> facet_verts(3u, 0u), facet_centre, facet_normal;
    std::vector<int> facet, perm;
    std::vector<double> angles;
    double min_angle;
    size_t min_idx;
    ArrayVector<double> all_normals = this->get_normals();
    for (size_t j=0; j<this->points.size(); ++j){
      facet = this->vertices_per_face[j];
      facet_normal = all_normals.extract(j);
      facet_verts.resize(facet.size());
      for (size_t i=0; i<facet.size(); ++i) facet_verts.set(i, this->vertices.extract(facet[i]));
      facet_centre = sum(facet_verts)/static_cast<double>(facet.size());
      facet_verts -= facet_centre; // these are now on-face vectors to each vertex
      facet_verts = facet_verts/norm(facet_verts); // and now just their directions;
      perm.resize(facet.size());
      perm[0] = 0; // always start with whichever vertex is first
      for (size_t i=1; i<facet.size(); ++i){
        angles = bare_winding_angles(facet_verts, perm[i-1], facet_normal);
        min_angle = 1e3;
        min_idx=facet.size()+1;
        for (size_t k=0; k<facet.size(); ++k)
          if (!approx_scalar(angles[k], 0.0) && angles[k] < min_angle){
            min_idx=k;
            min_angle = angles[k];
          }
        if (min_idx >= facet.size()) throw std::runtime_error("Error finding minimum winding angle polygon vertex");
        perm[i] = min_idx;
      }
      for (size_t i=0; i<facet.size(); ++i) sorted_vpp[j].push_back(facet[perm[i]]); // this could be part of the preceeding loop.
    }
    this->vertices_per_face = sorted_vpp;
  };
};


#endif // _POLYHEDRON_H_
