#ifndef _POLYHEDRON_H_
#define _POLYHEDRON_H_

#include <vector>
#include "arrayvector.h"
#include "debug.h"

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

template<class T> std::vector<T> reverse(const std::vector<T>& x){
  std::vector<T> r;
  for (size_t i = x.size(); i--;) r.push_back(x[i]);
  return r;
}
template<class T> std::vector<std::vector<T>> reverse_each(const std::vector<std::vector<T>>& x){
  std::vector<std::vector<T>> r;
  for (auto i: x) r.push_back(reverse(i));
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
                vertices_per_face(std::vector<std::vector<int>>()) {};
  // initalize from vertices, points, and three-plane intersection information
  Polyhedron(const ArrayVector<double>& v,
             const ArrayVector<double>& p,
             const ArrayVector<int>& fpv):
  vertices(v), points(p), normals(p/norm(p)) {
    this->keep_unique_vertices();
    this->find_all_faces_per_vertex();
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
    this->keep_unique_vertices();
    this->find_all_faces_per_vertex();
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
  Polyhedron mirror(void) const {
    // return Polyhedron(-1*this->vertices, -1*this->points, -1*this->normals, this->faces_per_vertex, reverse_each(this->vertices_per_face));
    return Polyhedron(-1*this->vertices, -1*this->points, -1*this->normals, this->faces_per_vertex, this->vertices_per_face);
  };
  template<class T> Polyhedron rotate(const std::array<T,9> rot) const {
    ArrayVector<double> newv(3u, this->vertices.size()), newp(3u, this->points.size()), newn(3u, this->normals.size());
    for (size_t i=0; i<this->vertices.size(); ++i)
      multiply_matrix_vector<double,T,double>(newv.datapointer(i), rot.data(), this->vertices.datapointer(i));
    for (size_t i=0; i<this->points.size(); ++i)
      multiply_matrix_vector<double,T,double>(newp.datapointer(i), rot.data(), this->points.datapointer(i));
    for (size_t i=0; i<this->normals.size(); ++i)
      multiply_matrix_vector<double,T,double>(newn.datapointer(i), rot.data(), this->normals.datapointer(i));
    return Polyhedron(newv, newp, newn, this->faces_per_vertex, this->vertices_per_face);
  };
  Polyhedron operator+(const Polyhedron& other) const {
    size_t d = this->vertices.numel();
    if (other.vertices.numel() != d) throw std::runtime_error("Only equal dimensionality polyhedra can be combined.");
    // combine all vertices, points, and normals; adjusting the fpv and vpf indexing
    size_t tvn = this->vertices.size(), ovn = other.vertices.size();
    size_t tfn = this->points.size(), ofn = other.points.size();
    ArrayVector<double> v(d, tvn+ovn);
    ArrayVector<double> p(d, tfn+ofn), n(d, tfn+ofn);
    for (size_t i=0; i<tvn; ++i) v.set(i,     this->vertices.extract(i));
    for (size_t i=0; i<ovn; ++i) v.set(tvn+i, other.vertices.extract(i));
    for (size_t i=0; i<tfn; ++i) p.set(i,     this->points.extract(i));
    for (size_t i=0; i<ofn; ++i) p.set(tfn+i, other.points.extract(i));
    for (size_t i=0; i<tfn; ++i) n.set(i,     this->normals.extract(i));
    for (size_t i=0; i<ofn; ++i) n.set(tfn+i, other.normals.extract(i));
    std::vector<std::vector<int>> fpv(this->faces_per_vertex), vpf(this->vertices_per_face);
    fpv.resize(tvn+ovn); vpf.resize(tfn+ofn);
    for (size_t i=0; i<ovn; ++i) for (auto j: other.faces_per_vertex[i]) fpv[tvn+i].push_back(tfn+j);
    for (size_t i=0; i<ofn; ++i) for (auto j: other.vertices_per_face[i]) vpf[tfn+i].push_back(tvn+j);

    // // look for equivalent vertices to remove
    // std::vector<bool> flg; std::vector<size_t> idx, map;
    // for (size_t i=0; i<v.size(); ++i){ flg.push_back(true); idx.push_back(i);}
    // int t = 3; //tolerance multiplier
    // for (size_t i=1; i<v.size(); ++i) for (size_t j=0; j<i; ++j){
    //   if (flg[i]&&flg[j]) flg[i]=!approx_vector(d, v.datapointer(i), v.datapointer(j), t);
    //   if (!flg[i]){ idx[i]=j; break; }
    // }
    // // combine the fpv for equivalent non-unique vertices
    // bool notpresent;
    // for (size_t i=0; i<v.size(); ++i) if (flg[i])
    //   for (size_t j=i+1; j<v.size(); ++j) if (i==idx[j])
    //     for (auto fj: fpv[j]){
    //       notpresent = true;
    //       for (auto fi: fpv[i]) if (fi==fj) notpresent = false;
    //       if (notpresent) fpv[i].push_back(fj);
    //     }
    // // we want map[i] to give the entry of the *extracted* unique vertices
    // size_t num=0;
    // for (size_t i=0; i<v.size(); ++i) map.push_back(flg[i]?num++:map[idx[i]]);
    // // extract the unique vertices
    // v = v.extract(flg);
    // // and the unique fpv:
    // std::vector<std::vector<int>> ufpv;
    // for (size_t i=0; i<v.size(); ++i) if (flg[i]) ufpv.push_back(fpv[i]);
    // // and replace the vertex indices in vpf using the map
    // for (size_t i=0; i<p.size(); ++i) for (size_t j=0; j<vpf[i].size(); ++i) if (!flg[vpf[i][j]]) vpf[i][j] = map[vpf[i][j]];
    //
    // // look for equivalent facets to remove
    // flg.resize(0); idx.resize(0); map.resize(0);
    // for (size_t i=0; i<p.size(); ++i){ flg.push_back(true); idx.push_back(i);}
    // for (size_t i=1; i<p.size(); ++i) for (size_t j=0; j<i; ++j){
    //   if (flg[i]&&flg[j]) flg[i] = !(approx_vector(d, p.datapointer(i), p.datapointer(j), t)&&approx_vector(d, n.datapointer(i), n.datapointer(j), t));
    //   if (!flg[i]){ idx[i]=j; break; }
    // }
    // // combine the vpf for equivalent non-unique facets
    // for (size_t i=0; i<p.size(); ++i) if (flg[i])
    // for (size_t j=i+1; j<p.size(); ++j) if (i==idx[j])
    // for (auto vj: vpf[j]){
    //   notpresent = true;
    //   for (auto vi: vpf[i]) if (vi==vj) notpresent = false;
    //   if (notpresent) vpf[i].push_back(vj);
    // }
    // num = 0; for (size_t i=0; i<p.size(); ++i) map.push_back(flg[i]?num++:map[idx[i]]);
    // p = p.extract(flg);
    // n = n.extract(flg);
    // std::vector<std::vector<int>> uvpf;
    // for (size_t i=0; i<p.size(); ++i) if (flg[i]) uvpf.push_back(vpf[i]);
    // for (size_t i=0; i<v.size(); ++i) for (size_t j=0; j<ufpv[i].size(); ++i) if (!flg[ufpv[i][j]]) ufpv[i][j] = map[ufpv[i][j]];
    //
    // // finally, double check that vpf and fpv only contain unique entries
    // vpf.resize(0); fpv.resize(0);
    // for (auto i: uvpf){
    //   vpf.push_back(std::vector<int>());
    //   for (auto j: i){
    //     notpresent = true;
    //     for (auto k: vpf.back()) if (j==k) notpresent = false;
    //     if (notpresent) vpf.back().push_back(j);
    //   }
    // }
    // for (auto i: ufpv){
    //   fpv.push_back(std::vector<int>());
    //   for (auto j: i){
    //     notpresent = true;
    //     for (auto k: fpv.back()) if (j==k) notpresent = false;
    //     if (notpresent) fpv.back().push_back(j);
    //   }
    // }

    return Polyhedron(v,p,n, fpv, vpf);
  }
  ArrayVector<double> get_vertices(void) const { return vertices; };
  ArrayVector<double> get_points(void) const { return points; };
  ArrayVector<double> get_normals(void) const { return normals; };
  std::vector<std::vector<int>> get_faces_per_vertex(void) const { return faces_per_vertex; };
  std::vector<std::vector<int>> get_vertices_per_face(void) const {return vertices_per_face; };
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
    if (found != nfv>>1)
      throw std::runtime_error("Problem finding all half-edge points.");
    return hep;
  };
  std::string string_repr(void) const {
    size_t nv = vertices.size(), nf=points.size();
    std::string repr = "Polyhedron with ";
    repr += std::to_string(nv) + " " + (1==nv?"vertex":"vertices") + " and ";
    repr += std::to_string(nf) + " " + (1==nf?"facet":"facets");
    debug_exec(repr += "; volume " +std::to_string(this->get_volume());)
    return repr;
  };
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
    std::array<int,3> tri;
    ArrayVector<double> a(3u, 1u), ba(3u, 1u), ca(3u, 1u);
    for (size_t f=0; f<normals.size(); ++f){
      a = this->vertices.extract(vertices_per_face[f][0]);
      for (int i=1; i<vertices_per_face[f].size()-1; ++i){ // loop over triangles
        ba = this->vertices.extract(vertices_per_face[f][ i ]) - a;
        ca = this->vertices.extract(vertices_per_face[f][i+1]) - a;
        vector_cross(n, ba.datapointer(0), ca.datapointer(0));
        subvol = vector_dot(a.datapointer(0), n);
        if (vector_dot(n, normals.datapointer(f)) < 0) subvol *= -1.0;
        volume += subvol;
      }
    }
    return volume/6.0; // not-forgetting the factor of 1/6
  };
protected:
  void keep_unique_vertices(void){
    // status_update(">");
    std::vector<bool> flg;
    for (size_t i=0; i<vertices.size(); ++i) flg.push_back(true);
    int t = 3; // a tolerance multiplier tuning parameter, 3 seems to work OK.
    size_t n = vertices.numel();
    for (size_t i=1; i<vertices.size(); ++i) for (size_t j=0; j<i; ++j)
      if (flg[i]&&flg[j]) flg[i]=!approx_vector(n, vertices.datapointer(i), vertices.datapointer(j), t);
    this->vertices = this->vertices.extract(flg);
    // status_update("<");
  }
  // This is the function which is bad for non-convex polyhedra, since it
  // assigns all faces with n⋅(v-p)=0 to a vertex irrespective of whether two
  // faces have opposite direction normals.
  void find_all_faces_per_vertex(void){
    // status_update(">");
    ArrayVector<double> vmp;
    std::vector<std::vector<int>> fpv(vertices.size());
    ArrayVector<bool> isonplane(1u, points.size());
    for (size_t i=0; i<vertices.size(); ++i){
      isonplane = dot(normals, vertices.extract(i) - points).is_approx("==",0.);
      for (size_t j=0; j<points.size(); ++j) if (isonplane.getvalue(j)) fpv[i].push_back(static_cast<int>(j));
    }
    this->faces_per_vertex = fpv;
    // status_update("<");
  }
  // void extend_unique_faces_per_vertex(const ArrayVector<int>& fpv){
  //   // status_update(">");
  //   // first look for vertices that are equivalent -- that is points where more than three planes interesect
  //   std::vector<bool> flg; std::vector<size_t> idx;
  //   for (size_t i=0; i<vertices.size(); ++i) flg.push_back(true);
  //   int t = 3; // a tolerance multiplier tuning parameter, 3 seems to work OK.
  //   size_t n = vertices.numel();
  //   for (size_t i=1; i<vertices.size(); ++i) for (size_t j=0; j<i; ++j){
  //     if (flg[i]&&flg[j]) flg[i]=!approx_vector(n, vertices.datapointer(i), vertices.datapointer(j), t);
  //     idx.push_back(flg[i] ? i : j); << This won't work as anticipated. it *adds* an entry for every j
  //     if (!flg[i]) break;
  //   }
  //   // extract the unique vertices
  //   this->vertices = this->vertices.extract(flg);
  //   // idx[i] gives the entry of vertices which is unique.
  //   // we want map[i] to give the entry of the *extracted* unique vertices
  //   std::vector<size_t> map; size_t num=0;
  //   for (size_t i=0; i<vertices.size(); ++i) map.push_back(flg[i]?num++:map[idx[i]]);
  //   // so that we can find all facets on wich each vertex sits:
  //   std::vector<std::vector<int>> ext(num);
  //   bool already_present; int tmp;
  //   for (size_t i=0; i<fpv.size(); ++i) for (size_t j=0; j<3u; ++j){
  //     already_present = false;
  //     tmp = fpv.getvalue(i,j);
  //     for (auto k: ext[map[i]]) if (k==tmp) already_present = true;
  //     if (!already_present) ext[map[i]].push_back(tmp);
  //   }
  //   // save the combined facets per vertex information:
  //   this->faces_per_vertex = ext;
  //   // status_update("<");
  // };
  void polygon_vertices_per_face(void) {
    // status_update(">");
    bool flag = true;
    // We have 3+ faces per vertex, so we can now find the vertices per face
    std::vector<std::vector<int>> vpf(this->points.size());
    for (size_t i=0; i<vertices.size(); ++i){
      for (auto facet: faces_per_vertex[i]){
        flag = true;
        for (auto vertex: vpf[facet]) if (vertex==i) flag = false;
        if (flag) vpf[facet].push_back(i);
      }
    }

    // additionally, we only want to keep faces which describe polygons
    std::vector<bool> is_polygon;
    for (size_t i=0; i<vpf.size(); ++i) is_polygon.push_back(vpf[i].size()>2);
    this->points = this->points.extract(is_polygon);
    this->normals = this->normals.extract(is_polygon);

    // we should modify faces_per_vertex here, to ensure its indexing is correct
    size_t count = 0, max=vpf.size(); std::vector<size_t> map;
    for (size_t i=0; i<max; ++i) map.push_back(is_polygon[i] ? count++ : max);
    std::vector<std::vector<int>> reduced_fpv(faces_per_vertex.size());
    for (size_t i=0; i<faces_per_vertex.size(); ++i)
    for (auto facet: faces_per_vertex[i]) if (is_polygon[facet]) reduced_fpv[i].push_back(map[facet]);
    this->faces_per_vertex = reduced_fpv;

    // plus cut-down the vertices_per_face vector
    std::vector<std::vector<int>> polygon_vpf;
    for (auto i: vpf) if (i.size()>2) polygon_vpf.push_back(i);
    this->vertices_per_face = polygon_vpf;
    // status_update("<");
  };
  void sort_polygons(void){
    // status_update(">");
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
    // status_update("<");
  };
};


#endif // _POLYHEDRON_H_
