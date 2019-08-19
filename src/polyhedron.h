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
                vertices_per_face(std::vector<std::vector<int>>()) {}
  // initalize from vertices, points, and three-plane intersection information
  Polyhedron(const ArrayVector<double>& v,
             const ArrayVector<double>& p
        )://,const ArrayVector<int>& fpv):
  vertices(v), points(p), normals(p/norm(p)) {
    this->keep_unique_vertices();
    this->find_all_faces_per_vertex();
    this->polygon_vertices_per_face();
    this->sort_polygons();
  }
  // initialize from vertices, points, and all relational information
  Polyhedron(const ArrayVector<double>& v,
             const ArrayVector<double>& p,
             const std::vector<std::vector<int>>& fpv,
             const std::vector<std::vector<int>>& vpf):
    vertices(v), points(p), normals(p/norm(p)), faces_per_vertex(fpv), vertices_per_face(vpf){}
  // initalize from vertices, points, normals, and three-plane intersection information
  Polyhedron(const ArrayVector<double>& v,
             const ArrayVector<double>& p,
             const ArrayVector<double>& n
        )://,const ArrayVector<int>& fpv):
  vertices(v), points(p), normals(n) {
    this->keep_unique_vertices();
    this->find_all_faces_per_vertex();
    this->polygon_vertices_per_face();
    this->sort_polygons();
  }
  // initialize from vertices, points, normals, and all relational information
  Polyhedron(const ArrayVector<double>& v,
             const ArrayVector<double>& p,
             const ArrayVector<double>& n,
             const std::vector<std::vector<int>>& fpv,
             const std::vector<std::vector<int>>& vpf):
    vertices(v), points(p), normals(n), faces_per_vertex(fpv), vertices_per_face(vpf){}
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
    // return Polyhedron(-1*this->vertices, -1*this->points, -1*this->normals, this->faces_per_vertex, reverse_each(this->vertices_per_face));
    return Polyhedron(-1*this->vertices, -1*this->points, -1*this->normals, this->faces_per_vertex, this->vertices_per_face);
  }
  template<class T> Polyhedron rotate(const std::array<T,9> rot) const {
    ArrayVector<double> newv(3u, this->vertices.size()), newp(3u, this->points.size()), newn(3u, this->normals.size());
    for (size_t i=0; i<this->vertices.size(); ++i)
      multiply_matrix_vector<double,T,double>(newv.datapointer(i), rot.data(), this->vertices.datapointer(i));
    for (size_t i=0; i<this->points.size(); ++i)
      multiply_matrix_vector<double,T,double>(newp.datapointer(i), rot.data(), this->points.datapointer(i));
    for (size_t i=0; i<this->normals.size(); ++i)
      multiply_matrix_vector<double,T,double>(newn.datapointer(i), rot.data(), this->normals.datapointer(i));
    return Polyhedron(newv, newp, newn, this->faces_per_vertex, this->vertices_per_face);
  }
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
    return Polyhedron(v,p,n, fpv, vpf);
  }
  ArrayVector<double> get_vertices(void) const { return vertices; }
  ArrayVector<double> get_points(void) const { return points; }
  ArrayVector<double> get_normals(void) const { return normals; }
  std::vector<std::vector<int>> get_faces_per_vertex(void) const { return faces_per_vertex; }
  std::vector<std::vector<int>> get_vertices_per_face(void) const {return vertices_per_face; }
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
        vector_cross(n, ba.datapointer(0), ca.datapointer(0));
        subvol = vector_dot(a.datapointer(0), n);
        if (vector_dot(n, normals.datapointer(f)) < 0) subvol *= -1.0;
        volume += subvol;
      }
    }
    return volume/6.0; // not-forgetting the factor of 1/6
  }
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
  void polygon_vertices_per_face(void) {
    // status_update(">");
    bool flag = true;
    // We have 3+ faces per vertex, so we can now find the vertices per face
    std::vector<std::vector<int>> vpf(this->points.size());
    for (size_t i=0; i<vertices.size(); ++i){
      for (auto facet: faces_per_vertex[i]){
        flag = true;
        for (auto vertex: vpf[facet]) if (static_cast<size_t>(vertex)==i) flag = false;
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
  }
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
  }
};

static bool ok_size(const size_t& a, const size_t& b){return (a==1u)||(a==b);}

template<class T> ArrayVector<bool> plane_outside_point(
  const ArrayVector<T>& n, const ArrayVector<T>& p, const ArrayVector<T>& x
){
  return dot(n, p-x).is_approx("<=", 0.); // true if x is closer to the origin
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
    ok.insert( plane_outside_point(n, p, at.extract(i)).all_true() , i);
  }
  return ok;
}
template<typename... L> bool one_intersection(L... args){
  return intersection(args...).getvalue(0);
}

/*! Find the polyhedron which results from slicing an existant polyhedron by
one or more plane passing through its volume. The part closest to the origin
is retained.*/
template<class T> Polyhedron bisect(const Polyhedron& pin, const ArrayVector<T>& n, const ArrayVector<T>& p) {
  if (n.numel()!=3u || p.numel()!=3u)
    throw std::runtime_error("Wrong number of vector elements");
  if (n.size()!=p.size())
    throw std::runtime_error("Equal number of vectors and points required");
  ArrayVector<bool> keep, valid, added(1u,0u);
  Polyhedron pout(pin);
  ArrayVector<double> pn, pp, pv=pout.get_vertices();
  std::vector<std::vector<int>> fpv, vpf;
  std::vector<int> del_vertices, cut_faces, face_vertices, new_vector;
  std::vector<int> vertex_map, face_map;
  ArrayVector<double> at(3u, 1u), ni(3u,1u), pi(3u,1u);
  size_t count;
  // copy the current vertices, normals, and relational information
  pv=pout.get_vertices();
  pn=pout.get_normals();
  pp=pout.get_points();
  fpv = pout.get_faces_per_vertex();
  vpf = pout.get_vertices_per_face();
  for (size_t i=0; i<n.size(); ++i){
    // check whether there's anything to do
    ni = n.extract(i);
    pi = p.extract(i);
    keep = plane_outside_point(ni, pi, pv);
    if (!keep.all_true()){
      // compile the list of to-be-deleted vertices
      del_vertices.empty();
      for (size_t j=0; j<keep.size(); ++j) if(!keep.getvalue(j))
        del_vertices.push_back(j);
      // and the list of facets which need to be cut or removed
      cut_faces.empty();
      for (auto x: del_vertices) for (auto f: fpv[x])
      if (std::find(cut_faces.begin(), cut_faces.end(), f)==cut_faces.end())
        cut_faces.push_back(f); // add f to the list if it's not present already
      // find the new intersection points of two neighbouring facets and the cutting plane
      count = 0;
      new_vector.empty();
      for (size_t j=0; j<cut_faces.size()-1; ++j)
      for (size_t k=j+1; k<cut_faces.size(); ++k)
      // check if the three planes intersect, and give a point not-outside of the polyhedron
      if (one_intersection(pn, pp, ni, pi, pn.extract(j), pp.extract(j), pn.extract(k), pp.extract(k), at)){
        // add the intersection point to all vertices
        pv = cat(pv, at);
        // add the new vertex to the list for each existing facet
        vpf[j].push_back(static_cast<int>(pv.size())-1);
        vpf[k].push_back(static_cast<int>(pv.size())-1);
        // and the yet-to-be-created facet
        new_vector.push_back(static_cast<int>(pv.size())-1);
        // keeping track of how many points we've added
        ++count;
        // plus add the face index to the faces_per_vertex list for this new vertex
        fpv.push_back({static_cast<int>(j), static_cast<int>(k), static_cast<int>(vpf.size())});
        // this will be a problem if two sets of facets give the same intersection point
      }
      // extend the vertices per face vector
      if (new_vector.size()) vpf.push_back(new_vector);
      // extend the keep ArrayVector<bool> to cover the new points
      added.resize(count);
      for (size_t j=0; j<count; ++j) added.insert(true, j);
      keep = cat(keep, added);
      // find the new indices for all vertices, and extract the kept vertices
      vertex_map.resize(pv.size());
      count = 0;
      for (size_t j=0; j<pv.size(); ++j)
        vertex_map[j]= (keep.getvalue(j) ? count++ : -1);
      pv = pv.extract(keep);
      // go through the vertices_per_face array, replacing vertex indicies
      // and making a face map to skip now-empty faces
      count = 0;
      for (size_t j=0; j<vpf.size(); ++j){
        new_vector.empty();
        for (auto x: vpf[j]) if (keep.getvalue(x)) new_vector.push_back(vertex_map[x]);
        vpf[j] = new_vector;
        face_map[j] = (vpf[j].size()>2 ? count++ : -1);
      }
      // make a face map, remove empty faces, and update the faces_per_vertex
      for (size_t j=0; j<fpv.size(); ++j){
        new_vector.empty();
        for (auto x: fpv[j]) if (face_map[x]>-1) new_vector.push_back(face_map[x]);
        fpv[j] = new_vector;
      }
      // we need to calculate the face centres
      pn = cat(pn, ni);
      pp = cat(pp, pi);
      for (size_t j=0; j<vpf.size(); ++j){
        at.resize(vpf[j].size());
        for (size_t k=0; k<vpf[j].size(); ++k) at.set(k, pv.extract(vpf[j][k]));
        pp.set(j, sum(at)/static_cast<double>(at.size()) );
      }

      pout = Polyhedron(pv, pp, pn);
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
template<class T>
Polyhedron polyhedron_box(std::array<T,3>& xmin, std::array<T,3>& xmax){
  ArrayVector<double> v(3u, 8u), p(3u, 6u);
  v.insert(xmin[0], 0, 0); v.insert(xmin[1], 0, 1); v.insert(xmin[2], 0, 2); // 000 0
  v.insert(xmin[0], 1, 0); v.insert(xmax[1], 1, 1); v.insert(xmin[2], 1, 2); // 010 1
  v.insert(xmin[0], 2, 0); v.insert(xmax[1], 2, 1); v.insert(xmax[2], 2, 2); // 011 2
  v.insert(xmin[0], 3, 0); v.insert(xmin[1], 3, 1); v.insert(xmax[2], 3, 2); // 001 3
  v.insert(xmax[0], 4, 0); v.insert(xmin[1], 4, 1); v.insert(xmin[2], 4, 2); // 100 4
  v.insert(xmax[0], 5, 0); v.insert(xmax[1], 5, 1); v.insert(xmin[2], 5, 2); // 110 5
  v.insert(xmax[0], 6, 0); v.insert(xmax[1], 6, 1); v.insert(xmax[2], 6, 2); // 111 6
  v.insert(xmax[0], 7, 0); v.insert(xmin[1], 7, 1); v.insert(xmax[2], 7, 2); // 101 7
  std::vector<std::vector<int>> vpf{{0,3,2,1},{4,5,6,7},{3,0,4,7},{1,2,6,5},{0,1,5,4},{3,7,6,2}};
  std::vector<std::vector<int>> fpv{{0,2,4},{0,3,4},{0,3,5},{0,2,5},{1,2,4},{1,3,4},{1,3,5},{1,2,5}};
  ArrayVector<T> tmp(3u, 4u);
  for (int i=0; i<6; ++i){
    for(int j=0; j<4; ++j) tmp.set(j, v.extract(vpf[i][j]));
    p.set(i, sum(tmp)/4.0);
  }
  return Polyhedron(v, p, p/norm(p), fpv, vpf);
}

#endif // _POLYHEDRON_H_
