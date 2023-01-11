#ifndef _BRILLE_POLYHEDRON_FACES_HPP_
#define _BRILLE_POLYHEDRON_FACES_HPP_

#include <atomic>
#include <random>
#include <utility>
#include <set>

#include "lattice_dual.hpp"
#include "array_.hpp"
#include "array_l_.hpp"
#include "geometry.hpp"
#include "approx_float.hpp"

namespace brille::polyhedron{
  template<class T, template<class> class A> std::enable_if_t<isBareArray<T,A>, A<T>>
  get_hkl(const A<T>& x){
    return x;
  }
  template<class T, template<class> class A> std::enable_if_t<isLatVec<T,A>, Array2<T>>
  get_hkl(const A<T>& x){
    return x.hkl();
  }


  // The vertex indices pointing into an external vertex array for a Polyhedron
  class Faces {
  public:
    using face_t = typename std::vector<ind_t>;
    using faces_t = typename std::vector<face_t>;
  protected:
    faces_t _faces;
  public:
    explicit Faces(): _faces() {}

    explicit Faces(faces_t faces): _faces(std::move(faces)) {}

    template<class T>
    explicit Faces(const std::vector<std::vector<T>>& faces){
      auto to_ind_t = [](const auto & i){return static_cast<ind_t>(i);};
      _faces.reserve(faces.size());
      for (const auto & f: faces){
        face_t _f; _f.reserve(f.size());
        std::transform(f.begin(), f.end(), std::back_inserter(_f), to_ind_t);
        _faces.push_back(_f);
      }
    }

//    template<class T, template<class> class A, std::enable_if_t<isArray<T,A>>* = nullptr> // this worked before moving isArray into brille::lattice::
//    template<class T, template<class> class A, std::enable_if_t<isArray<T,A>>* = nullptr> // doesn't work
    template<class T, template<class> class A> // hopefully doesn't capture anything non-brille::Array based
    explicit Faces(const A<T> & vertices, const T Ttol=T(0), const int tol=1){
      debug_update("Find convex hull of points\n", get_xyz(vertices).to_string());
      auto unique = vertices.is_unique(Ttol, tol);
      auto unique_vertices = vertices.extract(unique);
      auto[a, b, c] = find_convex_hull_planes(unique_vertices, Ttol, tol);
      auto fpv = find_planes_containing_point(a, b, c, unique_vertices, Ttol, tol);
      auto max = std::transform_reduce(fpv.begin(), fpv.end(), 0u, [](const auto & x, const auto & y){return x > y ? x : y;}, [](const auto & z){return z.size();});
      if (max > 2u) {
        _faces = polygon_faces(a, b, c, fpv, unique_vertices);
        debug_update("First-pass faces:\n", _faces);
        if (polygon_face_vertex_purge(unique_vertices, _faces)) {
          debug_update("Some vertices purged. Now vertices\nnp.array(\n",
                       get_xyz(unique_vertices).to_string(), "),\n", _faces);
          // the _faces vectors point into the smaller unique_vertices, so we
          // need to find *an* equivalent vertex in vertices for each retained
          // vertex in unique_vertices
          std::vector<ind_t> map(unique_vertices.size(0), vertices.size(0));
          for (ind_t i = 0; i < unique_vertices.size(0); ++i) {
            map[i] = vertices.first(cmp::eq, unique_vertices.view(i));
          }
          for (auto &face : _faces)
            for (auto &idx : face)
              idx = map[idx];
          debug_update("Face indexes updated to point into\nnp.array(\n",
                       vertices.to_string(), "),\n", _faces);
        }
      }
//      else {
//        info_update("The vertices \n", unique_vertices.to_string(), " are bounded by planes ",
//                    "\na:\n", a.to_string(), "\nb:\n", b.to_string(), "\nc:\n", c.to_string(),
//                    "but none are on three or more of the planes?!\n", fpv);
//      }
    }
    template<class T, template<class> class A>
    Faces(faces_t faces, const A<T> & vertices): _faces(std::move(faces)) {
      auto v = T(1) * vertices;
      std::tie(v, _faces) = remove_duplicate_points_and_update_face_indexing(v, _faces);
      _faces = polygon_faces(_faces, v);
      if (polygon_face_vertex_purge(v, _faces)){
        debug_update("Some vertices purged. Now vertices\nnp.array(\n", get_xyz(v).to_string(), "),\n", _faces);
        // the _faces vectors point into the smaller unique_vertices, so we
        // need to find *an* equivalent vertex in vertices for each retained
        // vertex in v
        std::vector<ind_t> map(v.size(0), vertices.size(0));
        for (ind_t i=0; i<v.size(0); ++i){
          map[i] = vertices.first(cmp::eq, v.view(i));
        }
        for (auto & face: _faces) for (auto & idx: face) idx = map[idx];
        debug_update("Face indexes updated to point into\nnp.array(\n",vertices.to_string(),"),\n",_faces);
      }
    }

    template<class T, template<class> class A>
    explicit Faces(faces_t faces, const A<T> & vertices, const A<T> & points, const A<T> & normals): _faces(std::move(faces)) {
      auto [a, b, c] = plane_points_from_normal(normals, points);
      auto faces_per_vertex = find_planes_containing_point(a, b, c, vertices);
      _faces = polygon_faces(a, b, c, faces_per_vertex, vertices);
    }

    Faces(const Faces & that) = default;
    Faces(Faces && that) noexcept : _faces(std::move(that._faces)) {}

    Faces& operator=(const Faces& that)= default;

    bool operator!=(const Faces&) const;
    bool operator==(const Faces& that) const {return !this->operator!=(that);}

    template<class I> Faces permute(const std::vector<I>& permutation) const {
      faces_t faces;
      faces.reserve(_faces.size());
      auto permute = [permutation](const auto & i) {return permutation[i];};
      for (const auto & face: _faces) {
        face_t perm;
        perm.reserve(face.size());
        std::transform(face.begin(), face.end(), std::back_inserter(perm), permute);
        faces.push_back(perm);
      }
      return Faces(faces);
    }

    // direct accessors
    [[nodiscard]] size_t size() const {return _faces.size();}
    [[nodiscard]] size_t face_count() const {return _faces.size();}
    [[nodiscard]] faces_t faces() const {return _faces;}
    [[nodiscard]] face_t face(const ind_t& i) const {
      assert(i < _faces.size());
      return _faces[i];
    }
    // calculated accessors
    [[nodiscard]] faces_t faces_per_vertex() const {return utils::invert_lists(_faces);}
    [[nodiscard]] face_t indexes() const {
      face_t all;
      for (const auto & f: _faces)
        for (const auto & x: f)
          if (std::find(all.begin(), all.end(), x)==all.end())
            all.push_back(x);
      return all;
    }
    [[nodiscard]] std::vector<std::vector<std::pair<ind_t, ind_t>>> edges_per_face() const;
    [[nodiscard]] Array2<ind_t> edges() const;
    [[nodiscard]] Array2<ind_t> planes() const;
    template<class T, template<class> class A>
    std::enable_if_t<isArray<T,A>, std::tuple<A<T>, A<T>, A<T>>>
    planes(const A<T>& x) const{
      auto p = this->planes();
      auto a = 0 * x.view(0); a.resize(p.size(0));
      auto b = 0 * x.view(0); b.resize(p.size(0));
      auto c = 0 * x.view(0); c.resize(p.size(0));
      for (ind_t i=0; i<p.size(0); ++i){
        auto i0 = p[{i, 0}];
        auto i1 = p[{i, 1}];
        auto i2 = p[{i, 2}];
        a.set(i, x.view(i0));
        b.set(i, x.view(i1));
        c.set(i, x.view(i2));
      }
      return std::make_tuple(a, b, c);
    }
    // calculated accessors requiring that vertices be provided
    // and returning the same type as the vertex array
    //! \brief Return the average of the vertices of a face, which *might* be the centroid
    template<class T, template<class> class A>
    std::enable_if_t<isArray<T,A>, A<T>> face_points(const A<T>&x) const {
      if (x.size(0) < 1){
        if (!_faces.empty()) throw std::runtime_error("No vertices provided for face points");
        return x;
      }
      auto p = 0 * x.view(0);
      p.resize(size());
      ind_t idx{0};
      for (const auto & face: _faces){
        p.set(idx++, x.extract(face).sum(0) / static_cast<T>(face.size()));
      }
      return p;
    }
    template<class T, template<class> class A>
    std::enable_if_t<isArray<T,A>, A<T>> face_normals(const A<T>& x) const {
      if (x.size(0) < 1){
        if (!_faces.empty()) throw std::runtime_error("No vertices provided for face normals");
        return x;
      }
      auto p = 0 * x.view(0);
      p.resize(size());
      ind_t idx{0};
      for (const auto & face: _faces){
        p.set(idx++, three_point_normal(x, face));
      }
      return p;
    }
    template<class T, template<class> class A>
    std::enable_if_t<isArray<T,A>, A<T>> half_edges(const A<T>& x) const {
      if (x.size(0) < 1){
        if (!_faces.empty()) throw std::runtime_error("No vertices provided for half edges");
        return x;
      }
      auto indexes = this->edges();
      auto half = 0 * x.view(0);
      half.resize(indexes.size(0));
      for (ind_t i=0; i<indexes.size(0); ++i){
        half.set(i, (x.view(indexes[{i, 0}]) + x.view(indexes[{i, 1}])) / T(2));
      }
      return half;
    }
    template<class T, template<class> class A>
    std::enable_if_t<isArray<T,A>, A<T>> centroid(const A<T>& x) const {
      if (x.size(0) < 1){
        if (!_faces.empty()) throw std::runtime_error("No vertices provided for centroid");
        return x;
      }
      auto cen = 0 * x.view(0);
      for (const auto & face: _faces) {
        const auto a{x.view(face[0])};
        for (ind_t i=1; i < face.size() - 1; ++i) {
          const auto b{x.view(face[i])};
          const auto c{x.view(face[i + 1])};
          auto cr = cross(b-a, c-a);
          auto sq = (((a+b)^T(2))+ ((b+c)^T(2)) + ((c+a)^T(2)));
          for (ind_t j=0; j<3; ++j) cen[j] += cr[j] * sq[j];
        }
      }
      return cen / (T(48) * this->volume(x));
    }
    template<class T, template<class> class A>
    std::enable_if_t<isArray<T,A>, A<T>> rand_rejection(const A<T>& x, const ind_t n, const unsigned int seed) const {
      auto tics = std::chrono::system_clock::now().time_since_epoch().count();
      std::default_random_engine gen(seed > 0 ? seed : static_cast<unsigned int>(tics));
      std::uniform_real_distribution<T> dst(T(0), T(1));
      auto direction = [&](){
        std::array<T,3> raw {{dst(gen), dst(gen), dst(gen)}};
        return from_std_like(x, raw);
      };
      auto min = x.min(0);
      auto delta = x.max(0) - min;
      auto points =  0 * x.view(0);
      points.resize(n);
      for (ind_t i=0; i < n; ){
        points.set(i, min + delta * direction());
        if (this->contains(x, points.view(i))[0]) ++i;
      }
      return points;
    }
    // or a scalar with type matching the array element type
    template<class T, template<class> class A>
    std::enable_if_t<isArray<T,A>, T> volume(const A<T>& x) const {
      /* per, e.g., http://wwwf.imperial.ac.uk/~rn/centroid.pdf
      For a polyhedron with N triangular faces, each with ordered vertices
      (aᵢ, bᵢ, cᵢ), one can define nᵢ = (bᵢ-aᵢ)×(cᵢ-aᵢ) for each face and then
      find that the volume of the polyhedron is V = 1/6 ∑ᵢ₌₁ᴺ aᵢ⋅ nᵢ
      */
      T volume{0};
      for (const auto & face: _faces){
        const auto a{x.view(face[0])};
        for (ind_t i=1; i < face.size() - 1; ++i){
          volume += dot(cross(x.view(face[i]) - a, x.view(face[i+1]) - a), a).sum();
        }
      }
      return volume / T(6);
    }
    template<class T, template<class> class A>
    std::enable_if_t<isArray<T,A>, T> minimum_distance(const A<T>& x) const {
      // For all pairs of (different) vertices, find the minimum distance
      // between any two:
      std::set<ind_t> index_set;
      for (const auto & face: _faces) for (const auto & idx: face) index_set.insert(idx);
      std::vector<ind_t> indexes;
      indexes.reserve(index_set.size());
      std::copy(index_set.cbegin(), index_set.cend(), std::back_inserter(indexes));
      if (indexes.size() < 2) return T(0);
      auto dist = (std::numeric_limits<T>::max)();
      for (ind_t i=0; i<x.size(0)-1; ++i){
        for (ind_t j=i+1; j<x.size(0); ++j){
          auto dij = norm(x.view(indexes[i]) - x.view(indexes[j])).sum();
          if (dij < dist) dist = dij;
        }
      }
      return dist;
    }
    template<class T, template<class> class A>
    std::enable_if_t<isArray<T,A>, T> circumsphere_radius(const A<T>& x) const{
      auto c2v = x - this->centroid(x);
      return norm(c2v).max(0).sum();
    }

    // return a modified copy of this Faces
    [[nodiscard]] Faces mirror() const {return Faces(reverse_each(_faces));}

    // geometric properties in relation to another point or polyhedron

//    template<class T, class R, template<class> class A, template<class> class B>
//    [[nodiscard]] std::enable_if_t<isArray<T,A> && isArray<R,B>, std::vector<bool>>
//    contains(const A<T>& v, const B<R>& x) const {
//      std::vector<bool> out(x.size(0));
//      // TODO Move n and p to per-thread variables instead of shared?
//      auto n = this->face_normals(v);
//      auto p = this->face_points(v);
//#pragma omp parallel for default(none) shared(out, n, p, x) schedule(dynamic)
//      for (ind_t i=0; i<x.size(0); ++i){
//        // FIXME, consider increasing the tolerance here!
//        out[i] = dot(n, x.view(i) - p).all(cmp::le, 0.);
//      }
//      return out;
//  }
    template<class T, class R, template<class> class A, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<T,A> && isArray<R,B>, std::vector<bool>>
    contains(const A<T>& v, const B<R>& x, const T t=T(0), const int n=1) const {
      std::vector<std::atomic<int>> tmp(x.size(0));
      A<T> pa, pb, pc;
      std::tie(pa, pb, pc) = this->planes(v);
//      const auto x_size = utils::u2s<long long>(x.size(0));
//      // making this parallel *also* makes it significantly slower!?
//#pragma omp parallel for default(none) shared(tmp, pa, pb, pc, x, x_size) schedule(dynamic)
//      for (long long i = 0; i < x_size; ++i) {
//        tmp[i] = point_inside_all_planes(pa, pb, pc, x.view(static_cast<ind_t>(i)), t, n) ? 1 : 0;
//      }
      for (ind_t i=0; i < x.size(0); ++i){
        tmp[i] = point_inside_all_planes(pa, pb, pc, x.view(i), t, n) ? 1 : 0;
      }
      std::vector<bool> out;
      out.reserve(tmp.size());
      std::transform(tmp.begin(), tmp.end(), std::back_inserter(out), [](const auto & i){return i>0;});
      return out;
    }

      /* point_inside_all_planes is now lattice-aware through pseudo_orient3d, no need to split its definition */
//    template<class T, class R, template<class> class A, template<class> class B>
//    [[nodiscard]] std::enable_if_t <isLatVec<T,A> && isLatVec<R,B>, std::vector<bool>>
//    contains(const A<T>& v, const B<R>& x) const {
//      assert(x.same_lattice(v) || x.star_lattice(v));
//      return this->contains(v.xyz(), x.xyz());
//    }


    template<class T, class R, template<class> class A, template<class> class B>
    [[nodiscard]] std::enable_if_t <isArray<T,A> &&  isArray<R,B>, size_t>
    face_index(const A<T>& v, const B<R>& a, const B<R>& b, const B<R>& c) const {
      auto match{_faces.size()};
      auto x = get_hkl(a);
      auto y = get_hkl(b);
      auto z = get_hkl(c);
      auto w = get_hkl(v);
      for (size_t i=0; i< _faces.size(); ++i) {
        const auto & face{_faces[i]};
        std::vector<double> o3d;
        o3d.reserve(face.size());
        for (const auto & vertex: face){
          o3d.push_back(orient3d(x.ptr(0), y.ptr(0), z.ptr(0), w.ptr(vertex)));
        }
        auto all_zero = std::all_of(o3d.begin(), o3d.end(), [](const auto & q){return approx_float::scalar(q, 0.);});
        if (all_zero) {
          all_zero = dot(three_point_normal(v, face[0], face[1], face[2]), three_point_normal(a, b, c)).all(cmp::gt, 0.);
        }
        if (all_zero) match = i;
      }
      return match;
    }
    template<class T, class R, template<class> class A, template<class> class B>
    [[nodiscard]] std::enable_if_t <isArray<T,A> &&  isArray<R,B>, size_t>
    has_face(const A<T>& v, const B<R>& a, const B<R>& b, const B<R>& c) const{
      return this->face_index(v,a,b,c) < _faces.size();
    }
    /*! \brief Check whether all vertices are behind a given plane
     * */
    template<class T, class R, template<class> class A, template<class> class B>
    [[nodiscard]] std::enable_if_t <isArray<T,A> &&  isArray<R,B>, bool>
    none_beyond(const A<T>& v, const B<R>& a, const B<R>& b, const B<R>& c) const {
      auto match = face_index(v, a, b, c); // look for a face which is the plane
      face_t f(v.size(0));
      std::iota(f.begin(), f.end(), 0);
      if (match < _faces.size()){
        // skip checking any vertices which are *on* the matched plane/face
        auto on_face = [s=_faces[match]](const auto & x){return std::find(s.begin(), s.end(), x) != s.end();};
        auto itr = std::remove_if(f.begin(), f.end(), on_face);
        f.erase(itr, f.end());
      }
      auto x = get_hkl(a);
      auto y = get_hkl(b);
      auto z = get_hkl(c);
      auto w = get_hkl(v);
      // using orient3d, w is 0 if its on the plane and > 0 if it's 'beneath' it
      auto is_negative = [&](const auto & i){
        return orient3d(x.ptr(0), y.ptr(0), z.ptr(0), w.ptr(i)) < 0;
      };
      // so all negative means no points are behind the plane and all
      // should be removed!
      return std::all_of(f.begin(), f.end(), is_negative);
    }

    template<class T, class R, template<class> class A, template<class> class B>
    [[nodiscard]] std::enable_if_t <isArray<T,A> &&  isArray<R,B>, std::tuple<A<T>, Faces>>
    divide(const A<T>& v, const B<R>& a, const B<R>& b, const B<R>& c) const {
      // this all seems unnecessary
      auto z = this->centroid(v);
      auto [cv, cf] = this->cut(v - z, a, b, c);
      return std::make_tuple(cv + z, cf);
    }

    template<class T, class R, template<class> class A, template<class> class B>
    [[nodiscard]] std::enable_if_t <isArray<T,A> &&  isArray<R,B>, std::tuple<A<T>, Faces>>
    one_cut(const A<T>& vin, const B<R>& a, const B<R>& b, const B<R>& c, const R Rtol=R(0), const int tol=1) const {
      assert(a.size(0) == b.size(0) && a.size(0) == c.size(0) && a.size(0) == 1);
      debug_update("plane passing through\nnp.array(\n", get_xyz(cat(0, a, b, c)).to_string(),")\n",
                     "(rlu) np.array(\n", cat(0, a, b, c).to_string(),")\n",
                     "with (rlu) normal ",three_point_normal(a,b,c).to_string());

      if (this->none_beyond(vin, a, b, c)) {
        return std::make_tuple(vin, Faces());
      }
      auto keep = point_inside_plane(a, b, c, vin);
//      verbose_update("keep ", keep, " of vertices\n", vin.to_string());
      if (std::find(keep.begin(), keep.end(), false) == keep.end()) return std::make_tuple(vin, *this);
      if (std::find(keep.begin(), keep.end(), true) == keep.end()) return std::make_tuple(vin, Faces());
      debug_update("keeping vertices ", keep, " of ", vin.to_string());
      auto edges_per_face = this->edges_per_face();
      auto cut_edges = keep_to_cut_edge_list(keep, edges_per_face);
      auto intersections = valid_edge_plane_intersections(a, b, c, vin, cut_edges);
      face_t new_face;
      ind_t vertex_count{0};
      auto v = T(1) * vin; // make a copy,  I hope
      auto f = this->faces(); // ditto
      auto pre_v_count = v.size(0);
      for (auto [first, second, edge, at]: intersections){
        auto index = v.first(cmp::eq, at, Rtol, Rtol, tol); // == v.size(0) if not found)
        if (index < keep.size()) {
          debug_update("Edge ", my_to_string(edge), " intersection ", at.to_string(0), " which is pre-existing vertex ", index, ": ", v.view(index).to_string(0));
          // ensure existing vertices which are 'cut' because they're on the plane are retained
          keep[index] = true;
        }
        if (pre_v_count <= index && index < v.size(0)) {
          // this should be impossible, since two edges which intersect a plane at the same point *must*
          // intersect it at *their* intersection point, which was a pre-existing vertex
          std::string msg =  "Duplicate intersection point for edge ";
          msg += my_to_string(edge) + " of\nnp.array(\n";
          msg += get_xyz(v).to_string() +  "),\n" +  lists_to_string(f) + "\n";
          msg += "intersection " + get_xyz(at).to_string(0);
          msg += " matches at index " + std::to_string(index);
          msg += "\n\nor in rlu:\nnp.array(\n" + v.to_string() + "),\n";
          msg += my_to_string(f) + "\n";
          msg += "intersection " + at.to_string(0) + " matches at index " + std::to_string(index);
          msg += "but there were only " + std::to_string(pre_v_count);
          msg += " vertices before the cut";
//          info_update(msg);
//          info_update("new face thus far: ", new_face);
//          info_update("comparison tolerance ", Rtol, " digits ", tol);
          msg += "\n\nconsider increasing tolerances from ";
          msg += my_to_string(Rtol) + " and " + std::to_string(tol);
          throw std::runtime_error(msg);
        }
        if (index == v.size(0)) {
//          debug_update("Edge ", my_to_string(edge), " intersection ", at.to_string(0), " which is a new vertex at ", index);
          v.append(0, at);
          ++vertex_count;
        }
        utils::add_between_if_missing(f[first], edge.first, edge.second, index);
        utils::add_between_if_missing(f[second], edge.second, edge.first, index);
        utils::add_if_missing(new_face, index);

      }
//      debug_update("new face: ", new_face);
      verbose_update_if(new_face.size()," with vertices \n", v.extract(new_face).to_string());
      // add the new face to the list of faces, if it is not present in the *input* faces
      // this might need to be more complex, we want to avoid doubling a face *or* partial face, but new_face has different indexing than the original
      if (new_face.size() > 2 && face_not_in_faces(new_face, this->faces())){
        f.push_back(sort_convex_polygon_face(a, b, c, new_face, v));
      }
      // We *may* have decided to keep all vertices and, in case of a rounding
      // error, also added new points. We should not add any points if we've
      // kept everything because keeping all vertices means no cut was performed
      if (std::find(keep.begin(), keep.end(), false) == keep.end()) return std::make_tuple(vin, *this);
      // extend to cover the newly found vertices
      for (ind_t z=0; z<vertex_count; ++z) keep.push_back(true);
//      debug_update("keeping vertices ", keep);

//      debug_update("np.array(", get_xyz(v).to_string(), ")\n,", f);

      // we *should* not sort faces again since we sorted the new face and inserted
      // any new vertices in the right spot already.
      std::tie(v, f) = remove_points_and_update_face_indexing(keep, v, f, false);
      // it's possible that some faces have extraneous edge points if rounding errors produced a 'new'
      // vertex partially along an existing edge
      f = remove_middle_colinear_points_from_faces(v, f, Rtol, tol);

//      debug_update("np.array(", get_xyz(v).to_string(), ")\n,", f);

      // remove any faces without three vertices
      f.erase(std::remove_if(f.begin(), f.end(), [](const auto & face){return unique(face).size() < 3;}), f.end());

      f = remove_dangling_faces(v.size(0), f);
      std::tie(v, f) = remove_faceless_points(v, f);
      return std::make_tuple(v, Faces(f));
    }

    template<class T, class R, template<class> class A, template<class> class B>
    [[nodiscard]] std::enable_if_t <isArray<T,A> && isArray<R,B>, std::tuple<A<T>, Faces>>
    cut(const A<T>& v, const B<R>&a, const B<R>&b, const B<R>&c, const R Rtol=R(0), const int tol=1) const {
      // this is, by far, the worst possible implementation of this algorithm
      assert(a.ndim()==2 && b.ndim()==2 && c.ndim() == 2);
      assert(a.size(1)==3 && b.size(1)==3 && c.size(1) == 3);
      assert(a.size(0) == b.size(0) && b.size(0) == c.size(0));
      auto ov = (T(1) * v).decouple();
      Faces of(*this);
      debug_update("Start cutting with ", a.size(0), " planes\nnp.array(", get_xyz(ov).to_string(),"),", of.python_string());
      for (ind_t i=0; i < a.size(0); ++i){
        std::tie(ov, of) = of.one_cut(ov, a.view(i), b.view(i), c.view(i), Rtol, tol);
        debug_update("After cut ", i, "\nnp.array(", get_xyz(ov).to_string(),"),", of.python_string());
      }
      return std::make_tuple(ov, of);
    }

#ifdef USE_HIGHFIVE
    template<class H> std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, bool> to_hdf(H& obj, const std::string& entry) const {
      bool ok{true};
      ok &= lists_to_hdf(_faces, obj, entry);
      return ok;
    }
    [[nodiscard]] bool to_hdf(const std::string& filename, const std::string& dataset, unsigned perm=HighFive::File::OpenOrCreate) const {
      HighFive::File file(filename, perm);
      return this->to_hdf(file, dataset);
    }
    template<class H> static std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, Faces> from_hdf(H& obj, const std::string& entry){
      faces_t faces = lists_from_hdf<ind_t>(obj, entry);
      return Faces(faces);
    }
    static Faces from_hdf(const std::string& filename, const std::string& dataset){
      HighFive::File file(filename, HighFive::File::ReadOnly);
      return Faces::from_hdf(file, dataset);
    }
#endif
    [[nodiscard]] Faces combine(const Faces& that, ind_t offset) const;

    [[nodiscard]] std::string python_string() const {
      return lists_to_string(_faces);
    }
  };


}
#endif
