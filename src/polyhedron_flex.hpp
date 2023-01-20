#ifndef _BRILLE_POLYHEDRON_FLEX_HPP_
#define _BRILLE_POLYHEDRON_FLEX_HPP_

#include <utility>
#include "lattice_dual.hpp"
#include "array_l_.hpp"
#include "polyhedron_faces.hpp"
#include "approx_float.hpp"

namespace brille::polyhedron{
  template<class T, template<class> class A>
  class Poly {
  public:
    using faces_t = Faces;
    using vertex_t = A<T>;
  protected:
    vertex_t _vertices;
    faces_t _faces;
  public:
    explicit Poly(): _vertices(), _faces() {}
    // Convex Hull constructors
    explicit Poly(const vertex_t & vertices, T Ttol=T(0), int tol=1): _vertices(vertices), _faces(vertices, Ttol, tol) {
      // This might be a bad idea, since we've been given a constant reference and might not want to make a copy
      this->finish_convex_hull(Ttol, tol);
    }
    explicit Poly(vertex_t && vertices, T Ttol=T(0), int tol=1) : _vertices(vertices),  _faces(vertices, Ttol, tol) {
      this->finish_convex_hull(Ttol, tol);
    }
    // Unverified constructors
    Poly(const vertex_t & vertices, faces_t faces) : _vertices(vertices), _faces(std::move(faces)) {}
    Poly(vertex_t && vertices, faces_t && faces) : _vertices(vertices), _faces(std::move(faces)) {}
    // Verified constructors
    Poly(const vertex_t & vertices, faces_t::faces_t proto_faces) : _vertices(vertices), _faces(proto_faces, vertices) {}
    Poly(const vertex_t & vertices, faces_t::faces_t && proto_faces) : _vertices(vertices), _faces(proto_faces, vertices) {}
    // Constructor from plane normals and points
    [[deprecated("Can you use the .cut property instead of this constructor?")]]
    Poly(const vertex_t & vertices, const vertex_t & points, const vertex_t & normals)
      :_vertices(vertices), _faces(vertices, points, normals) {
    }
    // copy constructors
    Poly(const Poly<T,A> & that): _vertices(that._vertices), _faces(that._faces) {}
    Poly(Poly<T,A> && that) noexcept : _vertices(that._vertices), _faces(std::move(that._faces)) {}
    // copy assignment
    Poly<T,A>& operator=(const Poly<T,A>& that){
      _vertices = that.vertices();
      _faces = that.faces();
      return *this;
    }
    Poly<T,A>& operator=(Poly<T,A>&& that) noexcept {
      _vertices = that._vertices;
      _faces = std::move(that._faces);
      return *this;
    }

  public:
    // direct accessors
    [[nodiscard]] ind_t vertex_count() const {return _vertices.size(0);}
    [[nodiscard]] size_t face_count() const {return _faces.size();}
    [[nodiscard]] vertex_t vertices() const {return _vertices;}
    [[nodiscard]] faces_t faces() const {return _faces;}
    // calculated accessors
    [[nodiscard]] faces_t::faces_t faces_per_vertex() const {return _faces.faces_per_vertex();}
    [[nodiscard]] vertex_t face_points() const {return _faces.face_points(_vertices);}
    [[nodiscard]] vertex_t face_normals() const {return _faces.face_normals(_vertices);};
    [[nodiscard]] Array2<ind_t> edges() const {return _faces.edges();};
    [[nodiscard]] vertex_t half_edges() const {return _faces.half_edges(_vertices);};
    [[nodiscard]] std::tuple<vertex_t, vertex_t, vertex_t> planes() const {return _faces.planes(_vertices);}
    [[nodiscard]] T volume() const {return _faces.volume(_vertices);}
    [[nodiscard]] vertex_t centroid() const {return _faces.centroid(_vertices);}
    [[nodiscard]] T circumsphere_radius() const {return _faces.circumsphere_radius(_vertices);}
    [[nodiscard]] vertex_t rand_rejection(const ind_t n, unsigned int seed=0) const {return _faces.rand_rejection(_vertices, n, seed);}

    [[nodiscard]] Poly<T,A> convex_hull() const {return Poly(_vertices);}

    bool is_not_approx(const Poly<T,A>& that, const T Ttol=T(0), const int tol=1) const {
      bool vertices_permuted{false};
      if (_vertices != that._vertices){
        vertices_permuted = _vertices.is_permutation(that._vertices, Ttol, Ttol, tol);
        if (!vertices_permuted) return true;
      }
      if (vertices_permuted) {
        auto permutation = _vertices.permutation_vector(that._vertices, Ttol, Ttol, tol);
        // auto permuted = that._faces.permute(permutation); // the permutation found permutes *our* vertices to theirs
        auto permuted = _faces.permute(permutation); // so we need to permute our face indices, not theirs
        return permuted != that._faces;
      }
      return _faces != that._faces;
    }
    bool is_approx(const Poly<T,A>& that, const T Ttol=T(0), const int tol=1) const {return !is_not_approx(that, Ttol, tol);}
//    template<class R, template<class> class B> std::enable_if_t<isArray<R,B>, bool>
    bool operator!=(const Poly<T,A>& that) const{ return is_not_approx(that);}
//    template<class R, template<class> class B> std::enable_if_t<isArray<R,B>, bool>
    bool operator==(const Poly<T,A>& that) const {return is_approx(that);}

//    template<class R, template<class> class B> std::enable_if_t<isArray<R,B>, Poly<T,A>>
    Poly<T,A>
    operator+(const Poly<T,A>& that) const {return this->combine(that);}

    Poly<T,A> combine(const Poly<T,A>& that, const T Ttol=T(0), const int tol=1) const {
      // combine vertices:
      auto vertices = cat(0, _vertices, that._vertices);
      // combine faces
      auto faces = _faces.combine(that._faces, _vertices.size(0)).faces();
      // check for duplicate vertices
      std::tie(vertices, faces) = remove_duplicate_points_and_update_face_indexing(vertices, faces, Ttol, tol);

      // check for overlapping faces
      /* such a check is hard, so we'll take a shortcut which is hopefully good enough:
       *    if any two faces contain three of the same vertices, they must be co-planar
       * we are, of course, discounting the possibility that two faces are coplanar and overlap without sharing points
       * but that seems *way* too difficult to handle at the moment. */
      std::vector<bool> needed(faces.size(), true);
      for (ind_t i=0; i < faces.size()-1;  ++i) if (needed[i]) {
        const auto & a{faces[i]};
        for (ind_t j = i + 1; j < faces.size(); ++j) if (needed[j]) {
          const auto & b{faces[j]};
          auto c = std::count_if(a.begin(), a.end(), [&b](const auto & v){return std::count(b.begin(), b.end(), v) > 0;});
          //      if (c == a.size()){
          //        // they share all of their points, they're either duplicates or face each other
          //        // in either case, we need to check for equal vertex ordering with rolling the starting index of one
          //        // e.g., [0, 1, 2] == [0, 1, 2] || [1, 2, 3] || [2, 3, 0] || [3, 0, 1] for 4-point duplicate check
          //        // or    [0, 2, 1] == [0, 1, 2] || [1, 2, 3] || [2, 3, 0] || [3, 0, 1] for 4-point opposites check
          //      }
          if (c > 2 && std::is_permutation(a.begin(), a.end(), b.begin())){
            // either these faces point in opposite directions and should be removed
            // or in the same direction and only one should be removed(?);
            needed[i] = is_positive_permutation(a, b); // keep the first face *if* they point the same way
            needed[j] = false;
            c = 0;
            // if neither are kept, this leaves behind a parting line. maybe we can handle this later?
          }
        }
      }
      if (std::find(needed.begin(), needed.end(), false) != needed.end()){
        for (ind_t i=0; i<faces.size(); ++i) if (!needed[i]) faces[i].clear();
        auto past = std::remove_if(faces.begin(), faces.end(), [](const auto& x){return x.empty();});
        faces.erase(past, faces.end()); // erase from past-last-good *to end*
        // figure out how to remove the unused vertices ...
//        std::remove_if();
      }
      return {vertices, Faces(faces)};
    }

    // return a modified copy of this Poly
    Poly<T,A> mirror() const {return {T(-1) * _vertices, _faces.mirror()};}
    Poly<T,A> centre() const {return {_vertices - this->centroid(), _faces};}
    template<class R> Poly<T,A> translate(const lattice::LVec<R>& vector) const {return {_vertices + vector, _faces};}

//    template<class R> Poly<T,A> rotate(const std::array<R,9>& rot) const {
//      // FIXME This can't compile since get_B_matrix only exists for Reciprocal
//      T B[9], invB[9], invB_R[9];
//      _vertices.get_lattice().get_B_matrix(B);
//      utils::matrix_inverse(invB, B);
//      utils::multiply_matrix_matrix(invB_R, invB, rot.data());
//      auto xyz = _vertices.xyz(); // this is B * (hkl)
//      for (ind_t i=0; i < xyz.size(0); ++i) {
//        utils::mul_mat_vec_inplace(3u, invB_R, xyz.ptr(i)); // rotate and convert back to (hkl)
//      }
//      return {xyz, _faces};
//    }
    Poly<T,A> apply(const PointSymmetry& ps, ind_t index) const {
      // FIXME This is only correct if A<T> == LQVec<T>!
      auto r_i = ps.get(index);
      if (LengthUnit::inverse_angstrom == _vertices.type()){
        // r_i is the rotation matrix for real-space vectors; but our vertices
        // are expressed in the reciprocal space so we must find its transposed
        r_i = transpose(r_i);
      }
      auto rotated =  0 * _vertices;
      for (ind_t i=0; i < _vertices.size(0); ++i){
        utils::multiply_matrix_vector(rotated.ptr(i), r_i.data(), _vertices.ptr(i));
      }
      // rotating has the possibility of reordering face vertices too...
      auto pre = this->face_normals();
      auto post = 0 * pre.view(0);
      typename faces_t::faces_t post_faces;
      for (ind_t i=0; i<pre.size(0); ++i){
        utils::multiply_matrix_vector(post.ptr(0), r_i.data(), pre.ptr(i));
        auto i_face = _faces.face(i);
        if (dot(post, three_point_normal(rotated, i_face)).sum() < 0){
          std::swap(i_face[1], i_face[2]);
          i_face = sort_convex_polygon_face(i_face, rotated);
        }
        post_faces.push_back(i_face);
      }
      auto faces = faces_t(post_faces);
      return {rotated, faces};
    }

    // geometric properties in relation to another point or polyhedron
    template<class R, template<class> class B, class... P>
    [[nodiscard]] std::vector<bool> contains(const B<R>& x, P... p) const {return _faces.contains(_vertices, x, p...);}
    template<class R, class... P>
    [[nodiscard]] std::vector<bool> contains(const std::vector<std::array<R,3>>& x, P... p) const{
      return this->contains(from_std_like(_vertices, x), p...);
    }

    //FIXME
    template<class R, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<R,B>, bool> intersects(const Poly<R,B>& that, const R Rtol=R(0), const int tol=1) const
    {
      auto overlap = this->intersection(that, Rtol, tol);
      if (!approx_float::scalar(overlap.volume() / (this->volume() + that.volume()), R(0), Rtol, Rtol, tol)){
//        info_update("this:\nP(",this->python_string(),")\nthat:\nP(",that.python_string(),")");
//        info_update("overlap:\nP(",overlap.python_string(),")");
//        overlap = this->intersection(that, Rtol, tol); // placed here to allow easier bad-state access to intersection
//        throw std::runtime_error("Overlap volume should be zero but is " + my_to_string(overlap.volume()));
        return true;
      }
      return false;
    }
    template<class R, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<R,B>, Poly<T,A>> intersection(const Poly<R,B>& that, const R Rtol=R(0), const int tol=1) const{
      auto [a, b, c] = that.planes();
      return this->cut(a, b, c, Rtol, tol);
    }
    template<class R, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<R,B>, Poly<T,A>> divide(const B<R>& a, const B<R>& b, const B<R>& c) const{
      auto [v, f] = _faces.divide(_vertices, a, b, c);
      return {v, f};
    }
    template<class R, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<R,B>, size_t> face_index(const B<R>& a, const B<R>& b, const B<R>& c) const{
      return _faces.face_index(_vertices, a, b, c);
    }
    template<class R, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<R,B>, bool> has_face(const B<R>& a, const B<R>& b, const B<R>& c) const {
      return _faces.has_face(_vertices, a, b, c);
    }
    template<class R, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<R,B>, bool> none_beyond(const B<R>& a, const B<R>& b, const B<R>& c) const{
      return _faces.none_beyond(_vertices, a, b, c);
    }

    template<class R, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<R,B>, Poly<T,A>>
    one_cut(const B<R>& a, const B<R>& b, const B<R>& c, const R Rtol=R(0), const int tol=1) const {
      auto [v, f] = _faces.one_cut(_vertices, a, b, c, Rtol, tol);
      return {v, f};
    }
    template<class R, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<R,B>, Poly<T,A>>
    cut(const B<R>& a, const B<R>& b, const B<R>& c, const R Rtol=R(0), const int tol=1) const {
      auto [v, f] = _faces.cut(_vertices, a, b, c, Rtol, tol);
//      verbose_update("Cut produced vertices\n", v.to_string(), "and faces\n", f.python_string());
      return {v, f};
    }

#ifdef USE_HIGHFIVE
    template<class H> std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, bool> to_hdf(H& obj, const std::string& entry) const {
      auto group = overwrite_group(obj, entry);
      bool ok{true};
      ok &= _vertices.to_hdf(group, "vertices");
      ok &= _faces.to_hdf(group, "faces");
      return ok;
    }
    [[nodiscard]] bool to_hdf(const std::string& filename, const std::string& dataset, unsigned perm=HighFive::File::OpenOrCreate) const {
      HighFive::File file(filename, perm);
      return this->to_hdf(file, dataset);
    }
    template<class H> static std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, Poly<T,A>> from_hdf(H& obj, const std::string& entry){
      auto group = obj.getGroup(entry);
      vertex_t v = vertex_t::from_hdf(group, "vertices");
      faces_t faces = faces_t::from_hdf(group, "faces");
      return Poly<T,A>(v, faces);
    }
    static Poly<T,A> from_hdf(const std::string& filename, const std::string& dataset){
      HighFive::File file(filename, HighFive::File::ReadOnly);
      return Poly<T,A>::from_hdf(file, dataset);
    }
#endif
    [[nodiscard]] std::string python_string() const {
      return "np.array("+get_xyz(_vertices).to_string()+"),"+_faces.python_string();
    }
  private:
    void finish_convex_hull(const T, const int) {
      // check for unused vertices and remove them
      auto present = _faces.indexes(); // unordered list of _vertices indexes present/used
      std::sort(present.begin(), present.end());
      //
      std::vector<ind_t> indexes(_vertices.size(0));
      std::iota(indexes.begin(), indexes.end(), 0u);
      // if all vertex indexes are in present, there's nothing to do
      if (std::includes(present.begin(), present.end(), indexes.begin(), indexes.end())) return;

      // make a map from old to new indexes, and an extraction list
      std::vector<bool> keep(indexes.size(), false);
      std::fill(indexes.begin(), indexes.end(), indexes.size());
      ind_t kept{0};
      for (const auto & x: present){
        indexes[x] = kept++;
        keep[x] = true;
      }
      // update faces vectors
      auto faces = _faces.faces();
      for (auto & face: faces) for (auto & x: face) x = indexes[x];
      _faces = Faces(faces);
      // and extract the used vertices
      _vertices = _vertices.extract(keep);
    }
	friend std::ostream& operator<<(std::ostream& os, const Poly<T,A> & p){
	  os << p.python_string();
		return os;
	}
  };

  template<class T, template<class> class A>
  std::enable_if_t<isArray<T,A>, Poly<T,A>>
  bounding_box(const A<T>& points){
//    verbose_update("bounding box of\n", points.to_string(), " in xyz frame\n", get_xyz(points).to_string());
    auto min = get_xyz(points).min(0);
    auto max = get_xyz(points).max(0);
//    verbose_update("Have extreme corners, ", min.to_string(0), " -- ", max.to_string(0));
    std::vector<std::array<T,3>> v{
        {min[{0, 0}], min[{0,1}], min[{0,2}]}, // 000 0
        {min[{0, 0}], max[{0,1}], min[{0,2}]}, // 010 1
        {min[{0, 0}], max[{0,1}], max[{0,2}]}, // 011 2
        {min[{0, 0}], min[{0,1}], max[{0,2}]}, // 001 3
        {max[{0, 0}], min[{0,1}], min[{0,2}]}, // 100 4
        {max[{0, 0}], max[{0,1}], min[{0,2}]}, // 110 5
        {max[{0, 0}], max[{0,1}], max[{0,2}]}, // 111 6
        {max[{0, 0}], min[{0,1}], max[{0,2}]}  // 101 7
    };
//    verbose_update("producing polygon vertices\n", v);
    auto faces = typename Poly<T,A>::faces_t({{3,0,4,7},{3,2,1,0},{0,1,5,4},{3,7,6,2},{7,4,5,6},{2,6,5,1}});
    auto hkl = from_xyz_like(points, bArray<T>::from_std(v));
//    verbose_update("which should match \n", hkl.xyz().to_string());
    return {hkl, faces};
  }


}
#endif
