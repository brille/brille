#ifndef _BRILLE_POLYHEDRON_FLEX_HPP_
#define _BRILLE_POLYHEDRON_FLEX_HPP_

#include "polyhedron.hpp"
#include "lattice_dual.hpp"

#include <utility>
#include "array_lvec.hpp"
#include "polyhedron_faces.hpp"

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
    explicit Poly(const vertex_t & vertices) : _vertices(vertices), _faces(vertices) {}
    explicit Poly(vertex_t && vertices) : _vertices(vertices),  _faces(vertices) {}
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

    template<class R, template<class> class B> std::enable_if_t<isArray<R,B>, bool>
      operator!=(const Poly<R,B>& that) const{
      bool vertices_permuted{false};
      if (_vertices != that._vertices){
        vertices_permuted = _vertices.is_permutation(that._vertices);
        if (!vertices_permuted) return true;
      }
      if (vertices_permuted) {
        auto permutation = _vertices.permutation_vector(that._vertices);
        auto permuted = that._faces.permute(permutation);
        return _faces != permuted;
      }
      return _faces != that._faces;
    }
    template<class R, template<class> class B> std::enable_if_t<isArray<R,B>, bool>
      operator==(const Poly<R,B>& that) const {return !this->operator!=(that);}
    template<class R, template<class> class B> std::enable_if_t<isArray<R,B>, Poly<T,A>>
      operator+(const Poly<R,B>& that) const{
      // combine vertices:
      auto vertices = cat(0, _vertices, that._vertices);
      // combine faces
      auto faces = _faces.combine(that._faces, _vertices.size(0)).faces();
      // check for duplicate vertices
      std::tie(vertices, faces) = remove_duplicate_points_and_update_face_indexing(vertices, faces);

      // check for overlapping faces
      /* such a check is hard, so we'll take a shortcut which is hopefully good enough:
       *    if any two faces contain three of the same vertices, they must be co-planar
       * we are, of course, discounting the possibility that two faces are coplanar and overlap without sharing points
       * but that seems *way* too difficult to handle at the moment. */
      for (ind_t i=0; i < faces.size()-1;  ++i) {
        const auto & a{faces[i]};
        for (ind_t j = i + 1; j < faces.size(); ++j) {
          const auto & b{faces[j]};
          auto c = std::count_if(a.begin(), a.end(), [&b](const auto & v){return std::count(b.begin(), b.end(), v) > 0;});
          //      if (c == a.size()){
          //        // they share all of their points, they're either duplicates or face each other
          //        // in either case, we need to check for equal vertex ordering with rolling the starting index of one
          //        // e.g., [0, 1, 2] == [0, 1, 2] || [1, 2, 3] || [2, 3, 0] || [3, 0, 1] for 4-point duplicate check
          //        // or    [0, 2, 1] == [0, 1, 2] || [1, 2, 3] || [2, 3, 0] || [3, 0, 1] for 4-point opposites check
          //      }
          if (c > 2) {
            std::string msg = "The faces " + std::to_string(i) + " and " + std::to_string(j);
            msg += " contain " + std::to_string(c) + " common vertex indexes, which indicates that they are coplanar!";
            throw std::runtime_error(msg);
            // If this gets thrown we should really think about handling this
          }
        }
      }
      return {vertices, Faces(faces)};
    }

    // return a modified copy of this Poly
    Poly<T,A> mirror() const {return {T(-1) * _vertices, _faces.mirror()};}
    Poly<T,A> centre() const {return {_vertices - this->centroid(), _faces};}
    template<class R> Poly<T,A> translate(const LQVec<R>& vector) const {return {_vertices + vector, _faces};}

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
      return {rotated, faces_t(post_faces)};
    }

    // geometric properties in relation to another point or polyhedron
    template<class R, template<class> class B>
    [[nodiscard]] std::vector<bool> contains(const B<R>& x) const {return _faces.contains(_vertices, x);}
    template<class R>
    [[nodiscard]] std::vector<bool> contains(const std::vector<std::array<R,3>>& x) const{
      return this->contains(from_std_like(_vertices, x));
    }

    //FIXME
    template<class R, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<R,B>, bool> intersects(const Poly<R,B>& that, const int tol=1) const
    {
      auto overlap = this->intersection(that);
      if (!approx::scalar(overlap.volume(), 0., tol)){
        info_update("this:\nP(",this->python_string(),")\nthat:\nP(",that.python_string(),")");
        info_update("overlap:\nP(",overlap.python_string(),")");
        overlap = this->intersection(that);
        throw std::runtime_error(std::to_string(overlap.volume()));
        return true;
      }
      return false;
    }
    template<class R, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<R,B>, Poly<T,A>> intersection(const Poly<R,B>& that) const{
      auto [a, b, c] = that.planes();
      return this->cut(a, b, c);
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
    one_cut(const B<R>& a, const B<R>& b, const B<R>& c, const int tol=1) const {
      auto [v, f] = _faces.one_cut(_vertices, a, b, c, tol);
      return {v, f};
    }
    template<class R, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<R,B>, Poly<T,A>>
    cut(const B<R>& a, const B<R>& b, const B<R>& c, const int tol=1) const {
      auto [v, f] = _faces.cut(_vertices, a, b, c, tol);
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
  };

  template<class T, template<class> class A>
  std::enable_if_t<isArray<T,A>, Poly<T,A>>
  bounding_box(const A<T>& points){
    auto min = points.min(0);
    auto max = points.max(0);
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
    auto faces = typename Poly<T,A>::faces_t({{3,0,4,7},{3,2,1,0},{0,1,5,4},{3,7,6,2},{7,4,5,6},{2,6,5,1}});
    auto hkl = Poly<T,A>::vertex_t::from_std(points.type(), points.lattice(), v);
    return {hkl, faces};
  }
}
#endif