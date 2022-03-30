#ifndef _BRILLE_POLYHEDRON_LATTICE_HPP_
#define _BRILLE_POLYHEDRON_LATTICE_HPP_

#include "lattice_dual.hpp"

#include <utility>
#include "array_l_.hpp"
namespace brille{
  template<class T>
  class LQPolyhedron {
  public:
    using face_t = typename std::vector<ind_t>;
    using faces_t = typename std::vector<face_t>;
    using vertex_t = lattice::LVec<T>;
  protected:
    vertex_t _vertices;
    faces_t _faces;
  public:
    explicit LQPolyhedron(): _vertices(), _faces() {}
    explicit LQPolyhedron(const vertex_t & vertices) : _vertices(vertices) { construct_convex_hull(); }
    explicit LQPolyhedron(vertex_t && vertices) : _vertices(vertices) { construct_convex_hull(); }
    LQPolyhedron(const vertex_t & vertices, faces_t faces) : _vertices(vertices), _faces(std::move(faces)) {
      construct_full_check();
    }
    LQPolyhedron(vertex_t && vertices, faces_t && faces) : _vertices(vertices), _faces(std::move(faces)) {
      construct_full_check();
    }
    [[deprecated("Can you use the .cut property instead of this constructor?")]]
    LQPolyhedron(const vertex_t & vertices, const vertex_t & points, const vertex_t & normals) : _vertices(vertices) {
      auto [b, c] = plane_points_from_normal(normals, points);
      auto faces_per_vertex = find_planes_containing_point(points, b, c, _vertices);
      _faces = polygon_faces(points, b, c, faces_per_vertex, _vertices);
      if (polygon_face_vetex_purge(_vertices, _faces)) debug_update("Some vertices purged");
    }
    LQPolyhedron(const LQPolyhedron<T> & that): _vertices(that._vertices), _faces(that._faces) {}
    LQPolyhedron(LQPolyhedron<T> && that) noexcept : _vertices(that._vertices), _faces(std::move(that._faces)) {}

    LQPolyhedron<T>& operator=(const LQPolyhedron<T>& that){
      _vertices = that.vertices();
      _faces = that.faces();
      return *this;
    }
    LQPolyhedron<T>& operator=(LQPolyhedron<T>&& that) noexcept {
      _vertices = that._vertices;
      _faces = std::move(that._faces);
      return *this;
    }

  private:
    void construct_convex_hull() {
      auto unique = _vertices.is_unique();
      _vertices = _vertices.extract(unique);
      auto[a, b, c] = find_convex_hull_planes(_vertices);
      auto fpv = find_planes_containing_point(a, b, c, _vertices);
      _faces = polygon_faces(a, b, c, fpv, _vertices);
      if (polygon_face_vertex_purge(_vertices, _faces)) {
        debug_update("Some vertices purged. Now vertices\n", _vertices.to_string(), "faces\n", _faces);
      }
      verbose_update("Finished constructing convex hull LQPolyhedron with", this->python_string());
    }
    void construct_full_check() {
      std::tie(_vertices, _faces) = remove_duplicate_points_and_update_face_indexing(_vertices, _faces);
      _faces = polygon_faces(_faces, _vertices);
      if (polygon_face_vertex_purge(_vertices, _faces)){
        debug_update("Some vertices purged. Now vertices\n", _vertices.to_string(), "faces\n", _faces);
      }
    }
  public:
    // direct accessors
    [[nodiscard]] ind_t vertex_count() const {return _vertices.size(0);}
    [[nodiscard]] size_t face_count() const {return _faces.size();}
    [[nodiscard]] vertex_t vertices() const {return _vertices;}
    [[nodiscard]] faces_t faces() const {return _faces;}
    // calculated accessors
    [[nodiscard]] faces_t faces_per_vertex() const {return utils::invert_lists(_faces);}
    [[nodiscard]] vertex_t face_points() const;
    [[nodiscard]] vertex_t face_normals() const;
    [[nodiscard]] Array2<ind_t> edges() const;
    [[nodiscard]] vertex_t half_edges() const;
    [[nodiscard]] T volume() const;
    [[nodiscard]] vertex_t centroid() const;
    [[nodiscard]] T circumsphere_radius() const;
    [[nodiscard]] vertex_t rand_rejection(ind_t, unsigned int seed=0) const;

    template<class R> bool operator!=(const LQPolyhedron<R>&) const;
    template<class R> bool operator==(const LQPolyhedron<R>& that) const {return !this->operator!=(that);}
    template<class R> LQPolyhedron<T> operator+(const LQPolyhedron<R>&) const;

    // return a modified copy of this LQPolyhedron
    LQPolyhedron<T> mirror() const {return {T(-1) * _vertices, reverse_each(_faces)};}
    LQPolyhedron<T> centre() const {return {_vertices - this->centroid(), _faces};}
    template<class R> LQPolyhedron<T> translate(const lattice::LVec<R>& vector) const {return {_vertices + vector, _faces};}
    template<class R> LQPolyhedron<T> rotate(const std::array<R,9>& rot) const;
    LQPolyhedron<T> apply(const PointSymmetry& ps, ind_t index) const;

    // geometric properties in relation to another point or polyhedron
    template<class R> [[nodiscard]] std::vector<bool> contains(const lattice::LVec<R>& x) const;
    template<class R> [[nodiscard]] std::vector<bool> contains(const std::vector<std::array<R,3>>& x) const{
      return this->contains(lattice::LVec<R>(_vertices.type(), _vertices.lattice(), Array2<R>::from_std(x)));
    }
    template<class R> [[nodiscard]] bool intersects(const LQPolyhedron<R>&) const;
    template<class R> [[nodiscard]] LQPolyhedron<T> intersection(const LQPolyhedron<R>&) const;
    template<class R> [[nodiscard]] LQPolyhedron<T> divide(const lattice::LVec<R>&, const lattice::LVec<R>&, const lattice::LVec<R>&) const;

    template<class R> [[nodiscard]] size_t face_index(const lattice::LVec<R>&, const lattice::LVec<R>&, const lattice::LVec<R>&) const;
    template<class R> [[nodiscard]] bool has_face(const lattice::LVec<R>& a, const lattice::LVec<R>& b, const lattice::LVec<R>& c) const {
      return this->face_index(a,b,c) < _faces.size();
    }
    template<class R> [[nodiscard]] bool none_beyond(const lattice::LVec<R>&, const lattice::LVec<R>&, const lattice::LVec<R>&) const;

    template<class R> LQPolyhedron<T> one_cut(const lattice::LVec<R>&, const lattice::LVec<R>&, const lattice::LVec<R>&) const;
    template<class R> LQPolyhedron<T> cut(const lattice::LVec<R>&, const lattice::LVec<R>&, const lattice::LVec<R>&) const;

#ifdef USE_HIGHFIVE
    template<class H> std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, bool> to_hdf(H& obj, const std::string& entry) const {
      auto group = overwrite_group(obj, entry);
      bool ok{true};
      ok &= _vertices.to_hdf(group, "vertices");
      ok &= lists_to_hdf(_faces, group, "faces");
      return ok;
    }
    [[nodiscard]] bool to_hdf(const std::string& filename, const std::string& dataset, unsigned perm=HighFive::File::OpenOrCreate) const {
      HighFive::File file(filename, perm);
      return this->to_hdf(file, dataset);
    }
    template<class H> static std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, LQPolyhedron<T>> from_hdf(H& obj, const std::string& entry){
      auto group = obj.getGroup(entry);
      vertex_t v = vertex_t::from_hdf(group, "vertices");
      faces_t faces = lists_from_hdf<ind_t>(group, "faces");
      return LQPolyhedron<T>(v, faces);
    }
    static LQPolyhedron<T> from_hdf(const std::string& filename, const std::string& dataset){
      HighFive::File file(filename, HighFive::File::ReadOnly);
      return LQPolyhedron<T>::from_hdf(file, dataset);
    }
#endif
  };
#include "polyhedron_lattice.tpp"

  template<class T> LQPolyhedron<T> bounding_box(const typename LQPolyhedron<T>::vertex_t&);
}
#endif