#ifndef BRILLE_TRELLIS_NODE_HPP_
#define BRILLE_TRELLIS_NODE_HPP_
/*! \file
    \author Greg Tucker
    \brief A class holding a hybrid grid of cuboid and triangulated tetrahedral
           cells and data for interpolation
*/
#include <queue>
#include <condition_variable>
#include <atomic>
#include <functional>
#include <utility>
#include "enums.hpp"
#include "triangulation_simple.hpp"
#include "interpolatordual.hpp"
#include "hdf_interface.hpp"
namespace brille {



  /*! \brief A base class for the differentiation of Node types

  The NullNode is not within the domain of the PolyhedronTrellis because it has
  a null intersection with the bounding Polyhedron.
  */
  class NullNode{
  public:
    //! Implicit construction of an empty NullNode
    NullNode() = default;
    //! Deconstruction of a NullNode
    virtual ~NullNode() = default;
    //! Return the type of this Node
    [[nodiscard]] virtual NodeType type() const {return NodeType::null;}
    //! Return the number of vertices this Node indexes
    [[nodiscard]] virtual ind_t vertex_count() const {return 0u;}
    //! Return the vertex indices of this Node
    [[nodiscard]] virtual std::vector<ind_t> vertices() const {return {};}
    //! Return the triangulated tetrahedra indices of this Node
    [[nodiscard]] virtual std::vector<std::array<ind_t,4>> vertices_per_tetrahedron() const {return {};}
    virtual //! Return the indices required and their weights for linear interpolation at a point
    bool indices_weights(const bArray<double>&,
                         const bArray<double>&,
                         std::vector<std::pair<ind_t,double>>&,
                         const bool
                        ) const {return false;}
    //! Return the volume of this Node
    [[nodiscard]] virtual double volume(const bArray<double>&) const {return 0.;}
#ifdef USE_HIGHFIVE
    //! Write to an HDF file
    template<class R>
    std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, bool>
    to_hdf(R&, const std::string&) const {
      return false;
    }
    //! Read from an HDF file
    template<class R>
    static std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, NullNode>
    from_hdf(R&, const std::string&) {return NullNode();}
    bool operator!=(const NullNode&) const { return true; }
#endif
  };
  /*! \brief A Node fully within the domain of the PolyhedronTrellis

  The eight vertices of this cuboid Node are all within the volume of the
  convex Polyhedron bounding the domain of the Polyhedrontrellis.
  It contains an ordered list of its vertex indices within the full list of
  all PolyhedronTrellis vertices.
  */
  class CubeNode: public NullNode {
  protected:
    /*!< \brief The eight vertex indices of the Node

    The vertex order is critical for the CubeNode and must be:
      (000), (100), (110), (010), (101), (001), (011), (111)
    */
    std::array<ind_t, 8> vertex_indices;
  public:
    bool operator!=(const CubeNode& other) const {
      return vertex_indices != other.vertex_indices;
    }
    bool operator==(const CubeNode& o) const {return !this->operator!=(o);}
    //! Implicit construction of a CubeNode with no volume
    CubeNode(): vertex_indices({{0,0,0,0,0,0,0,0}}) {}
    //! Construct from an array
    explicit CubeNode(const std::array<ind_t,8>& vi): vertex_indices(vi) {}
    //! Construct from a vector with 8 elements
    explicit CubeNode(const std::vector<ind_t>& vi): vertex_indices({{0,0,0,0,0,0,0,0}}) {
      if (vi.size() != 8) throw std::logic_error("CubeNode objects take 8 indices.");
      for (ind_t i=0; i<8u; ++i) vertex_indices[i] = vi[i];
    }
    //! Return the number of vertices this Node indexes
    [[nodiscard]] ind_t vertex_count() const override { return 8u;}
    //! Return the vertex indices of this Node
    [[nodiscard]] std::vector<ind_t> vertices() const override {
      std::vector<ind_t> out;
      for (auto v: vertex_indices) out.push_back(v);
      return out;
    }
    /*!\brief Return the indices required and their weights for linear
              interpolation at a point

    \param vertices all vertex positions of the PolyhedronTrellis
    \param        x the point at which linear interpolation is required
    \param[out]  iw the minimal number of vertex indices and their interpolation
                    weights required for linear interpoaltion within the Node
    \returns whether this Node contains `x` and `iw` has been set
    If the interpolation point is
      - one of the indexed PolyhedronTrellis points then that point index and
        a weight of 1 are set in `iw`.
      - on an edge of the Node the two vertices on the endpoints of the edge are
        set in `iw` along with weights proportional to the distance from the point
        to the *other* vertex -- that is, the line segment length with a vertex
        replaced by the point is proportional to the interpolation weight for that
        vertex.
      - on a face of the Node the four vertices at the corners of the face are
        set in `iw` with weights proportional to the area of the rectangle formed
        by the test point and the three other face vertices.
      - within the volume of the Node the eight corners indices are set in `iw`
        along with weights proportional to the volume of the cuboid formed by the
        point and the other seven Node corners.
      - outside of the volume no inidices or weights are set and the returned bool
        is false.
    */
    template<class T, class I, template<class> class A>
    std::enable_if_t<isArray<T,A>, bool>
    indices_weights(const A<T>& vertices, const A<T>& x, std::vector<std::pair<I,T>>& iw, const bool should_contain) const {
      // The CubeNode object contains the indices into `vertices` necessary to find
      // the 8 corners of the cube. Those indices should be ordered
      // (000) (100) (110) (010) (101) (001) (011) (111)
      // so that vertex_indices[i] and vertex_indices[7-i] are connected by a body diagonal
      auto node_verts = vertices.extract(vertex_indices);
      T node_volume = abs(node_verts.view(0)-node_verts.view(7)).prod(1)[0];
      auto w = abs(x - node_verts).prod(1)/node_volume; // the normalised volume of each sub-parallelpiped
      // If any normalised weights are greater than 1+eps() the point isn't in this node
      iw.clear();
      if (w.any(brille::cmp::gt, T(1)) && !should_contain) return false;
      auto needed = w.is(brille::cmp::gt, 0.);
      for (int i=0; i<8; ++i) if (needed[i]) {
        // the weight corresponds to the vertex opposite the one used to find the partial volume
        iw.emplace_back(vertex_indices[7-i],w[i]);
      }
      return true;
    }
    /*\brief Determine the volume of this Node

    \param vertices all vertex positions of the PolyhedronTrellis which this Node
                    indexes
    Since the vertex indices are in a defined order the body diagonal is easily
    extracted from two vertices and the Node volume is the product of the body
    diagonal elements.
    */
    template<class T, template<class> class A>
    [[nodiscard]]
    std::enable_if_t<isArray<T,A>, T>
    volume(const A<T>& vertices) const {
      // The CubeNode object contains the indices into `vertices` necessary to find
      // the 8 corners of the cube. Those indices should be ordered
      // (000) (100) (110) (010) (101) (001) (011) (111)
//      // so that vertex_indices[i] and vertex_indices[7-i] are connected by a body diagonal
//      return abs(vertices.view(vertex_indices[0])-vertices.view(vertex_indices[7])).prod(1)[0];
      // If the vertex coordinates are not absolut units, the body diagonal is not necessarily the node volume
      auto a = vertices.view(vertex_indices[1]);
      auto b = vertices.view(vertex_indices[0]);
      auto c = vertices.view(vertex_indices[3]);
      auto d = vertices.view(vertex_indices[5]);
      return pseudo_orient3d(a, b, c, d)[0];
    }
#ifdef USE_HIGHFIVE
    //! Write to an HDF file
    template<class R>
    std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, bool>
    to_hdf(R& obj, const std::string& name) const {
      if (obj.exist(name)) obj.unlink(name);
      auto group = obj.createGroup(name);
      group.createDataSet("vertex_indices", vertex_indices);
      return true;
    }
    //! Read from an HDF file
    template<class R>
    static
    std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, CubeNode>
    from_hdf(R& obj, const std::string& name){
      std::array<ind_t, 8> vi{};
      obj.getGroup(name).getDataSet("vertex_indices").read(vi);
      return CubeNode(vi);
    }
#endif
  };
  /*! \brief A Node at least partly within the domain of the PolyhedronTrellis

  If any of the eight vertices of Node are outside of the bounding Polyhedron of
  the PolyhedronTrellis then the surface of the Polyhedron passes through its
  volume. Since the surface of the Polyhedron are where degeneracies are allowed
  we do not want to interpolate across the surface and instead must find the
  intersection of the Node with the bounding Polyhedron. That intersection is
  itself another Polyhedron which can then be triangulated for use in linear
  interpolation.

  A PolyNode is always the intersection of a Trellis Node and the bounding
  Polyhedron and typically has lower volume than the Node.
  */
  class PolyNode: public NullNode {
  protected:
    std::vector<std::array<ind_t,4>> vi_t;  //!< vertex indices per triangulated tetrahedron
    std::vector<std::array<double,4>> ci_t; //!< circumsphere information per triangulated tetrahedra
    std::vector<double> vol_t;              //!< volume per triangulated tetrahedra
  public:
    bool operator!=(const PolyNode& other) const {
      if (vi_t != other.vi_t) return true;
      if (ci_t != other.ci_t) return true;
      if (vol_t != other.vol_t) return true;
      return false;
    }
    bool operator==(const PolyNode& o) const {return !this->operator!=(o);}
    //
    explicit PolyNode() = default;
    // actually constructing the tetrahedra from, e.g., a Polyhedron object will
    // need to be done elsewhere
    //! Construct with all parameters defined
    PolyNode(
      std::vector<std::array<ind_t,4>>  vit,
      std::vector<std::array<double,4>>  cit,
      std::vector<double>  volt
    ): vi_t(std::move(vit)), ci_t(std::move(cit)), vol_t(std::move(volt)) {}
    //! Return the number of triangulated tetrahedra in the PolyNode
    [[nodiscard]] ind_t tetrahedra_count() const {return static_cast<ind_t>(vi_t.size());}
    //! Return the number of unique vertex indices in the triangulated tetrahedra
    [[nodiscard]] ind_t vertex_count() const override { return static_cast<ind_t>(this->vertices().size());}
    //! Return the unique vertex indices from all triangulated tetrahedra
    [[nodiscard]] std::vector<ind_t> vertices() const override {
      std::vector<ind_t> out;
      for (auto tet: vi_t) for (auto idx: tet)
      if (std::find(out.begin(), out.end(), idx)==out.end()) out.push_back(idx);
      return out;
    }
    //! Return the vertex indices for each triangulated tetrahedra
    [[nodiscard]] std::vector<std::array<ind_t,4>> vertices_per_tetrahedron() const override {return vi_t;}
    /*!\brief Return the indices required and their weights for linear
              interpolation at a point

    \param vertices all vertex positions of the PolyhedronTrellis
    \param        x the point at which linear interpolation is required
    \param[out]  iw the minimal number of vertex indices and their interpolation
                    weights required for linear interpoaltion within the Node
    \returns whether a triangulated tetrahedra contains `x` and `iw` has been set

    If one of the triangulated tetrahedra contains the interpolation point and it
    is
      - one of the indexed PolyhedronTrellis points then that point index and
        a weight of 1 are set in `iw`.
      - on an edge of the containing tetrahedra the two vertices on the endpoints
        of the edge are set in `iw` along with weights proportional to the
        distance from the point to the *other* vertex -- that is, the line segment
        length with a vertex replaced by the point is proportional to the
        interpolation weight for that vertex.
      - on a face of the tetrahedra the three vertices at the corners of the face
        are set in `iw` with weights proportional to the area of the triangle
        formed by the test point and the two other face vertices.
      - within the volume of the tetrahedra the four corners indices are set in
        `iw` along with weights proportional to the volume of the tetrahedron
        formed by the point and the other three Node corners.

    If the point is not inside any of the triangulated tetrahedra no inidices
    or weights are set and the returned bool is false.
    */
    template<class T, class I, template<class> class A>
    std::enable_if_t<isArray<T,A>, bool>
    indices_weights(
      const A<T>& vertices,
      const A<T>& x,
      std::vector<std::pair<I,T>>& iw,
      const bool should_contain
    ) const {
      iw.clear();
      std::array<T,4> w{{0,0,0,0}};
      std::vector<T> most_neg(vi_t.size());
      for (ind_t i=0; i<vi_t.size(); ++i){
        most_neg[i] = tetrahedra_contains(i, vertices, x, w);
        if (most_neg[i] >= 0.){
          for (int j=0; j<4; ++j) if (!brille::approx_float::scalar(w[j], 0.))
            iw.emplace_back(static_cast<I>(vi_t[i][j]), w[j]);
          return true;
        }
      }
      // orient3d seems to be more precise than the inclusion checks done
      // thus far :/ so we need an escape-hatch in case a point 'should' be
      // in *one* of the contained tetrahedra but "isn't".
      if (should_contain){
        auto i = std::distance(most_neg.begin(), std::max_element(most_neg.begin(), most_neg.end()));
        verbose_update("PolyNode does not contain ",x.to_string(0)," but should\n\tUsing tetrahedron with ",most_neg[i]," weight for one vertex");
        bool allow_shortcut{false}; // just incase we chose a really bad tetrahedron
        auto val = tetrahedra_contains(static_cast<ind_t>(i), vertices, x, w, allow_shortcut);
        if (val > 0.){
          throw std::runtime_error("This shouldn't be possible");
        }
        for (int j=0; j<4; ++j) if (!brille::approx_float::scalar(w[j], 0.))
          iw.emplace_back(static_cast<I>(vi_t[i][j]), w[j]);
        return true;
      }
      return false;
    }
    //! Return the total triangulated volume of the PolyNode
    template<class T, template<class> class A>
    [[nodiscard]]
    std::enable_if_t<isArray<T,A>, T>
    volume(const A<T>&) const {
      return std::accumulate(vol_t.begin(), vol_t.end(), 0.);
    }
  private:
    template<template<class> class A>
    std::enable_if_t<isArray<double,A>, double>
    tetrahedra_contains(
      const ind_t t,
      const A<double>& v,
      const A<double>& x,
      std::array<double,4>& w,
      const bool allow_shortcut = true
    ) const {
      if (allow_shortcut){
        auto away = this->tetrahedra_might_contain(t,x);
        if (away < 0.) return away;
      }
      double vol6 = vol_t[t]*6.0;
      // use pseudo_orient3d which is lattice-aware
      w[0] = pseudo_orient3d(x, v.view(vi_t[t][1u]), v.view(vi_t[t][2u]), v.view(vi_t[t][3u]))[0]/vol6;
      w[1] = pseudo_orient3d(v.view(vi_t[t][0u]), x, v.view(vi_t[t][2u]), v.view(vi_t[t][3u]))[0]/vol6;
      w[2] = pseudo_orient3d(v.view(vi_t[t][0u]), v.view(vi_t[t][1u]), x, v.view(vi_t[t][3u]))[0]/vol6;
      w[3] = pseudo_orient3d(v.view(vi_t[t][0u]), v.view(vi_t[t][1u]), v.view(vi_t[t][2u]), x)[0]/vol6;
      if (std::any_of(w.begin(), w.end(), [](double z){return z < 0. && !brille::approx_float::scalar(z, 0.);}))
        return *std::min_element(w.begin(), w.end());
      return 0.;
    }
//    template<class T, template<class> class A>
//    [[nodiscard]] std::enable_if_t<isLatVec<T,A>, T>
//    tetrahedra_contains(const ind_t t, const A<T>& v, const A<T>& x, std::array<double,4>& w, const bool allow_shortcut=true) const {
//      // FIXME replace orient3d with a lattice-aware volume calculator
//      //    Then the two methods here could be combined and copying v avoided
//      return this->tetrahedra_contains(t, v.xyz(), x.xyz(), w, allow_shortcut);
//    }
    template<template<class> class A>
    [[nodiscard]] std::enable_if_t<isBareArray<double,A>, double>
    tetrahedra_might_contain(
      const ind_t t,
      const A<double>& x
    ) const {
      // find the vector from the circumsphere centre to x:
      double v[3];
      v[0]=ci_t[t][0]-x.val(0,0);
      v[1]=ci_t[t][1]-x.val(0,1);
      v[2]=ci_t[t][2]-x.val(0,2);
      // compute the squared length of v, and the circumsphere radius squared
      double d2{0}, r2 = ci_t[t][3]*ci_t[t][3];
      for (double i : v) d2 += i*i;
      // if the squared distance is no greater than the squared radius, x might be inside the tetrahedra
      //return d2 < r2 || brille::approx_float::scalar(d2, r2);
      return d2 < r2 || brille::approx_float::scalar(d2, r2) ? 0. : -d2;
    }
    template<class T, template<class> class A>
    [[nodiscard]] std::enable_if_t<isLatVec<T,A>, T>
    tetrahedra_might_contain(const ind_t t, const A<T>& x) const {
      return this->tetrahedra_might_contain(t, x.xyz());
    }

#ifdef USE_HIGHFIVE
  public:
    //! Write to an HDF file
    template<class R>
    std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, bool>
    to_hdf(R& obj, const std::string& name) const {
      using namespace HighFive;
      if (obj.exist(name)) obj.unlink(name);
      auto group = obj.createGroup(name);

      std::vector<std::vector<ind_t>> dv;
      std::vector<std::vector<double>> dc;
      for (const auto& vi: vi_t){
        std::vector<ind_t> i;
        for (const auto& v: vi) i.push_back(v);
        dv.push_back(i);
      }
      for (const auto& ci: ci_t){
        std::vector<double> i;
        for (const auto& v: ci) i.push_back(v);
        dc.push_back(i);
      }
      group.createDataSet("vi_t", dv);
      group.createDataSet("ci_t", dc);

      group.createDataSet("vol_t", vol_t);
      return true;
    }
    //! Read from an HDF file
    template<class R>
    static
        std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, PolyNode>
        from_hdf(R& obj, const std::string& name){
      std::vector<std::vector<ind_t>> dv;
      std::vector<std::vector<double>> dc;
      std::vector<std::array<ind_t, 4>> vi;
      std::vector<std::array<double, 4>> ci;
      std::vector<double> vol;
      auto group = obj.getGroup(name);
      group.getDataSet("vi_t").read(dv);
      group.getDataSet("ci_t").read(dc);
      group.getDataSet("vol_t").read(vol);

      for (const auto& v: dv){
        vi.push_back(std::array<ind_t,4>({v[0], v[1], v[2], v[3]}));
      }
      for (const auto& c: dc){
        ci.push_back(std::array<double,4>({c[0], c[1], c[2], c[3]}));
      }

      return {vi, ci, vol};
    }
#endif //USE_HIGHFIVE
  };

  /*! \brief A utility class to hold and index both CubeNode and PolyNode objects

  As both CubeNode and PolyNode are subclasses of NullNode they can not be easily
  distinguished in templates and some form of more complex differentiation is
  required.

  The NodeContainer contains a vector which contains the indices of all contained
  CubeNode and PolyNode objects within their respective vectors. The indexing
  vector holds at each element a pair of NodeType and index into the node type's
  vector.
  */
  class NodeContainer{
    using nodes_t = std::vector<std::pair<NodeType, ind_t>>;
    using cubes_t = std::vector<CubeNode>;
    using polys_t = std::vector<PolyNode>;
  protected:
    nodes_t nodes_;
    cubes_t cube_nodes_;
    polys_t poly_nodes_;
  public:
    explicit NodeContainer() = default;
    NodeContainer(nodes_t&& n, cubes_t&& c, polys_t&& p)
        : nodes_(std::move(n)), cube_nodes_(std::move(c)), poly_nodes_(std::move(p)) {}
    NodeContainer(const std::vector<NodeType>& types) {
      nodes_.reserve(types.size());
      ind_t n_cube{0}, n_poly{0}, null_idx{(std::numeric_limits<ind_t>::max)()};
      for (const auto & x: types){
        switch (x) {
        case NodeType::cube: nodes_.emplace_back(x, n_cube++); break;
        case NodeType::poly: nodes_.emplace_back(x, n_poly++); break;
        default: nodes_.emplace_back(x, null_idx);
        }
      }
      cube_nodes_.resize(n_cube);
      poly_nodes_.resize(n_poly);
    }
    bool operator!=(const NodeContainer& other) const {
      if (nodes_ != other.nodes_) return true;
      if (cube_nodes_ != other.cube_nodes_) return true;
      if (poly_nodes_ != other.poly_nodes_) return true;
      return false;
    }
    std::vector<NodeType> all_node_types() const {
      std::vector<NodeType> types;
      types.reserve(nodes_.size());
      for (const auto & n: nodes_) types.push_back(n.first);
      return types;
    }
    //
    //! Return the total number of index nodes
    [[nodiscard]] size_t size() const {return nodes_.size();}
    //! Return the number of indexed CubeNode objects
    [[nodiscard]] size_t cube_count() const {
      return std::count_if(nodes_.begin(),nodes_.end(),[](std::pair<NodeType,ind_t> n){return NodeType::cube == n.first;});
    }
    //! Return the nuber of indexed PolyNode objects
    [[nodiscard]] size_t poly_count() const {
      return std::count_if(nodes_.begin(),nodes_.end(),[](std::pair<NodeType,ind_t> n){return NodeType::poly == n.first;});
    }
    //! Return the nubmer of indexed NullNode objects
    [[nodiscard]] size_t null_count() const {
      return std::count_if(nodes_.begin(),nodes_.end(),[](std::pair<NodeType,ind_t> n){return NodeType::null == n.first;});
    }
    //! Push a CubeNode onto the back of the container
    void push_back(const CubeNode& n){
      nodes_.emplace_back(NodeType::cube, static_cast<ind_t>(cube_nodes_.size()));
      cube_nodes_.push_back(n);
    }
    //! Push a PolyNode onto the back of the container
    void push_back(const PolyNode& n){
      if (n.vertex_count() < 1)
        throw std::runtime_error("empty polynodes are not allowed!");
      nodes_.emplace_back(NodeType::poly, static_cast<ind_t>(poly_nodes_.size()));
      poly_nodes_.push_back(n);
    }
    //! Push a NullNode onto the back of the container
    void push_back(const NullNode&){
      nodes_.emplace_back(NodeType::null, (std::numeric_limits<ind_t>::max)());
    }
    void set(const ind_t i, CubeNode&& n){
      assert(i < nodes_.size());
      assert(nodes_[i].first == NodeType::cube);
      assert(nodes_[i].second < cube_nodes_.size());
      cube_nodes_[nodes_[i].second] = std::move(n);
    }
    void set(const ind_t i, PolyNode&& n){
      assert(i < nodes_.size());
      assert(nodes_[i].first == NodeType::poly);
      assert(nodes_[i].second < poly_nodes_.size());
      poly_nodes_[nodes_[i].second] = std::move(n);
    }
    void set(const ind_t, NullNode){
    }
    //! Return the NodeType of the indexed node
    [[nodiscard]] NodeType type(const ind_t i) const {
      return nodes_[i].first;
    }
    //! Return whether the indexed node is a CubeNode
    [[nodiscard]] bool is_cube(const ind_t i) const {return NodeType::cube == nodes_[i].first;}
    //! Return whether the indexed node is a PolyNode
    [[nodiscard]] bool is_poly(const ind_t i) const {return NodeType::poly == nodes_[i].first;}
    //! Return whether the indexed node is a NullNode
    [[nodiscard]] bool is_null(const ind_t i) const {
      const auto & f = nodes_[i].first;
      return NodeType::null == f || NodeType::assumed_null == f || NodeType::found_null == f;
    }
    [[nodiscard]] bool is_assumed_null(const ind_t i) const {
      return NodeType::assumed_null == nodes_[i].first;
    }
    [[nodiscard]] bool is_found_null(const ind_t i) const {
      return NodeType::found_null == nodes_[i].first;
    }
    //! Return the CubeNode at index i
    [[nodiscard]] const CubeNode& cube_at(const ind_t i) const {
      return cube_nodes_[nodes_[i].second];
    }
    //! Return the PolyNode at index i
    [[nodiscard]] const PolyNode& poly_at(const ind_t i) const {
      return poly_nodes_[nodes_[i].second];
    }
    //! Return the number of vertices indexed by the node at index i
    [[nodiscard]] ind_t vertex_count(const ind_t i) const {
      switch (nodes_[i].first){
        case NodeType::cube:
        return cube_nodes_[nodes_[i].second].vertex_count();
        case NodeType::poly:
        return poly_nodes_[nodes_[i].second].vertex_count();
        default:
        return 0;
      }
    }
    //! Return the unique vertex indices in the node at index i
    [[nodiscard]] std::vector<ind_t> vertices(const ind_t i) const{
      switch (nodes_[i].first){
        case NodeType::cube:
        return cube_nodes_[nodes_[i].second].vertices();
        case NodeType::poly:
        return poly_nodes_[nodes_[i].second].vertices();
        default:
        return {};
      }
    }
    //! Return the vertex indices for the triangulated tetrahedra held by the node at index i
    [[nodiscard]] std::vector<std::array<ind_t,4>> vertices_per_tetrahedron(const ind_t i) const{
      if (nodes_[i].first == NodeType::poly)
        return poly_nodes_[nodes_[i].second].vertices_per_tetrahedron();
      return {};
    }
    /*! Find the minimum number of vertex indices and their linear interpolation weights

    \param      i  the indexed node to interogate
    \param      v  all vertex positions of the PolyhedronTrellis
    \param      x  the interpolation point
    \param[out] iw the vertex indices and their linear interpolation weights
    \returns whether node at index `i` contains sufficient information to
             allow linear interpolation at `x`

    \see CubeNode::indices_weights, PolyNode::indices_weights
    */
    template<class T, class I, template<class> class A>
    std::enable_if_t<isArray<T,A>, bool>
    indices_weights(const ind_t i, const A<T>& v, const A<T>& x,
                    std::vector<std::pair<I,T>>& iw, const bool sc=false) const {
      switch (nodes_[i].first){
        case NodeType::cube:
        return cube_nodes_[nodes_[i].second].indices_weights(v,x,iw,sc);
        case NodeType::poly:
        return poly_nodes_[nodes_[i].second].indices_weights(v,x,iw,sc);
        case NodeType::null: {
        std::ostringstream oss;
        oss << "Attempting to access null node at index " << i;
        oss << " for point " << x.to_string(0);
        throw std::logic_error(oss.str());
        }
        case NodeType::assumed_null: {
        std::ostringstream oss;
        oss << "Attempting to access assumed-null node at index " << i;
        oss << " for point " << x.to_string(0);
        throw std::logic_error(oss.str());
        }
        case NodeType::found_null: {
        std::ostringstream oss;
        oss << "Attempting to access found-null node at index " << i;
        oss << " for point " << x.to_string(0);
        throw std::logic_error(oss.str());
        }
        default:
        return false;
      }
    }
    /*! Find the total interpolable volume of the node at index i

    \param verts all vertex positions of the PolyhedronTrellis
    \param i     the indexed node to interogate
    */
    template<class T, template<class> class A>
    [[nodiscard]] std::enable_if_t<isArray<T,A>, T>
    volume(const A<T>& verts, const ind_t i) const {
      switch (nodes_[i].first){
        case NodeType::cube:
        return cube_nodes_[nodes_[i].second].volume(verts);
        case NodeType::poly:
        return poly_nodes_[nodes_[i].second].volume(verts);
        default:
        return 0.;
      }
    }
    [[nodiscard]] NodeType node_type(const ind_t i) const {
      return nodes_[i].first;
    }
    [[nodiscard]] std::string node_type_string(const ind_t i) const {
      switch (nodes_[i].first){
      case NodeType::cube: return "cube";
      case NodeType::poly: return "polyhedron";
      case NodeType::null: return "null";
      case NodeType::assumed_null: return "assumed null";
      case NodeType::found_null: return "found null";
      }
      return "unknown node type";
    }
#ifdef USE_HIGHFIVE
  public:
    //! Write to an HDF file
    template<class R>
		std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, bool>
    to_hdf(R& obj, const std::string& name) const {
      if (obj.exist(name)) obj.unlink(name);
      auto group = obj.createGroup(name);
      std::vector<NodeType> nt;
      std::vector<ind_t> idx;
      for (const auto & [type, index]: nodes_){
        nt.push_back(type);
        idx.push_back(index);
      }
      group.createDataSet("types", nt);
      group.createDataSet("index", idx);
      auto cube = group.createGroup("cube_nodes");
      auto poly = group.createGroup("poly_nodes");
      cube.createAttribute("length", cube_nodes_.size());
      poly.createAttribute("length", poly_nodes_.size());
      size_t i{0};
      for (const auto & c: cube_nodes_) c.to_hdf(cube, std::to_string(i++));
      i = 0;
      for (const auto & p: poly_nodes_) p.to_hdf(poly, std::to_string(i++));
      return true;
    }
    //! Read from an HDF file
    template<class R>
    static std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, NodeContainer>
    from_hdf(R& obj, const std::string& name){
      auto group = obj.getGroup(name);
      std::vector<NodeType> nt;
      std::vector<ind_t> idx;
      group.getDataSet("types").read(nt);
      group.getDataSet("index").read(idx);
      if (nt.size() != idx.size())
        throw std::runtime_error("Error reading Node type-index pairs");
      std::vector<std::pair<NodeType, ind_t>> n;
      size_t nc{0}, np{0};
      for (size_t i=0; i<nt.size(); ++i){
        if (nt[i] == NodeType::cube) ++nc;
        if (nt[i] == NodeType::poly) ++np;
        n.emplace_back(nt[i], idx[i]);
      }
      std::vector<CubeNode> cubes(nc);
      std::vector<PolyNode> polys(np);
      auto cg = group.getGroup("cube_nodes");
      size_t res;
      cg.getAttribute("length").read(res);
      if (res != nc) throw std::runtime_error("Error with cube node count");
      for (size_t i=0; i<nc; ++i) cubes[i] = CubeNode::from_hdf(cg, std::to_string(i));
      auto pg = group.getGroup("poly_nodes");
      pg.getAttribute("length").read(res);
      if (res != np) throw std::runtime_error("Error with poly node count");
      for (size_t i=0; i<np; ++i) polys[i] = PolyNode::from_hdf(pg, std::to_string(i));
      return {std::move(n), std::move(cubes), std::move(polys)};
    }
#endif // USE_HIGHFIVE
  };

}
#endif
