/* This file is part of brille.

Copyright © 2019,2020 Greg Tucker <greg.tucker@stfc.ac.uk>

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

#ifndef BRILLE_TRELLIS_HPP_
#define BRILLE_TRELLIS_HPP_
/*! \file
    \author Greg Tucker
    \brief A class holding a hybrid grid of cuboid and triangulated tetrahedral
           cells and data for interpolation
*/
// #include <vector>
// #include <array>
#include <queue>
// #include <tuple>
// #include <mutex>
#include <condition_variable>
#include <atomic>
// #include <algorithm>
#include <functional>
// #include <omp.h>
// #include "array.hpp"
// #include "array2.hpp"
// #include "array_latvec.hpp" // defines bArray
#include "polyhedron.hpp"
// #include "utilities.hpp"
// #include "debug.hpp"
#include "triangulation_simple.hpp"
#include "interpolatordual.hpp"
// #include "permutation.hpp"
// #include "approx.hpp"
namespace brille {

/*
  Storing the lowest bin boundary (zero[3]), constant difference (step[3]),
  and number of bins (size[3]) and directly calculating the bin for a given x
  is a possible solution, but one which is no faster than storing bin bounds.
  Since storing the boundaries can enable non-uniform bins this seems like
  the better long-term solution.
*/
// template<class T, class I>
// static I find_bin(const T zero, const T step, const I size, const T x){
//   // we want to find i such that 0 <= i*step < x-zero < (i+1)*step
//   I i = x > zero ?  static_cast<I>(std::floor((x-zero)/static_cast<T>(step))) : 0;
//   return i < size ? i : size-1;
// }
// template<class T, class I>
// static int on_boundary(const T zero, const T step, const I size, const T x, const I i){
//   // if x is infinitesimally smaller than zero+step*(i+1)
//   if (i+1<size && brille::approx::scalar(zero+step*(i+1),x)) return  1;
//   // if x is infinitesimally larger than zero+step*i
//   if (i  >0    && brille::approx::scalar(zero+step*(i  ),x)) return -1;
//   return 0;
// }
template<class T>
static size_t find_bin(const std::vector<T>& bin_edges, const T x){
  auto found_at = std::find_if(bin_edges.begin(), bin_edges.end(), [x](double b){return b>x;});
  size_t d = std::distance(bin_edges.begin(), found_at);
  // if unfound make sure we return off-the-edge on the correct side:
  if (d > bin_edges.size()-1 && x < bin_edges.front()) d = 0;
  return d>0 ? d-1 : d;
}
template<class T>
static int on_boundary(const std::vector<T>& bin_edges, const T x, const size_t i){
  // if (i==0) then above d was *either* 0 or 1, otherwise d = i + 1;
  // if (i==0) we can't go lower in either case, so no problem.
  if (i+2<bin_edges.size() && brille::approx::scalar(bin_edges[i+1],x)) return  1;
  if (i  >0                && brille::approx::scalar(bin_edges[i  ],x)) return -1;
  return 0;
}

/*! \brief An enumeration to differentiate betwee Node types

\see NullNode, CubeNode, PolyNode
*/
enum class NodeType {null, cube, poly};

/*! \brief A base class for the differentiation of Node types

The NullNode is not within the domain of the PolyhedronTrellis because it has
a null intersection with the bounding Polyhedron.
*/
class NullNode{
public:
  //! Implicit construction of an empty NullNode
  NullNode() {}
  //! Deconstruction of a NullNode
  virtual ~NullNode() = default;
  //! Return the type of this Node
  virtual NodeType type(void) const {return NodeType::null;}
  //! Return the number of vertices this Node indexes
  virtual ind_t vertex_count() const {return 0u;}
  //! Return the vertex indices of this Node
  virtual std::vector<ind_t> vertices(void) const {return std::vector<ind_t>();}
  //! Return the triangulated tetrahedra indices of this Node
  virtual std::vector<std::array<ind_t,4>> vertices_per_tetrahedron(void) const {return std::vector<std::array<ind_t,4>>();}
  //! Return the indices required and their weights for linear interpolation at a point
  bool indices_weights(const bArray<double>&, const bArray<double>&, std::vector<std::pair<ind_t,double>>&) const {return false;}
  //! Return the volume of this Node
  double volume(const bArray<double>&) const {return 0.;}
};
/*! \brief A Node fully within the domain of the PolyhedronTrellis

The eight vertices of this cuboid Node are all within the volume of the
convex Polyhedron bounding the domain of the Polyhedrontrellis.
It contains an ordered list of its vertex indices within the full list of
all PolyhedronTrellis vertices.
*/
class CubeNode: public NullNode {
  /*!< \brief The eight vertex indices of the Node

  The vertex order is critical for the CubeNode and must be:
    (000), (100), (110), (010), (101), (001), (011), (111)
  */
  std::array<ind_t, 8> vertex_indices;
public:
  //! Implicit construction of a CubeNode with no volume
  CubeNode(): vertex_indices({{0,0,0,0,0,0,0,0}}) {}
  //! Construct from an array
  CubeNode(const std::array<ind_t,8>& vi): vertex_indices(vi) {}
  //! Construct from a vector with 8 elements
  CubeNode(const std::vector<ind_t>& vi): vertex_indices({{0,0,0,0,0,0,0,0}}) {
    if (vi.size() != 8) throw std::logic_error("CubeNode objects take 8 indices.");
    for (ind_t i=0; i<8u; ++i) vertex_indices[i] = vi[i];
  }
  //! Return the number of vertices this Node indexes
  ind_t vertex_count() const { return 8u;}
  //! Return the vertex indices of this Node
  std::vector<ind_t> vertices(void) const {
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
  bool indices_weights(
    const bArray<double>& vertices,
    const bArray<double>& x,
    std::vector<std::pair<ind_t,double>>& iw
  ) const {
    // The CubeNode object contains the indices into `vertices` necessary to find
    // the 8 corners of the cube. Those indices should be ordered
    // (000) (100) (110) (010) (101) (001) (011) (111)
    // so that vertex_indices[i] and vertex_indices[7-i] are connected by a body diagonal
    auto node_verts = vertices.extract(vertex_indices);
    double node_volume = abs(node_verts.view(0)-node_verts.view(7)).prod(1)[0];
    auto w = abs(x - node_verts).prod(1)/node_volume; // the normalised volume of each sub-parallelpiped
    // If any normalised weights are greater than 1+eps() the point isn't in this node
    iw.clear();
    if (w.any(brille::cmp::gt, 1.)) return false;
    auto needed = w.is(brille::cmp::gt, 0.);
    for (int i=0; i<8; ++i) if (needed[i]) {
      // the weight corresponds to the vertex opposite the one used to find the partial volume
      iw.push_back(std::make_pair(vertex_indices[7-i],w[i]));
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
  double volume(const bArray<double>& vertices) const {
    // The CubeNode object contains the indices into `vertices` necessary to find
    // the 8 corners of the cube. Those indices should be ordered
    // (000) (100) (110) (010) (101) (001) (011) (111)
    // so that vertex_indices[i] and vertex_indices[7-i] are connected by a body diagonal
    return abs(vertices.view(vertex_indices[0])-vertices.view(vertex_indices[7])).prod(1)[0];
  }
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
  std::vector<std::array<ind_t,4>> vi_t;  //!< vertex indices per triangulated tetrahedron
  std::vector<std::array<double,4>> ci_t; //!< circumsphere information per triangulated tetrahedra
  std::vector<double> vol_t;              //!< volume per triangulated tetrahedra
public:
  //! empty implicit constructor
  PolyNode() {};
  // actually constructing the tetrahedra from, e.g., a Polyhedron object will
  // need to be done elsewhere
  //! Construct with all parameters defined
  PolyNode(
    const std::vector<std::array<ind_t,4>>& vit,
    const std::vector<std::array<double,4>>& cit,
    const std::vector<double>& volt
  ): vi_t(vit), ci_t(cit), vol_t(volt) {}
  //! Return the number of triangulated tetrahedra in the PolyNode
  ind_t tetrahedra_count() const {return static_cast<ind_t>(vi_t.size());}
  //! Return the number of unique vertex indices in the triangulated tetrahedra
  ind_t vertex_count() const { return static_cast<ind_t>(this->vertices().size());}
  //! Return the unique vertex indices from all triangulated tetrahedra
  std::vector<ind_t> vertices(void) const {
    std::vector<ind_t> out;
    for (auto tet: vi_t) for (auto idx: tet)
    if (std::find(out.begin(), out.end(), idx)==out.end()) out.push_back(idx);
    return out;
  }
  //! Return the vertex indices for each triangulated tetrahedra
  std::vector<std::array<ind_t,4>> vertices_per_tetrahedron(void) const {return vi_t;}
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
  bool indices_weights(
    const bArray<double>& vertices,
    const bArray<double>& x,
    std::vector<std::pair<ind_t,double>>& iw
  ) const {
    iw.clear();
    std::array<double,4> w{{0,0,0,0}};
    for (ind_t i=0; i<vi_t.size(); ++i)
    if (this->tetrahedra_contains(i, vertices, x, w)){
      for (int j=0; j<4; ++j) if (!brille::approx::scalar(w[j],0.))
        iw.push_back(std::make_pair(vi_t[i][j],w[j]));
      return true;
    }
    return false;
  }
  //! Return the total triangulated volume of the PolyNode
  double volume(const bArray<double>&) const {
    return std::accumulate(vol_t.begin(), vol_t.end(), 0.);
  }
private:
  bool tetrahedra_contains(
    const ind_t t,
    const bArray<double>& v,
    const bArray<double>& x,
    std::array<double,4>& w
  ) const {
    if (!this->tetrahedra_might_contain(t,x)) return false;
    double vol6 = vol_t[t]*6.0;
    w[0] = orient3d( x.ptr(0),           v.ptr(vi_t[t][1u]), v.ptr(vi_t[t][2u]), v.ptr(vi_t[t][3u]) )/vol6;
    w[1] = orient3d( v.ptr(vi_t[t][0u]), x.ptr(0),           v.ptr(vi_t[t][2u]), v.ptr(vi_t[t][3u]) )/vol6;
    w[2] = orient3d( v.ptr(vi_t[t][0u]), v.ptr(vi_t[t][1u]), x.ptr(0),           v.ptr(vi_t[t][3u]) )/vol6;
    w[3] = orient3d( v.ptr(vi_t[t][0u]), v.ptr(vi_t[t][1u]), v.ptr(vi_t[t][2u]), x.ptr(0)           )/vol6;
    if (std::any_of(w.begin(), w.end(), [](double z){return z < 0. && !brille::approx::scalar(z, 0.);}))
      return false;
    return true;
  }
  bool tetrahedra_might_contain(
    const ind_t t,
    const bArray<double>& x
  ) const {
    // find the vector from the circumsphere centre to x:
    double v[3];
    v[0]=ci_t[t][0]-x.val(0,0);
    v[1]=ci_t[t][1]-x.val(0,1);
    v[2]=ci_t[t][2]-x.val(0,2);
    // compute the squared length of v, and the circumsphere radius squared
    double d2{0}, r2 = ci_t[t][3]*ci_t[t][3];
    for (int i=0; i<3; ++i) d2 += v[i]*v[i];
    // if the squared distance is no greater than the squared radius, x might be inside the tetrahedra
    return d2 < r2 || brille::approx::scalar(d2, r2);
  }
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
private:
  std::vector<std::pair<NodeType,ind_t>> nodes_;
  std::vector<CubeNode> cube_nodes_;
  std::vector<PolyNode> poly_nodes_;
public:
  //! Return the total number of index nodes
  size_t size(void) const {return nodes_.size();}
  //! Return the number of indexed CubeNode objects
  size_t cube_count() const {
    return std::count_if(nodes_.begin(),nodes_.end(),[](std::pair<NodeType,ind_t> n){return NodeType::cube == n.first;});
  }
  //! Return the nuber of indexed PolyNode objects
  size_t poly_count() const {
    return std::count_if(nodes_.begin(),nodes_.end(),[](std::pair<NodeType,ind_t> n){return NodeType::poly == n.first;});
  }
  //! Return the nubmer of indexed NullNode objects
  size_t null_count() const {
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
  //! Return the NodeType of the indexed node
  NodeType type(const ind_t i) const {
    return nodes_[i].first;
  }
  //! Return whether the indexed node is a CubeNode
  bool is_cube(const ind_t i) const {return NodeType::cube == nodes_[i].first;}
  //! Return whether the indexed node is a PolyNode
  bool is_poly(const ind_t i) const {return NodeType::poly == nodes_[i].first;}
  //! Return whether the indexed node is a NullNode
  bool is_null(const ind_t i) const {return NodeType::null == nodes_[i].first;}
  //! Return the CubeNode at index i
  const CubeNode& cube_at(const ind_t i) const {
    return cube_nodes_[nodes_[i].second];
  }
  //! Return the PolyNode at index i
  const PolyNode& poly_at(const ind_t i) const {
    return poly_nodes_[nodes_[i].second];
  }
  //! Return the number of vertices indexed by the node at index i
  ind_t vertex_count(const ind_t i) const {
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
  std::vector<ind_t> vertices(const ind_t i) const{
    switch (nodes_[i].first){
      case NodeType::cube:
      return cube_nodes_[nodes_[i].second].vertices();
      case NodeType::poly:
      return poly_nodes_[nodes_[i].second].vertices();
      default:
      return std::vector<ind_t>();
    }
  }
  //! Return the vertex indices for the triangulated tetrahedra held by the node at index i
  std::vector<std::array<ind_t,4>> vertices_per_tetrahedron(const ind_t i) const{
    if (nodes_[i].first == NodeType::poly)
      return poly_nodes_[nodes_[i].second].vertices_per_tetrahedron();
    return std::vector<std::array<ind_t,4>>();
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
  bool indices_weights(const ind_t i, const bArray<double>& v, const bArray<double>& x, std::vector<std::pair<ind_t,double>>& iw) const{
    switch (nodes_[i].first){
      case NodeType::cube:
      return cube_nodes_[nodes_[i].second].indices_weights(v,x,iw);
      case NodeType::poly:
      return poly_nodes_[nodes_[i].second].indices_weights(v,x,iw);
      case NodeType::null:
        throw std::logic_error("attempting to access null node!");
      default:
      return false;
    }
  }
  /*! Find the total interpolable volume of the node at index i

  \param verts all vertex positions of the PolyhedronTrellis
  \param i     the indexed node to interogate
  */
  double volume(const bArray<double>& verts, const ind_t i) const {
    switch (nodes_[i].first){
      case NodeType::cube:
      return cube_nodes_[nodes_[i].second].volume(verts);
      case NodeType::poly:
      return poly_nodes_[nodes_[i].second].volume(verts);
      default:
      return 0.;
    }
  }
};

/*! \brief A class implementing a hybrid Cartesian and n-simplex grid in 3 dimensions

The PolyhedronTrellis has a Polyhedron bounded domain over which it can linearly
interpolate arbitrary data.

For quick location of the vertices required for linear interpolation at an
arbitrary point the PolyhedronTrellis defines a Cartesian grid within the
bounding box of the Polyhedron, with intersection points defining a 'trellis'
and each set of eight intersections bounding a 'node'.
- If a node is wholey within the bounding Polyhedron then its vertices can be
  used directly for linear interpolation within that node, and it is a CubeNode.
- If a node intersects with the bounding Polyhedron surface then its
  intersection with that Polyhedron has less volume than the node, and it is a
  PolyNode. The vertices of a PolyNode are triangulated for use in linear
  interpolation.
- If a node does not intersect with the bounding Polyhedron then its vertices
  can never be used to interpolate within the domain of the PolyhedronTrellis
  and it is a NullNode.

The node containing an interpolation point can be calculated directly from the
its position and the spacing of the trellis. If that node is a CubeNode the
vertices necessary to perform the interpolation are trivial to determine.
If the node is a PolyNode a checks must be made to determine which of its
triangulated tetrahedra contains the point but, as the number of tetrahedra is
typically small, this is a relatively fast process; with the tetrahedra found
the vertices and their weights required for linear interpolation are again
trivial to determine.
*/
template<typename T, typename R>
class PolyhedronTrellis{
public:
  using data_t = DualInterpolator<T,R>; //!< the container which holds the data to interpolate and performs the interpolation
  using vert_t = bArray<double>;        //!< the container which holds the vertex positions of the PolyhedronTrellis
private:
  Polyhedron polyhedron_; //!< the Polyhedron bounding the domain of the PolyhedronTrellis
  data_t data_;           //!< data for interpolation stored for each indexed vertex of the PolyhedronT0rellis
  vert_t vertices_;       //!< the indexed vertices of the PolyhedronTrellis
  NodeContainer nodes_;   //!< the nodes of the trellis, indexing the vertices of the PolyhedronTrellis
  std::array<std::vector<double>,3> boundaries_; //!< The coordinates of the trellis intersections, which bound the nodes
public:
  /*! \brief Construct from a bounding Polyhedron

  \param polyhedron         the boundary of the PolyhedronTrellis domain
  \param max_volume         maximum node volume in the same units as the
                            Polyhedron volume
  \param always_triangulate control whether nodes fully within the domain of the
                            PolyhedronTrellis are CubeNode (true) or PolyNode
                            (false) objects
  */
  explicit PolyhedronTrellis(const Polyhedron& polyhedron, const double max_volume, const bool always_triangulate=false);
  // explicit PolyhedronTrellis(const Polyhedron& polyhedron, const double max_volume){
  //   this->construct(polyhedron, max_volume);
  // }
  // PolyhedronTrellis(const Polyhedron& polyhedron, const size_t number_density){
  //   double max_volume = polyhedron.get_volume()/static_cast<double>(number_density);
  //   this->construct(polyhedron, max_volume);
  // };
  //! Explicit empty constructor
  explicit PolyhedronTrellis(): vertices_(0,3) {}
  //! Return the number of trellis intersections
  ind_t expected_vertex_count() const {
    ind_t count = 1u;
    for (ind_t i=0; i<3u; ++i) count *= boundaries_[i].size();
    return count;
  }
  //! Return the number of indexed vertices
  ind_t vertex_count() const { return static_cast<ind_t>(vertices_.size(0)); }
  //! Return a constant reference to the indexed vertex positions
  const vert_t& vertices(void) const { return vertices_; }
  //! Replace the indexed vertex positions
  const vert_t& vertices(const bArray<double>& v){
    if (v.ndim()==2 && v.size(1)==3) vertices_ = v;
    return vertices_;
  }
  //! Return the vertex positions indexed by CubeNode objects
  vert_t cube_vertices(void) const {
    std::vector<bool> keep(vertices_.size(0), false);
    for (ind_t i=0; i<nodes_.size(); ++i)
    if (nodes_.is_cube(i))
    for (auto idx: nodes_.vertices(i)) keep[idx] = true;
    return vertices_.extract(keep);
  }
  //! Return the vertex positons indexed by PolyNode objects
  vert_t poly_vertices(void) const {
    std::vector<bool> keep(vertices_.size(0), false);
    for (ind_t i=0; i<nodes_.size(); ++i)
    if (nodes_.is_poly(i))
    for (auto idx: nodes_.vertices(i)) keep[idx] = true;
    return vertices_.extract(keep);
  }
  /*! Return the vertex indices of all PolyNode tetrahedra

  \returns indices into the full indexed vertices of the PolyhedronTrellis
           as returned by `PolyhedronTrellis::vertices`
  */
  std::vector<std::array<ind_t,4>> vertices_per_tetrahedron(void) const {
    std::vector<std::array<ind_t,4>> out;
    for (ind_t i=0; i<nodes_.size(); ++i)
    if (nodes_.is_poly(i))
    for (auto tet: nodes_.vertices_per_tetrahedron(i)) out.push_back(tet);
    return out;
  }
  /*! \brief Find the vertex indices and interpolation weights for a point

  \param x the point at which linear interpolation is to be perfomed
  \return The minimal list of vertex indices and weights to perform linear
          interpolation or an empty vector if the point is not in the domain of
          the PolyhedronTrellis.
  */
  std::vector<std::pair<ind_t,double>>
  indices_weights(const bArray<double>& x) const {
    std::vector<std::pair<ind_t,double>> iw;
    if (x.ndim()!=2 && x.size(0)!=1u && x.size(1)!=3u)
      throw std::runtime_error("The indices and weights can only be found for one point at a time.");
    nodes_.indices_weights(this->node_index(x), vertices_, x, iw);
    return iw;
  }
  //! Check that the held data can be used for linear interpolation
  template<class S>
  unsigned check_before_interpolating(const bArray<S>& x) const{
    unsigned int mask = 0u;
    if (data_.size()==0)
      throw std::runtime_error("The interpolation data must be filled before interpolating.");
    if (x.ndim()!=2 || x.size(1)!=3u)
      throw std::runtime_error("Only (n,3) two-dimensional Q vectors supported in interpolating.");
    if (x.stride().back()!=1)
      throw std::runtime_error("Contiguous vectors required for interpolation.");
    return mask;
  }
  /*! \brief Perform linear interpolation at one or more points

  \param x one or more points at which to perform linear interpolation of the
           stored data.
  \returns a tuple of the interpolated eigenvalues and eigenvectors for all
           points in `x`
  */
  std::tuple<brille::Array<T>, brille::Array<R>>
  interpolate_at(const bArray<double>& x) const {
    profile_update("Single thread interpolation at ",x.size(0)," points");
    this->check_before_interpolating(x);
    auto valsh = data_.values().shape();
    auto vecsh = data_.vectors().shape();
    valsh[0] = vecsh[0] = x.size(0);
    brille::Array<T> vals_out(valsh);
    brille::Array<R> vecs_out(vecsh);
    // vals and vecs are row-ordered contiguous by default, so we can create
    // mutable data-sharing Array2 objects for use with
    // Interpolator2::interpolate_at through the constructor:
    brille::Array2<T> vals2(vals_out);
    brille::Array2<R> vecs2(vecs_out);
    for (ind_t i=0; i<x.size(0); ++i){
      verbose_update("Locating ",x.to_string(i));
      auto indwghts = this->indices_weights(x.view(i));
      if (indwghts.size()<1)
        throw std::runtime_error("Point not found in PolyhedronTrellis");
      data_.interpolate_at(indwghts, vals2, vecs2, i);
    }
    return std::make_tuple(vals_out, vecs_out);
  }
  /*! \brief Perform linear interpolation in parallel at one or more points

  \param x       one or more points at which to perform linear interpolation of
                 the stored data.
  \param threads the number of OpenMP threads to use; the return value of
                 `omp_get_max_threads()` will be used if `threads` < 1.
  \returns a tuple of the interpolated eigenvalues and eigenvectors for all
           points in `x`
  */
  std::tuple<brille::Array<T>, brille::Array<R>>
  interpolate_at(const bArray<double>& x, const int threads) const {
    this->check_before_interpolating(x);
    omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
    profile_update("Parallel interpolation at ",x.size(0)," points with ",threads," threads");
    auto valsh = data_.values().shape();
    auto vecsh = data_.vectors().shape();
    valsh[0] = vecsh[0] = x.size(0);
    // shared between threads
    brille::Array<T> vals_out(valsh);
    brille::Array<R> vecs_out(vecsh);
    // vals and vecs are row-ordered contiguous by default, so we can create
    // mutable data-sharing Array2 objects for use with
    // Interpolator2::interpolate_at through the constructor:
    brille::Array2<T> vals2(vals_out);
    brille::Array2<R> vecs2(vecs_out);
    // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
    long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
    size_t n_unfound{0};
  #pragma omp parallel for default(none) shared(x,vals2,vecs2,xsize) reduction(+:n_unfound) schedule(dynamic)
    for (long long si=0; si<xsize; ++si){
      ind_t i = brille::utils::s2u<ind_t, long long>(si);
      auto indwghts = this->indices_weights(x.view(i));
      if (indwghts.size()>0) {
        data_.interpolate_at(indwghts, vals2, vecs2, i);
      } else {
        ++n_unfound;
      }
    }
    std::runtime_error("interpolate at failed to find "+std::to_string(n_unfound)+" point"+(n_unfound>1?"s.":"."));
    return std::make_tuple(vals_out, vecs_out);
  }
  //! Return the total number of nodes within the trellis
  ind_t node_count() {
    ind_t count = 1u;
    for (ind_t i=0; i<3u; ++i) count *= static_cast<ind_t>(boundaries_[i].size()-1);
    return count;
  }
  //! Return the number of nodes along each of the three dimensions of the trellis
  std::array<ind_t,3> size() const {
    std::array<ind_t,3> s;
    for (ind_t i=0; i<3u; ++i) s[i] = static_cast<ind_t>(boundaries_[i].size()-1);
    return s;
  }
  //! Return the span between neighbouring nodes along each of the three dimensions of the trellis
  std::array<ind_t,3> span() const {
    std::array<ind_t,3> s{{1,0,0}}, sz=this->size();
    for (ind_t i=1; i<3; ++i) s[i] = sz[i-1]*s[i-1];
    return s;
  }
  /*! \brief Find the trellis node subscript index containing an arbitrary point

  \param p a point within the bounding box of the Polyhedron
  \returns the subscript index of the trellis CubeNode or PolyNode containing
          the point
  \note Since a point can be on the surface of more than one trellis node and
        two adjacent nodes do not need to be of the same type it is possible
        that a simple indexing finds a NullNode (which can not be used for
        linear interpolation). To overcome this the method searches the
        neighbouring node(s) which the point is on the surface of to find a
        non-null node subscript index.
  */
  std::array<ind_t,3> node_subscript(const bArray<double>& p) const {
    std::array<ind_t,3> sub{{0,0,0}};
    for (ind_t dim=0; dim<3u; ++dim)
      sub[dim] = static_cast<ind_t>(find_bin(boundaries_[dim], p.val(0,dim)));
    // it's possible that a subscript could go beyond the last bin in any direction!
    bool bad = !subscript_ok_and_not_null(sub);
    if (bad){
      std::array<int,3> close{{0,0,0}};
      // determine if we are close to a boundary along any of the three binning
      // directions. if we are, on_boundary returns the direction in which we
      // can safely take a step without leaving the binned region
      for (ind_t i=0; i<3; ++i)
        close[i] = on_boundary(boundaries_[i], p.val(0,i), sub[i]);
      auto num_close = std::count_if(close.begin(), close.end(), [](int a){return a!=0;});
      // check one
      std::array<ind_t,3> newsub{sub};
      if (num_close > 0) for (int i=0; i<3 && bad; ++i) if (close[i]) {
        newsub = sub;
        newsub[i] += close[i];
        bad = !subscript_ok_and_not_null(newsub);
      }
      // check two
      if (bad && num_close>1)
      for (int i=0; i<3 && bad; ++i) if (close[i])
      for (int j=0; j<3 && bad; ++j) if (close[j]) {
        newsub = sub;
        newsub[i] += close[i];
        newsub[j] += close[j];
        bad = !subscript_ok_and_not_null(newsub);
      }
      // check all three
      if (bad && num_close>2){
        newsub = sub;
        for (int i=0; i<3; ++i) newsub[i] += close[i];
        bad = !subscript_ok_and_not_null(newsub);
      }
      if (!bad) sub = newsub;
    }
    info_update_if(bad,"The node subscript ",sub," for the point ",p.to_string()," is either invalid or points to a null node!");
    return sub;
  }
  //! Find the trellis node linear index for an arbitrary point
  template <class S> ind_t node_index(const S& p) const { return this->sub2idx(this->node_subscript(p)); }

  // return a list of non-null neighbouring nodes
  /*! \brief Return a list of all non-null nodes neighbouring a trellis node

  \param idx the trellis node linear index
  \returns the trellis node linear indices of all of the neighbouring nodes
           which do not contain NullNode objects, up to 26 in total.
  */
  std::vector<ind_t> node_neighbours(const ind_t idx) const {
    std::vector<ind_t> out;
    std::array<ind_t,3> sz{this->size()}, sp{this->span()}, sub;
    sub = this->idx2sub(idx, sp);
    for (ind_t i=0; i<3; ++i) if (sub[0]+i > 0 && sub[0]+i < sz[0]+1)
    for (ind_t j=0; j<3; ++j) if (sub[1]+j > 0 && sub[1]+j < sz[1]+1)
    for (ind_t k=0; k<3; ++k) if (sub[2]+k > 0 && sub[2]+k < sz[2]+1)
    if (!(1==i&&1==j&&1==k)) {
      std::array<ind_t,3> n_sub {{sub[0]+i-1, sub[1]+j-1, sub[2]+k-1}};
      ind_t n_idx = this->sub2idx(n_sub, sp);
      if (!nodes_.is_null(n_idx)) out.push_back(n_idx);
    }
    return out;
  }

  //! Return a constant reference to the CubeNode with a given trellis node linear index
  const CubeNode& cube_node(const ind_t idx) const {
    if (idx >= nodes_.size() || !nodes_.is_cube(idx))
      throw std::runtime_error("Out-of-bounds or non-cube node");
    return nodes_.cube_at(idx);
  }
  //! Return a constant reference to the PolyNode with a given trellis node linear index
  const PolyNode& poly_node(const ind_t idx) const {
    if (idx >= nodes_.size() || !nodes_.is_poly(idx))
      throw std::runtime_error("Out-of-bounds or non-polyhedron node");
    return nodes_.poly_at(idx);
  }
  //! Return a string representation of the size of the trellis
  std::string to_string(void) const {
    std::string str = "(";
    for (auto i: this->size()) str += " " + std::to_string(i);
    str += " )";
    return str;
  }
  //! Get a constant reference to the stored data
  const data_t& data(void) const {return data_;}
  //! Replace the data stored in the object
  template<typename... A> void replace_data(A... args) {data_.replace_data(args...);}
  //! Replace the eigenvalue data stored in the object
  template<typename... A> void replace_value_data(A... args) { data_.replace_value_data(args...); }
  //! Replace the eigenvector data stored in the object
  template<typename... A> void replace_vector_data(A... args) { data_.replace_vector_data(args...); }
  //! Replace the eigenvalue data cost information stored in the object
  template<typename... A> void set_value_cost_info(A... args) { data_.set_value_cost_info(args...); }
  //! Replace the eigenvector data cost information stored in the object
  template<typename... A> void set_vector_cost_info(A... args) {data_.set_vector_cost_info(args...);}
  //! Return the number of bytes used per Q point
  size_t bytes_per_point() const {return data_.bytes_per_point(); }
  //! Calculate the Debye-Waller factor for the provided Q points and ion masses
  template<template<class> class A>
  brille::Array<double> debye_waller(const A<double>& Qpts, const std::vector<double>& Masses, const double t_K) const{
    return data_.debye_waller(Qpts,Masses,t_K);
  }
  //! Determine the sorting permutation for every connected pair of vertices in the PolyhedronTrellis
  void sort(void){ data_.sort(); }
  //! Find the total volume of all trellis nodes
  double total_node_volume() const {
    double vol{0.};
    for (ind_t i=0; i<nodes_.size(); ++i)
      vol += nodes_.volume(vertices_, i);
    return vol;
  }
private:
  bool subscript_ok_and_not_null(const std::array<ind_t,3>& sub) const {
    return this->subscript_ok(sub) && !nodes_.is_null(this->sub2idx(sub));
  }
  bool subscript_ok(const std::array<ind_t,3>& sub) const {
    return this->subscript_ok(sub, this->size());
  }
  bool subscript_ok(const std::array<ind_t,3>& sub, const std::array<ind_t,3>& sz) const {
    for (ind_t dim=0; dim<3u; ++dim) if (sub[dim]>=sz[dim]) return false;
    return true;
  }
  ind_t sub2idx(const std::array<ind_t,3>& sub) const {
    return this->sub2idx(sub, this->span());
  }
  std::array<ind_t,3> idx2sub(const ind_t idx) const {
    return this->idx2sub(idx, this->span());
  }
  ind_t sub2idx(const std::array<ind_t,3>& sub, const std::array<ind_t,3>& sp) const {
    ind_t idx=0;
    for (ind_t dim=0; dim<3u; ++dim) idx += sp[dim]*sub[dim];
    // info_update("span ",sp," gives subscript ",sub," as linear ",idx);
    return idx;
  }
  std::array<ind_t,3> idx2sub(const ind_t idx, const std::array<ind_t,3>& sp) const {
    std::array<ind_t,3> sub{{0,0,0}};
    ind_t rem{idx};
    for (ind_t dim=3u; dim--;){
      sub[dim] = rem/sp[dim];
      rem -= sub[dim]*sp[dim];
    }
    return sub;
  }
  // template<typename S> void add_node(const S& node) {nodes_.push_back(node);}

  template<typename S>
  std::vector<ind_t> which_vertices_of_node(
    const std::vector<S>& t, const S value, const ind_t idx
  ) const {
    std::vector<ind_t> out;
    for (ind_t n: nodes_.vertices(idx)) if (value == t[n]) out.push_back(n);
    return out;
  }

  template<typename S>
  std::vector<ind_t> which_node_neighbours(
    const std::vector<S>& t, const S value, const ind_t idx
  ) const {
    std::vector<ind_t> out;
    for (ind_t n: this->node_neighbours(idx)) if (value == t[n]) out.push_back(n);
    return out;
  }
  template<typename S, typename Func>
  std::vector<ind_t> which_node_neighbours(
    const std::vector<S>& t, Func ufunc, const ind_t node
  ) const {
    std::vector<ind_t> out;
    for (ind_t n: this->node_neighbours(node)) if (ufunc(t[n])) out.push_back(n);
    return out;
  }
  std::set<size_t> collect_keys();
  std::set<size_t> collect_keys_node(const ind_t);

  std::vector<std::array<double,3>> trellis_centres() const {
    std::array<std::vector<double>,3> cents;
    for (size_t i=0; i<3; ++i) for (size_t j=0; j<boundaries_[i].size()-1; ++j)
      cents[i].push_back((boundaries_[i][j]+boundaries_[i][j+1])/2);
    std::vector<std::array<double,3>> centres;
    for (auto z: cents[2]) for (auto y: cents[1]) for (auto x: cents[0])
      centres.push_back({x,y,z});
    return centres;
  }
  std::array<size_t,3> trellis_centres_span() const {
    size_t cs0{boundaries_[0].size()-1}, cs1{boundaries_[1].size()-1};
    return std::array<size_t,3>({1, cs0, cs0*cs1});
  }
  std::vector<std::array<double,3>> trellis_intersections() const {
    std::vector<std::array<double,3>> intersections;
    for (auto z: boundaries_[2]) for (auto y: boundaries_[1]) for (auto x: boundaries_[0])
      intersections.push_back({x,y,z});
    return intersections;
  }
  std::array<ind_t,3> trellis_intersections_span() const {
    size_t bs0{boundaries_[0].size()}, bs1{boundaries_[1].size()};
    return std::array<ind_t,3>({1,static_cast<ind_t>(bs0),static_cast<ind_t>(bs0*bs1)});
  }
  std::vector<std::array<ind_t,3>> trellis_local_cube_indices() const {
    /* Each node with linear index idx has a subscripted index (i,j,k)
       and is surrounded by the trellis intersections of boundaries
       (i,j,k) + { (000), (100), (110), (010), (101), (001), (011), (111)};
    */
    // the order of the cube node intersections is paramount:
    std::vector<std::array<ind_t,3>> idx{{{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},{{1,0,1}},{{0,0,1}},{{0,1,1}},{{1,1,1}}};
    return idx;
  }
  Polyhedron trellis_local_cube() const {
    std::array<double,3> min_corner, max_corner;
    for (size_t i=0; i<3; ++i){
      double cen = (boundaries_[i][0] + boundaries_[i][1])/2;
      min_corner[i] = (boundaries_[i][0] < boundaries_[i][1] ? boundaries_[i][0] : boundaries_[i][1])-cen;
      max_corner[i] = (boundaries_[i][0] > boundaries_[i][1] ? boundaries_[i][0] : boundaries_[i][1])-cen;
    }
    return polyhedron_box(min_corner, max_corner);
  }
  double trellis_node_circumsphere_radius() const {
    // this will break if the boundaries_ are ever allowed to be non-uniform
    double a = boundaries_[0][1]-boundaries_[0][0];
    double b = boundaries_[1][1]-boundaries_[1][0];
    double c = boundaries_[2][1]-boundaries_[2][0];
    return 0.5*std::sqrt(a*a+b*b+c*c);
  }
};

#include "trellis.tpp"
} // end namespace brille
#endif
