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

#ifndef BRILLE_TRELLIS_POLY_HPP_
#define BRILLE_TRELLIS_POLY_HPP_
/*! \file
    \author Greg Tucker
    \brief A class holding a hybrid grid of cuboid and triangulated tetrahedral
           cells and data for interpolation
*/
#include "approx_config.hpp"
#include "approx_float.hpp"
#include "hdf_interface.hpp"
#include "interpolatordual.hpp"
#include "polyhedron_flex.hpp"
#include "trellis_node.hpp"
#include "triangulation_poly.hpp"
#include "vertex_map_set.h"
#include "vertex_index_map.h"
#include "vertex_map_knit.h"
#include <atomic>
#include <condition_variable>
#include <filesystem>
#include <functional>
#include <queue>
#include <utility>
#include "thread_exception.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <cassert>
#endif

#include <iostream>

namespace brille::polytrellis {

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
//   if (i+1<size && brille::approx_float::scalar(zero+step*(i+1),x)) return  1;
//   // if x is infinitesimally larger than zero+step*i
//   if (i  >0    && brille::approx_float::scalar(zero+step*(i  ),x)) return -1;
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
  if (i+2<bin_edges.size() && brille::approx_float::scalar(bin_edges[i+1],x)) return  1;
  if (i  >0                && brille::approx_float::scalar(bin_edges[i  ],x)) return -1;
  return 0;
}


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
template<class DataValues, class DataVectors, class VertexComponents, template<class> class VertexType>
class PolyTrellis{
public:
  using class_t = PolyTrellis<DataValues, DataVectors, VertexComponents, VertexType>;
  using boundary_t = polyhedron::Faces;
  using data_t = DualInterpolator<DataValues, DataVectors>;
  using vert_t = VertexType<VertexComponents>;
  using nodes_t = NodeContainer;
  using knots_t = std::array<std::vector<double>, 3u>;
  using approx_t = approx_float::Config;
  using poly_t = polyhedron::Poly<VertexComponents, VertexType>;
protected:
  boundary_t boundary_; //!< the Polyhedron bounding the domain of the PolyhedronTrellis
  data_t data_;         //!< data for interpolation stored for each indexed vertex of the PolyhedronTrellis
  vert_t vertices_;     //!< the indexed vertices of the PolyhedronTrellis
  nodes_t nodes_;       //!< the nodes of the trellis, indexing the vertices of the PolyhedronTrellis
  knots_t knots_;       //!< The coordinates of the trellis intersections
  approx_t approx_;     //!< Approximate comparisons configuration
public:
  bool operator!=(const class_t& other) const {
    if (boundary_ != other.boundary_) return true;
    if (data_ != other.data_) return true;
    if (vertices_ != other.vertices_) return true;
    if (nodes_ != other.nodes_) return true;
    if (knots_ != other.knots_) return true;
    return false;
  }
  /*!\brief Construct a PolyhedronTrellis from all required information
   *
   * */
  PolyTrellis(boundary_t p, const data_t& d, const vert_t& v, nodes_t n, knots_t b, approx_t cfg)
      : boundary_(std::move(p)), data_(d), vertices_(v), nodes_(std::move(n)), knots_(std::move(b)), approx_(cfg) {}
  //
  /*! \brief Construct from a bounding Polyhedron

  \param polyhedron         the boundary of the PolyTrellis domain
  \param max_volume         maximum node volume in the same units as the
                            boundary volume
  \param always_triangulate control whether nodes fully within the domain of the
                            PolyhedronTrellis are CubeNode (true) or PolyNode
                            (false) objects
  */
  PolyTrellis(const poly_t& polyhedron,
              double max_volume,
              bool always_triangulate=false
  ): boundary_(polyhedron.faces()), vertices_(polyhedron.vertices()){
    // the approximate configuration can not be a default function parameter
    // if we want to pick-up runtime changes in the namespace object
    this->construct(polyhedron, max_volume, always_triangulate, approx_float::config);
  }
  PolyTrellis(const poly_t& polyhedron,
              double max_volume,
              bool always_triangulate,
              approx_t cfg
  ): boundary_(polyhedron.faces()), vertices_(polyhedron.vertices())
  {
    this->construct(polyhedron, max_volume, always_triangulate, cfg);
  }
  void construct(const poly_t& poly,
                 double max_volume,
                 bool always_triangulate,
                 approx_t cfg);
private:
//  std::tuple<std::map<size_t, poly_t>, std::vector<std::vector<ind_t>>, std::vector<ind_t>, VertexType<VertexComponents>, ind_t>
//  part_one(const poly_t&, const VertexType<VertexComponents>&,
//           std::vector<NodeType>&, bool, VertexComponents, int);
//  std::tuple<std::map<size_t, poly_t>, std::vector<std::vector<ind_t>>>
//  part_one(const poly_t&, VertexMapSet<VertexComponents, VertexType>&,
//           std::vector<NodeType>&, bool, VertexComponents, int);
  std::tuple<
      std::map<size_t, poly_t>,
      VertexMapSet<VertexComponents, VertexType>,
      VertexIndexMap
      >
  part_one(const poly_t&, const VertexType<VertexComponents>&,
           std::vector<NodeType>&, bool, VertexComponents, int);

  void
  part_two(const std::map<size_t, poly_t>&, const std::vector<NodeType>&,
           ind_t, const VertexIndexMap&, VertexComponents, int);
public:
  //! Explicit empty constructor
  explicit PolyTrellis(): vertices_(0,3) {}
  //! Return the number of trellis intersections
  [[nodiscard]] ind_t expected_vertex_count() const {
    ind_t count = 1u;
    for (ind_t i=0; i<3u; ++i) count *= knots_[i].size();
    return count;
  }
  //! Return the number of indexed vertices
  [[nodiscard]] ind_t vertex_count() const { return static_cast<ind_t>(vertices_.size(0)); }
  //! Return a constant reference to the indexed vertex positions
  [[nodiscard]] const vert_t& vertices() const { return vertices_; }
  //! Replace the indexed vertex positions
  const vert_t& vertices(const bArray<double>& v){
    if (v.ndim()==2 && v.size(1)==3) vertices_ = v;
    return vertices_;
  }
  //! Return the vertex positions indexed by CubeNode objects
  [[nodiscard]] vert_t cube_vertices() const {
    std::vector<bool> keep(vertices_.size(0), false);
    for (ind_t i=0; i<nodes_.size(); ++i)
      if (nodes_.is_cube(i))
        for (auto idx: nodes_.vertices(i)) keep[idx] = true;
    return vertices_.extract(keep);
  }
  //! Return the vertex positons indexed by PolyNode objects
  [[nodiscard]] vert_t poly_vertices() const {
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
  [[nodiscard]] std::vector<std::array<ind_t,4>> vertices_per_tetrahedron() const {
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
  [[nodiscard]] std::vector<std::pair<ind_t,VertexComponents>>
  indices_weights(const vert_t& x) const {
    std::vector<std::pair<ind_t,VertexComponents>> iw{};
    if (x.ndim()!=2 && x.size(0)!=1u && x.size(1)!=3u)
      throw std::runtime_error("The indices and weights can only be found for one point at a time.");
    // if node_index does not throw an error then the indicated node *should*
    // contain the interpolation point (but rounding might be a problem)
    const bool should_contain = true;
    nodes_.indices_weights(this->node_index(x), vertices_, x, iw, should_contain);
    return iw;
  }
  //! Check that the held data can be used for linear interpolation
  unsigned check_before_interpolating(const vert_t& x) const{
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
  std::tuple<typename data_t::value_out_t, typename data_t::vector_out_t>
  interpolate_at(const vert_t& x) const {
    profile_update("Single thread interpolation at ",x.size(0)," points");
    this->check_before_interpolating(x);
    auto valsh = data_.values().shape();
    auto vecsh = data_.vectors().shape();
    valsh[0] = vecsh[0] = x.size(0);
    typename data_t::value_out_t vals_out(valsh);
    typename data_t::vector_out_t  vecs_out(vecsh);
    // vals and vecs are row-ordered contiguous by default, so we can create
    // mutable data-sharing Array2 objects for use with
    // Interpolator2::interpolate_at through the constructor:
    typename data_t::value_in_t vals2(vals_out);
    typename data_t::vector_in_t vecs2(vecs_out);
    for (ind_t i=0; i<x.size(0); ++i){
      verbose_update("Locating ",x.to_string(i));
      auto indwghts = this->indices_weights(x.view(i));
      if (indwghts.size()<1){
        std::string msg = "The point " + x.to_string(i) + " was not found";
        msg += " in the PolyhedronTrellis";
        throw std::runtime_error(msg);
      }
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
  std::tuple<typename data_t::value_out_t, typename data_t::vector_out_t>
  interpolate_at(const vert_t& x, const int threads) const {
    this->check_before_interpolating(x);
    profile_update("Parallel interpolation at ",x.size(0)," points with ",threads," threads");
    auto valsh = data_.values().shape();
    auto vecsh = data_.vectors().shape();
    valsh[0] = vecsh[0] = x.size(0);
    // shared between threads
    typename data_t::value_out_t vals_out(valsh);
    typename data_t::vector_out_t  vecs_out(vecsh);
    // vals and vecs are row-ordered contiguous by default, so we can create
    // mutable data-sharing Array2 objects for use with
    // Interpolator2::interpolate_at through the constructor:
    typename data_t::value_in_t vals2(vals_out);
    typename data_t::vector_in_t vecs2(vecs_out);
    // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
    size_t missing{0};
    ThreadException thread_ex;

    const auto pool = ThreadPool::getInstance();
    if (threads > 0) pool->resize(threads); else pool->resize();
    const auto workers = pool->size();
    auto task = [&](const size_t worker) {
      auto [f, l] = thread_slice(x.size(0), workers, worker);
      return [&,first=f,last=l]() {
        for (size_t i=first; i<last; ++i) {
          thread_ex.run([&]{
            if (auto i_w = indices_weights(x.view(i)); i_w.size()>0) {
              data_.interpolate_at(i_w, vals2, vecs2, i);
            } else {
              ++missing;
            }
          });
        }
      };
    };
    for (size_t thread=0; thread<workers; ++thread) pool->enqueue(task(thread));
    pool->wait();

    thread_ex.rethrow(); // only throws if error(s) were caught
    if (missing){
      std::ostringstream oss;
      oss << "interpolate_at failed to find " << missing << " point" << (missing > 1 ? "s." : ".");
      throw std::runtime_error(oss.str());
    }
    return std::make_tuple(vals_out, vecs_out);
  }

  //! Return the total number of nodes within the trellis
  ind_t node_count() {
    ind_t count = 1u;
    for (ind_t i=0; i<3u; ++i) count *= static_cast<ind_t>(knots_[i].size()-1);
    return count;
  }
  //! Return the number of nodes along each of the three dimensions of the trellis
  [[nodiscard]] std::array<ind_t,3> size() const {
    std::array<ind_t,3> s{};
    for (ind_t i=0; i<3u; ++i) s[i] = static_cast<ind_t>(knots_[i].size()-1);
    return s;
  }
  //! Return the span between neighbouring nodes along each of the three dimensions of the trellis
  [[nodiscard]] std::array<ind_t,3> span() const {
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
  [[nodiscard]] std::array<ind_t,3> node_subscript(const vert_t& p) const {
    std::array<ind_t,3> sub{{0,0,0}};
    auto pos = get_xyz(p.view(0)).to_std(); // the knots are in the absolute xyz frame
    for (ind_t dim=0; dim<3u; ++dim)
      sub[dim] = static_cast<ind_t>(find_bin(knots_[dim], pos[dim]));
    // it's possible that a subscript could go beyond the last bin in any direction!
    bool bad = !subscript_ok_and_not_null(sub);
    if (bad){
      std::array<int,3> close{{0,0,0}};
      // determine if we are close to a boundary along any of the three binning
      // directions. if we are, on_boundary returns the direction in which we
      // can safely take a step without leaving the binned region
      for (ind_t i=0; i<3; ++i)
        close[i] = on_boundary(knots_[i], pos[i], sub[i]);
      auto num_close = std::count_if(close.begin(), close.end(), [](int a){return a!=0;});
      // check one
      std::array<ind_t,3> newsub{sub};
      if (num_close > 0) for (int i=0; i<3 && bad; ++i) if (close[i]) {
        newsub = sub;
        newsub[i] += close[i];
        bad = !subscript_ok_and_not_null(newsub);
      }
      // check two
      if (bad && num_close>1) {
        for (int i = 0; i < 3 && bad; ++i) {
          if (close[i]) {
            for (int j = 0; j < 3 && bad; ++j) {
              if (close[j]) {
                newsub = sub;
                newsub[i] += close[i];
                newsub[j] += close[j];
                bad = !subscript_ok_and_not_null(newsub);
              }
            }
          }
        }
      }
      // check all three
      if (bad && num_close>2){
        newsub = sub;
        for (int i=0; i<3; ++i) newsub[i] += close[i];
        bad = !subscript_ok_and_not_null(newsub);
      }
      if (!bad) sub = newsub;
      else {
        auto node_type = nodes_.node_type_string(sub2idx(sub));
        info_update("The node subscript ", sub, " for the point\n\t(hkl) ",
                    p.to_string(0u), "\n\t(xyz) [[", pos, "]]\n",
                    " is either invalid or a ", node_type, " node");
      }
    }
    return sub;
  }
  //! Find the trellis node linear index for an arbitrary point
  template <class S> ind_t node_index(const S& p) const { return this->sub2idx(this->node_subscript(p)); }
  //! Return the (Cube) polyhedron representing an indexed node
  poly_t subscripted_node_poly(const std::array<ind_t, 3>& ijk) const {
    auto i0 = knots_[0][ijk[0]];
    auto i1 = knots_[0][ijk[0]+1];
    auto j0 = knots_[1][ijk[1]];
    auto j1 = knots_[1][ijk[1]+1];
    auto k0 = knots_[2][ijk[2]];
    auto k1 = knots_[2][ijk[2]+1];
    std::vector<std::array<VertexComponents,3>> v{
        {i0, j0, k0}, // 000 0
        {i0, j1, k0}, // 010 1
        {i0, j1, k1}, // 011 2
        {i0, j0, k1}, // 001 3
        {i1, j0, k0}, // 100 4
        {i1, j1, k0}, // 110 5
        {i1, j1, k1}, // 111 6
        {i1, j0, k1}, // 101 7
    };
    auto faces = typename poly_t::faces_t({{3,0,4,7},{3,2,1,0},{0,1,5,4},{3,7,6,2},{7,4,5,6},{2,6,5,1}});
    auto hkl = from_xyz_like(vertices_, bArray<VertexComponents>::from_std(v));
    return poly_t(hkl, faces);
  }
  //! Return the (Cube) polyhedron representing the node containing a point
  template <class S> poly_t point_in_node_poly(const S& p) const {
    return subscripted_node_poly(node_subscript(p));
  }
  [[nodiscard]] NodeType node_at_type(const std::array<ind_t, 3>& ijk) const {
    return nodes_.node_type(sub2idx(ijk));
  }
  template <class S>
  NodeType point_in_node_type(const S & p) const {
    return node_at_type(node_subscript(p));
  }

  [[nodiscard]] std::vector<NodeType> all_node_types() const { return nodes_.all_node_types(); }

  // return a list of non-null neighbouring nodes
  /*! \brief Return a list of all non-null nodes neighbouring a trellis node

  \param idx the trellis node linear index
  \returns the trellis node linear indices of all of the neighbouring nodes
           which do not contain NullNode objects, up to 26 in total.
  */
  [[nodiscard]] std::vector<ind_t> node_neighbours(const ind_t idx) const {
    std::vector<ind_t> out;
    std::array<ind_t,3> sz{this->size()}, sp{this->span()}, sub{};
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
  [[nodiscard]] const CubeNode& cube_node(const ind_t idx) const {
    if (idx >= nodes_.size() || !nodes_.is_cube(idx))
      throw std::runtime_error("Out-of-bounds or non-cube node");
    return nodes_.cube_at(idx);
  }
  //! Return a constant reference to the PolyNode with a given trellis node linear index
  [[nodiscard]] const PolyNode& poly_node(const ind_t idx) const {
    if (idx >= nodes_.size() || !nodes_.is_poly(idx))
      throw std::runtime_error("Out-of-bounds or non-polyhedron node");
    return nodes_.poly_at(idx);
  }
  //! Return a string representation of the size of the trellis
  [[nodiscard]] std::string to_string() const {
    std::string str = "(";
    for (auto i: this->size()) str += " " + std::to_string(i);
    str += " )";
    return str;
  }
  //! Get a constant reference to the stored data
  const data_t& data() const {return data_;}
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
  [[nodiscard]] size_t bytes_per_point() const {return data_.bytes_per_point(); }
  //! Determine the sorting permutation for every connected pair of vertices in the PolyhedronTrellis
  void sort(){ data_.sort(); }
  //! Find the total volume of all trellis nodes
  [[nodiscard]] double total_node_volume() const {
    double vol{0.};
    for (ind_t i=0; i<nodes_.size(); ++i) {
      auto v = nodes_.volume(vertices_, i);
      vol += v;
    }
    return vol;
  }
  [[nodiscard]] approx_t approx_config() const {return approx_;}
private:
  [[nodiscard]] bool subscript_ok_and_not_null(const std::array<ind_t,3>& sub) const {
    return this->subscript_ok(sub) && !nodes_.is_null(this->sub2idx(sub));
  }
  [[nodiscard]] bool subscript_ok(const std::array<ind_t,3>& sub) const {
    return this->subscript_ok(sub, this->size());
  }
  [[nodiscard]] bool subscript_ok(const std::array<ind_t,3>& sub, const std::array<ind_t,3>& sz) const {
    for (ind_t dim=0; dim<3u; ++dim) if (sub[dim]>=sz[dim]) return false;
    return true;
  }
  [[nodiscard]] ind_t sub2idx(const std::array<ind_t,3>& sub) const {
    return this->sub2idx(sub, this->span());
  }
  [[nodiscard]] std::array<ind_t,3> idx2sub(const ind_t idx) const {
    return this->idx2sub(idx, this->span());
  }
  [[nodiscard]] ind_t sub2idx(const std::array<ind_t,3>& sub, const std::array<ind_t,3>& sp) const {
    ind_t idx=0;
    for (ind_t dim=0; dim<3u; ++dim) idx += sp[dim]*sub[dim];
    // info_update("span ",sp," gives subscript ",sub," as linear ",idx);
    return idx;
  }
  [[nodiscard]] std::array<ind_t,3> idx2sub(const ind_t idx, const std::array<ind_t,3>& sp) const {
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
  std::set<size_t> collect_keys_node(ind_t);

  [[nodiscard]] std::vector<std::array<double,3>> trellis_centres() const {
    knots_t cents;
    for (size_t i=0; i<3; ++i) for (size_t j=0; j<knots_[i].size()-1; ++j)
      cents[i].push_back((knots_[i][j]+knots_[i][j+1])/2);
    std::vector<std::array<double,3>> centres;
    for (auto z: cents[2]) for (auto y: cents[1]) for (auto x: cents[0])
      centres.push_back({x,y,z});
    return centres;
  }
  [[nodiscard]] std::array<size_t,3> trellis_centres_span() const {
    // TODO: This is identical to span()?
    size_t cs0{knots_[0].size()-1}, cs1{knots_[1].size()-1};
    return std::array<size_t,3>({1, cs0, cs0*cs1});
  }
  [[nodiscard]] std::vector<std::array<double,3>> trellis_intersections() const {
    std::vector<std::array<double,3>> intersections;
    for (auto z: knots_[2]) for (auto y: knots_[1]) for (auto x: knots_[0])
      intersections.push_back({x,y,z});
    return intersections;
  }
  [[nodiscard]] std::array<ind_t,3> trellis_intersections_span() const {
    size_t bs0{knots_[0].size()}, bs1{knots_[1].size()};
    return std::array<ind_t,3>({1,static_cast<ind_t>(bs0),static_cast<ind_t>(bs0*bs1)});
  }
  [[nodiscard]] std::vector<std::array<ind_t,3>> trellis_local_cube_indices() const {
    /* Each node with linear index idx has a subscripted index (i,j,k)
       and is surrounded by the trellis intersections of boundaries (the knots)
       (i,j,k) + { (000), (100), (110), (010), (101), (001), (011), (111)};
    */
    // the order of the cube node intersections is paramount:
    std::vector<std::array<ind_t,3>> idx{{{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},{{1,0,1}},{{0,0,1}},{{0,1,1}},{{1,1,1}}};
    return idx;
  }
  [[nodiscard]] double maximum_node_circumsphere_radius() const {
    // this will break if the knots_ are ever allowed to be non-uniform
    double a = knots_[0][1]-knots_[0][0];
    double b = knots_[1][1]-knots_[1][0];
    double c = knots_[2][1]-knots_[2][0];
    return 0.5*std::sqrt(a*a+b*b+c*c);
  }
  //! Return the polyhedron face indexes into trellis_intersections for an indexed node
  [[nodiscard]] boundary_t trellis_node_faces(const ind_t index){
    auto ijk = idx2sub(index); // the node subscript
    auto ks = trellis_intersections_span();
    ind_t i0 = this->sub2idx(ijk, ks);     // 000 0
    ind_t i1 = i0         + ks[1]        ; // 010 1
    ind_t i2 = i0         + ks[1] + ks[2]; // 011 2
    ind_t i3 = i0                 + ks[2]; // 001 3
    ind_t i4 = i0 + ks[0]                ; // 100 4
    ind_t i5 = i0 + ks[0] + ks[1]        ; // 110 5
    ind_t i6 = i0 + ks[0] + ks[1] + ks[2]; // 111 6
    ind_t i7 = i0 + ks[0]         + ks[2]; // 101 7
    std::vector<std::vector<ind_t>> faces{
        {i3,i0,i4,i7}, // y == 0
        {i3,i2,i1,i0}, // x == 0
        {i0,i1,i5,i4}, // z == 0
        {i3,i7,i6,i2}, // z == 1
        {i7,i4,i5,i6}, // x == 1
        {i2,i6,i5,i1}  // y ==1
    };
    return boundary_t(faces);
  }

public:
  template<class HFObject>
  std::enable_if_t<std::is_base_of_v<HighFive::Object, HFObject>, bool>
  to_hdf(HFObject& obj, const std::string& entry) const {
    auto group = overwrite_group(obj, entry);
    bool ok{true};
    ok &= boundary_.to_hdf(group, "boundary");
    ok &= data_.to_hdf(group, "data");
    ok &= vertices_.to_hdf(group, "vertices");
    ok &= nodes_.to_hdf(group, "container");
    ok &= lists_to_hdf(knots_, group, "knots");
    ok &= approx_.to_hdf(group, "approx");
    return ok;
  }
  template<class HF>
  static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, class_t>
  from_hdf(HF& obj, const std::string& entry) {
    auto group = obj.getGroup(entry);
    auto p = boundary_t::from_hdf(group, "boundary");
    auto d = data_t::from_hdf(group, "data");
    auto v = vert_t::from_hdf(group, "vertices");
    auto n = nodes_t::from_hdf(group, "container");
    auto bl = lists_from_hdf<double>(group, "knots"); // returns a std::vector<std::vector<double>>
    if (bl.size() != 3) throw std::runtime_error("Error reading boundaries from file");
    knots_t b;
    for (size_t i=0; i<3u; ++i) b[i] = bl[i];
    //      return {p, d, v, n, b};
    auto cfg = approx_t::from_hdf(group, "approx");
    return PolyTrellis(p, d, v, n, b, cfg);
  }
  [[nodiscard]] bool to_hdf(const std::string& filename, const std::string& entry, const unsigned perm=HighFive::File::OpenOrCreate) const {
    HighFive::File file(filename, perm);
    return this->to_hdf(file, entry);
  }
  static class_t from_hdf(const std::string& filename, const std::string& entry){
    HighFive::File file(filename, HighFive::File::ReadOnly);
    return class_t::from_hdf(file, entry);
  }

};


template<class I>
static void add_to_maps(const I n_points, I& n_kept, std::vector<I>& map_idx, std::vector<I>& map, const I i)
{
  if (map_idx[i] > n_points) {
    map_idx[i] = n_kept++;
    debug_update("vertex ", i, " added to kept list; total kept now ", n_kept);
  }
  map.push_back(map_idx[i]);
}

template<class I, class T, template<class> class A>
static bool add_vertex(const I n_points,
                       I& n_kept,
                       std::vector<I>& map_idx,
                       const A<T>& points,
                       A<T>& extra,
                       I& n_extra,
                       const T s_tol,
                       const int d_tol,
                       const A<T>& vertex,
                       std::vector<I>& map
                       )
{
  auto equals = points.row_is(brille::cmp::eq, vertex, s_tol, s_tol, d_tol);
  auto no = equals.count();
  if (no) {
    if (no > 1) throw std::runtime_error("Too many matches to vertex");
    // modifies n_kept, map_index & map
    add_to_maps(n_points, n_kept, map_idx, map, equals.first());
    return false;
  }
  equals = extra.row_is(brille::cmp::eq, vertex, s_tol, s_tol, d_tol);
  no = equals.count(n_extra); // only count set points (protect against matching an uninitialized extra entry)
  if (no) {
    if (no > 1) throw std::runtime_error("Too many extra matches to vertex");
    map.push_back(n_points + equals.first(n_extra));
    return false;
  }
  if (extra.size(0) < n_extra + 1)
    extra.resize(2 * n_extra);
  extra.set(n_extra, vertex);
  map.push_back(n_points + n_extra++);
  return true;
}


template<class T, class R, class S, template<class> class A>
void PolyTrellis<T,R,S,A>::construct(const polyhedron::Poly<S,A>& poly,
                                  const double max_volume,
                                  const bool always_triangulate,
                                  const approx_t cfg)
{
  profile_update("Start of PolyTrellis construction");
  assert(poly.face_count() > 3 && poly.volume() > 0);

  approx_ = cfg;
  S s_tol = cfg.reciprocal<S>();
  int d_tol = cfg.digit();
  // find the extents of the polyhedron
  const auto v_hkl{poly.vertices()};
  const auto pv = get_xyz(poly.vertices());
  auto min = pv.min(0).to_std();
  auto max = pv.max(0).to_std();
  auto lengths = pv.max(0) - pv.min(0);
  lengths /= (lengths / std::cbrt(max_volume)).ceil();
  auto nl = lengths.to_std();

  // build-up the trellis intersection-point knots:
  for (int i = 0; i < 3; ++i) {
    knots_[i].reserve(static_cast<size_t>(std::ceil((max[i] - min[i]) / nl[i]) + 1));
    knots_[i].push_back(min[i]);
    while (knots_[i].back() < max[i])
      knots_[i].push_back(knots_[i].back() + nl[i]);
    debug_update("PolyTrellis has ", knots_[i].size() - 1, " bins along axis ",
                 i, ", with boundaries ", knots_[i]);
  }

  auto node_centres =
      from_xyz_like(v_hkl, from_std_like(pv, this->trellis_centres()));
  double max_dist =
      this->maximum_node_circumsphere_radius() + poly.circumsphere_radius();

  std::vector<NodeType> node_type;
  auto is_null = norm(node_centres - poly.centroid()).is(brille::cmp::gt, max_dist).to_std();
  std::transform(is_null.begin(), is_null.end(), std::back_inserter(node_type),
               [](const auto & b){return b ? NodeType::assumed_null : NodeType::cube;});

  // pull together the trellis knots
  auto all_points = from_xyz_like(v_hkl, from_std_like(pv, this->trellis_intersections()));

//  VertexMapSet vertex_set(all_points, s_tol, d_tol);

  // Go through all nodes and determine if they are null, cube, or poly
  // For cube and poly nodes, map their vertices to those in the knots
  // keeping track of 'extra' non-knot polyhedron vertices, and their number
  // the number of knots used in any node vertex mappings
  auto [stash, vertex_set, node_index_map] = part_one(poly, all_points, node_type, always_triangulate, s_tol, d_tol);

  // find the bounding polyhedra vertices in the knots or extra points
  std::vector<std::pair<MapVertexType, ind_t>> boundary_map;
  boundary_map.reserve(v_hkl.size(0));
  for (ind_t i = 0; i < v_hkl.size(0); ++i) {
    boundary_map.push_back(vertex_set.add(v_hkl.view(i)));
  }

  //auto consolidated_vertex_set = vertex_set.consolidate();
  vertices_ = vertex_set.consolidate().pristine();

  auto n_kept = vertex_set.preserved_count();
//  auto n_lost = vertex_set.pristine_count() - n_kept;

  // allocates and fills the NodeContainer
  part_two(stash, node_type, n_kept, node_index_map, s_tol, d_tol);
  // Now all non-null nodes have been populated with the indices of their vertices

  // create the Faces object too:
  // update the boundary map to account for point extraction:
  std::vector<ind_t> update_boundary_map;
  for (auto x: boundary_map){
//    update_boundary_map.push_back(idx < n_kept ? idx : idx - n_lost);
    update_boundary_map.push_back(x.first == MapVertexType::Pristine ? x.second : n_kept + x.second);
  }
  auto p_faces = poly.faces(); // Faces(std::vector<std::vector>>)
  auto pf_faces = p_faces.faces(); // std::vector<std::vector>>
  for (auto & face: pf_faces) for (auto & index: face) index = update_boundary_map[index];
  boundary_ = polyhedron::Faces(pf_faces);


  profile_update("  End of PolyTrellis Construction");
  // the data_ PermutationTable should be initialised now:
  data_.initialize_permutation_table(vertices_.size(0), this->collect_keys());
}

template<class T, class R, class S, template<class> class A>
std::tuple<
    std::map<size_t, typename PolyTrellis<T,R,S,A>::poly_t>,
    VertexMapSet<S,A>, VertexIndexMap
    >
PolyTrellis<T,R,S,A>::part_one(const poly_t& poly, const A<S>& all_points, std::vector<NodeType>& node_type, const bool always_triangulate, const S s_tol, const int d_tol) {
  profile_update("Starting PolyTrellis part_one");
  //  info_update("Using tolerance ", s_tol, " and digits ", d_tol);
  const ind_t nNodes = this->node_count();
  /*
  Find intersections corresponding to nodes that intersect with the polyhedron:
  The original approach was too simplistic as it *only* considered whether an
  intersection point is *inside* of the polyhedron. Such a criteria will not
  capture nodes where there is overlap but no vertices from the node are within
  the polyhedron.
  */
  std::vector<std::pair<VertexMapSet<S,A>, VertexIndexMap>> thread_pairs;

  auto Gamma = 0 * all_points.view(0);

//  VertexMapSet vertex_set(all_points, s_tol, d_tol);
//  auto vertex_index_map = VertexIndexMap();
  std::map<size_t, poly_t> poly_stash;
  std::mutex stash_mutex;

  debug_update_if(
    std::find(node_type.begin(), node_type.end(), NodeType::assumed_null) == node_type.end(),
    "No null nodes!");

  // Collect the cube node indexes to improve parallel loop work distribution
  size_t cubes{0};
  for (const auto nt : node_type) if (nt == NodeType::cube) ++cubes;
  std::vector<ind_t> cube_indexes;
  cube_indexes.reserve(cubes);
  for (ind_t i=0; i < nNodes; ++i) if (NodeType::cube == node_type[i]) cube_indexes.push_back(i);


  const auto pool = ThreadPool::getInstance();
  const auto workers = pool->size();
  // initialize the thread-pair data for each worker
  for (size_t worker=0; worker<workers; ++worker) {
    thread_pairs.push_back(std::make_pair(VertexMapSet<S,A>(all_points, s_tol, d_tol), VertexIndexMap()));
  }

  ThreadException thread_ex;
  auto task = [&](const size_t input_thread) {
    auto & [vms, vim] = thread_pairs[input_thread];
    auto [f, l] = thread_slice(cubes, workers, input_thread);
    // capturing 'input_thread' by reference would make all threads _share_ the same value!
    return [&,first=f,last=l,&vertex_map_set=vms,&vertex_index_map=vim]() {
      auto test0 = !always_triangulate;
      for (size_t cube=first; cube<last; ++cube) {
        auto node_index = cube_indexes[cube]; // the node index of this cube
        auto this_node_faces = trellis_node_faces(node_index);
        // this limits all_points to have the knots *first*
        auto this_node_poly = polyhedron::Poly(all_points, this_node_faces);
        auto intersection = poly.intersection(this_node_poly, s_tol, d_tol);
        /* FIXME A less-strict comparison might be useful but, at present, causes
         *       not-in-trellis errors */
        auto test1 = this_node_poly.volume() == intersection.volume();
        auto test2 = !this_node_poly.contains(Gamma)[0];
        node_type[node_index] = (test0 && test1 && test2) ? NodeType::cube : NodeType::poly;
        auto indexes = this_node_faces.indexes();
        if (NodeType::cube == node_type[node_index]) {
          // auto & node_map = thread_pairs[thread].second.get(i);
          auto & node_map = vertex_index_map.get(node_index);
          node_map.resize(indexes.size());
          for (ind_t j = 0; j < indexes.size(); j++) {
            // node_map[j] = thread_pairs[thread].first.preserve(indexes[j]);
            node_map[j] = vertex_map_set.preserve(indexes[j]);
          }
        } else {
          if (intersection.face_count() < 4 || intersection.volume() <= 0.) {
            node_type[node_index] = NodeType::found_null;
            continue; // we don't want to re-wrap the following code block
          }
          // find which of the vertices of the intersection are knots as well.
          const auto &iv{intersection.vertices()};

          // Determine which of the intersection vertices are from the pristine set
          std::map<ind_t, ind_t> i2p;
          for (auto idx : indexes) {
            if (auto ip = iv.row_is(cmp::eq, all_points.view(idx), s_tol, s_tol, d_tol); ip.count() == 1) {
              i2p.emplace(ip.first(), idx);
            }
          }
          // *SUPER IMPORTANT* Go through the vertices of the intersection IN
          // ORDER Check for the index in the already-constructed map
          // auto & node_map = thread_pairs[thread].second.get(i);
          auto & node_map = vertex_index_map.get(node_index);
          node_map.resize(iv.size(0));
          for (ind_t iv_index = 0; iv_index < iv.size(0); ++iv_index) {
            thread_ex.run([&]{
            if (auto search = i2p.find(iv_index); search != i2p.end()) {
              // node_map[iv_index] = thread_pairs[thread].first.preserve(search->second);
              node_map[iv_index] = vertex_map_set.preserve(search->second);
            } else {
					    // might-throw:
              // node_map[iv_index] = thread_pairs[thread].first.add(iv.view(iv_index), AddVertexType::Crafted);
              node_map[iv_index] = vertex_map_set.add(iv.view(iv_index), AddVertexType::Crafted);
            }
					  }); // end might-throw handler ThreadException runner
          }
        }
        // this was protected by a not-null check; but now if the node is null
        // there will have been a runtime error raised above, so there's no need
        // for a guard
        if (auto gamma_in = intersection.contains(Gamma); std::count(gamma_in.begin(), gamma_in.end(), true)) {
          auto int_vertices_are_gamma = intersection.vertices().row_is(cmp::eq, Gamma, s_tol, s_tol, d_tol);
          auto no = int_vertices_are_gamma.count();
          // if the gamma point *is* an intersection point
          // it is already in the node_index map.
          if (no == 0) {
            // thread_pairs[thread].second.append(i, thread_pairs[thread].first.origin_index());
            vertex_index_map.append(node_index, vertex_map_set.origin_index());
          }
        }
        {
          // Storing the intersection into the std::map is likely not thread safe
          std::unique_lock stash_lock(stash_mutex);
          // stash_mutex.lock();
          poly_stash.emplace(node_index, intersection);
          // stash_mutex.unlock();
        }
      }
    };
  };
  for (size_t i=0; i<workers; ++i) pool->enqueue(task(i));
  pool->wait();

  thread_ex.rethrow(); //Now handle any encountered errors if they exist

  /* Now combine the per-thread VertexMapSets and VertexIndexMaps */
  profile_update(" Start vertex maps reduction");
  auto comb = vertex_maps::parallel_reduce(thread_pairs);
  profile_update("  End of PolyTrellis part_one");
  return std::make_tuple(poly_stash, comb.first, comb.second);
}

template<class T, class R, class S, template<class> class A>
void
PolyTrellis<T,R,S,A>::part_two(
    const std::map<size_t, poly_t>& poly_stash,
    const std::vector<NodeType>& node_type,
    const ind_t n_kept,
    const VertexIndexMap& node_index_map,
    const S s_tol,
    const int d_tol
) {
  profile_update("Starting PolyTrellis part_two");
  // go through all cells again and construct the actual nodes:

  /* You might be tempted to do this in parallel, but the way NodeContainer
   * is implemented makes this impossible at present.
   *
   * To setup the NodeContainer for parallelisation we would need:
   *    - the total number of nodes (nNodes) ✅
   *    - the multiplicity of each node type (count using node_type) ✅
   *    - a pre-filling constructor for NodeType ✅
   *      + this would take the counts of each type and pre-allocate their
   *        (pointer|reference) storage ✅
   *      + it could take node_type as input and even fill-in the nodes_ array!✅
   *    - a thread-safe setter ❓ Maybe the implementation is ok
   * */
  // pre-allocate the NodeContainer storage and fill-in its mapping vector
  nodes_ = NodeContainer(node_type);
  auto count = this->node_count();

  // collect the cube and polyhedron node indexes to simplify the parallel loop
  // Collect the cube node indexes to improve parallel loop work distribution
  std::vector<ind_t> cube_indexes, poly_indexes;
  size_t cubes{0}, polys{0};
  for (const auto t: node_type){
    if (t == NodeType::cube) ++cubes;
    else if (t == NodeType::poly) ++polys;
  }
  cube_indexes.reserve(cubes);
  poly_indexes.reserve(polys);
  for (ind_t i=0; i < count; ++i){
    switch (node_type[i]){
    case NodeType::cube: cube_indexes.push_back(i); break;
    case NodeType::poly: poly_indexes.push_back(i); break;
    default: continue;
    }
  }

  profile_update("Cube and Poly node indexes collected");

  const auto pool = ThreadPool::getInstance();
  const auto workers = pool->size();

  auto cube_task = [&](const size_t thread) {
    auto [f, l] = thread_slice(cubes, workers, thread);
    return [&,first=f,last=l]() {
      for (size_t c_i=first; c_i<last; ++c_i) {
        const auto node_i = cube_indexes[c_i]; // linear index of the node to fill
        std::array<ind_t,8> fvi;
        for (ind_t j=0; j<8u; ++j) fvi[j] = node_index_map.decode(n_kept, node_i, j);
        // FIXME trellis_node_faces and CubeNode do not agree on vertex ordering!
        // CubeNode requires: [0,0,0], [1,0,0], [1,1,0], [0,1,0], [1,0,1], [0,0,1], [0,1,1], [1,1,1]
        //     faces assumes: [0,0,0], [0,1,0], [0,1,1], [0,0,1], [1,0,0], [1,1,0], [1,1,1], [1,0,1]
        // The Faces object has its vertices ordered to form a valid cube polyhedron via
        //    [3,0,4,7], [3,2,1,0], [0,1,5,4], [3,7,6,2], [7,4,5,6], [2,6,5,1]
        // and then the unique indexes are extracted in order [3,0,4,7,2,1,5,6] before being mapped into node_index_map
        // In the CubeNode vertex order indexing scheme this is [5,0,1,4,6,3,2,7] which requires the permutation
        // [1,2,6,5,3,0,4,7] to put the vertices in the 'correct' order:
        std::array<ind_t,8> cube_vert_idx{{fvi[1], fvi[2], fvi[6], fvi[5], fvi[3], fvi[0], fvi[4], fvi[7]}};
        nodes_.set(node_i, CubeNode(cube_vert_idx));
      }
    };
  };
  // Handle all polyhedron nodes (this almost certainly needs to be parallel)
  ind_t fatal_tri{0}, fatal_miss{0}, fatal_match{0}, hiccups{0};
  std::mutex tetgen_mutex;
  auto poly_task = [&](const size_t thread) {
    auto [f, l] = thread_slice(polys, workers, thread);
    return [&,first=f,last=l]() {
      for (size_t s_i=first; s_i<last; ++s_i) {
        auto i = poly_indexes[s_i];
        auto Gamma = 0 * vertices_.view(0);
        // combine the vertex indexes for preserved (kept pristine) and appended
        auto node_verts = vertices_.extract(node_index_map.decode(n_kept, i));

        // get a reference the stashed Polyhedron from before:
        const auto & this_node{poly_stash.at(i)};
        auto gin = this_node.contains(Gamma);
        bool contains_Gamma = std::count(gin.begin(), gin.end(), true) == 1;
        // Triangulate the node polyhedron into tetrahedra (the class name LQPolyTet is misleading...)
        polyhedron::LQPolyTet<S,A> tri_cut{};
        {
          // Parallel execution of tetgen somehow causes segfaults.
          std::unique_lock lock(tetgen_mutex);
          tri_cut = polyhedron::LQPolyTet(this_node, contains_Gamma);
        }
        if (tri_cut.get_vertices().size(0)<4){
          //something went wrong.
          /* A (somehow) likely culprit is that a face is missing from the cut
          cube and therefor is not a piecewise linear complex. try to re-form
          the input polyhedron and then re-triangulate.*/
          std::unique_lock lock(tetgen_mutex);
          tri_cut = polyhedron::LQPolyTet(this_node.convex_hull(), contains_Gamma);
          if (tri_cut.get_vertices().size(0)<4)
            fatal_tri += 1;
          else
            hiccups += 1;
        }
        int added_in_triangulation = tri_cut.any() ? tri_cut.count() : 0;
        // make sure we can match up the triangulated polyhedron vertices
        std::vector<ind_t> local_map;
        const auto& triverts{tri_cut.get_vertices()};
        for (ind_t j=0; j<triverts.size(0); ++j) {
          const auto trij{triverts.view(j)};
          auto known = node_verts.row_is(cmp::eq, trij, s_tol, s_tol, d_tol);
          auto no = known.count();
          if (no){
            if (no > 1){
              fatal_match += 1;
            }
            local_map.push_back(node_index_map.decode(n_kept, i, known.first()));
          } else {
            if (auto punt = vertices_.row_is(cmp::eq, trij, s_tol, s_tol, d_tol); punt.count() < 1 && added_in_triangulation > 0){
              // we've *ADDED* a new point to the triangulation, so we need to
              // add it to all points as well. hopefully this doesn't happen often
              local_map.push_back(vertices_.size(0));
              vertices_ = cat(0, vertices_, trij);
              --added_in_triangulation;
            } else {
              if (punt.count() < 1) {
                fatal_miss += 1;
              }
              local_map.push_back(punt.first());
            }
          }
        }
        // find the indices per tetrahedron, the tetrahedron circumsphere info,
        // and the tetrahedron volumes
        std::vector<std::array<ind_t,4>> idx_per_tet;
        const auto& local_ipt{tri_cut.get_vertices_per_tetrahedron()};
        for (ind_t j=0; j<local_ipt.size(0); ++j){
          std::array<ind_t,4> one_tet{0,0,0,0};
          for (ind_t k=0; k<4; ++k) one_tet[k] = local_map[local_ipt.val(j,k)];
          idx_per_tet.push_back(one_tet);
        }
        std::vector<std::array<double,4>> cci_per_tet;
        std::vector<double> vol_per_tet;
        for (ind_t j=0; j<tri_cut.number_of_tetrahedra(); ++j){
          cci_per_tet.push_back(tri_cut.circumsphere_info(j));
          vol_per_tet.push_back(tri_cut.volume(j));
        }
        if (idx_per_tet.size()<1){
          throw std::runtime_error("Triangulated node is actually Null!");
        }
        nodes_.set(i, PolyNode(idx_per_tet, cci_per_tet, vol_per_tet));
      }
    };
  };
  for (size_t i=0; i<workers; ++i) pool->enqueue(cube_task(i));
  for (size_t i=0; i<workers; ++i) pool->enqueue(poly_task(i));
  pool->wait();
  profile_update("Cube node vertices stored");

//
//   auto p_count = utils::u2s<long long>(poly_indexes.size());
// #pragma omp parallel for default(none) shared(n_kept, poly_stash, poly_indexes, p_count, node_index_map, s_tol, d_tol) reduction(+:fatal_tri, fatal_miss, fatal_match, hiccups)
//   for (long long s_i=0; s_i<p_count; ++s_i) {
//     auto i = poly_indexes[utils::s2u<ind_t>(s_i)];
//     //
//     auto Gamma = 0 * vertices_.view(0);
//     // combine the vertex indexes for preserved (kept pristine) and appended
//     auto node_verts = vertices_.extract(node_index_map.decode(n_kept, i));
//
//     // get a reference the stashed Polyhedron from before:
//     const auto & this_node{poly_stash.at(i)};
//     auto gin = this_node.contains(Gamma);
//     bool contains_Gamma = std::count(gin.begin(), gin.end(), true) == 1;
//     // Triangulate the node polyhedron into tetrahedra (the class name LQPolyTet is misleading...)
//     polyhedron::LQPolyTet<S,A> tri_cut{};
//     // this uses modified TetGen, which is now thread safe :)
//     /*
//      * ^ 2025-10-17 -- this statement appears to be wrong
//      *
//      * v For now, OMP critical sections _should_ let this work
//      */
// #pragma omp critical
//     {
//       tri_cut = polyhedron::LQPolyTet(this_node, contains_Gamma);
//     }
//     if (tri_cut.get_vertices().size(0)<4){
//       //something went wrong.
//       /* A (somehow) likely culprit is that a face is missing from the cut
//       cube and therefor is not a piecewise linear complex. try to re-form
//       the input polyhedron and then re-triangulate.*/
// #pragma omp critical
//       {
//         tri_cut = polyhedron::LQPolyTet(this_node.convex_hull(), contains_Gamma);
//       }
//       if (tri_cut.get_vertices().size(0)<4)
//         fatal_tri += 1;
//       else
//         hiccups += 1;
//     }
//     int added_in_triangulation = tri_cut.any() ? tri_cut.count() : 0;
//     // make sure we can match up the triangulated polyhedron vertices to the known ones:
//     std::vector<ind_t> local_map;
//     const auto& triverts{tri_cut.get_vertices()};
//     for (ind_t j=0; j<triverts.size(0); ++j){
//       const auto trij{triverts.view(j)};
//       auto known = node_verts.row_is(brille::cmp::eq, trij, s_tol, s_tol, d_tol);
//       auto no = known.count();
//       if (no){
//         if (no > 1){
//           // std::cout << trij << " matches\n" << cat(1, known, node_verts) << no << " times?!" << std::endl;
// //          for (const auto & x: node_index_map.get(i)) std::cout << " " << x.first << x.second;
//           // for (const auto & x: node_index_map.get(i)) std::cout << " " << x.second;
//           // std::cout << std::endl;
//           fatal_match += 1;
//         }
//         //local_map.push_back(poly_vert_idx[known.first()]);
//         local_map.push_back(node_index_map.decode(n_kept, i, known.first()));
//       } else {
//         auto punt = vertices_.row_is(brille::cmp::eq, trij, s_tol, s_tol, d_tol);
//         if (punt.count() < 1 && added_in_triangulation > 0){
//           // we've *ADDED* a new point to the triangulation, so we need to
//           // add it to all points as well. hopefully this doesn't happen often
//           local_map.push_back(vertices_.size(0));
//           vertices_ = cat(0, vertices_, trij);
//           --added_in_triangulation;
//         } else {
//           if (punt.count() < 1) {
//             fatal_miss += 1;
//           }
//           local_map.push_back(punt.first());
//         }
//       }
//     }
//
//     // find the indices per tetrahedron, the tetrahedron circumsphere info,
//     // and the tetrahedron volumes
//     std::vector<std::array<ind_t,4>> idx_per_tet;
//     const auto& local_ipt{tri_cut.get_vertices_per_tetrahedron()};
//     for (ind_t j=0; j<local_ipt.size(0); ++j){
//       std::array<ind_t,4> one_tet{0,0,0,0};
//       for (ind_t k=0; k<4; ++k) one_tet[k] = local_map[local_ipt.val(j,k)];
//       idx_per_tet.push_back(one_tet);
//     }
//     std::vector<std::array<double,4>> cci_per_tet;
//     std::vector<double> vol_per_tet;
//     for (ind_t j=0; j<tri_cut.number_of_tetrahedra(); ++j){
//       cci_per_tet.push_back(tri_cut.circumsphere_info(j));
//       vol_per_tet.push_back(tri_cut.volume(j));
//     }
//     if (idx_per_tet.size()<1){
//       throw std::runtime_error("Triangulated node is actually Null!");
//     }
//     nodes_.set(i, PolyNode(idx_per_tet, cci_per_tet, vol_per_tet));
//   }

  profile_update("Poly node vertexes stored");

  if (fatal_tri + fatal_miss + fatal_match){
    std::stringstream msg;
    if (fatal_tri) msg << fatal_tri << " Error(s) determining cut cube triangulation; ";
    if (fatal_match) msg << fatal_match << " Multiple matches of a triangulated vertex; ";
    if (fatal_miss) msg << fatal_miss << " Missing known vertex for triangulated vertex";
    throw std::runtime_error(msg.str());
  }
  debug_update_if(hiccups, "Bad vertex indexing occurred ", hiccups, " times");
  profile_update("  End of PolyTrellis part_two");
}


template<class T, class R, class S, template<class> class A>
std::set<size_t>
PolyTrellis<T,R,S,A>::collect_keys() {
  profile_update("Start of PolyTrellis permutation key collection");
  std::set<size_t> keys;
  std::mutex keys_mutex;

  const auto pool = ThreadPool::getInstance();
  const auto workers = pool->size();
  auto task = [&](const size_t worker) {
    auto [f, l] = thread_slice(nodes_.size(), workers, worker);
    return [&,first=f,last=l]() {
      for (size_t i=first; i<last; ++i) {
        if (const auto these = collect_keys_node(i); these.size()) {
          std::unique_lock lock(keys_mutex);
          keys.insert(these.begin(), these.end());
        }
      }
    };
  };
  for (size_t i=0; i<workers; ++i) pool->enqueue(task(i));
  pool->wait();
  profile_update("  End of PolyTrellis permutation key collection");
  return keys;
}

template<class T, class R, class S, template<class> class A>
std::set<size_t>
PolyTrellis<T,R,S,A>::collect_keys_node(const ind_t node){
  std::set<size_t> keys;
  if (!nodes_.is_null(node)){
    if (nodes_.is_poly(node)){
      // in the case of a polynode we must exploit the inner connectivity
      std::vector<std::array<ind_t,4>> tets = nodes_.vertices_per_tetrahedron(node);
      for (auto vt: tets){
        std::set<size_t> tmp = permutation_table_keys_from_indicies(vt.begin(), vt.end(), vertices_.size(0));
        keys.insert(tmp.begin(), tmp.end());
      }
    } else {
      std::vector<ind_t> vn = nodes_.vertices(node);
      keys = permutation_table_keys_from_indicies(vn.begin(), vn.end(), vertices_.size(0));
    }
  }
  return keys;
}


} // end namespace brille
#endif
