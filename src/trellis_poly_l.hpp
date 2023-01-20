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
#include <queue>
#include <condition_variable>
#include <atomic>
#include <functional>
#include <utility>
#include <filesystem>
#include "interpolatordual.hpp"
#include "hdf_interface.hpp"
#include "trellis_node.hpp"
#include "polyhedron_flex.hpp"
#include "triangulation_poly.hpp"
#include "approx_float.hpp"
#include "approx_config.hpp"

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
  std::tuple<std::map<size_t, poly_t>, std::vector<std::vector<ind_t>>, std::vector<ind_t>, VertexType<VertexComponents>, ind_t>
  part_one(const poly_t&, const VertexType<VertexComponents>&,
           std::vector<NodeType>&, bool, VertexComponents, int);
  void
  part_two(const std::map<size_t, poly_t>&, const std::vector<NodeType>& ,
           const std::vector<std::vector<ind_t>>&, ind_t, ind_t,
           VertexComponents, int);
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
    omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
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
    auto xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
    size_t n_unfound{0};
  #pragma omp parallel for default(none) shared(x,vals2,vecs2,xsize) reduction(+:n_unfound) schedule(dynamic)
    for (long long si=0; si<xsize; ++si){
      auto i = brille::utils::s2u<ind_t, long long>(si);
      auto indwghts = this->indices_weights(x.view(i));
      if (indwghts.size()>0) {
        data_.interpolate_at(indwghts, vals2, vecs2, i);
      } else {
        ++n_unfound;
      }
    }
    if (n_unfound){
      throw std::runtime_error("interpolate at failed to find "+std::to_string(n_unfound)+" point"+(n_unfound>1?"s.":"."));
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
      else
      info_update("The node subscript ",sub,
                  " for the point (hkl) ", p.to_string(0u), "(xyz) ", pos,
                  " is either invalid or points to a null node!");
    }
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

#ifdef USE_HIGHFIVE
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
  bool to_hdf(const std::string& filename, const std::string& entry, const unsigned perm=HighFive::File::OpenOrCreate) const {
    HighFive::File file(filename, perm);
    return this->to_hdf(file, entry);
  }
  static class_t from_hdf(const std::string& filename, const std::string& entry){
    HighFive::File file(filename, HighFive::File::ReadOnly);
    return class_t::from_hdf(file, entry);
  }
#endif // USE_HIGHFIVE
};

#include "trellis_poly_l.tpp"
} // end namespace brille
#endif
