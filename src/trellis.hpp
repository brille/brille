/* Copyright 2019 Greg Tucker
//
// This file is part of brille.
//
// brille is free software: you can redistribute it and/or modify it under the
// terms of the GNU Affero General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// brille is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with brille. If not, see <https://www.gnu.org/licenses/>.            */

#include <vector>
#include <array>
#include <algorithm>
#include <omp.h>
#include "arrayvector.hpp"
#include "latvec.hpp"
#include "polyhedron.hpp"
#include "utilities.hpp"
#include "debug.hpp"
#include "triangulation_simple.hpp"
#include "interpolation_data.hpp"

#ifndef _TRELLIS_H_
#define _TRELLIS_H_

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
//   if (i+1<size && approx_scalar(zero+step*(i+1),x)) return  1;
//   // if x is infinitesimally larger than zero+step*i
//   if (i  >0    && approx_scalar(zero+step*(i  ),x)) return -1;
//   return 0;
// }
template<class T>
static size_t find_bin(const std::vector<T>& bin_edges, const T x){
  auto found_at = std::find_if(bin_edges.begin(), bin_edges.end(), [x](double b){return b>x;});
  size_t d = std::distance(bin_edges.begin(), found_at);
  return d>0 ? d-1 : d;
}
template<class T>
static int on_boundary(const std::vector<T>& bin_edges, const T x, const size_t i){
  // if (i==0) then above d was *either* 0 or 1, otherwise d = i + 1;
  // if (i==0) we can't go lower in either case, so no problem.
  if (i+2<bin_edges.size() && approx_scalar(bin_edges[i+1],x)) return  1;
  if (i  >0                && approx_scalar(bin_edges[i  ],x)) return -1;
  return 0;
}

enum class NodeType {null, cube, poly};
// the number of nodes we might hold determines what type we need to store
// their indices:
//    bytes   maximum number    type
//      1                255    unsigned short
//      2             65,535    unsigned int
//      4      4,294,967,295    unsigned long
//      8            18×10¹⁸    unsigned long long (aka size_t)
// 65k is not enough. 4B *should* always be sufficient -- each node would
// occupy a fractional volume of ~2×10⁻¹⁰ of the polyhedron, which is overkill.
typedef unsigned long index_t;

class NullNode{
public:
  NullNode() {}
  virtual ~NullNode() = default;
  //
  virtual NodeType type(void) const {return NodeType::null;}
  virtual index_t vertex_count() const {return 0u;}
  virtual std::vector<index_t> vertices(void) const {return std::vector<index_t>();}
  virtual std::vector<std::array<index_t,4>> vertices_per_tetrahedron(void) const {return std::vector<std::array<index_t,4>>();}
  virtual bool indices_weights(const ArrayVector<double>&, const ArrayVector<double>&, std::vector<index_t>&, std::vector<double>&) const {return false;};
};
class CubeNode: public NullNode {
  std::array<index_t, 8> vertex_indices;
public:
  CubeNode(): vertex_indices({0,0,0,0,0,0,0,0}) {}
  CubeNode(const std::array<index_t,8>& vi): vertex_indices(vi) {}
  CubeNode(const std::vector<index_t>& vi): vertex_indices({0,0,0,0,0,0,0,0}) {
    if (vi.size() != 8) throw std::logic_error("CubeNode objects take 8 indices.");
    for (index_t i=0; i<8u; ++i) vertex_indices[i] = vi[i];
  }
  index_t vertex_count() const { return 8u;}
  std::vector<index_t> vertices(void) const {
    std::vector<index_t> out;
    for (auto v: vertex_indices) out.push_back(v);
    return out;
  }
  bool indices_weights(
    const ArrayVector<double>& vertices,
    const ArrayVector<double>& x,
    std::vector<index_t>& indices,
    std::vector<double>& weights
  ) const {
    // The CubeNode object contains the indices into `vertices` necessary to find
    // the 8 corners of the cube. Those indices should be ordered
    // (000) (100) (110) (010) (101) (001) (011) (111)
    // so that vertex_indices[i] and vertex_indices[7-i] are connected by a body diagonal
    ArrayVector<double> node_verts = vertices.extract(vertex_indices);
    double node_volume = abs(node_verts.extract(0)-node_verts.extract(7)).prod(1).getvalue(0,0);
    ArrayVector<double> w = abs(x - node_verts).prod(1)/node_volume; // the normalised volume of each sub-parallelpiped
    // If any normalised weights are greater than 1+eps() the point isn't in this node
    if (w.any_approx(Comp::gt,1.)) return false;
    ArrayVector<bool> needed = w.is_approx(Comp::gt, 0.);
    indices.clear();
    weights.clear();
    for (int i=0; i<8; ++i) if (needed.getvalue(i)) {
      // the weight corresponds to the vertex opposite the one used to find the partial volume
      indices.push_back(vertex_indices[7-i]);
      weights.push_back(w.getvalue(i));
    }
    return true;
  }
};
class PolyNode: public NullNode {
  std::vector<std::array<index_t,4>> vi_t;  //!< vertex indices per tetrahedra
  std::vector<std::array<double,4>> ci_t;   //!< circumsphere information per tetrahedra
  std::vector<double> vol_t;                //!< volume per tetrahedra
public:
  PolyNode() {};
  // actually constructing the tetrahedra from, e.g., a Polyhedron object will
  // need to be done elsewhere
  PolyNode(
    const std::vector<std::array<index_t,4>>& vit,
    const std::vector<std::array<double,4>>& cit,
    const std::vector<double>& volt
  ): vi_t(vit), ci_t(cit), vol_t(volt) {}
  // count-up the number fo unique vertices in the tetrahedra-triangulated polyhedron
  index_t tetrahedra_count() const {return vi_t.size();}
  index_t vertex_count() const { return this->vertices().size();}
  std::vector<index_t> vertices(void) const {
    std::vector<index_t> out;
    for (auto tet: vi_t) for (auto idx: tet)
    if (std::find(out.begin(), out.end(), idx)==out.end()) out.push_back(idx);
    return out;
  }
  std::vector<std::array<index_t,4>> vertices_per_tetrahedron(void) const {return vi_t;}
  bool indices_weights(
    const ArrayVector<double>& vertices,
    const ArrayVector<double>& x,
    std::vector<index_t>& indices,
    std::vector<double>& weights
  ) const {
    indices.clear();
    weights.clear();
    std::array<double,4> w{0,0,0,0};
    for (size_t i=0; i<vi_t.size(); ++i)
    if (this->tetrahedra_contains(i, vertices, x, w)){
      for (int j=0; j<4; ++j) if (!approx_scalar(w[j],0.)){
        indices.push_back(vi_t[i][j]);
        weights.push_back(w[j]);
      }
      return true;
    }
    return false;
  }
private:
  bool tetrahedra_contains(
    const index_t t,
    const ArrayVector<double>& v,
    const ArrayVector<double>& x,
    std::array<double,4>& w
  ) const {
    if (!this->tetrahedra_might_contain(t,x)) return false;
    double vol6 = vol_t[t]*6.0;
    w[0] = orient3d( x.data(),            v.data(vi_t[t][1u]), v.data(vi_t[t][2u]), v.data(vi_t[t][3u]) )/vol6;
    w[1] = orient3d( v.data(vi_t[t][0u]), x.data(),            v.data(vi_t[t][2u]), v.data(vi_t[t][3u]) )/vol6;
    w[2] = orient3d( v.data(vi_t[t][0u]), v.data(vi_t[t][1u]), x.data(),            v.data(vi_t[t][3u]) )/vol6;
    w[3] = orient3d( v.data(vi_t[t][0u]), v.data(vi_t[t][1u]), v.data(vi_t[t][2u]), x.data()            )/vol6;
    if (std::any_of(w.begin(), w.end(), [](double z){return z < 0. && !approx_scalar(z, 0.);}))
      return false;
    return true;
  }
  bool tetrahedra_might_contain(
    const index_t t,
    const ArrayVector<double>& x
  ) const {
    // find the vector from the circumsphere centre to x:
    double v[3];
    for (int i=0; i<3; ++i) v[i] = ci_t[t][i] - x.getvalue(0,i);
    // compute the squared length of v, and the circumsphere radius squared
    double d2{0}, r2 = ci_t[t][3]*ci_t[t][3];
    for (int i=0; i<3; ++i) d2 += v[i]*v[i];
    // if the squared distance is no greater than the squared radius, x might be inside the tetrahedra
    return d2 < r2 || approx_scalar(d2, r2);
  }
};

class NodeContainer{
  std::vector<std::pair<NodeType,index_t>> nodes_;
  std::vector<CubeNode> cube_nodes_;
  std::vector<PolyNode> poly_nodes_;
public:
  size_t size(void) const {return nodes_.size();}
  void push_back(const CubeNode& n){
    nodes_.emplace_back(NodeType::cube, cube_nodes_.size());
    cube_nodes_.push_back(n);
  }
  void push_back(const PolyNode& n){
    if (n.vertex_count() < 1)
      throw std::runtime_error("empty polynodes are not allowed!");
    nodes_.emplace_back(NodeType::poly, poly_nodes_.size());
    poly_nodes_.push_back(n);
  }
  void push_back(const NullNode&){
    nodes_.emplace_back(NodeType::null, (std::numeric_limits<index_t>::max)());
  }
  NodeType type(const index_t i) const {
    return nodes_[i].first;
  }
  bool is_cube(const index_t i) const {return NodeType::cube == nodes_[i].first;}
  bool is_poly(const index_t i) const {return NodeType::poly == nodes_[i].first;}
  bool is_null(const index_t i) const {return NodeType::null == nodes_[i].first;}
  const CubeNode& cube_at(const index_t i) const {
    return cube_nodes_[nodes_[i].second];
  }
  const PolyNode& poly_at(const index_t i) const {
    return poly_nodes_[nodes_[i].second];
  }
  index_t vertex_count(const index_t i) const {
    switch (nodes_[i].first){
      case NodeType::cube:
      return cube_nodes_[nodes_[i].second].vertex_count();
      case NodeType::poly:
      return poly_nodes_[nodes_[i].second].vertex_count();
      default:
      return 0;
    }
  }
  std::vector<index_t> vertices(const index_t i) const{
    switch (nodes_[i].first){
      case NodeType::cube:
      return cube_nodes_[nodes_[i].second].vertices();
      case NodeType::poly:
      return poly_nodes_[nodes_[i].second].vertices();
      default:
      return std::vector<index_t>();
    }
  }
  std::vector<std::array<index_t,4>> vertices_per_tetrahedron(const index_t i) const{
    if (nodes_[i].first == NodeType::poly)
      return poly_nodes_[nodes_[i].second].vertices_per_tetrahedron();
    return std::vector<std::array<index_t,4>>();
  }
  bool indices_weights(const index_t i, const ArrayVector<double>& v, const ArrayVector<double>& x, std::vector<index_t>& indices, std::vector<double>& weights) const{
    switch (nodes_[i].first){
      case NodeType::cube:
      return cube_nodes_[nodes_[i].second].indices_weights(v,x,indices,weights);
      case NodeType::poly:
      return poly_nodes_[nodes_[i].second].indices_weights(v,x,indices,weights);
      case NodeType::null:
        throw std::logic_error("attempting to access null node!");
      default:
      return false;
    }
  }
};

template<typename T> class PolyhedronTrellis{
  Polyhedron polyhedron_;                        //!< the Polyhedron bounding the Trellis
  InterpolationData<T> data_;                    //!< [optional] data stored at each Trellis vertex
  ArrayVector<double> vertices_;                 //!< The Trellis intersections inside the bounding Polyhedron
  NodeContainer nodes_;
  std::array<std::vector<double>,3> boundaries_; //!< The coordinates of the Trellis intersections, which bound the Trellis nodes
public:
  explicit PolyhedronTrellis(const Polyhedron& polyhedron, const double max_volume);
  // explicit PolyhedronTrellis(const Polyhedron& polyhedron, const double max_volume){
  //   this->construct(polyhedron, max_volume);
  // }
  // PolyhedronTrellis(const Polyhedron& polyhedron, const size_t number_density){
  //   double max_volume = polyhedron.get_volume()/static_cast<double>(number_density);
  //   this->construct(polyhedron, max_volume);
  // };
  PolyhedronTrellis(): vertices_({3,0}) {}
  index_t expected_vertex_count() const {
    index_t count = 1u;
    for (index_t i=0; i<3u; ++i) count *= boundaries_[i].size();
    return count;
  }
  index_t vertex_count() const { return vertices_.size(); }
  const ArrayVector<double>& vertices(void) const { return vertices_; }
  const ArrayVector<double>& vertices(const ArrayVector<double>& v){
    if (v.numel()==3) vertices_ = v;
    return vertices_;
  }
  ArrayVector<double> cube_vertices(void) const {
    std::vector<bool> keep(vertices_.size(), false);
    for (index_t i=0; i<nodes_.size(); ++i)
    if (nodes_.is_cube(i))
    for (auto idx: nodes_.vertices(i)) keep[idx] = true;
    return vertices_.extract(keep);
  }
  ArrayVector<double> poly_vertices(void) const {
    std::vector<bool> keep(vertices_.size(), false);
    for (index_t i=0; i<nodes_.size(); ++i)
    if (nodes_.is_poly(i))
    for (auto idx: nodes_.vertices(i)) keep[idx] = true;
    return vertices_.extract(keep);
  }
  std::vector<std::array<index_t,4>> vertices_per_tetrahedron(void) const {
    std::vector<std::array<index_t,4>> out;
    for (index_t i=0; i<nodes_.size(); ++i)
    if (nodes_.is_poly(i))
    for (auto tet: nodes_.vertices_per_tetrahedron(i)) out.push_back(tet);
    return out;
  }
  bool indices_weights(const ArrayVector<double>& x, std::vector<index_t>& indices, std::vector<double>& weights) const {
    if (x.size()!=1u || x.numel()!=3u)
      throw std::runtime_error("The indices and weights can only be found for one point at a time.");
    return nodes_.indices_weights(this->node_index(x), vertices_, x, indices, weights);
  }
  template<class R> unsigned check_before_interpolating(const ArrayVector<R>& x) const{
    unsigned int mask = 0u;
    if (this->data_.size()==0)
      throw std::runtime_error("The trellis must be filled before interpolating!");
    if (x.numel()!=3u)
      throw std::runtime_error("PolyhedronTrellis requires x values which are three-vectors.");
    return mask;
  }
  ArrayVector<T> interpolate_at(const ArrayVector<double>& x) const {
    verbose_update("Single thread interpolation at ",x.size()," points");
    this->check_before_interpolating(x);
    ArrayVector<T> out(data_.numel(), x.size());
    std::vector<index_t> indices;
    std::vector<double> weights;
    for (size_t i=0; i<x.size(); ++i){
      verbose_update("Locating ",x.to_string(i));
      if (!this->indices_weights(x.extract(i), indices, weights))
        throw std::runtime_error("Point not found in PolyhedronTrellis");
      verbose_update("Interpolate between vertices ", indices," with weights ",weights);
      data_.interpolate_at(indices, weights, out, i);
    }
    return out;
  }
  ArrayVector<T> interpolate_at(const ArrayVector<double>& x, const int threads) const {
    this->check_before_interpolating(x);
    omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
    verbose_update("Parallel interpolation at ",x.size()," points with ",threads," threads");
    // shared between threads
    ArrayVector<T> out(data_.numel(), x.size(), T(0)); // initialise to zero so that we can use a reduction
    // private to each thread
    std::vector<index_t> indices;
    std::vector<double> weights;
    // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
    long long xsize = unsigned_to_signed<long long, size_t>(x.size());
    size_t n_unfound{0};
  #pragma omp parallel for default(none) shared(x,out,xsize) private(indices, weights) reduction(+:n_unfound) schedule(dynamic)
    for (long long si=0; si<xsize; ++si){
      size_t i = signed_to_unsigned<size_t, long long>(si);
      if (this->indices_weights(x.extract(i), indices, weights)){
        data_.interpolate_at(indices, weights, out, i);
      } else {
        ++n_unfound;
      }
    }
    std::runtime_error("interpolate at failed to find "+std::to_string(n_unfound)+" point"+(n_unfound>1?"s.":"."));
    return out;
  }
  index_t node_count() {
    index_t count = 1u;
    for (index_t i=0; i<3u; ++i) count *= boundaries_[i].size()-1;
    return count;
  }
  std::array<index_t,3> size() const {
    std::array<index_t,3> s;
    for (index_t i=0; i<3u; ++i) s[i] = boundaries_[i].size()-1;
    return s;
  }
  std::array<index_t,3> span() const {
    std::array<index_t,3> s{1,0,0}, sz=this->size();
    for (index_t i=1; i<3; ++i) s[i] = sz[i-1]*s[i-1];
    return s;
  }
  // const std::array<std::vector<double>,3>& boundaries(void) const {return boundaries_;}
  //
  // Find the appropriate node for an arbitrary point:
  std::array<index_t,3> node_subscript(const ArrayVector<double>& p) const {
    std::array<index_t,3> sub{0,0,0};
    for (index_t dim=0; dim<3u; ++dim)
      sub[dim] = find_bin(boundaries_[dim], p.getvalue(0, dim));
    bool bad = nodes_.is_null(this->sub2idx(sub));
    if (bad){
      std::array<int,3> close{0,0,0};
      // determine if we are close to a boundary along any of the three binning
      // directions. if we are, on_boundary returns the direction in which we
      // can safely take a step without leaving the binned region
      for (int i=0; i<3; ++i) close[i] = on_boundary(boundaries_[i], p.getvalue(0,i), sub[i]);
      int num_close = std::count_if(close.begin(), close.end(), [](int a){return a!=0;});
      // check one
      std::array<index_t,3> newsub{sub};
      if (num_close > 0) for (int i=0; i<3 && bad; ++i) if (close[i]) {
        newsub = sub;
        newsub[i] += close[i];
        bad = nodes_.is_null(this->sub2idx(newsub));
      }
      // check two
      if (bad && num_close>1)
      for (int i=0; i<3 && bad; ++i) if (close[i])
      for (int j=0; j<3 && bad; ++j) if (close[j]) {
        newsub = sub;
        newsub[i] += close[i];
        newsub[j] += close[j];
        bad = nodes_.is_null(this->sub2idx(newsub));
      }
      // check all three
      if (bad && num_close>2){
        newsub = sub;
        for (int i=0; i<3; ++i) newsub[i] += close[i];
        bad = nodes_.is_null(this->sub2idx(newsub));
      }
      if (!bad) sub = newsub;
    }
    return sub;
  }
  // find the node linear index for a point
  template <class R> index_t node_index(const R& p) const { return this->sub2idx(this->node_subscript(p)); }

  const CubeNode& cube_node(const index_t idx) const {
    if (idx >= nodes_.size() || !nodes_.is_cube(idx))
      throw std::runtime_error("Out-of-bounds or non-cube node");
    return nodes_.cube_at(idx);
  }
  const PolyNode& poly_node(const index_t idx) const {
    if (idx >= nodes_.size() || !nodes_.is_poly(idx))
      throw std::runtime_error("Out-of-bounds or non-polyhedron node");
    return nodes_.poly_at(idx);
  }

  std::string to_string(void) const {
    std::string str = "(";
    for (auto i: this->size()) str += " " + std::to_string(i);
    str += " )";
    return str;
  }
  // Get a constant reference to the stored data
  const InterpolationData<T>& data(void) const {return data_;}
  // Replace the data stored in the object
  template<typename... A> void replace_data(A... args) { data_.replace_data(args...); }
  // Calculate the Debye-Waller factor for the provided Q points and ion masses
  template<template<class> class A>
  ArrayVector<double> debye_waller(const A<double>& Q, const std::vector<double>& M, const double t_K) const{
    return data_.debye_waller(Q,M,t_K);
  }
private:
  index_t sub2idx(const std::array<index_t,3>& sub) const {
    return this->sub2idx(sub, this->span());
  }
  std::array<index_t,3> idx2sub(const index_t idx) const {
    return this->idx2sub(idx, this->span());
  }
  index_t sub2idx(const std::array<index_t,3>& sub, const std::array<index_t,3>& sp) const {
    index_t idx=0;
    for (index_t dim=0; dim<3u; ++dim) idx += sp[dim]*sub[dim];
    // info_update("span ",sp," gives subscript ",sub," as linear ",idx);
    return idx;
  }
  std::array<index_t,3> idx2sub(const index_t idx, const std::array<index_t,3>& sp) const {
    std::array<index_t,3> sub{0,0,0};
    index_t rem{idx};
    for (index_t dim=3u; dim--;){
      sub[dim] = rem/sp[dim];
      rem -= sub[dim]*sp[dim];
    }
    return sub;
  }
  // template<typename R> void add_node(const R& node) {nodes_.push_back(node);}
};

#include "trellis.tpp"
#endif
