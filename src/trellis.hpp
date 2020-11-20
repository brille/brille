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
#include <vector>
#include <array>
#include <queue>
#include <tuple>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <algorithm>
#include <functional>
#include <omp.h>
#include "array.hpp"
#include "array2.hpp"
#include "array_latvec.hpp" // defines bArray
#include "polyhedron.hpp"
#include "utilities.hpp"
#include "debug.hpp"
#include "triangulation_simple.hpp"
#include "interpolatordual.hpp"
#include "permutation.hpp"
#include "approx.hpp"
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

enum class NodeType {null, cube, poly};


class NullNode{
public:
  using ind_t = brille::ind_t;
  NullNode() {}
  virtual ~NullNode() = default;
  //
  virtual NodeType type(void) const {return NodeType::null;}
  virtual ind_t vertex_count() const {return 0u;}
  virtual std::vector<ind_t> vertices(void) const {return std::vector<ind_t>();}
  virtual std::vector<std::array<ind_t,4>> vertices_per_tetrahedron(void) const {return std::vector<std::array<ind_t,4>>();}
  // bool indices_weights(const bArray<double>&, const bArray<double>&, std::vector<ind_t>&, std::vector<double>&) const {return false;}
  bool indices_weights(const bArray<double>&, const bArray<double>&, std::vector<std::pair<ind_t,double>>&) const {return false;}
  double volume(const bArray<double>&) const {return 0.;}
};
class CubeNode: public NullNode {
  std::array<ind_t, 8> vertex_indices;
public:
  CubeNode(): vertex_indices({{0,0,0,0,0,0,0,0}}) {}
  CubeNode(const std::array<ind_t,8>& vi): vertex_indices(vi) {}
  CubeNode(const std::vector<ind_t>& vi): vertex_indices({{0,0,0,0,0,0,0,0}}) {
    if (vi.size() != 8) throw std::logic_error("CubeNode objects take 8 indices.");
    for (ind_t i=0; i<8u; ++i) vertex_indices[i] = vi[i];
  }
  ind_t vertex_count() const { return 8u;}
  std::vector<ind_t> vertices(void) const {
    std::vector<ind_t> out;
    for (auto v: vertex_indices) out.push_back(v);
    return out;
  }
  // bool indices_weights(
  //   const bArray<double>& vertices,
  //   const bArray<double>& x,
  //   std::vector<ind_t>& indices,
  //   std::vector<double>& weights
  // ) const {
  //   // The CubeNode object contains the indices into `vertices` necessary to find
  //   // the 8 corners of the cube. Those indices should be ordered
  //   // (000) (100) (110) (010) (101) (001) (011) (111)
  //   // so that vertex_indices[i] and vertex_indices[7-i] are connected by a body diagonal
  //   auto node_verts = vertices.extract(vertex_indices);
  //   double node_volume = abs(node_verts.view(0)-node_verts.view(7)).prod(1)[0];
  //   auto w = abs(x - node_verts).prod(1)/node_volume; // the normalised volume of each sub-parallelpiped
  //   // If any normalised weights are greater than 1+eps() the point isn't in this node
  //   if (w.any(brille::cmp::gt, 1.)) return false;
  //   auto needed = w.is(brille::cmp::gt, 0.);
  //   indices.clear();
  //   weights.clear();
  //   for (int i=0; i<8; ++i) if (needed[i]) {
  //     // the weight corresponds to the vertex opposite the one used to find the partial volume
  //     indices.push_back(vertex_indices[7-i]);
  //     weights.push_back(w[i]);
  //   }
  //   return true;
  // }
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
  double volume(const bArray<double>& vertices) const {
    // The CubeNode object contains the indices into `vertices` necessary to find
    // the 8 corners of the cube. Those indices should be ordered
    // (000) (100) (110) (010) (101) (001) (011) (111)
    // so that vertex_indices[i] and vertex_indices[7-i] are connected by a body diagonal
    return abs(vertices.view(vertex_indices[0])-vertices.view(vertex_indices[7])).prod(1)[0];
  }
};
class PolyNode: public NullNode {
  std::vector<std::array<ind_t,4>> vi_t;  //!< vertex indices per tetrahedra
  std::vector<std::array<double,4>> ci_t;   //!< circumsphere information per tetrahedra
  std::vector<double> vol_t;                //!< volume per tetrahedra
public:
  PolyNode() {};
  // actually constructing the tetrahedra from, e.g., a Polyhedron object will
  // need to be done elsewhere
  PolyNode(
    const std::vector<std::array<ind_t,4>>& vit,
    const std::vector<std::array<double,4>>& cit,
    const std::vector<double>& volt
  ): vi_t(vit), ci_t(cit), vol_t(volt) {}
  // count-up the number fo unique vertices in the tetrahedra-triangulated polyhedron
  ind_t tetrahedra_count() const {return static_cast<ind_t>(vi_t.size());}
  ind_t vertex_count() const { return static_cast<ind_t>(this->vertices().size());}
  std::vector<ind_t> vertices(void) const {
    std::vector<ind_t> out;
    for (auto tet: vi_t) for (auto idx: tet)
    if (std::find(out.begin(), out.end(), idx)==out.end()) out.push_back(idx);
    return out;
  }
  std::vector<std::array<ind_t,4>> vertices_per_tetrahedron(void) const {return vi_t;}
  // bool indices_weights(
  //   const bArray<double>& vertices,
  //   const bArray<double>& x,
  //   std::vector<ind_t>& indices,
  //   std::vector<double>& weights
  // ) const {
  //   indices.clear();
  //   weights.clear();
  //   std::array<double,4> w{{0,0,0,0}};
  //   for (ind_t i=0; i<vi_t.size(); ++i)
  //   if (this->tetrahedra_contains(i, vertices, x, w)){
  //     for (int j=0; j<4; ++j) if (!brille::approx::scalar(w[j],0.)){
  //       indices.push_back(vi_t[i][j]);
  //       weights.push_back(w[j]);
  //     }
  //     return true;
  //   }
  //   return false;
  // }
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

class NodeContainer{
public:
  using ind_t = brille::ind_t;
private:
  std::vector<std::pair<NodeType,ind_t>> nodes_;
  std::vector<CubeNode> cube_nodes_;
  std::vector<PolyNode> poly_nodes_;
public:
  size_t size(void) const {return nodes_.size();}
  size_t cube_count() const {
    return std::count_if(nodes_.begin(),nodes_.end(),[](std::pair<NodeType,ind_t> n){return NodeType::cube == n.first;});
  }
  size_t poly_count() const {
    return std::count_if(nodes_.begin(),nodes_.end(),[](std::pair<NodeType,ind_t> n){return NodeType::poly == n.first;});
  }
  size_t null_count() const {
    return std::count_if(nodes_.begin(),nodes_.end(),[](std::pair<NodeType,ind_t> n){return NodeType::null == n.first;});
  }
  void push_back(const CubeNode& n){
    nodes_.emplace_back(NodeType::cube, static_cast<ind_t>(cube_nodes_.size()));
    cube_nodes_.push_back(n);
  }
  void push_back(const PolyNode& n){
    if (n.vertex_count() < 1)
      throw std::runtime_error("empty polynodes are not allowed!");
    nodes_.emplace_back(NodeType::poly, static_cast<ind_t>(poly_nodes_.size()));
    poly_nodes_.push_back(n);
  }
  void push_back(const NullNode&){
    nodes_.emplace_back(NodeType::null, (std::numeric_limits<ind_t>::max)());
  }
  NodeType type(const ind_t i) const {
    return nodes_[i].first;
  }
  bool is_cube(const ind_t i) const {return NodeType::cube == nodes_[i].first;}
  bool is_poly(const ind_t i) const {return NodeType::poly == nodes_[i].first;}
  bool is_null(const ind_t i) const {return NodeType::null == nodes_[i].first;}
  const CubeNode& cube_at(const ind_t i) const {
    return cube_nodes_[nodes_[i].second];
  }
  const PolyNode& poly_at(const ind_t i) const {
    return poly_nodes_[nodes_[i].second];
  }
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
  std::vector<std::array<ind_t,4>> vertices_per_tetrahedron(const ind_t i) const{
    if (nodes_[i].first == NodeType::poly)
      return poly_nodes_[nodes_[i].second].vertices_per_tetrahedron();
    return std::vector<std::array<ind_t,4>>();
  }
  // bool indices_weights(const ind_t i, const bArray<double>& v, const bArray<double>& x, std::vector<ind_t>& indices, std::vector<double>& weights) const{
  //   switch (nodes_[i].first){
  //     case NodeType::cube:
  //     return cube_nodes_[nodes_[i].second].indices_weights(v,x,indices,weights);
  //     case NodeType::poly:
  //     return poly_nodes_[nodes_[i].second].indices_weights(v,x,indices,weights);
  //     case NodeType::null:
  //       throw std::logic_error("attempting to access null node!");
  //     default:
  //     return false;
  //   }
  // }
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

template<typename T, typename R>
class PolyhedronTrellis{
public:
  using ind_t = brille::ind_t;
  using data_t = DualInterpolator<T,R>;
  using vert_t = bArray<double>;
private:
  Polyhedron polyhedron_;                        //!< the Polyhedron bounding the Trellis
  data_t data_;                  //!< [optional] data stored at each Trellis vertex
  vert_t vertices_;               //!< The Trellis intersections inside the bounding Polyhedron
  NodeContainer nodes_;
  std::array<std::vector<double>,3> boundaries_; //!< The coordinates of the Trellis intersections, which bound the Trellis nodes
public:
  explicit PolyhedronTrellis(const Polyhedron& polyhedron, const double max_volume, const bool always_triangulate=false);
  // explicit PolyhedronTrellis(const Polyhedron& polyhedron, const double max_volume){
  //   this->construct(polyhedron, max_volume);
  // }
  // PolyhedronTrellis(const Polyhedron& polyhedron, const size_t number_density){
  //   double max_volume = polyhedron.get_volume()/static_cast<double>(number_density);
  //   this->construct(polyhedron, max_volume);
  // };
  explicit PolyhedronTrellis(): vertices_(0,3) {}
  ind_t expected_vertex_count() const {
    ind_t count = 1u;
    for (ind_t i=0; i<3u; ++i) count *= boundaries_[i].size();
    return count;
  }
  ind_t vertex_count() const { return static_cast<ind_t>(vertices_.size(0)); }
  const vert_t& vertices(void) const { return vertices_; }
  const vert_t& vertices(const bArray<double>& v){
    if (v.ndim()==2 && v.size(1)==3) vertices_ = v;
    return vertices_;
  }
  vert_t cube_vertices(void) const {
    std::vector<bool> keep(vertices_.size(0), false);
    for (ind_t i=0; i<nodes_.size(); ++i)
    if (nodes_.is_cube(i))
    for (auto idx: nodes_.vertices(i)) keep[idx] = true;
    return vertices_.extract(keep);
  }
  vert_t poly_vertices(void) const {
    std::vector<bool> keep(vertices_.size(0), false);
    for (ind_t i=0; i<nodes_.size(); ++i)
    if (nodes_.is_poly(i))
    for (auto idx: nodes_.vertices(i)) keep[idx] = true;
    return vertices_.extract(keep);
  }
  std::vector<std::array<ind_t,4>> vertices_per_tetrahedron(void) const {
    std::vector<std::array<ind_t,4>> out;
    for (ind_t i=0; i<nodes_.size(); ++i)
    if (nodes_.is_poly(i))
    for (auto tet: nodes_.vertices_per_tetrahedron(i)) out.push_back(tet);
    return out;
  }
  // bool indices_weights(const bArray<double>& x, std::vector<ind_t>& indices, std::vector<double>& weights) const {
  //   if (x.ndim()!=2 && x.size(0)!=1u && x.size(1)!=3u)
  //     throw std::runtime_error("The indices and weights can only be found for one point at a time.");
  //   return nodes_.indices_weights(this->node_index(x), vertices_, x, indices, weights);
  // }
  std::vector<std::pair<ind_t,double>>
  indices_weights(const bArray<double>& x) const {
    std::vector<std::pair<ind_t,double>> iw;
    if (x.ndim()!=2 && x.size(0)!=1u && x.size(1)!=3u)
      throw std::runtime_error("The indices and weights can only be found for one point at a time.");
    nodes_.indices_weights(this->node_index(x), vertices_, x, iw);
    return iw;
  }
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
    for (size_t i=0; i<x.size(0); ++i){
      verbose_update("Locating ",x.to_string(i));
      auto indwghts = this->indices_weights(x.view(i));
      if (indwghts.size()<1)
        throw std::runtime_error("Point not found in PolyhedronTrellis");
      data_.interpolate_at(indwghts, vals2, vecs2, i);
    }
    return std::make_tuple(vals_out, vecs_out);
  }
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
    long long xsize = brille::utils::u2s<long long, size_t>(x.size(0));
    size_t n_unfound{0};
  #pragma omp parallel for default(none) shared(x,vals2,vecs2,xsize) reduction(+:n_unfound) schedule(dynamic)
    for (long long si=0; si<xsize; ++si){
      size_t i = brille::utils::s2u<size_t, long long>(si);
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
  ind_t node_count() {
    ind_t count = 1u;
    for (ind_t i=0; i<3u; ++i) count *= static_cast<ind_t>(boundaries_[i].size()-1);
    return count;
  }
  std::array<ind_t,3> size() const {
    std::array<ind_t,3> s;
    for (ind_t i=0; i<3u; ++i) s[i] = static_cast<ind_t>(boundaries_[i].size()-1);
    return s;
  }
  std::array<ind_t,3> span() const {
    std::array<ind_t,3> s{{1,0,0}}, sz=this->size();
    for (ind_t i=1; i<3; ++i) s[i] = sz[i-1]*s[i-1];
    return s;
  }
  // const std::array<std::vector<double>,3>& boundaries(void) const {return boundaries_;}
  //
  // Find the appropriate node for an arbitrary point:
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
  // find the node linear index for a point
  template <class S> ind_t node_index(const S& p) const { return this->sub2idx(this->node_subscript(p)); }

  // return a list of non-null neighbouring nodes
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

  const CubeNode& cube_node(const ind_t idx) const {
    if (idx >= nodes_.size() || !nodes_.is_cube(idx))
      throw std::runtime_error("Out-of-bounds or non-cube node");
    return nodes_.cube_at(idx);
  }
  const PolyNode& poly_node(const ind_t idx) const {
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
  //! Get a constant reference to the stored data
  const data_t& data(void) const {return data_;}
  //! Replace the data stored in the object
  template<typename... A> void replace_data(A... args) {data_.replace_data(args...);}
  template<typename... A> void replace_value_data(A... args) { data_.replace_value_data(args...); }
  template<typename... A> void replace_vector_data(A... args) { data_.replace_vector_data(args...); }
  template<typename... A> void set_value_cost_info(A... args) { data_.set_value_cost_info(args...); }
  template<typename... A> void set_vector_cost_info(A... args) {data_.set_vector_cost_info(args...);}
  //! Return the number of bytes used per Q point
  size_t bytes_per_point() const {return data_.bytes_per_point(); }
  //! Calculate the Debye-Waller factor for the provided Q points and ion masses
  template<template<class> class A>
  brille::Array<double> debye_waller(const A<double>& Qpts, const std::vector<double>& Masses, const double t_K) const{
    return data_.debye_waller(Qpts,Masses,t_K);
  }
  void sort(void){ data_.sort(); }
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
  std::array<size_t,3> trellis_intersections_span() const {
    size_t bs0{boundaries_[0].size()}, bs1{boundaries_[1].size()};
    return std::array<size_t,3>({1, bs0, bs0*bs1});
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
