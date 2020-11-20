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

#ifndef BRILLE_NEST_H_
#define BRILLE_NEST_H_
#include <set>
#include <vector>
#include <array>
#include <tuple>
#include <utility>
#include <algorithm>
#include <omp.h>
#include "array.hpp"
#include "array2.hpp"
#include "array_latvec.hpp"
#include "polyhedron.hpp"
#include "utilities.hpp"
#include "debug.hpp"
#include "triangulation_simple.hpp"
#include "interpolatordual.hpp"
#include "approx.hpp"
namespace brille {

template<typename T,size_t N>
static bool none_negative(const std::array<T,N>& x){
  return !std::any_of(x.begin(), x.end(), [](T z){return z<0 && !brille::approx::scalar(z,0.);});
}

class NestLeaf{
public:
  using ind_t = brille::ind_t;
private:
  std::array<ind_t,4> vi;
  std::array<double,4> centre_radius;
  double volume_;
public:
  explicit NestLeaf(): vi({{0,0,0,0}}), centre_radius({{0,0,0,0}}), volume_(0) {}
  // actually constructing the tetrahedra from, e.g., a Polyhedron object will
  // need to be done elsewhere
  NestLeaf(
    const std::array<ind_t,4>& vit,
    const std::array<double,4>& ci,
    const double vol
  ): vi(vit), centre_radius(ci), volume_(vol) {}
  //
  const std::array<ind_t,4>& vertices(void) const { return vi;}
  //
  double volume(void) const {return volume_;}
  //
  std::array<double,4> weights(const bArray<double>& v, const bArray<double>& x) const {
    std::array<double,4> w{{-1,-1,-1,-1}};
    if (this->might_contain(x)){
      double vol6 = volume_*6.0;
      w[0] = orient3d(x.ptr(0)    , v.ptr(vi[1]), v.ptr(vi[2]), v.ptr(vi[3])) / vol6;
      w[1] = orient3d(v.ptr(vi[0]), x.ptr(0)    , v.ptr(vi[2]), v.ptr(vi[3])) / vol6;
      w[2] = orient3d(v.ptr(vi[0]), v.ptr(vi[1]), x.ptr(0)    , v.ptr(vi[3])) / vol6;
      w[3] = orient3d(v.ptr(vi[0]), v.ptr(vi[1]), v.ptr(vi[2]), x.ptr(0)    ) / vol6;
    }
    return w;
  }
  bool contains(
    const bArray<double>& v,
    const bArray<double>& x,
    std::array<double,4>& w
  ) const {
    if (this->might_contain(x)){
      double vol6 = volume_*6.0;
      w[0] = orient3d(x.ptr(0)    , v.ptr(vi[1]), v.ptr(vi[2]), v.ptr(vi[3])) / vol6;
      w[1] = orient3d(v.ptr(vi[0]), x.ptr(0)    , v.ptr(vi[2]), v.ptr(vi[3])) / vol6;
      w[2] = orient3d(v.ptr(vi[0]), v.ptr(vi[1]), x.ptr(0)    , v.ptr(vi[3])) / vol6;
      w[3] = orient3d(v.ptr(vi[0]), v.ptr(vi[1]), v.ptr(vi[2]), x.ptr(0)    ) / vol6;
      return none_negative(w);
    }
    return false;
  }
  //
  std::string to_string(void) const {
    std::string msg = "[";
    for (auto i: vi) msg += " " + std::to_string(i);
    msg += " ]";
    return msg;
  }
private:
  bool might_contain(const bArray<double>& x) const {
    std::array<double,3> d;
    d[0] = x.val(0,0) - centre_radius[0];
    d[1] = x.val(0,1) - centre_radius[1];
    d[2] = x.val(0,2) - centre_radius[2];
    double d2{0}, r2 = centre_radius[3]*centre_radius[3];
    for (size_t i=0; i<3u; ++i) d2 += d[i]*d[i];
    return ( d2 < r2 || brille::approx::scalar(d2,r2) );
  }
};


class NestNode{
public:
  using ind_t = brille::ind_t;
private:
  bool is_root_;
  NestLeaf boundary_;
  std::vector<NestNode> branches_;
public:
  explicit NestNode(bool ir=false): is_root_(ir), boundary_() {}
  explicit NestNode(const NestLeaf& b): is_root_(false), boundary_(b) {}
  NestNode(
    const std::array<ind_t,4>& vit,
    const std::array<double,4>& ci,
    const double vol
  ): is_root_(false), boundary_(NestLeaf(vit,ci,vol)) {}
  bool is_root(void) const {return is_root_;}
  bool is_leaf(void) const {return !is_root_ && branches_.size()==0;}
  const NestLeaf& boundary(void) const {return boundary_;}
  const std::vector<NestNode>& branches(void) const {return branches_;}
  std::vector<NestNode>& branches(void) {return branches_;}
  double volume(void) const {return boundary_.volume();}
  template<typename... A> bool contains(A... args) {return boundary_.contains(args...);}
  template<typename... A> std::array<double,4> weights(A... args) {return boundary_.weights(args...);}
  std::vector<std::pair<ind_t,double>> indices_weights(
    const bArray<double>& v, const bArray<double>& x
  ) const
  {
    std::array<double,4> w;
    return __indices_weights(v,x,w);
  }
  std::vector<std::array<ind_t,4>> tetrahedra(void) const {
    std::vector<std::array<ind_t,4>> out;
    if (this->is_leaf()) out.push_back(boundary_.vertices());
    for (auto b: branches_) for (auto v: b.tetrahedra()) out.push_back(v);
    return out;
  }
  std::string to_string(const std::string& prefix, const bool not_last) const {
    std::string msg = prefix;
    msg += is_root_ ? "───┐" : not_last ? "├──" : "└──";
    if (!is_root_) msg += boundary_.to_string();
    msg += "\n";
    for (size_t i=0; i<branches_.size(); ++i)
      msg += branches_[i].to_string(prefix+(not_last?"|  ":"   "),i+1!=branches_.size());
    return msg;
  }
  std::set<size_t> collect_keys(const size_t nv) const {
    std::set<size_t> keys;
    if (this->is_leaf()){
      auto v = boundary_.vertices();
      keys = permutation_table_keys_from_indicies(v.begin(), v.end(), nv);
    } else {
      for (auto b: branches_) {
        std::set<size_t> t = b.collect_keys(nv);
        keys.insert(t.begin(), t.end());
      }
    }
    return keys;
  }
protected:
  std::vector<std::pair<ind_t,double>> __indices_weights(
    const bArray<double>& v, const bArray<double>& x, std::array<double,4>& w
  ) const
  {
    // This node is either the root (in which case it contains all tetrahedra)
    // or a tetrahedra below the root which definately contains the point x
    // In the second case w has been set for us by the calling function.
    if (this->is_leaf()){
      std::array<ind_t,4> vi = boundary_.vertices();
      std::vector<std::pair<ind_t,double>> iw;
      for (size_t i=0; i<4u; ++i) if (!brille::approx::scalar(w[i], 0.))
        iw.push_back(std::make_pair(vi[i], w[i]));
      return iw;
    }
    // This is not a leaf node. So continue down the tree
    for (auto b: branches_){
      w = b.weights(v,x);
      if (none_negative(w)) return b.__indices_weights(v,x,w);
    }
    std::vector<std::pair<ind_t,double>> empty;
    return empty;
  }
};
template<class T, class S>
class Nest{
public:
  using ind_t = brille::ind_t;
  using data_t = DualInterpolator<T,S>;
  using vert_t = bArray<double>;
private:
  NestNode root_;
  vert_t vertices_;
  data_t data_;
  // std::vector<size_t> map_; // vertices holds *all* vertices but data_ only holds information for terminal vertices!
public:
  std::string tree_string(void) const {
    std::string tree = root_.to_string("",false);
    return tree;
  }
  // Build using maximum leaf volume
  Nest(const Polyhedron& p, const double vol, const size_t nb=5u)
  : root_(true), vertices_(0u,3u)
  {
    this->construct(p, nb, vol);
    // this->make_all_to_terminal_map();
    size_t nvert = this->vertex_count();
    data_.initialize_permutation_table(nvert, root_.collect_keys(nvert));
  }
  // Build using desired leaf number density
  Nest(const Polyhedron& p, const size_t rho, const size_t nb=5u)
  : root_(true), vertices_(0u,3u)
  {
    this->construct(p, nb, p.get_volume()/static_cast<double>(rho));
    // this->make_all_to_terminal_map();
    size_t nvert = this->vertex_count();
    data_.initialize_permutation_table(nvert, root_.collect_keys(nvert));
  }
  std::vector<bool> vertex_is_leaf(void) const {
    std::vector<bool> vert_is_term(vertices_.size(0), false);
    for (auto tet: root_.tetrahedra()) for (auto idx: tet) vert_is_term[idx]=true;
    return vert_is_term;
  }
  const vert_t& all_vertices(void) const {return vertices_;}
  vert_t vertices(void) const{ return vertices_; }
  brille::ind_t vertex_count() const { return vertices_.size(0); }
  std::vector<std::array<ind_t,4>> tetrahedra(void) const {
    std::vector<std::array<ind_t,4>> all_tet = root_.tetrahedra();
    // we need to adjust indexing to be into vertices instead of all_vertices
    /* (do this later) */
    return all_tet;
  }
  std::vector<std::pair<ind_t,double>>
  indices_weights(const bArray<double> &x) const {
    if (x.ndim()!=2 || x.size(0) != 1u || x.size(1) != 3u)
      throw std::runtime_error("The indices and weights can only be found for one point at a time.");
    // return root_.indices_weights(vertices_, map_, x);
    return root_.indices_weights(vertices_, x);
  }
  template<class R>
  unsigned check_before_interpolating(const bArray<R>& x) const{
    unsigned int mask = 0u;
    if (data_.size()==0)
      throw std::runtime_error("The interpolation data must be filled before interpolating.");
    if (x.ndim()!=2 || x.size(1)!=3u)
      throw std::runtime_error("Only (n,3) two-dimensional Q vectors supported in interpolating.");
    if (x.stride().back()!=1)
      throw std::runtime_error("Contiguous vectors required for interpolation.");
    return mask;
  }
  std::tuple<brille::Array<T>, brille::Array<S>>
  interpolate_at(const bArray<double>& x) const {
    this->check_before_interpolating(x);
    auto valsh = data_.values().shape();
    auto vecsh = data_.values().shape();
    valsh[0] = x.size(0);
    vecsh[0] = x.size(0);
    brille::Array<T> vals(valsh);
    brille::Array<S> vecs(vecsh);
    // vals and vecs are row-ordered contiguous by default, so we can create
    // mutable data-sharing Array2 objects for use with
    // Interpolator2::interpolate_at through the constructor:
    brille::Array2<T> vals2(vals);
    brille::Array2<S> vecs2(vecs);
    for (ind_t i=0; i<x.size(0); ++i){
      // auto iw = root_.indices_weights(vertices_, map_, x.extract(i));
      auto iw = root_.indices_weights(vertices_, x.extract(i));
      data_.interpolate_at(iw, vals2, vecs2, i);
    }
    return std::make_tuple(vals, vecs);
  }
  std::tuple<brille::Array<T>, brille::Array<S>>
  interpolate_at(const bArray<double>& x, const int threads) const {
    this->check_before_interpolating(x);
    omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
    // not used in parallel region
    auto valsh = data_.values().shape();
    auto vecsh = data_.vectors().shape();
    valsh[0] = x.size(0);
    vecsh[0] = x.size(0);
    // shared between threads
    brille::Array<T> vals(valsh);
    brille::Array<S> vecs(vecsh);
    // vals and vecs are row-ordered contiguous by default, so we can create
    // mutable data-sharing Array2 objects for use with
    // Interpolator2::interpolate_at through the constructor:
    brille::Array2<T> vals2(vals);
    brille::Array2<S> vecs2(vecs);
    // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
    ind_t unfound=0;
    long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
  #pragma omp parallel for default(none) shared(x, vals2, vecs2) reduction(+:unfound) firstprivate(xsize) schedule(dynamic)
    for (long long si=0; si<xsize; ++si){
      ind_t i = brille::utils::s2u<ind_t, long long>(si);
      // auto iw = root_.indices_weights(vertices_, map_, x.extract(i));
      auto iw = root_.indices_weights(vertices_, x.extract(i));
      if (iw.size()){
        data_.interpolate_at(iw, vals2, vecs2, i);
      } else {
        ++unfound;
      }
    }
    if (unfound > 0){
      std::string msg = std::to_string(unfound) + " points not found in Nest";
      throw std::runtime_error(msg);
    }
    return std::make_tuple(vals, vecs);
  }
  const data_t& data(void) const {return data_;}
  template<typename... A> void replace_data(A... args) {data_.replace_data(args...);}
  template<typename... A> void replace_value_data(A... args) { data_.replace_value_data(args...); }
  template<typename... A> void replace_vector_data(A... args) { data_.replace_vector_data(args...); }
  template<typename... A> void set_value_cost_info(A... args) { data_.set_value_cost_info(args...); }
  template<typename... A> void set_vector_cost_info(A... args) {data_.set_vector_cost_info(args...);}
  //! Return the number of bytes used per Q point
  size_t bytes_per_point() const {return data_.bytes_per_point(); }
  template<template<class> class A>
  brille::Array<double> debye_waller(const A<double>& Qpts, const std::vector<double>& Masses, const double t_K) const{
    return data_.debye_waller(Qpts,Masses,t_K);
  }
  void sort() {data_.sort();}
private:
  void construct(const Polyhedron&, const size_t, const double);
  // void make_all_to_terminal_map(void) {
  //   std::vector<bool> vit = this->vertex_is_leaf();
  //   size_t nTerminal = std::count(vit.begin(), vit.end(), true);
  //   map_.clear();
  //   size_t idx{0};
  //   for (bool t: vit) map_.push_back(t ? idx++ : nTerminal);
  //   if (idx != nTerminal)
  //     throw std::runtime_error("This shouldn't happen");
  // }
  void subdivide(NestNode&, const size_t, const size_t, const double, const double, size_t&);
};

#include "nest.tpp"
} // namespace brille
#endif
