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

#include <set>
#include <vector>
#include <array>
#include <tuple>
#include <utility>
#include <algorithm>
#include <omp.h>
#include "arrayvector.hpp"
#include "latvec.hpp"
#include "polyhedron.hpp"
#include "utilities.hpp"
#include "debug.hpp"
#include "triangulation_simple.hpp"
#include "interpolation_data.hpp"

#ifndef _NEST_H_
#define _NEST_H_

template<typename T,size_t N>
static bool none_negative(const std::array<T,N>& x){
  return !std::any_of(x.begin(), x.end(), [](T z){return z<0 && !approx_scalar(z,0.);});
}

class NestLeaf{
  std::array<size_t,4> vi;
  std::array<double,4> centre_radius;
  double volume_;
public:
  NestLeaf(): vi({{0,0,0,0}}), centre_radius({{0,0,0,0}}), volume_(0) {}
  // actually constructing the tetrahedra from, e.g., a Polyhedron object will
  // need to be done elsewhere
  explicit NestLeaf(
    const std::array<size_t,4>& vit,
    const std::array<double,4>& ci,
    const double vol
  ): vi(vit), centre_radius(ci), volume_(vol) {}
  //
  const std::array<size_t,4>& vertices(void) const { return vi;}
  //
  double volume(void) const {return volume_;}
  // double volume(const ArrayVector<double>& v) const {
    // return orient3d(v.data(vi[0]), v.data(vi[1]), v.data(vi[2]), v.data(vi[3]))/6.0;
  // }
  //
  std::array<double,4> weights(const ArrayVector<double>& v, const ArrayVector<double>& x) const {
    std::array<double,4> w{{-1,-1,-1,-1}};
    if (this->might_contain(x)){
      // double vol6 = this->volume(v)*6.0;
      double vol6 = volume_*6.0;
      w[0] = orient3d(x.data()     , v.data(vi[1]), v.data(vi[2]), v.data(vi[3])) / vol6;
      w[1] = orient3d(v.data(vi[0]), x.data()     , v.data(vi[2]), v.data(vi[3])) / vol6;
      w[2] = orient3d(v.data(vi[0]), v.data(vi[1]), x.data()     , v.data(vi[3])) / vol6;
      w[3] = orient3d(v.data(vi[0]), v.data(vi[1]), v.data(vi[2]), x.data()     ) / vol6;
    }
    return w;
  }
  bool contains(
    const ArrayVector<double>& v,
    const ArrayVector<double>& x,
    std::array<double,4>& w
  ) const {
    if (this->might_contain(x)){
      // double vol6 = this->volume(v)*6.0;
      double vol6 = volume_*6.0;
      w[0] = orient3d(x.data()     , v.data(vi[1]), v.data(vi[2]), v.data(vi[3])) / vol6;
      w[1] = orient3d(v.data(vi[0]), x.data()     , v.data(vi[2]), v.data(vi[3])) / vol6;
      w[2] = orient3d(v.data(vi[0]), v.data(vi[1]), x.data()     , v.data(vi[3])) / vol6;
      w[3] = orient3d(v.data(vi[0]), v.data(vi[1]), v.data(vi[2]), x.data()     ) / vol6;
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
  bool might_contain(const ArrayVector<double>& x) const {
    std::array<double,3> d;
    for (size_t i=0; i<3u; ++i) d[i] = x.getvalue(0,i) - centre_radius[i];
    double d2{0}, r2 = centre_radius[3]*centre_radius[3];
    for (size_t i=0; i<3u; ++i) d2 += d[i]*d[i];
    return ( d2 < r2 || approx_scalar(d2,r2) );
  }
};


class NestNode;
class NestNode{
  bool is_root_;
  NestLeaf boundary_;
  std::vector<NestNode> branches_;
public:
  explicit NestNode(bool ir=false): is_root_(ir), boundary_() {}
  explicit NestNode(const NestLeaf& b): is_root_(false), boundary_(b) {}
  NestNode(
    const std::array<size_t,4>& vit,
    const std::array<double,4>& ci,
    const double vol
  ): is_root_(false), boundary_(NestLeaf(vit,ci,vol)) {}
  bool is_root(void) const {return is_root_;}
  bool is_leaf(void) const {return !is_root_ && branches_.size()==0;}
  const NestLeaf& boundary(void) const {return boundary_;}
  const std::vector<NestNode>& branches(void) const {return branches_;}
  std::vector<NestNode>& branches(void) {return branches_;}
  // double volume(const ArrayVector<double>& v) const {return boundary_.volume(v);}
  double volume(void) const {return boundary_.volume();}
  template<typename... A> bool contains(A... args) {return boundary_.contains(args...);}
  template<typename... A> std::array<double,4> weights(A... args) {return boundary_.weights(args...);}
  std::vector<std::pair<size_t,double>> indices_weights(
    const ArrayVector<double>& v,
    // const std::vector<size_t>& m,
    const ArrayVector<double>& x
  ) const {
    std::array<double,4> w;
    // return __indices_weights(v,m,x,w);
    return __indices_weights(v,x,w);
  }
  std::vector<std::array<size_t,4>> tetrahedra(void) const {
    std::vector<std::array<size_t,4>> out;
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
  std::vector<std::pair<size_t,double>> __indices_weights(
    const ArrayVector<double>& v,
    // const std::vector<size_t>& m,
    const ArrayVector<double>& x,
    std::array<double,4>& w
  ) const {
    // This node is either the root (in which case it contains all tetrahedra)
    // or a tetrahedra below the root which definately contains the point x
    // In the second case w has been set for us by the calling function.
    if (this->is_leaf()){
      std::array<size_t,4> vi = boundary_.vertices();
      std::vector<std::pair<size_t,double>> iw;
      for (size_t i=0; i<4u; ++i) if (!approx_scalar(w[i], 0.))
        // iw.push_back(std::make_pair(m[vi[i]], w[i]));
        iw.push_back(std::make_pair(vi[i], w[i]));
      return iw;
    }
    // This is not a leaf node. So continue down the tree
    for (auto b: branches_){
      w = b.weights(v,x);
      // if (none_negative(w)) return b.__indices_weights(v,m,x,w);
      if (none_negative(w)) return b.__indices_weights(v,x,w);
    }
    std::vector<std::pair<size_t,double>> empty;
    return empty;
  }
};
template<class T, class S>
class Nest{
  NestNode root_;
  ArrayVector<double> vertices_;
  InterpolationData<T,S> data_;
  // std::vector<size_t> map_; // vertices holds *all* vertices but data_ only holds information for terminal vertices!
public:
  std::string tree_string(void) const {
    std::string tree = root_.to_string("",false);
    return tree;
  }
  // Build using maximum leaf volume
  Nest(const Polyhedron& p, const double vol, const size_t nb=5u):
    root_(true), vertices_({3u,0u})
  {
    this->construct(p, nb, vol);
    // this->make_all_to_terminal_map();
    size_t nvert = this->vertex_count();
    data_.initialize_permutation_table(nvert, root_.collect_keys(nvert));
  }
  // Build using desired leaf number density
  Nest(const Polyhedron& p, const size_t rho, const size_t nb=5u):
    root_(true), vertices_({3u,0u})
  {
    this->construct(p, nb, p.get_volume()/static_cast<double>(rho));
    // this->make_all_to_terminal_map();
    size_t nvert = this->vertex_count();
    data_.initialize_permutation_table(nvert, root_.collect_keys(nvert));
  }
  std::vector<bool> vertex_is_leaf(void) const {
    std::vector<bool> vert_is_term(vertices_.size(), false);
    for (auto tet: root_.tetrahedra()) for (auto idx: tet) vert_is_term[idx]=true;
    return vert_is_term;
  }
  const ArrayVector<double>& all_vertices(void) const {return vertices_;}
  // ArrayVector<double> vertices(void) const{ return vertices_.extract(this->vertex_is_leaf()); }
  ArrayVector<double> vertices(void) const{ return vertices_; }
  size_t vertex_count() const { return vertices_.size(); }
  std::vector<std::array<size_t,4>> tetrahedra(void) const {
    std::vector<std::array<size_t,4>> all_tet = root_.tetrahedra();
    // we need to adjust indexing to be into vertices instead of all_vertices
    /* (do this later) */
    return all_tet;
  }
  std::vector<std::pair<size_t,double>> indices_weights(const ArrayVector<double> &x) const {
    if (x.size()!=1u || x.numel()!=3u)
      throw std::runtime_error("The indices and weights can only be found for one point at a time.");
    // return root_.indices_weights(vertices_, map_, x);
    return root_.indices_weights(vertices_, x);
  }
  template<class R> unsigned check_before_interpolating(const ArrayVector<R>& x) const{
    unsigned int mask = 0u;
    if (this->data_.size()==0)
      throw std::runtime_error("The trellis must be filled before interpolating!");
    if (x.numel()!=3u)
      throw std::runtime_error("Nest requires x values which are three-vectors.");
    return mask;
  }
  std::tuple<ArrayVector<T>, ArrayVector<S>>
  interpolate_at(const ArrayVector<double>& x) const {
    this->check_before_interpolating(x);
    ArrayVector<T> vals(data_.values().numel(), x.size());
    ArrayVector<S> vecs(data_.vectors().numel(), x.size());
    for (size_t i=0; i<x.size(); ++i){
      // auto iw = root_.indices_weights(vertices_, map_, x.extract(i));
      auto iw = root_.indices_weights(vertices_, x.extract(i));
      data_.interpolate_at(iw, vals, vecs, i);
    }
    return std::make_tuple(vals, vecs);
  }
  std::tuple<ArrayVector<T>, ArrayVector<S>>
  interpolate_at(const ArrayVector<double>& x, const int threads) const {
    this->check_before_interpolating(x);
    omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
    // shared between threads
    ArrayVector<T> vals(data_.values().numel(), x.size());
    ArrayVector<S> vecs(data_.vectors().numel(), x.size());
    // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
    size_t unfound=0;
    long xsize = unsigned_to_signed<long, size_t>(x.size());
  #pragma omp parallel for default(none) shared(x, vals, vecs) reduction(+:unfound) firstprivate(xsize) schedule(dynamic)
    for (long si=0; si<xsize; ++si){
      size_t i = signed_to_unsigned<size_t, long>(si);
      // auto iw = root_.indices_weights(vertices_, map_, x.extract(i));
      auto iw = root_.indices_weights(vertices_, x.extract(i));
      if (iw.size()){
        data_.interpolate_at(iw, vals, vecs, i);
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
  const InterpolationData<T,S>& data(void) const {return data_;}
  template<typename... A> void replace_value_data(A... args) { data_.replace_value_data(args...); }
  template<typename... A> void replace_vector_data(A... args) { data_.replace_vector_data(args...); }
  template<typename... A> void set_value_cost_info(A... args) { data_.set_value_cost_info(args...); }
  template<typename... A> void set_vector_cost_info(A... args) {data_.set_vector_cost_info(args...);}
  template<template<class> class A>
  ArrayVector<double> debye_waller(const A<double>& Q, const std::vector<double>& M, const double t_K) const{
    return data_.debye_waller(Q,M,t_K);
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
#endif
