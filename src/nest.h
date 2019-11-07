#include <vector>
#include <array>
#include <algorithm>
#include <omp.h>
#include "arrayvector.h"
#include "latvec.h"
#include "polyhedron.h"
#include "unsignedtosigned.h"
#include "debug.h"
#include "triangulation_simple.h"
#include "interpolation_data.h"

#ifndef _NEST_H_
#define _NEST_H_



class NestLeaf{
  std::array<size_t,4> vi;
public:
  NestLeaf(){}
  // actually constructing the tetrahedra from, e.g., a Polyhedron object will
  // need to be done elsewhere
  explicit NestLeaf(const std::array<size_t,4>& vit): vi(vit) {}
  const std::array<size_t,4>& vertices(void) const { return vi;}
  bool indices_weights(
    const ArrayVector<double>& v,
    const std::vector<size_t>& map,
    const ArrayVector<double>& x,
    std::vector<size_t>& indices,
    std::vector<double>& weights
  ) const {
    indices.clear();
    weights.clear();
    std::array<double,4> w{0,0,0,0};
    if (this->contains(v, x, w))
    for (size_t i=0; i<4u; ++i)
    if (!approx_scalar(w[i], 0.)){
      indices.push_back(map[vi[i]]);
      weights.push_back(w[i]);
    }
    return indices.size()>0;
  }
  double volume(const ArrayVector<double>& v) const {
    return orient3d(v.data(vi[0]), v.data(vi[1]), v.data(vi[2]), v.data(vi[3]))/6.0;
  }
  bool contains(
    const ArrayVector<double>& v,
    const ArrayVector<double>& x,
    std::array<double,4>& w
  ) const {
    if (!this->might_contain(v,x)) return false; // is this any faster (or right?)
    double vol6 = this->volume(v)*6.0;
    w[0] = orient3d(x.data()     , v.data(vi[1]), v.data(vi[2]), v.data(vi[3])) / vol6;
    w[1] = orient3d(v.data(vi[0]), x.data()     , v.data(vi[2]), v.data(vi[3])) / vol6;
    w[2] = orient3d(v.data(vi[0]), v.data(vi[1]), x.data()     , v.data(vi[3])) / vol6;
    w[3] = orient3d(v.data(vi[0]), v.data(vi[1]), v.data(vi[2]), x.data()     ) / vol6;
    if (std::any_of(w.begin(), w.end(), [](double z){return z < 0. && !approx_scalar(z, 0.);}))
      return false;
    return true;
  }
  bool contains(
    const ArrayVector<double>& v,
    const ArrayVector<double>& x
  ) const {
    std::array<double,4> w{0,0,0,0};
    return this->contains(v,x,w);
  }
  std::string to_string(void) const {
    std::string msg = "[";
    for (auto i: vi) msg += " " + std::to_string(i);
    msg += " ]";
    return msg;
  }
private:
  bool might_contain(const ArrayVector<double>& v, const ArrayVector<double>& x) const {
    double is = insphere(v.data(vi[0]),v.data(vi[1]),v.data(vi[2]),v.data(vi[3]),x.data());
    return is > 0 || approx_scalar(is, 0.);
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
  explicit NestNode(const std::array<size_t,4>& vit): is_root_(false), boundary_(NestLeaf(vit)) {}
  bool is_root(void) const {return is_root_;}
  bool is_leaf(void) const {return branches_.size()==0;}
  const NestLeaf& boundary(void) const {return boundary_;}
  const std::vector<NestNode>& branches(void) const {return branches_;}
  std::vector<NestNode>& branches(void) {return branches_;}
  double volume(const ArrayVector<double>& v) const {return boundary_.volume(v);}
  template<typename... A> bool contains(A... args) {return boundary_.contains(args...);}
  bool indices_weights(
    const ArrayVector<double>& v,
    const std::vector<size_t>& m,
    const ArrayVector<double>& x,
    std::vector<size_t>& i,
    std::vector<double>&w
  ) const {
    for (auto b: branches_) if (b.contains(v,x)) return b.indices_weights(v,m,x,i,w);
    if (!is_root_ && this->is_leaf()) return boundary_.indices_weights(v,m,x,i,w);
    return false;
  }
  std::vector<std::array<size_t,4>> tetrahedra(void) const {
    std::vector<std::array<size_t,4>> out;
    if (!is_root_ && this->is_leaf()) out.push_back(boundary_.vertices());
    else for (auto b: branches_) for (auto v: b.tetrahedra()) out.push_back(v);
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
};
template<typename T>
class Nest{
  NestNode root_;
  ArrayVector<double> vertices_;
  InterpolationData<T> data_;
  std::vector<size_t> map_; // vertices holds *all* vertices but data_ only holds information for terminal vertices!
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
    this->make_all_to_terminal_map();
  }
  // Build using desired leaf number density
  Nest(const Polyhedron& p, const size_t rho, const size_t nb=5u):
    root_(true), vertices_({3u,0u})
  {
    this->construct(p, nb, p.get_volume()/static_cast<double>(rho));
    this->make_all_to_terminal_map();
  }
  std::vector<bool> vertex_is_leaf(void) const {
    std::vector<bool> vert_is_term(vertices_.size(), false);
    for (auto tet: root_.tetrahedra()) for (auto idx: tet) vert_is_term[idx]=true;
    return vert_is_term;
  }
  const ArrayVector<double>& all_vertices(void) const {return vertices_;}
  ArrayVector<double> vertices(void) const{ return vertices_.extract(this->vertex_is_leaf()); }
  std::vector<std::array<size_t,4>> tetrahedra(void) const {
    std::vector<std::array<size_t,4>> all_tet = root_.tetrahedra();
    // we need to adjust indexing to be into vertices instead of all_vertices
    /* (do this later) */
    return all_tet;
  }
  bool indices_weights(const ArrayVector<double>& x, std::vector<size_t>& i, std::vector<double>& w) const {
    if (x.size()!=1u || x.numel()!=3u)
      throw std::runtime_error("The indices and weights can only be found for one point at a time.");
    return root_.indices_weights(vertices_, map_, x, i, w);
  }
  template<class R> unsigned check_before_interpolating(const ArrayVector<R>& x) const{
    unsigned int mask = 0u;
    if (this->data_.size()==0)
      throw std::runtime_error("The trellis must be filled before interpolating!");
    if (x.numel()!=3u)
      throw std::runtime_error("Nest requires x values which are three-vectors.");
    return mask;
  }
  ArrayVector<T> interpolate_at(const ArrayVector<double>& x, const int threads=1) const {
    this->check_before_interpolating(x);
    omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
    // shared between threads
    ArrayVector<T> out(data_.numel(), x.size());
    // private to each thread
    std::vector<size_t> indices;
    std::vector<double> weights;
    // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
    long xsize = unsigned_to_signed<long, size_t>(x.size());
  #pragma omp parallel for shared(x, out) private(indices, weights)
    for (long si=0; si<xsize; ++si){
      size_t i = signed_to_unsigned<size_t, long>(si);
      if (!root_.indices_weights(vertices_, map_, x.extract(i), indices, weights))
        throw std::runtime_error("Point not found in Nest");
      data_.interpolate_at(indices, weights, out, i);
    }
    return out;
  }
  const InterpolationData<T>& data(void) const {return data_;}
  template<typename... A> void replace_data(A... args){ data_.replace_data(args...);}
  template<template<class> class A>
  ArrayVector<double> debye_waller(const A<double>& Q, const std::vector<double>& M, const double t_K) const{
    return data_.debye_waller(Q,M,t_K);
  }
private:
  void construct(const Polyhedron&, const size_t, const double);
  void make_all_to_terminal_map(void) {
    std::vector<bool> vit = this->vertex_is_leaf();
    size_t nTerminal = std::count(vit.begin(), vit.end(), true);
    map_.clear();
    size_t idx{0};
    for (bool t: vit) map_.push_back(t ? idx++ : nTerminal);
    if (idx != nTerminal)
      throw std::runtime_error("This shouldn't happen");
  }
  void subdivide(NestNode&, const size_t, const size_t, const double, const double, size_t&);
};

#include "nest.hpp"
#endif
