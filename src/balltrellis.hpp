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

#ifndef BRILLE_BALLTRELLIS_H_
#define BRILLE_BALLTRELLIS_H_
#include <vector>
#include <array>
#include <algorithm>
#include "array_latvec.hpp"
#include "approx.hpp"
namespace brille {

class TrellisLeaf{
  std::array<double,3> _centre;
  double _squared_radius;
  size_t _index;
public:
  ~TrellisLeaf() = default;
  TrellisLeaf(): _squared_radius(0.), _index(0) {};
  TrellisLeaf(const std::array<double,3>& p, const double r, const size_t i): _centre(p), _squared_radius(r*r), _index(i) {};
  size_t index(const size_t idx){
    _index = idx;
    return _index;
  }
  size_t index() const { return _index; }
  double radius() const { return std::sqrt(_squared_radius); }
  const std::array<double,3>& centre() const {return _centre; }
  bool fuzzy_contains(const bArray<double>& x) const {
    assert(x.size() == 3u);
    double d=0;
    for (size_t i=0; i<3u; ++i){
      double tmp = x[i]-_centre[i];
      d += tmp*tmp;
    }
    return (d < _squared_radius || brille::approx::scalar(d, _squared_radius));
  }
  bool fuzzy_contains(const std::array<double,3>& x) const {
    double d=0;
    for (size_t i=0; i<3u; ++i) d += (x[i]-_centre[i])*(x[i]-_centre[i]);
    return (d < _squared_radius || brille::approx::scalar(d, _squared_radius));
  }
};


class Trellis{
  std::vector<std::vector<TrellisLeaf>> nodes_;
  std::array<double,9> xyz_;
  std::array<std::vector<double>,3> boundaries_;
  double max_leaf_radius_;
public:
  Trellis(): xyz_({1,0,0, 0,1,0, 0,0,1}), max_leaf_radius_(0.) {
    std::vector<double> everything;
    everything.push_back(std::numeric_limits<double>::lowest());
    everything.push_back((std::numeric_limits<double>::max)());
    this->boundaries(everything, everything, everything);
  };
  Trellis(const std::array<double,9>& abc, const std::array<std::vector<double>,3>& bounds): max_leaf_radius_(0.) {
    this->xyz(abc);
    this->boundaries(bounds);
    this->node_count(); /*resizes the vector*/
  }
  Trellis(const std::array<double,9>& abc,
          const std::array<std::vector<double>,3>& bounds,
          const std::vector<TrellisLeaf>& leaves): max_leaf_radius_(0.){
    this->xyz(abc);
    this->boundaries(bounds);
    this->node_count(); /*resizes the vector*/
    this->add_leaves(leaves);
  }
  size_t node_count() {
    size_t count = 1u;
    for (size_t i=0; i<3u; ++i) count *= boundaries_[i].size()-1;
    nodes_.resize(count);
    return count;
  }
  std::array<size_t,3> size() const {
    std::array<size_t,3> s;
    for (size_t i=0; i<3u; ++i) s[i] = boundaries_[i].size()-1;
    return s;
  }
  std::array<size_t,3> span() const {
    std::array<size_t,3> s{{1,0,0}}, sz=this->size();
    for (size_t i=1; i<3; ++i) s[i] = sz[i-1]*s[i-1];
    return s;
  }
  size_t boundaries(const std::array<std::vector<double>,3>& xyzb){
    if (std::all_of(xyzb.begin(), xyzb.end(), [](const std::vector<double>& a){ return a.size()>1; }))
      boundaries_ = xyzb;
    return this->node_count();
  }
  size_t boundaries(const std::vector<double>& xb, const std::vector<double>& yb, const std::vector<double>& zb){
    if (xb.size()>1) boundaries_[0] = xb;
    if (yb.size()>1) boundaries_[1] = yb;
    if (zb.size()>1) boundaries_[2] = zb;
    return this->node_count();
  }
  size_t xboundaries(const std::vector<double>& xb){ return this->boundaries(xb, boundaries_[1], boundaries_[2]);}
  size_t yboundaries(const std::vector<double>& yb){ return this->boundaries(boundaries_[0], yb, boundaries_[2]);}
  size_t zboundaries(const std::vector<double>& zb){ return this->boundaries(boundaries_[0], boundaries_[1], zb);}
  std::array<double,3> x() const { std::array<double,3> out{xyz_[0],xyz_[1],xyz_[2]}; return out; }
  std::array<double,3> y() const { std::array<double,3> out{xyz_[3],xyz_[4],xyz_[5]}; return out; }
  std::array<double,3> z() const { std::array<double,3> out{xyz_[6],xyz_[7],xyz_[8]}; return out; }
  std::array<double,3> x(const std::array<double,3>& n) { for (size_t i=0; i<3u; ++i) xyz_[i   ] = n[i]; return this->x();}
  std::array<double,3> y(const std::array<double,3>& n) { for (size_t i=0; i<3u; ++i) xyz_[i+3u] = n[i]; return this->y();}
  std::array<double,3> z(const std::array<double,3>& n) { for (size_t i=0; i<3u; ++i) xyz_[i+6u] = n[i]; return this->z();}
  const std::array<double,9>& xyz() const { return xyz_;}
  const std::array<double,9>& xyz(const std::array<double,9>& n){
    xyz_ = n;
    return xyz_;
  }
  const std::array<double,9>& xyz(const std::array<double,3>& nx, const std::array<double,3>& ny, const std::array<double,3>& nz){
    for (size_t i=0; i<3u; ++i){
      xyz_[i   ] = nx[i];
      xyz_[i+3u] = ny[i];
      xyz_[i+6u] = nz[i];
    }
    return xyz_;
  }
  const std::array<double,9>& xyz(const std::array<double,3>& nx, const std::array<double,3>& ny){
    for (size_t i=0; i<3u; ++i){
      xyz_[i   ] = nx[i];
      xyz_[i+3u] = ny[i];
    }
    brille::utils::vector_cross(xyz_.data()+6u, xyz_.data(), xyz_.data()+3u); // z = x × y
    return xyz_;
  }

  // Find the appropriate node for an arbitrary point:
  std::array<size_t,3> node_subscript(const std::array<double,3>& p) const {
    std::array<size_t,3> sub{{0,0,0}}, sz=this->size();
    for (size_t dim=0; dim<3u; ++dim){
      double p_dot_e = 0;
      for (size_t i=0; i<3u; ++i) p_dot_e += p[i]*xyz_[dim*3u + i];
      for (size_t i=0; i<sz[dim]; ++i) if ( p_dot_e < boundaries_[dim][i+1]) sub[dim] = i;
    }
    return sub;
  }
  std::array<size_t,3> node_subscript(const bArray<double>& p) const {
    assert(p.size() == 3u);
    std::array<size_t,3> sub{{0,0,0}}, sz=this->size();
    for (size_t dim=0; dim<3u; ++dim){
      double p_dot_e = 0;
      for (size_t i=0; i<3u; ++i) p_dot_e += p[i]*xyz_[dim*3u + i];
      for (size_t i=0; i<sz[dim]; ++i) if ( p_dot_e < boundaries_[dim][i+1]) sub[dim] = i;
    }
    return sub;
  }
  // or leaf, which has a centre and radius
  std::array<size_t,3> node_subscript(const TrellisLeaf& l) const {
    std::array<double,3> p = l.centre();
    double r = l.radius();
    std::array<size_t,3> sub{{0,0,0}}, sz=this->size();
    for (size_t dim=0; dim<3u; ++dim){
      double p_dot_e = 0;
      for (size_t i=0; i<3u; ++i) p_dot_e += p[i]*xyz_[3u*dim + i];
      for (size_t i=0; i<sz[dim]; ++i) if (p_dot_e+r < boundaries_[dim][i+1]) sub[dim] = i;
    }
    return sub;
  }
  // find the node linear index for a point or leaf
  template <class T> size_t node_index(const T& p) const { return this->sub2idx(this->node_subscript(p)); }

  // get the leaves located at a node
  const std::vector<TrellisLeaf>& node_leaves(const size_t idx) const {
    if (idx < nodes_.size()) return nodes_[idx];
    throw std::domain_error("Out of bounds index for Trellis' nodes");
  }
  const std::vector<TrellisLeaf>& node_leaves(const size_t idx, const std::vector<TrellisLeaf>& l) {
    if (idx < nodes_.size()){
      nodes_[idx] = l;
      return nodes_[idx];
    }
    throw std::domain_error("Out of bounds index for Trellis' nodes");
  }
  // get the leaves located at a node by point-in-the-node
  const std::vector<TrellisLeaf>& node_leaves(const std::array<double,3>& p) const {
    return this->node_leaves(this->node_index(p));
  }
  const std::vector<TrellisLeaf>& node_leaves(const bArray<double>& p) const {
    return this->node_leaves(this->node_index(p));
  }
  // add a leaf to the trellis:
  bool add_leaf(const TrellisLeaf& l){
    size_t idx = this->node_index(l);
    nodes_[idx].push_back(l);
    return true;
  }
  // add a number of leaves to the trellis:
  bool add_leaves(const std::vector<TrellisLeaf>& leaves){
    for (auto leaf: leaves){
      size_t idx = this->node_index(leaf);
      if (leaf.radius() > max_leaf_radius_) max_leaf_radius_ = leaf.radius(); // keep track of this on a per-node basis?
      nodes_[idx].push_back(leaf);
    }
    return true;
  }
  /* Get a vector of node indices to search when trying to find a point in a
     leaf. We only need to check against leaves with a node-distance no larger
     than 2√3 times the maximum leaf radius -- since nodes further away than
     this are incapable of holding leaves which reach points in this node.
     */
  template<class T> std::vector<size_t> nodes_to_search(const T& p) const {
    std::array<size_t,3> tsub, psub = this->node_subscript(p);
    std::array<size_t,3> xtnt = this->size();
    std::array<size_t,3> xspn = this->span();
    std::vector<size_t> to_search;
    // include all nodes which are no more than the maximum leaf diameter away
    // We only need to search in the increasing-subscripts direction due to how
    // the leaves were assigned to the trellis nodes to begin with.
    for (tsub[0]=psub[0]; tsub[0]<xtnt[0]; ++tsub[0])
    for (tsub[1]=psub[1]; tsub[1]<xtnt[1]; ++tsub[1])
    for (tsub[2]=psub[2]; tsub[2]<xtnt[2]; ++tsub[2]){
      // if (this->node_distance(psub, tsub) <= 3.5*max_leaf_radius_) // a bit more than sqrt(3)*d to account for body diagonal
      if (this->node_close_enough(psub, tsub))
        to_search.push_back(this->sub2idx(tsub, xspn));
    }
    // Assuming we're most likely to find the leaf at a node close to where
    // the point is located, sort the found nodes by their distance away
    std::sort(to_search.begin(), to_search.end(),
      [psub, xspn, this](const size_t i, const size_t j){
        return this->node_distance(psub,this->idx2sub(i,xspn)) < this->node_distance(psub,this->idx2sub(j,xspn));
      }
    );
    return to_search;
  }
  std::string to_string(void) const {
    std::string str = "(";
    for (auto i: this->size()) str += " " + std::to_string(i);
    str += " )";
    std::array<size_t,3> min_max_tot = this->node_leaves_min_max_tot();
    str += " {";
    str += std::to_string(min_max_tot[0]) + "--" + std::to_string(min_max_tot[1]);
    str += " : " + std::to_string(min_max_tot[2]/static_cast<double>(nodes_.size()));
    str += "}";
    return str;
  }
private:
  std::array<size_t,3> node_leaves_min_max_tot(void) const {
    std::array<size_t,3> mmt{(std::numeric_limits<size_t>::max)(), 0u, 0u};
    size_t leavescount;
    for (auto leaves: nodes_){
      leavescount = leaves.size();
      if (leavescount < mmt[0]) mmt[0] = leavescount;
      if (leavescount > mmt[1]) mmt[1] = leavescount;
      mmt[2] += leavescount;
    }
    return mmt;
  }
  size_t sub2idx(const std::array<size_t,3>& sub) const {
    return this->sub2idx(sub, this->span());
  }
  size_t sub2idx(const std::array<size_t,3>& sub, const std::array<size_t,3>& sp) const {
    size_t idx=0;
    for (size_t dim=0; dim<3u; ++dim) idx += sp[dim]*sub[dim];
    return idx;
  }
  std::array<size_t,3> idx2sub(const size_t idx, const std::array<size_t,3>& sp) const {
    std::array<size_t,3> sub{{0,0,0}};
    size_t rem{idx};
    for (size_t dim=3u; dim--;){
      sub[dim] = rem/sp[dim];
      rem -= sub[dim]*sp[dim];
    }
    return sub;
  }
  double node_distance(const std::array<size_t,3>& i, const std::array<size_t,3>& j) const {
    double d{0.};
    for (size_t dim=0; dim<3u; ++dim) d += (boundaries_[dim][i[dim]+1]-boundaries_[dim][j[dim]+1])*(boundaries_[dim][i[dim]+1]-boundaries_[dim][j[dim]+1]);
    return std::sqrt(d);
  }
  bool node_close_enough(const std::array<size_t,3>& i, const std::array<size_t,3>& j) const {
    double halfdist = this->node_distance(i,j)/2.0;
    if (halfdist < max_leaf_radius_) return true;
    return brille::approx::scalar(halfdist, max_leaf_radius_);
  }
};

Trellis construct_trellis(const std::vector<TrellisLeaf>& leaves, const double fraction=1.);

} // end namespace brille
#endif
