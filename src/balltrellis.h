#include <vector>
#include <array>
#include "arrayvector.h"

#ifndef _BALLTRELLIS_H_
#define _BALLTRELLIS_H_

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
  bool fuzzy_contains(const ArrayVector<double>& x) const {
    double d=0;
    for (size_t i=0; i<3u; ++i) d += (x.getvalue(0,i)-_centre[i])*(x.getvalue(0,i)-_centre[i]);
    return (d < _squared_radius || approx_scalar(d, _squared_radius));
  }
  bool fuzzy_contains(const std::array<double,3>& x) const {
    double d=0;
    for (size_t i=0; i<3u; ++i) d += (x[i]-_centre[i])*(x[i]-_centre[i]);
    return (d < _squared_radius || approx_scalar(d, _squared_radius));
  }
};


class Trellis{
  std::vector<std::vector<TrellisLeaf>> _nodes;
  std::array<double,9> _xyz;
  std::array<std::vector<double>,3> _boundaries;
public:
  Trellis(): _xyz({1,0,0, 0,1,0, 0,0,1}) {
    std::vector<double> everything;
    everything.push_back(std::numeric_limits<double>::lowest());
    everything.push_back((std::numeric_limits<double>::max)());
    this->boundaries(everything, everything, everything);
  };
  Trellis(const std::array<double,3>& abc, const std::array<std::vector<double>,3>& bounds){
    this->xyz(abs);
    this->boundaries(bounds);
    this->node_count(); /*resizes the vector*/
  }
  Trellis(const std::array<double,3>& abc,
          const std::array<std::vector<double>,3>& bounds,
          const std::vector<TrellisLeaf>& leaves){
    this->xyz(abs);
    this->boundaries(bounds);
    this->node_count(); /*resizes the vector*/
    this->add_leaves(leaves);
  }
  size_t node_count() const {
    size_t count = 1u;
    for (size_t i=0; i<3u; ++i) count *= _boundaries[i].size()-1;
    _nodes.resize(count);
    return count;
  }
  std::array<size_t,3> size() const {
    std::array<size_t,3> s;
    for (size_t i=0; i<3u; ++i) s[i] = _boundaries[i].size()-1;
    return s;
  }
  std::array<size_t,3> span() const {
    std::array<size_t,3> s{1,0,0}, z=this->size();
    for (size_t i=1; i<3; ++i) s[i] = z[i-1]*s[i-1];
    return s;
  }
  size_t boundaries(const std::vector<double>& xb, const std::vector<double>& yb, const std::vector<double>& zb){
    if (xb.size()>1) _boundaries[0] = xb;
    if (yb.size()>1) _boundaries[1] = yb;
    if (zb.size()>1) _boundaries[2] = zb;
    return this->node_count();
  }
  size_t xboundaries(const std::vector<double>& xb){ return this->boundaries(xb, _boundaries[1], _boundaries[2]);}
  size_t yboundaries(const std::vector<double>& yb){ return this->boundaries(_boundaries[0], yb, _boundaries[2]);}
  size_t zboundaries(const std::vector<double>& zb){ return this->boundaries(_boundaries[0], _boundaries[1], zb);}
  std::array<double,3> x() const { std::array<double,3> out{_xyz[0],_xyz[1],_xyz[2]}; return out; }
  std::array<double,3> y() const { std::array<double,3> out{_xyz[3],_xyz[4],_xyz[5]}; return out; }
  std::array<double,3> z() const { std::array<double,3> out{_xyz[6],_xyz[7],_xyz[8]}; return out; }
  std::array<double,3> x(const std::array<double,3>& n) { for (size_t i=0; i<3u; ++i) _xyz[i   ] = n[i]; return this->x();}
  std::array<double,3> y(const std::array<double,3>& n) { for (size_t i=0; i<3u; ++i) _xyz[i+3u] = n[i]; return this->y();}
  std::array<double,3> z(const std::array<double,3>& n) { for (size_t i=0; i<3u; ++i) _xyz[i+6u] = n[i]; return this->z();}
  const std::array<double,9>& xyz() const { return _xyz;}
  const std::array<double,9>& xyz(const std::array<double,9>& n){
    _xyz = n;
    return _xyz;
  }
  const std::array<double,9>& xyz(const std::array<double,3>& nx, const std::array<double,3>& ny, const std::array<double,3>& nz){
    for (size_t i=0; i<3u; ++i){
      _xyz[i   ] = nx[i];
      _xyz[i+3u] = ny[i];
      _xyz[i+6u] = nz[i];
    }
    return _xyz;
  }
  const std::array<double,9>& xyz(const std::array<double,3>& nx, const std::array<double,3>& ny){
    for (size_t i=0; i<3u; ++i){
      _xyz[i   ] = nx[i];
      _xyz[i+3u] = ny[i];
    }
    vector_cross(_xyz.data()+6u, _xyz.data(), _xyz.data()+3u); // z = x × y
    return _xyz;
  }

  // Find the appropriate node for an arbitrary point:
  std::array<size_t,3> node_subscript(const std::array<double,3>& p) const {
    std::array<size_t,3> sub{0,0,0}, z=this->size();
    double dot;
    for (size_t dim=0; dim<3u; ++dim){
      dot = 0;
      for (size_t i=0; i<3u; ++i) dot += p[i]*_xyz[dim*3u + i];
      for (size_t i=0; i<z[i]; ++i) if ( dot < _boundaries[dim][i+1]) sub[dim] = i;
    }
    return sub;
  }
  std::array<size_t,3> node_subscript(const ArrayVector<double>& p) const {
    std::array<size_t,3> sub{0,0,0}, z=this->size();
    double dot;
    for (size_t dim=0; dim<3u; ++dim){
      dot = 0;
      for (size_t i=0; i<3u; ++i) dot += p.getvalue(0,i)*_xyz[dim*3u + i];
      for (size_t i=0; i<z[i]; ++i) if ( dot < _boundaries[dim][i+1]) sub[dim] = i;
    }
    return sub;
  }
  // or leaf, which has a centre and radius
  std::array<size_t,3> node_subscript(const TrellisLeaf& l) const {
    std::array<double,3> p = l.centre();
    double dot, r = l.radius();
    std::array<size_t,3> sub{0,0,0}, z=this->size();
    for (size_t dim=0; dim<3u; ++dim){
      dot = 0;
      for (size_t i=0; i<3u; ++i) dot += p[i]*_xyz[3u*dim + i];
      for (size_t i=0; i<z[i]; ++i) if (dot+r < _boundaries[dim][i+1]) sub[dim] = i;
    }
    return sub;
  }
  // find the node linear index for a point or leaf
  template <class T>
  size_t node_index(const T& p) const {
    std::array<size_t,3> sp, sb;
    sp = this->span();
    sb = this->node_subscript(p);
    size_t idx=0;
    for (size_t dim=0; dim<3u; ++dim) idx += sp[dim]*sb[dim];
    return idx;
  }

  // get the leaves located at a node
  const std::vector<TrellisLeaf>& node_leaves(const size_t idx) const {
    if (idx < _nodes.size()) return _nodes[idx];
    throw std::domain_error("Out of bounds index for Trellis' nodes");
  }
  const std::vector<TrellisLeaf>& node_leaves(const size_t idx, const std::vector<TrellisLeaf>& l) {
    if (idx < _nodes.size()){
      _nodes[idx] = l;
      return _nodes[idx];
    }
    throw std::domain_error("Out of bounds index for Trellis' nodes");
  }
  // get the leaves located at a node by point-in-the-node
  const std::vector<TrellisLeaf>& node_leaves(const std::array<double,3>& p) const {
    return this->node_leaves(this->node_index(p));
  }
  const std::vector<TrellisLeaf>& node_leaves(const ArrayVector<double>& p) const {
    return this->node_leaves(this->node_index(p));
  }
  // add a leaf to the trellis:
  bool add_leaf(const TrellisLeaf& l){
    size_t idx = this->node_index(l);
    _nodes[idx].push_back(l);
    return true;
  }
  // add a number of leaves to the trellis:
  bool add_leaves(const std::vector<TrellisLeaf>& leaves){
    size_t idx;
    for (auto leaf: leaves){
      idx = this->node_index(leaf);
      _nodes[idx].push_back(leaf);
    }
    return true;
  }
};

Trellis construct_trellis(const std::vector<TrellisLeaf>& leaves, const size_t Nxyz=5);
Trellis construct_trellis(const std::vector<TrellisLeaf>& leaves, const std::array<size_t,3>& Nxyz);


#endif
