#include <vector>
#include <array>
#include <algorithm>
#include "arrayvector.h"

#ifndef _TRELLIS_H_
#define _TRELLIS_H_

class GeneralNode{
public:
  virtual size_t vertex_count() const;
  virtual bool indices_weights() const;
};
class CubeNode: public GeneralNode {
  std::array<size_t, 8> vertex_indices;
public:
  CubeNode(): vertex_indices({0,0,0,0,0,0,0,0}) {}
  CubeNode(const std::array<size_t,8>& vi): vertex_indices(vi) {}
  CubeNode(const std::vector<size_t>& vi): vertex_indices({0,0,0,0,0,0,0,0}) {
    if (vi.size() != 8) throw std::logic_error("TrellisNodeCube objects take 8 indices.");
    for (size_t i=0; i<8u; ++i) vertex_indices[i] = vi[i];
  }
  size_t vertex_count() const { return 8u;}
  // Calculate the relative interpolation weights for a point in the cube and
  // return the 1-8 finite weights and their associated vertex indices.
  // Return true if the point is in the cube, false otherwise.
  bool indices_weights(const ArrayVector<double>& vertices, const ArrayVector<double>& x, std::vector<size_t>& indices, std::vector<double>& weights);
};
class PolyNode: public GeneralNode {
  std::vector<std::array<size_t,4>> vertex_indices_per_tetrahedra;
public:
  PolyNode() {};
  // actually constructing the tetrahedra from, e.g., a Polyhedron object will
  // need to be done elsewhere
  PolyNode(const std::vector<std::array<size_t,4>>& vit): vertex_indices_per_tetrahedra(vit) {}
  // count-up the number fo unique vertices in the tetrahedra-triangulated polyhedron
  size_t vertex_count() const {
    std::vector<size_t> vi;
    for (auto tet: vertex_indices_per_tetrahedra)
    for (auto idx: tet)
    if (std::find(vi.begin(), vi.end(), idx)==vi.end()) vi.push_back(idx);
    return vi.size();
  }
  // this version needs to indentify *which* tetrahedron contains x first, then
  // calculate the relative interpolation weights for the 1-4 vertex indices.
  // Return true if the point is in one of the tetrahedra, otherwise false.
  bool indices_weights(const ArrayVector<double>& vertices, const ArrayVector<double>& x, std::vector<size_t>& indices, std::vector<double>& weights);
};
/* We might need an empty node, for now we can use a zero-tetrahedron PolyNode*/
// class EmptyNode: public GeneralNode {
// public:
//   EmptyNode() {};
//   size_t vertex_count() const { return 0; }
// };


class Trellis{
  ArrayVector<double> vertices_;
  std::vector<GeneralNode> nodes_;
  std::array<double,9> xyz_;
  std::array<std::vector<double>,3> boundaries_;
public:
  Trellis(): vertices_({3,0}), xyz_({1,0,0, 0,1,0, 0,0,1}) {
    std::vector<double> everything;
    everything.push_back(std::numeric_limits<double>::lowest());
    everything.push_back((std::numeric_limits<double>::max)());
    this->boundaries(everything, everything, everything);
  };
  Trellis(const std::array<double,9>& abc, const std::array<std::vector<double>,3>& bounds): vertices_({3,0}){
    this->xyz(abc);
    this->boundaries(bounds);
    this->node_count(); /*resizes the vector*/
  }
  size_t expected_vertex_count() const {
    size_t count = 1u;
    for (size_t i=0; i<3u; ++i) count *= boundaries_[i].size();
    return count;
  }
  size_t vertex_count() const {
    vertices_.size();
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
    std::array<size_t,3> s{1,0,0}, sz=this->size();
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
    vector_cross(xyz_.data()+6u, xyz_.data(), xyz_.data()+3u); // z = x × y
    return xyz_;
  }

  // Find the appropriate node for an arbitrary point:
  std::array<size_t,3> node_subscript(const std::array<double,3>& p) const {
    std::array<size_t,3> sub{0,0,0}, sz=this->size();
    for (size_t dim=0; dim<3u; ++dim){
      double p_dot_e = 0;
      for (size_t i=0; i<3u; ++i) p_dot_e += p[i]*xyz_[dim*3u + i];
      for (size_t i=0; i<sz[dim]; ++i) if ( p_dot_e < boundaries_[dim][i+1]) sub[dim] = i;
    }
    return sub;
  }
  std::array<size_t,3> node_subscript(const ArrayVector<double>& p) const {
    std::array<size_t,3> sub{0,0,0}, sz=this->size();
    for (size_t dim=0; dim<3u; ++dim){
      double p_dot_e = 0;
      for (size_t i=0; i<3u; ++i) p_dot_e += p.getvalue(0,i)*xyz_[dim*3u + i];
      for (size_t i=0; i<sz[dim]; ++i) if ( p_dot_e < boundaries_[dim][i+1]) sub[dim] = i;
    }
    return sub;
  }
  // find the node linear index for a point or leaf
  template <class T> size_t node_index(const T& p) const { return this->sub2idx(this->node_subscript(p)); }

  const GeneralNode& node(const size_t idx) const {
    if (idx < nodes_.size()) return nodes_[idx];
    throw std::domain_error("Out of bounds index for Trellis node");
  }
  GeneralNode& node(const size_t idx) {
    if (idx < nodes_.size()) return nodes_[idx];
    throw std::domain_error("Out of bounds index for Trellis node");
  }
  const GeneralNode& node(const size_t idx, const GeneralNode& n){
    if (idx < nodes_.size()){
      nodes_[idx] = n;
      return nodes_[idx];
    }
    throw std::domain_error("Out of bounds index for Trellis node");
  }
  // get the leaves located at a node by point-in-the-node
  const GeneralNode& node(const std::array<double,3>& p) const {
    return this->node(this->node_index(p));
  }
  const GeneralNode& node(const ArrayVector<double>& p) const {
    return this->node(this->node_index(p));
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

  size_t sub2idx(const std::array<size_t,3>& sub) const {
    return this->sub2idx(sub, this->span());
  }
  size_t sub2idx(const std::array<size_t,3>& sub, const std::array<size_t,3>& sp) const {
    size_t idx=0;
    for (size_t dim=0; dim<3u; ++dim) idx += sp[dim]*sub[dim];
    return idx;
  }
  std::array<size_t,3> idx2sub(const size_t idx, const std::array<size_t,3>& sp) const {
    std::array<size_t,3> sub{0,0,0};
    size_t rem{idx};
    for (size_t dim=3u; dim>0u; dim--){
      sub[dim] = rem/sp[dim];
      rem -= sub[dim]*sp[dim];
    }
    return sub;
  }

};

Trellis construct_trellis(const Polyhedron& poly, const double fraction=1.);


#endif
