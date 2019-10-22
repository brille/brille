#include <vector>
#include <array>
#include "arrayvector.h"

#ifndef _BALLTREE_H_
#define _BALLTREE_H_

class BallLeaf{
  std::array<double,3> _centre;
  double _radius;
  size_t _index;
public:
  ~BallLeaf() = default;
  BallLeaf(): _radius(0.), _index(0) {};
  BallLeaf(const std::array<double,3>& p, const double r, const size_t i): _centre(p), _radius(r), _index(i) {};
  BallLeaf(const BallLeaf& other){
    this->_centre = other.centre();
    this->_radius = other.radius();
    this->_index = other.index();
  }
  BallLeaf& operator=(const BallLeaf& other){
    this->_centre = other.centre();
    this->_radius = other.radius();
    this->_index = other.index();
    return *this;
  }
  size_t index(const size_t idx){
    _index = idx;
    return _index;
  }
  size_t index() const { return _index; }
  double radius() const { return _radius; }
  const std::array<double,3>& centre() const {return _centre; }
  bool fuzzy_contains(const std::array<double,3>& x) const {
    double d=0;
    for (size_t i=0; i<3u; ++i) d += (x[i]-_centre[i])*(x[i]-_centre[i]);
    d = std::sqrt(d) - _radius;
    return (d < 0 || approx_scalar(d, 0.));
  }
};

class BallNode{
  std::vector<BallNode> _children;
  std::vector<BallLeaf> _leaves;
  std::array<double,3> _centre;
  double _radius;
public:
  ~BallNode() = default;
  BallNode(): _radius(0) {};
  BallNode(const std::array<double,3>& p, const double r): _centre(p), _radius(r){};
  BallNode(const BallNode& other){
    this->_centre = other.centre();
    this->_radius = other.radius();
    // force deep copying
    for (auto child: other.children()) this->_children.push_back(BallNode(child));
    for (auto leaf: other.leaves()) this->_leaves.push_back(BallLeaf(leaf));
  }
  BallNode& operator=(const BallNode& other){
    this->_centre = other.centre();
    this->_radius = other.radius();
    // force deep copying
    for (auto child: other.children()) this->_children.push_back(BallNode(child));
    for (auto leaf: other.leaves()) this->_leaves.push_back(BallLeaf(leaf));
    return *this;
  }
  // child related methods
  size_t count_children() const { return _children.size(); }
  std::vector<BallNode> containing_children(const std::array<double,3>& x) const {
    std::vector<BallNode> out;
    for (auto child: _children) if (child.contains(x)) out.push_back(child);
    return out;
  }
  const std::vector<BallNode>& children() const { return _children; }
  void children(const std::vector<BallNode>& c){
    _children.clear();
    for (auto child: c) _children.push_back(child);
  }
  void addchild(const BallNode& c){ _children.push_back(c); }
  // this node related methods
  size_t count_leaves() const { return _leaves.size(); }
  const std::vector<BallLeaf>& leaves() const { return _leaves; }
  void leaves(const std::vector<BallLeaf>& l){
    _leaves.clear();
    for (auto leaf: l) _leaves.push_back(leaf);
  }
  void addleaf(const BallLeaf& l){ _leaves.push_back(l); }
  std::vector<BallLeaf> containing_leaves(const std::array<double,3>& x) const {
    std::vector<BallLeaf> out;
    for (auto leaf: _leaves) if (leaf.fuzzy_contains(x)) out.push_back(leaf);
    return out;
  }
  std::vector<size_t> containing_leaf_indexes(const std::array<double,3>& x) const {
    std::vector<size_t> out;
    for (auto leaf: _leaves) if (leaf.fuzzy_contains(x)) out.push_back(leaf.index());
    return out;
  }
  double radius() const { return _radius; }
  const std::array<double,3>& centre() const {return _centre; }
  bool contains(const std::array<double,3>& x) const {
    double d=0;
    for (size_t i=0; i<3u; ++i) d += (x[i]-_centre[i])*(x[i]-_centre[i]);
    return std::sqrt(d) <= _radius;
  }
  // (recursive) full tree related methods
  std::vector<BallLeaf> all_containing_leaves(const std::array<double,3>& x) const {
    std::vector<BallLeaf> acl, ccl;
    for (auto child: _children){
      ccl = child.all_containing_leaves(x);
      for (auto l: ccl) acl.push_back(l);
    }
    for (auto leaf: _leaves) if (leaf.fuzzy_contains(x)) acl.push_back(leaf);
    return acl;
  }
  std::vector<size_t> all_containing_leaf_indexes(const ArrayVector<double>& x) const {
    std::array<double,3> z;
    for (size_t i=0; i<3u; ++i) z[i] = x.getvalue(0,i);
    return this->all_containing_leaf_indexes(z);
  }
  std::vector<size_t> all_containing_leaf_indexes(const std::array<double,3>& x) const {
    std::vector<size_t> aci, cci;
    for (auto child: _children){
      cci = child.all_containing_leaf_indexes(x);
      for (auto i: cci) aci.push_back(i);
    }
    for (auto leaf: _leaves) if (leaf.fuzzy_contains(x)) aci.push_back(leaf.index());
    return aci;
  }
  std::string to_string() const {
    std::vector<std::string> levels;
    levels.push_back("level 0:");
    this->fill_strings(0, levels);
    std::string all_levels = levels[0];
    for (size_t i=1; i<levels.size(); ++i) all_levels += "\n"+levels[i];
    return all_levels;
  }
  void fill_strings(const size_t i, std::vector<std::string>& levels) const {
    if (_leaves.size()){
      for (auto leaf: _leaves) levels[i] += " " + std::to_string(leaf.index());
      levels[i] += " | ";
    }
    if (_children.size()){
      if (levels.size() <= i+1){
        levels.resize(i+2);
        levels[i+1] = "level " + std::to_string(i+1) +":";
      }
      for (auto child: _children) child.fill_strings(i+1, levels);
    }
  }
};

BallNode construct_balltree(const std::vector<BallLeaf>& leaves, const size_t max_depth=5);
BallNode construct_ballnode(const std::vector<BallLeaf>& leaves, std::array<double,3>& v);
bool bifurcate_balltree(BallNode& root, const std::vector<BallLeaf>& leaves, const size_t max_depth, const std::array<double,3>& at, const std::array<double,3>& along);


#endif
