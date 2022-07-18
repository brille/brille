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
/*! \file
    \author Greg Tucker
    \brief A class holding a triangulated tetrahedral mesh and data for interpolation
*/
#include <deque>
// #include <set>
// #include <vector>
// #include <array>
// #include <tuple>
// #include <utility>
// #include <algorithm>
// #include <omp.h>
// #include "array.hpp"
// #include "array2.hpp"
// #include "utilities.hpp"
// #include "debug.hpp"
#include "triangulation_simple.hpp"
#include "interpolatordual.hpp"
#include "approx_config.hpp"
#include "polyhedron_flex.hpp"

namespace brille {

template<typename T,size_t N>
static bool none_negative(const std::array<T,N>& x, const T t, const int n){
  auto is_negative = [t, n](T z){return z<0 && !brille::approx_float::scalar(z, 0., t, t, n);};
  auto any_negative = std::any_of(x.begin(), x.end(), is_negative);
  return !any_negative;
}

/*! \brief A single tetrahedron

The NestLeaf is a single tetrahedron which holds its center of mass position,
circumsphere radius, and volume along with the indices of its vertex positions.

If provided with the list of all vertex positions and a test point it can
determine whether the point is inside the tetrahedron and, if required, with
what weights the vertices should be combined to produce the linear interpolation
at the test point.
*/
class NestLeaf{
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
  [[nodiscard]] const std::array<ind_t,4>& vertices() const { return vi;}
  //
  [[nodiscard]] double volume() const {return volume_;}
  //
  [[nodiscard]] std::array<double,4> weights(const bArray<double>& v, const bArray<double>& x, const double t, const int n) const {
    std::array<double,4> w{{-1,-1,-1,-1}};
    if (this->might_contain(x,t,n)){
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
    std::array<double,4>& w,
    const double t,
    const int n
  ) const {
    if (this->might_contain(x,t,n)){
      double vol6 = volume_*6.0;
      w[0] = orient3d(x.ptr(0)    , v.ptr(vi[1]), v.ptr(vi[2]), v.ptr(vi[3])) / vol6;
      w[1] = orient3d(v.ptr(vi[0]), x.ptr(0)    , v.ptr(vi[2]), v.ptr(vi[3])) / vol6;
      w[2] = orient3d(v.ptr(vi[0]), v.ptr(vi[1]), x.ptr(0)    , v.ptr(vi[3])) / vol6;
      w[3] = orient3d(v.ptr(vi[0]), v.ptr(vi[1]), v.ptr(vi[2]), x.ptr(0)    ) / vol6;
      return none_negative(w, t, n);
    }
    return false;
  }
  //
  [[nodiscard]] std::string to_string() const {
    std::string msg = "[";
    for (auto i: vi) msg += " " + std::to_string(i);
    msg += " ]";
    return msg;
  }
private:
  [[nodiscard]] bool might_contain(const bArray<double>& x, const double t, const int n) const {
    std::array<double,3> d{0,0,0};
    d[0] = x.val(0,0) - centre_radius[0];
    d[1] = x.val(0,1) - centre_radius[1];
    d[2] = x.val(0,2) - centre_radius[2];
    double d2{0}, r2 = centre_radius[3]*centre_radius[3];
    for (size_t i=0; i<3u; ++i) d2 += d[i]*d[i];
    return ( d2 < r2 || brille::approx_float::scalar(d2,r2,t,t,n) );
  }
};

/*! \brief Part of a tetrahedron hierarchy

The NestNode may be any level of a tetrahedron hierarchy and always defines a
tetrahedral region of its domain.
If the NestNode is terminal it contains no further NetNodes, otherwise it must
hold a list of two or more smaller NestNodes which fill the space (or nest)
within it.
*/
class NestNode{
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
  //! Return whether this NestNode is the most coarse tetrahedron
  [[nodiscard]] bool is_root() const {return is_root_;}
  //! Return whether this NestNode is the most precise tetrahedron
  [[nodiscard]] bool is_leaf() const {return !is_root_ && branches_.size()==0;}
  //! Return a reference to the NestLeaf single tetrahedron which bounds this NestNode's domain
  [[nodiscard]] const NestLeaf& boundary() const {return boundary_;}
  //! Return a constant reference to the NestNodes which fill this NestNode's domain
  [[nodiscard]] const std::vector<NestNode>& branches() const {return branches_;}
  //! Return a reference to the NestNodes which fill this NestNode's domain
  std::vector<NestNode>& branches() {return branches_;}
  //! Return the volume of this NestNode's domain
  [[nodiscard]] double volume() const {return boundary_.volume();}
  //! Determine if a test point is inside of this NestNode's domain
  template<typename... A> bool contains(A... args) {return boundary_.contains(args...);}
  //! Determine the relative interpolation weights for the vertices of this NestNode if it contains a test point
  template<typename... A> std::array<double,4> weights(A... args) {return boundary_.weights(args...);}
  //! Return the indexed vertices in addition to the relative interpolation weights
  [[nodiscard]] std::vector<std::pair<ind_t,double>> indices_weights(
    const bArray<double>& v, const bArray<double>& x, const double t_tol, const int n_tol
  ) const
  {
    /* What if more than one branch 'contains' the point -- which happens any
     * time the point is on an internal face of the tetrahedron nest.
     * In this case we should follow all branches until we find all containing
     * leaves, then pick the 'best' one.
     * */
    std::deque<NestNode> work;
    for (const auto& b: branches_) work.emplace_back(b);
    std::vector<std::vector<std::pair<ind_t, double>>> solutions;
    while (!work.empty()) {
      auto w = work.front().weights(v, x, t_tol, n_tol);
      if (none_negative(w, t_tol, n_tol)) {
        if (work.front().is_leaf()) {
          auto vi = work.front().boundary().vertices();
          std::vector<std::pair<ind_t, double>> solution;
          for (size_t i = 0; i < 4u; ++i)
              solution.emplace_back(vi[i], w[i]);
          solutions.push_back(solution);
        } else {
          for (const auto &b : work.front().branches())
            work.emplace_back(b);
        }
      }
      work.pop_front();
    }
    // No solution found:
    if (solutions.empty()) return {};
    // More than one sol -- find the 'best' one:
    size_t best{0};
    if (solutions.size() > 1u) {
      for (size_t i=0; i<solutions.size(); ++i){
        bool better{true};
        for (auto & j : solutions[i]) if (j.second <= 0) better = false;
        if (better) best = i;
      }
    }
    // (Attempt to) deal with slightly-negative weights by reducing maximum
    auto is_zero = [t_tol, n_tol](double x){return approx_float::scalar(x, 0., t_tol, t_tol, n_tol);};
    auto sol = solutions[best];
    for (size_t i=0; i<4u; ++i){
      // check this in the loop incase sol[i] modified in previous loop
      if (is_zero(sol[i].second)) {
        size_t max_at{0};
        // same for max_at; previous loop may have modified maximum weight
        for (size_t j=0; j<4u; ++j){
          if (!is_zero(sol[j].second) && sol[j].second > sol[max_at].second) {
            max_at = j;
          }
        }
        sol[max_at].second += sol[i].second;
      }
    }
    // Extract just the non-zero weight point(s)
    std::vector<std::pair<ind_t, double>> out;
    for (const auto & p: sol) if (!is_zero(p.second)) out.push_back(p);
    return out;
  }
  //! Pull together all tetrahedra vertex indices at or below this level of the hierarchy
  [[nodiscard]] std::vector<std::array<ind_t,4>> tetrahedra() const {
    std::vector<std::array<ind_t,4>> out;
    if (this->is_leaf()) out.push_back(boundary_.vertices());
    for (auto b: branches_) for (auto v: b.tetrahedra()) out.push_back(v);
    return out;
  }
  [[nodiscard]] std::string to_string(const std::string& prefix, const bool not_last) const {
    std::string msg = prefix;
    msg += is_root_ ? "───┐" : not_last ? "├──" : "└──";
    if (!is_root_) msg += boundary_.to_string();
    msg += "\n";
    for (size_t i=0; i<branches_.size(); ++i)
      msg += branches_[i].to_string(prefix+(not_last?"|  ":"   "),i+1!=branches_.size());
    return msg;
  }
  //! Pull together all PermutationTable keys for the connected vertices at or below this level of the hierarchy
  [[nodiscard]] std::set<size_t> collect_keys(const size_t nv) const {
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
};

/*!
\brief A nested triangulated tetrahedral mesh with eigenvalue and eigenvector data

One way of dividing three dimensional space is to fill it with a tiling of
tetrehedra. Each tetrahedra has four vertices and four triangular faces.
As each face can be shared with one other tetrahedra, each tetrahedron has up
to four neighbours (tetrahedra on the surface of a space will have one fewer
neighbour per surface). There is no limit to how many tetrahedra any vertex
can contribute to.

If one or more values are defined for every vertex in the tetrahedral mesh then
it can be used to perform linear interpolation for any arbitrary point within
the bound space.

With no guaranteed ordering of the tethrahedra finding which tetrahedra contains
the interpolation point can require testing all tetrahedra for inclusion.
As such a simple tetrahedral mesh is not well suited for *fast* interpolation.
In order to overcome this limitation, this class uses a hierarchy of nested
tetrahedra arranged in a tree to limit the number of inclusion tests required
during interpolation.
If a point is within a tetrahedron at a given depth of the hierarchy then a list
of tetrahedra that it might be in at the next lower level is available.
These next-lower tetrahedra have the property that they fill the higher-level
tetrahedra.
*/
template<class DataValues, class DataVectors, class VertexComponents, template<class> class VertexType>
class Nest{
public:
  using class_t = Nest<DataValues, DataVectors, VertexComponents, VertexType>;
  using data_t = DualInterpolator<DataValues, DataVectors>;
  using vert_t = VertexType<VertexComponents>;
  using root_t = NestNode;
  using approx_t = approx_float::Config;
  using poly_t = polyhedron::Poly<VertexComponents, VertexType>;
private:
  root_t root_;
  vert_t vertices_;
  data_t data_;
  approx_t approx_;
  // std::vector<size_t> map_; // vertices holds *all* vertices but data_ only holds information for terminal vertices!
public:
  [[nodiscard]] approx_t approx_config() const {return approx_;}
  [[nodiscard]] std::string tree_string() const {
    std::string tree = root_.to_string("",false);
    return tree;
  }

  template<class Z> Nest(const poly_t& p, const Z x, const ind_t nb=5u):
  root_(true), vertices_(0u, 3u)
  {
    this->pre_construct(p, x, nb, approx_float::config);
  }
  template<class Z> Nest(const poly_t& p, const Z x, const ind_t nb, const approx_t a):
  root_(true), vertices_(0u, 3u)
  {
    this->pre_construct(p, x, nb, a);
  }
  void pre_construct(const poly_t & p, const ind_t rho, const ind_t nb, const approx_t a){
    approx_ = a;
    this->construct(p, nb, p.volume() / static_cast<double>(rho));
    auto n_vert = vertex_count();
    data_.initialize_permutation_table(n_vert, root_.collect_keys(n_vert));
  }
  void pre_construct(const poly_t & p, const double volume, const ind_t nb, const approx_t a){
    approx_ = a;
    this->construct(p, nb, volume);
    auto n_vert = vertex_count();
    data_.initialize_permutation_table(n_vert, root_.collect_keys(n_vert));
  }

  [[nodiscard]] std::vector<bool> vertex_is_leaf() const {
    std::vector<bool> vert_is_term(vertices_.size(0), false);
    for (auto tet: root_.tetrahedra()) for (auto idx: tet) vert_is_term[idx]=true;
    return vert_is_term;
  }
  const vert_t& all_vertices() const {return vertices_;}
  vert_t vertices() const{ return vertices_; }
  [[nodiscard]] brille::ind_t vertex_count() const { return vertices_.size(0); }
  [[nodiscard]] std::vector<std::array<ind_t,4>> tetrahedra() const {
    std::vector<std::array<ind_t,4>> all_tet = root_.tetrahedra();
    // we need to adjust indexing to be into vertices instead of all_vertices
    /* (do this later) */
    return all_tet;
  }
  std::vector<std::pair<ind_t,double>>
  indices_weights(const vert_t &x) const {
    if (x.ndim()!=2 || x.size(0) != 1u || x.size(1) != 3u)
      throw std::runtime_error("The indices and weights can only be found for one point at a time.");
    auto t = approx_.reciprocal<double>();
    auto n = approx_.digit();
    return root_.indices_weights(vertices_, x, t, n);
  }
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
  std::tuple<brille::Array<DataValues>, brille::Array<DataVectors>>
  interpolate_at(const vert_t& x) const {
    this->check_before_interpolating(x);
    auto valsh = data_.values().shape();
    auto vecsh = data_.values().shape();
    valsh[0] = x.size(0);
    vecsh[0] = x.size(0);
    brille::Array<DataValues> vals(valsh);
    brille::Array<DataVectors> vecs(vecsh);
    // vals and vecs are row-ordered contiguous by default, so we can create
    // mutable data-sharing Array2 objects for use with
    // Interpolator2::interpolate_at through the constructor:
    brille::Array2<DataValues> vals2(vals);
    brille::Array2<DataVectors> vecs2(vecs);
    auto t = approx_.reciprocal<double>();
    auto n = approx_.digit();
    for (ind_t i=0; i<x.size(0); ++i){
      // auto iw = root_.indices_weights(vertices_, map_, x.extract(i));
      auto iw = root_.indices_weights(vertices_, x.extract(i), t, n);
      data_.interpolate_at(iw, vals2, vecs2, i);
    }
    return std::make_tuple(vals, vecs);
  }
  std::tuple<brille::Array<DataValues>, brille::Array<DataVectors>>
  interpolate_at(const vert_t& x, const int threads) const {
    this->check_before_interpolating(x);
    omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
    // not used in parallel region
    auto valsh = data_.values().shape();
    auto vecsh = data_.vectors().shape();
    valsh[0] = x.size(0);
    vecsh[0] = x.size(0);
    // shared between threads
    brille::Array<DataValues> vals(valsh);
    brille::Array<DataVectors> vecs(vecsh);
    // vals and vecs are row-ordered contiguous by default, so we can create
    // mutable data-sharing Array2 objects for use with
    // Interpolator2::interpolate_at through the constructor:
    brille::Array2<DataValues> vals2(vals);
    brille::Array2<DataVectors> vecs2(vecs);
    // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
    ind_t unfound=0;
    auto t = approx_.reciprocal<double>();
    auto n = approx_.digit();
    auto xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
  #pragma omp parallel for default(none) shared(x, vals2, vecs2) reduction(+:unfound) firstprivate(xsize, t, n) schedule(dynamic)
    for (long long si=0; si<xsize; ++si){
      auto i = brille::utils::s2u<ind_t, long long>(si);
      auto iw = root_.indices_weights(vertices_, x.extract(i), t, n);
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
  const data_t& data() const {return data_;}
  template<typename... A> void replace_data(A... args) {data_.replace_data(args...);}
  template<typename... A> void replace_value_data(A... args) { data_.replace_value_data(args...); }
  template<typename... A> void replace_vector_data(A... args) { data_.replace_vector_data(args...); }
  template<typename... A> void set_value_cost_info(A... args) { data_.set_value_cost_info(args...); }
  template<typename... A> void set_vector_cost_info(A... args) {data_.set_vector_cost_info(args...);}
  //! Return the number of bytes used per Q point
  [[nodiscard]] size_t bytes_per_point() const {return data_.bytes_per_point(); }
  void sort() {data_.sort();}
private:
  void construct(const poly_t&, size_t, double);
  // void make_all_to_terminal_map(void) {
  //   std::vector<bool> vit = this->vertex_is_leaf();
  //   size_t nTerminal = std::count(vit.begin(), vit.end(), true);
  //   map_.clear();
  //   size_t idx{0};
  //   for (bool t: vit) map_.push_back(t ? idx++ : nTerminal);
  //   if (idx != nTerminal)
  //     throw std::runtime_error("This shouldn't happen");
  // }
  void subdivide(root_t&, size_t, size_t, double, double, ind_t&);
};

#include "nest.tpp"
} // namespace brille
#endif
