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
#ifndef BRILLE_INTERPOLATOR_HPP_
#define BRILLE_INTERPOLATOR_HPP_
#include <vector>
#include <array>
#include <utility>
#include <mutex>
#include <cassert>
#include <functional>
#include <omp.h>
#include "array_latvec.hpp" // defines bArray
#include "phonon.hpp"
#include "permutation.hpp"
#include "permutation_table.hpp"
#include "approx.hpp"
#include "utilities.hpp"
namespace brille {

//template<class T> using CostFunction = std::function<typename CostTraits<T>::type(ind_t, T*, T*)>;
template<class T>
using CostFunction = std::function<double(brille::ind_t, const T*, const T*)>;

template<class T> struct is_complex {enum{value = false};};
template<class T> struct is_complex<std::complex<T>> {enum {value=true};};
// template<bool C, typename T> using enable_if_t = typename std::enable_if<C,T>::type;

enum class RotatesLike {
  Real, Reciprocal, Axial, Gamma
};

template<class T>
class Interpolator{
public:
  using ind_t = brille::ind_t;
  template<class Z> using element_t =std::array<Z,3>;
  using costfun_t = CostFunction<T>;
  using shape_t = std::vector<ind_t>;
private:
  bArray<T> data_;      //!< The stored Array2 of points indexed like the holding-Object's vertices
  shape_t shape_;       //!< The shape of the input Array (or Array2)
  element_t<ind_t> _elements; //!< The number of each element type per point and per mode
  RotatesLike rotlike_;   //!< How the elements of `data_` rotate
  element_t<double> _costmult; //!< The relative (multiplicative) cost for differences in each element type
  costfun_t _scalarfun; //!< A function to calculate differences between the scalars at two stored points
  costfun_t _vectorfun; //!< A function to calculate differences between the vectors at two stored points
  //costfun_t _matrixfun; //!< A function to calculate the differences between matrices at two stored points
public:
  explicit Interpolator(size_t scf_type=0, size_t vcf_type=0)
  : data_(0,0), _elements({{0,0,0}}), rotlike_{RotatesLike::Real}, _costmult({{1,1,1}})
  {
    this->set_cost_info(scf_type, vcf_type);
  }
  Interpolator(costfun_t scf, costfun_t vcf)
  : data_(0,0), _elements({{0,0,0}}), rotlike_{RotatesLike::Real},
    _costmult({{1,1,1}}), _scalarfun(scf), _vectorfun(vcf)
  {}
  Interpolator(bArray<T>& d, shape_t sh, element_t<ind_t> el, RotatesLike rl)
  : data_(d), shape_(sh), _elements(el), rotlike_{rl}, _costmult({{1,1,1}})
  {
    this->set_cost_info(0,0);
    this->check_elements();
  }
  Interpolator(bArray<T>& d, shape_t sh, element_t<ind_t> el, RotatesLike rl, size_t csf, size_t cvf, element_t<double> wg)
  : data_(d), shape_(sh), _elements(el), rotlike_{rl}, _costmult(wg)
  {
    this->set_cost_info(csf, cvf);
    this->check_elements();
  }
  // use the Array2<T>(const Array<T>&) constructor
  Interpolator(brille::Array<T>& d, element_t<ind_t> el, RotatesLike rl)
  : data_(d), shape_(d.shape()), _elements(el), rotlike_{rl}, _costmult({{1,1,1}})
  {
    this->set_cost_info(0,0);
    this->check_elements();
  }
  // use the Array2<T>(const Array<T>&) constructor
  Interpolator(brille::Array<T>& d, element_t<ind_t> el, RotatesLike rl, size_t csf, size_t cvf, element_t<double> wg)
  : data_(d), shape_(d.shape()), _elements(el), rotlike_{rl}, _costmult(wg)
  {
    this->set_cost_info(csf, cvf);
    this->check_elements();
  }
  //
  void setup_fake(const ind_t sz, const ind_t br){
    data_ = bArray<T>(sz, br);
    shape_ = {sz, br};
    _elements = {1u,0u,0u};
  }
  //
  void set_cost_info(const int scf, const int vcf){
    switch (scf){
      default:
      this->_scalarfun = [](ind_t n, const T* i, const T* j){
        double s{0};
        for (ind_t z=0; z<n; ++z) s += brille::utils::magnitude(i[z]-j[z]);
        return s;
      };
    }
    switch (vcf){
      case 1:
      debug_update("selecting brille::utils::vector_distance");
      this->_vectorfun = [](ind_t n, const T* i, const T* j){
        return brille::utils::vector_distance(n, i, j);
      };
      break;
      case 2:
      debug_update("selecting 1-brille::utils::vector_product");
      this->_vectorfun = [](ind_t n, const T* i, const T* j){
        return 1-brille::utils::vector_product(n, i, j);
      };
      break;
      case 3:
      debug_update("selecting brille::utils::vector_angle");
      this->_vectorfun = [](ind_t n, const T* i, const T* j){
        return brille::utils::vector_angle(n, i, j);
      };
      break;
      case 4:
      debug_update("selecting brille::utils::hermitian_angle");
      this->_vectorfun = [](ind_t n, const T* i, const T* j){
         return brille::utils::hermitian_angle(n,i,j);
      };
      break;
      default:
      debug_update("selecting sin**2(brille::utils::hermitian_angle)");
      // this->_vectorfun = [](ind_t n, T* i, T* j){return std::abs(std::sin(brille::utils::hermitian_angle(n, i, j)));};
      this->_vectorfun = [](ind_t n, const T* i, const T* j){
        auto sin_theta_H = std::sin(brille::utils::hermitian_angle(n, i, j));
        return sin_theta_H*sin_theta_H;
      };
    }
  }
  void set_cost_info(const int scf, const int vcf, const element_t<double>& elcost){
    _costmult = elcost;
    this->set_cost_info(scf, vcf);
  }
  //
  ind_t size(void) const {return data_.size(0);}
  ind_t branches(void) const {
    /*
    The data Array is allowed to be anywhere from 1 to 5 dimensional.
    Possible shapes are:
      (P,)            - one branch with one scalar per point
      (P, X)          - one branch with X elements per point
      (P, B, Y)       - B branches with Y elements per branch per point
      (P, B, V, 3)    - B branches with V 3-vectors per branch per point
      (P, B, M, 3, 3) - B branches with M (3,3) matrices per branch per point
    In the two-dimensional case there is a potential ambiguity in that X counts
    both the number of branches and the number of elements.
      X = B*Y
    where Y is the sum of scalar, vector, and matrix elements per branch
      Y = S + 3*V + 9*M
    */
    ind_t nd = shape_.size();
    ind_t b = nd>1 ? shape_[1] : 1u;
    if (2 == nd){
      ind_t y = std::accumulate(_elements.begin(), _elements.end(), ind_t(0));
      // zero-elements is the special (initial) case → 1 scalar per branch
      if (y > 0) b /= y;
    }
    return b;
  }
  bool only_vector_or_matrix(void) const {
    // if V or M can not be deduced directly from the shape of data_
    // then this Interpolator represents mixed data
    // purely-scalar data is also classed as 'mixed' for our purposes
    ind_t nd = shape_.size();
    if (5u == nd && 3u == shape_[4] && 3u == shape_[3]) return true;
    if (4u == nd && 3u == shape_[2]) return true;
    if (nd < 4u) return false; // (P,), (P,X), (P,B,Y)
    std::string msg = "Interpolator can not handle a {";
    for (auto x: shape_) msg += " " + std::to_string(x);
    msg += " } data array";
    throw std::runtime_error(msg);
  }
  const bArray<T>& data(void) const {return data_;}
  shape_t shape(void) const {return shape_;};
  brille::Array<T> array(void) const {return brille::Array<T>(data_,shape_);}
  element_t<ind_t> elements(void) const {return _elements;}
  //
  template<class... Args>
  void interpolate_at(Args... args) const {
      this->interpolate_at_mix(args...);
  }
  //
  template<class R, class RotT>
  bool rotate_in_place(bArray<T>& x,
                       const LQVec<R>& q,
                       const RotT& rt,
                       const PointSymmetry& ps,
                       const std::vector<size_t>& r,
                       const std::vector<size_t>& invr,
                       const int nth) const {
    switch (rotlike_){
      case RotatesLike::Real:       return this->rip_real(x,ps,r,invr,nth);
      case RotatesLike::Axial:      return this->rip_axial(x,ps,r,invr,nth);
      case RotatesLike::Reciprocal: return this->rip_recip(x,ps,r,invr,nth);
      case RotatesLike::Gamma:      return this->rip_gamma(x,q,rt,ps,r,invr,nth);
      default: throw std::runtime_error("Impossible RotatesLike value!");
    }
  }
  //
  RotatesLike rotateslike() const { return rotlike_; }
  RotatesLike rotateslike(const RotatesLike a) {
    rotlike_ = a;
    return rotlike_;
  }
  // Replace the data within this object.
  template<class I>
  void replace_data(
      const bArray<T>& nd,
      const shape_t sh,
      const std::array<I,3>& ne,
      const RotatesLike rl = RotatesLike::Real)
  {
    data_ = nd;
    shape_ = sh;
    rotlike_ = rl;
    // convert the elements datatype as necessary
    if (ne[1]%3)
      throw std::logic_error("Vectors must have 3N elements per branch");
    if (ne[2]%9)
      throw std::logic_error("Matrices must have 9N elements per branch");
    for (size_t i=0; i<3u; ++i) _elements[i] = static_cast<ind_t>(ne[i]);
    this->check_elements();
  }
  template<class I>
  void replace_data(
      const brille::Array<T>& nd,
      const std::array<I,3>& ne,
      const RotatesLike rl = RotatesLike::Real)
  {
    data_ = bArray<T>(nd);
    shape_ = nd.shape();
    rotlike_ = rl;
    // convert the elements datatype as necessary
    if (ne[1]%3)
      throw std::logic_error("Vectors must have 3N elements per branch");
    if (ne[2]%9)
      throw std::logic_error("Matrices must have 9N elements per branch");
    for (size_t i=0; i<3u; ++i) _elements[i] = static_cast<ind_t>(ne[i]);
    this->check_elements();
  }
  // Replace the data in this object without specifying the data shape or its elements
  // this variant is necessary since the template specialization above can not have a default value for the elements
  template<template<class> class A>
  void replace_data(const A<T>& nd){
    return this->replace_data(nd, element_t<ind_t>({{0,0,0}}));
  }
  ind_t branch_span() const { return this->branch_span(_elements);}
  //
  std::string to_string() const {
    std::string str= "{ ";
    for (auto s: shape_) str += std::to_string(s) + " ";
    str += "} data";
    auto b = this->branches();
    if (b){
      str += " with " + std::to_string(b) + " mode";
      if (b>1) str += "s";
    }
    auto n = std::count_if(_elements.begin(), _elements.end(), [](ind_t a){return a>0;});
    if (n){
      str += " of ";
      std::array<std::string,3> types{"scalar", "vector", "matrix"};
      for (size_t i=0; i<3u; ++i) if (_elements[i]) {
        str += std::to_string(_elements[i]) + " " + types[i];
        if (--n>1) str += ", ";
        if (1==n) str += " and ";
      }
      str += " element";
      if (this->branch_span()>1) str += "s";
    }
    return str;
  }

  template<class S>
  void add_cost(const ind_t, const ind_t, std::vector<S>&, bool) const;

  template<typename I>
  bool any_equal_modes(const I idx) const {
    return this->any_equal_modes_(static_cast<ind_t>(idx), this->branches(), this->branch_span());
  }
  size_t bytes_per_point() const {
    size_t n_elements = data_.numel()/data_.size(0);
    return n_elements * sizeof(T);
  }
private:
  void check_elements(void){
    // check the input for correctness
    ind_t x = this->branch_span(_elements);
    switch (shape_.size()) {
      case 1u: // 1 scalar per branch per point
        if (0u == x) x = _elements[0] = 1u;
        if (x > 1u) throw std::runtime_error("1-D data must represent one scalar per point!") ;
        break;
      case 2u: // (P, B*Y)
        if (0u == x) x = _elements[0] = shape_[1]; // one branch with y scalars per point
        if (shape_[1] % x)
          throw std::runtime_error("2-D data requires an integer number of branches!");
        break;
      case 3u: // (P, B, Y)
        if (0u == x) x = _elements[0] = shape_[2];
        if (shape_[2] != x)
          throw std::runtime_error("3-D data requires that the last dimension matches the specified number of elements!");
        break;
      case 4u: // (P, B, V, 3)
        if (3u != shape_[3])
          throw std::runtime_error("4-D data can only be 3-vectors");
        if (0u == x) x = _elements[1] = shape_[2]*3u;
        if (shape_[2]*3u != x)
          throw std::runtime_error("4-D data requires that the last two dimensions match the specified number of vector elements!");
        break;
      case 5u: // (P, B, M, 3, 3)
        if (3u != shape_[3] || 3u != shape_[4])
          throw std::runtime_error("5-D data can only be matrices");
        if (0u == x) x = _elements[2] = shape_[2]*9u;
        if (shape_[2]*9u != x)
          throw std::runtime_error("5-D data requires the last three dimensions match the specified number of matrix elements!");
        break;
      default: // higher dimensions not (yet) supported
        throw std::runtime_error("Interpolator data is expected to be 1- to 5-D");
    }
  }
  bool any_equal_modes_(const ind_t idx, const ind_t b_, const ind_t s_) {
    // since we're probably only using this when the data is provided and
    // most eigenproblem solvers sort their output by eigenvalue magnitude it is
    // most-likely for mode i and mode i+1 to be equal.
    // ∴ search (i,j), (i+1,j+1), (i+2,j+2), ..., i ∈ (0,N], j ∈ (1,N]
    // for each j = i+1, i+2, i+3, ..., i+N-1
    if (b_ < 2) return false;
    // data_ is always 2D: (N,1), (N,B), or (N,Y)
    for (ind_t offset=1; offset < b_; ++offset)
    for (ind_t i=0, j=offset; j < b_; ++i, ++j)
    if (brille::approx::vector(s_, data_.ptr(idx, i*s_), data_.ptr(idx, j*s_))) return true;
    // no matches
    return false;
  }
  template<typename I> ind_t branch_span(const std::array<I,3>& e) const {
    return static_cast<ind_t>(e[0])+static_cast<ind_t>(e[1])+static_cast<ind_t>(e[2]);
  }
  element_t<ind_t> count_scalars_vectors_matrices(void) const {
    element_t<ind_t> no{_elements[0], _elements[1]/3u, _elements[2]/9u};
    return no;
  }
  // the 'mixed' variants of the rotate_in_place implementations
  bool rip_real(bArray<T>&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  bool rip_recip(bArray<T>&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  bool rip_axial(bArray<T>&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  template<class R>
  bool rip_gamma_complex(bArray<T>&, const LQVec<R>&, const GammaTable&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  template<class R, class S=T>
  enable_if_t<is_complex<S>::value, bool>
  rip_gamma(bArray<T>& x, const LQVec<R>& q, const GammaTable& gt, const PointSymmetry& ps, const std::vector<size_t>& r, const std::vector<size_t>& ir, const int nth) const{
    return rip_gamma_complex(x, q, gt, ps, r, ir, nth);
  }
  template<class R, class S=T>
  enable_if_t<!is_complex<S>::value, bool>
  rip_gamma(bArray<T>&, const LQVec<R>&, const GammaTable&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const{
    throw std::runtime_error("RotatesLike == Gamma requires complex valued data!");
  }

  // interpolate_at_*
  void interpolate_at_mix(const std::vector<std::vector<ind_t>>&, const std::vector<ind_t>&, const std::vector<double>&, bArray<T>&, const ind_t, const bool) const;
  void interpolate_at_mix(const std::vector<std::vector<ind_t>>&, const std::vector<std::pair<ind_t,double>>&, bArray<T>&, const ind_t, const bool) const;
};

#include "interpolator2_at.tpp"
#include "interpolator2_axial.tpp"
#include "interpolator2_cost.tpp"
#include "interpolator2_gamma.tpp"
#include "interpolator2_real.tpp"
#include "interpolator2_recip.tpp"

} // namespace brille
#endif
