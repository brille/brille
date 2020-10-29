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
#include <vector>
#include <array>
#include <utility>
#include <mutex>
#include <cassert>
#include <functional>
#include <omp.h>
#include "array.hpp"
#include "array_latvec.hpp" // defines bArray
#include "phonon.hpp"
#include "permutation.hpp"
#include "permutation_table.hpp"
#include "approx.hpp"

#ifndef _INNER_INTERPOLATION_DATA_HPP_
#define _INNER_INTERPOLATION_DATA_HPP_

//template<class T> using CostFunction = std::function<typename CostTraits<T>::type(ind_t, T*, T*)>;
template<class T>
using CostFunction = std::function<double(brille::ind_t, const T*, const T*)>;

template<class T> struct is_complex {enum{value = false};};
template<class T> struct is_complex<std::complex<T>> {enum {value=true};};
// template<bool C, typename T> using enable_if_t = typename std::enable_if<C,T>::type;

enum class RotatesLike {
  Real, Reciprocal, Axial, Gamma
};

template<class T, class P>
class InnerInterpolationData{
public:
  using ind_t = brille::ind_t;
  template<class R, class S> using data_t = brille::Array<R>;
  using shape_t = typename data_t<T,P>::shape_t;
  using ElementsType = std::array<brille::ind_t,3>;
  using ElementsCost = std::array<double,3>;
private:
  data_t<T,P> data_;      //!< The stored eigenvalue-like Array indexed like the holding-Object's vertices
  ElementsType elements_; //!< The number of scalars, vector elements, and matrix elements per `data_` array
  RotatesLike rotlike_;   //!< How the elements of `data_` rotate
  ElementsCost costs_;    //!< The cost assigned to each type for equivalent mode assignment
  CostFunction<T> scalar_cost_function;
  CostFunction<T> vector_cost_function;
public:
  explicit InnerInterpolationData(size_t scf_type=0, size_t vcf_type=0)
  : data_(0,0), elements_({{0,0,0}}), rotlike_{RotatesLike::Real}, costs_({{1,1,1}})
  {
    this->set_cost_info(scf_type, vcf_type);
  }
  InnerInterpolationData(CostFunction<T> scf, CostFunction<T> vcf)
  : data_(0,0), elements_({{0,0,0}}), rotlike_{RotatesLike::Real},
    costs_({{1,1,1}}), scalar_cost_function(scf), vector_cost_function(vcf)
  {}
  //
  void setup_fake(const ind_t sz, const ind_t br){
    data_ = data_t<T,P>(sz, br);
    elements_ = {1u,0u,0u};
  }
  //
  void set_cost_info(const int scf, const int vcf){
    switch (scf){
      default:
      this->scalar_cost_function = [](ind_t n, const T* i, const T* j){
        double s{0};
        for (ind_t z=0; z<n; ++z) s += brille::utils::magnitude(i[z]-j[z]);
        return s;
      };
    }
    switch (vcf){
      case 1:
      debug_update("selecting brille::utils::vector_distance");
      this->vector_cost_function = [](ind_t n, const T* i, const T* j){
        return brille::utils::vector_distance(n, i, j);
      };
      break;
      case 2:
      debug_update("selecting 1-brille::utils::vector_product");
      this->vector_cost_function = [](ind_t n, const T* i, const T* j){
        return 1-brille::utils::vector_product(n, i, j);
      };
      break;
      case 3:
      debug_update("selecting brille::utils::vector_angle");
      this->vector_cost_function = [](ind_t n, const T* i, const T* j){
        return brille::utils::vector_angle(n, i, j);
      };
      break;
      case 4:
      debug_update("selecting brille::utils::hermitian_angle");
      this->vector_cost_function = [](ind_t n, const T* i, const T* j){
         return brille::utils::hermitian_angle(n,i,j);
      };
      break;
      default:
      debug_update("selecting sin**2(brille::utils::hermitian_angle)");
      // this->vector_cost_function = [](ind_t n, T* i, T* j){return std::abs(std::sin(brille::utils::hermitian_angle(n, i, j)));};
      this->vector_cost_function = [](ind_t n, const T* i, const T* j){
        auto sin_theta_H = std::sin(brille::utils::hermitian_angle(n, i, j));
        return sin_theta_H*sin_theta_H;
      };
    }
  }
  void set_cost_info(const int scf, const int vcf, const ElementsCost& elcost){
    costs_ = elcost;
    this->set_cost_info(scf, vcf);
  }
  //
  ind_t size(void) const {return data_.size(0);}
  ind_t branches(void) const {
    ind_t nd = data_.ndim();
    if (nd < 2)
      throw std::runtime_error("data must have size (nPts, nBranches) or (nPts, nBranches, nElements), not (nPts,)");
    ind_t num = data_.size(1);
    if (2 == nd){
      ind_t el_per_branch = std::accumulate(elements_.begin(), elements_.end(), ind_t(0));
      // zero-elements is the special (initial) case
      // otherwise we can safely divide *AND* n_branches*el_per_branch ≡ num
      if (el_per_branch > 0) num /= el_per_branch;
    }
    return num;
  }
  const data_t<T,P>& data(void) const {return data_;}
  const ElementsType& elements(void) const {return elements_;}
  //
  template<class S>
  void interpolate_at(const std::vector<std::vector<ind_t>>&, const std::vector<ind_t>&, const std::vector<double>&, data_t<T,S>&, const ind_t, const bool) const;
  //
  template<class S>
  void interpolate_at(const std::vector<std::vector<ind_t>>&, const std::vector<std::pair<ind_t,double>>&, data_t<T,S>&, const ind_t, const bool) const;
  //
  template<class S, class R, class Q, class RotT>
  bool rotate_in_place(data_t<T,S>& x,
                       const LQVec<R>& q,
                       const RotT& rt,
                       const PointSymmetry& ps,
                       const std::vector<size_t>& r,
                       const std::vector<size_t>& invr,
                       const int nth) const {
    switch (rotlike_){
      case RotatesLike::Real: return this->rip_real(x, ps, r, invr, nth);
      case RotatesLike::Axial: return this->rip_axial(x, ps, r, invr, nth);
      case RotatesLike::Reciprocal: return this->rip_recip(x, ps, r, invr, nth);
      case RotatesLike::Gamma: return this->rip_gamma(x, q, rt, ps, r, invr, nth);
      default:
        throw std::runtime_error("Impossible RotatesLike value!");
    }
  }
  //
  RotatesLike rotateslike() const { return rotlike_; }
  RotatesLike rotateslike(const RotatesLike a) {
    rotlike_ = a;
    return rotlike_;
  }
  // Replace the data within this object.
  template<class S, class I>
  void replace_data(
      const data_t<T,S>& nd,
      const std::array<I,3>& ne,
      const RotatesLike rl = RotatesLike::Real)
  {
    data_ = nd;
    rotlike_ = rl;
    // convert the elements datatype as necessary
    for (size_t i=0; i<3u; ++i) elements_[i] = static_cast<ind_t>(ne[i]);
    if (ne[1]%3)
      throw std::logic_error("Vectors must have 3N elements per branch");
    if (ne[2]%9)
      throw std::logic_error("Matrices must have 9N elements per branch");
    // check the input for correctness
    ind_t total_elements = 1u;
    // scalar + eigenvector + vector + matrix elements
    ind_t known_elements = this->branch_span(ne);
    // no matter what, shape[0] should be the number of gridded points
    if (data_.ndim()>2){
      // if the number of dimensions of the array is greater than two,
      // the second element is the number of modes per point
      for (size_t i=2u; i<data_.ndim(); ++i) total_elements *= data_.size(i);
    } else {
      // shape is [n_points, n_branches*n_elements] or [n_points,], so there is only one mode
      total_elements = data_.ndim() > 1 ? data_.size(1) : 1u;
      // in case of n_branches*n_elements along the second dimesnion:
      if (total_elements % known_elements == 0){
        ind_t n_branches = total_elements/known_elements;
        total_elements /= n_branches;
        // or just total_elements = known_elements;
      }
    }
    if (0 == known_elements) elements_[0] = total_elements;
    if (known_elements && known_elements != total_elements){
      std::string msg;
      msg = "Inconsistent elements: " + std::to_string(known_elements) + " = ";
      msg += std::to_string(elements_[0]) + "+" + std::to_string(elements_[1]) + "+";
      msg += std::to_string(elements_[2]) + " ≠ ";
      msg += std::to_string(total_elements);
      throw std::runtime_error(msg);
    }
  }
  // Replace the data in this object without specifying the data shape or its elements
  // this variant is necessary since the template specialization above can not have a default value for the elements
  void replace_data(const data_t<T,P>& nd){
    return this->replace_data(nd, ElementsType({{0,0,0}}));
  }
  ind_t branch_span() const { return this->branch_span(elements_);}
  //
  std::string to_string() const {
    std::string str= "( ";
    for (auto s: data_.shape()) str += std::to_string(s) + " ";
    str += ") data";
    auto branches_ = this->branches();
    if (branches_){
      str += " with " + std::to_string(branches_) + " mode";
      if (branches_>1) str += "s";
    }
    auto neltypes = std::count_if(elements_.begin(), elements_.end(), [](ind_t a){return a>0;});
    if (neltypes){
      str += " of ";
      std::array<std::string,3> types{"scalar", "vector", "matrix"};
      for (size_t i=0; i<3u; ++i) if (elements_[i]) {
        str += std::to_string(elements_[i]) + " " + types[i];
        if (--neltypes>1) str += ", ";
        if (1==neltypes) str += " and ";
      }
      str += " element";
      if (this->branch_span()>1) str += "s";
    }
    return str;
  }
  template<typename S>
  void add_cost(const ind_t i0, const ind_t i1, std::vector<S>& cost, const bool arbitrary_phase_allowed) const {
    // can this be extended to account for the arbitrary phase too?
    S s_cost{0}, v_cost{0}, m_cost{0};
    ind_t nel2 = std::sqrt(elements_[2]);
    ind_t span = this->branch_span();
    ind_t b_ = this->branches();
    ind_t s_ = (data_.ndim() > 2) ? 1u : span;
    auto e_ = elements_;
    brille::shape_t x0{i0,0}, x1{i1,0};
    if (arbitrary_phase_allowed){ // if the vector_cost_function uses the Hermitian angle, e^iθ *never* matters.
      auto phased = std::unique_ptr<T[]>(new T[span]);
      for (ind_t i=0; i<b_; ++i){
        x0[1] = i*s_;
        const T* d0_i = data_.ptr(x0);
        for (ind_t j=0; j<b_; ++j){
          x1[1] = j*s_;
          const T* d1_j = data_.ptr(x1);
          T eith = brille::utils::antiphase(span, d0_i, d1_j);
          for (size_t s=0; s<span; ++s) phased[s] = eith*d1_j[s];
          if (e_[0]) s_cost = this->scalar_cost_function(e_[0], d0_i, phased.get());
          if (e_[1]) v_cost = this->vector_cost_function(e_[1], d0_i+e_[0], phased.get()+e_[0]);
          if (e_[2]) m_cost = brille::utils::frobenius_distance(nel2, d0_i+e_[0]+e_[1], phased.get()+e_[0]+e_[1]);
          cost[i*b_+j] += costs_[0]*s_cost + costs_[1]*v_cost + costs_[2]*m_cost;
        }
      }
    } else {
      for (ind_t i=0; i<b_; ++i){
        x0[1] = i*s_;
        const T* d0_i = data_.ptr(x0);
        for (ind_t j=0; j<b_; ++j){
          x1[1] = j*s_;
          const T* d1_j = data_.ptr(x1);
          if (e_[0]) s_cost = this->scalar_cost_function(e_[0], d0_i, d1_j);
          if (e_[1]) v_cost = this->vector_cost_function(e_[1], d0_i+e_[0], d1_j+e_[0]);
          if (e_[2]) m_cost = brille::utils::frobenius_distance(nel2, d0_i+e_[0]+e_[1], d1_j+e_[0]+e_[1]);
          cost[i*b_+j] += costs_[0]*s_cost + costs_[1]*v_cost + costs_[2]*m_cost;
        }
      }
    }
  }
  // template<typename I>
  // void permute_modes(const I idx, const std::vector<ind_t>& p){
  //   std::vector<ind_t> perm;
  //   std::copy(p.begin(), p.end(), std::back_inserter(perm));
  //   // perm is used as scratch space by Array::permute_modes
  //   data_.permute_modes(idx, perm);
  // }
  // template<typename I>
  // void inverse_permute_modes(const I idx, const std::vector<ind_t>& invp){
  //   std::vector<ind_t> invperm;
  //   std::copy(invp.begin(), invp.end(), std::back_inserter(invperm));
  //   // invperm is used as scratch space by Array::inverse_permute_modes
  //   data_.inverse_permute_modes(idx, invperm);
  // }
  template<typename I>
  bool any_equal_modes(const I idx) {
    size_t span = static_cast<size_t>(this->branch_span());
    // since we're probably only using this when the data is provided and
    // most eigenproblem solvers sort their output by eigenvalue magnitude it is
    // most-likely for mode i and mode i+1 to be equal.
    // ∴ search (i,j), (i+1,j+1), (i+2,j+2), ..., i ∈ (0,N], j ∈ (1,N]
    // for each j = i+1, i+2, i+3, ..., i+N-1
    ind_t b_ = this->branches();
    ind_t s_ = (data_.ndim() > 2) ? 1u : span;
    brille::shape_t xi{static_cast<ind_t>(idx),0}, xj{static_cast<ind_t>(idx),0};
    if (b_ > 1)
    for (ind_t offset=1; offset < b_; ++offset)
    for (ind_t i=0, j=offset; j < b_; ++i, ++j){
      xi[1] = i*s_;
      xj[1] = j*s_;
      if (brille::approx::vector(span, data_.ptr(xi), data_.ptr(xj))) return true;
    }
    // no matches
    return false;
  }
  size_t bytes_per_point() const {
    size_t n_elements = data_.numel()/data_.size(0);
    return n_elements * sizeof(T);
  }
private:
  template<typename I> ind_t branch_span(const std::array<I,3>& e) const {
    return static_cast<ind_t>(e[0])+static_cast<ind_t>(e[1])+static_cast<ind_t>(e[2]);
  }
  ElementsType count_scalars_vectors_matrices(void) const {
    ElementsType no{elements_[0], elements_[1]/3u, elements_[2]/9u};
    return no;
  }
  template<class S>
  bool rip_real(data_t<T,S>&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  template<class S>
  bool rip_recip(data_t<T,S>&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  template<class S>
  bool rip_axial(data_t<T,S>&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  template<class S, class R, class Q>
  bool rip_gamma_complex(data_t<T,S>&, const LQVec<R>&, const GammaTable&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  template<class P0, class R, class P1, class S=T>
  enable_if_t<is_complex<S>::value, bool>
  rip_gamma(data_t<T,P0>& x, const LQVec<R>& q, const GammaTable& gt, const PointSymmetry& ps, const std::vector<size_t>& r, const std::vector<size_t>& ir, const int nth) const{
    return rip_gamma_complex(x, q, gt, ps, r, ir, nth);
  }
  template<class P0, class R, class P1, class S=T>
  enable_if_t<!is_complex<S>::value, bool>
  rip_gamma(data_t<T,P0>&, const LQVec<R>&, const GammaTable&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const{
    throw std::runtime_error("RotatesLike == Gamma requires complex valued data!");
  }
};

template<class T, class P> template<class Q>
void InnerInterpolationData<T,P>::interpolate_at(
  const std::vector<std::vector<ind_t>>& permutations,
  const std::vector<ind_t>& indices,
  const std::vector<double>& weights,
  data_t<T,Q>& out,
  const ind_t to,
  const bool arbitrary_phase_allowed
) const {
  if (indices.size()==0 || weights.size()==0)
    throw std::logic_error("Interpolation requires input data!");
  ind_t span = this->branch_span();
  ind_t b_ = this->branches();
  ind_t s_ = (data_.ndim() > 2) ? 1u : span;
  verbose_update("Combining\n",data_.extract(indices).to_string(),"with weights ", weights);
  shape_t xidx{0,0}, oidx{to,0}, zidx{indices[0],0};
  if (arbitrary_phase_allowed){
    for (size_t x=0; x<indices.size(); ++x){
      xidx[0] = indices[x];
      for (ind_t b=0; b < b_; ++b){
        xidx[1] = (permutations[x][b])*s_;
        zidx[1] = oidx[1] = b*s_;
        const T* ptr0 = data_.ptr(zidx);
        const T* ptrX = data_.ptr(xidx);
        T eith = brille::utils::antiphase(span, ptr0, ptrX);
        T* ptrO = out.ptr(oidx);
        for (ind_t s=0; s<span; ++s) ptrO[s] += weights[x]*eith*ptrX[s];
      }
    }
  } else {
    for (size_t x=0; x<indices.size(); ++x){
      xidx[0] = indices[x];
      for (ind_t b=0; b < b_; ++b){
        xidx[1] = (permutations[x][b])*s_;
        oidx[1] = b*s_;
        const T* ptrX = data_.ptr(xidx);
        T* ptrO = out.ptr(oidx);
        for (ind_t s=0; s<span; ++s) ptrO[s] += weights[x]*ptrX[s];
      }
    }
  }
}

template<class T, class P> template<class Q>
void InnerInterpolationData<T,P>::interpolate_at(
  const std::vector<std::vector<ind_t>>& permutations,
  const std::vector<std::pair<ind_t,double>>& indices_weights,
  data_t<T,Q>& out,
  const ind_t to,
  const bool arbitrary_phase_allowed
) const {
  if (indices_weights.size()==0)
    throw std::logic_error("Interpolation requires input data!");
  ind_t span = this->branch_span();
  std::vector<int> dummy;
  ind_t b_ = this->branches();
  ind_t s_ = (data_.ndim() > 2) ? 1u : span;
  if (arbitrary_phase_allowed){
    std::transform(permutations.begin(), permutations.end(), indices_weights.begin(), std::back_inserter(dummy),
    [&](const std::vector<ind_t>& perm, const std::pair<ind_t,double>& iw){
      for (ind_t b=0; b<b_; ++b)
      {
        shape_t oidx{to, b*s_}, xidx{iw.first, perm[b]*s_}, zidx{indices_weights[0].first, b*s_};
        const T *ptrX = data_.ptr(xidx);
        const T *ptr0 = data_.ptr(zidx);
        T eith = brille::utils::antiphase(span, ptr0, ptrX);
        T *ptrO = out.ptr(oidx);
        for (ind_t s=0; s<span; ++s) ptrO[s] += iw.second*eith*ptrX[s];
      }
      return 1;
    });
  } else {
    std::transform(permutations.begin(), permutations.end(), indices_weights.begin(), std::back_inserter(dummy),
    [&](const std::vector<ind_t>& perm, const std::pair<ind_t,double>& iw){
      for (ind_t b=0; b<b_; ++b)
      {
        shape_t oidx{to, b*s_}, xidx{iw.first, perm[b]*s_};
        T *ptrO = out.ptr(oidx);
        const T *ptrX = data_.ptr(xidx);
        for (ind_t s=0; s<span; ++s) ptrO[s] += iw.second*ptrX[s];
      }
      return 1;
    });
  }
}

template<class T, class P> template<class Q>
bool InnerInterpolationData<T,P>::rip_recip(
  data_t<T,Q>& x, const PointSymmetry& ptsym, const std::vector<size_t>& r, const std::vector<size_t>& invR, const int nthreads
) const {
  profile_update("Start InnerInterpolationData::rip_recip method");
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  ElementsType no = this->count_scalars_vectors_matrices();
  if (!std::any_of(no.begin()+1, no.end(), [](ind_t n){return n>0;}))
    return false;
  T tmp_v[3], tmp_m[9];
  const ind_t b_{this->branches()};
  const ind_t s_{(data_.ndim() > 2) ? 1u : this->branch_span()};
  shape_t xidx{0,0};
  if (data_.ndim() > 2) xidx.resize(3);
  std::array<int,9> ident = {1,0,0, 0,1,0, 0,0,1};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
// otherwise gcc complains that const b_ and s_ are *already* shared
#pragma omp parallel for default(none) shared(x,ptsym,r,invR) private(tmp_v,tmp_m) firstprivate(xidx,ident,no,xsize) schedule(static)
#else
#pragma omp parallel for default(none) shared(b_,s_,x,ptsym,r,invR) private(tmp_v,tmp_m) firstprivate(xidx,ident,no,xsize) schedule(static)
#endif
  for (long long si=0; si<xsize; ++si){
    T* xptr;
    xidx[0] = brille::utils::s2u<ind_t, long long>(si);
    if (!brille::approx::matrix(3, ident.data(), ptsym.get(r[xidx[0]]).data()))
    for (ind_t b=0; b<b_; ++b){
      xidx.back() = 0u; // either xidx[1] or xidx[2]
      xidx[1] = b*s_;  // if ndim()==2 we *must* include the branch span offset
      // scalar elements do not need to be rotated, so skip them
      xidx.back() += no[0]; // either xidx[1] or xidx[2]
      // we chose R such that Q = Rᵀq + τ.
      for (ind_t v=0; v<no[1]; ++v){
        xptr = x.ptr(xidx);
        brille::utils::mul_mat_vec(tmp_v, 3u, transpose(ptsym.get(r[xidx[0]])).data(), xptr);
        for (int j=0; j<3; ++j) xptr[j] = tmp_v[j];
        xidx.back() += 3u; // shift 3 for each vector
      }
      for (ind_t m=0; m<no[2]; ++m){
        xptr = x.ptr(xidx);
        // Calculate R*M*R⁻¹ in two steps
        // first calculate M*R⁻¹, storing in tmp_m
        brille::utils::mul_mat_mat(tmp_m, 3u, xptr, transpose(ptsym.get(invR[xidx[0]])).data());
        // next calculate R*tmp_m, storing back in the x array
        brille::utils::mul_mat_mat(xptr, 3u, transpose(ptsym.get(r[xidx[0]])).data(), tmp_m);
        xidx.back() += 9u; // shift 9 for each matrix
      }
    }
  }
  profile_update("  End InnerInterpolationData::rip_recip method");
  return true;
}

template<class T, class P> template<class Q>
bool InnerInterpolationData<T,P>::rip_real(
  data_t<T,Q>& x, const PointSymmetry& ptsym, const std::vector<size_t>& r, const std::vector<size_t>& invR, const int nthreads
) const {
  profile_update("Start InnerInterpolationData::rip_real method");
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  ElementsType no = this->count_scalars_vectors_matrices();
  if (!std::any_of(no.begin()+1, no.end(), [](ind_t n){return n>0;}))
    return false;
  T tmp_v[3], tmp_m[9];
  const ind_t b_{this->branches()};
  const ind_t s_{(data_.ndim() > 2) ? 1u : this->branch_span()};
  shape_t xidx{0,0};
  if (data_.ndim() > 2) xidx.resize(3);
  std::array<int,9> ident = {1,0,0, 0,1,0, 0,0,1};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
#pragma omp parallel for default(none) shared(x,ptsym,r,invR) private(tmp_v,tmp_m) firstprivate(xidx,ident,no,xsize) schedule(static)
#else
#pragma omp parallel for default(none) shared(b_,s_,x,ptsym,r,invR) private(tmp_v,tmp_m) firstprivate(xidx,ident,no,xsize) schedule(static)
#endif
  for (long long si=0; si<xsize; ++si){
    xidx[0] = brille::utils::s2u<ind_t, long long>(si);
    T* xptr;
    if (!brille::approx::matrix(3, ident.data(), ptsym.get(r[xidx[0]]).data()))
    for (ind_t b=0; b<b_; ++b){
      xidx.back() = 0;
      xidx[1] = b*s_;
      // scalar elements do not need to be rotated, so skip them
      xidx.back() += no[0];
      // rotate real vectors: since Q = Rᵀq + τ → Rv
      for (ind_t v=0; v<no[1]; ++v){
        xptr = x.ptr(xidx);
        brille::utils::mul_mat_vec(tmp_v, 3u, ptsym.get(r[xidx[0]]).data(), xptr);
        for (int j=0; j<3; ++j) xptr[j] = tmp_v[j];
        xidx.back() += 3u;
      }
      for (ind_t m=0; m<no[2]; ++m){
        xptr = x.ptr(xidx);
        // Calculate R⁻¹*M*R in two steps
        // first calculate M*R, storing in tmp_m
        brille::utils::mul_mat_mat(tmp_m, 3u, xptr, ptsym.get(invR[xidx[0]]).data());
        // next calculate R⁻¹*tmp_m, storing back in the x array
        brille::utils::mul_mat_mat(xptr, 3u, ptsym.get(r[xidx[0]]).data(), tmp_m);
        xidx.back() += 9u;
      }
    }
  }
  profile_update("  End InnerInterpolationData::rip_real method");
  return true;
}

template<class T, class P> template<class Q>
bool InnerInterpolationData<T,P>::rip_axial(
  data_t<T,Q>& x, const PointSymmetry& ptsym, const std::vector<size_t>& r, const std::vector<size_t>& invR, const int nthreads
) const {
  profile_update("Start InnerInterpolationData::rip_axial method");
  omp_set_num_threads((nthreads > 0) ? nthreads : omp_get_max_threads());
  ElementsType no = this->count_scalars_vectors_matrices();
  if (!std::any_of(no.begin() + 1, no.end(), [](ind_t n) {return n > 0; }))
      return false;
  T tmp_v[3], tmp_m[9];
  std::vector<T> detR;
  if (no[1])
      std::transform(r.begin(), r.end(), std::back_inserter(detR),
          [ptsym](const size_t z) { return static_cast<T>(brille::utils::matrix_determinant(ptsym.get(z).data())); });
  const ind_t b_{this->branches()};
  const ind_t s_{(data_.ndim() > 2) ? 1u : this->branch_span()};
  shape_t xidx{0,0};
  if (data_.ndim() > 2) xidx.resize(3);
  std::array<int,9> ident = {1,0,0, 0,1,0, 0,0,1};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
#pragma omp parallel for default(none) shared(x,ptsym,r,invR,detR) private(tmp_v,tmp_m) firstprivate(xidx,ident,no,xsize) schedule(static)
#else
#pragma omp parallel for default(none) shared(b_,s_,x,ptsym,r,invR,detR) private(tmp_v,tmp_m) firstprivate(xidx,ident,no,xsize) schedule(static)
#endif
  for (long long si = 0; si < xsize; ++si) {
    xidx[0]= brille::utils::s2u<ind_t, long long>(si);
    T* xptr;
    if (!brille::approx::matrix(3, ident.data(), ptsym.get(r[xidx[0]]).data()))
    for (ind_t b = 0; b < b_; ++b) {
      xidx.back() = 0u;
      xidx[1] = b*s_;
      // scalar elements do not need to be rotated, so skip them
      xidx.back() += no[0];
      // rotate real vectors: since Q = Rᵀq + τ → R⁻¹*v
      for (ind_t v = 0; v < no[1]; ++v) {
        xptr = x.ptr(xidx);
        brille::utils::mul_mat_vec(tmp_v, 3u, ptsym.get(invR[xidx[0]]).data(), xptr);
        for (int j = 0; j < 3; ++j) xptr[j] = detR[xidx[0]]*tmp_v[j];
        xidx.back() += 3u;
      }
      for (ind_t m = 0; m < no[2]; ++m) {
        xptr = x.ptr(xidx);
        // Calculate R⁻¹*M*R in two steps
        // first calculate M*R, storing in tmp_m
        brille::utils::mul_mat_mat(tmp_m, 3u, xptr, ptsym.get(r[xidx[0]]).data());
        // next calculate R⁻¹*tmp_m, storing back in the x array
        brille::utils::mul_mat_mat(xptr, 3u, ptsym.get(invR[xidx[0]]).data(), tmp_m);
        xidx.back() += 9u;
      }
    }
  }
  profile_update("  End InnerInterpolationData::rip_axial method");
  return true;
}

template<class T, class P, class R, class Q,
class S = typename std::common_type<T,R>::type,
class I = typename brille::ind_t
>
static std::complex<S>
e_iqd(const LQVec<T>& q, const I i, const bArray<R>& d, const size_t j)
{
  // const double pi = 3.14159265358979323846;
  S dotqd{0};
  I jj = static_cast<I>(j);
  for (I k=0; k<3u; ++k)
    dotqd += q.val(i,k)*d.val(jj,k);
  std::complex<S> i_2pi_x(0, 2*brille::pi*dotqd);
  return std::exp(i_2pi_x);
}

template<class T, class P>
template<class P0, class R, class P1>
bool InnerInterpolationData<T,P>::rip_gamma_complex(
  data_t<T,P0>& x, const LQVec<R>& q, const GammaTable& pgt,
  const PointSymmetry& ptsym, const std::vector<size_t>& ridx, const std::vector<size_t>& invRidx,
  const int nthreads
) const {
  profile_update("Start InnerInterpolationData::rip_gamma_complex method");
  // construct a lambda to calculate the phase for qᵢ, and [R⁻¹xₖ - xᵥ]
  auto e_iqd_gt = [q,pgt](ind_t i, ind_t k, size_t r){
    return e_iqd(q, i, pgt.vectors(), pgt.vector_index(k,r));
  };
  /* For speed purposes the GammaTable vectors do not contain lattice information
  when accessed by .vector, as below. There is another accessor .ldvector which
  adds the lattice, but this significantly slows down the dot product due to
  equivalent/inverse lattice checks reconstructing the Direct/Reciprocal objects
  multiple times *for every q point*(!).
  To avoid that check we should do it once, now. Then the dot product of
  q = (h,k,l) and d = (a,b,c) is q⋅d = 2π (h*a + k*b + l*c).
  */
  if (! pgt.lattice().isstar(q.get_lattice()))
    throw std::runtime_error("The q points and GammaTable must be in mutually reciprocal lattices");
  verbose_update("InnerInterpolationData::rip_gamma_complex called with ",nthreads," threads");
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  ElementsType no = this->count_scalars_vectors_matrices();
  if (!std::any_of(no.begin()+1, no.end(), [](ind_t n){return n>0;}))
    return false;
  // for Γ transformations of tensors, there *must* be N×N in total:
  ind_t Nmat = static_cast<ind_t>(std::sqrt(no[2]))/3;
  if (no[2] != 9*Nmat*Nmat){
    std::cout << "Atomic displacement Gamma transformation requires NxN 3x3 tensors!" << std::endl;
    return false;
  }
  const ind_t b_{this->branches()};
  const ind_t s_{(data_.ndim() > 2) ? 1u : this->branch_span()};
  shape_t xidx{0,0};
  if (data_.ndim() > 2) xidx.resize(3);
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
#pragma omp parallel for default(none) shared(x,q,pgt,ptsym,ridx,invRidx,e_iqd_gt) firstprivate(xidx,no,Nmat,xsize) schedule(static)
#else
#pragma omp parallel for default(none) shared(b_,s_,x,q,pgt,ptsym,ridx,invRidx,e_iqd_gt) firstprivate(xidx,no,Nmat,xsize) schedule(static)
#endif
  for (long long si=0; si<xsize; ++si){
    xidx[0] = brille::utils::s2u<ind_t, long long>(si);
    ind_t offset{0};
    T* xptr;
    for (ind_t b=0; b<b_; ++b){
      xidx.back() = 0u;
      xidx[1] = b*s_;
      // scalar elements do not need to be rotated, so skip them
      xidx.back() += no[0];
      shape_t yidx{xidx};
      // the *full 3N* eigenvector for mode j transforms as
      //  ϵʲ(Rq)  =  Γ(q;R) eʲ(q)
      // where Γ(q;R) is a 3N×3N matrix with elements given by
      // the submatrices Γₖᵥᵅᵝ(q;R) = Rᵅᵝ δ(v,F₀⁻¹(k,R)) exp{i q⋅[R⁻¹xₖ - xᵥ]}
      // -- I *think* that for a given k and R ∑ᵥ δ(v,F₀⁻¹(k,R)) ≡ 1
      //    that is, all equivalent-atom mappings are singular. THIS SHOULD BE VERIFIED!
      // The GammaTable contains predetermined F₀(k,R) and R⁻¹xₖ - xᵥ
      // to make this calculation more straightforward.
      //
      // Rotate, permute, and apply the phase factor simultaneously into a temporary array
      if (no[1]>0){
        T tmp_v[3];
        std::vector<T> tmpvecs(no[1]*3u);
        for (ind_t k=0; k<no[1]; ++k){
          // is q already expressed in the right lattice? hopefully!
          // if the vector rotates by Γ is *must* be complex, so T is a complex type
          // use its constructor to make i q⋅[R⁻¹xₖ - xᵥ] and calculate the (k,R) phase
          T phase = e_iqd_gt(xidx[0], k, invRidx[xidx[0]]);
          xptr = x.ptr(xidx);
          brille::utils::mul_mat_vec(tmp_v, 3u, ptsym.get(invRidx[xidx[0]]).data(), xptr);
          auto v0idx = 3u*pgt.F0(k, invRidx[xidx[0]]);
          for (int j=0; j<3; ++j) tmpvecs[v0idx+j] = phase*tmp_v[j];
          xidx.back() += 3u;
        }
        xptr = x.ptr(yidx);
        for (ind_t j=0; j<no[1]*3u; ++j) xptr[j] = tmpvecs[j];
        offset += no[1]*3u;
      }

      if (no[2]>0){
        ind_t xback = xidx.back();
        T tmp_m[9], tmp_m2[9];
        std::vector<T> tmpmats(no[2]*9u);
        for (ind_t n=0; n<Nmat; ++n){
          T rph = e_iqd_gt(xidx[0], n, ridx[xidx[0]]);
          ind_t v = static_cast<ind_t>(pgt.F0(n, ridx[xidx[0]]));
          for (ind_t m=0; m<Nmat; ++m){
            T invRph = e_iqd_gt(xidx[0], m, invRidx[xidx[0]]);
            ind_t k = static_cast<ind_t>(pgt.F0(m, invRidx[xidx[0]]));
            // Calculate R⁻¹*M*R in two steps
            // first calculate M*R, storing in tmp_m
            xidx.back() = xback + (n*Nmat+m)*9u;
            brille::utils::mul_mat_mat(tmp_m, 3u, x.ptr(xidx), ptsym.get(ridx[xidx[0]]).data());
            // next calculate R⁻¹*tmp_m, storing in the temporary all matrix array
            brille::utils::mul_mat_mat(tmp_m2, 3u, ptsym.get(invRidx[xidx[0]]).data(), tmp_m);
            // include the R R⁻¹ phase factor
            for (int j=0; j<9; ++j) tmpmats[(v*Nmat+k)*9u+j] = rph*invRph*tmp_m2[j];
          }
        }
        xidx.back() = xback; // move back to the start of all matrices
        xptr = x.ptr(xidx);
        for (ind_t j=0; j<no[2]*9u; ++j) xptr[j] = tmpmats[j];
      }
    }
  }
  profile_update("  End InnerInterpolationData::rip_gamma_complex method");
  return true;
}



#endif
