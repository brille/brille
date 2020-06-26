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
#include "arrayvector.hpp"
#include "latvec.hpp"
#include "phonon.hpp"
#include "permutation.hpp"
#include "permutation_table.hpp"

#ifndef _INTERPOLATION_DATA_H_
#define _INTERPOLATION_DATA_H_

typedef unsigned long element_t;
typedef std::vector<size_t> ShapeType;
typedef std::array<element_t,3> ElementsType;
typedef std::array<double,3> ElementsCost;

using int_t = PermutationTable::int_t;


//template<class T> using CostFunction = std::function<typename CostTraits<T>::type(element_t, T*, T*)>;
template<class T> using CostFunction = std::function<double(element_t, T*, T*)>;

template<class T> struct is_complex {enum{value = false};};
template<class T> struct is_complex<std::complex<T>> {enum {value=true};};
template<bool C, typename T> using enable_if_t = typename std::enable_if<C,T>::type;

enum class RotatesLike {
  Real, Reciprocal, Axial, Gamma
};

template<class T> class InnerInterpolationData{
  ArrayVector<T> data_;   //!< The stored eigenvalue-like ArrayVector indexed like the holding-Object's vertices
  ShapeType shape_;       //!< A std::vector to indicate a possible higher-dimensional shape of each `data_` array
  ElementsType elements_; //!< The number of scalars, vector elements, and matrix elements per `data_` array
  element_t branches_;    //!< The number of branches contained per data array
  RotatesLike rotlike_;   //!< How the elements of `data_` rotate
  ElementsCost costs_;    //!< The cost assigned to each type for equivalent mode assignment
  CostFunction<T> scalar_cost_function;
  CostFunction<T> vector_cost_function;
public:
  explicit InnerInterpolationData(size_t scf_type=0, size_t vcf_type=0):
  data_({0,0}), shape_({0,0}), elements_({{0,0,0}}), branches_(0),
  rotlike_{RotatesLike::Real}, costs_({{1,1,1}}) {
    this->set_cost_info(scf_type, vcf_type);
  };
  InnerInterpolationData(CostFunction<T> scf, CostFunction<T> vcf):
    data_({0,0}), shape_({0,0}), elements_({{0,0,0}}), branches_(0),
    rotlike_{RotatesLike::Real}, costs_({{1,1,1}}), scalar_cost_function(scf), vector_cost_function(vcf){};
  //
  void setup_fake(const size_t sz, const element_t br){
    data_.refresh(br, sz);
    shape_ = {sz, static_cast<size_t>(br)};
    elements_ = {1u,0u,0u};
    branches_ = br;
  }
  //
  void set_cost_info(const int scf, const int vcf){
    switch (scf){
      default:
      this->scalar_cost_function = [](element_t n, T* i, T* j){
        double s{0};
        for (element_t z=0; z<n; ++z) s += magnitude(i[z]-j[z]);
        return s;
      };
    }
    switch (vcf){
      case 1:
      debug_update("selecting vector_distance");
      this->vector_cost_function = [](element_t n, T* i, T* j){return vector_distance(n, i, j);};
      break;
      case 2:
      debug_update("selecting 1-vector_product");
      this->vector_cost_function = [](element_t n, T* i, T* j){return 1-vector_product(n, i, j);};
      break;
      case 3:
      debug_update("selecting vector_angle");
      this->vector_cost_function = [](element_t n, T* i, T* j){return vector_angle(n, i, j);};
      break;
      case 4:
      debug_update("selecting hermitian_angle");
      this->vector_cost_function = [](element_t n, T* i, T* j){ return hermitian_angle(n,i,j);};
      break;
      default:
      debug_update("selecting sin**2(hermitian_angle)");
      // this->vector_cost_function = [](element_t n, T* i, T* j){return std::abs(std::sin(hermitian_angle(n, i, j)));};
      this->vector_cost_function = [](element_t n, T* i, T* j){
        auto sin_theta_H = std::sin(hermitian_angle(n, i, j));
        return sin_theta_H*sin_theta_H;
      };
    }
  }
  void set_cost_info(const int scf, const int vcf, const ElementsCost& elcost){
    costs_ = elcost;
    this->set_cost_info(scf, vcf);
  }
  //
  size_t size(void) const {return data_.size();}
  size_t numel(void) const {return data_.numel();}
  const ArrayVector<T>& data(void) const {return data_;}
  const ShapeType& shape(void) const {return shape_;}
  const ElementsType& elements(void) const {return elements_;}
  element_t branches(void) const {return branches_;}
  //
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  void interpolate_at(const std::vector<std::vector<int_t>>&, const std::vector<I>&, const std::vector<double>&, ArrayVector<T>&, const size_t, const bool) const;
  //
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  void interpolate_at(const std::vector<std::vector<int_t>>&, const std::vector<std::pair<I,double>>&, ArrayVector<T>&, const size_t, const bool) const;
  //
  template<class R, class RotT>
  bool rotate_in_place(ArrayVector<T>& x,
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
  template<typename I> void replace_data(
      const ArrayVector<T>& nd,
      const ShapeType& ns,
      const std::array<I,3>& ne,
      const RotatesLike rl = RotatesLike::Real
    ){
    data_ = nd;
    shape_ = ns;
    rotlike_ = rl;
    // convert the elements datatype as necessary
    for (size_t i=0; i<3u; ++i) elements_[i] = static_cast<element_t>(ne[i]);
    if (ne[1]%3)
      throw std::logic_error("Vectors must have 3N elements per branch");
    if (ne[2]%9)
      throw std::logic_error("Matrices must have 9N elements per branch");
    // check the input for correctness
    element_t total_elements = 1u;
    // scalar + eigenvector + vector + matrix*matrix elements
    element_t known_elements = this->branch_span(ne);
    // no matter what, shape[0] should be the number of gridded points
    if (ns.size()>2){
      // if the number of dimensions of the shape array is greater than two,
      // the second element is the number of modes per point                    */
      branches_ = static_cast<element_t>(ns[1]);
      for (size_t i=2u; i<ns.size(); ++i) total_elements *= static_cast<element_t>(ns[i]);
    } else {
      // shape is [n_points, n_elements] or [n_points,], so there is only one mode
      branches_ = 1u;
      total_elements = static_cast<element_t>( ns.size() > 1 ? ns[1] : 1u );
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
  // Replace the data in this object without specifying the data shape
  template<typename I> void replace_data(const ArrayVector<T>& nd, const std::array<I,3>& ne){
    ShapeType ns{nd.size(), nd.numel()};
    return this->replace_data(nd, ns, ne);
  }
  // Replace the data in this object without specifying the data shape or its elements
  // this variant is necessary since the template specialization above can not have a default value for the elements
  void replace_data(const ArrayVector<T>& nd){
    return this->replace_data(nd, ElementsType({{0,0,0}}));
  }
  element_t branch_span() const { return this->branch_span(elements_);}
  //
  std::string to_string() const {
    std::string str= "( ";
    for (auto s: shape_) str += std::to_string(s) + " ";
    str += ") data";
    if (branches_){
      str += " with " + std::to_string(branches_) + " mode";
      if (branches_>1) str += "s";
    }
    auto neltypes = std::count_if(elements_.begin(), elements_.end(), [](element_t a){return a>0;});
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
  template<typename I, typename S>
  void add_cost(const I i0, const I i1, std::vector<S>& cost, const bool arbitrary_phase_allowed) const {
    // can this be extended to account for the arbitrary phase too?
    S s_cost{0}, v_cost{0}, m_cost{0};
    T *d0_i, *d1_j;
    element_t nel2 = std::sqrt(elements_[2]);
    element_t span = this->branch_span();
    if (arbitrary_phase_allowed){ // if the vector_cost_function uses the Hermitian angle, e^iθ *never* matters.
      auto phased = std::unique_ptr<T[]>(new T[span]);
      for (size_t i=0; i<branches_; ++i) for (size_t j=0; j<branches_; ++j){
        d0_i = data_.data(i0,i*span);
        d1_j = data_.data(i1,j*span);
        T eith = antiphase(span, d0_i, d1_j);
        for (size_t s=0; s<span; ++s) phased[s] = eith*d1_j[s];
        if (elements_[0]) s_cost = this->scalar_cost_function(elements_[0], d0_i, phased.get());
        if (elements_[1]) v_cost = this->vector_cost_function(elements_[1], d0_i+elements_[0], phased.get()+elements_[0]);
        if (elements_[2]) m_cost = frobenius_distance(nel2, d0_i+elements_[0]+elements_[1], phased.get()+elements_[0]+elements_[1]);
        cost[i*branches_+j] += costs_[0]*s_cost + costs_[1]*v_cost + costs_[2]*m_cost;
      }
    } else {
      for (size_t i=0; i<branches_; ++i) for (size_t j=0; j<branches_; ++j){
        d0_i = data_.data(i0,i*span);
        d1_j = data_.data(i1,j*span);
        if (elements_[0]) s_cost = this->scalar_cost_function(elements_[0], d0_i, d1_j);
        if (elements_[1]) v_cost = this->vector_cost_function(elements_[1], d0_i+elements_[0], d1_j+elements_[0]);
        if (elements_[2]) m_cost = frobenius_distance(nel2, d0_i+elements_[0]+elements_[1], d1_j+elements_[0]+elements_[1]);
        cost[i*branches_+j] += costs_[0]*s_cost + costs_[1]*v_cost + costs_[2]*m_cost;
      }
    }
    // info_update_if(arbitrary_phase_allowed, i0,'-',i1,'\n',cost);
  }
  template<typename I>
  void permute_modes(const I idx, const std::vector<int_t>& p){
    std::vector<int_t> perm;
    std::copy(p.begin(), p.end(), std::back_inserter(perm));
    // perm is used as scratch space by ArrayVector::permute_modes
    data_.permute_modes(idx, perm);
  }
  template<typename I>
  void inverse_permute_modes(const I idx, const std::vector<int_t>& invp){
    std::vector<int_t> invperm;
    std::copy(invp.begin(), invp.end(), std::back_inserter(invperm));
    // invperm is used as scratch space by ArrayVector::inverse_permute_modes
    data_.inverse_permute_modes(idx, invperm);
  }
  template<typename I>
  bool any_equal_modes(const I idx) const {
    size_t span = static_cast<size_t>(this->branch_span());
    // since we're probably only using this when the data is provided and
    // most eigenproblem solvers sort their output by eigenvalue magnitude it is
    // most-likely for mode i and mode i+1 to be equal.
    // ∴ search (i,j), (i+1,j+1), (i+2,j+2), ..., i ∈ (0,N], j ∈ (1,N]
    // for each j = i+1, i+2, i+3, ..., i+N-1
    for (element_t offset=1; offset < branches_; ++offset)
    for (element_t i=0, j=offset; j < branches_; ++i, ++j)
    if (approx_vector(span, data_.data(idx, i*span), data_.data(idx, j*span)))
      return true;
    // no matches
    return false;
  }
private:
  template<typename I> element_t branch_span(const std::array<I,3>& e) const {
    return static_cast<element_t>(e[0])+static_cast<element_t>(e[1])+static_cast<element_t>(e[2]);
  }
  ElementsType count_scalars_vectors_matrices(void) const {
    ElementsType no{elements_[0], elements_[1]/3u, elements_[2]/9u};
    return no;
  }
  bool rip_real(ArrayVector<T>&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  bool rip_recip(ArrayVector<T>&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  bool rip_axial(ArrayVector<T>&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  template<class R>
  bool rip_gamma_complex(ArrayVector<T>&, const LQVec<R>&, const GammaTable&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  template<class R, class S=T>
  enable_if_t<is_complex<S>::value, bool>
  rip_gamma(ArrayVector<T>& x, const LQVec<R>& q, const GammaTable& gt, const PointSymmetry& ps, const std::vector<size_t>& r, const std::vector<size_t>& ir, const int nth) const{
    return rip_gamma_complex(x, q, gt, ps, r, ir, nth);
  }
  template<class R, class S=T>
  enable_if_t<!is_complex<S>::value, bool>
  rip_gamma(ArrayVector<T>&, const LQVec<R>&, const GammaTable&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const{
    throw std::runtime_error("RotatesLike == Gamma requires complex valued data!");
  }
};

template<typename T> template<typename I, typename>
void InnerInterpolationData<T>::interpolate_at(
  const std::vector<std::vector<int_t>>& permutations,
  const std::vector<I>& indices,
  const std::vector<double>& weights,
  ArrayVector<T>& out,
  const size_t to,
  const bool arbitrary_phase_allowed
) const {
  if (indices.size()==0 || weights.size()==0)
    throw std::logic_error("Interpolation requires input data!");
  T *out_to = out.data(to), *ptr0 = data_.data(indices[0]);
  element_t span = this->branch_span();
  verbose_update("Combining\n",data_.extract(indices).to_string(),"with weights ", weights);
  if (arbitrary_phase_allowed){
    for (size_t x=0; x<indices.size(); ++x){
      T *ptrX = data_.data(indices[x]);
      for (element_t b=0; b < branches_; ++b){
        element_t p = static_cast<element_t>(permutations[x][b]);
        T eith = antiphase(span, ptr0+b*span, ptrX+p*span);
        for (size_t s=0; s<span; ++s) out_to[b*span+s] += weights[x]*eith*ptrX[p*span + s];
      }
    }
  } else {
    for (size_t x=0; x<indices.size(); ++x){
      T *ptrX = data_.data(indices[x]);
      for (element_t b=0; b < branches_; ++b){
        element_t p = static_cast<element_t>(permutations[x][b]);
        for (size_t s=0; s<span; ++s) out_to[b*span+s] += weights[x]*ptrX[p*span+s];
      }
    }
  }
}

template<typename T> template<typename I, typename>
void InnerInterpolationData<T>::interpolate_at(
  const std::vector<std::vector<int_t>>& permutations,
  const std::vector<std::pair<I,double>>& indices_weights,
  ArrayVector<T>& out,
  const size_t to,
  const bool arbitrary_phase_allowed
) const {
  if (indices_weights.size()==0)
    throw std::logic_error("Interpolation requires input data!");
  T *out_to = out.data(to), *ptr0 = data_.data(indices_weights[0].first);
  element_t span = this->branch_span();
  std::vector<int> dummy;
  if (arbitrary_phase_allowed){
    std::transform(permutations.begin(), permutations.end(), indices_weights.begin(), std::back_inserter(dummy),
    [&](const std::vector<int_t>& perm, const std::pair<I,double>& iw){
      T *ptrX = data_.data(iw.first);
      for (element_t b=0; b<branches_; ++b){
        element_t p = static_cast<element_t>(perm[b]);
        T eith = antiphase(span, ptr0+b*span, ptrX+p*span);
        for (size_t s=0; s<span; ++s) out_to[b*span+s] += iw.second*eith*ptrX[p*span+s];
      }
      return 1;
    });
  } else {
    std::transform(permutations.begin(), permutations.end(), indices_weights.begin(), std::back_inserter(dummy),
    [&](const std::vector<int_t>& perm, const std::pair<I,double>& iw){
      T *ptrX = data_.data(iw.first);
      for (element_t b=0; b<branches_; ++b){
        element_t p = static_cast<element_t>(perm[b]);
        for (size_t s=0; s<span; ++s) out_to[b*span+s] += iw.second*ptrX[p*span+s];
      }
      return 1;
    });
  }
}

template<typename T>
bool InnerInterpolationData<T>::rip_recip(
  ArrayVector<T>& x, const PointSymmetry& ptsym, const std::vector<size_t>& r, const std::vector<size_t>& invR, const int nthreads
) const {
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  ElementsType no = this->count_scalars_vectors_matrices();
  if (!std::any_of(no.begin()+1, no.end(), [](element_t n){return n>0;}))
    return false;
  T tmp_v[3], tmp_m[9];
  element_t offset, sp = this->branch_span();
  std::array<int,9> ident = {1,0,0, 0,1,0, 0,0,1};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = unsigned_to_signed<long long, size_t>(x.size());
#pragma omp parallel for default(none) shared(x,ptsym,r,invR) private(offset, tmp_v, tmp_m) firstprivate(ident,no, sp, xsize) schedule(static)
  for (long long si=0; si<xsize; ++si){
    size_t i = signed_to_unsigned<size_t, long long>(si);
    if (!approx_matrix(3, ident.data(), ptsym.get(r[i]).data()))
    for (element_t b=0; b<branches_; ++b){
      // scalar elements do not need to be rotated, so skip them
      offset = b*sp + no[0];
      // we chose R such that Q = Rᵀq + τ.
      for (element_t v=0; v<no[1]; ++v){
        mul_mat_vec(tmp_v, 3u, transpose(ptsym.get(r[i])).data(), x.data(i, offset+v*3u));
        for (int j=0; j<3; ++j) x.insert(tmp_v[j], i, offset+v*3u+j);
      }
      offset += no[1]*3u;
      for (element_t m=0; m<no[2]; ++m){
        // Calculate R*M*R⁻¹ in two steps
        // first calculate M*R⁻¹, storing in tmp_m
        mul_mat_mat(tmp_m, 3u, x.data(i, offset+m*9u), transpose(ptsym.get(invR[i])).data());
        // next calculate R*tmp_m, storing back in the x array
        mul_mat_mat(x.data(i, offset+m*9u), 3u, transpose(ptsym.get(r[i])).data(), tmp_m);
      }
    }
  }
  return true;
}

template<typename T>
bool InnerInterpolationData<T>::rip_real(
  ArrayVector<T>& x, const PointSymmetry& ptsym, const std::vector<size_t>& r, const std::vector<size_t>& invR, const int nthreads
) const {
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  ElementsType no = this->count_scalars_vectors_matrices();
  if (!std::any_of(no.begin()+1, no.end(), [](element_t n){return n>0;}))
    return false;
  T tmp_v[3], tmp_m[9];
  element_t offset, sp = this->branch_span();
  std::array<int,9> ident = {1,0,0, 0,1,0, 0,0,1};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = unsigned_to_signed<long long, size_t>(x.size());
#pragma omp parallel for default(none) shared(x,ptsym,r,invR) private(offset, tmp_v, tmp_m) firstprivate(ident, no, sp, xsize) schedule(static)
  for (long long si=0; si<xsize; ++si){
    size_t i = signed_to_unsigned<size_t, long long>(si);
    if (!approx_matrix(3, ident.data(), ptsym.get(r[i]).data()))
    for (element_t b=0; b<branches_; ++b){
      // scalar elements do not need to be rotated, so skip them
      offset = b*sp + no[0];
      // rotate real vectors: since Q = Rᵀq + τ → Rv
      for (element_t v=0; v<no[1]; ++v){
        mul_mat_vec(tmp_v, 3u, ptsym.get(r[i]).data(), x.data(i, offset+v*3u));
        for (int j=0; j<3; ++j) x.insert(tmp_v[j], i, offset+v*3u+j);
      }
      offset += no[1]*3u;
      for (element_t m=0; m<no[2]; ++m){
        // Calculate R⁻¹*M*R in two steps
        // first calculate M*R, storing in tmp_m
        mul_mat_mat(tmp_m, 3u, x.data(i, offset+m*9u), ptsym.get(invR[i]).data());
        // next calculate R⁻¹*tmp_m, storing back in the x array
        mul_mat_mat(x.data(i, offset+m*9u), 3u, ptsym.get(r[i]).data(), tmp_m);
      }
    }
  }
  return true;
}

template<typename T>
bool InnerInterpolationData<T>::rip_axial(
  ArrayVector<T>& x, const PointSymmetry& ptsym, const std::vector<size_t>& r, const std::vector<size_t>& invR, const int nthreads
) const {
  omp_set_num_threads((nthreads > 0) ? nthreads : omp_get_max_threads());
  ElementsType no = this->count_scalars_vectors_matrices();
  if (!std::any_of(no.begin() + 1, no.end(), [](element_t n) {return n > 0; }))
      return false;
  T tmp_v[3], tmp_m[9];
  std::vector<T> detR;
  if (no[1])
      std::transform(r.begin(), r.end(), std::back_inserter(detR),
          [ptsym](const size_t z) { return static_cast<T>(matrix_determinant(ptsym.get(z).data())); });
  element_t offset, sp = this->branch_span();
  std::array<int,9> ident = {1,0,0, 0,1,0, 0,0,1};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = unsigned_to_signed<long long, size_t>(x.size());
#pragma omp parallel for default(none) shared(x,ptsym,r,invR,detR) private(offset, tmp_v, tmp_m) firstprivate(ident, no, sp, xsize) schedule(static)
  for (long long si = 0; si < xsize; ++si) {
    size_t i = signed_to_unsigned<size_t, long long>(si);
    if (!approx_matrix(3, ident.data(), ptsym.get(r[i]).data()))
    for (element_t b = 0; b < branches_; ++b) {
        // scalar elements do not need to be rotated, so skip them
        offset = b * sp + no[0];
        // rotate real vectors: since Q = Rᵀq + τ → R⁻¹*v
        for (element_t v = 0; v < no[1]; ++v) {
            mul_mat_vec(tmp_v, 3u, ptsym.get(invR[i]).data(), x.data(i, offset + v * 3u));
            for (int j = 0; j < 3; ++j) x.insert(detR[i]*tmp_v[j], i, offset + v * 3u + j);
        }
        offset += no[1] * 3u;
        for (element_t m = 0; m < no[2]; ++m) {
            // Calculate R⁻¹*M*R in two steps
            // first calculate M*R, storing in tmp_m
            mul_mat_mat(tmp_m, 3u, x.data(i, offset + m * 9u), ptsym.get(r[i]).data());
            // next calculate R⁻¹*tmp_m, storing back in the x array
            mul_mat_mat(x.data(i, offset + m * 9u), 3u, ptsym.get(invR[i]).data(), tmp_m);
        }
    }
  }
  return true;
}

template<class T, class R, class S = typename std::common_type<T,R>::type>
static std::complex<S> e_iqd(const LQVec<T>& q, const size_t i, const ArrayVector<R>& d, const size_t j){
  const double pi = 3.14159265358979323846;
  S dotqd{0};
  for (size_t k=0; k<3u; ++k) dotqd += q.getvalue(i,k) * d.getvalue(j,k);
  std::complex<S> i_2pi_x(0, 2*pi*dotqd);
  return std::exp(i_2pi_x);
}

template<typename T> template<typename R>
bool InnerInterpolationData<T>::rip_gamma_complex(
  ArrayVector<T>& x, const LQVec<R>& q, const GammaTable& pgt,
  const PointSymmetry& ptsym, const std::vector<size_t>& ridx, const std::vector<size_t>& invRidx,
  const int nthreads
) const {
  // construct a lambda to calculate the phase for qᵢ, and [R⁻¹xₖ - xᵥ]
  auto e_iqd_gt = [q,pgt](size_t i, element_t k, size_t r){
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
  verbose_update("InnerInterpolationData::rip_gamma_complex called with ",threads," threads");
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  ElementsType no = this->count_scalars_vectors_matrices();
  if (!std::any_of(no.begin()+1, no.end(), [](element_t n){return n>0;}))
    return false;
  // for Γ transformations of tensors, there *must* be N×N in total:
  size_t Nmat = static_cast<size_t>(std::sqrt(no[2]));
  if (no[2] != Nmat*Nmat){
    std::cout << "Atomic displacement Gamma transformation requires NxN 3x3 tensors!" << std::endl;
    return false;
  }
  element_t sp = this->branch_span();
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = unsigned_to_signed<long long, size_t>(x.size());
#pragma omp parallel for default(none) \
                         shared(x, q, pgt, ptsym, ridx, invRidx, e_iqd_gt) \
                         firstprivate(no, Nmat, sp, xsize) \
                         schedule(static)
  for (long long si=0; si<xsize; ++si){
    size_t i = signed_to_unsigned<size_t, long long>(si);
    element_t offset{0};
    for (element_t b=0; b<branches_; ++b){
      // scalar elements do not need to be rotated, so skip them
      offset = b*sp + no[0];
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
        for (element_t k=0; k<no[1]; ++k){
          // is q already expressed in the right lattice? hopefully!
          // if the vector rotates by Γ is *must* be complex, so T is a complex type
          // use its constructor to make i q⋅[R⁻¹xₖ - xᵥ] and calculate the (k,R) phase
          T phase = e_iqd_gt(i, k, invRidx[i]);
          mul_mat_vec(tmp_v, 3u, ptsym.get(invRidx[i]).data(), x.data(i, offset+k*3u));
          auto v0idx = 3u*pgt.F0(k, invRidx[i]);
          for (int j=0; j<3; ++j) tmpvecs[v0idx+j] = phase*tmp_v[j];
        }
        for (element_t j=0; j<no[1]*3u; ++j) x.insert(tmpvecs[j], i, offset+j);
        offset += no[1]*3u;
      }

      if (no[2]>0){
        T tmp_m[9], tmp_m2[9];
        std::vector<T> tmpmats(no[2]*9u);
        for (element_t n=0; n<Nmat; ++n){
          T rph = e_iqd_gt(i, n, ridx[i]);
          element_t v = static_cast<element_t>(pgt.F0(n, ridx[i]));
          for (element_t m=0; m<Nmat; ++m){

            T invRph = e_iqd_gt(i, m, invRidx[i]);
            element_t k = static_cast<element_t>(pgt.F0(m, invRidx[i]));
            // Calculate R⁻¹*M*R in two steps
            // first calculate M*R, storing in tmp_m
            mul_mat_mat(tmp_m, 3u, x.data(i, offset+(n*3u+m)*9u), ptsym.get(ridx[i]).data());
            // next calculate R⁻¹*tmp_m, storing in the temporary all matrix array
            mul_mat_mat(tmp_m2, 3u, ptsym.get(invRidx[i]).data(), tmp_m);
            // include the R R⁻¹ phase factor
            for (int j=0; j<9; ++j) tmpmats[(v*3u+k)*9u+j] = rph*invRph*tmp_m2[j];
          }
        }
        for (element_t j=0; j<no[2]*9u; ++j) x.insert(tmpmats[j], i, offset+j);
      }
    }
  }
  return true;
}


template<class T, class R> class InterpolationData{
  InnerInterpolationData<T> values_;
  InnerInterpolationData<R> vectors_;
  PermutationTable permutation_table_;
public:
  InterpolationData(): values_(), vectors_(), permutation_table_(0,0) {};
  //
  void validate_values() {
    if (values_.size()!=vectors_.size() || values_.branches()!=vectors_.branches())
      values_.setup_fake(vectors_.size(), vectors_.branches());
  }
  void validate_vectors() {
    if (values_.size()!=vectors_.size() || values_.branches()!=vectors_.branches())
      vectors_.setup_fake(values_.size(), values_.branches());
  }
  //
  size_t size() const {
    assert(values_.size() == vectors_.size());
    return values_.size();
  }
  const InnerInterpolationData<T>& values() const {return this->values_;}
  const InnerInterpolationData<R>& vectors() const {return this->vectors_;}
  element_t branches() const {
    assert(values_.branches() == vectors_.branches());
    return values_.branches();
  }
  void rotateslike(const RotatesLike values, const RotatesLike vectors) {
    values_.rotateslike(values);
    vectors_.rotateslike(vectors);
  }
  RotatesLike values_rotate_like(const RotatesLike a){ return values_.rotateslike(a); }
  RotatesLike vectors_rotate_like(const RotatesLike a){ return vectors_.rotateslike(a); }
  //
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  void interpolate_at(const std::vector<I>&, const std::vector<double>&, ArrayVector<T>&, ArrayVector<R>&, const size_t) const;
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  void interpolate_at(const std::vector<std::pair<I,double>>&, ArrayVector<T>&, ArrayVector<R>&, const size_t) const;
  //
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  std::vector<int_t> get_permutation(const I, const I) const;
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  std::vector<std::vector<int_t>> get_permutations(const std::vector<I>&) const;
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  std::vector<std::vector<int_t>> get_permutations(const std::vector<std::pair<I,double>>&) const;
  //
//  bool rotate_in_place(ArrayVector<T>& vals, ArrayVector<R>& vecs, const std::vector<std::array<int,9>>& r) const {
//    return values_.rotate_in_place(vals, r) && vectors_.rotate_in_place(vecs, r);
//  }
//  bool rotate_in_place(ArrayVector<T>& vals, ArrayVector<R>& vecs, const std::vector<std::array<int,9>>& r, const int n) const {
//    return values_.rotate_in_place(vals, r, n) && vectors_.rotate_in_place(vecs, r, n);
//  }
  //
  // Replace the data within this object.
  template<typename... A> void replace_value_data(A... args) {
    values_.replace_data(args...);
    this->validate_vectors();
    this->update_permutation_table();
  }
  template<typename... A> void replace_vector_data(A... args) {
    vectors_.replace_data(args...);
    this->validate_values();
    this->update_permutation_table();
  }
  // used during holding-object initialisation
  void initialize_permutation_table(const size_t nverts, const std::set<size_t>& keys){
    this->permutation_table_ = PermutationTable(nverts, this->branches(), keys);
  }
  void update_permutation_table(){
    // preserve the keys in the permutation table, if possible
    this->permutation_table_.refresh(this->size(), this->branches());
  }
  //
  void set_value_cost_info(const int csf, const int cvf, const ElementsCost& elcost){
    values_.set_cost_info(csf, cvf, elcost);
  }
  void set_vector_cost_info(const int csf, const int cvf, const ElementsCost& elcost){
    vectors_.set_cost_info(csf, cvf, elcost);
  }
  // create a string representation of the values and vectors
  std::string to_string() const {
    std::string str = "value " + values_.to_string() + " vector " + vectors_.to_string();
    return str;
  }
  // Calculate the Debye-Waller factor for the provided Q points and ion masses
  template<template<class> class A>
  ArrayVector<double> debye_waller(const A<double>& Q, const std::vector<double>& M, const double t_K) const;
  //
  template<typename I, typename S=typename CostTraits<T>::type, typename=std::enable_if_t<std::is_integral<I>::value> >
  std::vector<S> cost_matrix(const I i0, const I i1) const;
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value>>
  void permute_modes(const I i, const std::vector<int_t>& p){
    values_.permute_modes(i, p);
    vectors_.permute_modes(i, p);
  }
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value>>
  void inverse_permute_modes(const I i, const std::vector<int_t>& p){
    values_.inverse_permute_modes(i, p);
    vectors_.inverse_permute_modes(i, p);
  }
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value>>
  bool is_degenerate(const I idx) const {
    return values_.any_equal_modes(idx);
    // we could try and do something fancier, but it's probaby not useful.
  }
  void sort(void);
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value>>
  bool determine_permutation_ij(const I i, const I j, std::mutex& map_mutex);
private:
  ArrayVector<double> debye_waller_sum(const ArrayVector<double>& Q, const double t_K) const;
  ArrayVector<double> debye_waller_sum(const LQVec<double>& Q, const double beta) const{ return this->debye_waller_sum(Q.get_xyz(), beta); }
};


template<typename T, class R>
ArrayVector<double>
InterpolationData<T,R>::debye_waller_sum(const ArrayVector<double>& Q, const double t_K) const {
  const double hbar = 6.582119569E-13; // meV⋅s
  const double kB   = 8.617333252E-2; // meV⋅K⁻¹
  size_t nQ = Q.size();
  ElementsType vector_elements = vectors_.elements();
  size_t nIons = vector_elements[1] / 3u; // already checked to be correct
  ArrayVector<double> WdQ(nIons,nQ); // Wᵈ(Q) has nIons entries per Q point
  double coth_en, Q_dot_e_2;
  size_t values_span = values_.branch_span();
  size_t vector_span = vectors_.branch_span();
  size_t vector_nq = vectors_.size();
  element_t nbr = vectors_.branches();
  const double beta = kB*t_K; // meV
  const double pref{hbar*hbar/static_cast<double>(2*vector_nq)}; // meV²⋅s²
  // for each input Q point
  for (size_t Qidx=0; Qidx<nQ; ++Qidx){
    // and each ion
    for (size_t d=0; d<nIons; ++d){
      double qj_sum{0};
      // sum over all reduced q in the first Brillouin zone
      for (size_t q=0; q<vector_nq; ++q){
        // and over all 3*nIon branches at each q
        for (size_t j=0; j<nbr; ++j){
          // for each branch energy, find <2nₛ+1>/ħωₛ ≡ coth(2ħωₛβ)/ħωₛ
          coth_en = coth_over_en(values_.data().getvalue(q,j*values_span), beta);
          // and find |Q⋅ϵₛ|². Note: vector_product(x,y) *is* |x⋅y|²
          Q_dot_e_2 = vector_product(3u, Q.data(Qidx), vectors_.data().data(q,j*vector_span+1u+3u*d));
          // adding |Q⋅ϵₛ|²coth(2ħωₛβ)/ħωₛ to the sum over s for [Qidx, d]
          qj_sum += Q_dot_e_2 * coth_en;
        }
      }
      // with the sum over s complete, normalize by ħ²/2 divided by the number
      // of points in the Brillouin zone and store the result at W[Qidx, d];
      WdQ.insert(qj_sum*pref, Qidx, d);
    }
  }
  return WdQ;
}

template<typename T, class R> template<template<class> class A>
ArrayVector<double>
InterpolationData<T,R>::debye_waller(const A<double>& Q, const std::vector<double>& M, const double t_K) const {
  ElementsType vector_elements = vectors_.elements();
  size_t nIons = vector_elements[1] / 3u;
  if (0 == nIons || vector_elements[1] != nIons*3u)
    throw std::runtime_error("Debye-Waller factor requires 3-vector eigenvector(s).");
  if (M.size() != nIons)
    throw std::runtime_error("Debye-Waller factor requires an equal number of ions and masses.");
  ArrayVector<double> WdQ = this->debye_waller_sum(Q, t_K);
  ArrayVector<double> factor(1u, Q.size());
  double d_sum;
  for (size_t Qidx=0; Qidx<Q.size(); ++Qidx){
    d_sum = double(0);
    for (size_t d=0; d<nIons; ++d){
      d_sum += std::exp(WdQ.getvalue(Qidx, d)/M[d]);
    }
    factor.insert(d_sum*d_sum, Qidx);
  }
  return factor;
}

template<class T, class R> template<typename I, typename>
void
InterpolationData<T,R>::interpolate_at(
  const std::vector<I>& indices,
  const std::vector<double>& weights,
  ArrayVector<T>& values_out,
  ArrayVector<R>& vectors_out,
  const size_t to
) const {
  std::vector<std::vector<int_t>> permutations = this->get_permutations(indices);
  values_.interpolate_at(permutations, indices, weights, values_out, to, false);
  vectors_.interpolate_at(permutations, indices, weights, vectors_out, to, true);
}

template<class T, class R> template<typename I, typename>
void
InterpolationData<T,R>::interpolate_at(
  const std::vector<std::pair<I,double>>& indices_weights,
  ArrayVector<T>& values_out,
  ArrayVector<R>& vectors_out,
  const size_t to
) const {
  std::vector<std::vector<int_t>> permutations = this->get_permutations(indices_weights);
  values_.interpolate_at(permutations, indices_weights, values_out, to, false);
  vectors_.interpolate_at(permutations, indices_weights, vectors_out, to, true);
}

// template<class T, class R> template<typename I, typename>
// std::vector<int_t>
// InterpolationData<T,R>::get_permutation(const I i, const I j) const {
//   std::vector<int_t> perm;
//   /*
//   Since missing permutations are added, and might be the same for two threads,
//   we can not allow multiple threads to perform their get call at the same time
//   otherwise the table might end up with Nthread copies of every unique
//   permutation vector.
//   This is probably an unacceptable performance hit.
//   */
//   #pragma omp critical
//   {
//     perm = permutation_table_.safe_get(i, j);
//     // perm is empty if i:j is not present in the table
//     if (perm.empty()){
//       jv_permutation_fill(this->cost_matrix(i, j), perm);
//       permutation_table_.set(i, j, perm);
//     }
//   }
//   // jv_permutation_fill(this->cost_matrix(i, j), perm);
//   // perm = jv_permutation(this->cost_matrix(i, j));
//   // and return it
//   return perm;
// }
template<class T, class R> template<typename I, typename>
std::vector<int_t>
InterpolationData<T,R>::get_permutation(const I i, const I j) const {
  return permutation_table_.safe_get(i, j);
}

template<class T, class R> template<typename I, typename>
std::vector<std::vector<int_t>>
InterpolationData<T,R>::get_permutations(const std::vector<I>& indices) const {
  std::vector<std::vector<int_t>> perms;
  // find the minimum index so that permutation(pvt,idx) is always ordered
  I pvt{indices[0]};
  //for (const I idx: indices) if (idx < pvt) pvt = idx;
  for (const I idx: indices) perms.push_back(this->get_permutation(pvt, idx));
  return perms;
}
template<class T, class R> template<typename I, typename>
std::vector<std::vector<int_t>>
InterpolationData<T,R>::get_permutations(const std::vector<std::pair<I,double>>& iw) const {
  std::vector<std::vector<int_t>> perms;
  I pvt{iw[0].first};
  //for (const auto piw: iw) if (piw.first < pvt) pvt = piw.first;
  for (const auto piw: iw) perms.push_back(this->get_permutation(pvt, piw.first));
  return perms;
}

template<class T, class R> template<typename I, typename S, typename>
std::vector<S>
InterpolationData<T,R>::cost_matrix(const I i0, const I i1) const {
  element_t Nbr{this->branches()};
  std::vector<S> cost(Nbr*Nbr, S(0));
  if (i0==i1){
    for (element_t j=0; j<Nbr*Nbr; j+=Nbr+1) cost[j] = S(-1);
  } else {
    values_.add_cost(i0, i1, cost, false);
    vectors_.add_cost(i0, i1, cost, true);
  }
  return cost;
}

template<class T, class R>
void InterpolationData<T,R>::sort(void){
  std::set<size_t> keys = permutation_table_.keys();
  // find the keys corresponding to one triangular part of the matrix (i<j)
  std::vector<std::array<size_t,2>> tri_ij;
  tri_ij.reserve(keys.size()/2);
  size_t no = this->size();
  for (const auto & key: keys){
    size_t i = key/no;
    if (i*(no+1) < key) tri_ij.push_back({i, key-i*no});
  }
  info_update("Finding permutations for ",keys.size()," connections between the ",no," vertices");
  // now find the permutations in parallel
  std::mutex m;
  long long nok = unsigned_to_signed<long long, size_t>(tri_ij.size());
  #pragma omp parallel for default(none) shared(tri_ij, m, nok)
  for (long long sk=0; sk<nok; ++sk){
    size_t k = signed_to_unsigned<size_t, long long>(sk);
    this->determine_permutation_ij(tri_ij[k][0], tri_ij[k][1], m);
  }
  info_update("Done");
}

template<class T, class R> template<typename I, typename>
bool
InterpolationData<T,R>::determine_permutation_ij(const I i, const I j, std::mutex& map_mutex){
  // if (!permutation_table_.value_needed(i,j)) return false;
  std::vector<int_t> row, col;
  jv_permutation_fill(this->cost_matrix(i,j), row, col);
  std::unique_lock<std::mutex> map_lock(map_mutex);
  permutation_table_.overwrite(i, j, row);
  permutation_table_.overwrite(j, i, col);
  map_lock.unlock();
  return true;
}

#endif
