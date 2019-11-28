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
#include "arrayvector.hpp"

#ifndef _INTERPOLATION_DATA_H_
#define _INTERPOLATION_DATA_H_

template<class T> class InterpolationData{
  typedef unsigned long element_t;
  typedef std::vector<size_t> ShapeType;
  typedef std::array<element_t,4> ElementsType;
  ArrayVector<T> data_;   //!< The stored ArrayVector indexed like the holding-PolyhedronTrellis' vertices
  ShapeType shape_;       //!< A std::vector to indicate a possible higher-dimensional shape of each `data` array
  ElementsType elements_; //!< The number of scalars, normalised eigenvector elements, vector elements, and matrix elements per data array
  element_t branches_;    //!< The number of branches contained per data array
public:
  InterpolationData(): data_({0,0}), shape_({0,0}), elements_({0,0,0,0}), branches_(0){};
  size_t size(void) const {return data_.size();}
  size_t numel(void) const {return data_.numel();}
  const ArrayVector<T>& data(void) const {return data_;}
  const ShapeType& shape(void) const {return shape_;}
  const ElementsType& elements(void) const {return elements_;}
  element_t branches(void) const {return branches_;}
  //
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  void interpolate_at(const std::vector<I>&, const std::vector<double>&, ArrayVector<T>&, const size_t) const;
  //
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  void interpolate_at(const std::vector<std::pair<I,double>>&, ArrayVector<T>&,const size_t) const;
  //
  bool rotate_in_place(ArrayVector<T>&, const std::vector<std::array<int,9>>&) const;
  bool rotate_in_place(ArrayVector<T>&, const std::vector<std::array<int,9>>&, const int) const;
  //
  void replace_data(const ArrayVector<T>&, const ShapeType&, const ElementsType&);
  void replace_data(const ArrayVector<T>& nd, const ElementsType& ne=ElementsType({0,0,0,0})){
    ShapeType ns{nd.size(), nd.numel()};
    return this->replace_data(nd, ns, ne);
  }
  // Calculate the Debye-Waller factor for the provided Q points and ion masses
  template<template<class> class A>
  ArrayVector<double> debye_waller(const A<double>& Q, const std::vector<double>& M, const double t_K) const;
private:
  ArrayVector<double> debye_waller_sum(const LQVec<double>& Q, const double beta) const;
  ArrayVector<double> debye_waller_sum(const ArrayVector<double>& Q, const double t_K) const;
  element_t branch_span(void) const { return this->branch_span(elements_);}
  element_t branch_span(const ElementsType& e) const { return e[0]+e[1]+e[2]+e[3]*e[3]; }
  ElementsType count_scalars_eigenvectors_vectors_matrices(void) const {
    ElementsType no{elements_[0],0,0,0};
    for (int i=1; i<4; ++i) no[i] = elements_[i]/3u;
    return no;
  }
};

template<typename T> template<typename I, typename>
void InterpolationData<T>::interpolate_at(
  const std::vector<I>& indices,
  const std::vector<double>& weights,
  ArrayVector<T>& out,
  const size_t to
) const {
  if (indices.size()==0 || weights.size()==0)
    throw std::logic_error("Interpolation requires input data!");
  T *out_to = out.data(to), *ptr0 = data_.data(indices[0]);
  element_t span = this->branch_span();
  for (size_t x=0; x<indices.size(); ++x){
    T *ptrX = data_.data(indices[x]);
    // loop over the branches:
    for (element_t b=0; b < branches_; ++b){
      // find the weighted sum of each scalar:
      for (element_t s=0; s < elements_[0]; ++s) out_to[b*span+s] += weights[x]*ptrX[b*span+s];
      // handle any eigenvectors
      if (elements_[1]){
        element_t o1 = b*span + elements_[0];
        // find the arbitrary phase eⁱᶿ between different point eigenvectors
        T eith = antiphase(hermitian_product(elements_[1], ptr0+o1, ptrX+o1));
        // remove the arbitrary phase while adding the weighted value
        for (element_t e=0; e<elements_[1]; ++e) out_to[o1+e] += weights[x]*(eith*ptrX[o1+e]);
      }
      // handle any remaining standard vector or matrix elements
      for (element_t r=elements_[0]+elements_[1]; r<span; ++r)
        out_to[b*span+r] += weights[x]*ptrX[b*span+r];
    }
  }
  // ensure that any eigenvectors are still normalised
  if (elements_[1]) for (element_t b=0; b<branches_; ++b){
    element_t o2 = b*span + elements_[0];
    auto normI = std::sqrt(inner_product(elements_[1], out_to+o2, out_to+o2));
    for (element_t e=0; e<elements_[1]; ++e) out_to[o2+e]/=normI;
  }
}

template<typename T> template<typename I, typename>
void InterpolationData<T>::interpolate_at(
  const std::vector<std::pair<I,double>>& indices_weights,
  ArrayVector<T>& out,
  const size_t to
) const {
  if (indices_weights.size()==0)
    throw std::logic_error("Interpolation requires input data!");
  T *out_to = out.data(to), *ptr0 = data_.data(indices_weights[0].first);
  element_t span = this->branch_span();
  for (auto iw: indices_weights){
    T *ptrX = data_.data(iw.first);
    // loop over the branches:
    for (element_t b=0; b < branches_; ++b){
      // find the weighted sum of each scalar:
      for (element_t s=0; s < elements_[0]; ++s) out_to[b*span+s] += iw.second*ptrX[b*span+s];
      // handle any eigenvectors
      if (elements_[1]){
        element_t o1 = b*span + elements_[0];
        // find the arbitrary phase eⁱᶿ between different point eigenvectors
        T eith = antiphase(hermitian_product(elements_[1], ptr0+o1, ptrX+o1));
        // remove the arbitrary phase while adding the weighted value
        for (element_t e=0; e<elements_[1]; ++e) out_to[o1+e] += iw.second*(eith*ptrX[o1+e]);
      }
      // handle any remaining standard vector or matrix elements
      for (element_t r=elements_[0]+elements_[1]; r<span; ++r)
        out_to[b*span+r] += iw.second*ptrX[b*span+r];
    }
  }
  // ensure that any eigenvectors are still normalised
  if (elements_[1]) for (element_t b=0; b<branches_; ++b){
    element_t o2 = b*span + elements_[0];
    auto normI = std::sqrt(inner_product(elements_[1], out_to+o2, out_to+o2));
    for (element_t e=0; e<elements_[1]; ++e) out_to[o2+e]/=normI;
  }
}

template<typename T>
void InterpolationData<T>::replace_data(
  const ArrayVector<T>& nd,
  const ShapeType& ns,
  const ElementsType& ne
){
  data_ = nd;
  shape_ = ns;
  elements_ = ne;
  if (ne[1]%3)
    throw std::logic_error("Eigenvectors must have 3 elements per atom per branch");
  if (ne[2]%3)
    throw std::logic_error("Vectors must consist of 3N elements per branch");
  if (ne[3] && ne[3]!=3u)
    throw std::logic_error("Matrices must be 3x3 for each branch");
  // check the input for correctness
  element_t total_elements = 1u;
  // scalar + eigenvector + vector + matrix*matrix elements
  element_t known_elements = this->branch_span(ne);
  // no matter what, shape[0] should be the number of gridded points
  if (ns.size()>2){
    // if the number of dimensions of the shape array is greater than two,
    // the second element is the number of modes per point                    */
    branches_ = ns[1];
    for (size_t i=2u; i<ns.size(); ++i) total_elements *= ns[i];
  } else {
    // shape is [n_points, n_elements] or [n_points,], so there is only one mode
    branches_ = 1u;
    total_elements = ns.size() > 1 ? ns[1] : 1u;
  }
  if (0 == known_elements) elements_[0] = total_elements;
  if (known_elements && known_elements != total_elements){
    std::string msg;
    msg = "Inconsistent elements: " + std::to_string(known_elements) + " = ";
    msg += std::to_string(elements_[0]) + "+" + std::to_string(elements_[1]) + "+";
    msg += std::to_string(elements_[2]) + "+" + std::to_string(elements_[3]) + "² ≠ ";
    msg += std::to_string(total_elements);
    throw std::runtime_error(msg);
  }
}


template<typename T>
ArrayVector<double> InterpolationData<T>::debye_waller_sum(
  const LQVec<double>& Q,
  const double t_K
) const{
  return this->debye_waller_sum(Q.get_xyz(), t_K);
}

template<typename T>
ArrayVector<double> InterpolationData<T>::debye_waller_sum(
  const ArrayVector<double>& Q,
  const double t_K
) const {
  const double hbar = 6.582119569E-13; // meV⋅s
  const double kB   = 8.617333252E-2; // meV⋅K⁻¹
  size_t nQ = Q.size();
  size_t nIons = elements_[1] / 3u; // already checked to be correct
  ArrayVector<double> WdQ(nIons,nQ); // Wᵈ(Q) has nIons entries per Q point
  double coth_en, Q_dot_e_2;
  size_t span = 1u + nIons*3u + elements_[2] + elements_[3]*elements_[3];
  size_t nq = shape_[0];

  const double beta = kB*t_K; // meV
  const double pref{hbar*hbar/static_cast<double>(2*nq)}; // meV²⋅s²

  double qj_sum;
  // for each input Q point
  for (size_t Qidx=0; Qidx<nQ; ++Qidx){
    // and each ion
    for (size_t d=0; d<nIons; ++d){
      qj_sum = double(0);
      // sum over all reduced q in the first Brillouin zone
      for (size_t q=0; q<nq; ++q){
        // and over all 3*nIon branches at each q
        for (size_t j=0; j<branches_; ++j){
          // for each branch energy, find <2nₛ+1>/ħωₛ ≡ coth(2ħωₛβ)/ħωₛ
          coth_en = coth_over_en(data_.getvalue(q,j*span), beta);
          // and find |Q⋅ϵₛ|². Note: vector_product(x,y) *is* |x⋅y|²
          Q_dot_e_2 = vector_product(3u, Q.data(Qidx), data_.data(q,j*span+1u+3u*d));
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

template<typename T>
template<template<class> class A>
ArrayVector<double> InterpolationData<T>::debye_waller(
  const A<double>& Q,
  const std::vector<double>& M,
  const double t_K
) const {
  size_t nIons = elements_[1] / 3u;
  if (0 == nIons || elements_[1]*3u != nIons)
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

template<typename T>
bool InterpolationData<T>::rotate_in_place(
  ArrayVector<T>& x, const std::vector<std::array<int,9>>& r
) const {
  ElementsType no = this->count_scalars_eigenvectors_vectors_matrices();
  if (!std::any_of(no.begin()+1, no.end(), [](element_t n){return n>0;}))
    return false;
  T tmp_v[3], tmp_m[9];
  std::vector<std::array<int,9>> invR;
  if (no[3]){
    invR.resize(r.size());
    for (size_t i=0; i<r.size(); ++i) matrix_inverse(invR[i].data(), r[i].data());
  }
  element_t offset, nvec = no[1]+no[2], sp = this->branch_span();
  for (size_t i=0; i<x.size(); ++i)
  if (r[i][0] + r[i][4] + r[i][8] != 3)
  for (element_t b=0; b<branches_; ++b){
    // scalar elements do not need to be rotated, so skip them
    offset = b*sp + no[0];
    // eigenvectors and regular vectors rotate the same way
    for (element_t v=0; v<nvec; ++v){
      mul_mat_vec(tmp_v, 3u, r[i].data(), x.data(i, offset+v*3u));
      for (int j=0; j<3; ++j) x.insert(tmp_v[j], i, offset+v*3u+j);
    }
    offset += nvec*3u;
    for (element_t m=0; m<no[3]; ++m){
      // Calculate R*M*R⁻¹ in two steps
      // first calculate M*R⁻¹, storing in tmp_m
      mul_mat_mat(tmp_m, 3u, x.data(i, offset+m*9u), invR[i].data());
      // next calculate R*tmp_m, storing back in the x array
      mul_mat_mat(x.data(i, offset+m*9u), 3u, r[i].data(), tmp_m);
    }
  }
  return true;
}

template<typename T>
bool InterpolationData<T>::rotate_in_place(
  ArrayVector<T>& x, const std::vector<std::array<int,9>>& r, const int nthreads
) const {
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  ElementsType no = this->count_scalars_eigenvectors_vectors_matrices();
  if (!std::any_of(no.begin()+1, no.end(), [](element_t n){return n>0;}))
    return false;
  T tmp_v[3], tmp_m[9];
  std::vector<std::array<int,9>> invR;
  if (no[3]){
    long long rsize = unsigned_to_signed<long long, size_t>(r.size());
    invR.resize(r.size());
    #pragma omp parallel for default(none) shared(invR, r) firstprivate(rsize) schedule(dynamic)
    for (long long si=0; si<rsize; ++si){
      size_t i = signed_to_unsigned<size_t, long long>(si);
      matrix_inverse(invR[i].data(), r[i].data());
    }
  }
  element_t offset, nvec = no[1]+no[2], sp = this->branch_span();
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = unsigned_to_signed<long long, size_t>(x.size());
#pragma omp parallel for default(none) shared(x,r,invR) private(offset, tmp_v, tmp_m) firstprivate(nvec, no, sp, xsize) schedule(static)
  for (long long si=0; si<xsize; ++si){
    size_t i = signed_to_unsigned<size_t, long long>(si);
    if (r[i][0] + r[i][4] + r[i][8] != 3)
    for (element_t b=0; b<branches_; ++b){
      // scalar elements do not need to be rotated, so skip them
      offset = b*sp + no[0];
      // eigenvectors and regular vectors rotate the same way
      for (element_t v=0; v<nvec; ++v){
        mul_mat_vec(tmp_v, 3u, r[i].data(), x.data(i, offset+v*3u));
        for (int j=0; j<3; ++j) x.insert(tmp_v[j], i, offset+v*3u+j);
      }
      offset += nvec*3u;
      for (element_t m=0; m<no[3]; ++m){
        // Calculate R*M*R⁻¹ in two steps
        // first calculate M*R⁻¹, storing in tmp_m
        mul_mat_mat(tmp_m, 3u, x.data(i, offset+m*9u), invR[i].data());
        // next calculate R*tmp_m, storing back in the x array
        mul_mat_mat(x.data(i, offset+m*9u), 3u, r[i].data(), tmp_m);
      }
    }
  }
  return true;
}
 #endif
