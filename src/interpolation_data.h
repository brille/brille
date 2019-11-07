#include <vector>
#include <array>
#include <utility>
#include "arrayvector.h"

#ifndef _INTERPOLATION_DATA_H_
#define _INTERPOLATION_DATA_H_

template<class T> class InterpolationData{
  typedef std::vector<size_t> ShapeType;
  typedef std::array<unsigned,4> ElementsType;
  ArrayVector<T> data_;   //!< The stored ArrayVector indexed like the holding-PolyhedronTrellis' vertices
  ShapeType shape_;       //!< A std::vector to indicate a possible higher-dimensional shape of each `data` array
  ElementsType elements_; //!< The number of scalars, normalised eigenvector elements, vector elements, and matrix elements per data array
  size_t branches_;       //!< The number of branches contained per data array
public:
  InterpolationData(): data_({0,0}), shape_({0,0}), elements_({0,0,0,0}), branches_(0){};
  size_t size(void) const {return data_.size();}
  size_t numel(void) const {return data_.numel();}
  const ArrayVector<T>& data(void) const {return data_;}
  const ShapeType& shape(void) const {return shape_;}
  const ElementsType& elements(void) const {return elements_;}
  size_t branches(void) const {return branches_;}
  //
  void interpolate_at(const std::vector<size_t>&, const std::vector<double>&, ArrayVector<T>&, const size_t) const;
  void interpolate_at(const std::vector<std::pair<size_t,double>>&, ArrayVector<T>&,const size_t) const;
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
};

template<typename T>
void InterpolationData<T>::interpolate_at(
  const std::vector<size_t>& indices,
  const std::vector<double>& weights,
  ArrayVector<T>& out,
  const size_t to
) const {
  if (indices.size()==0 || weights.size()==0)
    throw std::logic_error("Interpolation requires input data!");
  new_unsafe_interpolate_to(data_, elements_, branches_, indices, weights, out, to);
}
template<typename T>
void InterpolationData<T>::interpolate_at(
  const std::vector<std::pair<size_t,double>>& indices_weights,
  ArrayVector<T>& out,
  const size_t to
) const {
  if (indices_weights.size()==0)
    throw std::logic_error("Interpolation requires input data!");
  std::vector<size_t> indices;
  std::vector<double> weights;
  for (auto iw: indices_weights){
    indices.push_back(iw.first);
    weights.push_back(iw.second);
  }
  new_unsafe_interpolate_to(data_, elements_, branches_, indices, weights, out, to);
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
  // check the input for correctness
  size_t total_elements = 1u;
  // scalar + eigenvector + vector + matrix*matrix elements
  size_t known_elements = static_cast<size_t>(ne[0])+static_cast<size_t>(ne[1])+static_cast<size_t>(ne[2])
                        + static_cast<size_t>(ne[3])*static_cast<size_t>(ne[3]);
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


 #endif
