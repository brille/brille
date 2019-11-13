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

template<class T>
template<class R, class S>
ArrayVector<S> MapGrid3<T>::debye_waller_sum(const LQVec<R>& Q, const R t_K) const {
  return this->debye_waller_sum(Q.get_xyz(), t_K);
}
// T is likely complex, so we need CostTraits<T> to determine the common type S.
template<class T>
template<class R, class S>
ArrayVector<S> MapGrid3<T>::debye_waller_sum(const ArrayVector<R>& Q, const R t_K) const {
  const S hbar = 6.582119569E-13; // meV⋅s
  const S kB   = 8.617333252E-2; // meV⋅K⁻¹
  if (Q.numel() != 3)
    throw std::runtime_error("Debye-Waller factor requires 3-vector Q.");
  if (this->elements[0] != 1u)
    throw std::runtime_error("Debye-Waller factor requires one scalar (energy) per mode.");
  size_t nIons = this->elements[1] / 3u;
  if (0 == nIons || this->elements[1]*3u != nIons)
    throw std::runtime_error("Debye-Waller factor requires 3-vector eigenvector(s).");
  size_t nQ = Q.size();
  ArrayVector<S> WdQ(nIons,nQ); // Wᵈ(Q) has nIons entries per Q point

  S coth_en, Q_dot_e_2;
  size_t stride = 1u + nIons*3u + this->elements[2] + this->elements[3]*this->elements[3];
  size_t nq = this->shape.getvalue(0u);

  const S beta = kB*t_K; // meV
  const S pref{hbar*hbar/static_cast<S>(2*nq)}; // meV²⋅s²

  S qj_sum;
  // for each input Q point
  for (size_t Qidx=0; Qidx<nQ; ++Qidx){
    // and each ion
    for (size_t d=0; d<nIons; ++d){
      qj_sum = S(0);
      // sum over all reduced q in the first Brillouin zone
      for (size_t q=0; q<nq; ++q){
        // and over all 3*nIon branches at each q
        for (size_t j=0; j<this->branches; ++j){
          // for each branch energy, find <2nₛ+1>/ħωₛ ≡ coth(2ħωₛβ)/ħωₛ
          coth_en = coth_over_en(this->data.getvalue(q,j*stride), beta);
          // and find |Q⋅ϵₛ|². Note: vector_product(x,y) *is* |x⋅y|²
          Q_dot_e_2 = vector_product(3u, Q.data(Qidx), this->data.data(q,j*stride+1u+3u*d));
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

/*! \brief Calculate the Debye-Waller factor for one or more Q points.

@param Q An array of N 3-vectors expressed in reciprocal lattice units or
         inverse angstrom. If rlu are provided Q must be a LQVec such that Å⁻¹
         can be calculated.
@param M The ion masses in meV⋅s²⋅Å⁻² // or should we take amu and convert? 1 amu == 1.03642688E-25 meV s² Å⁻²
@param T The temperature in K.
@returns The Debye-Waller factor at each of the N Q points, which is unitless.
*/
template<class T>
template<class R, template<class> class A, class S>
ArrayVector<S> MapGrid3<T>::debye_waller(const A<R>& Q, const std::vector<R>& M, const R t_K) const {
  size_t nIons = this->elements[1] / 3u;
  if (0 == nIons || this->elements[1]*3u != nIons)
    throw std::runtime_error("Debye-Waller factor requires 3-vector eigenvector(s).");
  if (M.size() != nIons)
    throw std::runtime_error("Debye-Waller factor requires an equal number of ions and masses.");
  ArrayVector<S> WdQ = this->debye_waller_sum(Q, t_K);
  ArrayVector<S> factor(1u, Q.size());
  S d_sum;
  for (size_t Qidx=0; Qidx<Q.size(); ++Qidx){
    d_sum = S(0);
    for (size_t d=0; d<nIons; ++d){
      d_sum += std::exp(WdQ.getvalue(Qidx, d)/M[d]);
    }
    factor.insert(d_sum*d_sum, Qidx);
  }
  return factor;
}
