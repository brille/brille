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
#ifndef BRILLE_DUALINTERPOLATOR_H_
#define BRILLE_DUALINTERPOLATOR_H_
#include "interpolator2.hpp"
namespace brille {

template<class T, class R>
class DualInterpolator{
public:
  using ind_t = brille::ind_t;
  template<class Z> using element_t = std::array<Z,3>;
private:
  Interpolator<T> values_;
  Interpolator<R> vectors_;
  PermutationTable permutation_table_;
public:
  DualInterpolator(): values_(), vectors_(), permutation_table_(0,0) {};
  DualInterpolator(Interpolator<T>& val, Interpolator<R>& vec)
  : values_(val), vectors_(vec), permutation_table_(0,0)
  {
    if (values_.size() != vectors_.size() || values_.branches() != vectors_.branches())
      throw std::runtime_error("Inconsistent values and vectors provided to DualInterpolator");
  }
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
  const Interpolator<T>& values() const {return this->values_;}
  const Interpolator<R>& vectors() const {return this->vectors_;}
  ind_t branches() const {
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
  // template<template<class> class A>
  // void
  // interpolate_at(const std::vector<ind_t>& indices, const std::vector<double>& weights, A<T>& values_out, A<R>& vectors_out, const ind_t to) const {
  //   auto permutations = this->get_permutations(indices);
  //   values_.interpolate_at(permutations, indices, weights, values_out, to, false);
  //   vectors_.interpolate_at(permutations, indices, weights, vectors_out, to, true);
  // }
  template<template<class> class A>
  void
  interpolate_at(const std::vector<std::pair<ind_t,double>>& indices_weights, A<T>& values_out, A<R>& vectors_out, const ind_t to) const {
    auto permutations = this->get_permutations(indices_weights);
    values_.interpolate_at(permutations, indices_weights, values_out, to, false);
    vectors_.interpolate_at(permutations, indices_weights, vectors_out, to, true);
  }
  //
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  std::vector<typename PermutationTable::ind_t>
  get_permutation(const I, const I) const;
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  std::vector<std::vector<typename PermutationTable::ind_t>>
  get_permutations(const std::vector<I>&) const;
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  std::vector<std::vector<typename PermutationTable::ind_t>>
  get_permutations(const std::vector<std::pair<I,double>>&) const;
  //
  // Replace the data within this object.
  void replace_data(Interpolator<T>& val, Interpolator<R>& vec){
    if (vec.size() != val.size() || vec.branches() != val.branches())
      throw std::runtime_error("The values and vectors must have matching number of points and branches");
    values_ = val;
    vectors_ = vec;
    this->update_permutation_table();
  }
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
  void set_value_cost_info(const int csf, const int cvf, const element_t<double>& elcost){
    values_.set_cost_info(csf, cvf, elcost);
  }
  void set_vector_cost_info(const int csf, const int cvf, const element_t<double>& elcost){
    vectors_.set_cost_info(csf, cvf, elcost);
  }
  // create a string representation of the values and vectors
  std::string to_string() const {
    std::string str = "value " + values_.to_string() + " vector " + vectors_.to_string();
    return str;
  }
  // Calculate the Debye-Waller factor for the provided Q points and ion masses
  // This returns a 1-D brille::Array (so it can't be a brille::Array2)
  template<template<class> class A>
  brille::Array<double> debye_waller(const A<double>& Qpts, const std::vector<double>& Masses, const double t_K) const;
  //
  template<typename I, typename S=typename CostTraits<T>::type, typename=std::enable_if_t<std::is_integral<I>::value> >
  std::vector<S> cost_matrix(const I i0, const I i1) const;
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value>>
  void permute_modes(const I i, const std::vector<ind_t>& p){
    values_.permute_modes(i, p);
    vectors_.permute_modes(i, p);
  }
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value>>
  void inverse_permute_modes(const I i, const std::vector<ind_t>& p){
    values_.inverse_permute_modes(i, p);
    vectors_.inverse_permute_modes(i, p);
  }
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value>>
  bool is_degenerate(const I idx) {
    return values_.any_equal_modes(idx);
    // we could try and do something fancier, but it's probaby not useful.
  }
  void sort(void);
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value>>
  bool determine_permutation_ij(const I i, const I j, std::mutex& map_mutex);
  size_t bytes_per_point() const {
    return values_.bytes_per_point() + vectors_.bytes_per_point();
  }
private:
  bArray<double> debye_waller_sum(const bArray<double>& Qpts, const double t_K) const;
  bArray<double> debye_waller_sum(const LQVec<double>& Qpts, const double beta) const{ return this->debye_waller_sum(Qpts.get_xyz(), beta); }
};


template<class T, class R>
bArray<double>
DualInterpolator<T,R>::debye_waller_sum(const bArray<double>& Qpts, const double t_K) const {
  const double hbar = 6.582119569E-13; // meV⋅s
  const double kB   = 8.617333252E-2; // meV⋅K⁻¹
  size_t nQ = Qpts.size(0);
  element_t<ind_t> vector_elements = vectors_.elements();
  size_t nIons = vector_elements[1] / 3u; // already checked to be correct
  bArray<double> WdQ(nQ,nIons); // Wᵈ(Q) has nIons entries per Q point
  double coth_en, Q_dot_e_2;
  size_t vector_nq = vectors_.size();
  ind_t nbr = vectors_.branches();
  const double beta = kB*t_K; // meV
  const double pref{hbar*hbar/static_cast<double>(2*vector_nq)}; // meV²⋅s²
  const auto& val{ values_.data()};
  const auto& vec{vectors_.data()};
  ind_t vecsp = vectors_.branch_span();
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
          coth_en = brille::utils::coth_over_en(val.val(q,j), beta);
          // and find |Q⋅ϵₛ|². Note: brille::utils::vector_product(x,y) *is* |x⋅y|²
          Q_dot_e_2 = brille::utils::vector_product(3u, Qpts.ptr(Qidx), vec.ptr(q,j*vecsp+3u*d));
          // adding |Q⋅ϵₛ|²coth(2ħωₛβ)/ħωₛ to the sum over s for [Qidx, d]
          qj_sum += Q_dot_e_2 * coth_en;
        }
      }
      // with the sum over s complete, normalize by ħ²/2 divided by the number
      // of points in the Brillouin zone and store the result at W[Qidx, d];
      WdQ.val(Qidx,d) = qj_sum*pref;
    }
  }
  return WdQ;
}

template<class T, class R>
template<template<class> class A>
brille::Array<double>
DualInterpolator<T,R>::debye_waller(const A<double>& Qpts, const std::vector<double>& Masses, const double t_K) const {
  element_t<ind_t> vector_elements = vectors_.elements();
  size_t nIons = vector_elements[1] / 3u;
  if (0 == nIons || vector_elements[1] != nIons*3u)
    throw std::runtime_error("Debye-Waller factor requires 3-vector eigenvector(s).");
  if (Masses.size() != nIons)
    throw std::runtime_error("Debye-Waller factor requires an equal number of ions and masses.");
  auto WdQ = this->debye_waller_sum(Qpts, t_K); // {nQ, nAtoms}
  brille::shape_t fshape{Qpts.size(0)}; // (nQ,)
  brille::Array<double> factor(fshape); // we don't want an Array2 here!
  double d_sum;
  for (size_t Qidx=0; Qidx<Qpts.size(0); ++Qidx){
    d_sum = double(0);
    for (size_t d=0; d<nIons; ++d)
      d_sum += std::exp(WdQ.val(Qidx,d)/Masses[d]);
    factor[Qidx] = d_sum*d_sum;
  }
  return factor;
}



// template<class T, class R> template<typename I, typename>
// std::vector<PermutationTable::ind_t>
// DualInterpolator<T,R>::get_permutation(const I i, const I j) const {
//   std::vector<PermutationTable::ind_t> perm;
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
std::vector<typename PermutationTable::ind_t>
DualInterpolator<T,R>::get_permutation(const I i, const I j) const {
  return permutation_table_.safe_get(i, j);
}

template<class T, class R> template<typename I, typename>
std::vector<std::vector<typename PermutationTable::ind_t>>
DualInterpolator<T,R>::get_permutations(const std::vector<I>& indices) const {
  std::vector<std::vector<PermutationTable::ind_t>> perms;
  // find the minimum index so that permutation(pvt,idx) is always ordered
  I pvt{indices[0]};
  //for (const I idx: indices) if (idx < pvt) pvt = idx;
  for (const I & idx: indices) perms.push_back(this->get_permutation(pvt, idx));
  return perms;
}
template<class T, class R> template<typename I, typename>
std::vector<std::vector<typename PermutationTable::ind_t>>
DualInterpolator<T,R>::get_permutations(const std::vector<std::pair<I,double>>& iw) const {
  std::vector<std::vector<typename PermutationTable::ind_t>> perms;
  I pvt{iw[0].first};
  //for (const auto piw: iw) if (piw.first < pvt) pvt = piw.first;
  for (const auto & piw: iw) perms.push_back(this->get_permutation(pvt, piw.first));
  return perms;
}

template<class T, class R> template<typename I, typename S, typename>
std::vector<S>
DualInterpolator<T,R>::cost_matrix(const I i0, const I i1) const {
  ind_t Nbr{this->branches()};
  std::vector<S> cost(Nbr*Nbr, S(0));
  if (i0==i1){
    for (ind_t j=0; j<Nbr*Nbr; j+=Nbr+1) cost[j] = S(-1);
  } else {
     values_.add_cost(i0, i1, cost, false);
    vectors_.add_cost(i0, i1, cost, true);
  }
  return cost;
}

template<class T, class R>
void
DualInterpolator<T,R>::sort(void){
  std::set<size_t> keys = permutation_table_.keys();
  // find the keys corresponding to one triangular part of the matrix (i<j)
  std::vector<std::array<size_t,2>> tri_ij;
  tri_ij.reserve(keys.size()/2);
  size_t no = this->size();
  for (const auto & key: keys){
    size_t i = key/no;
    if (i*(no+1) < key) tri_ij.push_back({i, key-i*no});
  }
  debug_update("Finding permutations for ",keys.size()," connections between the ",no," vertices");
  // now find the permutations in parallel
  std::mutex m;
  long long nok = brille::utils::u2s<long long, size_t>(tri_ij.size());
  #pragma omp parallel for default(none) shared(tri_ij, m, nok)
  for (long long sk=0; sk<nok; ++sk){
    size_t k = brille::utils::s2u<size_t, long long>(sk);
    this->determine_permutation_ij(tri_ij[k][0], tri_ij[k][1], m);
  }
  debug_update("Done");
}

template<class T, class R> template<typename I, typename>
bool
DualInterpolator<T,R>::determine_permutation_ij(const I i, const I j, std::mutex& map_mutex){
  // if (!permutation_table_.value_needed(i,j)) return false;
  std::vector<int> row, col; // jv_permutation has difficulty with unsigned integers
  jv_permutation_fill(this->cost_matrix(i,j), row, col);
  std::unique_lock<std::mutex> map_lock(map_mutex);
  permutation_table_.overwrite(i, j, row);
  permutation_table_.overwrite(j, i, col);
  map_lock.unlock();
  return true;
}

} // namespace brille
#endif
