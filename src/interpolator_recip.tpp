/* This file is part of brille.

Copyright © 2020 Greg Tucker <greg.tucker@stfc.ac.uk>

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

template<class T>
bool Interpolator<T>::rip_recip(
  bArray<T>& x, const PointSymmetry& ptsym, const std::vector<size_t>& r, const std::vector<size_t>& invR, const int nthreads
) const {
  profile_update("Start Interpolator::rip_recip method");
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  auto no = this->count_scalars_vectors_matrices();
  if (!std::any_of(no.begin()+1, no.end(), [](ind_t n){return n>0;}))
    return false;
  T tmp_v[3], tmp_m[9];
  const ind_t b_{this->branches()}, s_{this->branch_span()};
  std::array<int,9> ident = {1,0,0, 0,1,0, 0,0,1};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
// otherwise gcc complains that const b_ is *already* shared
#pragma omp parallel for default(none) shared(x,ptsym,r,invR) private(tmp_v,tmp_m) firstprivate(ident,no,xsize) schedule(static)
#else
#pragma omp parallel for default(none) shared(b_,s_,x,ptsym,r,invR) private(tmp_v,tmp_m) firstprivate(ident,no,xsize) schedule(static)
#endif
  for (long long si=0; si<xsize; ++si){
    ind_t i = brille::utils::s2u<ind_t, long long>(si);
    T * xi = x.ptr(i);
    if (!brille::approx::matrix(3, ident.data(), ptsym.get(r[i]).data())){
      for (ind_t b=0; b<b_; ++b){
        // scalar elements do not need to be rotated, so skip them
        ind_t o = b*s_ + no[0];
        // we chose R such that Q = Rᵀq + τ.
        for (ind_t v=0; v<no[1]; ++v){
          brille::utils::mul_mat_vec(tmp_v, 3u, transpose(ptsym.get(r[i])).data(), xi+o);
          for (int j=0; j<3; ++j) x.val(i,o+j) = tmp_v[j];
          o += 3u; // shift 3 for each vector
        }
        for (ind_t m=0; m<no[2]; ++m){
          // Calculate R*M*R⁻¹ in two steps
          // first calculate M*R⁻¹, storing in tmp_m
          brille::utils::mul_mat_mat(tmp_m, 3u, xi+o, transpose(ptsym.get(invR[i])).data());
          // next calculate R*tmp_m, storing back in the x array
          brille::utils::mul_mat_mat(xi+o, 3u, transpose(ptsym.get(r[i])).data(), tmp_m);
          o += 9u; // shift 9 for each matrix
        }
      }
    }
  }
  profile_update("  End Interpolator::rip_recip method");
  return true;
}
