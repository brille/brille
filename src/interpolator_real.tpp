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
bool Interpolator<T>::rip_real(
  bArray<T>& x, const PointSymmetry& ptsym, const std::vector<size_t>& r, const std::vector<size_t>& invR, const int nthreads
) const {
  profile_update("Start Interpolator::rip_real method");
  auto no = this->count_scalars_vectors_matrices();
  if (!std::any_of(no.begin()+1, no.end(), [](ind_t n){return n>0;}))
    return false;

  const auto pool = ThreadPool::getInstance();
  if (nthreads) pool->resize(nthreads); else pool->resize();
  const auto workers = pool->size();
  auto task = [&](const size_t worker) {
    auto [f, l] = thread_slice(x.size(0), workers, worker);
    return [&,first=f,last=l,b_=branches(),s_=branch_span()]() {
      const std::array<int,9> ident = {1,0,0, 0,1,0, 0,0,1};
      for (size_t i=first; i<last; ++i) {
        auto xi = x.ptr(i);
        if (!approx_float::matrix(3, ident.data(), ptsym.get(r[i]).data())){
          for (ind_t b=0; b<b_; ++b){
            T tmp_v[3];
            // scalar elements do not need to be rotated, so skip them
            ind_t o = b*s_ + no[0];
            // rotate real vectors: since Q = Rᵀq + τ → Rv
            for (ind_t v=0; v<no[1]; ++v){
              utils::mul_mat_vec(tmp_v, 3u, ptsym.get(r[i]).data(), xi+o);
              for (int j=0; j<3; ++j) x.val(i,o+j) = tmp_v[j];
              o += 3u; // shift 3 for each vector
            }
            for (ind_t m=0; m<no[2]; ++m){
              T tmp_m[9];
              // Calculate R*M*R⁻¹ in two steps
              // first calculate M*R⁻¹, storing in tmp_m
              utils::mul_mat_mat(tmp_m, 3u, xi+o, ptsym.get(invR[i]).data());
              // next calculate R*tmp_m, storing back in the x array
              utils::mul_mat_mat(xi+o, 3u, ptsym.get(r[i]).data(), tmp_m);
              o += 9u; // shift 9 for each matrix
            }
          }
        }
      }
    };
  };
  for (size_t thread=0; thread<workers; ++thread) pool->enqueue(task(thread));
  pool->wait();

  profile_update("  End Interpolator::rip_real method");
  return true;
}
