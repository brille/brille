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

// Perform sanity checks before attempting to interpolate
template<class T, class S> template<typename R>
unsigned int
Mesh3<T,S>::check_before_interpolating(const bArray<R>& x) const{
  unsigned int mask = 0u;
  if (data_.size()==0)
    throw std::runtime_error("The interpolation data must be filled before interpolating.");
  if (x.ndim()!=2 || x.size(1)!=3u)
    throw std::runtime_error("Only (n,3) two-dimensional Q vectors supported in interpolating.");
  if (x.stride().back()!=1)
    throw std::runtime_error("Contiguous vectors required for interpolation.");
  return mask;
}
//! Perform linear interpolating at the specified points in the mesh's orthonormal frame
template<class T, class S> template<typename R>
std::tuple<brille::Array<T>,brille::Array<S>>
Mesh3<T,S>::interpolate_at(const bArray<R>& x) const {
  this->check_before_interpolating(x);
  auto valsh = data_.values().shape();
  auto vecsh = data_.vectors().shape();
  valsh[0] = x.size(0);
  vecsh[0] = x.size(0);
  brille::Array<T> vals(valsh);
  brille::Array<S> vecs(vecsh);
  // vals and vecs are row-ordered contiguous by default, so we can create
  // mutable data-sharing Array2 objects for use with
  // Interpolator2::interpolate_at through the constructor:
  brille::Array2<T> vals2(vals);
  brille::Array2<S> vecs2(vecs);
  for (ind_t i=0; i<x.size(0); ++i){
    verbose_update("Locating ",x.to_string(i));
    auto verts_weights = this->mesh.locate(x.view(i));
    if (verts_weights.size()<1){
      debug_update("Point ",x.to_string(i)," not found in tetrahedra!");
      throw std::runtime_error("Point not found in tetrahedral mesh");
    }
    data_.interpolate_at(verts_weights, vals2, vecs2, i);
  }
  return std::make_tuple(vals, vecs);
}
template<class T, class S> template<typename R>
std::tuple<brille::Array<T>,brille::Array<S>>
Mesh3<T,S>::parallel_interpolate_at(const bArray<R>& x, const int threads) const {
  omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
  this->check_before_interpolating(x);
  // not used in parallel region
  auto valsh = data_.values().shape();
  auto vecsh = data_.vectors().shape();
  valsh[0] = x.size(0);
  vecsh[0] = x.size(0);
  // shared between threads
  brille::Array<T> vals(valsh);
  brille::Array<S> vecs(vecsh);
  // vals and vecs are row-ordered contiguous by default, so we can create
  // mutable data-sharing Array2 objects for use with
  // Interpolator2::interpolate_at through the constructor:
  brille::Array2<T> vals2(vals);
  brille::Array2<S> vecs2(vecs);
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long xsize = brille::utils::u2s<long, ind_t>(x.size(0));
#pragma omp parallel for default(none) shared(x, vals2, vecs2, xsize) schedule(dynamic)
  for (long si=0; si<xsize; ++si){
    ind_t i = brille::utils::s2u<ind_t, long>(si);
    auto verts_weights = this->mesh.locate(x.view(i));
    data_.interpolate_at(verts_weights, vals2, vecs2, i);
  }
  return std::make_tuple(vals, vecs);
}

// template<class T, class S> template<typename R>
// std::vector<ind_t>
// Mesh3<T,S>::which_neighbours(const std::vector<R>& t, const R value, const ind_t v) const{
//   std::vector<ind_t> out;
//   for (ind_t n: this->mesh.neighbours(v)) if (t[n] == value) out.push_back(n);
//   return out;
// }
