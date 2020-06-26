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

// Perform sanity checks before attempting to interpolate
template<class T, class S> template<typename R>
unsigned int
Mesh3<T,S>::check_before_interpolating(const ArrayVector<R>& x) const{
  unsigned int mask = 0u;
  if (data_.size()==0)
    throw std::runtime_error("The mesh must be filled before interpolating!");
  if (x.numel()!=3u)
    throw std::runtime_error("Mesh3 requires x values which are three-vectors.");
  return mask;
}
//! Perform linear interpolating at the specified points in the mesh's orthonormal frame
template<class T, class S> template<typename R>
std::tuple<ArrayVector<T>,ArrayVector<S>>
Mesh3<T,S>::interpolate_at(const ArrayVector<R>& x) const {
  this->check_before_interpolating(x);
  ArrayVector<T> vals(data_.values().numel(), x.size());
  ArrayVector<S> vecs(data_.vectors().numel(), x.size());
  std::vector<size_t> vertices;
  std::vector<double> weights;
  size_t found_tet, max_valid_tet = this->mesh.number_of_tetrahedra()-1;
  for (size_t i=0; i<x.size(); ++i){
    verbose_update("Locating ",x.to_string(i));
    found_tet = this->mesh.locate(x.extract(i), vertices, weights);
    debug_update_if(found_tet > max_valid_tet,"Point ",x.to_string(i)," not found in tetrahedra!");
    if (found_tet > max_valid_tet)
      throw std::runtime_error("Point not found in tetrahedral mesh");
    verbose_update("Interpolate between vertices ", vertices," with weights ",weights);
    data_.interpolate_at(vertices, weights, vals, vecs, i);
  }
  return std::make_tuple(vals, vecs);
}
template<class T, class S> template<typename R>
std::tuple<ArrayVector<T>,ArrayVector<S>>
Mesh3<T,S>::parallel_interpolate_at(const ArrayVector<R>& x, const int threads) const {
  omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
  this->check_before_interpolating(x);
  // shared between threads
  ArrayVector<T> vals(data_.values().numel(), x.size());
  ArrayVector<S> vecs(data_.vectors().numel(), x.size());
  // private to each thread
  std::vector<size_t> indexes;
  std::vector<double> weights;
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long xsize = unsigned_to_signed<long, size_t>(x.size());
#pragma omp parallel for default(none) shared(x, vals, vecs, xsize) private(indexes, weights) schedule(dynamic)
  for (long si=0; si<xsize; ++si){
    size_t i = signed_to_unsigned<size_t, long>(si);
    this->mesh.locate(x.extract(i), indexes, weights);
    data_.interpolate_at(indexes, weights, vals, vecs, i);
  }
  return std::make_tuple(vals, vecs);
}

template<class T, class S> template<typename R>
std::vector<size_t>
Mesh3<T,S>::which_neighbours(const std::vector<R>& t, const R value, const size_t v) const{
  std::vector<size_t> out;
  for (size_t n: this->mesh.neighbours(v)) if (t[n] == value) out.push_back(n);
  return out;
}
