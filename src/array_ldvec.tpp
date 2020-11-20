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
template<class T> template<typename... A>
LDVec<T>
LDVec<T>::view(A... args) const {
  return LDVec<T>(this->get_lattice(), this->bArray<T>::view(args...));
}
template<class T> template<typename... A>
LDVec<T>
LDVec<T>::extract(A... args) const {
  return LDVec<T>(this->get_lattice(), this->bArray<T>::extract(args...));
}
template<class T>
bArray<T>
LDVec<T>::get_hkl() const {
  return bArray<T>(*this); // strip off the Lattice information
}
template<class T>
bArray<double>
LDVec<T>::get_xyz() const {
  assert(this->is_row_ordered() && this->is_contiguous() && this->size(this->ndim()-1) == 3);
  bArray<double> xyz(this->shape(), this->stride());
  std::vector<double> toxyz = this->get_lattice().get_xyz_transform();
  auto tshape = this->shape();
  tshape.back() = 0u;
  for (auto x: this->subItr(tshape))
    brille::utils::multiply_matrix_vector(xyz.ptr(x), toxyz.data(), this->ptr(x));
  return xyz;
}
template<class T>
LQVec<double>
LDVec<T>::star() const {
  assert(this->is_row_ordered() && this->is_contiguous() && this->size(this->ndim()-1) == 3);
  std::vector<double> cvmt = this->get_lattice().get_covariant_metric_tensor();
  // shape_t ost(this->ndim(), 0);
  LQVec<double> slv(this->get_lattice().star(), this->shape(), this->stride());
  auto fx = this->shape(); fx.back() = 0;
  for (auto x: this->subItr(fx))
    brille::utils::multiply_matrix_vector(slv.ptr(x), cvmt.data(), this->ptr(x));
  slv /= 2.0*brille::pi; // ai= gij/2/pi * ai_star
  return slv;
}
template<class T>
LDVec<double>
LDVec<T>::cross(const size_t i, const size_t j) const {
  assert(this->is_row_ordered() && this->is_contiguous() && this->ndim()==2 && this->size(this->ndim()-1) == 3);
  bool bothok = (i<this->size(0) && j<this->size(0));
  LDVec<double> out(this->get_lattice(), bothok? 1u : 0u, 3u);
  if (bothok){
    Direct dlat = this->get_lattice();
    LQVec<double> lqv(dlat.star(), 1u);
    brille::utils::vector_cross<double,T,T,3>(lqv.ptr(0), this->ptr(i), this->ptr(j));
    lqv *= dlat.get_volume()/2.0/brille::pi;
    out =  lqv.star();
  }
  return out;
}

template<class T>
double
LDVec<T>::dot(const size_t i, const size_t j) const {
  assert(this->is_row_ordered() && this->is_contiguous() && this->ndim()==2 && this->size(this->ndim()-1) == 3);
  if (i>=this->size(0) || j>=this->size(0))
    throw std::out_of_range("attempted out of bounds access by dot");
  Direct lat = this->get_lattice();
  std::vector<double> len{lat.get_a(), lat.get_b(), lat.get_c()};
  std::vector<double> ang{lat.get_alpha(), lat.get_beta(), lat.get_gamma()};
  return same_lattice_dot(this->view(i),this->view(j),len,ang);
}

template<class T>
void
LDVec<T>::check_array(){
  auto last = this->ndim()-1;
  // the last dimension must cover 3-vectors
  if(this->size(last) != 3) throw std::runtime_error("LDVec objects must have a last dimension of size 3");
  // which must be contiguous in memory for most operations
  auto st = this->stride();
  if(st[last] != 1) throw std::runtime_error("LDVec objects must have a contiguous last dimension");
}
