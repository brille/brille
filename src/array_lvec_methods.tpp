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

namespace brille::lattice {

  template<class T>
  template<typename... A>
  LVec <T>
  LVec<T>::view(A... args) const {
    return LVec<T>(_type, _lattice, this->bArray<T>::view(args...));
  }

  template<class T>
  template<typename... A>
  LVec <T>
  LVec<T>::extract(A... args) const {
    return LVec<T>(_type, _lattice, this->bArray<T>::extract(args...));
  }

  template<class T>
  bArray <T>
  LVec<T>::hkl() const {
    return bArray<T>(*this); // strip off the Lattice information
  }

  template<class T>
  bArray<double>
  LVec<T>::xyz() const {
    assert(this->is_row_ordered() && this->is_contiguous() && this->size(this->ndim() - 1) == 3);
    bArray<double> xyz(this->shape(), this->stride());
    auto to_xyz = _lattice.to_xyz(_type);
    auto t_shape = this->shape();
    t_shape.back() = 0u;
    for (auto x: this->subItr(t_shape))
      brille::utils::multiply_matrix_vector(xyz.ptr(x), to_xyz.data(), this->ptr(x));
    return xyz;
  }

  template<class T>
  LVec<double>
  LVec<T>::star() const {
    assert(this->is_row_ordered() && this->is_contiguous() && this->size(this->ndim() - 1) == 3);
    auto tensor = _lattice.metric(_type);
    auto lu = LengthUnit::angstrom == _type ? LengthUnit::inverse_angstrom : LengthUnit::angstrom;
    LVec<double> slv(lu, _lattice, this->shape(), this->stride());
    auto fx = this->shape();
    fx.back() = 0;
    for (auto x: this->subItr(fx))
      brille::utils::multiply_matrix_vector(slv.ptr(x), tensor.data(), this->ptr(x));
    slv /= math::two_pi; // ai= gij/2/pi * ai_star
    return slv;
  }

  template<class T>
  LVec<double>
  LVec<T>::cross(const ind_t i, const ind_t j) const {
    assert(this->is_row_ordered() && this->is_contiguous() && this->ndim() == 2 && this->size(this->ndim() - 1) == 3);
    bool both_ok = (i < this->size(0) && j < this->size(0));
    LVec<double> out(_type, _lattice, both_ok ? 1u : 0u, 3u);
    if (both_ok) {
      auto lu = LengthUnit::angstrom == _type ? LengthUnit::inverse_angstrom : LengthUnit::angstrom;
      LVec<double> ldv(lu, _lattice, 1u);
      brille::utils::vector_cross<double, T, T, 3>(ldv.ptr(0), this->ptr(i), this->ptr(j));
      ldv *= _lattice.volume(_type) / brille::math::two_pi;
      out = ldv.star();
    }
    return out;
  }

  template<class T>
  double
  LVec<T>::dot(const ind_t i, const ind_t j) const {
    assert(this->is_row_ordered() && this->is_contiguous() && this->ndim() == 2 && this->size(this->ndim() - 1) == 3);
    if (i >= this->size(0) || j >= this->size(0))
      throw std::out_of_range("attempted out of bounds access by dot");
    return same_lattice_dot(this->view(i), this->view(j), _lattice.metric(_type));
  }

  template<class T>
  void
  LVec<T>::check_array() {
    auto last = this->ndim() - 1;
    // the last dimension must cover 3-vectors
    if (this->size(last) != 3) throw std::runtime_error("LVec objects must have a last dimension of size 3");
    // which must be contiguous in memory for most operations
    auto st = this->stride();
    if (st[last] != 1) throw std::runtime_error("LVec objects must have a contiguous last dimension");
  }

}