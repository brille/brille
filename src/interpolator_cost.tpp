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
template<typename S>
void
Interpolator<T>::add_cost(const ind_t i0, const ind_t i1, std::vector<S>& cost, const bool arbitrary_phase_allowed) const {
  const T * x0 = data_.ptr(i0), * x1 = data_.ptr(i1);
  S s_cost{0}, v_cost{0}, m_cost{0};
  auto e_ = _elements;
  const ind_t s_{this->branch_span()}, b_{this->branches()}, mo_{e_[0]+e_[1]};
  if (arbitrary_phase_allowed){ // if the _vectorfun uses the Hermitian angle, e^iθ *never* matters.
    auto phased = std::unique_ptr<T[]>(new T[s_]);
    for (ind_t i=0; i<b_; ++i){
      const T * x0i = x0+i*s_;
      for (ind_t j=0; j<b_; ++j){
        brille::utils::inplace_antiphase(s_, x0i, x1+j*s_, phased.get());
        if (e_[0]) s_cost = this->_scalarfun(e_[0], x0i, phased.get());
        if (e_[1]) v_cost = this->_vectorfun(e_[1], x0i+e_[0], phased.get()+e_[0]);
        if (e_[2]){
          m_cost = 0;
          for (ind_t m=0; m<e_[2]/9; ++m)
            m_cost += brille::utils::frobenius_distance(3u, x0i+mo_+9u*m, phased.get()+mo_+9u*m);
        }
        cost[i*b_+j] += _costmult[0]*s_cost + _costmult[1]*v_cost + _costmult[2]*m_cost;
      }
    }
  } else {
    for (ind_t i=0; i<b_; ++i){
      const T * x0i = x0+i*s_;
      for (ind_t j=0; j<b_; ++j){
        const T * x1j = x1+j*s_;
        if (e_[0]) s_cost = this->_scalarfun(e_[0], x0i, x1j);
        if (e_[1]) v_cost = this->_vectorfun(e_[1], x0i+e_[0], x1j+e_[0]);
        if (e_[2]){
          m_cost = 0;
          for (ind_t m=0; m<e_[2]/9; ++m)
            m_cost += brille::utils::frobenius_distance(3u, x0i+mo_+9u*m, x1j+mo_+9u*m);
        }
        cost[i*b_+j] += _costmult[0]*s_cost + _costmult[1]*v_cost + _costmult[2]*m_cost;
      }
    }
  }
}
