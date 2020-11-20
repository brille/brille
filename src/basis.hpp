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
#ifndef BASIS_HPP
#define BASIS_HPP
#include <vector>
#include <array>
#include <tuple>
#include "symmetry.hpp"
#include "utilities.hpp"
#include "approx.hpp"
namespace brille {

class Basis{
public:
  using point = std::array<double,3>;
  using index = unsigned long;
private:
  std::vector<point> positions_;
  std::vector<index> types_;
public:
  explicit Basis(){};
  Basis(const std::vector<point>& pos): positions_(pos){
    types_.resize(positions_.size());
    std::iota(types_.begin(), types_.end(), 0u);
  }
  Basis(const std::vector<point>& pos, const std::vector<index>& typ): positions_(pos){
    if (pos.size() == typ.size()){
      types_ = typ;
    } else {
      types_.resize(positions_.size());
      std::iota(types_.begin(), types_.end(), 0u);
    }
  }
  size_t size() const {return this->positions_.size();}
  std::vector<point> positions() const {return this->positions_; }
  point position(const size_t i) const {
    if (i<positions_.size())
      return this->positions_[i];
    throw std::runtime_error("Provided index is out of bounds");
  }
  std::vector<index> types() const {return types_;}
  // Check if there is an atom in the basis with position ⃗k' =  ⃗K + ⃗G
  // if so, return its index and ⃗G
  // otherwise, return the number of atoms and ⃗K - ⃗G (which is not an existant atom position)
  std::tuple<bool, size_t> equivalent_to(const point& Kappa) const {
    // find κ' equivalent to K in the first unit cell, all elements ∈ [0,1)
    // We need to protect against mapping Kᵢ ≈ -0 to 1; which you might introduce
    // by looking for an equivalent Κ' = Κ%1. This discontinuity near 0 and 1
    // is exactly what we don't want. Instead map numbers near 0 and near 1 to
    // near zero before checking for point equivalency
    auto checker = [Kappa](const point& p){
      point d, z{{0,0,0}};
      // find the difference vector % 1, with the discontinuity moved to 0.5
      for (int i=0; i<3; ++i) d[i] = Kappa[i]-p[i]+0.5;
      for (int i=0; i<3; ++i) d[i] = std::abs(d[i]-std::floor(d[i]))-0.5;
      return brille::approx::vector(d.data(), z.data());
    };
    // now search for κ'
    auto kp_itr = std::find_if(positions_.begin(), positions_.end(), checker);
    bool found = kp_itr != positions_.end();
    // κ' (or size(positions_) if no equivalent position found)
    size_t kp = std::distance(positions_.begin(), kp_itr);
    return std::make_tuple(found, kp);
  }
  // apply the Motion op (a Spacegroup symmetry operation) to atom position k
  // then check if the resulting position has an equivalent index k'
  template<class T, class R> std::tuple<bool, size_t> equivalent_after_operation(const size_t k, const Motion<T,R>& op){
    if (k>=positions_.size()) throw std::runtime_error("invalid atom position index");
    point K_pos = op.move_point(positions_[k]);
    return this->equivalent_to(K_pos);
  }
  // apply the matrix op (a Pointgroup symmetry operation) to atom position k
  // then check if the resulting position ahs an equivalent index k'
  template<class T> std::tuple<bool, size_t> equivalent_after_operation(const size_t k, const std::array<T,9>& op){
    if (k>=positions_.size()) throw std::runtime_error("invalid atom positon index");
    point K_pos;
    brille::utils::multiply_matrix_vector(K_pos.data(), op.data(), positions_[k].data());
    return this->equivalent_to(K_pos);
  }
  //
  std::string to_string() const {
    std::string repr;
    for (size_t i=0; i<this->size(); ++i)
      repr += std::to_string(types_[i]) + " : " + "( "
            + std::to_string(positions_[i][0]) + " "
            + std::to_string(positions_[i][1]) + " "
            + std::to_string(positions_[i][2])  + " )\n";
    return repr;
  }
};

} // end namespace brille
#endif // BASIS_HPP
