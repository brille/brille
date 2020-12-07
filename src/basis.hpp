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
/*! \file
    \author Greg Tucker
    \brief The atoms and their positions within a lattice is the Basis, class definition.
*/
#ifndef BASIS_HPP
#define BASIS_HPP
// #include <vector>
// #include <array>
// #include <tuple>
#include "symmetry.hpp"
// #include "utilities.hpp"
// #include "approx.hpp"
#include "types.hpp"
namespace brille {

/*! \brief The atom basis of a lattice

The positions and types of atoms within the basis of a lattice
*/
class Basis{
public:
  using point = std::array<double,3>; //!< The container for an atom position
  // using index = unsigned long;
private:
  std::vector<point> positions_; //!< all atom positions
  std::vector<ind_t> types_;     //!< all atom types
public:
  //! Explicit empty constructor
  explicit Basis(){};
  //! Construct from atom positions only, all of which are of a unique type
  Basis(const std::vector<point>& pos): positions_(pos){
    types_.resize(positions_.size());
    std::iota(types_.begin(), types_.end(), 0u);
  }
  /*! \brief Construct from atom positions and types

  \param pos the postions of all atoms
  \param typ the atom types of all atoms
  \note if `typ` is not the same size as `pos` it is not used and instead all
        atom are assumed to be of unique type
  */
  Basis(const std::vector<point>& pos, const std::vector<ind_t>& typ): positions_(pos){
    if (pos.size() == typ.size()){
      types_ = typ;
    } else {
      types_.resize(positions_.size());
      std::iota(types_.begin(), types_.end(), 0u);
    }
  }
  //! Return the number of atoms in the basis
  size_t size() const {return this->positions_.size();}
  //! Return the atom positions in the basis
  std::vector<point> positions() const {return this->positions_; }
  //! Return the position of an indexed atom
  point position(const size_t i) const {
    if (i<positions_.size())
      return this->positions_[i];
    throw std::runtime_error("Provided index is out of bounds");
  }
  //! Return all atom types in the basis
  std::vector<ind_t> types() const {return types_;}

  /*! \brief Determine if an atom exists in the basis

  \param Kappa a position which may not be in the first unit cell
  \returns a tuple containing
            - whether an atom is in the basis at position  ⃗k' =  ⃗K + ⃗G
              for a lattice vector ⃗G
            - that atom's index if it exists, or the total number of atoms
              in the basis if it does not
  */
  std::tuple<bool, ind_t> equivalent_to(const point& Kappa) const {
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
    ind_t kp = static_cast<ind_t>(std::distance(positions_.begin(), kp_itr));
    return std::make_tuple(found, kp);
  }
  /*! \brief Determine the equivalent atom index after a Symmetry operation

  Apply the Motion to the position of atom with index `k`
  then check if the resulting position has an equivalent index `k'`.

  \param k  the atom index to apply the operation to
  \param op the Spacegroup Symmetry operation
  \return the output of `equivalent_to` after applying the operation
  */
  template<class T, class R>
  std::tuple<bool, ind_t> equivalent_after_operation(const size_t k, const Motion<T,R>& op){
    if (k>=positions_.size()) throw std::runtime_error("invalid atom position index");
    point K_pos = op.move_point(positions_[k]);
    return this->equivalent_to(K_pos);
  }
  /*! \brief Determine the equivalent atom index after a PointSymmetry operation

  Apply the matrix to the position of the atom with index `k`
  then check if the resulting position has an equivalent index `k'`.

  \param k  the atom index to apply the operation to
  \param op the Pointgroup PointSymmetry operation
  \return the output of `equivalent_to` after applying the operation
  */
  template<class T>
  std::tuple<bool, ind_t> equivalent_after_operation(const size_t k, const std::array<T,9>& op){
    if (k>=positions_.size()) throw std::runtime_error("invalid atom positon index");
    point K_pos;
    brille::utils::multiply_matrix_vector(K_pos.data(), op.data(), positions_[k].data());
    return this->equivalent_to(K_pos);
  }
  //! Return a string representation of the atom types and positions
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
