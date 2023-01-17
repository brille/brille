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
#include <utility>

#include "symmetry.hpp"
#include "types.hpp"
#include "approx_float.hpp"

namespace brille {

/*! \brief The atom basis of a lattice

The positions and types of atoms within the basis of a lattice
*/
class Basis{
public:
  using element_t = double;
  using point = std::array<element_t,3>; //!< The container for an atom position
  // using index = unsigned long;
private:
  std::vector<point> positions_; //!< all atom positions
  std::vector<ind_t> types_;     //!< all atom types
public:
  //! Explicit empty constructor
  explicit Basis()= default;
  //! Construct from atom positions only, all of which are of a unique type
  explicit Basis(std::vector<point> pos): positions_(std::move(pos)){
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
  [[nodiscard]] size_t size() const {return this->positions_.size();}
  //! Return the atom positions in the basis
  [[nodiscard]] std::vector<point> positions() const {return this->positions_; }
  //! Return the position of an indexed atom
  [[nodiscard]] point position(const size_t i) const {
    if (i<positions_.size())
      return this->positions_[i];
    throw std::runtime_error("Provided index is out of bounds");
  }
  //! Return all atom types in the basis
  [[nodiscard]] std::vector<ind_t> types() const {return types_;}

  /*! \brief Determine the closest atom in the basis to a given point

  \param Kappa a position which may not be in the first unit cell
  \returns a tuple containing
            - whether an atom is in the basis at position  ⃗k' =  ⃗K + ⃗G
              for a lattice vector ⃗G
            - that atom's index if it exists, or the closest atom's index if not
            - the difference vector between ⃗k and  ⃗K + ⃗G
  */
  [[nodiscard]] std::tuple<bool, ind_t, point> closest_to(const ind_t type, const point& Kappa, element_t e_tol = element_t(0), int n_tol = 0) const {
    // find κ' equivalent to K in the first unit cell, all elements ∈ [0,1)
    // We need to protect against mapping Kᵢ ≈ -0 to 1; which you might introduce
    // by looking for an equivalent Κ' = Κ%1. This discontinuity near 0 and 1
    // is exactly what we don't want. Instead, map numbers near 0 and near 1 to
    // near zero before checking for point equivalency
    auto check = [type,Kappa,e_tol,n_tol](const point& p, const ind_t& t){
      point d{{0,0,0}}, z{{0,0,0}};
      if (t != type) return std::make_tuple(false, false, z);
      // find the difference vector % 1, with the discontinuity moved to 0.5
      for (int i=0; i<3; ++i) d[i] = Kappa[i]-p[i]+0.5;
      for (int i=0; i<3; ++i) d[i] = std::abs(d[i]-std::floor(d[i]))-0.5;
      auto ok = approx_float::equal(d, z, e_tol, e_tol, n_tol);
      //      info_update_if(ok,  "   match ", Kappa, " to ", p);
      //      info_update_if(!ok, "no match ", Kappa, " to ", p, " since ", d);
      return std::make_tuple(true, ok, d);
    };
    // now search for κ'
    std::vector<std::tuple<bool, bool, point>> v;
    v.reserve(positions_.size());
    for (size_t i=0; i<size(); ++i) v.push_back(check(positions_[i], types_[i]));
    // look for one with matching type and position
    auto itr = std::find_if(v.begin(), v.end(), [](const auto& a){return std::get<0>(a) && std::get<1>(a);});
    bool found = itr != v.end();
    // κ'
    ind_t kp;
    if (found) {
      kp = static_cast<ind_t>(std::distance(v.begin(), itr));
    } else {
      // no match within tolerance, so we must find the closest κ'
      std::vector<double> distances;
      distances.reserve(v.size());
      for (const auto& ed: v) if (std::get<0>(ed)) {
        // if the type matches, append the distance to the list
        double distance{0.};
        for (const auto & x: std::get<point>(ed)) distance += x * x;
        distances.push_back(distance);
      } else {
        // otherwise the distance is, effectively, infinite
        distances.push_back((std::numeric_limits<double>::max)());
      }
      auto min_at = std::min_element(distances.begin(), distances.end());
      kp = static_cast<ind_t>(std::distance(distances.begin(), min_at));
    }
    return std::make_tuple(found, kp, std::get<point>(v[kp]));
  }

  /*! \brief Determine if an atom exists in the basis

   \param type the atom type searched for
  \param Kappa a position which may not be in the first unit cell
   \param e_tol options floating point approximate comparison tolerance
   \param n_tol optional, ~10^[number of floating point digits] to compare
  \returns a tuple containing
            - whether an atom is in the basis at position  ⃗k' =  ⃗K + ⃗G
              for a lattice vector ⃗G
            - that atom's index if it exists, or the total number of atoms
              in the basis if it does not
  */
  [[nodiscard]] std::tuple<bool, ind_t> equivalent_to(const ind_t type, const point& Kappa, element_t e_tol = element_t(0), int n_tol = 0) const {
    // find κ' equivalent to K in the first unit cell, all elements ∈ [0,1)
    // We need to protect against mapping Kᵢ ≈ -0 to 1; which you might introduce
    // by looking for an equivalent Κ' = Κ%1. This discontinuity near 0 and 1
    // is exactly what we don't want. Instead map numbers near 0 and near 1 to
    // near zero before checking for point equivalency
    auto checker = [type,Kappa,e_tol,n_tol](const point& p, const ind_t& t){
      // if the atom type is different, the positions "can't" be the same
      if (t != type) return false;
      point d{{0,0,0}}, z{{0,0,0}};
      // find the difference vector % 1, with the discontinuity moved to 0.5
      for (int i=0; i<3; ++i) d[i] = Kappa[i]-p[i]+0.5;
      for (int i=0; i<3; ++i) d[i] = std::abs(d[i]-std::floor(d[i]))-0.5;
      auto ok = approx_float::equal(d, z, e_tol, e_tol, n_tol);
//      info_update_if(ok,  "   match ", Kappa, " to ", p);
//      info_update_if(!ok, "no match ", Kappa, " to ", p, " since ", d);
      return ok;
    };
    // now search for κ'
    for (size_t i=0; i<positions_.size(); ++i)
      if (checker(positions_[i], types_[i]))
        return std::make_tuple(true, static_cast<ind_t>(i));
    return std::make_tuple(false, static_cast<ind_t>(positions_.size()));
//    auto kp_itr = std::find_if(positions_.begin(), positions_.end(), checker);
//    bool found = kp_itr != positions_.end();
//    // κ' (or size(positions_) if no equivalent position found)
//    auto kp = static_cast<ind_t>(std::distance(positions_.begin(), kp_itr));
//    return std::make_tuple(found, kp);
  }
  /*! \brief Determine the equivalent atom index after a Symmetry operation

  Apply the Motion to the position of atom with index `k`
  then check if the resulting position has an equivalent index `k'`.

  \param k  the atom index to apply the operation to
  \param op the Spacegroup Symmetry operation
  \return the output of `equivalent_to` after applying the operation
  */
  template<class T, class R>
  std::tuple<bool, ind_t> equivalent_after_operation(const size_t k, const Motion<T,R>& op, element_t e_tol = element_t(0), int n_tol = 0){
    if (k>=positions_.size()) throw std::runtime_error("invalid atom position index");
    point K_pos = op.move_point(positions_[k]);
    return this->equivalent_to(types_[k], K_pos, e_tol, n_tol);
  }
  template<class T, class R>
  bool snap_to(const std::vector<Motion<T,R>> ops, element_t e_tol = element_t(0), int n_tol = 0){
    auto all_ok = [](const auto & a){
      return std::all_of(a.begin(), a.end(), [](const auto& b){return b;});
    };
    // every atom must have an equivalent atom for every operation
    std::vector<bool> is_ok(positions_.size(), true);
    for (size_t k=0; k<positions_.size(); ++k){
      for (const auto & op: ops){
        is_ok[k] = std::get<bool>(equivalent_after_operation(k, op, e_tol, n_tol));
        if (!is_ok[k]) break;
      }
      // if (!is_ok[k]) break; // we can not break early in case multiple atoms need to be modified
    }
    if (all_ok(is_ok)) return true;
    // not all ok; figure out *which* need to be snapped
    for (size_t k=0; k<positions_.size(); ++k){
      if (!is_ok[k]){
        // we know the k-th atom position is (slightly) wrong, find how much
        // it needs to be moved to be 'right':
        point delta{{0,0,0}};
        ind_t count{0};
        for (const auto & op: ops){
          point K_pos = op.move_point(positions_[k]);
          auto c = closest_to(types_[k], K_pos, e_tol, n_tol); // (found_or_not, closest_atom_index, difference_vector)
          if (std::get<ind_t>(c) == k){
            count++;
            for (size_t xi=0; xi<3; ++xi) delta[xi] += std::get<point>(c)[xi];
          }
        }
        // calculate the average difference vector for the closest atom
        if (count) for (auto & d: delta) d /= static_cast<double>(count);
        // and add it to the k-th position, thereby fixing its wrong position
        for (size_t xi=0; xi<3; ++xi) positions_[k][xi] += delta[xi];
      }
    }
    // double check that we've snapped appropriately
    for (size_t k=0; k<positions_.size(); ++k)
      for (const auto & op: ops)
        if (!std::get<bool>(equivalent_after_operation(k, op, e_tol, n_tol)))
          return false;
    return true;
  }
  /*! \brief Determine the equivalent atom index after a PointSymmetry operation

  Apply the matrix to the position of the atom with index `k`
  then check if the resulting position has an equivalent index `k'`.

  \param k  the atom index to apply the operation to
  \param op the Pointgroup PointSymmetry operation
  \return the output of `equivalent_to` after applying the operation
  */
  template<class T>
  std::tuple<bool, ind_t> equivalent_after_operation(const size_t k, const std::array<T,9>& op, element_t e_tol = element_t(0), int n_tol = 0){
    if (k>=positions_.size()) throw std::runtime_error("invalid atom positon index");
    point K_pos{{0,0,0}};
    brille::utils::multiply_matrix_vector(K_pos.data(), op.data(), positions_[k].data());
    return this->equivalent_to(types_[k], K_pos, e_tol, n_tol);
  }
  //! Return a string representation of the atom types and positions
  [[nodiscard]] std::string to_string() const {
    std::string repr;
    for (size_t i=0; i<this->size(); ++i)
      repr += std::to_string(types_[i]) + " : " + "( "
            + my_to_string(positions_[i][0]) + " "
            + my_to_string(positions_[i][1]) + " "
            + my_to_string(positions_[i][2])  + " )\n";
    return repr;
  }
#ifdef USE_HIGHFIVE
    template<class HF>
    std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, bool>
    to_hdf(HF& obj, const std::string& entry) const{
        auto group = overwrite_group(obj, entry);
        group.createAttribute("size", size());
        if (size()){
          // HighFive can't handle std::vector<std::array<T,3>> ?!?
          std::vector<std::vector<double>> p;
          for (const auto & pos: positions_){
            std::vector<double> in;
            for (const auto& pin: pos) in.push_back(pin);
            p.push_back(in);
          }
          group.createDataSet("positions", p);
          group.createDataSet("types", types_);
        }
        return true;
    }
    // Input from HDF5 file/object
    template<class HF>
    static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, Basis>
    from_hdf(HF& obj, const std::string& entry){
        auto group = obj.getGroup(entry);
        std::vector<std::array<double,3>> p;
        std::vector<ind_t> t;
        size_t num{0};
        group.getAttribute("size").read(num);
        if (num) {
          auto pos = group.getDataSet("positions");
          std::vector<size_t> pos_shape = pos.getDimensions();
          if (pos_shape.size() != 2u || pos_shape[1] != 3u || pos_shape[0] != num) {
            throw std::runtime_error("Position should be (N,3) in shape!");
          }
          // or can we only read a std::vector<std::vector>?
          auto tot = pos.getElementCount();
          auto *pos_data = new double[tot]();
          pos.read(pos_data);
          for (size_t i = 0; i < tot; i += 3u) {
            std::array<double, 3> x{
                {pos_data[i], pos_data[i + 1], pos_data[i + 2]}};
            p.push_back(x);
          }
          delete[] pos_data;

          group.getDataSet("types").read(t);
        }
        return {p,t};
    }
#endif //USE_HIGHFIVE
};

} // end namespace brille
#endif // BASIS_HPP
