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
#ifndef BRILLE_PHONON_HPP_
#define BRILLE_PHONON_HPP_
#include <array>
#include <tuple>
#include "array_latvec.hpp" // defines bArray
namespace brille {

/*! \brief A superclass for all rotation-required tabulated information

When information is interpolated in the irreducible part of reciprocal space
it may or may not be simple to apply symmetry operation(s) to 'rotate' the
information to another part of reciprocal space. In simple cases there is no
need to tabulate information, but in more complex cases such tabulation is
worthwhile. To simplify the calling structure of the implemented rotations it
will be useful to have a dummy empty table.
*/
class RotateTable{};

/*! \brief A convenient store for (F₀(k,R), ⃗rₖ - R ⃗rₗ) pairs

The application of a symmetry operation to a crystal strcture leaves its total
energy unchanged. The same is true for the total potential energy of the crystal.
The force constant matrix, howerever, which is the combination of all forces
acting on atom k for a displacement of atom l, is not invariant but instead
the atom labels are permuted and individual-pairs are transformed as tensors.

The dynamical matrix is the mass noramalised form of the Fourier transform of the
force constant matrix. In addition to the atom label permutation application of
a symmetry operation to the dynamical matrix appends a per-equivalent-atom-pair
phase factor due to the possibility that the operation moves equivalent atoms
in one unit cell into different unit cells.

Phonons are the solutions to the equations of motion described by the dynamical
matrix. The solutions make use of the dynamical matrix eigenvectors_ which also
acquire a per-atom phase factor upon operation of a symmetry element.
The phase factor is
    exp[i ⃗Q⋅( ⃗rₖ - R ⃗rₗ )] = exp[i R ⃗q ⋅ ( ⃗rₗ - R ⃗rₖ)]
                          = exp[i ⃗q ⋅ (R⁻¹ ⃗rₗ - ⃗rₖ)]
where l is the equivalent atom which R transforms k to
    l = F₀(k,R)

brille uses only the point group operations to find equivalent qᵢᵣ to q = Q - τ
A pointgroup never contains more than 48 elements, and for each the permutation
and vector are unique.
This class is designed to hold the Nₚ × Nₐ values of F₀(k,R) and (R⁻¹ ⃗rₖ - ⃗rₗ)
where there are Nₚ pointgroup operations and Nₐ atoms in the Basis.

Rather than storing a 2D array of F₀(k,R) and a 2D array of (R⁻¹ ⃗rₖ - ⃗rₗ),
let the standard template library do the heavy lifting.
If k ∈ {0, Nₐ-1} and r ∈ {0, Nₚ-1} index a std::map with keys that are
the pairs (k,r); and use the values to store (N₀(k,Rᵣ), v)
where v is an index into a LDVec of the (R⁻¹ ⃗rₗ - ⃗rₖ).

The LDVec *could* contain only unique (R⁻¹ ⃗rₗ - ⃗rₖ) values, but this is left for
a future exercise if it becomes necessary.
*/
class GammaTable: public RotateTable {
public:
  using ind_t = unsigned;
private:
  std::vector<ind_t> point2space_; //! maps Rᵣ to Sᵣ
  // use std::vectors instead of std::map for the mapping since we know the
  // total number of keys and how to calculate their positions in the vector
  ind_t n_atoms;
  ind_t n_sym_ops;
  std::vector<ind_t> l_mapping; //! maps (κ,r) to l=Nₒ(κ,Sᵣ)
  std::vector<ind_t> v_mapping; //! maps (κ,r) to v
  Direct lattice_;
  bArray<double> vectors_; //! element v is (Rᵣ⁻¹ ⃗rₖ - ⃗rₗ)
public:
  explicit GammaTable(): n_atoms(0), n_sym_ops(0) {
    l_mapping.resize(0);
    v_mapping.resize(0);
  }
  GammaTable(const Direct& dlat, const int time_reversal=0){
    this->construct(dlat, time_reversal);
  }
  bool construct(const Direct& dlat, const int time_reversal=0){
    lattice_ = dlat;
    PointSymmetry ps = dlat.get_pointgroup_symmetry(time_reversal);
    Symmetry spgsym = dlat.get_spacegroup_symmetry(time_reversal);
    Basis bs = dlat.get_basis();
    // resize all vectors/arrays
    n_atoms = bs.size();
    n_sym_ops = ps.size();
    point2space_.resize(n_sym_ops);
    l_mapping.resize(n_atoms*n_sym_ops);
    v_mapping.resize(n_atoms*n_sym_ops);
    vectors_ = bArray<double>({n_atoms*n_sym_ops+1u, 3u}, 0.); // always put (0,0,0) first
    // construct a mapping of pointgroup indices to spacegroup indices
    // -- this mapping is likely not invertable, but it shouldn't (doesn't?)
    //    matter. I think.
    for (ind_t i=0; i<ps.size(); ++i){
      point2space_[i] = spgsym.find_matrix_index(ps.get(i));
      if (point2space_[i]>=spgsym.size()){
        info_update("The point group operation\n",ps.get(i),"was not found in the spacegroup!");
        throw std::runtime_error("Something has gone wrong with the correspondence of spacegroup to pointgroup");
      }
    }
    ind_t count{1u}; // for (0,0,0)
    // fill in the mappings
    for (ind_t k=0; k<bs.size(); ++k) for (ind_t r=0; r<ps.size(); ++r){
      bool found;
      ind_t l;
      auto motion = spgsym.getm(point2space_[r]);
      std::tie(found,l) = bs.equivalent_after_operation(k, motion);
      if (!found){
        info_update(bs.to_string(),"\ndoes not have an equivalent atom to ",k," for symmetry operation\n",motion.getr(),"+",motion.gett());
        throw std::runtime_error("All atoms in the basis *must* be mapped to an equivalent atom for *all* symmetry operations");
      }
      // calculate (R⁻¹ ⃗rₗ - ⃗rₖ)
      std::array<double,3> vec = motion.inverse().move_point(bs.position(l));
      auto rk = bs.position(k);
      for (int i=0; i<3; ++i) vec[i] -= rk[i];
      // check if this vector is in vectors_
      // look for an equal vector within the first count vectors_ -- return its index, or count if none match
      // count >= 1, so view is fine:
      ind_t v = norm(vectors_.view(0,count) - vec).is(brille::cmp::eq, 0.).first();
      // store the vector if its not already present
      if (count == v) vectors_.set(count++, vec);
      ind_t key = this->calc_key(k, r);
      l_mapping[key] = l;
      v_mapping[key] = v;
    }
    vectors_.resize(count); // not really necessary memory copy?
    return true;
  }
  template<class Ik, class Ir>
  ind_t F0(Ik k, Ir r) const {
    return l_mapping[this->calc_key(k,r)];
  }
  const bArray<double>& vectors() const {return vectors_;}
  template<class Ik, class Ir>
  ind_t vector_index(Ik k, Ir r) const {
    return v_mapping[this->calc_key(k,r)];
  }
  template<class Ik, class Ir>
  bArray<double> vector(Ik k, Ir r) const {
    return vectors_.extract(this->vector_index(k,r));
    // aternatively
    // return vectors_.view(this->vector_index(k,r));
  }
  template<class Ik, class Ir>
  LDVec<double> ldvector(Ik k, Ir r) const {
    return LDVec<double>(lattice_, this->vector(k,r));
  }
  const Direct& lattice() const {return lattice_;}
private:
  template<class Ik, class Ir> ind_t calc_key(Ik k, Ir r) const {
    if (k<n_atoms && r<n_sym_ops)
      return static_cast<ind_t>(k)*n_sym_ops + static_cast<ind_t>(r);
    throw std::runtime_error("Attempting to access out of bounds mapping!");
  }
};

} // namespace brille
#endif
