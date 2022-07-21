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
/*! \file
    \author Greg Tucker
    \brief Classes for the effect of symmetry on phonon eigenvectors
*/
#include "array_.hpp" // defines bArray
#include "array_l_.hpp"
//#include "lattice_dual.hpp"

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

/*! \brief A convenient store for
           \f$\left(F_0(k,R), \mathbf{r}_k - R \mathbf{r}_l\right)\f$ pairs

The application of a symmetry operation to a crystal strcture leaves its total
energy unchanged. The same is true for the total potential energy of the crystal.
The force constant matrix, howerever, which is the combination of all forces
acting on atom \f$k\f$ for a displacement of atom \f$l\f$, is not invariant but
instead the atom labels are permuted and individual-pairs are transformed
as tensors.

The dynamical matrix is the mass noramalised form of the Fourier transform of the
force constant matrix. In addition to the atom label permutation application of
a symmetry operation to the dynamical matrix appends a per-equivalent-atom-pair
phase factor due to the possibility that the operation moves equivalent atoms
in one unit cell into different unit cells.

Phonons are the solutions to the equations of motion described by the dynamical
matrix. The solutions make use of the dynamical matrix eigenvectors_ which also
acquire a per-atom phase factor upon operation of a symmetry element.
The phase factor is
@f[
e^{i \mathbf{Q}\cdot\left(\mathbf{r}_k - R \mathbf{r}_l\right)}
= e^{i R \mathbf{q} \cdot \left(\mathbf{r}_l - R \mathbf{r}_k\right)}
= e^{i \mathbf{q} \cdot \left( R^{-1} \mathbf{r}_l - \mathbf{r}_k\right)}
@f]
where \f$l\f$ is the equivalent atom which \f$R\f$ transforms \f$k\f$ to
\f$l = F_0(k,R)\f$

brille uses only the point group operations to find equivalent
\f$\mathbf{q}_i\f$ to \f$\mathbf{q} = \mathbf{Q} - \mathbf{\tau}\f$.
A pointgroup never contains more than \f$48\f$ elements, and for each
the permutation and vector are unique.
This class is designed to hold the \f$N_p \times N_a\f$ values of \f$F_0(k,R)\f$
and \f$\left(R^{-1}\mathbf{r}_k - \mathbf{r}_l\right)\f$ where there are
\f$N_p\f$ pointgroup operations and \f$N_a\f$ atoms in the Basis.


Rather than storing a 2D array of \f$F_0(k,R)\f$ and a 2D array of
\f$\left(R^{-1}\mathbf{r}_k - \mathbf{r}_l\right)\f$,
the standard template library does the heavy lifting in GammaTable.
If \f$ k \in (0, N_a-1)\f$ and \f$r \in (0,N_p-1)\f$ index a `std::map`
with keys that are the pairs \f$(k,r)\f$ then the mapped values are
\f$\left(N_0(k,R_r), v\right)\f$ where \f$v\f$ is an index into an
LDVec of all \f$\left(R^{-1} \mathbf{r}_l - \mathbf{r}_k\right)\f$.


The LDVec *could* contain only unique
\f$\left(R^{-1} \mathbf{r}_l - \mathbf{r}_k\right)\f$
but this is left for a future exercise if it becomes necessary.

*/
class GammaTable: public RotateTable {
 public:
  using lattice_t = lattice::Lattice<double>;
//   using ind_t = unsigned;
private:
  std::vector<ind_t> point2space_; //! maps Rᵣ to Sᵣ
  // use std::vectors instead of std::map for the mapping since we know the
  // total number of keys and how to calculate their positions in the vector
  ind_t n_atoms;
  ind_t n_sym_ops;
  std::vector<ind_t> l_mapping; //! maps (κ,r) to l=Nₒ(κ,Sᵣ)
  std::vector<ind_t> v_mapping; //! maps (κ,r) to v
  lattice_t lattice_;
  bArray<double> vectors_; //! element v is (Rᵣ⁻¹ ⃗rₖ - ⃗rₗ)
public:
  GammaTable(bool init, const lattice_t& dlat, const int time_reversal=0, double e_tol=0., int n_tol=1): lattice_(dlat) {
    if (init){
      this->construct(dlat, time_reversal, e_tol, n_tol);
    } else {
      l_mapping.resize(0);
      v_mapping.resize(0);
    }
  }
  bool construct(const lattice_t& dlat, const int time_reversal=0, double e_tol=0., int n_tol=1){
    lattice_ = dlat;
    auto ps = dlat.pointgroup_symmetry();
    auto spgsym = dlat.spacegroup_symmetry();
    if (time_reversal){
      ps = ps.add_space_inversion();
      spgsym = spgsym.add_space_inversion();
    }
    Basis bs = dlat.basis();
    // resize all vectors/arrays
    n_atoms = static_cast<ind_t>(bs.size());
    n_sym_ops = static_cast<ind_t>(ps.size());
    point2space_.resize(n_sym_ops);
    l_mapping.resize(n_atoms*n_sym_ops);
    v_mapping.resize(n_atoms*n_sym_ops);
    vectors_ = bArray<double>({n_atoms*n_sym_ops+1u, 3u}, 0.); // always put (0,0,0) first
    // construct a mapping of pointgroup indices to spacegroup indices
    // -- this mapping is likely not invertable, but it shouldn't (doesn't?)
    //    matter. I think.
    for (ind_t i=0; i<ps.size(); ++i){
      point2space_[i] = static_cast<ind_t>(spgsym.find_matrix_index(ps.get(i)));
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
      std::tie(found,l) = bs.equivalent_after_operation(k, motion, e_tol, n_tol);
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
      ind_t v = norm(vectors_.view(0,count) - vec).is(brille::cmp::eq, 0., e_tol, n_tol).first();
      // store the vector if it is not already present
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
  [[nodiscard]] const bArray<double>& vectors() const {return vectors_;}
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
  lattice::LVec<double> ldvector(Ik k, Ir r) const {
    return lattice::LVec<double>(LengthUnit::angstrom, lattice_, this->vector(k,r));
  }
  [[nodiscard]] const lattice_t& lattice() const {return lattice_;}
private:
  template<class Ik, class Ir> [[nodiscard]] ind_t calc_key(Ik k, Ir r) const {
    if (k<n_atoms && r<n_sym_ops)
      return static_cast<ind_t>(k)*n_sym_ops + static_cast<ind_t>(r);
    throw std::runtime_error("Attempting to access out of bounds mapping!");
  }
};

} // namespace brille
#endif
