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

/*! \file */

#ifndef BRILLE_SYMMETRY_H_
#define BRILLE_SYMMETRY_H_
/*! \file
    \author Greg Tucker
    \brief Classes for a lattice spacegroup Symmetry operations
*/
#include "symmetry_common.hpp"
#include "motion.hpp"
#include "bravais.hpp"
namespace brille {


/*****************************************************************************\
| Symmetry class                                                              |
|-----------------------------------------------------------------------------|
|   Holds N 3x3 matrices R and N 3 vectors T which are the motions of a       |
|   space group. A motion is the combination of a rotation and a translation. |
|-----------------------------------------------------------------------------|
| Member functions:                                                           |
|   set, setrot, settran   copy a given matrix and/or vector into R and/or T  |
|                          for the motion i.                                  |
|   get, getrot, gettran   copy the rotation and/or translation of motion i   |
|                          into the provided matrix and/or vector.            |
|   getrot, gettran        return a pointer to the first element of the       |
|                          rotation or translation of motion i.               |
|   resize                 grow or shrink the number of motions that the      |
|                          object can hold -- this necessitates a memory copy |
|                          if the old and new sizes are non-zero.             |
|   size                   return the number of motions the object can store. |
\*****************************************************************************/
/*! \brief The spacegroup symmetry operations of a lattice

The spacegroup of a crystallographic lattice is comprised of a set of symmetry
operations, each of which is a Motion.
The operations form a group, that is for any two elements Mᵢ and Mⱼ ∈ G their
product Mᵢⱼ = MᵢMⱼ is also a member of the group, Mᵢⱼ ∈ G.

Since the group is complete and all Mᵢⱼ ∈ G we can reduce the Motions used to
describe a group to a subset which are able to *generate* the group through
repeated multiplication. Any such subset are the generators of the group and
these can be obtained from a Symmetry operation via a class method.

If only the generators are specified when a Symmetry object is constructed then
the full group can be obtained via a class method.

*/
class Symmetry{
public:
  using Motions = std::vector<Motion<int,double>>;
private:
  Motions M;
public:
  //! Empty constructor
  Symmetry(size_t n=0) { M.resize(n); }
  //! Construct from a list of symmetry operation Motions
  Symmetry(const Motions& m): M(m) {};
  //! Construct from a CIF xyz encoded string of one or more Motions
  Symmetry(const std::string& strrep){ this->from_ascii(strrep); };
  //! Return the centring type of the symmetry operations
  Bravais                getcentring() const;
  //! Return a reference to the held Motions
  const Motions&         getallm() const {return this->M;}
  //! Return just the matrix part of all held Motions
  Matrices<int>          getallr() const;
  //! Return just the translation part of all held motions
  Vectors<double>        getallt() const;
  //! Return the number of symmetry operations held
  size_t                 size()    const { return M.size(); }
  //! Increase or decrease the allocated space for Motions
  size_t                 resize(const size_t i)                                ;
  /*! \brief Add one Motion to the list of Motions from its matrix and vector parts

  \param r A pointer to the first element of a 3x3 row-ordered contiguous array
  \param t A pointer to the first element of a 3-element contiguous array
  */
  size_t                 add(const int *r, const double *t)                    ;
  /*! \brief Add one Motion to the list of Motions from its Matrix and Vector parts

  \param W A reference to the Matrix part of the Motion to add
  \param w A reference to the Vector part of the Motion to add
  */
  size_t                 add(const Matrix<int>& W, const Vector<double>& w)    ;
  /*! \brief Add one Motion to the list of Motions

  \param M A reference to the Motion to add
  */
  size_t                 add(const Motion<int,double>& M)                      ;
  /*! \brief Add one CIF xyz encoded Motion to the list of Motions

  \param motion A comma separated list of the combined matrix and vector parts
                of a Motion, possibly terminated by a semicolon.

  See, e.g., the IUCr definition of
  [_space_group_symop_operation_xyz](https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Ispace_group_symop_operation_xyz.html)
  */
  size_t                 add(const std::string& motion)                        ;
  //! Replace the held Motions by those described in a CIF xyz encoded string
  bool                   from_ascii(const std::string& s)                      ;
  //! Remove the `i`ᵗʰ element from the Motions
  size_t                 erase(const size_t i)                                 ;
  //! Return the `i`ᵗʰ element then remove it from the Motions
  Motion<int,double>     pop(const size_t i)                                   ;
  //! Check for the presence of a specified Motion in the Motions
  bool                   has(const Motion<int,double>&)                   const;
  //! Return the Matrix part of the `i`ᵗʰ Motion
  Matrix<int>            getr(const size_t i)                             const;
  //! Return the Vector partof the `i`ᵗʰ Motion
  Vector<double>         gett(const size_t i)                             const;
  //! Return the `i`ᵗʰ Motion
  Motion<int,double>     getm(const size_t i)                             const;
  /*! \brief Find and return the order of the `i`ᵗʰ Motion

  \param i The index of the Motion for which to find its order

  The order of a Motion, o, is given by Mᵒ ≡ (𝟙| ⃗0). This method finds the number of
  successive applications of the `i`ᵗʰ Motion required to equal the identity
  operation.
  */
  int order(const size_t i) const;
  //! Return the order of all Motions
  std::vector<int> orders(void) const;
  /*! \brief Generate all group operations from the known subset

  Assuming that the held Motions represet a subset of the full Symmetry group
  examine all pairwise products of the known Motions to find missing group
  members, adding the unique missing members to the known Motions.
  Repeat this until no new Motions are found.
  */
  Symmetry generate(void) const;
  /*! \brief Find a subset of Motions which can generate all stored Motions

  If the Symmetry object represents a complete group, it must contain Motions
  G = {M₁, M₂, …, Mₙ} which have the property for any i ∈ (1,n), j ∈ (1,n)
  MᵢMⱼ = Mₖ ∈ G.

  This algorithm traverses through G finding pairwise products Mₖ=MᵢMⱼ which are
  members of the group and eliminates the Mₖ.
  Once all pairwise product of the elements of the remaining Motions G' do not
  yield members of the subset, this subset is comprised soley of (one set of)
  generators of the group G.

  This returned Symmetry object contains the subgroup G'.
  */
  Symmetry generators(void) const;
  bool operator==(const Symmetry& other) const;
  bool operator!=(const Symmetry& other) const { return !this->operator==(other);}
  //! Check whether the Motions contain (̄𝟙, ⃗0)
  bool has_space_inversion() const;
  /*! \brief Find a Matrix in the Motions and return its index

  If none of the Motions have a matching Matrix part the total number of
  Motions is instead returned.
  */
  size_t  find_matrix_index(const Matrix<int>&) const;

#ifdef USE_HIGHFIVE
  // Output to HDF5 file/object
  template<class HF>
  std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, bool>
  to_hdf(HF& obj, const std::string& entry) const{
    std::vector<HF_Motion<int, double>> hfm;
    for (const auto & x: M) hfm.push_back(HF_Motion(x.getr(), x.gett()));
    auto ds = overwrite_data(obj, entry, hfm);
    return true;
  }
  // Input from HDF5 file/object
  template<class HF>
  static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, Symmetry>
  from_hdf(HF& obj, const std::string& entry){
    std::vector<HF_Motion<int, double>> hfm;
    obj.getDataSet(entry).read(hfm);
    Motions m;
    for(const auto& x: hfm){
      auto [W, w] = x.tuple();
      Motion<int, double> mot(W, w);
      m.push_back(mot);
    }
    return Symmetry(m);
  }
#endif
};

} // end namespace brille
#endif
