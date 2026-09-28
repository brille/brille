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

#ifndef BRILLE_PRIMITIVE_HPP_
#define BRILLE_PRIMITIVE_HPP_
/*! \file
    \author Greg Tucker
    \brief Defines a class to hold transformation matrices between conventional
           and primitive lattices for seven centring types.
*/
#include "spg_database.hpp"
#include "linear_algebra.hpp"
namespace brille {
/*! \brief A class to hold transformation matricies and their inverse, with the matrix
stored in an object determined by a provided Bravais centring type.

For each centred real-space-filling lattice there is a transformation matrix P
which converts its basis vectors into those of an equivalent primitive
space-filling lattice via,
    (aₚ,bₚ,cₚ) = (aₛ,bₛ,cₛ)P

The components of a real-space vector expressed in units of the conventional
cell vectors transforms to the primitive cell expression by
    xₚ = P⁻¹xₛ
this can be shown by considering Aₛ = (aₛ,bₛ,cₛ) and Aₚ = (aₚ,bₚ,cₚ) then
    Aₚ = Aₛ P
and noting that the vector expressed in either lattice remains unchanged
    Aₚ xₚ = Aₛ xₛ
    (Aₛ P) xₚ = Aₛ xₛ
    P xₚ = xₛ
    xₚ = P⁻¹ xₛ

If we express the reciprocal lattice vectors as the columns of a matrix B, then
    Aₛ Bₛᵀ ≡ 2π 𝟙
from which, we can derive the transformation of the reciprocal lattice vectors
    Aₚ Bₚᵀ = Aₛ Bₛᵀ
    Aₛ P Bₚᵀ = Aₛ Bₛᵀ
    Bₚᵀ = P⁻¹ Bₛᵀ
    Bₚ = Bₛ (P⁻¹)ᵀ
and a reciprocal lattice vector qₛ can be expressed in the primitive reciprocal
lattice since Bₚ qₚ = Bₛ qₛ
    Bₛ (P⁻¹)ᵀ qₚ = Bₛ qₛ
    (Pᵀ)⁻¹ qₚ = qₛ
    qₚ = Pᵀ qₛ
where the commutability of the inverse and transpose, (Pᵀ)⁻¹ ≡ (P⁻¹)ᵀ, has been
utilized.

This class holds P and P⁻¹ for a given centering type.
*/
class PrimitiveTransform{
private:
  Bravais bravais;    //!< The Bravais enum value
public:
  //! Construct from the centring type of the conventional lattice
  explicit PrimitiveTransform(const Bravais c): bravais{c} {}
//  /*! \brief Return the transformation matrix P
//
//  P converts the conventional basis vectors into those of an equivalent
//  primitive space-filling lattice via
//
//      (aₚ,bₚ,cₚ) = (aₛ,bₛ,cₛ)P
//
//  for (x,y,z) a matrix with its columns the indicated vectors.
//
//  \returns a flattened-3x3 representation of the matrix
//  */
//  std::array<double,9> get_P(void) const {
//    switch (bravais){
//      case Bravais::_: throw std::runtime_error("Invalid Bravais centring");
//      case Bravais::I: return {-0.5, 0.5, 0.5,  0.5,-0.5, 0.5,  0.5, 0.5,-0.5};
//      case Bravais::F: return {  0 , 0.5, 0.5,  0.5,  0 , 0.5,  0.5, 0.5,  0 };
//      case Bravais::A: return {  1 ,  0 ,  0 ,   0 , 0.5,-0.5,   0 , 0.5, 0.5};
//      case Bravais::B: return { 0.5,  0 ,-0.5,   0 ,  1 ,  0 ,  0.5,  0 , 0.5};
//      case Bravais::C: return { 0.5, 0.5,  0 , -0.5, 0.5,  0 ,   0 ,  0 ,  1 };
//      case Bravais::R: return {2./3.,-1./3.,-1./3.,  1./3., 1./3.,-2./3.,  1./3., 1./3., 1./3.};
//      default: return {  1 ,  1 ,  1,    1 ,  1 ,  1 ,   1 ,  1 ,  1 };
//    }
//  }
  /*! \brief Return the transformation matrix P

  P converts the conventional basis vectors into those of an equivalent
  primitive space-filling lattice via

      (aₚ,bₚ,cₚ) = (aₛ,bₛ,cₛ)P

  for (x,y,z) a matrix with its columns the indicated vectors.
  This method returns six times the matrix since P can only contain the values
  0, ±1/3, ±1/2, ±2/3, 1 which are 0, ±2, ±3, ±4, 6 sixths, respectively.

  \returns a flattened-3x3 representation of the matrix            */
  [[nodiscard]] std::array<int,9> get_6P() const {
    switch (bravais){
    case Bravais::_: throw std::runtime_error("Invalid Bravais centring");
    case Bravais::I: return {-3, 3, 3,  3,-3, 3,  3, 3,-3};
    case Bravais::F: return { 0, 3, 3,  3, 0, 3,  3, 3, 0};
    case Bravais::A: return { 6, 0, 0,  0, 3,-3,  0, 3, 3};
    case Bravais::B: return { 3, 0,-3,  0, 6, 0,  3, 0, 3};
    case Bravais::C: return { 3, 3, 0, -3, 3, 0,  0, 0, 6};
    case Bravais::R: return { 4,-2,-2,  2, 2,-4,  2, 2, 2};
    default: return {6,6,6, 6,6,6, 6,6,6};
    }
  }
  /*! \brief Return the inverse of the transformation matrix P

  The components of a real-space vector expressed in units of the conventional
  cell vectors transforms to the primitive cell expression by

      xₚ = P⁻¹xₛ

  where the vectors xᵢ are column vectors.

  \returns a flattened-3x3 representation of the matrix
  */
  [[nodiscard]] std::array<int,9> get_invP() const {
    switch (bravais){
      case Bravais::_: throw std::runtime_error("Invalid Bravais centring");
      case Bravais::I: return { 0, 1, 1,  1, 0, 1,  1, 1, 0};
      case Bravais::F: return {-1, 1, 1,  1,-1, 1,  1, 1,-1};
      case Bravais::A: return { 1, 0, 0,  0, 1, 1,  0,-1, 1};
      case Bravais::B: return { 1, 0, 1,  0, 1, 0, -1, 0, 1};
      case Bravais::C: return { 1,-1, 0,  1, 1, 0,  0, 0, 1};
      case Bravais::R: return { 1, 0, 1, -1, 1, 1,  0,-1, 1};
      default: return { 1, 0, 0,  0, 1, 0,  0, 0, 1};
    }
  }
//  /*! \brief Return the transpose of the transformation matrix P
//
//  A vector expressed in the conventional reciprocal lattice units can be
//  expressed in the primitive reciprocal lattice by means of the transformation
//
//    qₚ = Pᵀ qₛ
//
//  where the qᵢ are column vectors.
//
//  \returns a flattened-3x3 representation of the matrix
//  */
//  std::array<double,9> get_Pt(void) const { return transpose(this->get_P()); }
  /*! \brief Return the transpose of the transformation matrix 6*P

  A vector expressed in the conventional reciprocal lattice units can be
  expressed in the primitive reciprocal lattice by means of the transformation

          qₚ = Pᵀ qₛ

  where the qᵢ are column vectors.

  \returns a flattened-3x3 representation of the matrix
                   */
  [[nodiscard]] std::array<int,9> get_6Pt() const { return transpose(this->get_6P()); }
  /*! \brief Return the transpose of the inverse of the transformation matrix P

  If the conventional reciprocal lattice basis vectors form the rows of a matrix
  Bₛ then the primitive reciprocal lattice basis vectors are related via P as

    Bₚ = Bₛ (P⁻¹)ᵀ

  \returns a flattened-3x3 representation of the matrix
  */
  [[nodiscard]] std::array<int,9> get_invPt() const { return transpose(this->get_invP());}
  void print() const{
    auto sixM = this->get_6P();
    auto invM = this->get_invP();
    for (int i=0; i<3; ++i){
      printf("%3s", i==1 ? "to" :"");
      for (int j=0; j<3; ++j) printf(" % 4.2f", sixM[i*3+j]/6.);
      printf("%5s", i==1 ? "from" : "");
      for (int j=0; j<3; ++j) printf(" % 2d",invM[i*3+j]);
      printf("\n");
    }
  }
  //! A check if the Matrices *should* do anything, if the centring *is not* Primitive
  [[nodiscard]] bool does_anything() const { return (bravais!=Bravais::P);}
  //! A check if the Matrices *shouldn't* do anything, if the centring *is* Primitive
  [[nodiscard]] bool does_nothing() const { return (bravais==Bravais::P);}
  //! Return a basic string representation of the object
  [[nodiscard]] std::string string_repr() const {
    std::string repr = "<" + bravais_string(bravais) + " PrimitiveTransform object>";
    return repr;
  }
};


/*! \brief Type information for templated transformation functions

Since the transformation matrices stored in the PrimitiveTransform type do not
have the same internal data types, templated functions benefit from having a
single definition of the internal types (especially in the case of changing one
or both in the future).
*/
struct PrimitiveTraits{
  using P = double;
  using sixP = int;
  using invP = int;
};
} // end namespace brille
#endif
