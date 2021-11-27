#ifndef BRILLE_ROTATES_HPP_
#define BRILLE_ROTATES_HPP_
/*! \file
    \author Greg Tucker
    \brief An enumeration to represent the effect of a lattice symmetry operation
*/
namespace brille {
/*! \brief Represents how a lattice symmetry operation effects a vector or tensor

The symmetry operations of a lattice contain matrices which are either rotations
or rotoinversions. Within `brille` these matrices are expressed in the relative
units of the real space lattice.

How a vector transforms under application of the symmetry operation then depends
on whehter it is a real space vector, a reciprocal space vector, or an axial
vector. Similarly real space and reciprocal space tensors transform differently
under application of a symmetry operation.

An extra enumerated value, `Gamma`, is used to indicate that a quantity
represents an eigenvector of the grand dynamical matrix, which undergoes a
permutation and has an additional complex phase applied when transformed by a
spacegroup symmetry operation.
*/
enum class RotatesLike {Real, Reciprocal, Axial, Gamma};
}

#endif
