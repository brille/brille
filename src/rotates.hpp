#ifndef BRILLE_ROTATES_HPP_
#define BRILLE_ROTATES_HPP_
/*! \file
    \author Greg Tucker
    \brief Enumeration describing how the value rotates in a lattice symmetry operation
*/
namespace brille {
/*! \brief Represents how the value rotates in a lattice symmetry operation

The symmetry operations of a lattice contain matrices which are either rotations
or rotoinversions. Within `brille` these matrices are expressed in the relative
units of the real space lattice.

How a tensor transforms under application of the symmetry operation then depends
on whether it is a vector or pseudovector (axial).

Similarly the units of the tensors affect the transformation under application
of a symmetry operation (e.g. whether they are real or reciprocal space, or
or are in Cartesian or fractional units). The units are represented with the
LengthUnit enumeration, which should be used alongside this RotatesLike enumeration.

An extra enumerated value, `Gamma`, is used here to indicate that a quantity
represents an eigenvector of the grand dynamical matrix, which undergoes a
permutation and has an additional complex phase applied when transformed by a
spacegroup symmetry operation.
*/
enum class RotatesLike {vector, pseudovector, Gamma};
}

#endif
