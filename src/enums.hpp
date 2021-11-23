#ifndef BRILLE_ENUMS_HPP_
#define BRILLE_ENUMS_HPP_

#include "rotates.hpp"
#include "bravais.hpp"

namespace brille {
/*! \brief The units of anglular quantities supplied to some Lattice methods

The value of a measured angle can be expressed in a variety of different units
which prevents unambiguous calculations involving the magnitude of angles if
their unit is not also specified.

Common units used to measure angles include the radian and the degree which,
for most crystallographic systems, do not have overlapping magnitude ranges.
The radian is a measurement of the arclength swept by an angle divided by the
radius of the arc and a full circle, having circumference 2πr, has an internal
angle of 2π radian.
The degree is defined as 1/360ᵗʰ of a circle's arc so that a circle has an
internal angle of 360°.
Most lattices have basis vectors with relative angles more than 2π° ≈ 6.283°
and less than π radian ≈ 3.141 radian. Therefore, if the magnitude of a lattice's
angles are known but not whether they are expressed in radian or degrees we can
infer the unit based on which section of the positive numberline the three
magnitudes fall in: [(0,3.14159): radian, (6.283, 180): degree]

This enum exists in case an ambiguous case is found. In such a case the unit
can be provided and no inferrence will be used.

An third angle unit is understood -- multiples of π radian -- 'pi'.
*/
enum class AngleUnit { not_provided, radian, degree, pi };

/*! \brief The units of length quantities supplied to some Lattice methods

*/
enum class LengthUnit { none, angstrom, inverse_angstrom};

/*! \brief An enumeration to differentiate betwee Node types

\see NullNode, CubeNode, PolyNode
*/
enum class NodeType {null, cube, poly};




}


#endif