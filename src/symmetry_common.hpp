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
#ifndef COMMON_SYMMETRY_HPP
#define COMMON_SYMMETRY_HPP
/*! \file
    \author Greg Tucker
    \brief Utility definitions for PointSymmetry and Symmetry
*/
// #include <array>
// #include <vector>
#include <algorithm>
#include "utilities.hpp"

namespace brille {
//! A 3x3 matrix in flattened-form
template<class T> using Matrix = std::array<T,9>;
//! A 3 vector
template<class T> using Vector = std::array<T,3>;
//! A list of 3x3 flattened matrices
template<class T> using Matrices = std::vector<Matrix<T>>;
//! A list of 3 vectors
template<class T> using Vectors = std::vector<Vector<T>>;

} //end namespace brille
#endif // COMMON_SYMMETRY_HPP
