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
#ifndef TYPES_HPP
#define TYPES_HPP

#include <vector>

// the number of nodes we might hold determines what type we need to store
// their indices:
//    bytes   maximum number    type
//      1                255    unsigned short
//      2             65,535    unsigned int
//      4      4,294,967,295    unsigned long
//      8            18×10¹⁸    unsigned long long (aka size_t)
// 65k is not enough. 4B *should* always be sufficient -- each node would
// occupy a fractional volume of ~2×10⁻¹⁰ of the polyhedron, which is overkill.
namespace brille {
  using ind_t = unsigned;
  using shape_t = std::vector<ind_t>;
}

#endif
