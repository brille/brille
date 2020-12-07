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

/* This file has evolved from pointgroup.h distributed as part of spglib.
   Changes have been made to introduce C++ style classes as well as other
   modifications to suit the needs of brille.
   spglib was licensed under the following BSD 3-clause license:

 Copyright (C) 2008 Atsushi Togo
 All rights reserved.

 This file is part of spglib. https://github.com/atztogo/spglib

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in
   the documentation and/or other materials provided with the
   distribution.

 * Neither the name of the phonopy project nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE. */

#ifndef BRILLE_POINTGROUP_HPP_
#define BRILLE_POINTGROUP_HPP_
/*! \file
    \author Atsushi Togo
    \author Greg Tucker
    \brief Lattice spacegroup information class and utility functions

    The Pointgroup class is a C++ rewrite of the C struct `Pointgroup` from
    [pointgroup.h](https://github.com/spglib/spglib/blob/develop/src/pointgroup.h)
    which is part of [spglib](https://github.com/atztogo/spglib).

*/
// #include <string>
#include "symmetry.hpp"
#include "pointsymmetry.hpp"
namespace brille {

/*! \brief The geometric crystal classes

There are seven geometric crystal classes with matching pointgroup and lattice

| holohedry | enumeration |
|-----------|-------------|
| ̄1         | triclinic   |
| 2/m       | monoclinic  |
| mmm       | orthogonal  |
| ̄3m        | trigonal    |
| 4/mmm     | tetragonal  |
| 6/mmm     | hexagonal   |
| m3m       | cubic       |

See, e.g., [IUCr Online Dictionary of Crystallography](https://dictionary.iucr.org/Holohedry)
*/
enum class Holohedry {_, triclinic, monoclinic, orthogonal, tetragonal, trigonal, hexagonal, cubic};
//! Return a string representation of the Holohendry enumerated value
std::string my_to_string(const Holohedry& h);
/*! \brief The eleven Laue classes

The eleven geometric crystal classes
| Laue class | enumeration |
|------------|-------------|
| ̄1          | _1          |
| 2/m        | _2m         |
| mmm        | _mmm        |
| ̄3          | _3          |
| ̄3m         | _3m         |
| 4/m        | _4m         |
| 4/mmm      | _4mmm       |
| 6/m        | _6m         |
| 6/mmm      | _6mmm       |
| m3         | _m3         |
| m3m        | _m3m        |

See, e.g., [IUCr Online Dictionary of Crystallography](https://dictionary.iucr.org/Laue_class)
*/
enum class Laue {_, _1, _2m, _mmm, _4m, _4mmm, _3, _3m, _6m, _6mmm, _m3, _m3m};
//! Return a string representation of the Laue enumerated value
std::string my_to_string(const Laue& l);

/*! \brief Information about a lattice pointgroup

There are 32 crystallographic pointgroups. The information held in this class
is initialized as part of the static namespace global `pointgroup_data`.
*/
class Pointgroup{
public:
  int number;              //!< A serial number
  std::string symbol;      //!< The Hermann-Mauguin notation symbol
  std::string schoenflies; //!< The Schoenflies notation symbol
  Holohedry holohedry;     //!< The geometric class of the pointgroup
  Laue laue;               //!< The Laue class of the pointgroup
  // Initializers:
  //! empty constructor
  Pointgroup(): number(0), holohedry(Holohedry::_), laue(Laue::_){}
  //! Construct from a 'serial number'
  Pointgroup(const int no): number(no) {setup();}
  //! Construct from all required information
  Pointgroup(const int no, const std::string& sym, const std::string& sch, const Holohedry& h, const Laue& l):
    number(no), symbol(sym), schoenflies(sch), holohedry(h), laue(l) {}
  //! Return a string representation of the contained information
  std::string to_string(void) const {
    std::string str = "<Pointgroup";
    str += " " + symbol;
    str += " " + schoenflies;
    str += " " + my_to_string(holohedry);
    str += " " + my_to_string(laue);
    return str+">";
  }
  //! Return the serial number
  int get_number(void) const {return number;}
  //! Return the Hermann-Mauguin symbol
  std::string get_symbol(void) const { return symbol; }
  //! Return the Schoenflies symbol
  std::string get_schoenflies(void) const { return schoenflies; }
  //! Return the Holohedry enumerated geometric class
  Holohedry get_holohedry(void) const { return holohedry; }
  //! Return the Laue enumerated class
  Laue get_laue(void) const { return laue; }
  //! Return a string representation of the geometric class
  std::string get_holohedry_string(void) const { return my_to_string(holohedry); }
  //! Return a string representation of the Laue class
  std::string get_laue_string(void) const { return my_to_string(laue); }
protected:
  void setup(void);
};
//! Set the transformation matrix from the set of pointgroup symmetry operations
Pointgroup ptg_get_transformation_matrix(int *transform_mat, const int *rotations, const int num_rotations);
//! Return a pointgroup PointSymmetry object from the provided symmetry operations
PointSymmetry ptg_get_pointsymmetry(const int *rotations, const int num_rotations);
//! Set the pointgroup symmetry operations from the 'Hall number'
int get_pointgroup_rotations_hall_number(int *rotations, const int max_size, const int hall_number, const int is_time_reversal);
//! Determine the isometry of a pointgroup operation in matrix form
int isometry_value(const int *rot);
//! Determine the order of a pointgroup operation in matrix form
int rotation_order(const int *rot);
//! Return three 3-vectors: the characteristic axis of a pointgroup operation and two perpendicular axes
std::array<std::array<int,3>,3> rotation_axis_and_perpendicular_vectors(const int* rot);
//! Remove duplicates from a list of pointgroup operations in matrix form
std::vector<std::array<int,9>> get_unique_rotations(const std::vector<std::array<int,9>>&, const int);
} // namespace brille
#endif
