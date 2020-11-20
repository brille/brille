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
#include <string>
#include "symmetry.hpp"
#include "pointsymmetry.hpp"
namespace brille {

enum class Holohedry {_, triclinic, monoclinic, orthogonal, tetragonal, trigonal, hexagonal, cubic};
std::string my_to_string(const Holohedry& h);
enum class Laue {_, _1, _2m, _mmm, _4m, _4mmm, _3, _3m, _6m, _6mmm, _m3, _m3m};
std::string my_to_string(const Laue& l);

class Pointgroup{
public:
  int number;
  std::string symbol;
  std::string schoenflies;
  Holohedry holohedry;
  Laue laue;
  // Initializers:
  Pointgroup(): number(0), holohedry(Holohedry::_), laue(Laue::_){}
  Pointgroup(const int no): number(no) {setup();}
  Pointgroup(const int no, const std::string& sym, const std::string& sch, const Holohedry& h, const Laue& l):
    number(no), symbol(sym), schoenflies(sch), holohedry(h), laue(l) {}
  std::string to_string(void) const {
    std::string str = "<Pointgroup";
    str += " " + symbol;
    str += " " + schoenflies;
    str += " " + my_to_string(holohedry);
    str += " " + my_to_string(laue);
    return str+">";
  }
  int get_number(void) const {return number;}
  std::string get_symbol(void) const { return symbol; }
  std::string get_schoenflies(void) const { return schoenflies; }
  Holohedry get_holohedry(void) const { return holohedry; }
  Laue get_laue(void) const { return laue; }
  std::string get_holohedry_string(void) const { return my_to_string(holohedry); }
  std::string get_laue_string(void) const { return my_to_string(laue); }
protected:
  void setup(void);
};

Pointgroup ptg_get_transformation_matrix(int *transform_mat, const int *rotations, const int num_rotations);
PointSymmetry ptg_get_pointsymmetry(const int *rotations, const int num_rotations);

int get_pointgroup_rotations_hall_number(int *rotations, const int max_size, const int hall_number, const int is_time_reversal);

int isometry_value(const int *rot);
int rotation_order(const int *rot);
std::array<std::array<int,3>,3> rotation_axis_and_perpendicular_vectors(const int* rot);
std::vector<std::array<int,9>> get_unique_rotations(const std::vector<std::array<int,9>>&, const int);
} // namespace brille
#endif
