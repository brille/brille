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

/* This file has evolved from spg_database.h distributed as part of spglib.
   Changes have been made to introduce C++ style classes as well as other
   modifications to suit the needs of brille.
   spglib was licensed under the following BSD 3-clause license:

 Copyright (C) 2010 Atsushi Togo
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

#ifndef BRILLE_SPACEGROUP_DATABASE_HPP_
#define BRILLE_SPACEGROUP_DATABASE_HPP_
#include <cstring>
#include "bravais.hpp"
#include "pointgroup.hpp"
namespace brille {

class Spacegroup{
public:
  int number;
  std::string schoenflies;
  std::string hall_symbol;
  std::string international;
  std::string international_full;
  std::string international_short;
  std::string choice;
  Bravais bravais;
  int pointgroup_number;
  int hall_number;
  // Initializers
  Spacegroup(): number(0), bravais(Bravais::_), pointgroup_number(0), hall_number(0) {}
  Spacegroup(const int no, const std::string& sf, const std::string& hs,
             const std::string& its, const std::string& itf,
             const std::string& ith, const std::string& ch,
             const Bravais br, const int pno, const int hno):
    number(no), schoenflies(sf), hall_symbol(hs), international(its),
    international_full(itf), international_short(ith), choice(ch), bravais(br),
    pointgroup_number(pno), hall_number(hno) {};
  Spacegroup(int _hall_number) { set_from_hall_number(_hall_number); }
  //
  int get_hall_number(void) const { return this->hall_number; }
  int get_international_table_number(void) const { return this->number; }
  int get_pointgroup_number(void) const { return this->pointgroup_number; }
  std::string get_schoenflies_symbol(void) const {return this->schoenflies; }
  std::string get_hall_symbol(void) const {return this->hall_symbol; }
  std::string get_international_table_symbol(void) const {return this->international; }
  std::string get_international_table_full(void) const {return this->international_full; }
  std::string get_international_table_short(void) const {return this->international_short; }
  std::string get_choice(void) const {return this->choice; }
  Bravais get_bravais_type() const { return this->bravais; }
  int set_hall_number(const int nh) { this->hall_number = nh; return this->hall_number; }
  int set_international_table_number(const int itn) { this->number=itn; return this->number; }
  int set_pointgroup_number(const int pgn) { this->pointgroup_number=pgn; return this->pointgroup_number; }
  std::string set_schoenflies_symbol(const std::string& ns) { this->schoenflies=ns; return this->schoenflies; }
  std::string set_hall_symbol(const std::string& ns) { this->hall_symbol=ns; return this->hall_symbol; }
  std::string set_international_table_symbol(const std::string& ns) {this->international=ns; return this->international; }
  std::string set_international_table_full(const std::string& ns) {this->international_full=ns; return this->international_full; }
  std::string set_international_table_short(const std::string& ns) {this->international_short=ns; return this->international_short; }
  std::string set_choice(const std::string& ns) {this->choice=ns; return this->choice; }
  Bravais set_bravais_type(const Bravais b) {this->bravais = b; return this->bravais; }
  std::string string_repr(void) const {
    std::string repr;
    repr += " IT(" + std::to_string(this->number) + "): "
          + this->international;
    if (this->choice.size()>0)
      repr += " [" + this->choice + "]";
    repr += ";";
    repr += " Hall(" + std::to_string(this->hall_number) + "): "
          + this->hall_symbol;
    std::string left("< Spacegroup "), right(" >");
    return left + repr + right;
  }
  Pointgroup get_pointgroup(void) const {return Pointgroup(this->pointgroup_number);}
  Symmetry get_spacegroup_symmetry() const;
  PointSymmetry get_pointgroup_symmetry(const int time_reversal=0) const;
private:
  // void set_hall_number(void);
  void set_from_hall_number(const int);
};

bool hall_number_ok(const int hall_number);
// PointSymmetry make_pointgroup_symmetry_object(const int hall_number, const int time_reversal=0);

int international_number_to_hall_number(const int number, const std::string& choice="");
int international_string_to_hall_number(const std::string& itname, const std::string& choice="");
int hall_symbol_to_hall_number(const std::string& hsymbol);
int string_to_hall_number(const std::string&, const std::string& choice="");

} // end namespace brille
#endif
