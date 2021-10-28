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
/*! \file
    \author Atsushi Togo
    \author Greg Tucker
    \brief Lattice spacegroup information class and utility functions

    The Spacegroup class is a C++ rewrite of the C struct `SpacegroupType` from
    [spg_database.h](https://github.com/spglib/spglib/blob/develop/src/spg_database.h)
    which is part of [spglib](https://github.com/atztogo/spglib).

*/
#include <cstring>
// #include "bravais.hpp"
#include "pointgroup.hpp"
namespace brille {

/*! \brief Information about a lattice spacegroup

There are many possible descriptions of a lattice; this class holds a subset
of this information for a given lattice. The information is sourced from the
`ALL_SPACEGROUPS` constant table.
*/
class Spacegroup{
public:
  int number;                        //!< International Tables spacegroup number
  std::string schoenflies;           //!< Schoenflies symbol
  std::string hall_symbol;           //!< Hall symbol
  std::string international;         //!< International Table spacegroup name
  std::string international_full;    //!< Hermann-Mauguin spacegroup name
  std::string international_short;   //!< Compact International Table spacegroup name
  std::string choice;                //!< Unique axis, setting, or cell choice
  Bravais bravais;                   //!< Lattice centring type
  int pointgroup_number;             //!< The pointgroup number
  int hall_number;                   //!< The 'serial no.' index from [Seto's Home Page](http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en)
  // Initializers
  //! Empty initializer
  Spacegroup(): number(0), bravais(Bravais::_), pointgroup_number(0), hall_number(0) {}
  //! Construct with all required information
  Spacegroup(const int no, const std::string& sf, const std::string& hs,
             const std::string& its, const std::string& itf,
             const std::string& ith, const std::string& ch,
             const Bravais br, const int pno, const int hno):
    number(no), schoenflies(sf), hall_symbol(hs), international(its),
    international_full(itf), international_short(ith), choice(ch), bravais(br),
    pointgroup_number(pno), hall_number(hno) {};
  //! Construct from a serial 'Hall number'
  Spacegroup(int _hall_number) { set_from_hall_number(_hall_number); }
  //! Return the serial 'Hall number' for this Spacegroup
  int get_hall_number(void) const { return this->hall_number; }
  //! Return the spacegroup number from the International Tables of Crystallography
  int get_international_table_number(void) const { return this->number; }
  //! Return the pointgroup serial index
  int get_pointgroup_number(void) const { return this->pointgroup_number; }
  //! Return the Schoenflies notation symbol
  std::string get_schoenflies_symbol(void) const {return this->schoenflies; }
  //! Return the generators of the Spacegroup encoded in a Hall symbol
  std::string get_hall_symbol(void) const {return this->hall_symbol; }
  //! Return the International Tables of Crystallography spacegroup name
  std::string get_international_table_symbol(void) const {return this->international; }
  //! Return the Herman-Mauguin spacegroup name
  std::string get_international_table_full(void) const {return this->international_full; }
  //! Return a compact version of th International Tables of Crystallography spacegroup name
  std::string get_international_table_short(void) const {return this->international_short; }
  //! Return the unique axis, setting, or cell choice
  std::string get_choice(void) const {return this->choice; }
  //! Return the lattice centring type Bravais enumerated value
  Bravais get_bravais_type() const { return this->bravais; }
  //! Set the serial 'Hall number'
  int set_hall_number(const int nh) { this->hall_number = nh; return this->hall_number; }
  //! Set the International Tables of Crystallography spacegroup number
  int set_international_table_number(const int itn) { this->number=itn; return this->number; }
  //! Set the pointgroup serial index
  int set_pointgroup_number(const int pgn) { this->pointgroup_number=pgn; return this->pointgroup_number; }
  //! Set the Schoenflies notation symbol
  std::string set_schoenflies_symbol(const std::string& ns) { this->schoenflies=ns; return this->schoenflies; }
  //! Set the generators of the spacegroup in Hall notation
  std::string set_hall_symbol(const std::string& ns) { this->hall_symbol=ns; return this->hall_symbol; }
  //! Set the International Tables of Crystallography spacegroup name
  std::string set_international_table_symbol(const std::string& ns) {this->international=ns; return this->international; }
  //! Set the Hermann-Mauguin spacegroup name
  std::string set_international_table_full(const std::string& ns) {this->international_full=ns; return this->international_full; }
  //! Set the compact version of the Internation Tables of Crystallography spacegroup name
  std::string set_international_table_short(const std::string& ns) {this->international_short=ns; return this->international_short; }
  //! Set the unique axis, setting, or cell choice
  std::string set_choice(const std::string& ns) {this->choice=ns; return this->choice; }
  //! Set the lattice centring type
  Bravais set_bravais_type(const Bravais b) {this->bravais = b; return this->bravais; }
  //! Return a string representation of the information contained in this object
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
  //! Return the Pointgroup associated with this object's pointgroup serial index
  Pointgroup get_pointgroup(void) const {return Pointgroup(this->pointgroup_number);}
  //! Return the Spacegroup Symmetry object associated with this object
  Symmetry get_spacegroup_symmetry() const;
  //! Return the Pointgroup PointSymmetry obect associated with this object
  PointSymmetry get_pointgroup_symmetry(const int time_reversal=0) const;
private:
  // void set_hall_number(void);
  void set_from_hall_number(const int);
};
//! Verify that the provided serial 'Hall number' is valid
bool hall_number_ok(const int hall_number);
// PointSymmetry make_pointgroup_symmetry_object(const int hall_number, const int time_reversal=0);
//! Find the serial 'Hall number' for a specified International Tables of Crystallography spacegroup number
int international_number_to_hall_number(const int number, const std::string& choice="");
//! Find the serial 'Hall number' for a specified International Tables of Crystallography spacegroup name
int international_string_to_hall_number(const std::string& itname, const std::string& choice="");
//! Find the serial 'Hall number' for one of Hall's proposed 530 spacegroup generator symbols
int hall_symbol_to_hall_number(const std::string& hsymbol);
/*! \brief Locate a serial 'Hall number' from a string

Uses both `international_string_to_hall_number` and `hall_symbol_to_hall_number`
to find the serial 'Hall number' for a provided string.

\param symbol The International Tables of Crystallography spacegroup name
              (or its shortened variant), Hermann-Mauguin symbol, or one of
              Hall's 530 proposed spacegroup symbols.
\param choice The unique axis, setting, or cell choice; if required for the
              differentiation between otherwise identical symbols.
\returns A valid serial 'Hall number' or zero
*/
int string_to_hall_number(const std::string& symbol, const std::string& choice="");

} // end namespace brille
#endif
