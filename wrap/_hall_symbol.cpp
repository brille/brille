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
#include <pybind11/pybind11.h>
#include "hall_symbol.hpp"

void wrap_hallsymbol(pybind11::module &m){
  using namespace pybind11::literals;
  using namespace brille;
  pybind11::class_<HallSymbol> cls(m, "HallSymbol", R"pbdoc(
    A crystallographic spacegroup's symmetries encoded in Hall's notation

    Hall proposed a compact unambiguous notation for the representation of the
    generators of a spacegroup. Within his notation each motion is comprised of
    a character with one or more subscripts and superscripts which describe its
    order, unique axis, and translation. The notation specifies that, depending
    on the position of a motion and details of any preceeding motion, some or
    all of the sub- and superscripts can be omitted. The :class:`HallSymbol` has
    been written to handle the logic necessary to decode a Hall symbol into its
    equivalent motions.
    An added compilcation arises when the Hall symbol is encoded as an ASCII
    string. Namely, there are no sub- or superscript glyphs and some scheme must
    be enacted to represent them.

  )pbdoc");

  cls.def(pybind11::init([](std::string symbol){
    HallSymbol out;
    out.from_ascii(symbol);
    return out;
  }),"Hall_symbol"_a);

  cls.def("__repr__",&HallSymbol::to_ascii);

  cls.def_property_readonly("generators",&HallSymbol::get_generators);
}
