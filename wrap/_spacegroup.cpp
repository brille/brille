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
#include <stdexcept>
#include <string>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "spg_database.hpp"

void wrap_spacegroup(pybind11::module & m){
  using namespace pybind11::literals;
  using namespace brille;
  pybind11::class_<Spacegroup> cls(m,"Spacegroup",
                                   R"pbdoc(The space group information as in :py:mod:`spglib:

  Equivalent to the struct `SpacegroupType` from
  `spg_database.h <https://github.com/spglib/spglib/blob/develop/src/spg_database.h>`_
)pbdoc");

  cls.def(pybind11::init([](const std::string& symbol, const std::string& choice){
    const auto n = string_to_hall_number(symbol, choice);
    if (n <= 0 || n >= 531)
      throw std::invalid_argument("'" + symbol + "'" + (choice.empty() ? "" : " with choice '" + choice + "'")
        + " is not a space group's Hall symbol, or its Hermann-Mauguin symbol or International Tables name"
          " (with a valid setting choice)");
    return Spacegroup(n);
  }), "symbol"_a, "choice"_a="", R"pbdoc(
    The space group setting named by a Hall symbol, or by a Hermann-Mauguin symbol or
    International Tables name with an optional setting choice, as
    :py:func:`brille.Lattice` accepts them.
  )pbdoc");

  cls.def_static("all", [](){
    std::vector<Spacegroup> out;
    for (int n = 1; n < 531; ++n) out.emplace_back(n);
    return out;
  }, R"pbdoc(
    Every space group setting brille knows, in the order of its table.
  )pbdoc");

  cls.def_property_readonly("international_table_number", &Spacegroup::get_international_table_number);

  cls.def_property_readonly("pointgroup_number", &Spacegroup::get_pointgroup_number);

  cls.def_property_readonly("schoenflies_symbol", &Spacegroup::get_schoenflies_symbol);

  cls.def_property_readonly("hall_symbol", &Spacegroup::get_hall_symbol);

  cls.def_property_readonly("international_table_symbol", &Spacegroup::get_international_table_symbol);

  cls.def_property_readonly("international_table_full", &Spacegroup::get_international_table_full);

  cls.def_property_readonly("international_table_short", &Spacegroup::get_international_table_short);

  cls.def_property_readonly("choice", &Spacegroup::get_choice);

  cls.def("__repr__",&Spacegroup::string_repr);
}
