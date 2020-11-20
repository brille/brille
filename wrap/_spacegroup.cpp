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
#include "spg_database.hpp"

void wrap_spacegroup(pybind11::module & m){
  using namespace pybind11::literals;
  using namespace brille;
  pybind11::class_<Spacegroup> cls(m,"Spacegroup");

  cls.def(pybind11::init<int>(),"Hall number"_a);

  cls.def_property_readonly("hall_number", &Spacegroup::get_hall_number);

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
