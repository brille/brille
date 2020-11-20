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
#include "pointgroup.hpp"

void wrap_pointgroup(pybind11::module & m){
  using namespace brille;
  pybind11::class_<Pointgroup> cls(m, "Pointgroup");
  cls.def(pybind11::init<int>(), pybind11::arg("Pointgroup number"));
  cls.def_property_readonly("number",&Pointgroup::get_number);
  cls.def_property_readonly("symbol",&Pointgroup::get_symbol);
  cls.def_property_readonly("holohedry",&Pointgroup::get_holohedry_string);
  cls.def_property_readonly("laue",&Pointgroup::get_laue_string);
  cls.def("__repr__",&Pointgroup::to_string);
}
