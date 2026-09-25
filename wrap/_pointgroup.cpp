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
  pybind11::class_<Pointgroup> cls(m, "Pointgroup",
                                   R"pbdoc(Point group information as in :py:mod:`spglib`

  Wrapped access to the struct originally in
  `pointgroup.h <https://github.com/spglib/spglib/blob/develop/src/pointgroup.h>`_
)pbdoc");
  cls.def(pybind11::init<int>(), pybind11::arg("Pointgroup number"),
          R"pbdoc(Initialize from the serial index in the static `pointgroup_data` array)pbdoc");
  cls.def_property_readonly("number",&Pointgroup::get_number);
  cls.def_property_readonly("symbol",&Pointgroup::get_symbol);
  cls.def_property_readonly("holohedry",&Pointgroup::get_holohedry_string,
                            R"pbdoc(Return a string representation of the Holohedry value

  One of `triclinic`, `monoclinic`, `orthogonal`, `tetragonal`, `trigonal`, `hexagonal`, or `cubic`.
)pbdoc");
  cls.def_property_readonly("laue",&Pointgroup::get_laue_string,
                            R"pbdoc(Return a string representation of the Laue class

  One of `1`, `2m`, `mmm`, `4m`, `4mmm`, `3`, `3m`, `6m`, `6mmm`, `m3`, or `m3m`.
)pbdoc");
  cls.def("__repr__",&Pointgroup::to_string);
}
