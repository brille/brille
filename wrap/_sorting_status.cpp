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
#include "sorting_status.hpp"

void wrap_sortingstatus(pybind11::module &m){
  using namespace pybind11::literals;
  using namespace brille;
  pybind11::class_<SortingStatus> cls(m, "SortingStatus", R"pbdoc(
    An object representing the status of a single object under sorting.

    Internally may be represented as a single unsigned integer where two or
    more bits are reserved for various flags, or as a set of boolean values
    and an integer.

  )pbdoc");

  cls.def(pybind11::init<bool, bool, unsigned int>(),"sorted"_a,"locked"_a,"visits"_a);
  cls.def_property_readonly("sorted",(bool (SortingStatus::*)(void) const) &SortingStatus::sorted,"Return the sorted flag");
  cls.def_property_readonly("locked",(bool (SortingStatus::*)(void) const) &SortingStatus::locked,"Return the locked flag");
  cls.def_property_readonly("visits",(SortingStatus::status_t (SortingStatus::*)(void) const) &SortingStatus::visits,"Return the visit count");
  cls.def("__repr__",&SortingStatus::to_string);
}
