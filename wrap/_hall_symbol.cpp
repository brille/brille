/* Copyright 2019 Greg Tucker
//
// This file is part of brille.
//
// brille is free software: you can redistribute it and/or modify it under the
// terms of the GNU Affero General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// brille is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with brille. If not, see <https://www.gnu.org/licenses/>.            */
#include <pybind11/pybind11.h>
#include "hall_symbol.hpp"

void wrap_hallsymbol(pybind11::module &m){
  using namespace pybind11::literals;
  pybind11::class_<HallSymbol> cls(m, "HallSymbol");

  cls.def(pybind11::init([](std::string symbol){
    HallSymbol out;
    out.from_ascii(symbol);
    return out;
  }),"Hall_symbol"_a);

  cls.def("__repr__",&HallSymbol::to_ascii);

  cls.def_property_readonly("generators",&HallSymbol::get_generators);
}
