/* This file is part of brille.

Copyright © 2020 Greg Tucker <greg.tucker@stfc.ac.uk>

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
#include "bravais.hpp"

namespace py = pybind11;

void wrap_bravais(py::module &m){
  using namespace brille;
  py::enum_<Bravais>(m, "Bravais")
    .value("invalid", Bravais::_)
    .value("P", Bravais::P)
    .value("A", Bravais::A)
    .value("B", Bravais::B)
    .value("C", Bravais::C)
    .value("I", Bravais::I)
    .value("F", Bravais::F)
    .value("R", Bravais::R);
}
