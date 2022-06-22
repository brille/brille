/* This file is part of brille.

Copyright © 2022 Greg Tucker <gregory.tucker@ess.eu>

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
#include "enums.hpp"

namespace py = pybind11;

void wrap_enums(py::module &m){
  using namespace brille;
  py::enum_<AngleUnit>(m, "AngleUnit", R"pbdoc(
  The units of a number representing an angle
  )pbdoc")
    .value("not_provided", AngleUnit::not_provided, R"pbdoc(unknown, will be inferred)pbdoc")
    .value("radian", AngleUnit::radian, R"pbdoc(radian)pbdoc")
    .value("degree", AngleUnit::degree, R"pbdoc(degree)pbdoc")
    .value("pi", AngleUnit::pi, R"pbdoc(radian divided by pi)pbdoc");

  py::enum_<LengthUnit>(m, "LengthUnit", R"pbdoc(
  The units of a number representing a length
  )pbdoc")
    .value("none", LengthUnit::none)
    .value("angstrom", LengthUnit::angstrom, R"pbdoc(10⁻¹⁰ meter)pbdoc")
    .value("inverse_angstrom", LengthUnit::inverse_angstrom, R"pbdoc(1 / angstrom)pbdoc");
}
