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
#include <pybind11/operators.h>

#include "_c_to_python.hpp"
#include "symmetry.hpp"
#include "spg_database.hpp"

/*!
The symmetry operations of a spacegroup
*/
void wrap_symmetry(pybind11::module & m){
  using namespace pybind11::literals;
  pybind11::class_<Symmetry> cls(m, "Symmetry");

  cls.def(pybind11::init([](int hall){return Spacegroup(hall).get_spacegroup_symmetry();}),"Hall number"_a);

  cls.def_property_readonly("size",&Symmetry::size);

  cls.def_property_readonly("W",[](Symmetry& ps){
    std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size()), 3, 3};
    return sva2np(sz, ps.getallr());
  });

  cls.def_property_readonly("w",[](Symmetry& ps){
    std::vector<ssize_t> sz={static_cast<ssize_t>(ps.size()), 3};
    return sva2np(sz, ps.getallt());
  });

  cls.def("generate",&Symmetry::generate);

  cls.def("generators",&Symmetry::generators);

  cls.def(pybind11::self == pybind11::self);
}
