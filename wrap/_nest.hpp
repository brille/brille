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

#include "_array.hpp"
#include "_common_grid.hpp"

#include "nest.hpp"
#include "bz_nest.hpp"

#ifndef __NEST_H
#define __NEST_H

namespace py = pybind11;
typedef long slong; // ssize_t is only defined for gcc?
typedef unsigned long element_t;

template<class T, class R>
void declare_bznestq(py::module &m, const std::string &typestr){
  using namespace pybind11::literals;
  using Class = BrillouinZoneNest3<T,R>;
  std::string pyclass_name = std::string("BZNestQ")+typestr;
  py::class_<Class> cls(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr());
  // Initializer (BrillouinZone, maximum node volume fraction)
  cls.def(py::init<BrillouinZone,double,size_t>(), "brillouinzone"_a, "max_volume"_a, "max_branchings"_a=5);
  cls.def(py::init<BrillouinZone,size_t,size_t>(), "brillouinzone"_a, "number_density"_a, "max_branchings"_a=5);

  cls.def_property_readonly("BrillouinZone",[](const Class& cobj){return cobj.get_brillouinzone();});

  cls.def_property_readonly("invA",[](const Class& cobj){
    return brille::a2py(cobj.get_xyz());
  });
  cls.def_property_readonly("all_invA",[](const Class& cobj){
    return brille::a2py(cobj.get_all_xyz());
  });

  cls.def_property_readonly("rlu",[](const Class& cobj){
    return brille::a2py(cobj.get_hkl());
  });
  cls.def_property_readonly("all_rlu",[](const Class& cobj){
    return brille::a2py(cobj.get_all_hkl());
  });

  cls.def_property_readonly("tetrahedra",[](const Class& cobj){
    return cobj.tetrahedra();
  });

  cls.def_property_readonly("bytes_per_point", &Class::bytes_per_point);

  def_grid_fill(cls);
  def_grid_ir_interpolate(cls);
  def_grid_sort(cls);
  def_grid_debye_waller(cls);
}

#endif
