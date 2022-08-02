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
#include "_array.hpp"
#include "_common_grid.hpp"
#include "nest.hpp"
#include "bz_nest.hpp"

#ifndef WRAP_BRILLE_NEST_HPP_
#define WRAP_BRILLE_NEST_HPP_

namespace py = pybind11;

template<class T, class R, class S>
void declare_bznestq(py::module &m, const std::string &typestr){
  using namespace pybind11::literals;
  using namespace brille;
  using Class = BrillouinZoneNest3<T,R,S>;
  std::string pyclass_name = std::string("BZNestQ")+typestr;
  py::class_<Class> cls(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr());
  // Initializer (BrillouinZone, maximum node volume fraction)
  cls.def(py::init<BrillouinZone,double,ind_t>(), "brillouin_zone"_a, "max_volume"_a, "max_branchings"_a=5);
  cls.def(py::init<BrillouinZone,ind_t,ind_t>(), "brillouin_zone"_a, "number_density"_a, "max_branchings"_a=5);

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
//  def_grid_debye_waller(cls);

#ifdef USE_HIGHFIVE
  def_grid_hdf_interface(cls, pyclass_name);
#endif
}

#endif
