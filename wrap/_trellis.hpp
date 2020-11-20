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

#include <iostream>
#include <pybind11/pybind11.h>
#include "_common_grid.hpp"
#include "trellis.hpp"
#include "bz_trellis.hpp"

#ifndef WRAP_BRILLE_TRELLIS_HPP_
#define WRAP_BRILLE_TRELLIS_HPP_

namespace py = pybind11;

template<class T,class R>
void declare_bztrellisq(py::module &m, const std::string &typestr){
  using namespace pybind11::literals;
  using namespace brille;
  using Class = BrillouinZoneTrellis3<T,R>;
  std::string pyclass_name = std::string("BZTrellisQ")+typestr;
  py::class_<Class> cls(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr());
  // Initializer (BrillouinZone, maximum node volume fraction, always_triangulate)
  cls.def(py::init<BrillouinZone,double,bool>(), "brillouinzone"_a, "node_volume_fraction"_a=0.1, "always_triangulate"_a=false);

  cls.def_property_readonly("BrillouinZone",[](const Class& cobj){return cobj.get_brillouinzone();});

  cls.def_property_readonly("invA",[](const Class& cobj){
    return brille::a2py(cobj.get_xyz());
  });
  cls.def_property_readonly("inner_invA",[](const Class& cobj){
    return brille::a2py(cobj.get_inner_xyz());
  });
  cls.def_property_readonly("outer_invA",[](const Class& cobj){
    return brille::a2py(cobj.get_outer_xyz());
  });

  cls.def_property_readonly("rlu",[](const Class& cobj){
    return brille::a2py(cobj.get_hkl());
  });
  cls.def_property_readonly("inner_rlu",[](const Class& cobj){
    return brille::a2py(cobj.get_inner_hkl());
  });
  cls.def_property_readonly("outer_rlu",[](const Class& cobj){
    return brille::a2py(cobj.get_outer_hkl());
  });

  cls.def_property_readonly("tetrahedra",[](const Class& cobj){
    return cobj.get_vertices_per_tetrahedron();
  });

  // cls.def("__repr__",&Class::to_string);
  def_grid_fill(cls);
  def_grid_interpolate(cls);
  def_grid_ir_interpolate(cls);
  def_grid_sort(cls);
  def_grid_debye_waller(cls);
}

#endif
