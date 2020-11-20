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
#include "bz_mesh.hpp"

#ifndef WRAP_BRILLE_MESH_HPP_
#define WRAP_BRILLE_MESH_HPP_

namespace py = pybind11;
typedef long slong; // ssize_t is only defined for gcc?

template<class T,class R>
void declare_bzmeshq(py::module &m, const std::string &typestr){
  using namespace pybind11::literals;
  using namespace brille;
  using Class = BrillouinZoneMesh3<T,R>;
  std::string pyclass_name = std::string("BZMeshQ")+typestr;
  py::class_<Class> cls(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr());
  // Initializer (BrillouinZone, max-volume, is-volume-rlu)
  cls.def(py::init<BrillouinZone,double,int,int>(), "brillouinzone"_a, "max_size"_a=-1., "num_levels"_a=3, "max_points"_a=-1);
  cls.def_property_readonly("BrillouinZone",[](const Class& cobj){return cobj.get_brillouinzone();});
  cls.def_property_readonly("rlu",[](const Class& cobj){
    return brille::a2py(cobj.get_mesh_hkl());
  });
  cls.def_property_readonly("invA",[](const Class& cobj){
    return brille::a2py(cobj.get_mesh_xyz());
  });
  cls.def_property_readonly("tetrahedra",[](const Class& cobj){
    return brille::a2py(cobj.get_mesh_tetrehedra());
  });
  cls.def("__repr__",&Class::to_string);

  def_grid_fill(cls);
  def_grid_ir_interpolate(cls);
  def_grid_sort(cls);
  def_grid_debye_waller(cls);
}

#endif
