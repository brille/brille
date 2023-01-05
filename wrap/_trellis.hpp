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
#include "bz_trellis.hpp"
#include "approx_config.hpp"

#ifndef WRAP_BRILLE_TRELLIS_HPP_
#define WRAP_BRILLE_TRELLIS_HPP_

namespace py = pybind11;

template<class T,class R,class S>
void declare_bztrellisq(py::module &m, const std::string &typestr){
  using namespace pybind11::literals;
  using namespace brille;
  using Class = BrillouinZoneTrellis3<T,R,S>;
  std::string pyclass_name = std::string("BZTrellisQ")+typestr;
  py::class_<Class> cls(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr());
  // Initializer (BrillouinZone, maximum node volume fraction, always_triangulate)
  cls.def(py::init<BrillouinZone,double,bool>(), "brillouin_zone"_a, "node_volume_fraction"_a=0.1, "always_triangulate"_a=false);
  cls.def(py::init<BrillouinZone,double,bool,approx_float::Config>(), "brillouin_zone"_a, "node_volume_fraction"_a, "always_triangulate"_a, "approx_config"_a);

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

  cls.def("node_at", [](const Class& cobj, const std::array<ind_t, 3>& sub){
    return cobj.subscripted_node_poly(sub);
  }, "subscript"_a);

  cls.def("node_containing", [](const Class & c, py::array_t<double> pyX){
    using namespace brille;
    using namespace brille::lattice;
    brille::Array2<double> bX = brille::py2a2(pyX);
    auto qv = LQVec<double>(c.get_brillouinzone().get_lattice(), bX);
    return c.point_in_node_poly(qv);
  }, "Q"_a);

  cls.def("node_at_type", [](const Class & c, const std::array<ind_t, 3>& s){
    return c.node_at_type(s);
  }, "subscript"_a);

  cls.def("node_containing_type", [](const Class & c, py::array_t<double> pyX){
    using namespace brille;
    using namespace brille::lattice;
    brille::Array2<double> bX = brille::py2a2(pyX);
    auto qv = LQVec<double>(c.get_brillouinzone().get_lattice(), bX);
    return c.point_in_node_type(qv);
  }, "Q"_a);

  cls.def("all_node_types", [](const Class & c){
    return c.all_node_types();
  });

#ifdef USE_HIGHFIVE
  def_grid_hdf_interface(cls, pyclass_name);
#endif
}

#endif
