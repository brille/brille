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
#include "_c_to_python.hpp"
#include "polyhedron.hpp"
#include "utilities.hpp"

void wrap_polyhedron(pybind11::module &m){
  using namespace pybind11::literals;
  using namespace brille;
  pybind11::class_<Polyhedron> cls(m, "Polyhedron");

  cls.def(pybind11::init([](py::array_t<double> pyv){
    return Polyhedron(brille::py2a2(pyv));
  }),"vertices"_a);

  cls.def(pybind11::init([](py::array_t<double> pyv, const std::vector<std::vector<int>>& vpf){
    return Polyhedron(brille::py2a2(pyv), vpf);
  }),"vertices"_a, "vertices_per_face"_a);

  cls.def_property_readonly("vertices",[](const Polyhedron& o){return brille::a2py(o.get_vertices());});

  cls.def_property_readonly("points",[](const Polyhedron& o){return brille::a2py(o.get_points());});

  cls.def_property_readonly("normals",[](const Polyhedron& o){return brille::a2py(o.get_normals());});

  cls.def_property_readonly("vertices_per_face",&Polyhedron::get_vertices_per_face);

  cls.def_property_readonly("faces_per_vertex",&Polyhedron::get_faces_per_vertex);

  cls.def_property_readonly("volume",&Polyhedron::get_volume);

  cls.def_property_readonly("mirror",&Polyhedron::mirror);

  cls.def_property_readonly("centre",&Polyhedron::centre);

  cls.def("intersection",&Polyhedron::intersection);

  cls.def("rotate",[](const Polyhedron& o, pybind11::array_t<double> rot){
    pybind11::buffer_info info = rot.request();
    if (info.ndim != 2)
      throw std::runtime_error("Number of dimensions of rotation matrix must be two");
    if (info.shape[0]!=3 || info.shape[1]!=3)
      throw std::runtime_error("Rotation matrix must be 3x3");
    std::array<double, 9> crot;
    double* ptr = (double*) info.ptr;
    auto s0 = info.strides[0]/sizeof(double);
    auto s1 = info.strides[1]/sizeof(double);
    for (size_t i=0; i<3u; ++i) for (size_t j=0; j<3u; ++j)
    crot[i*3u+j] = ptr[i*s0 + j*s1];
    return o.rotate(crot);
  });
}
