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
#include "_lattice.hpp"

namespace py = pybind11;

void wrap_lattice(py::module &m){
  using namespace pybind11::literals;
  using namespace brille;
  // Declare the interface to the superclass Lattice
  py::class_<Lattice> cls(m,"Lattice");

  cls.def(py::init<double,double,double,double,double,double,int>(),
          "a"_a,"b"_a,"c"_a,"alpha"_a=brille::halfpi,"beta"_a=brille::halfpi,"gamma"_a=brille::halfpi,
          "HallNumber"_a=1);
  //
  cls.def(py::init([](py::array_t<double> lens, py::array_t<double> angs, int hall) {
    py::buffer_info linfo = lens.request(), ainfo = angs.request();
    if ( linfo.ndim!=1 || ainfo.ndim!=1)
      throw std::runtime_error("Number of dimensions must be one");
    if ( linfo.shape[0] < 3 || ainfo.shape[0] < 3 )
      throw std::runtime_error("(At least) three lengths and angles required.");
    return Lattice((double*) linfo.ptr, linfo.strides,
                   (double*) ainfo.ptr, ainfo.strides, hall);
  }),"lengths"_a,"angles"_a,"HallNumber"_a=1);
  //
  cls.def(py::init([](py::array_t<double> vecs, int hall) {
    py::buffer_info info = vecs.request();
    if ( info.ndim!=2 )
      throw std::runtime_error("Number of dimensions must be two");
    if ( info.shape[0] != 3 || info.shape[1] != 3 )
      throw std::runtime_error("Three three-vectors required.");
    return Lattice((double*) info.ptr, info.strides, hall);
  }),"vectors"_a,"HallNumber"_a=1);

  // accessors
  cls.def_property_readonly("a",     &Lattice::get_a);
  cls.def_property_readonly("b",     &Lattice::get_b);
  cls.def_property_readonly("c",     &Lattice::get_c);
  cls.def_property_readonly("alpha", &Lattice::get_alpha);
  cls.def_property_readonly("beta",  &Lattice::get_beta);
  cls.def_property_readonly("gamma", &Lattice::get_gamma);
  cls.def_property_readonly("volume",&Lattice::get_volume);
  cls.def_property_readonly("bravais",&Lattice::get_bravais_type);
  //cls.def_property("spacegroup",&Lattice::get_spacegroup_symmetry,&Lattice::set_spacegroup_symmetry);
  cls.def_property("spacegroup",
    [](const Lattice& l){
    return l.get_spacegroup_symmetry();
    },
    [](Lattice& l, const Symmetry& s){
    return l.set_spacegroup_symmetry(s);
  });
  //cls.def_property_readonly("pointgroup",&Lattice::get_pointgroup_symmetry);
  cls.def_property_readonly("pointgroup",
    [](const Lattice& l){
    return l.get_pointgroup_symmetry();
    });
  // The next line would require that the Basis type is exposed to Python
  //cls.def_property_readonly("basis",&Lattice::get_basis);

  cls.def("get_covariant_metric_tensor",[](Lattice &l){
    auto result = py::array_t<double, py::array::c_style >({3,3});
    py::buffer_info bi = result.request();
    double *cmt = (double *) bi.ptr;
    l.get_covariant_metric_tensor( cmt );
    return result;
  });

  cls.def("get_contravariant_metric_tensor",[](Lattice &l){
    auto result = py::array_t<double, py::array::c_style >({3,3});
    py::buffer_info bi = result.request();
    l.get_contravariant_metric_tensor((double *) bi.ptr );
    return result;
  });

  cls.def_property_readonly("star",[](Lattice &){throw std::runtime_error("Bare Lattices do not have a reciprocal!");});
  cls.def("issame",&Lattice::issame);
  cls.def("__eq__",&Lattice::issame);
  cls.def("isapprox",&Lattice::isapprox);
  cls.def("isstar",[](Lattice &, Lattice ){throw std::runtime_error("Bare Lattices do not have a reciprocal!");});
  cls.def("__repr__",&Lattice::string_repr);

  // Inform pybind11 that the specializations Direct and Reciprocal exist
  py::class_<Direct,Lattice> direct(m,"Direct");
  py::class_<Reciprocal,Lattice> reciprocal(m,"Reciprocal");
  // and define their common extensions
  declare_lattice_methods<Direct>(direct,"Å");
  declare_lattice_methods<Reciprocal>(reciprocal,"Å⁻¹");

  reciprocal.def_property_readonly("B",[](Reciprocal &r){
    auto result = py::array_t<double, py::array::c_style>({3,3});
    py::buffer_info bi = result.request();
    r.get_B_matrix((double *)bi.ptr);
    return result;
  });
}
