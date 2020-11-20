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
#include "_c_to_python.hpp"
#include "primitive.hpp"

void wrap_primitivetransform(pybind11::module & m){
  using namespace pybind11::literals;
  using namespace brille;
  pybind11::class_<PrimitiveTransform> cls(m,"PrimitiveTransform");

  cls.def(pybind11::init<int>(),"Hall number"_a);

  cls.def_property_readonly("P",[](const PrimitiveTransform &p){
    std::vector<ssize_t> sz={3,3};
    return sa2np(sz,p.get_P());
  });

  cls.def_property_readonly("invP",[](const PrimitiveTransform &p){
    std::vector<ssize_t> sz={3,3};
    return sa2np(sz,p.get_invP());
  });

  cls.def_property_readonly("Pt",[](const PrimitiveTransform &p){
    std::vector<ssize_t> sz={3,3};
    return sa2np(sz,p.get_Pt());
  });

  cls.def_property_readonly("invPt",[](const PrimitiveTransform &p){
    std::vector<ssize_t> sz={3,3};
    return sa2np(sz,p.get_invPt());
  });

  cls.def_property_readonly("does_anything",&PrimitiveTransform::does_anything);

  cls.def_property_readonly("is_primitive",&PrimitiveTransform::does_nothing);

  cls.def("__repr__",&PrimitiveTransform::string_repr);
}
