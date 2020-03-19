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
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <thread>

#include "_c_to_python.hpp"
#include "lattice.hpp"
#include "utilities.hpp"

#ifndef __BINDING_H
#define __BINDING_H

namespace py = pybind11;
typedef long slong; // ssize_t is only defined for gcc?

template<class R, class T>
void declare_lattice_scl_init(py::class_<T,Lattice> &pclass, const std::string &lenunit, const std::string &argname, const R defarg){
  pclass.def(py::init<double,double,double, double,double,double, R>(),
    py::arg( ("a/"+lenunit).c_str() ),
    py::arg( ("b/"+lenunit).c_str() ),
    py::arg( ("c/"+lenunit).c_str() ),
    py::arg("alpha/radian")=PI/2,
    py::arg("beta/radian")=PI/2,
    py::arg("gamma/radian")=PI/2,
    py::arg(argname.c_str())=defarg );
}
template<class R, class T>
void declare_lattice_vec_init(py::class_<T,Lattice> &pclass, const std::string &lenunit, const std::string &argname, const R defarg){
  pclass.def(py::init( [](py::array_t<double> lens, py::array_t<double> angs, const R groupid){
    py::buffer_info linfo = lens.request(), ainfo = angs.request();
    if ( linfo.ndim!=1 || ainfo.ndim!=1)
      throw std::runtime_error("Number of dimensions must be one");
    if ( linfo.shape[0] < 3 || ainfo.shape[0] < 3 )
      throw std::runtime_error("(At least) three lengths and angles required.");
    return T((double*) linfo.ptr, linfo.strides,
             (double*) ainfo.ptr, ainfo.strides, groupid);
  }), py::arg( ("lengths/"+lenunit).c_str() ),
      py::arg("angles/radian"), py::arg(argname.c_str())=defarg);
}
template<class R, class T>
void declare_lattice_mat_init(py::class_<T,Lattice> &pclass, const std::string &argname, const R defarg){
  pclass.def(py::init( [](py::array_t<double> vecs, const R groupid){
    py::buffer_info info = vecs.request();
    if ( info.ndim!=2 )
      throw std::runtime_error("Number of dimensions must be two");
    if ( info.shape[0] != 3 || info.shape[0] != 3 )
      throw std::runtime_error("Three three-vectors required.");
    return T((double *) info.ptr, info.strides, groupid);
  }), py::arg( "lattice vectors" ), py::arg(argname.c_str())=defarg);
}

template<class T>
void declare_lattice_methods(py::class_<T,Lattice> &pclass, const std::string &lenunit) {
    declare_lattice_scl_init(pclass,lenunit,"hall",1);
    declare_lattice_scl_init(pclass,lenunit,"ITname","P_1");
    declare_lattice_vec_init(pclass,lenunit,"hall",1);
    declare_lattice_vec_init(pclass,lenunit,"ITname","P_1");
    declare_lattice_mat_init(pclass,"hall",1);
    declare_lattice_mat_init(pclass,"ITname","P_1");
    pclass.def_property_readonly("star",&T::star,"Return the dual lattice");
    pclass.def_property_readonly("xyz_transform",[](T &d){
      auto result = py::array_t<double, py::array::c_style>({3,3});
      py::buffer_info bi = result.request();
      size_t c = static_cast<size_t>(bi.strides[0]/sizeof(double));
      size_t r = static_cast<size_t>(bi.strides[1]/sizeof(double));
      d.get_xyz_transform( (double *) bi.ptr, c, r);
      return result;
    });
    pclass.def_property_readonly("lattice_matrix",[](T &d){
      auto result = py::array_t<double, py::array::c_style>({3,3});
      py::buffer_info bi = result.request();
      size_t c = static_cast<size_t>(bi.strides[0]/sizeof(double));
      size_t r = static_cast<size_t>(bi.strides[1]/sizeof(double));
      d.get_lattice_matrix( (double *) bi.ptr, c, r);
      return result;
    });
    pclass.def("isstar",(bool (T::*)(const Direct&    ) const) &T::isstar);
    pclass.def("isstar",(bool (T::*)(const Reciprocal&) const) &T::isstar);
    pclass.def_property_readonly("primitive",&T::primitive,"Return the primitive lattice");
}

#endif
