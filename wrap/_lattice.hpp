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
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <thread>

#include "_array.hpp"
#include "_c_to_python.hpp"
#include "lattice.hpp"
#include "array.hpp"
#include "utilities.hpp"


#ifndef WRAP_BRILLE_LATTICE_HPP_
#define WRAP_BRILLE_LATTICE_HPP_
namespace py = pybind11;
namespace br = brille;

template<class R, class T>
void declare_lattice_scl_init(py::class_<T,br::Lattice> &pclass, const std::string &lenunit, const std::string &argname, const R defarg){
  pclass.def(py::init<double,double,double, double,double,double, R>(),
    py::arg( ("a/"+lenunit).c_str() ),
    py::arg( ("b/"+lenunit).c_str() ),
    py::arg( ("c/"+lenunit).c_str() ),
    py::arg("alpha")=90.0,
    py::arg("beta")=90.0,
    py::arg("gamma")=90.0,
    py::arg(argname.c_str())=defarg );
}
template<class R, class T>
void declare_lattice_vec_init(py::class_<T,br::Lattice> &pclass, const std::string &lenunit, const std::string &argname, const R defarg){
  pclass.def(py::init( [](py::array_t<double> lens, py::array_t<double> angs, const R groupid){
    py::buffer_info linfo = lens.request(), ainfo = angs.request();
    if ( linfo.ndim!=1 || ainfo.ndim!=1)
      throw std::runtime_error("Number of dimensions must be one");
    if ( linfo.shape[0] < 3 || ainfo.shape[0] < 3 )
      throw std::runtime_error("(At least) three lengths and angles required.");
    return T((double*) linfo.ptr, linfo.strides,
             (double*) ainfo.ptr, ainfo.strides, groupid);
  }), py::arg( ("lengths/"+lenunit).c_str() ),
      py::arg("angles"), py::arg(argname.c_str())=defarg);
}
template<class R, class T>
void declare_lattice_mat_init(py::class_<T,br::Lattice> &pclass, const std::string &argname, const R defarg){
  pclass.def(py::init( [](py::array_t<double> vecs, const R groupid){
    py::buffer_info info = vecs.request();
    if ( info.ndim!=2 )
      throw std::runtime_error("Number of dimensions must be two");
    if ( info.shape[0] != 3 || info.shape[1] != 3 )
      throw std::runtime_error("Three three-vectors required.");
    return T((double *) info.ptr, info.strides, groupid);
  }), py::arg( "lattice vectors" ), py::arg(argname.c_str())=defarg);
}

template<class R, class T>
void declare_lattice_mat_basis_init(py::class_<T,br::Lattice> &pclass, const std::string &argname, const R defarg){
  pclass.def(py::init( [](py::array_t<double> vecs, py::array_t<double> pypos, std::vector<unsigned long> typ, const R groupid){
    py::buffer_info info = vecs.request();
    if ( info.ndim!=2 )
      throw std::runtime_error("Number of dimensions must be two");
    if ( info.shape[0] != 3 || info.shape[1] != 3 )
      throw std::runtime_error("Three three-vectors required.");
    T lattice((double *) info.ptr, info.strides, groupid);
    // pull in the basis information
    // br::Array<double> avpos = br::py2a(pypos);
    br::Array2<double> avpos = br::py2a2(pypos);
    std::vector<std::array<double,3>> pos;
    for (br::ind_t i=0; i<avpos.size(0); ++i){
      std::array<double,3> one;
      for (br::ind_t j=0; j<3u; ++j){
        one[j] = avpos.val(i,j);
      }
      pos.push_back(one);
    }
    lattice.set_basis(pos,typ);
    return lattice;
  }), py::arg( "lattice vectors" ), py::arg("atom positions"), py::arg("atom types"), py::arg(argname.c_str())=defarg);
}

template<class T>
void declare_lattice_mat_basis_sym_init(py::class_<T,br::Lattice> &pclass){
  pclass.def(py::init( [](py::array_t<double> vecs, py::array_t<double> pypos, std::vector<unsigned long> typ, const br::Symmetry& sym){
    py::buffer_info info = vecs.request();
    if (info.ndim !=2 || info.shape[0] != 3 || info.shape[1] != 3)
      throw std::runtime_error("Basis vectors must be a 3x3 array");
    br::Array2<double> latmat = br::py2a2(vecs);
    br::Array2<double> avpos = br::py2a2(pypos);
    std::vector<std::array<double,3>> pos;
    for (br::ind_t i=0; i<avpos.size(0); ++i){
      std::array<double,3> one;
      for (br::ind_t j=0; j<3u; ++j) one[j] = avpos.val(i,j);
      pos.push_back(one);
    }
    return T(latmat, pos, typ, sym);
  }), py::arg("lattice vectors"), py::arg("atom positions"), py::arg("atom types"), py::arg("spacegroup operations"));
}

template<class T>
void declare_lattice_methods(py::class_<T,br::Lattice> &pclass, const std::string &lenunit) {
  using namespace brille;
  std::string hnum = "Hall Number";
  int hnumdef = 1;
  std::string hstr = "IT Name | Hall Symbol | Seitz notation symmetry";
  std::string hstrdef = "P_1";
    declare_lattice_scl_init(pclass,lenunit,hnum,hnumdef);
    declare_lattice_scl_init(pclass,lenunit,hstr,hstrdef);
    declare_lattice_vec_init(pclass,lenunit,hnum,hnumdef);
    declare_lattice_vec_init(pclass,lenunit,hstr,hstrdef);
    declare_lattice_mat_init(pclass,hnum,hnumdef);
    declare_lattice_mat_init(pclass,hstr,hstrdef);
    declare_lattice_mat_basis_init(pclass,hnum,hnumdef);
    declare_lattice_mat_basis_init(pclass,hstr,hstrdef);
    declare_lattice_mat_basis_sym_init(pclass);
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
