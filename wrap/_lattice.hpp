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
#include <utility>

#include "_array.hpp"
#include "_c_to_python.hpp"
#include "lattice_dual.hpp"
#include "utilities.hpp"


#ifndef WRAP_BRILLE_LATTICE_HPP_
#define WRAP_BRILLE_LATTICE_HPP_
namespace py = pybind11;
namespace br = brille;
namespace lt = brille::lattice;

template<class R, class T>
void declare_lattice_vec_init(py::class_<T> &pclass, const std::string &lenunit, const std::string &argname, const R defarg){
  pclass.def(py::init( [](const py::array_t<double>& lens, const py::array_t<double>& angs, const R groupid){
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
void declare_lattice_mat_init(py::class_<T> &pclass, const std::string &argname, const R defarg){
  pclass.def(py::init( [](const py::array_t<double>& vecs, const R groupid){
    py::buffer_info info = vecs.request();
    if ( info.ndim!=2 )
      throw std::runtime_error("Number of dimensions must be two");
    if ( info.shape[0] != 3 || info.shape[1] != 3 )
      throw std::runtime_error("Three three-vectors required.");
    return T((double *) info.ptr, info.strides, groupid);
  }), py::arg( "lattice vectors" ), py::arg(argname.c_str())=defarg);
}

template<class R, class T>
void declare_lattice_mat_basis_init(py::class_<T> &pclass, const std::string &argname, const R defarg){
  pclass.def(py::init( [](const py::array_t<double>& vecs, py::array_t<double> pypos, std::vector<br::ind_t> typ, const R groupid){
    py::buffer_info info = vecs.request();
    if ( info.ndim!=2 )
      throw std::runtime_error("Number of dimensions must be two");
    if ( info.shape[0] != 3 || info.shape[1] != 3 )
      throw std::runtime_error("Three three-vectors required.");
    T lattice((double *) info.ptr, info.strides, groupid);
    // pull in the basis information
    // br::Array<double> avpos = br::py2a(pypos);
    br::Array2<double> avpos = br::py2a2(std::move(pypos));
    std::vector<std::array<double,3>> pos;
    for (br::ind_t i=0; i<avpos.size(0); ++i){
      std::array<double,3> one{{0,0,0}};
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
void declare_lattice_mat_basis_sym_init(py::class_<T> &pclass){
  pclass.def(py::init( [](const py::array_t<double>& vecs, py::array_t<double> pypos, std::vector<br::ind_t> typ, const br::Symmetry& sym){
    py::buffer_info info = vecs.request();
    if (info.ndim !=2 || info.shape[0] != 3 || info.shape[1] != 3)
      throw std::runtime_error("Basis vectors must be a 3x3 array");
    br::Array2<double> latmat = br::py2a2(vecs);
    br::Array2<double> avpos = br::py2a2(std::move(pypos));
    std::vector<std::array<double,3>> pos;
    for (br::ind_t i=0; i<avpos.size(0); ++i){
      std::array<double,3> one{{0,0,0}};
      for (br::ind_t j=0; j<3u; ++j) one[j] = avpos.val(i,j);
      pos.push_back(one);
    }
    return T(latmat, pos, typ, sym);
  }), py::arg("lattice vectors"), py::arg("atom positions"), py::arg("atom types"), py::arg("spacegroup operations"));
}

template<class T>
void declare_lattice_methods(py::class_<T> &pclass, const std::string &lenunit) {
  using namespace brille;
  std::string hstr = "IT Name | Hall Symbol | Seitz notation symmetry";
  std::string hstrdef = "P_1";
    declare_lattice_scl_init(pclass,lenunit,hstr,hstrdef);
    declare_lattice_vec_init(pclass,lenunit,hstr,hstrdef);
    declare_lattice_mat_init(pclass,hstr,hstrdef);
    declare_lattice_mat_basis_init(pclass,hstr,hstrdef);
    declare_lattice_mat_basis_sym_init(pclass);
    pclass.def_property_readonly("star",&T::star,R"pbdoc(
    Return the dual lattice

    A lattice described by the basis vectors :math:`\mathbf{a}_1`,
    :math:`\mathbf{a}_2`, and :math:`\mathbf{a}_3` has a unit cell volume
    given by :math:`V = \mathbf{a}_1 \cdot (\mathbf{a}_2\times\mathbf{a}_3)`.
    It's dual lattice is defined by the basis vectors :math:`\mathbf{b}_1`,
    :math:`\mathbf{b}_2`, and :math:`\mathbf{b}_3` with the relationship

    .. math::
        \mathbf{b}_i = \frac{2\pi}{V}\mathbf{a}_j\times \mathbf{a}_k

    for :math:`(i,j,k) \in \left\{(1,2,3),(2,3,1),(3,1,2)\right\}`
    and has unit cell volume :math:`V* = 2\pi/V`

    )pbdoc");
    pclass.def_property_readonly("real_space_vectors",[](T &d){
      auto xyz = d.to_xyz();
      return py::array_t<double>({3,3}, {1, 3}, transpose(xyz.data(LengthUnit::angstrom)));
    });
    pclass.def_property_readonly("reciprocal_space_vectors",[](T &d){
      auto xyz = d.to_xyz();
      return py::array_t<double>({3,3}, {1, 3}, transpose(xyz.data(LengthUnit::inverse_angstrom)));
    });
    pclass.def("isstar",(bool (T::*)(const T&) const) &T::is_same);
    pclass.def("issame",(bool (T::*)(const T&) const) &T::is_same);
    pclass.def_property_readonly("primitive",&T::primitive,R"pbdoc(
    Return a primitive (non-centred) lattice equivalent to this one

    Every centred lattice has any number of non-centred lattices which are
    equally valid in that they span the same space. This method returns the
    primitive lattice related to a centred lattice by the translation vectors
    described in :py:class:`~brille._brille.Bravais`.
    )pbdoc");
}

#endif
