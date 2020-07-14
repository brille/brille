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
#ifndef __MESH_H
#define __MESH_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <thread>
#include <tuple>

#include "_c_to_python.hpp"
#include "_interpolation_data.hpp"
#include "bz_mesh.hpp"
#include "utilities.hpp"


namespace py = pybind11;
typedef long slong; // ssize_t is only defined for gcc?

template<class T,class R>
void declare_bzmeshq(py::module &m, const std::string &typestr){
  using namespace pybind11::literals;
  using Class = BrillouinZoneMesh3<T,R>;
  std::string pyclass_name = std::string("BZMeshQ")+typestr;
  py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
  // Initializer (BrillouinZone, max-volume, is-volume-rlu)
  .def(py::init<BrillouinZone,double,int,int>(), "brillouinzone"_a, "max_size"_a=-1., "num_levels"_a=3, "max_points"_a=-1)
  .def_property_readonly("BrillouinZone",[](const Class& cobj){return cobj.get_brillouinzone();})
  .def_property_readonly("rlu",[](const Class& cobj){return av2np(cobj.get_mesh_hkl());})
  .def_property_readonly("invA",[](const Class& cobj){return av2np(cobj.get_mesh_xyz());})
  .def_property_readonly("tetrahedra",[](const Class& cobj){return av2np(cobj.get_mesh_tetrehedra());})

  .def("fill",[](Class& cobj,
    py::array_t<T> pyvals, py::array_t<int, py::array::c_style> pyvalelrl,
    py::array_t<R> pyvecs, py::array_t<int, py::array::c_style> pyvecelrl,
    bool sort
  ){
    ArrayVector<T> vals;
    ArrayVector<R> vecs;
    std::vector<size_t> val_sh, vec_sh;
    std::array<element_t, 3> val_el{{0,0,0}}, vec_el{{0,0,0}};
    RotatesLike val_rl, vec_rl;
    size_t count = cobj.size();
    std::tie(vals,val_sh,val_el,val_rl)=fill_check(pyvals,pyvalelrl,count);
    std::tie(vecs,vec_sh,vec_el,vec_rl)=fill_check(pyvecs,pyvecelrl,count);

    cobj.replace_value_data(vals, val_sh, val_el, val_rl);
    cobj.replace_vector_data(vecs, vec_sh, vec_el, vec_rl);
    if (sort) cobj.sort();
  }, "values_data"_a, "values_elements"_a, "vectors_data"_a, "vectors_elements"_a, "sort"_a=false)

  .def("fill",[](Class& cobj,
    // py::array_t<T> pyvals, py::array_t<int, py::array::c_style> pyvalelrl,
    // py::array_t<R> pyvecs, py::array_t<int, py::array::c_style> pyvecelrl
    py::array_t<T> pyvals, py::array_t<int> pyvalel, py::array_t<double> pyvalwght,
    py::array_t<R> pyvecs, py::array_t<int> pyvecel, py::array_t<double> pyvecwght,
    bool sort
  ){
    ArrayVector<T> vals;
    ArrayVector<R> vecs;
    std::vector<size_t> val_sh, vec_sh;
    std::array<element_t, 3> val_el{{0,0,0}}, vec_el{{0,0,0}};
    std::array<double,3> val_wght{{1,1,1}}, vec_wght{{1,1,1}};
    RotatesLike val_rl, vec_rl;
    int val_sf{0}, val_vf{0}, vec_sf{0}, vec_vf{0};
    size_t count = cobj.size();
    std::tie(vals,val_sh,val_el,val_rl,val_sf,val_vf,val_wght)=fill_check(pyvals,pyvalel,pyvalwght,count);
    std::tie(vecs,vec_sh,vec_el,vec_rl,vec_sf,vec_vf,vec_wght)=fill_check(pyvecs,pyvecel,pyvecwght,count);

    cobj.replace_value_data(vals, val_sh, val_el, val_rl);
    cobj.replace_vector_data(vecs, vec_sh, vec_el, vec_rl);
    cobj.set_value_cost_info(val_sf, val_vf, val_wght);
    cobj.set_vector_cost_info(vec_sf, vec_vf, vec_wght);
    if (sort) cobj.sort();
  }, "values_data"_a, "values_elements"_a, "values_weights"_a,"vectors_data"_a, "vectors_elements"_a, "vectors_weights"_a, "sort"_a=false)


  .def_property_readonly("values",[](Class& cobj){
    const auto & v{cobj.data().values()};
    return av2np_shape(v.data(), v.shape(), false);
  })
  .def_property_readonly("vectors",[](Class& cobj){
    const auto & v{cobj.data().vectors()};
    return av2np_shape(v.data(), v.shape(), false);
  })

  .def("ir_interpolate_at",[](Class& cobj,
                           py::array_t<double> pyX,
                           const bool& useparallel,
                           const int& threads, const bool& no_move){
    py::buffer_info bi = pyX.request();
    if ( bi.shape[bi.ndim-1] !=3 )
      throw std::runtime_error("Interpolation requires one or more 3-vectors");
    // store shape of X before three-vector dimension for shaping output
    std::vector<ssize_t> preshape;
    for (ssize_t i=0; i < bi.ndim-1; ++i) preshape.push_back(bi.shape[i]);
    // copy the Python X array
    BrillouinZone b = cobj.get_brillouinzone();
    Reciprocal lat = b.get_lattice();
    LQVec<double> qv(lat,(double*)bi.ptr, bi.shape, bi.strides); //memcopy
    // perform the interpolation and rotate and vectors/tensors afterwards
    const int maxth(static_cast<int>(std::thread::hardware_concurrency()));
    int nthreads = (useparallel) ? ((threads < 1) ? maxth : threads) : 1;
    ArrayVector<T> valres;
    ArrayVector<R> vecres;
    std::tie(valres, vecres) = cobj.ir_interpolate_at(qv, nthreads, no_move);
    // copy results to Python arrays and return
    py::array_t<T, py::array::c_style> valout = iid2np(valres, cobj.data().values(),  preshape);
    py::array_t<R, py::array::c_style> vecout = iid2np(vecres, cobj.data().vectors(), preshape);
    return std::make_tuple(valout, vecout);
  },"Q"_a,"useparallel"_a=false,"threads"_a=-1,"do_not_move_points"_a=false)

  .def("ir_interpolate_at_dw",[](Class& cobj,
                           py::array_t<double> pyX,
                           py::array_t<double> pyM, double temp_k,
                           const bool& useparallel,
                           const int& threads, const bool& no_move){
    py::buffer_info bi = pyX.request();
    if ( bi.shape[bi.ndim-1] !=3 )
      throw std::runtime_error("Interpolation requires one or more 3-vectors");
    // store shape of X before three-vector dimension for shaping output
    std::vector<ssize_t> preshape;
    for (ssize_t i=0; i < bi.ndim-1; ++i) preshape.push_back(bi.shape[i]);
    // copy the Python X array
    BrillouinZone b = cobj.get_brillouinzone();
    Reciprocal lat = b.get_lattice();
    LQVec<double> qv(lat,(double*)bi.ptr, bi.shape, bi.strides); //memcopy
    // perform the interpolation and rotate and vectors/tensors afterwards
    const int maxth(static_cast<int>(std::thread::hardware_concurrency()));
    int nthreads = (useparallel) ? ((threads < 1) ? maxth : threads) : 1;
    ArrayVector<T> valres;
    ArrayVector<R> vecres;
    std::tie(valres, vecres) = cobj.ir_interpolate_at(qv, nthreads, no_move);
    // copy results to Python arrays and return
    py::array_t<T, py::array::c_style> valout = iid2np(valres, cobj.data().values(),  preshape);
    py::array_t<R, py::array::c_style> vecout = iid2np(vecres, cobj.data().vectors(), preshape);
    // calculate the Debye-Waller factor
    auto Wd = av2np_squeeze(cobj.debye_waller(qv, np2vec(pyM), temp_k));
    return std::make_tuple(valout, vecout, Wd);
  },"Q"_a,"M/amu"_a,"temperature/K"_a,"useparallel"_a=false,"threads"_a=-1,"do_not_move_points"_a=false)


  .def("debye_waller",[](Class& cobj, py::array_t<double> pyQ, py::array_t<double> pyM, double temp_k){
    // handle Q
    py::buffer_info bi = pyQ.request();
    if ( bi.shape[bi.ndim-1] !=3 )
      throw std::runtime_error("debye_waller requires one or more 3-vectors");
    BrillouinZone b = cobj.get_brillouinzone();
    Reciprocal lat = b.get_lattice();
    LQVec<double> cQ(lat, (double*)bi.ptr, bi.shape, bi.strides); //memcopy
    return av2np_squeeze(cobj.debye_waller(cQ, np2vec(pyM), temp_k));
  }, "Q"_a, "masses"_a, "Temperature_in_K"_a)
  .def("__repr__",&Class::to_string)

  .def("sort",&Class::sort)
  ;
}

#endif
