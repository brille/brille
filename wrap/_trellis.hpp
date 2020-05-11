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
#ifndef __TRELLIS_H
#define __TRELLIS_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <thread>

#include "_c_to_python.hpp"
#include "_interpolation_data.hpp"
#include "trellis.hpp"
#include "bz_trellis.hpp"
#include "utilities.hpp"


namespace py = pybind11;

typedef long slong; // ssize_t is only defined for gcc?
typedef unsigned long element_t;

template<class T,class R>
void declare_bztrellisq(py::module &m, const std::string &typestr){
  using namespace pybind11::literals;
  using Class = BrillouinZoneTrellis3<T,R>;
  std::string pyclass_name = std::string("BZTrellisQ")+typestr;
  py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
  // Initializer (BrillouinZone, maximum node volume fraction)
  .def(py::init<BrillouinZone,double>(), "brillouinzone"_a, "node_volume_fraction"_a=0.1)

  .def_property_readonly("BrillouinZone",[](const Class& cobj){return cobj.get_brillouinzone();})

  .def_property_readonly("invA",[](const Class& cobj){return av2np(cobj.get_xyz());})
  .def_property_readonly("inner_invA",[](const Class& cobj){return av2np(cobj.get_inner_xyz());})
  .def_property_readonly("outer_invA",[](const Class& cobj){return av2np(cobj.get_outer_xyz());})

  .def_property_readonly("rlu",[](const Class& cobj){return av2np(cobj.get_hkl());})
  .def_property_readonly("inner_rlu",[](const Class& cobj){return av2np(cobj.get_inner_hkl());})
  .def_property_readonly("outer_rlu",[](const Class& cobj){return av2np(cobj.get_outer_hkl());})

  .def_property_readonly("tetrahedra",[](const Class& cobj){return cobj.get_vertices_per_tetrahedron();})

  .def("fill",[](Class& cobj,
    // py::array_t<T> pyvals, py::array_t<int, py::array::c_style> pyvalelrl,
    // py::array_t<R> pyvecs, py::array_t<int, py::array::c_style> pyvecelrl
    py::array_t<T> pyvals, py::array_t<int> pyvalel,
    py::array_t<R> pyvecs, py::array_t<int> pyvecel
  ){
    ArrayVector<T> vals;
    ArrayVector<R> vecs;
    std::vector<size_t> val_sh, vec_sh;
    std::array<element_t, 3> val_el{{0,0,0}}, vec_el{{0,0,0}};
    RotatesLike val_rl, vec_rl;
    size_t count = cobj.vertex_count();
    std::tie(vals,val_sh,val_el,val_rl)=fill_check(pyvals,pyvalel,count);
    std::tie(vecs,vec_sh,vec_el,vec_rl)=fill_check(pyvecs,pyvecel,count);

    cobj.replace_value_data(vals, val_sh, val_el, val_rl);
    cobj.replace_vector_data(vecs, vec_sh, vec_el, vec_rl);
  }, "values_data"_a, "values_elements"_a, "vectors_data"_a, "vectors_elements"_a)

  //.def_property_readonly("data", /*get data*/ [](Class& cobj){ return av2np_shape(cobj.data().data(), cobj.data().shape(), false);})
  .def_property_readonly("values",[](Class& cobj){
    const auto & v{cobj.data().values()};
    return av2np_shape(v.data(), v.shape(), false);
  })
  .def_property_readonly("vectors",[](Class& cobj){
    const auto & v{cobj.data().vectors()};
    return av2np_shape(v.data(), v.shape(), false);
  })

  .def("interpolate_at",[](Class& cobj,
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
    std::tie(valres, vecres) = cobj.interpolate_at(qv, nthreads, no_move);
    // copy the results to Python arrays and return
    auto valout = iid2np(valres, cobj.data().values(),  preshape);
    auto vecout = iid2np(vecres, cobj.data().vectors(), preshape);
    return std::make_tuple(valout, vecout);
  },"Q"_a,"useparallel"_a=false,"threads"_a=-1,"do_not_move_points"_a=false)

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
    auto valout = iid2np(valres, cobj.data().values(),  preshape);
    auto vecout = iid2np(vecres, cobj.data().vectors(), preshape);
    return std::make_tuple(valout, vecout);
  },"Q"_a,"useparallel"_a=false,"threads"_a=-1,"do_not_move_points"_a=false)

  .def("debye_waller",[](Class& cobj, py::array_t<double> pyQ, py::array_t<double> pyM, double temp_k){
    // handle Q
    py::buffer_info bi = pyQ.request();
    if ( bi.shape[bi.ndim-1] !=3 )
      throw std::runtime_error("debye_waller requires one or more 3-vectors");
    // ssize_t npts = 1;
    // if (bi.ndim > 1) for (ssize_t i=0; i<bi.ndim-1; i++) npts *= bi.shape[i];
    BrillouinZone b = cobj.get_brillouinzone();
    Reciprocal lat = b.get_lattice();
    LQVec<double> cQ(lat, (double*)bi.ptr, bi.shape, bi.strides); //memcopy
    // handle the masses
    py::buffer_info mi = pyM.request();
    if ( mi.ndim != 1u )
      throw std::runtime_error("debey_waller requires masses as a 1-D vector.");
    size_t span = static_cast<size_t>(mi.strides[0])/sizeof(double);
    std::vector<double> masses(mi.shape[0]);
    double * mass_ptr = (double*) mi.ptr;
    for (size_t i=0; i<static_cast<size_t>(mi.shape[0]); ++i) masses.push_back(mass_ptr[i*span]);
    return av2np_squeeze(cobj.debye_waller(cQ, masses, temp_k));
  }, "Q"_a, "masses"_a, "Temperature_in_K"_a)

  // .def("__repr__",&Class::to_string)
  //
  .def("multi_sort_perm",
    [](Class& cobj, const double wS, const double wV,
                    const double wM, const int vwf){
    return av2np(cobj.multi_sort_perm(wS,wV,wM,vwf));
  }, "scalar_cost_weight"_a=1,
     "vector_cost_weight"_a=1,
     "matrix_cost_weight"_a=1,
     "vector_weight_function"_a=0
  )
  ;
}

#endif
