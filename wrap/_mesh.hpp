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

#include "_c_to_python.hpp"
#include "bz_mesh.hpp"
#include "utilities.hpp"


namespace py = pybind11;
typedef long slong; // ssize_t is only defined for gcc?

template<class T>
void declare_bzmeshq(py::module &m, const std::string &typestr){
  using namespace pybind11::literals;
  using Class = BrillouinZoneMesh3<T>;
  std::string pyclass_name = std::string("BZMeshQ")+typestr;
  py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
  // Initializer (BrillouinZone, max-volume, is-volume-rlu)
  .def(py::init<BrillouinZone,double,int,int>(), "brillouinzone"_a, "max_size"_a=-1., "num_levels"_a=3, "max_points"_a=-1)
  .def_property_readonly("BrillouinZone",[](const Class& cobj){return cobj.get_brillouinzone();})
  .def_property_readonly("rlu",[](const Class& cobj){return av2np(cobj.get_mesh_hkl());})
  .def_property_readonly("invA",[](const Class& cobj){return av2np(cobj.get_mesh_xyz());})
  .def_property_readonly("tetrahedra",[](const Class& cobj){return av2np(cobj.get_mesh_tetrehedra());})
  .def("fill",[](Class& cobj, py::array_t<T> pydata, py::array_t<int, py::array::c_style> pyel){
    py::buffer_info bi;
    // copy-in the elements array
    bi = pyel.request();
    if (bi.ndim != 1) throw std::runtime_error("elements must be a 1-D array");
    std::array<size_t, 3> el{{0,0,0}};
    int* intel = (int*)bi.ptr;
    for (ssize_t i=0; i<bi.shape[0] && i<3; ++i) el[i] = static_cast<size_t>(intel[i]);
    // copy-in the data ArrayVector
    bi = pydata.request();
    ArrayVector<T> data((T*)bi.ptr, bi.shape, bi.strides);
    //
    if (cobj.size() != data.size()){
      std::string msg = "Provided " + std::to_string(data.size())
                      + " data inputs but expected "
                      + std::to_string(cobj.size()) + "!";
      throw std::runtime_error(msg);
    }
    std::vector<size_t> shape;
    for (auto s: bi.shape) shape.push_back(static_cast<size_t>(s));
    cobj.replace_data(data, shape, el);
  }, "data"_a, "elements"_a)
  .def_property_readonly("data", /*get data*/ [](Class& cobj){
    auto & id = cobj.data(); // the InterpolationData object
    return av2np_shape(id.data(), id.shape(), false);})
  .def("multi_sort_perm",
    [](Class& cobj, const double wS, const double wV,
                    const double wM, const int vwf){
    return av2np(cobj.multi_sort_perm(wS,wV,wM,vwf));
  }, "scalar_cost_weight"_a=1,
     "vector_cost_weight"_a=1,
     "matrix_cost_weight"_a=1,
     "vector_weight_function"_a=0
  )
  .def("ir_interpolate_at",[](Class& cobj,
                           py::array_t<double> pyX,
                           const bool& useparallel,
                           const int& threads, const bool& no_move){
    py::buffer_info bi = pyX.request();
    if ( bi.shape[bi.ndim-1] !=3 )
      throw std::runtime_error("Interpolation requires one or more 3-vectors");
    // ssize_t npts = 1;
    // if (bi.ndim > 1) for (ssize_t i=0; i<bi.ndim-1; i++) npts *= bi.shape[i];
    BrillouinZone b = cobj.get_brillouinzone();
    Reciprocal lat = b.get_lattice();
    LQVec<double> qv(lat,(double*)bi.ptr, bi.shape, bi.strides); //memcopy

    // perform the interpolation and rotate and vectors/tensors afterwards
    int nthreads = (useparallel) ? ((threads < 1) ? static_cast<int>(std::thread::hardware_concurrency()) : threads) : 1;
    ArrayVector<T> lires = cobj.interpolate_at(qv, nthreads, no_move);
    // and then make sure we return an numpy array of appropriate size:
    std::vector<ssize_t> outshape;
    for (ssize_t i=0; i < bi.ndim-1; ++i) outshape.push_back(bi.shape[i]);
    if (cobj.data().shape().size() > 1){
      auto data_shape = cobj.data().shape();
      // the shape of each element is data_shape[1,...,data_ndim-1]
      for (size_t i=1; i<data_shape.size(); ++i) outshape.push_back(data_shape[i]);
    }
    // size_t total_elements = 1;
    // for (ssize_t osi : outshape) total_elements *= osi;
    size_t total_elements = signed_to_unsigned<size_t>(std::accumulate(outshape.begin(), outshape.end(), 1, std::multiplies<ssize_t>()));
    if (total_elements != lires.numel()*lires.size() ){
      std::cout << "outshape: (";
      for (ssize_t osi : outshape) std::cout << osi << "," ;
      std::cout << "\b)" << std::endl;
      std::string msg = "Why do we expect " + std::to_string(total_elements)
                      + " but only get back " + std::to_string(lires.numel()*lires.size())
                      + " (" + std::to_string(lires.numel())
                      + " × " + std::to_string(lires.size()) + ")";
      throw std::runtime_error(msg);
    }
    auto liout = py::array_t<T,py::array::c_style>(outshape);
    T *rptr = (T*) liout.request().ptr;
    for (size_t i =0; i< lires.size(); i++)
      for (size_t j=0; j< lires.numel(); j++)
        rptr[i*lires.numel()+j] = lires.getvalue(i,j);
    return liout;
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
  .def("__repr__",&Class::to_string);
}

#endif
