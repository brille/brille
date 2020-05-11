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
#ifndef __GRID_H
#define __GRID_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <thread>

#include "_c_to_python.hpp"
#include "_interpolation_data.hpp"
#include "bz_grid.hpp"
#include "utilities.hpp"


namespace py = pybind11;
typedef long slong; // ssize_t is only defined for gcc?

template<class T, class R>
void declare_bzgridq(py::module &m, const std::string &typestr) {
    using namespace pybind11::literals;
    using Class = BrillouinZoneGrid3<T,R>;
    std::string pyclass_name = std::string("BZGridQ") + typestr;
    py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
    // Initializer (BrillouinZone, [half-]Number_of_steps vector)
    .def(py::init([](BrillouinZone &b, py::array_t<size_t> pyN){
      py::buffer_info bi = pyN.request();
      if (bi.ndim != 1) throw std::runtime_error("halfN must be a 1-D array");
      if (bi.shape[0] < 3) throw std::runtime_error("halfN must have three elements");
      if (bi.strides[0]!=sizeof(size_t)) throw std::runtime_error("halfN must be contiguous");
      Class cobj( b, (size_t*)bi.ptr );
      return cobj;
    }),"brillouinzone"_a,"halfN"_a)
    // Initializer (BrillouinZone, step_size vector, flag_for_whether_step_size_is_in_rlu_or_inverse_angstrom)
    .def(py::init([](BrillouinZone &b, py::array_t<double> pyD, const bool& isrlu){
      py::buffer_info bi = pyD.request();
      if (bi.ndim != 1) throw std::runtime_error("stepsize must be a 1-D array");
      if (bi.shape[0] < 3) throw std::runtime_error("stepsize must have three elements");
      if (bi.strides[0]!=sizeof(double)) throw std::runtime_error("stepsize must be contiguous");
      return Class( b, (double*)bi.ptr, isrlu ? 1 : 0 );
    }),"brillouinzone"_a,"step"_a,"rlu"_a=true)
    // // Initializer (BrillouinZone, tetrahedron volume, flag_for_whether_volume_is_in_rlu_or_inverse_angstrom)
    // .def(py::init([](BrillouinZone &b, const double& vol, const bool& isrlu){
    //   return Class( b, vol, isrlu ? 1 : 0 );
    // }),"brillouinzone"_a,"volume"_a,"rlu"_a=true)
    .def_property_readonly("N",[](const Class& cobj){ return av2np_squeeze(cobj.get_N());})
    .def_property_readonly("halfN",[](const Class& cobj){ return av2np_squeeze((cobj.get_N()-1)/2);})
    .def_property_readonly("BrillouinZone",[](const Class& cobj){ return cobj.get_brillouinzone();} )
    .def_property_readonly("grid_rlu",[](const Class& cobj){ return av2np(cobj.get_grid_hkl());} )
    .def_property_readonly("grid_invA",[](const Class& cobj){ return av2np(cobj.get_grid_xyz());} )
    .def_property_readonly("rlu",[](const Class& cobj){ return av2np(cobj.get_mapped_hkl());} )
    .def_property_readonly("invA",[](const Class& cobj){ return av2np(cobj.get_mapped_xyz());} )
    .def("fill",[](Class& cobj,
      py::array_t<T> pyvals, py::array_t<int, py::array::c_style> pyvalelrl,
      py::array_t<R> pyvecs, py::array_t<int, py::array::c_style> pyvecelrl
    ){
      ArrayVector<T> vals;
      ArrayVector<R> vecs;
      std::vector<size_t> val_sh, vec_sh;
      std::array<element_t, 3> val_el{{0,0,0}}, vec_el{{0,0,0}};
      RotatesLike val_rl, vec_rl;
      size_t count = cobj.valid_mapping_count();
      std::tie(vals,val_sh,val_el,val_rl)=fill_check(pyvals,pyvalelrl,count);
      std::tie(vecs,vec_sh,vec_el,vec_rl)=fill_check(pyvecs,pyvecelrl,count);

      cobj.replace_value_data(vals, val_sh, val_el, val_rl);
      cobj.replace_vector_data(vecs, vec_sh, vec_el, vec_rl);
    }, "values_data"_a, "values_elements"_a, "vectors_data"_a, "vectors_elements"_a)
    .def_property_readonly("values",[](Class& cobj){
      const auto & v{cobj.data().values()};
      return av2np_shape(v.data(), v.shape(), false);
    })
    .def_property_readonly("vectors",[](Class& cobj){
      const auto & v{cobj.data().vectors()};
      return av2np_shape(v.data(), v.shape(), false);
    })
    .def("multi_sort_perm",
      [](Class& cobj, const double wS, const double wV,
                      const double wM, const int vwf){
      return av2np(cobj.multi_sort_perm(wS,wV,wM,vwf));
    }, "scalar_cost_weight"_a=1,
       "vector_cost_weight"_a=1,
       "matrix_cost_weight"_a=1,
       "vector_weight_function"_a=0
    )
    .def("centre_sort_perm",
      [](Class& cobj, const double wS, const double wV,
                      const double wM, const int vwf){
      return av2np(cobj.centre_sort_perm(wS,wV,wM,vwf));
    }, "scalar_cost_weight"_a=1,
       "vector_cost_weight"_a=1,
       "matrix_cost_weight"_a=1,
       "vector_weight_function"_a=0
    )
    .def_property("map",
      /*get map*/ [](Class& cobj){
        std::vector<ssize_t> shape(3); // the map is 3D
        size_t expected = 1u;
        for (size_t i=0; i<3u; ++i){
          shape[i] = cobj.size(i);
          expected *= cobj.size(i);
        }
        auto np = py::array_t<slong,py::array::c_style>(shape);
        size_t nret = cobj.unsafe_get_map( (slong*)np.request().ptr );
        if (nret != expected)
          // I guess nret is smaller, otherwise we probably already encountered a segfault
          throw std::runtime_error("Something has gone horribly wrong with getting the map.");
        return np;
      },
      /*set map*/ [](Class& cobj, py::array_t<slong> pymap){
        py::buffer_info bi = pymap.request();
        if (bi.ndim != 3) throw std::runtime_error("The mapping must be a 3 dimensional array");
        for (size_t i=0; i<3; ++i) if (bi.shape[i]<0 || static_cast<size_t>(bi.shape[i])!=cobj.size(i))
          throw std::runtime_error("The new map shape must match the old map"); // or we could resize it, but that is more difficult
        for (size_t i=0; i<3; ++i) if (bi.strides[i]<0 || static_cast<size_t>(bi.strides[i])!=cobj.stride(i))
          throw std::runtime_error("The new map strides must match the old map");
        if (cobj.maximum_mapping( (slong*)bi.ptr ) > cobj.data().size() )
          throw std::runtime_error("The largest integer in the new mapping exceeds the number of data elements.");
        cobj.unsafe_set_map( (slong*)bi.ptr ); //no error, so this works.
    })
    .def("interpolate_at",[](Class& cobj,
                             py::array_t<double> pyX,
                             const bool& moveinto,
                             const bool& useparallel,
                             const int& threads){
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
      if (moveinto){
        LQVec<double> Qv(qv); // second memcopy
        LQVec<int>  tauv(lat, qv.size()); // filled by moveinto
        bool success = b.moveinto(Qv,qv,tauv);
        if (!success)
          throw std::runtime_error("failed to move all Q into the first Brillouin Zone");
      }
      // perform the interpolation and rotate and vectors/tensors afterwards
      const int maxth(static_cast<int>(std::thread::hardware_concurrency()));
      int nthreads = (useparallel) ? ((threads < 1) ? maxth : threads) : 1;
      ArrayVector<T> valres;
      ArrayVector<R> vecres;
      std::tie(valres, vecres) = useparallel ?  cobj.parallel_linear_interpolate_at(qv, nthreads) : cobj.linear_interpolate_at(qv);
      // copy the results to Python arrays and return
      auto valout = iid2np(valres, cobj.data().values(),  preshape);
      auto vecout = iid2np(vecres, cobj.data().vectors(), preshape);
      return std::make_tuple(valout, vecout);
    },"Q"_a,"moveinto"_a=true,"useparallel"_a=false,"threads"_a=-1)
    //
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
    //
    // .def("sum_data",[](Class& cobj, const int axis, const bool squeeze){
    //   return av2np_shape( cobj.sum_data(axis), cobj.data_shape(), squeeze);
    // },"axis"_a,"squeeze"_a=true)
    .def("nearest_map_index",[](Class& cobj, py::array_t<double> pyX, const bool isrlu){
      py::buffer_info xinfo = pyX.request();
      if (xinfo.ndim!=1)
        throw std::runtime_error("x must be a 1-D array");
      if (xinfo.shape[0]<3)
        throw std::runtime_error("x must have three elements");
      // Copy the python x position into an ArrayVector
      ArrayVector<double> x_invA((double*)xinfo.ptr, xinfo.shape, xinfo.strides);
      if (isrlu){
        // if x was provided in relative lattice units, conver it to Å⁻¹
        BrillouinZone bz = cobj.get_brillouinzone();
        Reciprocal rlat = bz.get_lattice();
        LQVec<double> x_rlu(rlat, x_invA);
        x_invA = x_rlu.get_xyz();
      }
      size_t subidx[3];
      unsigned int flg = cobj.nearest_index(x_invA.data(0,0), subidx);
      py::tuple ret = py::make_tuple(flg, cobj.sub2map(subidx));
      return ret;
    },"x"_a,"isrlu"_a=true)
    .def("floor_map_index",[](Class& cobj, py::array_t<double> pyX, const bool isrlu){
      py::buffer_info xinfo = pyX.request();
      if (xinfo.ndim!=1)
        throw std::runtime_error("x must be a 1-D array");
      if (xinfo.shape[0]<3)
        throw std::runtime_error("x must have three elements");
      // Copy the python x position into an ArrayVector
      ArrayVector<double> x_invA((double*)xinfo.ptr, xinfo.shape, xinfo.strides);
      if (isrlu){
        // if x was provided in relative lattice units, conver it to Å⁻¹
        BrillouinZone bz = cobj.get_brillouinzone();
        Reciprocal rlat = bz.get_lattice();
        LQVec<double> x_rlu(rlat,x_invA);
        x_invA = x_rlu.get_xyz();
      }
      size_t subidx[3];
      unsigned int flg = cobj.floor_index(x_invA.data(0,0), subidx);
      py::tuple ret = py::make_tuple(flg, cobj.sub2map(subidx));
      return ret;
    },"x"_a,"isrlu"_a=true)
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
    }, "Q"_a, "masses"_a, "Temperature_in_K"_a);
}

// template<class T>
// void declare_bzgridqe(py::module &m, const std::string &typestr) {
//     using namespace pybind11::literals;
//     using Class = BrillouinZoneGrid4<T>;
//     std::string pyclass_name = std::string("BZGridQE") + typestr;
//     py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
//     // Initializer (BrillouinZone, [half-]Number_of_steps vector)
//     .def(py::init([](BrillouinZone &b, py::array_t<double> pySpec, py::array_t<size_t> pyN){
//       py::buffer_info b0 = pySpec.request();
//       if (b0.ndim!=1) throw std::runtime_error("spec must be a 1-D array");
//       if (b0.shape[0]<3) throw std::runtime_error("spec must have three elements");
//       if (b0.strides[0]!=sizeof(double)) throw std::runtime_error("spec must be contiguous");
//       py::buffer_info bi = pyN.request();
//       if (bi.ndim != 1) throw std::runtime_error("halfN must be a 1-D array");
//       if (bi.shape[0] < 3) throw std::runtime_error("halfN must have three elements");
//       if (bi.strides[0]!=sizeof(size_t)) throw std::runtime_error("halfN must be contiguous");
//       Class cobj( b, (double*)b0.ptr, (size_t*)bi.ptr );
//       return cobj;
//     }),"brillouinzone"_a,"spec"_a,"halfN"_a)
//     // Initializer (BrillouinZone, step_size vector, flag_for_whether_step_size_is_in_rlu_or_inverse_angstrom)
//     .def(py::init([](BrillouinZone &b, py::array_t<double> pySpec, py::array_t<double> pyD, const bool& isrlu){
//       py::buffer_info b0 = pySpec.request();
//       if (b0.ndim!=1) throw std::runtime_error("spec must be a 1-D array");
//       if (b0.shape[0]<3) throw std::runtime_error("spec must have three elements");
//       if (b0.strides[0]!=sizeof(double)) throw std::runtime_error("spec must be contiguous");
//       py::buffer_info bi = pyD.request();
//       if (bi.ndim != 1) throw std::runtime_error("stepsize must be a 1-D array");
//       if (bi.shape[0] < 3) throw std::runtime_error("stepsize must have three elements");
//       if (bi.strides[0]!=sizeof(double)) throw std::runtime_error("stepsize must be contiguous");
//       return Class( b, (double*)b0.ptr, (double*)bi.ptr, isrlu ? 1 : 0 );
//     }),"brillouinzone"_a,"spec"_a,"step"_a,"rlu"_a=true)
//     .def_property_readonly("N",    [](const Class& cobj){return av2np_squeeze(cobj.get_N());} )
//     .def_property_readonly("halfN",[](const Class& cobj){return av2np_squeeze(cobj.get_halfN());} )
//     .def_property_readonly("spec", [](const Class& cobj){return av2np_squeeze(cobj.get_spec());} )
//     .def_property_readonly("BrillouinZone",[](const Class& cobj){ return cobj.get_brillouinzone();} )
//     .def_property_readonly("grid_rlu_Q", [](const Class& cobj){ return av2np(cobj.get_grid_hkl());} )
//     .def_property_readonly("grid_invA_Q",[](const Class& cobj){ return av2np(cobj.get_grid_xyz());} )
//     .def_property_readonly("rlu_Q",      [](const Class& cobj){ return av2np(cobj.get_mapped_hkl());} )
//     .def_property_readonly("invA_Q",     [](const Class& cobj){ return av2np(cobj.get_mapped_xyz());} )
//     .def_property_readonly("grid_rlu",   [](const Class& cobj){ return av2np(cobj.get_grid_hkle());} )
//     .def_property_readonly("grid_invA",  [](const Class& cobj){ return av2np(cobj.get_grid_xyzw());} )
//     .def_property_readonly("rlu",        [](const Class& cobj){ return av2np(cobj.get_mapped_hkle());} )
//     .def_property_readonly("invA",       [](const Class& cobj){ return av2np(cobj.get_mapped_xyzw());} )
//     .def("fill",[](Class& cobj, py::array_t<T> pydata){
//       py::buffer_info bi = pydata.request();
//       ArrayVector<T> data((T*)bi.ptr, bi.shape, bi.strides);
//       std::vector<size_t> shape;
//       for (auto s: bi.shape) shape.push_back(static_cast<size_t>(s));
//       int mapExceedsNewData = cobj.check_map(data);
//       if (mapExceedsNewData){
//         std::string msg = "Provided " + std::to_string(data.size())
//                         + " data inputs but expected "
//                         + std::to_string(cobj.maximum_mapping()) + "!";
//         throw std::runtime_error(msg);
//       }
//       std::array<size_t, 3> el{{0,0,0}};
//       cobj.replace_data(data,shape,el); // no error, so this will work for sure
//       // return mapExceedsNewData; // let the calling function deal with this?
//     })
//     .def_property("map",
//       /*get map*/ [](Class& cobj){
//         std::vector<ssize_t> shape(4); // the map is 4D
//         size_t expected = 1u;
//         for (size_t i=0; i<4u; ++i){
//           shape[i] = cobj.size(i);
//           expected *= cobj.size(i);
//         }
//         auto np = py::array_t<slong,py::array::c_style>(shape);
//         size_t nret = cobj.unsafe_get_map( (slong*)np.request().ptr );
//         if (nret != expected)
//           // I guess nret is smaller, otherwise we probably already encountered a segfault
//           throw std::runtime_error("Something has gone horribly wrong with getting the map.");
//         return np;
//       },
//       /*set map*/ [](Class& cobj, py::array_t<slong,py::array::c_style> pymap){
//         py::buffer_info bi = pymap.request();
//         if (bi.ndim != 4) throw std::runtime_error("The mapping must be a 4 dimensional array");
//         for (size_t i=0; i<4; ++i) if (bi.shape[i]<0 || static_cast<size_t>(bi.shape[i])!=cobj.size(i))
//           throw std::runtime_error("The new map shape must match the old map"); // or we could resize it, but that is more difficult
//         for (size_t i=0; i<4; ++i) if (bi.strides[i]<0 || static_cast<size_t>(bi.strides[i])!=cobj.stride(i))
//           throw std::runtime_error("The new map strides must match the old map");
//         if (cobj.maximum_mapping( (slong*)bi.ptr ) > cobj.data().size() )
//           throw std::runtime_error("The largest integer in the new mapping exceeds the number of data elements.");
//         cobj.unsafe_set_map( (slong*)bi.ptr ); //no error, so this works.
//     })
//     .def("interpolate_at",[](Class& cobj, py::array_t<double> pyX, const bool& moveinto, const bool& useparallel, const int& threads){
//       py::buffer_info bi = pyX.request();
//       if ( bi.shape[bi.ndim-1] !=4 )
//         throw std::runtime_error("Interpolation requires one or more 4-vectors");
//       ssize_t npts = 1;
//       if (bi.ndim > 1) for (ssize_t i=0; i<bi.ndim-1; i++) npts *= bi.shape[i];
//       BrillouinZone b = cobj.get_brillouinzone();
//       Reciprocal lat = b.get_lattice();
//       ArrayVector<double> qEv((double*)bi.ptr, bi.shape, bi.strides); //memcopy
//       if (moveinto){
//         LQVec<double> Qv(lat, qEv, 0); // 0 ==> truncate qEv such that it has numel()==3.
//         LQVec<double> qv(lat,npts); // filled by moveinto
//         LQVec<int>  tauv(lat,npts); // filled by moveinto
//         bool success = b.moveinto(Qv,qv,tauv);
//         if (!success)
//           throw std::runtime_error("failed to move all Q into the first Brillouin Zone");
//         // replace the first three elements of qEv with qv in absolute units.
//         ArrayVector<double> qv_invA = qv.get_xyz();
//         for (ssize_t i=0; i<npts; ++i) for(size_t j=0; j<3u; ++j) qEv.insert(qv_invA.getvalue(i,j), i,j);
//       }
//       // do the interpolation for each point in qv
//       ArrayVector<T> lires;
//       if (useparallel){
//         lires = cobj.parallel_linear_interpolate_at(qEv,threads);
//       } else {
//         lires = cobj.linear_interpolate_at(qEv);
//       }
//       // and then make sure we return an numpy array of appropriate size:
//       std::vector<ssize_t> outshape;
//       for (ssize_t i=0; i < bi.ndim-1; ++i) outshape.push_back(bi.shape[i]);
//       if (cobj.data().shape().size() > 1){
//         auto data_shape = cobj.data().shape();
//         // the shape of each element is data_shape[1,...,data_ndim-1]
//         for (size_t i=1; i<data_shape.size(); ++i) outshape.push_back(data_shape[i]);
//       }
//       // size_t total_elements = 1;
//       // for (ssize_t osi : outshape) total_elements *= osi;
//       size_t total_elements = signed_to_unsigned<size_t>(std::accumulate(outshape.begin(), outshape.end(), 1, std::multiplies<ssize_t>()));
//       if (total_elements != lires.numel()*lires.size() ){
//         std::cout << "outshape: (";
//         for (ssize_t osi : outshape) std::cout << osi << "," ;
//         std::cout << "\b)" << std::endl;
//         std::string msg = "Why do we expect " + std::to_string(total_elements)
//                         + " but only get back " + std::to_string(lires.numel()*lires.size())
//                         + " (" + std::to_string(lires.numel())
//                         + " × " + std::to_string(lires.size()) + ")";
//         throw std::runtime_error(msg);
//       }
//       auto liout = py::array_t<T,py::array::c_style>(outshape);
//       T *rptr = (T*) liout.request().ptr;
//       for (size_t i =0; i< lires.size(); i++)
//         for (size_t j=0; j< lires.numel(); j++)
//           rptr[i*lires.numel()+j] = lires.getvalue(i,j);
//       return liout;
//     },"QE"_a,"moveinto"_a=true,"useparallel"_a=false,"threads"_a=-1)
//     // .def("sum_data",[](Class& cobj, const int axis, const bool squeeze){
//     //   return av2np_shape( cobj.sum_data(axis), cobj.data_shape(), squeeze);
//     // },"axis"_a,"squeeze"_a=true)
//     ;
// }

#endif
