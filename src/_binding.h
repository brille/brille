/*! \file */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <thread>

#include "bz_grid.h"
#include "bz_mesh.h"

#ifndef __BINDING_H
#define __BINDING_H

namespace py = pybind11;
using namespace pybind11::literals; // bring in "[name]"_a to be interpreted as py::arg("[name]")

template<typename T, size_t N> py::array_t<T> sa2np(const std::vector<ssize_t>& sz, const std::array<T,N>& sv){
  size_t numel = 1;
  for (ssize_t i: sz) numel *= i;
  if (N != numel){
    std::string msg = "Inconsistent required shape ( ";
    for (ssize_t i: sz) msg += std::to_string(i) + " ";
    msg += ") and array size " + std::to_string(N);
    throw std::runtime_error(msg);
  }
  auto np = py::array_t<T,py::array::c_style>(sz);
  T *ptr = (T*) np.request().ptr;
  for (size_t i=0; i<numel; ++i) ptr[i] = sv[i];
  return np;
}
template<typename T> py::array_t<T> sv2np(const std::vector<ssize_t>& sz, const std::vector<T>& sv){
  size_t numel = 1;
  for (ssize_t i: sz) numel *= i;
  if (sv.size() != numel){
    std::string msg = "Inconsistent required shape ( ";
    for (ssize_t i: sz) msg += std::to_string(i) + " ";
    msg += ") and vector size " + std::to_string(sv.size());
    throw std::runtime_error(msg);
  }
  auto np = py::array_t<T,py::array::c_style>(sz);
  T *ptr = (T*) np.request().ptr;
  for (size_t i=0; i<numel; ++i) ptr[i] = sv[i];
  return np;
}
template<typename T, size_t N>
py::array_t<T> sva2np(const std::vector<ssize_t>&sz,
                      const std::vector<std::array<T,N>>& sva)
{
  size_t numel = 1;
  for (ssize_t i: sz) numel *= i;
  if (sva.size()*N != numel){
    std::string msg = "Inconsistent required shape ( ";
    for (ssize_t i: sz) msg += std::to_string(i) + " ";
    msg += ") and vector<array<" + std::to_string(N) + "> size ";
    msg += std::to_string(sva.size()) + ".";
    throw std::runtime_error(msg);
  }
  auto np = py::array_t<T, py::array::c_style>(sz);
  T *ptr = (T*) np.request().ptr;
  size_t i=0;
  for (std::array<T,N> array: sva) for (T value: array) ptr[i++]=value;
  return np;
}

template<typename T> py::array_t<T> av2np(const ArrayVector<T>& av){
  std::vector<ssize_t> shape(2); // ArrayVectors are 2D by default
  shape[0] = av.size();
  shape[1] = av.numel();
  auto np = py::array_t<T,py::array::c_style>(shape);
  T *rptr = (T*) np.request().ptr;
  for (size_t i =0; i< av.size(); i++)
    for (size_t j=0; j< av.numel(); j++)
      rptr[i*av.numel()+j] = av.getvalue(i,j);
  return np;
}
template<typename T> py::array_t<T> av2np_squeeze(const ArrayVector<T>& av){
  std::vector<ssize_t> shape(1); // assume we'll be able to squeeze.
  if (av.size()==1u)
    shape[0] = av.numel();
  else{
    if (av.numel()==1u)
      shape[0] = av.size();
    else
      return av2np(av); // no squeezing possible
  }
  auto np = py::array_t<T,py::array::c_style>(shape);
  T *rptr = (T*) np.request().ptr;
  if (av.size()==1u)
    for (size_t j=0; j<av.numel(); ++j)
      rptr[j] = av.getvalue(0,j);
  else
    for (size_t i=0; i<av.size(); ++i)
      rptr[i] = av.getvalue(i,0);
  return np;
}
/*! \brief Convert an ArrayVector to a numpy.ndarray, using a defined shape

The ArrayVector type is always two-dimensional, but is used by MapGrid3 and
MapGrid4 to represent higher-dimensional data stored at each mapped grid point.
When passing back gridded information to Python, we can re-create the shape
information in the created numpy.ndarray.
@param av The (numel(),size()) ArrayVector
@param inshape A (1,1+array_ndim) ArrayVector consisting of
              (n_arrays, array_dim₀, array_dim₁, …, array_dimₙ)
@param squeeze A flag to indicate if, for av.size()==1, whether to return a
               (1, array_dim₀, array_dim₁, ..., array_dimₙ) or a
               (array_dim₀, array_dim₁, ..., array_dimₙ) numpy.ndarray.
               Or to remove all singleton dimensions from the output.
@returns A pybind wrapped numpy.ndarray with shape
        (av.size(), array_dim₀, array_dim₁, …, array_dimₙ)
*/
template<typename T>
py::array_t<T> av2np_shape(const ArrayVector<T>& av,
                           const ArrayVector<size_t>& inshape,
                           const bool squeeze = false){
  std::vector<ssize_t> outshape;
  if (!(squeeze && av.size()==1))
    outshape.push_back( (ssize_t)av.size() );
  for (size_t i=1; i<inshape.size(); ++i)
    if (!(squeeze && inshape.getvalue(i)==1))
      outshape.push_back( (ssize_t)inshape.getvalue(i) );
  size_t numel = 1;
  for (ssize_t osi : outshape) numel *= (size_t)osi;
  if (numel != av.size()*av.numel()){
    return squeeze ? av2np_squeeze(av) : av2np(av);
    std::string msg = "Expected " + std::to_string(numel)
                    + " but only have " + std::to_string(av.numel()*av.size())
                    + " (" + std::to_string(av.numel())
                    + " × " + std::to_string(av.size()) + ")";
    throw std::runtime_error(msg);
  }
  auto out = py::array_t<T,py::array::c_style>(outshape);
  T *ptr = (T*)out.request().ptr;
  for (size_t i=0; i<av.size(); i++)
    for (size_t j=0; j<av.numel(); j++)
      ptr[i*av.numel()+j] = av.getvalue(i,j);
  return out;
}

std::string long_version(){
  using namespace symbz::version;
  std::string v = version_number;
  if (!std::string(git_revision).empty()){
    v += "-" + std::string(git_revision).substr(0,7)
       + "@" + git_branch;
  }
  return v;
}

typedef long slong; // ssize_t is only defined for gcc?

template<class T>
void declare_bzmeshq(py::module &m, const std::string &typestr){
  using Class = BrillouinZoneMesh3<T>;
  std::string pyclass_name = std::string("BZMeshQ")+typestr;
  py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
  // Initializer (BrillouinZone, max-volume, is-volume-rlu)
  .def(py::init<BrillouinZone,double,double,double,double,int>(),
       py::arg("brillouinzone"),
       py::arg("max_size")=-1.,
       py::arg("min_angle")=20.,
       py::arg("max_angle")=-1.,
       py::arg("min_ratio")=-1.,
       py::arg("max_points")=-1
      )
  .def_property_readonly("BrillouinZone",[](const Class& cobj){return cobj.get_brillouinzone();})
  .def_property_readonly("rlu",[](const Class& cobj){return av2np(cobj.get_mesh_hkl());})
  .def_property_readonly("invA",[](const Class& cobj){return av2np(cobj.get_mesh_xyz());})
  .def_property_readonly("tetrahedra",[](const Class& cobj){return av2np(cobj.get_mesh_tetrehedra());})
  .def("fill",[](Class& cobj, py::array_t<T> pydata, py::array_t<int, py::array::c_style> pyel){
    py::buffer_info bi;
    // copy-in the elements array
    bi = pyel.request();
    if (bi.ndim != 1) throw std::runtime_error("elements must be a 1-D array");
    std::array<unsigned, 4> el{0,0,0,0};
    int* intel = (int*)bi.ptr;
    for (ssize_t i=0; i<bi.shape[0] && i<4; ++i) el[i] = static_cast<unsigned>(intel[i]);
    // copy-in the data ArrayVector
    bi = pydata.request();
    ssize_t ndim = bi.ndim;
    ArrayVector<T> data((T*)bi.ptr, bi.shape, bi.strides);
    //
    if (cobj.size() != data.size()){
      std::string msg = "Provided " + std::to_string(data.size())
                      + " data inputs but expected "
                      + std::to_string(cobj.size()) + "!";
      throw std::runtime_error(msg);
    }
    ArrayVector<size_t> shape(1, ndim);
    for (ssize_t i=0; i<ndim; ++i) shape.insert(bi.shape[i], static_cast<size_t>(i));
    cobj.replace_data(data, shape, el);
  }, py::arg("data"), py::arg("elements"))
  .def_property_readonly("data", /*get data*/ [](Class& cobj){ return av2np_shape(cobj.get_data(), cobj.data_shape(), false);})
  .def("interpolate_at",[](Class& cobj,
                           py::array_t<double> pyX,
                           const bool& useparallel,
                           const int& threads){
    py::buffer_info bi = pyX.request();
    if ( bi.shape[bi.ndim-1] !=3 )
      throw std::runtime_error("Interpolation requires one or more 3-vectors");
    ssize_t npts = 1;
    if (bi.ndim > 1) for (ssize_t i=0; i<bi.ndim-1; i++) npts *= bi.shape[i];
    BrillouinZone b = cobj.get_brillouinzone();
    Reciprocal lat = b.get_lattice();
    LQVec<double> qv(lat,(double*)bi.ptr, bi.shape, bi.strides); //memcopy

    int nthreads = (useparallel) ? ((threads < 1) ? static_cast<int>(std::thread::hardware_concurrency()) : threads) : 1;
    // perform the interpolation and rotate and vectors/tensors afterwards
    ArrayVector<T> lires = cobj.interpolate_at(qv, nthreads);
    // and then make sure we return an numpy array of appropriate size:
    std::vector<ssize_t> outshape;
    for (ssize_t i=0; i < bi.ndim-1; ++i) outshape.push_back(bi.shape[i]);
    if (cobj.data_ndim() > 1){
      ArrayVector<size_t> data_shape = cobj.data_shape();
      // the shape of each element is data_shape[1,...,data_ndim-1]
      for (size_t i=1; i<data_shape.size(); ++i) outshape.push_back( data_shape.getvalue(i) );
    }
    size_t total_elements = 1;
    for (ssize_t osi : outshape) total_elements *= osi;
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
  },py::arg("Q"),py::arg("useparallel")=false,py::arg("threads")=-1)
  .def("debye_waller",[](Class& cobj, py::array_t<double> pyQ, py::array_t<double> pyM, double temp_k){
    // handle Q
    py::buffer_info bi = pyQ.request();
    if ( bi.shape[bi.ndim-1] !=3 )
      throw std::runtime_error("debye_waller requires one or more 3-vectors");
    ssize_t npts = 1;
    if (bi.ndim > 1) for (ssize_t i=0; i<bi.ndim-1; i++) npts *= bi.shape[i];
    BrillouinZone b = cobj.get_brillouinzone();
    Reciprocal lat = b.get_lattice();
    LQVec<double> cQ(lat, (double*)bi.ptr, bi.shape, bi.strides); //memcopy
    // handle the masses
    py::buffer_info mi = pyM.request();
    if ( mi.ndim != 1u )
      throw std::runtime_error("debey_waller requires masses as a 1-D vector.");
    size_t span = mi.strides[0]/sizeof(double);
    std::vector<double> masses(mi.shape[0]);
    double * mass_ptr = (double*) mi.ptr;
    for (size_t i=0; i<mi.shape[0]; ++i) masses.push_back(mass_ptr[i*span]);
    return av2np_squeeze(cobj.debye_waller(cQ, masses, temp_k));
  }, py::arg("Q"), py::arg("masses"), py::arg("Temperature_in_K"));
}

template<class T>
void declare_bzgridq(py::module &m, const std::string &typestr) {
    using Class = BrillouinZoneGrid3<T>;
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
    }),py::arg("brillouinzone"),py::arg("halfN"))
    // Initializer (BrillouinZone, step_size vector, flag_for_whether_step_size_is_in_rlu_or_inverse_angstrom)
    .def(py::init([](BrillouinZone &b, py::array_t<double> pyD, const bool& isrlu){
      py::buffer_info bi = pyD.request();
      if (bi.ndim != 1) throw std::runtime_error("stepsize must be a 1-D array");
      if (bi.shape[0] < 3) throw std::runtime_error("stepsize must have three elements");
      if (bi.strides[0]!=sizeof(double)) throw std::runtime_error("stepsize must be contiguous");
      return Class( b, (double*)bi.ptr, isrlu ? 1 : 0 );
    }),py::arg("brillouinzone"),py::arg("step"),py::arg("rlu")=true)
    // // Initializer (BrillouinZone, tetrahedron volume, flag_for_whether_volume_is_in_rlu_or_inverse_angstrom)
    // .def(py::init([](BrillouinZone &b, const double& vol, const bool& isrlu){
    //   return Class( b, vol, isrlu ? 1 : 0 );
    // }),py::arg("brillouinzone"),py::arg("volume"),py::arg("rlu")=true)
    .def_property_readonly("N",[](const Class& cobj){ return av2np_squeeze(cobj.get_N());})
    .def_property_readonly("halfN",[](const Class& cobj){ return av2np_squeeze((cobj.get_N()-1)/2);})
    .def_property_readonly("BrillouinZone",[](const Class& cobj){ return cobj.get_brillouinzone();} )
    .def_property_readonly("rlu",[](const Class& cobj){ return av2np(cobj.get_grid_hkl());} )
    .def_property_readonly("invA",[](const Class& cobj){ return av2np(cobj.get_grid_xyz());} )
    .def_property_readonly("mapped_rlu",[](const Class& cobj){ return av2np(cobj.get_mapped_hkl());} )
    .def_property_readonly("mapped_invA",[](const Class& cobj){ return av2np(cobj.get_mapped_xyz());} )
    .def("fill",[](Class& cobj, py::array_t<T> pydata, py::array_t<int, py::array::c_style> pyel){
      py::buffer_info bi;
      // copy-in the elements array
      bi = pyel.request();
      if (bi.ndim != 1) throw std::runtime_error("elements must be a 1-D array");
      std::array<unsigned, 4> el{0,0,0,0};
      int* intel = (int*)bi.ptr;
      for (ssize_t i=0; i<bi.shape[0] && i<4; ++i) el[i] = static_cast<unsigned>(intel[i]);
      // copy-in the data ArrayVector
      bi = pydata.request();
      ssize_t ndim = bi.ndim;
      ArrayVector<T> data((T*)bi.ptr, bi.shape, bi.strides);
      if (cobj.check_map(data) /*non-zero indicates too-many data*/){
        std::string msg = "Provided " + std::to_string(data.size())
                        + " data inputs but expected "
                        + std::to_string(cobj.maximum_mapping()) + "!";
        throw std::runtime_error(msg);
      }
      ArrayVector<size_t> shape(1, ndim);
      for (ssize_t i=0; i<ndim; ++i) shape.insert(bi.shape[i], static_cast<size_t>(i));
      cobj.replace_data(data, shape, el);
    }, py::arg("data"), py::arg("elements"))
    .def_property_readonly("data",
      /*get data*/ [](Class& cobj){
        return av2np_shape(cobj.get_data(), cobj.data_shape(), false);
      }
      // ,
      // /*set data*/[](Class& cobj, py::array_t<T,py::array::c_style> pydata){
      //   py::buffer_info bi = pydata.request();
      //   ssize_t ndim = bi.ndim;
      //   size_t numel=1, numarrays=bi.shape[0];
      //   if (ndim > 1) for (ssize_t i=1; i<ndim; ++i) numel *= bi.shape[i];
      //   ArrayVector<T> data(numel, numarrays, (T*)bi.ptr);
      //   ArrayVector<size_t> shape(1,ndim);
      //   for (ssize_t i=0; i<ndim; ++i) shape.insert(bi.shape[i], (size_t)i );
      //   int mapExceedsNewData = cobj.check_map(data);
      //   if (mapExceedsNewData) {
      //     std::string msg = "Provided " + std::to_string(data.size())
      //                     + " data inputs but expected "
      //                     + std::to_string(cobj.maximum_mapping()) + "!";
      //     throw std::runtime_error(msg);
      //   }
      //   cobj.replace_data(data,shape); // no error, so this will work for sure
      // }
    )
    .def("multi_sort_perm",
      [](Class& cobj, const double wS, const double wE, const double wV,
                      const double wM, const int ewf, const int vwf){
      return av2np(cobj.multi_sort_perm(wS,wE,wV,wM,ewf,vwf));
    }, py::arg("scalar_cost_weight")=1,
       py::arg("eigenvector_cost_weight")=1,
       py::arg("vector_cost_weight")=1,
       py::arg("matrix_cost_weight")=1,
       py::arg("eigenvector_weight_function")=0,
       py::arg("vector_weight_function")=0
    )
    .def("centre_sort_perm",
      [](Class& cobj, const double wS, const double wE, const double wV,
                      const double wM, const int ewf, const int vwf){
      return av2np(cobj.centre_sort_perm(wS,wE,wV,wM,ewf,vwf));
    }, py::arg("scalar_cost_weight")=1,
       py::arg("eigenvector_cost_weight")=1,
       py::arg("vector_cost_weight")=1,
       py::arg("matrix_cost_weight")=1,
       py::arg("eigenvector_weight_function")=0,
       py::arg("vector_weight_function")=0
    )
    .def_property("map",
      /*get map*/ [](Class& cobj){
        std::vector<ssize_t> shape(3); // the map is 3D
        shape[0] = cobj.size(0);
        shape[1] = cobj.size(1);
        shape[2] = cobj.size(2);
        auto np = py::array_t<slong,py::array::c_style>(shape);
        size_t nret = cobj.unsafe_get_map( (slong*)np.request().ptr );
        if (nret != shape[0]*shape[1]*shape[2])
          // I guess nret is smaller, otherwise we probably already encountered a segfault
          throw std::runtime_error("Something has gone horribly wrong with getting the map.");
        return np;
      },
      /*set map*/ [](Class& cobj, py::array_t<slong> pymap){
        py::buffer_info bi = pymap.request();
        if (bi.ndim != 3) throw std::runtime_error("The mapping must be a 3 dimensional array");
        for (size_t i=0; i<3; ++i) if (bi.shape[i]!=cobj.size(i))
          throw std::runtime_error("The new map shape must match the old map"); // or we could resize it, but that is more difficult
        for (size_t i=0; i<3; ++i) if (bi.strides[i]!=cobj.stride(i))
          throw std::runtime_error("The new map strides must match the old map");
        if (cobj.maximum_mapping( (slong*)bi.ptr ) > cobj.num_data() )
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
      ssize_t npts = 1;
      if (bi.ndim > 1) for (ssize_t i=0; i<bi.ndim-1; i++) npts *= bi.shape[i];
      BrillouinZone b = cobj.get_brillouinzone();
      Reciprocal lat = b.get_lattice();
      LQVec<double> qv(lat,(double*)bi.ptr, bi.shape, bi.strides); //memcopy
      if (moveinto){
        LQVec<double> Qv(qv); // second memcopy
        LQVec<int>  tauv(lat,npts); // filled by moveinto
        bool success = b.moveinto(Qv,qv,tauv);
        if (!success)
          throw std::runtime_error("failed to move all Q into the first Brillouin Zone");
      }
      // do the interpolation for each point in qv
      ArrayVector<T> lires;
      if (useparallel){
        lires = cobj.parallel_linear_interpolate_at(qv,threads);
      } else {
        lires = cobj.linear_interpolate_at(qv);
      }
      // and then make sure we return an numpy array of appropriate size:
      std::vector<ssize_t> outshape;
      for (ssize_t i=0; i < bi.ndim-1; ++i) outshape.push_back(bi.shape[i]);
      if (cobj.data_ndim() > 1){
        ArrayVector<size_t> data_shape = cobj.data_shape();
        // the shape of each element is data_shape[1,...,data_ndim-1]
        for (size_t i=1; i<data_shape.size(); ++i) outshape.push_back( data_shape.getvalue(i) );
      }
      size_t total_elements = 1;
      for (ssize_t osi : outshape) total_elements *= osi;
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
    },py::arg("Q"),py::arg("moveinto")=true,py::arg("useparallel")=false,py::arg("threads")=-1)
    //
    .def("ir_interpolate_at",[](Class& cobj,
                             py::array_t<double> pyX,
                             const bool& useparallel,
                             const int& threads, const bool& no_move){
      py::buffer_info bi = pyX.request();
      if ( bi.shape[bi.ndim-1] !=3 )
        throw std::runtime_error("Interpolation requires one or more 3-vectors");
      ssize_t npts = 1;
      if (bi.ndim > 1) for (ssize_t i=0; i<bi.ndim-1; i++) npts *= bi.shape[i];
      BrillouinZone b = cobj.get_brillouinzone();
      Reciprocal lat = b.get_lattice();
      LQVec<double> qv(lat,(double*)bi.ptr, bi.shape, bi.strides); //memcopy

      int nthreads = (useparallel) ? ((threads < 1) ? static_cast<int>(std::thread::hardware_concurrency()) : threads) : 1;
      // perform the interpolation and rotate and vectors/tensors afterwards
      ArrayVector<T> lires = cobj.ir_interpolate_at(qv, nthreads, no_move);
      // and then make sure we return an numpy array of appropriate size:
      std::vector<ssize_t> outshape;
      for (ssize_t i=0; i < bi.ndim-1; ++i) outshape.push_back(bi.shape[i]);
      if (cobj.data_ndim() > 1){
        ArrayVector<size_t> data_shape = cobj.data_shape();
        // the shape of each element is data_shape[1,...,data_ndim-1]
        for (size_t i=1; i<data_shape.size(); ++i) outshape.push_back( data_shape.getvalue(i) );
      }
      size_t total_elements = 1;
      for (ssize_t osi : outshape) total_elements *= osi;
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
    },py::arg("Q"),py::arg("useparallel")=false,py::arg("threads")=-1,py::arg("do_not_move_points")=false)
    //
    .def("sum_data",[](Class& cobj, const int axis, const bool squeeze){
      return av2np_shape( cobj.sum_data(axis), cobj.data_shape(), squeeze);
    },py::arg("axis"),py::arg("squeeze")=true)
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
      unsigned int flg = cobj.nearest_index(x_invA.datapointer(0,0), subidx);
      py::tuple ret = py::make_tuple(flg, cobj.sub2map(subidx));
      return ret;
    },py::arg("x"),py::arg("isrlu")=true)
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
      unsigned int flg = cobj.floor_index(x_invA.datapointer(0,0), subidx);
      py::tuple ret = py::make_tuple(flg, cobj.sub2map(subidx));
      return ret;
    },py::arg("x"),py::arg("isrlu")=true)
    .def("debye_waller",[](Class& cobj, py::array_t<double> pyQ, py::array_t<double> pyM, double temp_k){
      // handle Q
      py::buffer_info bi = pyQ.request();
      if ( bi.shape[bi.ndim-1] !=3 )
        throw std::runtime_error("debye_waller requires one or more 3-vectors");
      ssize_t npts = 1;
      if (bi.ndim > 1) for (ssize_t i=0; i<bi.ndim-1; i++) npts *= bi.shape[i];
      BrillouinZone b = cobj.get_brillouinzone();
      Reciprocal lat = b.get_lattice();
      LQVec<double> cQ(lat, (double*)bi.ptr, bi.shape, bi.strides); //memcopy
      // handle the masses
      py::buffer_info mi = pyM.request();
      if ( mi.ndim != 1u )
        throw std::runtime_error("debey_waller requires masses as a 1-D vector.");
      size_t span = mi.strides[0]/sizeof(double);
      std::vector<double> masses(mi.shape[0]);
      double * mass_ptr = (double*) mi.ptr;
      for (size_t i=0; i<mi.shape[0]; ++i) masses.push_back(mass_ptr[i*span]);
      return av2np_squeeze(cobj.debye_waller(cQ, masses, temp_k));
    }, py::arg("Q"), py::arg("masses"), py::arg("Temperature_in_K"));
}

template<class T>
void declare_bzgridqe(py::module &m, const std::string &typestr) {
    using Class = BrillouinZoneGrid4<T>;
    std::string pyclass_name = std::string("BZGridQE") + typestr;
    py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
    // Initializer (BrillouinZone, [half-]Number_of_steps vector)
    .def(py::init([](BrillouinZone &b, py::array_t<double> pySpec, py::array_t<size_t> pyN){
      py::buffer_info b0 = pySpec.request();
      if (b0.ndim!=1) throw std::runtime_error("spec must be a 1-D array");
      if (b0.shape[0]<3) throw std::runtime_error("spec must have three elements");
      if (b0.strides[0]!=sizeof(double)) throw std::runtime_error("spec must be contiguous");
      py::buffer_info bi = pyN.request();
      if (bi.ndim != 1) throw std::runtime_error("halfN must be a 1-D array");
      if (bi.shape[0] < 3) throw std::runtime_error("halfN must have three elements");
      if (bi.strides[0]!=sizeof(size_t)) throw std::runtime_error("halfN must be contiguous");
      Class cobj( b, (double*)b0.ptr, (size_t*)bi.ptr );
      return cobj;
    }),py::arg("brillouinzone"),py::arg("spec"),py::arg("halfN"))
    // Initializer (BrillouinZone, step_size vector, flag_for_whether_step_size_is_in_rlu_or_inverse_angstrom)
    .def(py::init([](BrillouinZone &b, py::array_t<double> pySpec, py::array_t<double> pyD, const bool& isrlu){
      py::buffer_info b0 = pySpec.request();
      if (b0.ndim!=1) throw std::runtime_error("spec must be a 1-D array");
      if (b0.shape[0]<3) throw std::runtime_error("spec must have three elements");
      if (b0.strides[0]!=sizeof(double)) throw std::runtime_error("spec must be contiguous");
      py::buffer_info bi = pyD.request();
      if (bi.ndim != 1) throw std::runtime_error("stepsize must be a 1-D array");
      if (bi.shape[0] < 3) throw std::runtime_error("stepsize must have three elements");
      if (bi.strides[0]!=sizeof(double)) throw std::runtime_error("stepsize must be contiguous");
      return Class( b, (double*)b0.ptr, (double*)bi.ptr, isrlu ? 1 : 0 );
    }),py::arg("brillouinzone"),py::arg("spec"),py::arg("step"),py::arg("rlu")=true)
    .def_property_readonly("N",    [](const Class& cobj){return av2np_squeeze(cobj.get_N());} )
    .def_property_readonly("halfN",[](const Class& cobj){return av2np_squeeze(cobj.get_halfN());} )
    .def_property_readonly("spec", [](const Class& cobj){return av2np_squeeze(cobj.get_spec());} )
    .def_property_readonly("BrillouinZone",[](const Class& cobj){ return cobj.get_brillouinzone();} )
    .def_property_readonly("rlu_Q",[](const Class& cobj){ return av2np(cobj.get_grid_hkl());} )
    .def_property_readonly("invA_Q",[](const Class& cobj){ return av2np(cobj.get_grid_xyz());} )
    .def_property_readonly("mapped_rlu_Q",[](const Class& cobj){ return av2np(cobj.get_mapped_hkl());} )
    .def_property_readonly("mapped_invA_Q",[](const Class& cobj){ return av2np(cobj.get_mapped_xyz());} )
    .def_property_readonly("rlu",[](const Class& cobj){ return av2np(cobj.get_grid_hkle());} )
    .def_property_readonly("invA",[](const Class& cobj){ return av2np(cobj.get_grid_xyzw());} )
    .def_property_readonly("mapped_rlu",[](const Class& cobj){ return av2np(cobj.get_mapped_hkle());} )
    .def_property_readonly("mapped_invA",[](const Class& cobj){ return av2np(cobj.get_mapped_xyzw());} )
    .def("fill",[](Class& cobj, py::array_t<T> pydata){
      py::buffer_info bi = pydata.request();
      ssize_t ndim = bi.ndim;
      ArrayVector<T> data((T*)bi.ptr, bi.shape, bi.strides);
      ArrayVector<size_t> shape(1,ndim);
      for (ssize_t i=0; i<ndim; ++i) shape.insert(bi.shape[i], (size_t)i );
      int mapExceedsNewData = cobj.check_map(data);
      if (mapExceedsNewData){
        std::string msg = "Provided " + std::to_string(data.size())
                        + " data inputs but expected "
                        + std::to_string(cobj.maximum_mapping()) + "!";
        throw std::runtime_error(msg);
      }
      cobj.replace_data(data,shape); // no error, so this will work for sure
      // return mapExceedsNewData; // let the calling function deal with this?
    })
    .def_property("map",
      /*get map*/ [](Class& cobj){
        std::vector<ssize_t> shape(4); // the map is 4D
        for (int i=0; i<4; i++) shape[i] = cobj.size(i);
        auto np = py::array_t<slong,py::array::c_style>(shape);
        size_t nret = cobj.unsafe_get_map( (slong*)np.request().ptr );
        if (nret != shape[0]*shape[1]*shape[2]*shape[3])
          // I guess nret is smaller, otherwise we probably already encountered a segfault
          throw std::runtime_error("Something has gone horribly wrong with getting the map.");
        return np;
      },
      /*set map*/ [](Class& cobj, py::array_t<slong,py::array::c_style> pymap){
        py::buffer_info bi = pymap.request();
        if (bi.ndim != 4) throw std::runtime_error("The mapping must be a 4 dimensional array");
        for (size_t i=0; i<4; ++i) if (bi.shape[i]!=cobj.size(i))
          throw std::runtime_error("The new map shape must match the old map"); // or we could resize it, but that is more difficult
        for (size_t i=0; i<4; ++i) if (bi.strides[i]!=cobj.stride(i))
          throw std::runtime_error("The new map strides must match the old map");
        if (cobj.maximum_mapping( (slong*)bi.ptr ) > cobj.num_data() )
          throw std::runtime_error("The largest integer in the new mapping exceeds the number of data elements.");
        cobj.unsafe_set_map( (slong*)bi.ptr ); //no error, so this works.
    })
    .def("interpolate_at",[](Class& cobj, py::array_t<double> pyX, const bool& moveinto, const bool& useparallel, const int& threads){
      py::buffer_info bi = pyX.request();
      if ( bi.shape[bi.ndim-1] !=4 )
        throw std::runtime_error("Interpolation requires one or more 4-vectors");
      ssize_t npts = 1;
      if (bi.ndim > 1) for (ssize_t i=0; i<bi.ndim-1; i++) npts *= bi.shape[i];
      BrillouinZone b = cobj.get_brillouinzone();
      Reciprocal lat = b.get_lattice();
      ArrayVector<double> qEv((double*)bi.ptr, bi.shape, bi.strides); //memcopy
      if (moveinto){
        LQVec<double> Qv(lat, qEv, 0); // 0 ==> truncate qEv such that it has numel()==3.
        LQVec<double> qv(lat,npts); // filled by moveinto
        LQVec<int>  tauv(lat,npts); // filled by moveinto
        bool success = b.moveinto(Qv,qv,tauv);
        if (!success)
          throw std::runtime_error("failed to move all Q into the first Brillouin Zone");
        // replace the first three elements of qEv with qv in absolute units.
        ArrayVector<double> qv_invA = qv.get_xyz();
        for (ssize_t i=0; i<npts; ++i) for(size_t j=0; j<3u; ++j) qEv.insert(qv_invA.getvalue(i,j), i,j);
      }
      // do the interpolation for each point in qv
      ArrayVector<T> lires;
      if (useparallel){
        lires = cobj.parallel_linear_interpolate_at(qEv,threads);
      } else {
        lires = cobj.linear_interpolate_at(qEv);
      }
      // and then make sure we return an numpy array of appropriate size:
      std::vector<ssize_t> outshape;
      for (ssize_t i=0; i < bi.ndim-1; ++i) outshape.push_back(bi.shape[i]);
      if (cobj.data_ndim() > 1){
        ArrayVector<size_t> data_shape = cobj.data_shape();
        // the shape of each element is data_shape[1,...,data_ndim-1]
        for (size_t i=1; i<data_shape.size(); ++i) outshape.push_back( data_shape.getvalue(i) );
      }
      size_t total_elements = 1;
      for (ssize_t osi : outshape) total_elements *= osi;
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
    },py::arg("QE"),py::arg("moveinto")=true,py::arg("useparallel")=false,py::arg("threads")=-1)
    .def("sum_data",[](Class& cobj, const int axis, const bool squeeze){
      return av2np_shape( cobj.sum_data(axis), cobj.data_shape(), squeeze);
    },py::arg("axis"),py::arg("squeeze")=true);
}

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
    pclass.def("isstar",(bool (T::*)(const Direct    ) const) &T::isstar);
    pclass.def("isstar",(bool (T::*)(const Reciprocal) const) &T::isstar);
    pclass.def_property_readonly("primitive",&T::primitive,"Return the primitive lattice");
}

template<class T>
void declare_polyhedron(py::module &m, const std::string &typestr) {
    std::string pyclass_name = std::string("Polyhedron") + typestr;
    py::class_<T>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
    .def_property_readonly("vertices",[](const T& o){return av2np(o.get_vertices());})
    .def_property_readonly("points",[](const T& o){return av2np(o.get_points());})
    .def_property_readonly("normals",[](const T& o){return av2np(o.get_normals());})
    .def_property_readonly("vertices_per_face",&T::get_vertices_per_face)
    .def_property_readonly("faces_per_vertex",&T::get_faces_per_vertex)
    .def_property_readonly("volume",&T::get_volume)
    .def_property_readonly("mirror",&T::mirror)
    .def("rotate",[](const T& o, py::array_t<double> rot){
      py::buffer_info info = rot.request();
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
#endif
