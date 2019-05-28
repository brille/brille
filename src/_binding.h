/*! \file */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

#include "bz_grid.h"

#ifndef __BINDING_H
#define __BINDING_H

namespace py = pybind11;
using namespace pybind11::literals; // bring in "[name]"_a to be interpreted as py::arg("[name]")

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
void declare_bzgridq(py::module &m, const std::string &typestr) {
    using Class = BrillouinZoneGrid3<T>;
    std::string pyclass_name = std::string("BZGridQ") + typestr;
    py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
    // Initializer (BrillouinZone, [half-]Number_of_steps vector)
    .def(py::init([](BrillouinZone &b, py::array_t<size_t,py::array::c_style> pyN){
      py::buffer_info bi = pyN.request();
      if (bi.ndim != 1) throw std::runtime_error("halfN must be a 1-D array");
      if (bi.shape[0] < 3) throw std::runtime_error("halfN must have three elements");
      Class cobj( b, (size_t*)bi.ptr );
      return cobj;
    }),py::arg("brillouinzone"),py::arg("halfN"))
    // Initializer (BrillouinZone, step_size vector, flag_for_whether_step_size_is_in_rlu_or_inverse_angstrom)
    .def(py::init([](BrillouinZone &b, py::array_t<double,py::array::c_style> pyD, const bool& isrlu){
      py::buffer_info bi = pyD.request();
      if (bi.ndim != 1) throw std::runtime_error("stepsize must be a 1-D array");
      if (bi.shape[0] < 3) throw std::runtime_error("stepsize must have three elements");
      return Class( b, (double*)bi.ptr, isrlu ? 1 : 0 );
    }),py::arg("brillouinzone"),py::arg("step"),py::arg("rlu")=true)
    .def_property_readonly("N",[](const Class& cobj){ return av2np_squeeze(cobj.get_N());})
    .def_property_readonly("halfN",[](const Class& cobj){ return av2np_squeeze((cobj.get_N()-1)/2);})
    .def_property_readonly("BrillouinZone",[](const Class& cobj){ return cobj.get_brillouinzone();} )
    .def_property_readonly("rlu",[](const Class& cobj){ return av2np(cobj.get_grid_hkl());} )
    .def_property_readonly("invA",[](const Class& cobj){ return av2np(cobj.get_grid_xyz());} )
    .def_property_readonly("mapped_rlu",[](const Class& cobj){ return av2np(cobj.get_mapped_hkl());} )
    .def_property_readonly("mapped_invA",[](const Class& cobj){ return av2np(cobj.get_mapped_xyz());} )
    .def("fill",[](Class& cobj, py::array_t<T,py::array::c_style> pydata){
      py::buffer_info bi = pydata.request();
      ssize_t ndim = bi.ndim;
      /* ndim  assumed interpolation data type  numel /
      /   1              scalar                   1   /
      /   2              vector                shape[1]  /
      /   3       matrix / rank 2 tensor       shape[1]*shape[2]       /
      /   N         rank N-1 tensor            prod(shape,1,N-1)      */
      size_t numel=1, numarrays=bi.shape[0];
      if (ndim > 1) for (ssize_t i=1; i<ndim; ++i) numel *= bi.shape[i];
      ArrayVector<T> data(numel, numarrays, (T*)bi.ptr);
      ArrayVector<size_t> shape(1,ndim);
      for (ssize_t i=0; i<ndim; ++i) shape.insert(bi.shape[i], (size_t)i );
      int mapExceedsNewData = cobj.check_map(data);
      if (mapExceedsNewData) {
        std::string msg = "Provided " + std::to_string(data.size())
                        + " data inputs but expected "
                        + std::to_string(cobj.maximum_mapping()) + "!";
        throw std::runtime_error(msg);
      }
      cobj.replace_data(data,shape); // no error, so this will work for sure
      // return mapExceedsNewData; // let the calling function deal with this?
    })
    .def_property_readonly("sort_perm",[](Class& cobj){  return av2np(cobj.sort_perm()); })
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
      /*set map*/ [](Class& cobj, py::array_t<slong,py::array::c_style> pymap){
        py::buffer_info bi = pymap.request();
        if (bi.ndim != 3) throw std::runtime_error("The mapping must be a 3 dimensional array");
        for (size_t i=0; i<3; ++i) if (bi.shape[i]!=cobj.size(i))
          throw std::runtime_error("The new map shape must match the old map"); // or we could resize it, but that is more difficult
        if (cobj.maximum_mapping( (slong*)bi.ptr ) > cobj.num_data() )
          throw std::runtime_error("The largest integer in the new mapping exceeds the number of data elements.");
        cobj.unsafe_set_map( (slong*)bi.ptr ); //no error, so this works.
    })
    .def("interpolate_at",[](Class& cobj, py::array_t<double,py::array::c_style> pyX, const bool& moveinto, const bool& useparallel, const int& threads){
      py::buffer_info bi = pyX.request();
      if ( bi.shape[bi.ndim-1] !=3 )
        throw std::runtime_error("Interpolation requires one or more 3-vectors");
      ssize_t npts = 1;
      if (bi.ndim > 1) for (ssize_t i=0; i<bi.ndim-1; i++) npts *= bi.shape[i];
      BrillouinZone b = cobj.get_brillouinzone();
      Reciprocal lat = b.get_lattice();
      LQVec<double> qv(lat,npts, (double*)bi.ptr ); //memcopy
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
        for (ssize_t i=1; i<data_shape.size(); ++i) outshape.push_back( data_shape.getvalue(i) );
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
    },py::arg("Q"),py::arg("moveinto")=true,py::arg("useparallel")=false,py::arg("threads")=-1);
}

template<class T>
void declare_bzgridqe(py::module &m, const std::string &typestr) {
    using Class = BrillouinZoneGrid4<T>;
    std::string pyclass_name = std::string("BZGridQE") + typestr;
    py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
    // Initializer (BrillouinZone, [half-]Number_of_steps vector)
    .def(py::init([](BrillouinZone &b, py::array_t<double,py::array::c_style> pySpec, py::array_t<size_t,py::array::c_style> pyN){
      py::buffer_info b0 = pySpec.request();
      if (b0.ndim!=1) throw std::runtime_error("spec must be a 1-D array");
      if (b0.shape[0]<3) throw std::runtime_error("spec must have three elements");
      py::buffer_info bi = pyN.request();
      if (bi.ndim != 1) throw std::runtime_error("halfN must be a 1-D array");
      if (bi.shape[0] < 3) throw std::runtime_error("halfN must have three elements");
      Class cobj( b, (double*)b0.ptr, (size_t*)bi.ptr );
      return cobj;
    }),py::arg("brillouinzone"),py::arg("spec"),py::arg("halfN"))
    // Initializer (BrillouinZone, step_size vector, flag_for_whether_step_size_is_in_rlu_or_inverse_angstrom)
    .def(py::init([](BrillouinZone &b, py::array_t<double,py::array::c_style> pySpec, py::array_t<double,py::array::c_style> pyD, const bool& isrlu){
      py::buffer_info b0 = pySpec.request();
      if (b0.ndim!=1) throw std::runtime_error("spec must be a 1-D array");
      if (b0.shape[0]<3) throw std::runtime_error("spec must have three elements");
      py::buffer_info bi = pyD.request();
      if (bi.ndim != 1) throw std::runtime_error("stepsize must be a 1-D array");
      if (bi.shape[0] < 3) throw std::runtime_error("stepsize must have three elements");
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
    .def("fill",[](Class& cobj, py::array_t<T,py::array::c_style> pydata){
      py::buffer_info bi = pydata.request();
      ssize_t ndim = bi.ndim;
      size_t numel=1, numarrays=bi.shape[0];
      if (ndim > 1) for (ssize_t i=1; i<ndim; ++i) numel *= bi.shape[i];
      ArrayVector<T> data(numel, numarrays, (T*)bi.ptr);
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
        if (cobj.maximum_mapping( (slong*)bi.ptr ) > cobj.num_data() )
          throw std::runtime_error("The largest integer in the new mapping exceeds the number of data elements.");
        cobj.unsafe_set_map( (slong*)bi.ptr ); //no error, so this works.
    })
    .def("interpolate_at",[](Class& cobj, py::array_t<double,py::array::c_style> pyX, const bool& moveinto, const bool& useparallel, const int& threads){
      py::buffer_info bi = pyX.request();
      if ( bi.shape[bi.ndim-1] !=4 )
        throw std::runtime_error("Interpolation requires one or more 4-vectors");
      ssize_t npts = 1;
      if (bi.ndim > 1) for (ssize_t i=0; i<bi.ndim-1; i++) npts *= bi.shape[i];
      BrillouinZone b = cobj.get_brillouinzone();
      Reciprocal lat = b.get_lattice();
      ArrayVector<double> qEv(4u,npts, (double*)bi.ptr ); //memcopy
      if (moveinto){
        LQVec<double> Qv(lat, qEv, 0); // 0 ==> truncate qEv such that it has numel()==3.
        LQVec<double> qv(lat,npts); // filled by moveinto
        LQVec<int>  tauv(lat,npts); // filled by moveinto
        bool success = b.moveinto(Qv,qv,tauv);
        if (!success)
          throw std::runtime_error("failed to move all Q into the first Brillouin Zone");
        // replace the first three elements of qEv with qv in absolute units.
        ArrayVector<double> qv_invA = qv.get_xyz();
        for (size_t i=0; i<npts; ++i) for(size_t j=0; j<3u; ++j) qEv.insert(qv_invA.getvalue(i,j), i,j);
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
        for (ssize_t i=1; i<data_shape.size(); ++i) outshape.push_back( data_shape.getvalue(i) );
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
    },py::arg("QE"),py::arg("moveinto")=true,py::arg("useparallel")=false,py::arg("threads")=-1);
  }

template<class R, class T>
void declare_lattice_init(py::class_<T,Lattice> &pclass, const std::string &lenunit, const std::string &argname, const R defarg){
  pclass.def(py::init( [](py::array_t<double> lens, py::array_t<double> angs, const R groupid){
    py::buffer_info linfo = lens.request(), ainfo = angs.request();
    if ( linfo.ndim!=1 || ainfo.ndim!=1)
      throw std::runtime_error("Number of dimensions must be one");
    if ( linfo.shape[0] < 3 || ainfo.shape[0] < 3 )
      throw std::runtime_error("(At least) three lengths and angles required.");
    double *lengths = (double *) linfo.ptr, *angles = (double *) ainfo.ptr;
    return T(lengths,angles,groupid);
  }), py::arg( ("lengths/"+lenunit).c_str() ),
      py::arg("angles/radian"), py::arg(argname.c_str())=defarg);
}

template<class T>
void declare_lattice_methods(py::class_<T,Lattice> &pclass, const std::string &lenunit) {
    declare_lattice_init(pclass,lenunit,"hall",1);
    declare_lattice_init(pclass,lenunit,"ITname","P_1");
    pclass.def("star",&T::star,"Return the dual lattice");
    // pclass.def(py::init( [](py::array_t<double> lens, py::array_t<double> angs, const int hall){
    //   py::buffer_info linfo = lens.request(), ainfo = angs.request();
    //   if ( linfo.ndim!=1 || ainfo.ndim!=1)
    //     throw std::runtime_error("Number of dimensions must be one");
    //   if ( linfo.shape[0] < 3 || ainfo.shape[0] < 3 )
    //     throw std::runtime_error("(At least) three lengths and angles required.");
    //   double *lengths = (double *) linfo.ptr, *angles = (double *) ainfo.ptr;
    //   return T(lengths,angles,hall);
    // }), py::arg( ("lengths/"+lenunit).c_str() ),
    //     py::arg("angles/radian"), py::arg("hall")=1);
    pclass.def(py::init<double,double,double, double,double,double, int>(),
               py::arg( ("a/"+lenunit).c_str() ),
               py::arg( ("b/"+lenunit).c_str() ),
               py::arg( ("c/"+lenunit).c_str() ),
               py::arg("alpha/radian")=PI/2,
               py::arg("beta/radian")=PI/2,
               py::arg("gamma/radian")=PI/2,
               py::arg("HallNumber")=1 );
    pclass.def("xyz_transform",[](T &d){
      auto result = py::array_t<double, py::array::c_style>({3,3});
      py::buffer_info bi = result.request();
      d.get_xyz_transform( (double *) bi.ptr);
      return result;
    });
    pclass.def("lattice_matrix",[](T &d){
      auto result = py::array_t<double, py::array::c_style>({3,3});
      py::buffer_info bi = result.request();
      d.get_lattice_matrix( (double *) bi.ptr);
      return result;
    });
    pclass.def("isstar",(bool (T::*)(const Direct    ) const) &T::isstar);
    pclass.def("isstar",(bool (T::*)(const Reciprocal) const) &T::isstar);
    pclass.def("primitive",&T::primitive,"Return the primitive lattice");
}

#endif
