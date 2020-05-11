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
#include "bz.hpp"
#include "_c_to_python.hpp"

namespace py = pybind11;

void wrap_brillouinzone(py::module & m){
  using namespace pybind11::literals; // bring in "[name]"_a to be interpreted as py::arg("[name]")
  using CLS = BrillouinZone;
  py::class_<CLS> cls(m,"BrillouinZone");
  cls.def(py::init<Reciprocal,bool,int,bool,bool>(),
          "lattice"_a, "use_primitive"_a=true, "search_length"_a=1,
          "time_reversal_symmetry"_a=false, "wedge_search"_a=true);

  // return the internal Lattice object
  cls.def_property_readonly("lattice", [](const CLS &b){ return b.get_lattice();} );

  // access the polyhedra directly
  cls.def_property_readonly("polyhedron",&CLS::get_polyhedron);
  cls.def_property_readonly("ir_polyhedron",[](const CLS &b){return b.get_ir_polyhedron(true);});
  cls.def_property_readonly("ir_polyhedron_generated",[](const CLS &b){return b.get_ir_polyhedron(false);});

  // first Brillouin zone polyhedron
  cls.def_property_readonly("normals",                  [](const CLS &b){return av2np(b.get_normals().get_hkl());});
  cls.def_property_readonly("normals_invA",             [](const CLS &b){return av2np(b.get_normals().get_xyz());});
  cls.def_property_readonly("normals_primitive",        [](const CLS &b){return av2np(b.get_primitive_normals().get_hkl());});
  cls.def_property_readonly("points",                   [](const CLS &b){return av2np(b.get_points().get_hkl());});
  cls.def_property_readonly("points_invA",              [](const CLS &b){return av2np(b.get_points().get_xyz());});
  cls.def_property_readonly("points_primitive",         [](const CLS &b){return av2np(b.get_primitive_points().get_hkl());});
  cls.def_property_readonly("vertices",                 [](const CLS &b){return av2np(b.get_vertices().get_hkl());});
  cls.def_property_readonly("vertices_invA",            [](const CLS &b){return av2np(b.get_vertices().get_xyz());});
  cls.def_property_readonly("half_edge_points",         [](const CLS &b){return av2np(b.get_half_edges().get_hkl());});
  cls.def_property_readonly("half_edge_points_invA",    [](const CLS &b){return av2np(b.get_half_edges().get_xyz());});
  cls.def_property_readonly("vertices_primitive",       [](const CLS &b){return av2np(b.get_primitive_vertices().get_hkl());});
  cls.def_property_readonly("faces_per_vertex",  &CLS::get_faces_per_vertex);
  cls.def_property_readonly("vertices_per_face", &CLS::get_vertices_per_face);

  // irreducible first Brillouin zone polyhedron
  cls.def_property_readonly("ir_normals",                  [](const CLS &b){return av2np(b.get_ir_normals().get_hkl());});
  cls.def_property_readonly("ir_normals_invA",             [](const CLS &b){return av2np(b.get_ir_normals().get_xyz());});
  cls.def_property_readonly("ir_normals_primitive",        [](const CLS &b){return av2np(b.get_ir_primitive_normals().get_hkl());});
  cls.def_property_readonly("ir_points",                   [](const CLS &b){return av2np(b.get_ir_points().get_hkl());});
  cls.def_property_readonly("ir_points_invA",              [](const CLS &b){return av2np(b.get_ir_points().get_xyz());});
  cls.def_property_readonly("ir_points_primitive",         [](const CLS &b){return av2np(b.get_ir_primitive_points().get_hkl());});
  cls.def_property_readonly("ir_vertices",                 [](const CLS &b){return av2np(b.get_ir_vertices().get_hkl());});
  cls.def_property_readonly("ir_vertices_invA",            [](const CLS &b){return av2np(b.get_ir_vertices().get_xyz());});
  cls.def_property_readonly("ir_vertices_primitive",       [](const CLS &b){return av2np(b.get_ir_primitive_vertices().get_hkl());});
  cls.def_property_readonly("ir_faces_per_vertex",  &CLS::get_ir_faces_per_vertex);
  cls.def_property_readonly("ir_vertices_per_face", &CLS::get_ir_vertices_per_face);

  // irreducible reciprocal space wedge
  cls.def_property_readonly("wedge_normals",            [](const CLS &b){return av2np(b.get_ir_wedge_normals().get_hkl());});
  cls.def_property_readonly("wedge_normals_invA",       [](const CLS &b){return av2np(b.get_ir_wedge_normals().get_xyz());});
  cls.def_property_readonly("wedge_normals_primitive",  [](const CLS &b){return av2np(b.get_primitive_ir_wedge_normals().get_hkl());});

  // check whether one or more points are inside
  cls.def("isinside",[](CLS &b, py::array_t<double> p){
    py::buffer_info bi = p.request();
    ssize_t ndim = bi.ndim;
    if ( bi.shape[ndim-1] !=3) throw std::runtime_error("one or more 3-dimensional points is required");
    ssize_t npts = 1;
    if (ndim > 1)  for (ssize_t i=0; i<ndim-1; i++) npts *= bi.shape[i];
    LQVec<double> pv( b.get_lattice(), (double*) bi.ptr, bi.shape, bi.strides); // this is a copy :(
    ArrayVector<bool> resultv = b.isinside(pv);
    std::vector<ssize_t> outshape(ndim>1 ? ndim-1 : 1);
    if (ndim > 1){
      for (ssize_t i=0; i<ndim-1; i++) outshape[i] = bi.shape[i];
    } else {
      outshape[0] = 1;
    }
    auto result = py::array_t<bool, py::array::c_style>(outshape);
    bool *rptr = (bool *) result.request().ptr;
    for (ssize_t i=0; i<npts; i++) rptr[i] = resultv.getvalue(i);
    return result;
  },"points"_a);

  cls.def("moveinto",[](CLS &b, py::array_t<double> Q){
    py::buffer_info bi = Q.request();
    ssize_t ndim=bi.ndim;
    if (bi.shape[ndim-1] !=3) throw std::runtime_error("one or more 3-dimensional Q points is required");
    ssize_t npts = 1;
    if (ndim > 1) for (ssize_t i=0; i<ndim-1; i++) npts *= bi.shape[i];
    LQVec<double> Qv( b.get_lattice(),(double*) bi.ptr, bi.shape, bi.strides); // memcopy
    LQVec<double> qv(b.get_lattice(), npts); // output
    LQVec<int>  tauv(b.get_lattice(), npts); // output
    bool success = b.moveinto(Qv,qv,tauv);
    if (!success) throw std::runtime_error("failed to move all Q into the first Brillouin Zone");
    auto qout = py::array_t<double, py::array::c_style>(bi.shape);
    auto tout = py::array_t<int, py::array::c_style>(bi.shape);
    double *qptr = (double *) qout.request().ptr;
    int *tptr = (int *) tout.request().ptr;
    for (ssize_t i=0; i<npts; ++i){
      for (size_t j=0; j<3u; ++j){
        qptr[3u*i+j] = qv.getvalue(i,j);
        tptr[3u*i+j] = tauv.getvalue(i,j);
      }
    }
    return py::make_tuple(qout,tout);
  }, "Q"_a);

  cls.def("ir_moveinto",[](CLS &b, py::array_t<double> Q){
    py::buffer_info bi = Q.request();
    ssize_t ndim = bi.ndim;
    if (bi.shape[ndim-1] != 3)
      throw std::runtime_error("One or more 3-dimensional Q points are required.");
    ssize_t npts = 1;
    if (ndim > 1) for (ssize_t i=0; i<ndim-1; ++i) npts *= bi.shape[i];
    LQVec<double> Qv(b.get_lattice(), (double*) bi.ptr, bi.shape, bi.strides);
    // prepare intermediate outputs
    LQVec<double> qv(b.get_lattice(), npts);
    LQVec<int>  tauv(b.get_lattice(), npts);
    std::vector<size_t> rotidx(npts), invrotidx(npts);
    if (!b.ir_moveinto(Qv, qv, tauv, rotidx, invrotidx))
      throw std::runtime_error("Moving points into irreducible zone failed.");
    // get the pointgroup symmetry operations indexed by rotidx and invrotidx
    PointSymmetry ptsym = b.get_pointgroup_symmetry();
    // prepare Python outputs
    auto qout = py::array_t<double, py::array::c_style>(bi.shape);
    auto tout = py::array_t<int,    py::array::c_style>(bi.shape);
    // The rotations array has an extra dimension compared to q and tau
    bi.shape.push_back(3);
    auto rout = py::array_t<int,    py::array::c_style>(bi.shape);
    auto invrout = py::array_t<int, py::array::c_style>(bi.shape);
    // grab pointers to the underlying data blocks
    double *qptr = (double *) qout.request().ptr;
    int *tptr = (int *) tout.request().ptr;
    int *rptr = (int *) rout.request().ptr;
    int *iptr = (int *) invrout.request().ptr;
    for (ssize_t i=0; i<npts; ++i)
    for (size_t j=0; j<3u; ++j){
      qptr[3u*i+j] = qv.getvalue(i,j);
      tptr[3u*i+j] = tauv.getvalue(i,j);
      for (size_t k=0; k<3u; ++k) {
        rptr[9u*i+3u*j+k] = ptsym.get(rotidx[i])[3u*j+k];
        iptr[9u*i+3u*j+k] = ptsym.get(invrotidx[i])[3u*j+k];
      }
    }
    return py::make_tuple(qout, tout, rout, invrout);
  }, "Q"_a);

  cls.def("ir_moveinto_wedge",[](CLS &b, py::array_t<double> Q){
    py::buffer_info bi = Q.request();
    ssize_t ndim = bi.ndim;
    if (bi.shape[ndim-1] != 3)
      throw std::runtime_error("One or more 3-dimensional Q points are required.");
    ssize_t npts = 1;
    if (ndim > 1) for (ssize_t i=0; i<ndim-1; ++i) npts *= bi.shape[i];
    LQVec<double> Qv(b.get_lattice(), (double*) bi.ptr, bi.shape, bi.strides);
    // prepare intermediate outputs
    LQVec<double> qv(b.get_lattice(), npts);
    std::vector<std::array<int,9>> rots(npts);
    if (!b.ir_moveinto_wedge(Qv, qv, rots))
      throw std::runtime_error("Moving points into irreducible zone failed.");
    // prepare Python outputs
    auto qout = py::array_t<double, py::array::c_style>(bi.shape);
    // The rotations array has an extra dimension compared to q and tau
    bi.shape.push_back(3);
    auto rout = py::array_t<int,    py::array::c_style>(bi.shape);
    // grab pointers to the underlying data blocks
    double *qptr = (double *) qout.request().ptr;
    int *rptr = (int *) rout.request().ptr;
    for (ssize_t i=0; i<npts; ++i)
    for (size_t j=0; j<3u; ++j){
      qptr[3u*i+j] = qv.getvalue(i,j);
      for (size_t k=0; k<3u; ++k) rptr[9u*i+3u*j+k] = rots[i][3u*j+k];
    }
    return py::make_tuple(qout, rout);
  }, "Q"_a);
}
