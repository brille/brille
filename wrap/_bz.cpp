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
#include "bz.hpp"
#include "_c_to_python.hpp"
#include "_array.hpp"

namespace py = pybind11;

void wrap_brillouinzone(py::module & m){
  using namespace pybind11::literals; // bring in "[name]"_a to be interpreted as py::arg("[name]")
  using namespace brille;
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
  cls.def_property_readonly("normals",                  [](const CLS &b){return brille::a2py(b.get_normals().get_hkl());});
  cls.def_property_readonly("normals_invA",             [](const CLS &b){return brille::a2py(b.get_normals().get_xyz());});
  cls.def_property_readonly("normals_primitive",        [](const CLS &b){return brille::a2py(b.get_primitive_normals().get_hkl());});
  cls.def_property_readonly("points",                   [](const CLS &b){return brille::a2py(b.get_points().get_hkl());});
  cls.def_property_readonly("points_invA",              [](const CLS &b){return brille::a2py(b.get_points().get_xyz());});
  cls.def_property_readonly("points_primitive",         [](const CLS &b){return brille::a2py(b.get_primitive_points().get_hkl());});
  cls.def_property_readonly("vertices",                 [](const CLS &b){return brille::a2py(b.get_vertices().get_hkl());});
  cls.def_property_readonly("vertices_invA",            [](const CLS &b){return brille::a2py(b.get_vertices().get_xyz());});
  cls.def_property_readonly("half_edge_points",         [](const CLS &b){return brille::a2py(b.get_half_edges().get_hkl());});
  cls.def_property_readonly("half_edge_points_invA",    [](const CLS &b){return brille::a2py(b.get_half_edges().get_xyz());});
  cls.def_property_readonly("vertices_primitive",       [](const CLS &b){return brille::a2py(b.get_primitive_vertices().get_hkl());});
  cls.def_property_readonly("faces_per_vertex",  &CLS::get_faces_per_vertex);
  cls.def_property_readonly("vertices_per_face", &CLS::get_vertices_per_face);

  // irreducible first Brillouin zone polyhedron
  cls.def_property_readonly("ir_normals",                  [](const CLS &b){return brille::a2py(b.get_ir_normals().get_hkl());});
  cls.def_property_readonly("ir_normals_invA",             [](const CLS &b){return brille::a2py(b.get_ir_normals().get_xyz());});
  cls.def_property_readonly("ir_normals_primitive",        [](const CLS &b){return brille::a2py(b.get_ir_primitive_normals().get_hkl());});
  cls.def_property_readonly("ir_points",                   [](const CLS &b){return brille::a2py(b.get_ir_points().get_hkl());});
  cls.def_property_readonly("ir_points_invA",              [](const CLS &b){return brille::a2py(b.get_ir_points().get_xyz());});
  cls.def_property_readonly("ir_points_primitive",         [](const CLS &b){return brille::a2py(b.get_ir_primitive_points().get_hkl());});
  cls.def_property_readonly("ir_vertices",                 [](const CLS &b){return brille::a2py(b.get_ir_vertices().get_hkl());});
  cls.def_property_readonly("ir_vertices_invA",            [](const CLS &b){return brille::a2py(b.get_ir_vertices().get_xyz());});
  cls.def_property_readonly("ir_vertices_primitive",       [](const CLS &b){return brille::a2py(b.get_ir_primitive_vertices().get_hkl());});
  cls.def_property_readonly("ir_faces_per_vertex",  &CLS::get_ir_faces_per_vertex);
  cls.def_property_readonly("ir_vertices_per_face", &CLS::get_ir_vertices_per_face);

  // irreducible reciprocal space wedge
  cls.def_property_readonly("wedge_normals",            [](const CLS &b){return brille::a2py(b.get_ir_wedge_normals().get_hkl());});
  cls.def_property_readonly("wedge_normals_invA",       [](const CLS &b){return brille::a2py(b.get_ir_wedge_normals().get_xyz());});
  cls.def_property_readonly("wedge_normals_primitive",  [](const CLS &b){return brille::a2py(b.get_primitive_ir_wedge_normals().get_hkl());});

  // check whether one or more points are inside
  cls.def("isinside",[](CLS &b, py::array_t<double> p){
    // brille::Array<double> sp = brille::py2a(p);
    brille::Array2<double> sp = brille::py2a2(p);
    if (sp.size(sp.ndim()-1) != 3)
      throw std::runtime_error("The last dimension must have size 3");
    LQVec<double> pv( b.get_lattice(), sp); // no copy :)
    return b.isinside(pv);
  },"points"_a, R"pbdoc(
    Determine whether each of the provided reciprocal lattice points is located
    within the first Brillouin zone

    Parameters
    ----------
    Q : :py:class:`numpy.ndarray`
      A 2+ dimensional array of three-vectors (`Q.shape[-1]==3`) expressed in
      units of the reciprocal lattice.

    Returns
    -------
    :py:class:`numpy.ndarray`
      Logical array with one less dimension than `Q` and shape `Q.shape[0:-1]`
  )pbdoc");

  cls.def("moveinto",[](CLS &b, py::array_t<double> Q, int threads){
    // brille::Array<double> sp = brille::py2a(Q);
    brille::Array2<double> sp = brille::py2a2(Q);
    if (sp.size(sp.ndim()-1) != 3)
      throw std::runtime_error("The last dimension must have size 3");
    LQVec<double> Qv( b.get_lattice(),  sp); // view
    LQVec<double> qv(b.get_lattice(), sp.shape(), sp.stride()); // output
    LQVec<int>  tauv(b.get_lattice(), sp.shape(), sp.stride()); // output
    bool success = b.moveinto(Qv,qv,tauv,threads);
    if (!success) throw std::runtime_error("failed to move all Q into the first Brillouin Zone");
    return py::make_tuple(brille::a2py(qv), brille::a2py(tauv));
  }, "Q"_a, "threads"_a=0, R"pbdoc(
    Find points equivalent to those provided within the first Brillouin zone.

    Parameters
    ----------
    Q : :py:class:`numpy.ndarray`
      A 2+ dimensional array of three-vectors (`Q.shape[-1]==3`) expressed in
      units of the reciprocal lattice.
    threads : integer, optional (default 0)
      The number of parallel threads that should be used. If this value is less
      than one the maximum number of OpenMP threads will be used -- this value
      can be controled by the environment variable `OMP_NUM_THREADS` and is
      typically the number of logical cores if not explicitly set.

    Returns
    -------
    :py:class:`numpy.ndarray`
      The array of equivalent irreducible points for all Q.
  )pbdoc");

  cls.def("ir_moveinto",[](CLS &b, py::array_t<double> Q, int threads){
    // brille::Array<double> sp = brille::py2a(Q);
    brille::Array2<double> sp = brille::py2a2(Q);
    if (sp.size(sp.ndim()-1) != 3)
      throw std::runtime_error("The last dimension must have size 3");
    LQVec<double> Qv( b.get_lattice(),  sp); // view
    // prepare intermediate outputs
    LQVec<double> qv(b.get_lattice(), sp.shape(), sp.stride()); // output
    LQVec<int>  tauv(b.get_lattice(), sp.shape(), sp.stride()); // output
    std::vector<size_t> rotidx(Qv.numel()/3), invrotidx(Qv.numel()/3);
    if (!b.ir_moveinto(Qv, qv, tauv, rotidx, invrotidx, threads))
      throw std::runtime_error("Moving points into irreducible zone failed.");
    // get the pointgroup symmetry operations indexed by rotidx and invrotidx
    PointSymmetry ptsym = b.get_pointgroup_symmetry();
    // prepare Python outputs
    // The rotations array has an extra dimension compared to q and tau
    std::vector<ssize_t> sh;
    for (auto s: Qv.shape()) sh.push_back(static_cast<ssize_t>(s));
    sh.push_back(3);
    auto rout = py::array_t<int,    py::array::c_style>(sh);
    auto invrout = py::array_t<int, py::array::c_style>(sh);
    // grab pointers to the underlying data blocks
    int *rptr = (int *) rout.request().ptr;
    int *iptr = (int *) invrout.request().ptr;
    for (ssize_t i=0; i<Qv.numel()/3; ++i)
    for (size_t j=0; j<3u; ++j) for (size_t k=0; k<3u; ++k) {
      rptr[9u*i+3u*j+k] = ptsym.get(rotidx[i])[3u*j+k];
      iptr[9u*i+3u*j+k] = ptsym.get(invrotidx[i])[3u*j+k];
    }
    return py::make_tuple(brille::a2py(qv), brille::a2py(tauv), rout, invrout);
  }, "Q"_a, "threads"_a=0, R"pbdoc(
    Find points equivalent to those provided within the irreducible Brillouin zone.

    The BrillouinZone object defines a volume of reciprocal space which contains
    an irreducible part of the full reciprocal-space. This method will find
    points equivalent under the operations of the lattice which fall within this
    irreducible volume.

    Parameters
    ----------
    Q : :py:class:`numpy.ndarray`
      A 2+ dimensional array of three-vectors (`Q.shape[-1]==3`) expressed in
      units of the reciprocal lattice.
    threads : integer, optional (default 0)
      The number of parallel threads that should be used. If this value is less
      than one the maximum number of OpenMP threads will be used -- this value
      can be controled by the environment variable `OMP_NUM_THREADS` and is
      typically the number of logical cores if not explicitly set.

    Returns
    -------
    :py:class:`numpy.ndarray`
      The array of equivalent irreducible points for all Q.
  )pbdoc");

  cls.def("ir_moveinto_wedge",[](CLS &b, py::array_t<double> Q, int threads){
    // brille::Array<double> sp = brille::py2a(Q);
    brille::Array2<double> sp = brille::py2a2(Q);
    if (sp.size(sp.ndim()-1) != 3)
      throw std::runtime_error("The last dimension must have size 3");
    LQVec<double> Qv( b.get_lattice(),  sp); // view
    // prepare intermediate outputs
    LQVec<double> qv(b.get_lattice(), sp.shape(), sp.stride()); // output
    std::vector<std::array<int,9>> rots(Qv.numel()/3);
    if (!b.ir_moveinto_wedge(Qv, qv, rots, threads))
      throw std::runtime_error("Moving points into irreducible zone failed.");
    // prepare Python outputs
    // The rotations array has an extra dimension compared to q and tau
    std::vector<ssize_t> sh;
    for (auto s: Qv.shape()) sh.push_back(static_cast<ssize_t>(s));
    sh.push_back(3);
    auto rout = py::array_t<int,    py::array::c_style>(sh);
    // grab pointers to the underlying data blocks
    int *rptr = (int *) rout.request().ptr;
    for (ssize_t i=0; i<Qv.numel()/3; ++i)
    for (size_t j=0; j<3u; ++j) for (size_t k=0; k<3u; ++k)
      rptr[9u*i+3u*j+k] = rots[i][3u*j+k];
    return py::make_tuple(brille::a2py(qv), rout);
  }, "Q"_a, "threads"_a=0, R"pbdoc(
    Find points equivalent to those provided within the irreducible wedge.

    The BrillouinZone object defines a wedge of reciprocal space which contains
    an irreducible part of the full-space 4π steradian solid angle. This method
    will find points equivalent under the pointgroup operations of the lattice
    which fall within this irreducible solid angle and maintain their absolute
    magntiude.

    Parameters
    ----------
    Q : :py:class:`numpy.ndarray`
      A 2+ dimensional array of three-vectors (`Q.shape[-1]==3`) expressed in
      units of the reciprocal lattice.
    threads : integer, optional (default 0)
      The number of parallel threads that should be used. If this value is less
      than one the maximum number of OpenMP threads will be used -- this value
      can be controled by the environment variable `OMP_NUM_THREADS` and is
      typically the number of logical cores if not explicitly set.

    Returns
    -------
    :py:class:`numpy.ndarray`
      The array of equivalent in-wedge points for all Q
  )pbdoc");
}
