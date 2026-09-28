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

#include <numeric>
#if defined(__GLIBC__)
#include <malloc.h>
#endif
#include <optional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "_array.hpp"
#include "_common_grid.hpp"
#include "bz_mesh.hpp"

#ifndef WRAP_BRILLE_MESH_HPP_
#define WRAP_BRILLE_MESH_HPP_

namespace py = pybind11;

template<class T,class R,class S>
void declare_bzmeshq(py::module &m, const std::string &typestr){
  using namespace pybind11::literals;
  using namespace brille;
  using Class = BrillouinZoneMesh3<T,R,S>;
  std::string pyclass_name = std::string("BZMeshQ")+typestr;
  py::class_<Class> cls(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr());
  // Initializer (BrillouinZone, max-volume, is-volume-rlu)
  cls.def(py::init([](const BrillouinZone& bz, const double max_size, const int num_levels, const int max_points){
    auto mesh = [&]{
      py::gil_scoped_release release;
      return Class(bz, max_size, num_levels, max_points);
    }();
    if (mesh.refinement_limited()) {
      const std::string msg = "max_points made the mesh grid coarser than max_size asked for; its tetrahedra"
        " are larger than requested.";
      if (PyErr_WarnEx(PyExc_RuntimeWarning, msg.c_str(), 1) < 0) throw py::error_already_set();
    }
    return mesh;
  }), "brillouin_zone"_a, "max_size"_a=-1., "num_levels"_a=3, "max_points"_a=-1);
  cls.def_property_readonly("refinement_limited", [](const Class& cobj){return cobj.refinement_limited();},
    "Whether max_points made the mesh coarser than max_size asked for");
  cls.def_property_readonly("BrillouinZone",[](const Class& cobj){return cobj.get_brillouinzone();});
  cls.def_property_readonly("rlu",[](const Class& cobj){
    return brille::a2py(cobj.get_mesh_hkl());
  });
  cls.def_property_readonly("invA",[](const Class& cobj){
    return brille::a2py(cobj.get_mesh_xyz());
  });
  cls.def_property_readonly("tetrahedra",[](const Class& cobj){
    return brille::a2py(cobj.get_mesh_tetrehedra());
  });
  // `where` as tetrahedron indices: None for all, a boolean mask, or indices
  auto tetrahedra_of = [](const Class& cobj, const py::object& where){
    const auto count = cobj.get_mesh_tetrehedra().size(0);
    std::vector<ind_t> out;
    if (where.is_none()) {
      out.resize(count);
      std::iota(out.begin(), out.end(), ind_t(0));
      return out;
    }
    auto numpy = py::module_::import("numpy");
    py::array a = numpy.attr("asarray")(where);
    if (a.ndim() != 1) throw py::value_error("where must be a 1-D boolean mask or array of tetrahedron indices");
    if (py::str(a.dtype().attr("kind")).cast<std::string>() == "b") {
      if (static_cast<ind_t>(a.shape(0)) != count)
        throw py::value_error("a boolean where must have one entry per tetrahedron (" + std::to_string(count) + ")");
      const py::array_t<bool> mask = py::cast<py::array_t<bool>>(a);
      auto m = mask.unchecked<1>();
      for (py::ssize_t i = 0; i < m.shape(0); ++i) if (m(i)) out.push_back(static_cast<ind_t>(i));
      return out;
    }
    // a declared type, not auto: GCC 13 takes auto here as dependent and wants `template` before unchecked
    const py::array_t<long long> indices = py::cast<py::array_t<long long>>(numpy.attr("asarray")(a, "dtype"_a = "int64"));
    auto idx = indices.unchecked<1>();
    for (py::ssize_t i = 0; i < idx.shape(0); ++i) {
      if (idx(i) < 0 || static_cast<ind_t>(idx(i)) >= count)
        throw py::index_error("tetrahedron index " + std::to_string(idx(i)) + " is not in the mesh (" + std::to_string(count) + ")");
      out.push_back(static_cast<ind_t>(idx(i)));
    }
    return out;
  };
  auto edge_limit = [](const std::optional<double>& resolution, const double per){
    if (!resolution) return 0.0;
    if (*resolution <= 0 || per <= 0) throw py::value_error("resolution and points_per_resolution must be positive");
    return *resolution / per;
  };
  cls.def_property_readonly("refinable", [](const Class& cobj){return cobj.refinable();},
    "Whether the mesh can be refined; a mesh read from a file written before refinement existed can't be");
  cls.def_property_readonly("holds_triangulation", [](const Class& cobj){return cobj.holds_triangulation();},
    "Whether the triangulation that refinement works on is in memory (see release_triangulation)");
  cls.def("release_triangulation", [](Class& cobj){
    cobj.release_triangulation();
#if defined(__GLIBC__)
    malloc_trim(0);
#endif
  },
R"pbdoc(
Free the memory refinement holds between refinements.

After :py:meth:`refine` (or :py:meth:`refinement_points`) the mesh keeps the
triangulation refinement works on, several times the memory of the mesh itself.
This frees it. The mesh is unchanged and can still be refined: the triangulation
is rebuilt when next needed, which costs a build of the mesh plus a replay of the
refinements made so far.
)pbdoc");
  cls.def("refinement_points", [tetrahedra_of, edge_limit](const Class& cobj, const py::object& where,
                                                           const std::optional<double>& resolution, const double per){
    const auto tets = tetrahedra_of(cobj, where);
    const auto min_edge = edge_limit(resolution, per);
    auto points = [&]{ py::gil_scoped_release release; return cobj.refinement_points_hkl(tets, min_edge); }();
    return brille::a2py(points);
  }, "where"_a=py::none(), "resolution"_a=py::none(), "points_per_resolution"_a=2.0,
R"pbdoc(
The points that :py:meth:`refine` would add, without changing the mesh.

Evaluate your model at these points and compare with
:py:meth:`ir_interpolate_at` there to decide whether refining is worth it; then
pass the model's values to :py:meth:`refine` with the same arguments.

Parameters
----------
where : None, bool array or int array
  The tetrahedra to split: all of them (None), those where a boolean mask with
  one entry per tetrahedron is true, or those with the given indices.
  Neighbouring tetrahedra are split as needed to keep the mesh conforming, and
  split edges on the zone boundary are split with their symmetry equivalents, so
  that equivalent zone faces keep matching.
resolution : float, optional
  The resolution limit, in inverse Angstrom. No edge is split to below
  ``resolution / points_per_resolution``: a tetrahedron whose longest edge is at
  most twice that is left whole.
points_per_resolution : float, optional (default: 2)
  How finely to resolve ``resolution``.

Returns
-------
numpy.ndarray
  The new vertices, shape (N, 3), in relative lattice units like :py:attr:`rlu`.
  :py:meth:`refine` appends them to the vertices in this order.
)pbdoc");
  cls.def("refine", [tetrahedra_of, edge_limit](Class& cobj, const py::object& where,
                                                const py::object& values, const py::object& vectors,
                                                const std::optional<double>& resolution, const double per){
    const auto tets = tetrahedra_of(cobj, where);
    const auto min_edge = edge_limit(resolution, per);
    std::optional<brille::Array<T>> vals;
    std::optional<brille::Array<R>> vecs;
    py::array_t<T> pyvals;
    py::array_t<R> pyvecs;
    if (!values.is_none()) { pyvals = values.cast<py::array_t<T>>(); vals = brille::py2a(pyvals); }
    if (!vectors.is_none()) { pyvecs = vectors.cast<py::array_t<R>>(); vecs = brille::py2a(pyvecs); }
    // only planning runs without the GIL: applying replaces the held data, which can
    // release a numpy buffer
    auto plan = [&]{ py::gil_scoped_release release; return cobj.plan_refinement(tets, min_edge); }();
    return brille::a2py(cobj.apply_refinement_hkl(std::move(plan), vals, vecs));
  }, "where"_a=py::none(), "values"_a=py::none(), "vectors"_a=py::none(),
     "resolution"_a=py::none(), "points_per_resolution"_a=2.0,
R"pbdoc(
Refine the mesh by bisecting tetrahedra; existing vertices keep their indices and data.

Parameters
----------
where, resolution, points_per_resolution
  As for :py:meth:`refinement_points`, which gives the points this adds.
values, vectors : numpy.ndarray, optional
  If the mesh holds data (after :py:meth:`fill`), the data for the new vertices,
  laid out per point as the filled data and for exactly the points
  :py:meth:`refinement_points` returns, in that order. Not allowed before the
  mesh is filled.

Returns
-------
numpy.ndarray
  The new vertices, shape (N, 3), in relative lattice units, appended to
  :py:attr:`rlu` in this order.

Note
----
The mode permutations found by :py:meth:`sort` are reset; sort again after
refining if needed.
)pbdoc");
  cls.def("__repr__",&Class::to_string);

  def_grid_fill(cls);
  def_grid_ir_interpolate(cls);
  def_grid_sort(cls);
//  def_grid_debye_waller(cls);

  def_grid_hdf_interface(cls, pyclass_name);
}

#endif
