/* This file is part of brille.

Copyright © 2026 Greg Tucker <gregory.tucker@ess.eu>

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
#include <pybind11/stl.h>
#include "lattice_grid.hpp"

namespace py = pybind11;

namespace {
py::array_t<long long> tetrahedra_array(const std::vector<brille::latticetri::tetrahedron> & tets) {
  py::array_t<long long> out({static_cast<py::ssize_t>(tets.size()), py::ssize_t(4), py::ssize_t(3)});
  auto r = out.mutable_unchecked<3>();
  for (size_t i = 0; i < tets.size(); ++i)
    for (int k = 0; k < 4; ++k)
      for (int j = 0; j < 3; ++j) r(i, k, j) = tets[i][k][j];
  return out;
}
}

/* The grid of the structured mesh that is to replace TetGen, exposed privately so
 * that tests can compare it with the Python reference implementation. Not API. */
void wrap_lattice_grid(py::module & m) {
  using namespace brille::latticetri;
  using namespace pybind11::literals;
  py::class_<Grid> cls(m, "_LatticeGrid", R"pbdoc(
    A periodic triangulation of a lattice, invariant under its point group.

    Internal: part of the structured mesh under development; for tests only.
    Grid points are integer coordinates in the lattice basis, times `scale`.
  )pbdoc");
  cls.def(py::init([](const py::array_t<double> & basis_rows, const std::vector<py::array_t<long long>> & operations) {
    auto b = basis_rows.unchecked<2>();
    std::array<double, 9> rows{};
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) rows[3 * i + k] = b(i, k);
    std::vector<mat3i> ops;
    for (const auto & o: operations) {
      auto a = o.unchecked<2>();
      mat3i r{};
      for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) r[3 * i + k] = a(i, k);
      ops.push_back(r);
    }
    return Grid(rows, ops);
  }), "basis_rows"_a, "operations"_a);
  cls.def_property_readonly_static("scale", [](const py::object &) { return scale; });
  cls.def_property_readonly("degenerate", &Grid::degenerate);
  cls.def_property_readonly("pattern", [](const Grid & g) { return tetrahedra_array(g.pattern()); });
  cls.def("patch", [](const Grid & g, double radius) { return tetrahedra_array(g.patch(radius)); }, "radius"_a);
  cls.def("invariant", &Grid::invariant, "radius"_a);
  cls.def("locate", [](const Grid & g, const std::array<double, 3> & x) {
    tetrahedron t{};
    std::array<double, 4> w{};
    if (!g.locate(x, t, w)) throw std::runtime_error("point not located");
    return py::make_tuple(tetrahedra_array({t}), w);
  }, "x"_a);
}
