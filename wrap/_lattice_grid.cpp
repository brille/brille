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
#include "lattice_boundary.hpp"
#include "lattice_tri.hpp"
#include "lattice_mesh.hpp"

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

namespace {
py::array_t<double> coordinates(const brille::exact::Geometry & geom, const std::vector<brille::exact::Point> & pts) {
  py::array_t<double> out({static_cast<py::ssize_t>(pts.size()), py::ssize_t(3)});
  auto r = out.mutable_unchecked<2>();
  for (size_t i = 0; i < pts.size(); ++i) {
    const auto x = geom.coordinates(pts[i]);
    for (int k = 0; k < 3; ++k) r(i, k) = x[k];
  }
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
  py::class_<Boundary> bnd(m, "_LatticeBoundary", R"pbdoc(
    The irreducible zone's boundary, exactly: faces, pairing cells, special points.

    Internal: part of the structured mesh under development; for tests only.
    Coordinates are in the primitive reciprocal lattice basis.
  )pbdoc");
  bnd.def(py::init([](const py::array_t<double> & metric, const std::vector<py::array_t<long long>> & operations,
                      const std::optional<std::vector<std::array<long long, 3>>> & cone) {
    auto g = metric.unchecked<2>();
    std::array<double, 9> G{};
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) G[3 * i + k] = g(i, k);
    std::vector<mat3i> ops;
    for (const auto & o: operations) {
      auto a = o.unchecked<2>();
      mat3i r{};
      for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) r[3 * i + k] = a(i, k);
      ops.push_back(r);
    }
    py::gil_scoped_release release;
    return Boundary(G, ops, cone);
  }), "metric"_a, "operations"_a, "cone"_a = py::none());
  bnd.def_property_readonly("faces", [](const Boundary & b) {
    py::list out;
    for (const auto & f: b.faces()) out.append(coordinates(b.geometry(), f.vertices));
    return out;
  });
  bnd.def_property_readonly("cells", [](const Boundary & b) {
    py::list out;
    for (const auto & cs: b.cells()) {
      py::list face;
      for (const auto & c: cs) face.append(coordinates(b.geometry(), c.vertices));
      out.append(face);
    }
    return out;
  });
  bnd.def_property_readonly("special_points", [](const Boundary & b) { return coordinates(b.geometry(), b.special_points()); });
  bnd.def_property_readonly("map_count", [](const Boundary & b) { return b.maps().size(); });

  py::class_<LatticeTri> tri(m, "_LatticeTri", R"pbdoc(
    The structured mesh of the irreducible zone: the grid clipped to the zone.

    Internal: under development; for tests only. Vertices are in the primitive
    reciprocal lattice basis; the grid lattice is that lattice divided by `n`.
  )pbdoc");
  tri.def(py::init([](const py::array_t<double> & metric, const std::vector<py::array_t<long long>> & operations, long long n,
                      const std::optional<std::vector<std::array<long long, 3>>> & cone) {
    auto g = metric.unchecked<2>();
    std::array<double, 9> G{};
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) G[3 * i + k] = g(i, k);
    std::vector<mat3i> ops;
    for (const auto & o: operations) {
      auto a = o.unchecked<2>();
      mat3i r{};
      for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) r[3 * i + k] = a(i, k);
      ops.push_back(r);
    }
    py::gil_scoped_release release;
    return LatticeTri(G, ops, n, cone);
  }), "metric"_a, "operations"_a, "n"_a, "cone"_a = py::none());
  tri.def_property_readonly("vertices", [](const LatticeTri & t) {
    py::array_t<double> out({static_cast<py::ssize_t>(t.vertices().size()), py::ssize_t(3)});
    auto r = out.mutable_unchecked<2>();
    for (size_t i = 0; i < t.vertices().size(); ++i) for (int k = 0; k < 3; ++k) r(i, k) = t.vertices()[i][k];
    return out;
  });
  tri.def_property_readonly("tetrahedra", [](const LatticeTri & t) {
    py::array_t<long long> out({static_cast<py::ssize_t>(t.tetrahedra().size()), py::ssize_t(4)});
    auto r = out.mutable_unchecked<2>();
    for (size_t i = 0; i < t.tetrahedra().size(); ++i) for (int k = 0; k < 4; ++k) r(i, k) = static_cast<long long>(t.tetrahedra()[i][k]);
    return out;
  });
  tri.def_property_readonly("clipped", &LatticeTri::clipped);
  tri.def_property_readonly("self_paired_ties", &LatticeTri::self_paired_ties);
  tri.def("planes_of", [](const LatticeTri & t, const size_t v) { const auto & s = t.planes_of(v); return std::vector<int>(s.begin(), s.end()); });
  tri.def("refine", [](LatticeTri & t, const py::array_t<long long> & marked, const double min_edge) {
    auto r = marked.unchecked<2>();
    std::vector<std::array<size_t, 4>> tets;
    for (py::ssize_t i = 0; i < r.shape(0); ++i) tets.push_back({static_cast<size_t>(r(i, 0)), static_cast<size_t>(r(i, 1)), static_cast<size_t>(r(i, 2)), static_cast<size_t>(r(i, 3))});
    py::gil_scoped_release release;
    t.refine(tets, min_edge);
  }, "marked"_a, "min_edge"_a = 0.0, "Refine the marked tetrahedra (rows of vertex indices) once, with closure");
  tri.def_property_readonly("faces", [](const LatticeTri & t) {
    py::list out;
    for (const auto & f: t.boundary().faces()) out.append(coordinates(t.geometry(), f.vertices));
    return out;
  });

  cls.def("locate", [](const Grid & g, const std::array<double, 3> & x) {
    tetrahedron t{};
    std::array<double, 4> w{};
    if (!g.locate(x, t, w)) throw std::runtime_error("point not located");
    return py::make_tuple(tetrahedra_array({t}), w);
  }, "x"_a);

  m.def("_lattice_mesh_inputs", [](const brille::BrillouinZone & bz) {
    const auto in = brille::LatticeMesh::inputs(bz);
    py::array_t<double> metric({3, 3}), basis({3, 3});
    auto g = metric.mutable_unchecked<2>();
    auto b = basis.mutable_unchecked<2>();
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) { g(i, k) = in.metric[3 * i + k]; b(i, k) = in.basis[3 * i + k]; }
    py::list ops;
    for (const auto & r: in.ops) {
      py::array_t<long long> o({3, 3});
      auto a = o.mutable_unchecked<2>();
      for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) a(i, k) = r[3 * i + k];
      ops.append(o);
    }
    return py::make_tuple(metric, ops, in.cone, basis);
  }, "brillouin_zone"_a, R"pbdoc(
    (metric, operations, cone, basis) for meshing a zone's irreducible part: the
    primitive reciprocal metric, the point group on primitive reciprocal
    coordinates, the zone's wedge as integer normals (c·x >= 0 inside), and the
    primitive reciprocal vectors as columns. Internal; for tests only.
  )pbdoc");
}
