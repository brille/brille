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
#ifndef BRILLE_LATTICE_BOUNDARY_HPP_
#define BRILLE_LATTICE_BOUNDARY_HPP_
/*! \file
\brief The irreducible zone's boundary for the structured mesh, exactly

The irreducible polyhedron (first zone ∩ the point group's Dirichlet cone) with
named vertices, its faces, the maps that pair them, and the pairing cells: each
face F_j split into the convex pieces F_j ∩ g(F_i) of positive area. Cells pair
whole, so meshing each cell consistently makes paired faces match (see the design
note). All decisions use exact predicates (exact_geometry.hpp).

Coordinates are in the primitive reciprocal lattice basis; the metric G and the
point group operations (acting on those coordinates) are the inputs.
*/
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <tuple>
#include "exact_polytope.hpp"

namespace brille::latticetri {
using exact::Expansion;
using exact::Geometry;
using exact::Metric;
using exact::Plane;
using exact::Point;
using int3 = exact::int3;
using mat3i = exact::mat3i;

//! A map x -> g x + t sending one face's plane onto another's
struct FaceMap {
  mat3i g;
  int3 t;
  int from;
  int to;
};

using exact::Polygon;
using exact::Polytope;

class Boundary {
  Geometry geom_;
  std::vector<mat3i> ops_;
  std::vector<Plane> planes_;            // face planes, oriented inside n·x <= d
  std::vector<Polygon> faces_;           // one per face plane
  std::vector<std::vector<Polygon>> cells_;
  std::vector<std::vector<FaceMap>> cell_maps_;
  std::vector<FaceMap> maps_;
  std::vector<Point> special_;
  std::optional<std::vector<int3>> cone_;

public:
  /*! \param metric the reciprocal metric
      \param ops the point group acting on reciprocal coordinates
      \param cone the wedge, as integer c with c·x >= 0 inside; the Dirichlet cone of
             `ops` about (7, 3, 1) if not given
      \param reach the largest lattice translation, per coordinate, searched for face maps */
  Boundary(const std::array<double, 9> & metric, std::vector<mat3i> ops,
           std::optional<std::vector<int3>> cone = std::nullopt, const int reach = 3)
      : geom_(Metric(symmetrized(metric, ops))), ops_(std::move(ops)), cone_(std::move(cone)) {
    for (const auto & g: ops_)
      if (!geom_.metric().invariant_under(g))
        throw std::invalid_argument("the metric is not exactly invariant under the point group");
    polyhedron();
    face_maps(reach);
    pairing_cells();
    special_points(reach);
  }
  [[nodiscard]] const Geometry & geometry() const { return geom_; }
  [[nodiscard]] const std::vector<Plane> & planes() const { return planes_; }
  [[nodiscard]] const std::vector<Polygon> & faces() const { return faces_; }
  [[nodiscard]] const std::vector<std::vector<Polygon>> & cells() const { return cells_; }
  [[nodiscard]] const std::vector<FaceMap> & maps() const { return maps_; }
  [[nodiscard]] const std::vector<Point> & special_points() const { return special_; }

private:
  /*! \brief The metric nearest G that the point group leaves exactly invariant, in doubles

  The invariant metrics are spanned by integer matrices: the group sums of the six
  elementary symmetric matrices. G is written in an independent set of them, and
  its coordinates are rounded to a common grid 2⁻⁴⁰ below the largest, so every
  entry, a sum of small integer multiples of the coordinates, is exact. Rounding
  each entry of the group average instead leaves relations such as G₂₂ = G₀₂ + G₁₂
  (body-centred lattices in a primitive basis) true only to round-off.
  */
  static std::array<double, 9> symmetrized(const std::array<double, 9> & g, const std::vector<mat3i> & ops) {
    using six = std::array<double, 6>;
    constexpr std::array<std::array<int, 2>, 6> entries{{{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}}};
    // the group sum of the elementary symmetric matrix with ones at (k, l) and (l, k)
    auto group_sum = [&](const int k, const int l) {
      six out{};
      for (const auto & w: ops)
        for (size_t e = 0; e < 6; ++e) {
          const auto [i, j] = entries[e];
          // (wᵀ S w)_ij = w_ki w_lj + w_li w_kj, or w_ki w_kj when k == l
          long long v = w[3 * k + i] * w[3 * l + j];
          if (k != l) v += w[3 * l + i] * w[3 * k + j];
          out[e] += static_cast<double>(v);
        }
      return out;
    };
    // an independent set, by Gram-Schmidt on copies
    std::vector<six> basis, orthogonal;
    for (const auto & [k, l]: entries) {
      const auto b = group_sum(k, l);
      auto r = b;
      for (const auto & q: orthogonal) {
        double d{0}, n{0};
        for (size_t e = 0; e < 6; ++e) { d += r[e] * q[e]; n += q[e] * q[e]; }
        for (size_t e = 0; e < 6; ++e) r[e] -= d / n * q[e];
      }
      double rn{0}, bn{0};
      for (size_t e = 0; e < 6; ++e) { rn += r[e] * r[e]; bn += b[e] * b[e]; }
      if (bn > 0 && rn > 1e-18 * bn) { basis.push_back(b); orthogonal.push_back(r); }
    }
    // G's coordinates in that basis, by least squares: (EᵀE) t = Eᵀ g
    const size_t n = basis.size();
    six gv{};
    for (size_t e = 0; e < 6; ++e) gv[e] = g[3 * entries[e][0] + entries[e][1]];
    std::vector<std::vector<double>> A(n, std::vector<double>(n + 1, 0.0));
    for (size_t r = 0; r < n; ++r) {
      for (size_t c = 0; c < n; ++c) for (size_t e = 0; e < 6; ++e) A[r][c] += basis[r][e] * basis[c][e];
      for (size_t e = 0; e < 6; ++e) A[r][n] += basis[r][e] * gv[e];
    }
    for (size_t c = 0; c < n; ++c) {
      size_t pivot = c;
      for (size_t r = c + 1; r < n; ++r) if (std::abs(A[r][c]) > std::abs(A[pivot][c])) pivot = r;
      std::swap(A[c], A[pivot]);
      for (size_t r = 0; r < n; ++r) {
        if (r == c) continue;
        const double f = A[r][c] / A[c][c];
        for (size_t k = c; k <= n; ++k) A[r][k] -= f * A[c][k];
      }
    }
    std::vector<double> t(n);
    double largest{0}, coefficient{0};
    for (size_t r = 0; r < n; ++r) { t[r] = A[r][n] / A[r][r]; largest = std::max(largest, std::abs(t[r])); }
    for (const auto & b: basis) for (const auto x: b) coefficient = std::max(coefficient, std::abs(x));
    if (largest == 0) throw std::invalid_argument("the metric has no invariant part");
    // a grid fine enough to move G by ~1e-12, coarse enough that six terms with
    // integer factors up to `coefficient` sum exactly
    const double quantum = std::ldexp(1.0, std::ilogb(largest * coefficient) - 40);
    for (auto & x: t) x = std::round(x / quantum) * quantum;
    std::array<double, 9> out{};
    for (size_t e = 0; e < 6; ++e) {
      double v{0};
      for (size_t r = 0; r < n; ++r) v += t[r] * basis[r][e];
      const auto [i, j] = entries[e];
      out[3 * i + j] = out[3 * j + i] = v;
    }
    return out;
  }

  /*! The irreducible polyhedron: a box cut by the zone planes and the cone planes,
  with each vertex named by three planes and its exact set of incident planes. */
  void polyhedron() {
    std::vector<Plane> box;
    for (int i = 0; i < 3; ++i) {
      int3 e{0, 0, 0};
      e[i] = 1;
      box.push_back(Plane::integer_plane(e, 3));
      e[i] = -1;
      box.push_back(Plane::integer_plane(e, 3));
    }
    std::vector<Polytope::Vertex> corners;
    for (int a: {0, 1}) for (int b: {2, 3}) for (int c: {4, 5}) corners.push_back({{{box[a], box[b], box[c]}}, {a, b, c}});
    Polytope P(geom_, box, corners);
    // zone planes x·Gτ <= τᵀGτ/2, shortest first
    std::vector<int3> taus;
    for (int a = -2; a <= 2; ++a) for (int b = -2; b <= 2; ++b) for (int c = -2; c <= 2; ++c) if (a || b || c) taus.push_back({a, b, c});
    std::sort(taus.begin(), taus.end(), [&](const int3 & x, const int3 & y) {
      return geom_.metric().form(x, x).estimate() < geom_.metric().form(y, y).estimate();
    });
    for (const auto & t: taus) P.cut(Plane::metric_plane(t, t));
    // the Dirichlet cone of the point group under M = Σ gᵀg: x·M(p - g p) >= 0
    std::vector<int3> normals = cone_ ? *cone_ : std::vector<int3>{};
    if (!cone_) {
      std::array<long long, 9> M{};
      for (const auto & g: ops_) for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) M[3 * i + j] += g[3 * k + i] * g[3 * k + j];
      const int3 p{7, 3, 1};
      for (const auto & g: ops_) {
        const auto gp = Plane::apply(g, p);
        if (gp == p) continue;
        const int3 d{p[0] - gp[0], p[1] - gp[1], p[2] - gp[2]};
        int3 c{M[0] * d[0] + M[1] * d[1] + M[2] * d[2], M[3] * d[0] + M[4] * d[1] + M[5] * d[2], M[6] * d[0] + M[7] * d[1] + M[8] * d[2]};
        const long long h = std::gcd(std::gcd(std::llabs(c[0]), std::llabs(c[1])), std::llabs(c[2]));
        for (auto & x: c) x /= h;
        if (std::find(normals.begin(), normals.end(), c) == normals.end()) normals.push_back(c);
      }
    }
    for (const auto & c: normals) P.cut(Plane::integer_plane({-c[0], -c[1], -c[2]}, 0));
    for (const auto & v: P.vertices())
      for (int i = 0; i < 6; ++i)
        if (v.incident.count(i)) throw std::runtime_error("the box bounds the zone: enlarge it");
    for (const int f: P.faces(6)) {
      planes_.push_back(P.planes()[static_cast<size_t>(f)]);
      faces_.push_back(P.polygon(f));
    }
  }

  //! whether a and b are the same plane with the same inside
  [[nodiscard]] bool parallel_names(const Plane & a, const Plane & b) const {
    if (a.kind == b.kind) {
      const auto & x = a.a;
      const auto & y = b.a;
      return x[0] * y[1] - x[1] * y[0] == 0 && x[1] * y[2] - x[2] * y[1] == 0 && x[0] * y[2] - x[2] * y[0] == 0;
    }
    return true;   // a metric and an integer plane: decide exactly
  }

  /*! Maps x -> g x + t, for every operation and t within `reach`, sending a face's
  plane onto a face plane (other than the identity on a face) */
  void face_maps(const int reach) {
    for (size_t i = 0; i < planes_.size(); ++i)
      for (const auto & g: ops_)
        for (size_t j = 0; j < planes_.size(); ++j) {
          if (!parallel_names(planes_[i].mapped(g, {0, 0, 0}), planes_[j])) continue;
          for (long long a = -reach; a <= reach; ++a)
            for (long long b = -reach; b <= reach; ++b)
              for (long long c = -reach; c <= reach; ++c) {
                const int3 t{a, b, c};
                if (i == j && t == int3{0, 0, 0} && g == identity()) continue;
                if (geom_.same(planes_[i].mapped(g, t), planes_[j])) maps_.push_back({g, t, static_cast<int>(i), static_cast<int>(j)});
              }
        }
  }
  static mat3i identity() { return {1, 0, 0, 0, 1, 0, 0, 0, 1}; }

  /*! Each face F_j split into the convex cells F_j ∩ g(F_i) of positive area */
  void pairing_cells() {
    cells_.assign(planes_.size(), {});
    cell_maps_.assign(planes_.size(), {});
    for (const auto & m: maps_) {
      const auto & src = faces_[static_cast<size_t>(m.from)];
      Polygon cell = faces_[static_cast<size_t>(m.to)];
      for (const auto & e: src.edges) {
        cell = exact::clip(geom_, cell, planes_[static_cast<size_t>(m.to)], e.mapped(m.g, m.t));
        if (cell.vertices.empty()) break;
      }
      if (cell.vertices.empty()) continue;
      cells_[static_cast<size_t>(m.to)].push_back(cell);
      cell_maps_[static_cast<size_t>(m.to)].push_back(m);
    }
  }

  //! Whether the point lies strictly inside an edge of one of the cells
  [[nodiscard]] bool on_cell_edge(const Point & x) const {
    for (size_t j = 0; j < cells_.size(); ++j) {
      if (!geom_.on(x, planes_[j])) continue;
      for (const auto & cell: cells_[j])
        for (size_t k = 0; k < cell.vertices.size(); ++k) {
          if (!geom_.on(x, cell.edges[k])) continue;
          // inside the cell: no other edge has x outside
          bool inside{true};
          for (size_t l = 0; l < cell.edges.size() && inside; ++l) if (l != k) inside = geom_.side(x, cell.edges[l]) <= 0;
          if (inside) return true;
        }
    }
    return false;
  }

  /*! Cell corners, and every image of one (operation and translation) that lies on
  a cell edge: paired cells then share their vertices along edges. */
  void special_points(const int reach) {
    auto known = [&](const Point & p) { return std::any_of(special_.begin(), special_.end(), [&](const Point & q) { return geom_.same(p, q); }); };
    for (const auto & cs: cells_) for (const auto & c: cs) for (const auto & p: c.vertices) if (!known(p)) special_.push_back(p);
    std::vector<Point> todo = special_;
    // approximate face data, to rule out most images before exact tests
    while (!todo.empty()) {
      const Point x = todo.back();
      todo.pop_back();
      const auto xc = geom_.coordinates(x);
      for (const auto & g: ops_) {
        std::array<double, 3> gx{0, 0, 0};
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) gx[i] += static_cast<double>(g[3 * i + k]) * xc[k];
        for (long long a = -reach; a <= reach; ++a)
          for (long long b = -reach; b <= reach; ++b)
            for (long long c = -reach; c <= reach; ++c) {
              const std::array<double, 3> y{gx[0] + a, gx[1] + b, gx[2] + c};
              if (!near_boundary(y)) continue;
              const Point q = geom_.mapped(x, g, {a, b, c});
              if (known(q) || !inside(q) || !on_cell_edge(q)) continue;
              special_.push_back(q);
              todo.push_back(q);
            }
      }
    }
  }
  [[nodiscard]] bool inside(const Point & q) const {
    return std::all_of(planes_.begin(), planes_.end(), [&](const Plane & p) { return geom_.side(q, p) <= 0; });
  }
  [[nodiscard]] bool near_boundary(const std::array<double, 3> & y) const {
    bool on_any{false};
    for (const auto & p: planes_) {
      const auto c = p.coefficients(geom_.metric());
      const double v = c[0].estimate() * y[0] + c[1].estimate() * y[1] + c[2].estimate() * y[2] - c[3].estimate();
      const double scale = std::abs(c[0].estimate()) + std::abs(c[1].estimate()) + std::abs(c[2].estimate()) + std::abs(c[3].estimate());
      if (v > 1e-9 * scale) return false;
      on_any |= std::abs(v) <= 1e-9 * scale;
    }
    return on_any;
  }
};
}
#endif
