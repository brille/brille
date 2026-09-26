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
#include <set>
#include <tuple>
#include "exact_geometry.hpp"

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

//! A convex polygon in a face plane: its vertices, and the plane of each edge (vertex k to k+1)
struct Polygon {
  std::vector<Point> vertices;
  std::vector<Plane> edges;
};

class Boundary {
  Geometry geom_;
  std::vector<mat3i> ops_;
  std::vector<Plane> planes_;            // face planes, oriented inside n·x <= d
  std::vector<Polygon> faces_;           // one per face plane
  std::vector<std::vector<Polygon>> cells_;
  std::vector<std::vector<FaceMap>> cell_maps_;
  std::vector<FaceMap> maps_;
  std::vector<Point> special_;

public:
  Boundary(const std::array<double, 9> & metric, std::vector<mat3i> ops, const int reach = 3)
      : geom_(Metric(symmetrized(metric, ops))), ops_(std::move(ops)) {
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
  /*! (1/N) Σ gᵀ G g, each entry a fixed integer combination of G's upper triangle
  summed in a fixed order, so that entries equal by symmetry are identical (as in
  brille's Lattice). */
  static std::array<double, 9> symmetrized(const std::array<double, 9> & g, const std::vector<mat3i> & ops) {
    std::array<double, 9> out{};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        std::array<long long, 9> c{};
        for (const auto & w: ops)
          for (int k = 0; k < 3; ++k)
            for (int l = 0; l < 3; ++l) c[3 * std::min(k, l) + std::max(k, l)] += w[3 * k + i] * w[3 * l + j];
        double sum{0};
        for (int k = 0; k < 3; ++k) for (int l = k; l < 3; ++l) if (c[3 * k + l]) sum += static_cast<double>(c[3 * k + l]) * g[3 * k + l];
        out[3 * i + j] = sum / static_cast<double>(ops.size());
      }
    return out;
  }

  struct Vertex { Point point; std::set<int> incident; };

  /*! The irreducible polyhedron: a box cut by the zone planes and the cone planes,
  with each vertex named by three planes and its exact set of incident planes. */
  void polyhedron() {
    std::vector<Plane> all;
    for (int i = 0; i < 3; ++i) {
      int3 e{0, 0, 0};
      e[i] = 1;
      all.push_back(Plane::integer_plane(e, 3));
      e[i] = -1;
      all.push_back(Plane::integer_plane(e, 3));
    }
    std::vector<Vertex> vs;
    for (int a: {0, 1}) for (int b: {2, 3}) for (int c: {4, 5}) vs.push_back({{{all[a], all[b], all[c]}}, {a, b, c}});
    // zone planes x·Gτ <= τᵀGτ/2, shortest first
    std::vector<int3> taus;
    for (int a = -2; a <= 2; ++a) for (int b = -2; b <= 2; ++b) for (int c = -2; c <= 2; ++c) if (a || b || c) taus.push_back({a, b, c});
    std::sort(taus.begin(), taus.end(), [&](const int3 & x, const int3 & y) {
      return geom_.metric().form(x, x).estimate() < geom_.metric().form(y, y).estimate();
    });
    std::vector<Plane> cuts;
    for (const auto & t: taus) cuts.push_back(Plane::metric_plane(t, t));
    // the Dirichlet cone of the point group under M = Σ gᵀg: x·M(p - g p) >= 0
    std::array<long long, 9> M{};
    for (const auto & g: ops_) for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) M[3 * i + j] += g[3 * k + i] * g[3 * k + j];
    const int3 p{7, 3, 1};
    std::vector<int3> normals;
    for (const auto & g: ops_) {
      const auto gp = Plane::apply(g, p);
      if (gp == p) continue;
      const int3 d{p[0] - gp[0], p[1] - gp[1], p[2] - gp[2]};
      int3 c{M[0] * d[0] + M[1] * d[1] + M[2] * d[2], M[3] * d[0] + M[4] * d[1] + M[5] * d[2], M[6] * d[0] + M[7] * d[1] + M[8] * d[2]};
      const long long h = std::gcd(std::gcd(std::llabs(c[0]), std::llabs(c[1])), std::llabs(c[2]));
      for (auto & x: c) x /= h;
      if (std::find(normals.begin(), normals.end(), c) == normals.end()) normals.push_back(c);
    }
    for (const auto & c: normals) cuts.push_back(Plane::integer_plane({-c[0], -c[1], -c[2]}, 0));
    for (const auto & h: cuts) cut(all, vs, h);
    // faces: planes (not the box) with at least three vertices
    for (size_t i = 6; i < all.size(); ++i) {
      std::vector<Point> fv;
      for (const auto & v: vs) if (v.incident.count(static_cast<int>(i))) fv.push_back(v.point);
      if (fv.size() >= 3) {
        planes_.push_back(all[i]);
        faces_.push_back(order(fv, all[i], all, vs, static_cast<int>(i)));
      }
    }
    for (size_t i = 0; i < 6; ++i)
      for (const auto & v: vs)
        if (v.incident.count(static_cast<int>(i))) throw std::runtime_error("the box bounds the zone: enlarge it");
  }

  void cut(std::vector<Plane> & all, std::vector<Vertex> & vs, const Plane & h) {
    std::vector<int> side(vs.size());
    bool any_out{false};
    for (size_t i = 0; i < vs.size(); ++i) { side[i] = geom_.side(vs[i].point, h); any_out |= side[i] > 0; }
    const int k = static_cast<int>(all.size());
    all.push_back(h);
    if (!any_out) {
      for (size_t i = 0; i < vs.size(); ++i) if (side[i] == 0) vs[i].incident.insert(k);
      return;
    }
    std::vector<Vertex> out;
    for (size_t u = 0; u < vs.size(); ++u) {
      if (side[u] <= 0) continue;
      for (size_t v = 0; v < vs.size(); ++v) {
        if (side[v] >= 0 || !adjacent(vs, u, v)) continue;
        std::vector<int> common;
        std::set_intersection(vs[u].incident.begin(), vs[u].incident.end(), vs[v].incident.begin(), vs[v].incident.end(), std::back_inserter(common));
        // two of the common planes that meet h in a point
        bool made{false};
        for (size_t a = 0; a < common.size() && !made; ++a)
          for (size_t b = a + 1; b < common.size() && !made; ++b) {
            Point x{{all[common[a]], all[common[b]], h}};
            if (geom_.independent(x) == 0) continue;
            std::set<int> inc(common.begin(), common.end());
            inc.insert(k);
            out.push_back({x, inc});
            made = true;
          }
      }
    }
    std::vector<Vertex> kept;
    for (size_t i = 0; i < vs.size(); ++i) {
      if (side[i] > 0) continue;
      if (side[i] == 0) vs[i].incident.insert(k);
      kept.push_back(vs[i]);
    }
    for (auto & v: out) kept.push_back(v);
    vs.swap(kept);
  }

  static bool adjacent(const std::vector<Vertex> & vs, const size_t u, const size_t v) {
    std::vector<int> common;
    std::set_intersection(vs[u].incident.begin(), vs[u].incident.end(), vs[v].incident.begin(), vs[v].incident.end(), std::back_inserter(common));
    if (common.size() < 2) return false;
    for (size_t w = 0; w < vs.size(); ++w) {
      if (w == u || w == v) continue;
      if (std::includes(vs[w].incident.begin(), vs[w].incident.end(), common.begin(), common.end())) return false;
    }
    return true;
  }

  //! The face's vertices in cyclic order, with the plane of each edge
  Polygon order(const std::vector<Point> & fv, const Plane & plane, const std::vector<Plane> & all,
                const std::vector<Vertex> & vs, const int face) const {
    const auto c = plane.coefficients(geom_.metric());
    std::array<double, 3> n{c[0].estimate(), c[1].estimate(), c[2].estimate()};
    std::vector<std::array<double, 3>> x;
    for (const auto & p: fv) x.push_back(geom_.coordinates(p));
    std::array<double, 3> m{0, 0, 0};
    for (const auto & y: x) for (int i = 0; i < 3; ++i) m[i] += y[i] / static_cast<double>(x.size());
    std::array<double, 3> u{x[0][0] - m[0], x[0][1] - m[1], x[0][2] - m[2]};
    const std::array<double, 3> w{n[1] * u[2] - n[2] * u[1], n[2] * u[0] - n[0] * u[2], n[0] * u[1] - n[1] * u[0]};
    std::vector<std::pair<double, size_t>> angle;
    for (size_t i = 0; i < x.size(); ++i) {
      const std::array<double, 3> d{x[i][0] - m[0], x[i][1] - m[1], x[i][2] - m[2]};
      angle.emplace_back(std::atan2(d[0] * w[0] + d[1] * w[1] + d[2] * w[2], d[0] * u[0] + d[1] * u[1] + d[2] * u[2]), i);
    }
    std::sort(angle.begin(), angle.end());
    Polygon out;
    for (const auto & [a, i]: angle) out.vertices.push_back(fv[i]);
    // the edge plane between consecutive vertices: another plane both are incident to
    auto incident_of = [&](const Point & p) -> const std::set<int> & {
      for (const auto & v: vs) if (geom_.same(v.point, p)) return v.incident;
      throw std::runtime_error("face vertex not found");
    };
    for (size_t k = 0; k < out.vertices.size(); ++k) {
      const auto & a = incident_of(out.vertices[k]);
      const auto & b = incident_of(out.vertices[(k + 1) % out.vertices.size()]);
      int edge{-1};
      for (int i: a) if (i != face && i >= 6 && b.count(i)) { edge = i; break; }
      if (edge < 0) throw std::runtime_error("no plane shared by consecutive face vertices");
      out.edges.push_back(all[edge]);
    }
    return out;
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

  /*! polygon ∩ {n·x <= d}, exactly; empty if the result has no area */
  [[nodiscard]] Polygon clip(const Polygon & poly, const Plane & face, const Plane & h) const {
    const size_t n = poly.vertices.size();
    std::vector<int> s(n);
    bool out_any{false}, in_any{false};
    for (size_t k = 0; k < n; ++k) { s[k] = geom_.side(poly.vertices[k], h); out_any |= s[k] > 0; in_any |= s[k] < 0; }
    if (!out_any) return poly;
    if (!in_any) return {};
    Polygon out;
    for (size_t k = 0; k < n; ++k) {
      const size_t l = (k + 1) % n;
      const auto & P = poly.vertices[k];
      const auto & E = poly.edges[k];
      if (s[k] < 0 || (s[k] == 0 && s[l] <= 0)) { out.vertices.push_back(P); out.edges.push_back(E); }
      else if (s[k] == 0) { out.vertices.push_back(P); out.edges.push_back(h); }
      if ((s[k] < 0 && s[l] > 0) || (s[k] > 0 && s[l] < 0)) {
        out.vertices.push_back({{face, E, h}});
        out.edges.push_back(s[k] < 0 ? h : E);
      }
    }
    if (out.vertices.size() < 3) return {};
    return out;
  }

  /*! Each face F_j split into the convex cells F_j ∩ g(F_i) of positive area */
  void pairing_cells() {
    cells_.assign(planes_.size(), {});
    cell_maps_.assign(planes_.size(), {});
    for (const auto & m: maps_) {
      const auto & src = faces_[static_cast<size_t>(m.from)];
      Polygon cell = faces_[static_cast<size_t>(m.to)];
      for (const auto & e: src.edges) {
        cell = clip(cell, planes_[static_cast<size_t>(m.to)], e.mapped(m.g, m.t));
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
