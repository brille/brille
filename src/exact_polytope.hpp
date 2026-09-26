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
#ifndef BRILLE_EXACT_POLYTOPE_HPP_
#define BRILLE_EXACT_POLYTOPE_HPP_
#include <algorithm>
#include <cmath>
#include <iterator>
#include <set>
#include <vector>
#include "exact_geometry.hpp"

namespace brille::exact {

//! A convex polygon in a plane: its vertices in order, and the plane of each edge (vertex k to k+1)
struct Polygon {
  std::vector<Point> vertices;
  std::vector<Plane> edges;
};

/*! \brief A convex polytope, exactly: half-spaces n·x <= d and named vertices

Each vertex is named by three of the planes and knows every plane it lies on, so
cutting, adjacency and faces are decided without tolerances.
*/
class Polytope {
public:
  struct Vertex { Point point; std::set<int> incident; };
private:
  const Geometry * geom_;
  std::vector<Plane> planes_;
  std::vector<Vertex> vertices_;
public:
  //! A polytope from its planes and its vertices (with their incident plane indices)
  Polytope(const Geometry & geom, std::vector<Plane> planes, std::vector<Vertex> vertices)
      : geom_(&geom), planes_(std::move(planes)), vertices_(std::move(vertices)) {}
  //! The tetrahedron bounded by four planes, each oriented with the tetrahedron inside
  static Polytope tetrahedron(const Geometry & geom, const std::array<Plane, 4> & p) {
    std::vector<Vertex> vs;
    for (int skip = 3; skip >= 0; --skip) {
      std::array<Plane, 3> t;
      std::set<int> inc;
      int j{0};
      for (int k = 0; k < 4; ++k) if (k != skip) { t[j++] = p[k]; inc.insert(k); }
      vs.push_back({{t}, inc});
    }
    return {geom, {p.begin(), p.end()}, vs};
  }
  [[nodiscard]] const std::vector<Plane> & planes() const { return planes_; }
  [[nodiscard]] const std::vector<Vertex> & vertices() const { return vertices_; }
  [[nodiscard]] const Geometry & geometry() const { return *geom_; }

  //! Intersect with the half-space h (n·x <= d); returns whether anything was cut off
  bool cut(const Plane & h) {
    std::vector<int> side(vertices_.size());
    bool out_any{false};
    for (size_t i = 0; i < vertices_.size(); ++i) { side[i] = geom_->side(vertices_[i].point, h); out_any |= side[i] > 0; }
    const int k = static_cast<int>(planes_.size());
    planes_.push_back(h);
    if (!out_any) {
      for (size_t i = 0; i < vertices_.size(); ++i) if (side[i] == 0) vertices_[i].incident.insert(k);
      return false;
    }
    std::vector<Vertex> created;
    for (size_t u = 0; u < vertices_.size(); ++u) {
      if (side[u] <= 0) continue;
      for (size_t v = 0; v < vertices_.size(); ++v) {
        if (side[v] >= 0 || !adjacent(u, v)) continue;
        const auto common = shared(vertices_[u], vertices_[v]);
        bool made{false};
        for (size_t a = 0; a < common.size() && !made; ++a)
          for (size_t b = a + 1; b < common.size() && !made; ++b) {
            Point x{{planes_[static_cast<size_t>(common[a])], planes_[static_cast<size_t>(common[b])], h}};
            if (geom_->independent(x) == 0) continue;
            std::set<int> inc(common.begin(), common.end());
            inc.insert(k);
            created.push_back({x, inc});
            made = true;
          }
      }
    }
    std::vector<Vertex> kept;
    for (size_t i = 0; i < vertices_.size(); ++i) {
      if (side[i] > 0) continue;
      if (side[i] == 0) vertices_[i].incident.insert(k);
      kept.push_back(vertices_[i]);
    }
    for (auto & v: created) kept.push_back(v);
    vertices_.swap(kept);
    return true;
  }

  //! Whether all vertices lie on one plane: no volume
  [[nodiscard]] bool flat() const {
    if (vertices_.size() < 4) return true;
    for (size_t p = 0; p < planes_.size(); ++p)
      if (std::all_of(vertices_.begin(), vertices_.end(), [&](const Vertex & v) { return v.incident.count(static_cast<int>(p)); })) return true;
    return false;
  }

  //! Indices of the planes holding a face (three or more vertices), skipping the first `skip`
  [[nodiscard]] std::vector<int> faces(const int skip = 0) const {
    std::vector<int> out;
    for (int p = skip; p < static_cast<int>(planes_.size()); ++p) {
      int n{0};
      for (const auto & v: vertices_) n += v.incident.count(p) ? 1 : 0;
      if (n >= 3) out.push_back(p);
    }
    return out;
  }

  //! The face on plane `face`: its vertices in cyclic order, with the plane of each edge
  [[nodiscard]] Polygon polygon(const int face) const {
    std::vector<size_t> ids;
    for (size_t i = 0; i < vertices_.size(); ++i) if (vertices_[i].incident.count(face)) ids.push_back(i);
    const auto c = planes_[static_cast<size_t>(face)].coefficients(geom_->metric());
    const std::array<double, 3> n{c[0].estimate(), c[1].estimate(), c[2].estimate()};
    std::vector<std::array<double, 3>> x;
    for (const auto i: ids) x.push_back(geom_->coordinates(vertices_[i].point));
    std::array<double, 3> m{0, 0, 0};
    for (const auto & y: x) for (int i = 0; i < 3; ++i) m[i] += y[i] / static_cast<double>(x.size());
    const std::array<double, 3> u{x[0][0] - m[0], x[0][1] - m[1], x[0][2] - m[2]};
    const std::array<double, 3> w{n[1] * u[2] - n[2] * u[1], n[2] * u[0] - n[0] * u[2], n[0] * u[1] - n[1] * u[0]};
    std::vector<std::pair<double, size_t>> angle;
    for (size_t i = 0; i < x.size(); ++i) {
      const std::array<double, 3> d{x[i][0] - m[0], x[i][1] - m[1], x[i][2] - m[2]};
      angle.emplace_back(std::atan2(d[0] * w[0] + d[1] * w[1] + d[2] * w[2], d[0] * u[0] + d[1] * u[1] + d[2] * u[2]), ids[i]);
    }
    std::sort(angle.begin(), angle.end());
    Polygon out;
    std::vector<size_t> order;
    for (const auto & [a, i]: angle) { order.push_back(i); out.vertices.push_back(vertices_[i].point); }
    for (size_t k = 0; k < order.size(); ++k) {
      const auto common = shared(vertices_[order[k]], vertices_[order[(k + 1) % order.size()]]);
      int edge{-1};
      for (const int p: common) if (p != face) { edge = p; break; }
      if (edge < 0) throw std::runtime_error("no plane is shared by consecutive face vertices");
      out.edges.push_back(planes_[static_cast<size_t>(edge)]);
    }
    return out;
  }

private:
  static std::vector<int> shared(const Vertex & a, const Vertex & b) {
    std::vector<int> common;
    std::set_intersection(a.incident.begin(), a.incident.end(), b.incident.begin(), b.incident.end(), std::back_inserter(common));
    return common;
  }
  [[nodiscard]] bool adjacent(const size_t u, const size_t v) const {
    const auto common = shared(vertices_[u], vertices_[v]);
    if (common.size() < 2) return false;
    for (size_t w = 0; w < vertices_.size(); ++w) {
      if (w == u || w == v) continue;
      if (std::includes(vertices_[w].incident.begin(), vertices_[w].incident.end(), common.begin(), common.end())) return false;
    }
    return true;
  }
};

/*! polygon ∩ {n·x <= d}, exactly; empty if the result has no area. The polygon
lies in the plane `face`. */
inline Polygon clip(const Geometry & geom, const Polygon & poly, const Plane & face, const Plane & h) {
  const size_t n = poly.vertices.size();
  std::vector<int> s(n);
  bool out_any{false}, in_any{false};
  for (size_t k = 0; k < n; ++k) { s[k] = geom.side(poly.vertices[k], h); out_any |= s[k] > 0; in_any |= s[k] < 0; }
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
}
#endif
