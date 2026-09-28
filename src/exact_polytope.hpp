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
#include <optional>
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
  /*! The tetrahedron bounded by four planes, each oriented with the tetrahedron inside.
  Corner k, opposite plane k, is named by the other three planes unless `corners`
  gives it another name. */
  static Polytope tetrahedron(const Geometry & geom, const std::array<Plane, 4> & p,
                              const std::optional<std::array<Point, 4>> & corners = std::nullopt) {
    std::vector<Vertex> vs;
    for (int skip = 3; skip >= 0; --skip) {
      std::array<Plane, 3> t;
      std::set<int> inc;
      int j{0};
      for (int k = 0; k < 4; ++k) if (k != skip) { t[j++] = p[k]; inc.insert(k); }
      vs.push_back({corners ? (*corners)[static_cast<size_t>(skip)] : Point{t}, inc});
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
    // Where three or more planes share a line, a new vertex can lie on planes its edge's
    // ends don't share, and two edges can give the same point: complete each new
    // vertex's planes, and keep one vertex per point
    std::vector<Vertex> unique;
    for (auto & v: created) {
      for (int p = 0; p < k; ++p) if (!v.incident.count(p) && geom_->on(v.point, planes_[static_cast<size_t>(p)])) v.incident.insert(p);
      auto same = std::find_if(unique.begin(), unique.end(), [&](const Vertex & u) { return geom_->same(u.point, v.point); });
      if (same == unique.end()) unique.push_back(v);
      else same->incident.insert(v.incident.begin(), v.incident.end());
    }
    std::vector<Vertex> kept;
    for (size_t i = 0; i < vertices_.size(); ++i) {
      if (side[i] > 0) continue;
      if (side[i] == 0) vertices_[i].incident.insert(k);
      kept.push_back(vertices_[i]);
    }
    for (auto & v: unique) kept.push_back(v);
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

  /*! \brief The face on plane `face`: its vertices in cyclic order, with the plane of each edge

  The order comes from the planes the vertices lie on, not from their coordinates,
  so vertices far closer together than round-off in their coordinates still come
  in the right order. It is counter-clockwise seen from outside, i.e. from the
  side the face plane's normal points to.
  */
  [[nodiscard]] Polygon polygon(const int face) const {
    std::vector<size_t> ids;
    for (size_t i = 0; i < vertices_.size(); ++i) if (vertices_[i].incident.count(face)) ids.push_back(i);
    std::vector<std::array<double, 3>> x;
    for (const auto i: ids) x.push_back(geom_->coordinates(vertices_[i].point));
    // the edges: consecutive face vertices along the line where the face meets another plane
    std::vector<std::vector<std::pair<size_t, int>>> next(ids.size());   // (position in ids, edge plane)
    auto linked = [&](const size_t a, const size_t b) {
      return std::any_of(next[a].begin(), next[a].end(), [&](const auto & e) { return e.first == b; });
    };
    for (int p = 0; p < static_cast<int>(planes_.size()); ++p) {
      if (p == face) continue;
      std::vector<size_t> on;
      for (size_t k = 0; k < ids.size(); ++k) if (vertices_[ids[k]].incident.count(p)) on.push_back(k);
      if (on.size() < 2) continue;
      if (on.size() > 2) {
        // collinear: order along the line
        std::array<double, 3> d{};
        double far{-1};
        for (const auto k: on) {
          const std::array<double, 3> e{x[k][0] - x[on[0]][0], x[k][1] - x[on[0]][1], x[k][2] - x[on[0]][2]};
          const double n = e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
          if (n > far) { far = n; d = e; }
        }
        auto along = [&](const size_t k) { return (x[k][0] - x[on[0]][0]) * d[0] + (x[k][1] - x[on[0]][1]) * d[1] + (x[k][2] - x[on[0]][2]) * d[2]; };
        std::sort(on.begin(), on.end(), [&](const size_t a, const size_t b) { return along(a) < along(b); });
      }
      for (size_t k = 0; k + 1 < on.size(); ++k) {
        if (linked(on[k], on[k + 1])) continue;
        next[on[k]].emplace_back(on[k + 1], p);
        next[on[k + 1]].emplace_back(on[k], p);
      }
    }
    // walk the cycle
    std::vector<size_t> order{0};
    std::vector<int> edges;
    std::vector<bool> used(ids.size(), false);
    used[0] = true;
    for (size_t step = 0; step < ids.size(); ++step) {
      const auto here = order.back();
      auto it = std::find_if(next[here].begin(), next[here].end(), [&](const auto & e) {
        return !used[e.first] || (e.first == 0 && order.size() == ids.size() && order.size() > 2);
      });
      if (it == next[here].end()) throw std::runtime_error("the face's vertices do not form a cycle");
      edges.push_back(it->second);
      if (it->first == 0) break;
      order.push_back(it->first);
      used[it->first] = true;
    }
    if (order.size() != ids.size() || edges.size() != ids.size()) throw std::runtime_error("the face's vertices do not form a cycle");
    // counter-clockwise about the plane's normal
    const auto c = planes_[static_cast<size_t>(face)].coefficients(geom_->metric());
    const std::array<double, 3> n{c[0].estimate(), c[1].estimate(), c[2].estimate()};
    double turn{0};
    for (size_t k = 0; k < order.size(); ++k) {
      const auto & a = x[order[k]];
      const auto & b = x[order[(k + 1) % order.size()]];
      turn += n[0] * (a[1] * b[2] - a[2] * b[1]) + n[1] * (a[2] * b[0] - a[0] * b[2]) + n[2] * (a[0] * b[1] - a[1] * b[0]);
    }
    if (turn < 0) {
      std::reverse(order.begin() + 1, order.end());
      // edge k joins order[k] and order[k+1]: reversing the walk reverses the edges too
      std::vector<int> reversed;
      for (size_t k = 0; k < edges.size(); ++k) reversed.push_back(edges[(edges.size() - 1 - k) % edges.size()]);
      edges = reversed;
    }
    Polygon out;
    for (size_t k = 0; k < order.size(); ++k) {
      out.vertices.push_back(vertices_[ids[order[k]]].point);
      out.edges.push_back(planes_[static_cast<size_t>(edges[k])]);
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
