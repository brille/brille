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
#ifndef BRILLE_LATTICE_TRI_HPP_
#define BRILLE_LATTICE_TRI_HPP_
/*! \file
\brief The structured mesh of the irreducible zone: the grid clipped to the zone

Grid tetrahedra inside the irreducible polyhedron are kept. Those crossing its
boundary are clipped exactly, and each clipped piece (convex) is split into
tetrahedra from its centroid. Faces on the zone boundary are first split into
their pairing cells, and every edge gets the special points and cell-edge
crossings lying on it, so that paired faces are meshed identically (see the
design note). Every vertex but the centroids is a named point, and vertices are
shared between pieces by exact comparison; a centroid is shared by the set of
vertices it averages.
*/
#include <map>
#include <optional>
#include "lattice_boundary.hpp"
#include "lattice_grid.hpp"

namespace brille::latticetri {

class LatticeTri {
  Boundary boundary_;
  long long denominator_;              // Λ* coordinates are grid points / denominator_
  std::array<double, 9> cholesky_{};   // rows: a basis with the metric G
  std::vector<std::array<double, 3>> vertices_;          // Λ* coordinates
  std::vector<std::optional<Point>> named_;
  std::vector<std::array<size_t, 4>> tetrahedra_;
  std::map<std::array<long long, 3>, std::vector<size_t>> buckets_;
  std::map<std::vector<size_t>, size_t> centroids_;
  size_t clipped_{0};

public:
  /*! \param metric the reciprocal metric in the primitive reciprocal basis
      \param ops the point group, acting on those coordinates
      \param n the grid lattice is Λ* divided by n */
  LatticeTri(const std::array<double, 9> & metric, const std::vector<mat3i> & ops, const long long n)
      : boundary_(metric, ops), denominator_(scale * n) {
    factor();
    std::array<double, 9> rows{};
    for (int i = 0; i < 9; ++i) rows[i] = cholesky_[i] / static_cast<double>(n);
    const Grid grid(rows, ops);
    double radius{0};
    for (const auto & f: boundary_.faces())
      for (const auto & p: f.vertices) {
        const auto x = cartesian(geometry().coordinates(p));
        radius = std::max(radius, std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));
      }
    double longest{0};
    for (int i = 0; i < 3; ++i) longest = std::max(longest, std::sqrt(rows[3 * i] * rows[3 * i] + rows[3 * i + 1] * rows[3 * i + 1] + rows[3 * i + 2] * rows[3 * i + 2]));
    for (const auto & t: grid.patch(radius + 3 * longest)) process(t);
  }
  [[nodiscard]] const Geometry & geometry() const { return boundary_.geometry(); }
  [[nodiscard]] const Boundary & boundary() const { return boundary_; }
  [[nodiscard]] const std::vector<std::array<double, 3>> & vertices() const { return vertices_; }
  [[nodiscard]] const std::vector<std::array<size_t, 4>> & tetrahedra() const { return tetrahedra_; }
  [[nodiscard]] size_t clipped() const { return clipped_; }

private:
  void factor() {
    // G = L Lᵀ, L lower triangular; the rows of L have the metric G
    const auto & G = geometry().metric().values();
    auto & L = cholesky_;
    L[0] = std::sqrt(G[0]);
    L[3] = G[3] / L[0];
    L[4] = std::sqrt(G[4] - L[3] * L[3]);
    L[6] = G[6] / L[0];
    L[7] = (G[7] - L[6] * L[3]) / L[4];
    L[8] = std::sqrt(G[8] - L[6] * L[6] - L[7] * L[7]);
  }
  [[nodiscard]] std::array<double, 3> cartesian(const std::array<double, 3> & x) const {
    std::array<double, 3> out{0, 0, 0};
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) out[k] += x[i] * cholesky_[3 * i + k];
    return out;
  }

  //! A grid point as a named point: x_i = y_i / D is the integer plane D x_i = y_i
  [[nodiscard]] Point grid_point(const point & y) const {
    return {{Plane::integer_plane({denominator_, 0, 0}, y[0]), Plane::integer_plane({0, denominator_, 0}, y[1]),
             Plane::integer_plane({0, 0, denominator_}, y[2])}};
  }
  //! The plane through three grid points, oriented with the fourth inside
  [[nodiscard]] Plane face_plane(const point & a, const point & b, const point & c, const point & d) const {
    const int3 u{b[0] - a[0], b[1] - a[1], b[2] - a[2]}, v{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    int3 N{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
    const long long g = std::gcd(std::gcd(std::llabs(N[0]), std::llabs(N[1])), std::llabs(N[2]));
    for (auto & x: N) x /= g;
    long long k = N[0] * a[0] + N[1] * a[1] + N[2] * a[2];
    const long long kd = N[0] * d[0] + N[1] * d[1] + N[2] * d[2];
    if (kd > k) { for (auto & x: N) x = -x; k = -k; }
    // N·y = k with y = D x
    return Plane::integer_plane({N[0] * denominator_, N[1] * denominator_, N[2] * denominator_}, k);
  }

  size_t named(const Point & p) {
    const auto x = geometry().coordinates(p);
    std::array<long long, 3> key{};
    for (int i = 0; i < 3; ++i) key[i] = std::llround(x[i] * 1e9);
    for (long long a = -1; a <= 1; ++a)
      for (long long b = -1; b <= 1; ++b)
        for (long long c = -1; c <= 1; ++c) {
          auto it = buckets_.find({key[0] + a, key[1] + b, key[2] + c});
          if (it == buckets_.end()) continue;
          for (const auto id: it->second) if (named_[id] && geometry().same(*named_[id], p)) return id;
        }
    const size_t id = vertices_.size();
    vertices_.push_back(x);
    named_.emplace_back(p);
    buckets_[key].push_back(id);
    return id;
  }
  size_t centroid(std::vector<size_t> ids) {
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    auto it = centroids_.find(ids);
    if (it != centroids_.end()) return it->second;
    std::array<double, 3> x{0, 0, 0};
    for (const auto i: ids) for (int k = 0; k < 3; ++k) x[k] += vertices_[i][k] / static_cast<double>(ids.size());
    const size_t id = vertices_.size();
    vertices_.push_back(x);
    named_.emplace_back(std::nullopt);
    centroids_[ids] = id;
    return id;
  }

  void process(const tetrahedron & t) {
    std::array<Point, 4> corners;
    for (int k = 0; k < 4; ++k) corners[k] = grid_point(t[k]);
    const auto & planes = boundary_.planes();
    std::vector<size_t> crossing;
    bool touching{false};
    for (size_t j = 0; j < planes.size(); ++j) {
      int out{0}, on{0};
      for (const auto & c: corners) { const int s = geometry().side(c, planes[j]); out += s > 0; on += s == 0; }
      if (out == 4 || (out + on == 4 && out > 0)) return;        // outside, or touching only from outside
      if (out > 0) crossing.push_back(j);
      touching |= on > 0;
    }
    if (crossing.empty() && !touching) {
      tetrahedra_.push_back({named(corners[0]), named(corners[1]), named(corners[2]), named(corners[3])});
      return;
    }
    std::array<Plane, 4> fp;
    for (int k = 0; k < 4; ++k) {
      std::array<point, 3> f;
      int j{0};
      for (int l = 0; l < 4; ++l) if (l != k) f[j++] = t[l];
      fp[k] = face_plane(f[0], f[1], f[2], t[k]);
    }
    auto P = Polytope::tetrahedron(geometry(), fp);
    for (const auto j: crossing) P.cut(planes[j]);
    if (P.flat()) return;
    if (!crossing.empty()) ++clipped_;
    triangulate(P);
  }

  //! whether x lies on the polytope (inside or on every plane)
  [[nodiscard]] bool within(const Polytope & P, const Point & x) const {
    return std::all_of(P.planes().begin(), P.planes().end(), [&](const Plane & p) { return geometry().side(x, p) <= 0; });
  }

  //! Insert into each edge of the polygon the candidate points lying strictly inside it
  [[nodiscard]] std::vector<Point> with_points(const Polytope & P, const Polygon & poly, const Plane & face, const std::vector<Point> & candidates) const {
    std::vector<Point> out;
    const size_t n = poly.vertices.size();
    for (size_t k = 0; k < n; ++k) {
      const auto & a = poly.vertices[k];
      const auto & b = poly.vertices[(k + 1) % n];
      out.push_back(a);
      std::vector<std::pair<double, Point>> mids;
      const auto xa = geometry().coordinates(a), xb = geometry().coordinates(b);
      for (const auto & x: candidates) {
        if (!geometry().on(x, face) || !geometry().on(x, poly.edges[k])) continue;
        if (geometry().same(x, a) || geometry().same(x, b)) continue;
        // strictly between a and b: inside every other edge of the polygon, and in the piece
        bool inside{true};
        for (size_t l = 0; l < n && inside; ++l) if (l != k) inside = geometry().side(x, poly.edges[l]) <= 0;
        if (!inside || !within(P, x)) continue;
        const auto xx = geometry().coordinates(x);
        double d{0};
        for (int i = 0; i < 3; ++i) d += (xx[i] - xa[i]) * (xb[i] - xa[i]);
        mids.emplace_back(d, x);
      }
      std::sort(mids.begin(), mids.end(), [](const auto & u, const auto & v) { return u.first < v.first; });
      // the same point may be found twice (e.g. as a special point and a crossing)
      for (const auto & [d, x]: mids)
        if (!geometry().same(x, out.back())) out.push_back(x);
    }
    return out;
  }

  //! Fan a polygon (vertex ids) from its centroid; a triangle is kept
  void fan(const std::vector<size_t> & poly, std::vector<std::array<size_t, 3>> & triangles) {
    if (poly.size() == 3) { triangles.push_back({poly[0], poly[1], poly[2]}); return; }
    const size_t c = centroid(poly);
    for (size_t k = 0; k < poly.size(); ++k) triangles.push_back({poly[k], poly[(k + 1) % poly.size()], c});
  }

  void triangulate(const Polytope & P) {
    const auto & planes = boundary_.planes();
    const auto & cells = boundary_.cells();
    // candidate points: special points near the piece, cell corners and cell-edge crossings
    std::vector<Point> candidates;
    std::array<double, 3> lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const auto & v: P.vertices()) {
      const auto x = geometry().coordinates(v.point);
      for (int i = 0; i < 3; ++i) { lo[i] = std::min(lo[i], x[i]); hi[i] = std::max(hi[i], x[i]); }
    }
    auto near = [&](const Point & p) {
      const auto x = geometry().coordinates(p);
      for (int i = 0; i < 3; ++i) if (x[i] < lo[i] - 1e-9 || x[i] > hi[i] + 1e-9) return false;
      return true;
    };
    for (const auto & s: boundary_.special_points()) if (near(s)) candidates.push_back(s);
    // faces of the piece, with the boundary plane each lies on, if any
    std::vector<std::pair<int, int>> faces;   // (piece plane, boundary plane or -1)
    for (const int f: P.faces()) {
      int bi{-1};
      for (size_t j = 0; j < planes.size() && bi < 0; ++j) if (geometry().same(P.planes()[static_cast<size_t>(f)], planes[j])) bi = static_cast<int>(j);
      faces.emplace_back(f, bi);
    }
    // cell edges crossing the piece's edges that lie in a boundary plane
    for (size_t j = 0; j < planes.size(); ++j) {
      for (const auto & cell: cells[j])
        for (const auto & ce: cell.edges)
          for (const auto & [f, bi]: faces) {
            if (bi != static_cast<int>(j)) continue;
            const auto poly = P.polygon(f);
            for (const auto & e: poly.edges) {
              Point x{{planes[j], e, ce}};
              if (geometry().independent(x) == 0 || !within(P, x)) continue;
              bool in_cell{true};
              for (const auto & other: cell.edges) in_cell &= geometry().side(x, other) <= 0;
              if (in_cell) candidates.push_back(x);
            }
          }
    }
    // pieces touching a boundary plane only along an edge: that edge still needs the crossings
    for (size_t j = 0; j < planes.size(); ++j) {
      bool has_face{false};
      for (const auto & [f, bi]: faces) has_face |= bi == static_cast<int>(j);
      if (has_face) continue;
      std::vector<size_t> on;
      for (size_t i = 0; i < P.vertices().size(); ++i) if (geometry().on(P.vertices()[i].point, planes[j])) on.push_back(i);
      if (on.size() != 2) continue;
      const auto & va = P.vertices()[on[0]];
      const auto & vb = P.vertices()[on[1]];
      std::vector<int> common;
      std::set_intersection(va.incident.begin(), va.incident.end(), vb.incident.begin(), vb.incident.end(), std::back_inserter(common));
      if (common.size() < 2) continue;
      for (const auto & cell: cells[j])
        for (const auto & ce: cell.edges) {
          Point x{{P.planes()[static_cast<size_t>(common[0])], P.planes()[static_cast<size_t>(common[1])], ce}};
          if (geometry().independent(x) == 0 || !within(P, x) || !geometry().on(x, planes[j])) continue;
          bool in_cell{true};
          for (const auto & other: cell.edges) in_cell &= geometry().side(x, other) <= 0;
          if (in_cell) candidates.push_back(x);
        }
    }
    // the faces, split into cells on the boundary; then, with every candidate known,
    // insert the points on each part's edges and fan it
    struct Part { Plane face; Polygon polygon; size_t original; };
    std::vector<Part> parts;
    bool changed = P.vertices().size() != 4;
    for (const auto & [f, bi]: faces) {
      const auto & face = P.planes()[static_cast<size_t>(f)];
      const auto poly = P.polygon(f);
      if (bi < 0) { parts.push_back({face, poly, poly.vertices.size()}); continue; }
      size_t count{0};
      for (const auto & cell: cells[static_cast<size_t>(bi)]) {
        Polygon q = poly;
        for (const auto & ce: cell.edges) {
          q = exact::clip(geometry(), q, face, ce);
          if (q.vertices.empty()) break;
        }
        if (q.vertices.empty()) continue;
        for (const auto & v: q.vertices) candidates.push_back(v);
        parts.push_back({face, q, poly.vertices.size()});
        ++count;
      }
      if (count > 1) changed = true;
    }
    std::vector<std::array<size_t, 3>> triangles;
    for (const auto & part: parts) {
      const auto pts = with_points(P, part.polygon, part.face, candidates);
      if (pts.size() != part.original) changed = true;
      std::vector<size_t> ids;
      for (const auto & p: pts) ids.push_back(named(p));
      fan(ids, triangles);
    }
    std::vector<size_t> corner_ids;
    for (const auto & v: P.vertices()) corner_ids.push_back(named(v.point));
    if (!changed) {
      tetrahedra_.push_back({corner_ids[0], corner_ids[1], corner_ids[2], corner_ids[3]});
      return;
    }
    const size_t c = centroid(corner_ids);
    for (const auto & tri: triangles) tetrahedra_.push_back({tri[0], tri[1], tri[2], c});
  }
};
}
#endif
