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
#include <iterator>
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
  std::vector<std::set<int>> on_planes_;              // boundary planes each vertex lies on
  std::vector<std::array<size_t, 4>> tetrahedra_;
  // refinement
  using edge_t = std::pair<size_t, size_t>;
  using tet_t = std::array<size_t, 4>;
  std::set<tet_t> tets_;
  std::map<edge_t, std::set<tet_t>> edge_tets_;
  std::map<std::array<long long, 3>, std::vector<size_t>> positions_;
  using key_t = std::array<long long, 6>;
  std::map<edge_t, key_t> class_cache_;
  double min_edge2_{0};
  double tie_{1e-9};
  size_t self_paired_ties_{0};
  std::map<std::array<long long, 3>, std::vector<size_t>> buckets_;
  std::map<std::vector<size_t>, size_t> centroids_;
  size_t clipped_{0};
  std::vector<std::array<double, 4>> plane_values_;   // boundary planes n·x <= d as doubles, for a filter
  std::map<point, size_t> grid_ids_;                  // the vertex at each grid point already used
  std::vector<mat3i> ops_;
  std::array<double, 3> lo_{}, hi_{};                 // bounding box of the vertices

public:
  /*! \param metric the reciprocal metric in the primitive reciprocal basis
      \param ops the point group, acting on those coordinates
      \param n the grid lattice is Λ* divided by n
      \param cone the wedge (see Boundary) */
  LatticeTri(const std::array<double, 9> & metric, const std::vector<mat3i> & ops, const long long n,
             std::optional<std::vector<int3>> cone = std::nullopt)
      : boundary_(metric, ops, std::move(cone)), denominator_(scale * n) {
    factor();
    std::array<double, 9> rows{};
    for (int i = 0; i < 9; ++i) rows[i] = cholesky_[i] / static_cast<double>(n);
    const Grid grid(rows, ops);
    for (const auto & q: boundary_.planes()) {
      const auto c = q.coefficients(geometry().metric());
      plane_values_.push_back({c[0].estimate(), c[1].estimate(), c[2].estimate(), c[3].estimate()});
    }
    // Only grid tetrahedra that can meet the zone: the pattern's cells over the zone's
    // bounding box, in grid coordinates (n x), widened by the pattern's reach and one cell
    // for round-off. They are made and processed one at a time.
    std::array<double, 3> zlo{}, zhi{};
    bool first{true};
    for (const auto & f: boundary_.faces())
      for (const auto & p: f.vertices) {
        const auto x = geometry().coordinates(p);
        for (int i = 0; i < 3; ++i) {
          const double y = x[i] * static_cast<double>(n);
          zlo[i] = first ? y : std::min(zlo[i], y);
          zhi[i] = first ? y : std::max(zhi[i], y);
        }
        first = false;
      }
    int3 reach_lo{0, 0, 0}, reach_hi{0, 0, 0};
    for (const auto & t: grid.pattern())
      for (const auto & p: t)
        for (int i = 0; i < 3; ++i) {
          reach_lo[i] = std::min(reach_lo[i], floor_div(p[i], scale));
          reach_hi[i] = std::max(reach_hi[i], -floor_div(-p[i], scale));
        }
    std::array<long long, 3> from{}, to{};
    for (int i = 0; i < 3; ++i) {
      from[i] = static_cast<long long>(std::floor(zlo[i])) - reach_hi[i] - 1;
      to[i] = static_cast<long long>(std::ceil(zhi[i])) - reach_lo[i] + 1;
    }
    for (long long a = from[0]; a <= to[0]; ++a)
      for (long long b = from[1]; b <= to[1]; ++b)
        for (long long c = from[2]; c <= to[2]; ++c)
          for (const auto & t: grid.pattern()) {
            tetrahedron s{};
            for (int k = 0; k < 4; ++k) s[k] = Grid::shifted(t[k], {a, b, c});
            process(s);
          }
    ops_ = ops;
    lo_ = hi_ = vertices_.front();
    for (const auto & v: vertices_) for (int i = 0; i < 3; ++i) { lo_[i] = std::min(lo_[i], v[i]); hi_[i] = std::max(hi_[i], v[i]); }
  }
  [[nodiscard]] const Geometry & geometry() const { return boundary_.geometry(); }
  [[nodiscard]] const Boundary & boundary() const { return boundary_; }
  [[nodiscard]] const std::vector<std::array<double, 3>> & vertices() const { return vertices_; }
  [[nodiscard]] const std::vector<std::array<size_t, 4>> & tetrahedra() const { return tetrahedra_; }
  [[nodiscard]] size_t clipped() const { return clipped_; }
  [[nodiscard]] size_t self_paired_ties() const { return self_paired_ties_; }
  //! The boundary planes a vertex lies on
  [[nodiscard]] const std::set<int> & planes_of(const size_t v) const { return on_planes_[v]; }

  /*! \brief Refine the marked tetrahedra once each (with closure)

  Inside the zone this is longest-edge bisection with longest-edge propagation.
  On its boundary, a triangle splits only along its own longest edge (ties broken
  by the edge's equivalence class), and splitting a boundary edge splits every
  equivalent boundary edge too, so paired faces keep matching. A tetrahedron whose
  longest edge is no longer than twice `min_edge` (Cartesian, Å⁻¹) is not split, so
  refinement never makes a tetrahedron smaller than the resolution limit.
  */
  void refine(const std::vector<tet_t> & marked, const double min_edge = 0) {
    // the refinement indexes, built on first use: a mesh that is never refined doesn't need them
    if (tets_.empty()) {
      for (size_t i = 0; i < vertices_.size(); ++i) index_position(i);
      for (auto t: tetrahedra_) add(t);
    }
    min_edge2_ = min_edge * min_edge;
    std::vector<tet_t> todo;
    for (auto t: marked) { std::sort(t.begin(), t.end()); if (tets_.count(t) && refinable(t)) todo.push_back(t); }
    std::sort(todo.begin(), todo.end());
    for (const auto & t: todo) if (tets_.count(t)) split_edge(chosen(t));
    tetrahedra_.assign(tets_.begin(), tets_.end());
  }

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
    std::set<int> on;
    for (size_t j = 0; j < boundary_.planes().size(); ++j) if (geometry().on(p, boundary_.planes()[j])) on.insert(static_cast<int>(j));
    on_planes_.push_back(on);
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
    on_planes_.push_back(common_planes(ids));
    centroids_[ids] = id;
    return id;
  }

  /*! The side of grid point y (named `exact`) of boundary plane j: decided in doubles
  when the value is far from zero compared with its round-off, else exactly */
  [[nodiscard]] int grid_side(const point & y, const Point & exact, const size_t j) const {
    const auto & c = plane_values_[j];
    const auto D = static_cast<double>(denominator_);
    double value{-c[3]}, size{std::abs(c[3])};
    for (int i = 0; i < 3; ++i) {
      const double term = c[i] * (static_cast<double>(y[i]) / D);
      value += term;
      size += std::abs(term);
    }
    if (std::abs(value) > 1e-12 * size) return value > 0 ? 1 : -1;
    return geometry().side(exact, boundary_.planes()[j]);
  }

  //! The vertex at grid point y (named `p`): found by its integer coordinates once it is known
  size_t grid_vertex(const point & y, const Point & p) {
    if (const auto it = grid_ids_.find(y); it != grid_ids_.end()) return it->second;
    const auto id = named(p);
    grid_ids_.emplace(y, id);
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
      for (int k = 0; k < 4; ++k) { const int s = grid_side(t[k], corners[k], j); out += s > 0; on += s == 0; }
      if (out == 4 || (out + on == 4 && out > 0)) return;        // outside, or touching only from outside
      if (out > 0) crossing.push_back(j);
      touching |= on > 0;
    }
    if (crossing.empty() && !touching) {
      tetrahedra_.push_back({grid_vertex(t[0], corners[0]), grid_vertex(t[1], corners[1]), grid_vertex(t[2], corners[2]),
                             grid_vertex(t[3], corners[3])});
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
    // the candidates on the face, found once for all its edges
    std::vector<std::pair<const Point *, std::array<double, 3>>> on_face;
    for (const auto & x: candidates) if (geometry().on(x, face)) on_face.emplace_back(&x, geometry().coordinates(x));
    for (size_t k = 0; k < n; ++k) {
      const auto & a = poly.vertices[k];
      const auto & b = poly.vertices[(k + 1) % n];
      out.push_back(a);
      std::vector<std::pair<double, Point>> mids;
      const auto xa = geometry().coordinates(a), xb = geometry().coordinates(b);
      for (const auto & [candidate, xx]: on_face) {
        const auto & x = *candidate;
        if (!geometry().on(x, poly.edges[k])) continue;
        if (geometry().same(x, a) || geometry().same(x, b)) continue;
        // strictly between a and b: inside every other edge of the polygon, and in the piece
        bool inside{true};
        for (size_t l = 0; l < n && inside; ++l) if (l != k) inside = geometry().side(x, poly.edges[l]) <= 0;
        if (!inside || !within(P, x)) continue;
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

  // --- refinement -------------------------------------------------------------
  [[nodiscard]] std::set<int> common_planes(const std::vector<size_t> & ids) const {
    std::set<int> out = on_planes_[ids[0]];
    for (size_t k = 1; k < ids.size(); ++k) {
      std::set<int> both;
      std::set_intersection(out.begin(), out.end(), on_planes_[ids[k]].begin(), on_planes_[ids[k]].end(), std::inserter(both, both.begin()));
      out.swap(both);
    }
    return out;
  }
  static std::array<long long, 3> position_key(const std::array<double, 3> & x) {
    return {std::llround(x[0] * 1e7), std::llround(x[1] * 1e7), std::llround(x[2] * 1e7)};
  }
  void index_position(const size_t id) { positions_[position_key(vertices_[id])].push_back(id); }
  //! The vertex at x (Λ* coordinates), if any: vertices are far apart compared with round-off
  [[nodiscard]] std::optional<size_t> vertex_at(const std::array<double, 3> & x) const {
    const auto k = position_key(x);
    for (long long a = -1; a <= 1; ++a)
      for (long long b = -1; b <= 1; ++b)
        for (long long c = -1; c <= 1; ++c) {
          auto it = positions_.find({k[0] + a, k[1] + b, k[2] + c});
          if (it == positions_.end()) continue;
          for (const auto id: it->second) {
            double d{0};
            for (int i = 0; i < 3; ++i) d = std::max(d, std::abs(vertices_[id][i] - x[i]));
            if (d <= 1e-9) return id;
          }
        }
    return std::nullopt;
  }
  static edge_t edge(const size_t a, const size_t b) { return a < b ? edge_t{a, b} : edge_t{b, a}; }
  static std::array<edge_t, 6> edges(const tet_t & t) {
    return {edge(t[0], t[1]), edge(t[0], t[2]), edge(t[0], t[3]), edge(t[1], t[2]), edge(t[1], t[3]), edge(t[2], t[3])};
  }
  void add(tet_t t) {
    std::sort(t.begin(), t.end());
    tets_.insert(t);
    for (const auto & e: edges(t)) edge_tets_[e].insert(t);
  }
  void remove(const tet_t & t) {
    tets_.erase(t);
    for (const auto & e: edges(t)) {
      auto it = edge_tets_.find(e);
      it->second.erase(t);
      if (it->second.empty()) edge_tets_.erase(it);
    }
  }
  [[nodiscard]] double length2(const edge_t & e) const {
    std::array<double, 3> d{};
    for (int i = 0; i < 3; ++i) d[i] = vertices_[e.first][i] - vertices_[e.second][i];
    const auto & G = geometry().metric().values();
    double out{0};
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) out += d[i] * G[3 * i + j] * d[j];
    return out;
  }
  template <size_t N>
  [[nodiscard]] std::vector<edge_t> longest(const std::array<edge_t, N> & es) const {
    double top{0};
    for (const auto & e: es) top = std::max(top, length2(e));
    std::vector<edge_t> out;
    for (const auto & e: es) if (length2(e) >= top * (1 - tie_)) out.push_back(e);
    return out;
  }
  [[nodiscard]] bool refinable(const tet_t & t) const {
    double top{0};
    for (const auto & e: edges(t)) top = std::max(top, length2(e));
    return top > 4 * min_edge2_;
  }
  [[nodiscard]] bool boundary_edge(const edge_t & e) const { return !common_planes({e.first, e.second}).empty(); }
  //! a face of t lying in a boundary plane, used by t alone
  [[nodiscard]] bool boundary_triangle(const std::array<size_t, 3> & f) const {
    if (common_planes({f[0], f[1], f[2]}).empty()) return false;
    const auto it = edge_tets_.find(edge(f[0], f[1]));
    if (it == edge_tets_.end()) return false;
    int n{0};
    for (const auto & t: it->second) if (std::count(t.begin(), t.end(), f[2])) ++n;
    return n == 1;
  }
  [[nodiscard]] std::vector<std::array<size_t, 3>> boundary_triangles(const tet_t & t) const {
    std::vector<std::array<size_t, 3>> out;
    for (int skip = 0; skip < 4; ++skip) {
      std::array<size_t, 3> f{};
      int j{0};
      for (int k = 0; k < 4; ++k) if (k != skip) f[j++] = t[k];
      if (boundary_triangle(f)) out.push_back(f);
    }
    return out;
  }
  /*! The boundary edges equivalent to e (e included): its images x -> g x + t, for
  every operation and lattice translation, that are boundary edges of the mesh */
  [[nodiscard]] std::vector<edge_t> orbit(const edge_t & e) const {
    std::set<edge_t> seen{e};
    const auto & a = vertices_[e.first];
    const auto & b = vertices_[e.second];
    for (const auto & g: ops_) {
      std::array<double, 3> ga{}, gb{};
      for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k) {
          ga[i] += static_cast<double>(g[3 * i + k]) * a[k];
          gb[i] += static_cast<double>(g[3 * i + k]) * b[k];
        }
      // translations bringing both images inside the bounding box of the vertices
      std::array<long long, 3> lo{}, hi{};
      for (int i = 0; i < 3; ++i) {
        lo[i] = static_cast<long long>(std::ceil(std::max(lo_[i] - ga[i], lo_[i] - gb[i]) - 1e-9));
        hi[i] = static_cast<long long>(std::floor(std::min(hi_[i] - ga[i], hi_[i] - gb[i]) + 1e-9));
      }
      for (long long t0 = lo[0]; t0 <= hi[0]; ++t0)
        for (long long t1 = lo[1]; t1 <= hi[1]; ++t1)
          for (long long t2 = lo[2]; t2 <= hi[2]; ++t2) {
            const auto p = vertex_at({ga[0] + t0, ga[1] + t1, ga[2] + t2});
            if (!p) continue;
            const auto q = vertex_at({gb[0] + t0, gb[1] + t1, gb[2] + t2});
            if (!q) continue;
            const auto y = edge(*p, *q);
            if (edge_tets_.count(y) && boundary_edge(y)) seen.insert(y);
          }
    }
    return {seen.begin(), seen.end()};
  }
  /*! An invariant label for e's equivalence class, from its position alone: the
  smallest, over the operations, of the image's direction (up to sign) and its
  midpoint modulo the lattice, both on a fine grid. Equivalent edges share it
  whether or not their images are (still) edges of the mesh. */
  [[nodiscard]] key_t class_key(const edge_t & e) {
    auto it = class_cache_.find(e);
    if (it != class_cache_.end()) return it->second;
    constexpr double fine{1e7};
    const auto q = static_cast<long long>(fine);
    const auto & a = vertices_[e.first];
    const auto & b = vertices_[e.second];
    std::optional<key_t> best;
    for (const auto & g: ops_) {
      std::array<long long, 3> d{}, c{};
      for (int i = 0; i < 3; ++i) {
        double ga{0}, gb{0};
        for (int k = 0; k < 3; ++k) { ga += static_cast<double>(g[3 * i + k]) * a[k]; gb += static_cast<double>(g[3 * i + k]) * b[k]; }
        d[i] = std::llround((gb - ga) * fine);
        c[i] = ((std::llround((ga + gb) * fine / 2) % q) + q) % q;
      }
      const std::array<long long, 3> nd{-d[0], -d[1], -d[2]};
      if (nd > d) d = nd;
      const key_t k{d[0], d[1], d[2], c[0], c[1], c[2]};
      if (!best || k < *best) best = k;
    }
    return class_cache_[e] = *best;
  }
  [[nodiscard]] edge_t chosen_2d(const std::array<size_t, 3> & f) {
    const std::array<edge_t, 3> es{edge(f[0], f[1]), edge(f[0], f[2]), edge(f[1], f[2])};
    auto cands = longest(es);
    if (cands.size() > 1) {
      std::vector<std::pair<key_t, edge_t>> keyed;
      for (const auto & e: cands) keyed.emplace_back(class_key(e), e);
      std::sort(keyed.begin(), keyed.end());
      if (keyed[0].first == keyed[1].first) ++self_paired_ties_;
      return keyed[0].second;
    }
    return cands[0];
  }
  [[nodiscard]] edge_t chosen(const tet_t & t) {
    const auto cands = longest(edges(t));
    const auto faces = boundary_triangles(t);
    std::vector<edge_t> ok;
    for (const auto & e: cands) {
      bool good{true};
      for (const auto & f: faces)
        if (std::count(f.begin(), f.end(), e.first) && std::count(f.begin(), f.end(), e.second)) good &= chosen_2d(f) == e;
      if (good) ok.push_back(e);
    }
    return ok.empty() ? *std::min_element(cands.begin(), cands.end()) : *std::min_element(ok.begin(), ok.end());
  }
  void split_edge(const edge_t & e, const int depth = 0, const bool synchronized = false) {
    if (depth > 1000) throw std::runtime_error("refinement propagated too far");
    if (!edge_tets_.count(e)) return;
    const bool on_boundary = boundary_edge(e);
    // e must be the refinement edge of its boundary triangles and of every tetrahedron
    // around it. Splitting another edge can change either, so recheck both each time.
    while (edge_tets_.count(e)) {
      std::optional<edge_t> other;
      if (on_boundary)
        for (const auto & t: edge_tets_.at(e)) {
          for (const auto & f: boundary_triangles(t))
            if (std::count(f.begin(), f.end(), e.first) && std::count(f.begin(), f.end(), e.second)) {
              const auto c = chosen_2d(f);
              if (c != e) { other = c; break; }
            }
          if (other) break;
        }
      if (!other)
        for (const auto & t: edge_tets_.at(e)) { const auto c = chosen(t); if (c != e) { other = c; break; } }
      if (!other) break;
      split_edge(*other, depth + 1);
    }
    if (!edge_tets_.count(e)) return;
    // the equivalent boundary edges, found before e is split
    std::vector<edge_t> partners;
    if (on_boundary && !synchronized) {
      for (const auto & f: orbit(e)) if (f != e) partners.push_back(f);
    }
    // the midpoint
    std::array<double, 3> x{};
    for (int i = 0; i < 3; ++i) x[i] = (vertices_[e.first][i] + vertices_[e.second][i]) / 2;
    const size_t m = vertices_.size();
    vertices_.push_back(x);
    named_.emplace_back(std::nullopt);
    on_planes_.push_back(common_planes({e.first, e.second}));
    index_position(m);
    const auto around = edge_tets_.at(e);
    for (const auto & t: around) {
      remove(t);
      tet_t a = t, b = t;
      for (auto & v: a) if (v == e.second) v = m;
      for (auto & v: b) if (v == e.first) v = m;
      add(a);
      add(b);
    }
    for (const auto & f: partners) split_edge(f, depth + 1, true);
  }
};
}
#endif
