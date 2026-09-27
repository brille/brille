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
#include <atomic>
#include <iterator>
#include <map>
#include <optional>
#include <unordered_map>
#include "lattice_boundary.hpp"
#include "lattice_grid.hpp"
#include "thread_pool.h"

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
  std::map<std::vector<size_t>, size_t> centroids_;
  size_t clipped_{0};
  std::vector<std::array<double, 4>> plane_values_;   // boundary planes n·x <= d as doubles, for a filter
  struct PointHash {
    size_t operator()(const std::array<long long, 3> & p) const {
      return std::hash<long long>()(p[0]) ^ (std::hash<long long>()(p[1]) * 0x9e3779b97f4a7c15ULL) ^ (std::hash<long long>()(p[2]) * 0xc2b2ae3d27d4eb4fULL);
    }
  };
  struct NameHash {
    size_t operator()(const std::array<long long, 21> & k) const {
      size_t h{0xcbf29ce484222325ULL};
      for (const auto v: k) h ^= static_cast<size_t>(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
      return h;
    }
  };
  std::unordered_map<std::array<long long, 3>, std::vector<size_t>, PointHash> buckets_;
  std::unordered_map<point, size_t, PointHash> grid_ids_;   // the vertex at each grid point already used
  std::unordered_map<std::array<long long, 21>, size_t, NameHash> names_;   // the vertex for each point name already used
  std::vector<mat3i> ops_;
  std::array<double, 3> lo_{}, hi_{};                 // bounding box of the vertices

public:
  /*! \param metric the reciprocal metric in the primitive reciprocal basis
      \param ops the point group, acting on those coordinates
      \param n the grid lattice is Λ* divided by n
      \param cone the wedge (see Boundary) */
  LatticeTri(const std::array<double, 9> & metric, const std::vector<mat3i> & ops, const long long n,
             std::optional<std::vector<int3>> cone = std::nullopt)
      : boundary_(metric, ops, std::move(cone)), denominator_(scale * n), ops_(ops) {
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
    pending_t pending;
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
    const auto pool = brille::ThreadPool::getInstance();
    const auto workers = std::max<size_t>(1, pool->size());
    {
      // classified in parallel, one slab of cells at a time, then kept in slab order,
      // so the vertex numbering doesn't depend on the thread count
      const auto slabs = static_cast<size_t>(to[0] - from[0] + 1);
      std::vector<std::vector<tetrahedron>> inside(slabs);
      std::vector<pending_t> queued(slabs);
      std::atomic<size_t> next{0};
      for (size_t w = 0; w < workers; ++w)
        pool->enqueue([&]() {
          std::vector<size_t> crossing;
          for (size_t i = next++; i < slabs; i = next++) {
            const long long a = from[0] + static_cast<long long>(i);
            for (long long b = from[1]; b <= to[1]; ++b)
              for (long long c = from[2]; c <= to[2]; ++c)
                for (const auto & t: grid.pattern()) {
                  tetrahedron s{};
                  for (int k = 0; k < 4; ++k) s[k] = Grid::shifted(t[k], {a, b, c});
                  switch (classify(s, crossing)) {
                    case Kind::inside: inside[i].push_back(s); break;
                    case Kind::boundary: queued[i].emplace_back(s, crossing); break;
                    case Kind::outside: break;
                  }
                }
          }
        });
      pool->wait();
      for (size_t i = 0; i < slabs; ++i) {
        for (const auto & s: inside[i]) keep(s);
        std::vector<tetrahedron>().swap(inside[i]);
        for (auto & q: queued[i]) pending.push_back(std::move(q));
      }
    }
    // The pieces at the boundary are clipped and prepared in parallel, touching
    // nothing shared, then named and split in order, so the result doesn't depend on
    // the thread count.
    std::vector<std::optional<Piece>> pieces(pending.size());
    {
      std::atomic<size_t> next{0};
      for (size_t w = 0; w < workers; ++w)
        pool->enqueue([&]() {
          for (size_t i = next++; i < pending.size(); i = next++) pieces[i] = clip_piece(pending[i].first, pending[i].second);
        });
      pool->wait();
    }
    for (const auto & piece: pieces)
      if (piece) {
        if (piece->clipped) ++clipped_;
        emit(*piece);
      }
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

  //! A point's name up to the order of its planes and their orientations
  using name_t = std::array<long long, 21>;
  static name_t name_of(const Point & p) {
    std::array<std::array<long long, 7>, 3> planes;
    for (int r = 0; r < 3; ++r) {
      const auto & q = p.planes[static_cast<size_t>(r)];
      auto & k = planes[static_cast<size_t>(r)];
      k = {q.kind == Plane::Kind::metric ? 1 : 0, q.a[0], q.a[1], q.a[2], q.b[0], q.b[1], q.b[2]};
      if (q.kind == Plane::Kind::integer) {
        // c·x = k and -c·x = -k are the same plane
        const auto lead = q.a[0] ? q.a[0] : (q.a[1] ? q.a[1] : q.a[2]);
        if (lead < 0) for (int i = 1; i <= 4; ++i) k[static_cast<size_t>(i)] = -k[static_cast<size_t>(i)];
      } else {
        // x·Gσ = uᵀGσ and x·G(-σ) = uᵀG(-σ) are the same plane
        const auto lead = q.a[0] ? q.a[0] : (q.a[1] ? q.a[1] : q.a[2]);
        if (lead < 0) for (int i = 1; i <= 3; ++i) k[static_cast<size_t>(i)] = -k[static_cast<size_t>(i)];
      }
    }
    std::sort(planes.begin(), planes.end());
    name_t out{};
    for (size_t r = 0; r < 3; ++r) for (size_t i = 0; i < 7; ++i) out[7 * r + i] = planes[r][i];
    return out;
  }
  /*! The vertex at p. A point named by the same planes as one already found is that
  vertex, with no arithmetic; otherwise the vertices near p are compared exactly. */
  size_t named(const Point & p) {
    auto name = name_of(p);
    if (const auto it = names_.find(name); it != names_.end()) return it->second;
    const auto id = named_exactly(p);
    names_.emplace(std::move(name), id);
    return id;
  }
  size_t named_exactly(const Point & p) {
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

  //! A face of a piece, with the points inserted on its polygon's edges
  struct PieceFace {
    std::vector<Point> points;                 //!< in cyclic order
    std::vector<size_t> segment_edge;          //!< the polygon edge holding the segment from point k to k+1
    std::vector<Plane> edge_planes;            //!< the plane of each polygon edge
    int plane{-1};                             //!< the piece's plane holding the face
    bool boundary{false};                      //!< whether that plane is a zone boundary plane
  };
  //! A piece of a grid tetrahedron, ready to triangulate
  struct Piece {
    std::vector<Point> corners;
    std::vector<std::set<int>> incident;       //!< the piece planes each corner lies on
    std::vector<PieceFace> faces;
    bool changed{false};                       //!< not simply the grid tetrahedron
    bool clipped{false};
  };

  using pending_t = std::vector<std::pair<tetrahedron, std::vector<size_t>>>;
  enum class Kind { outside, inside, boundary };
  /*! Whether a grid tetrahedron is outside the zone, inside it, or crosses or touches
  its boundary; then `crossing` holds the boundary planes it crosses */
  [[nodiscard]] Kind classify(const tetrahedron & t, std::vector<size_t> & crossing) const {
    std::array<Point, 4> corners;
    for (int k = 0; k < 4; ++k) corners[k] = grid_point(t[k]);
    const auto & planes = boundary_.planes();
    crossing.clear();
    bool touching{false};
    for (size_t j = 0; j < planes.size(); ++j) {
      int out{0}, on{0};
      for (int k = 0; k < 4; ++k) { const int s = grid_side(t[k], corners[k], j); out += s > 0; on += s == 0; }
      if (out == 4 || (out + on == 4 && out > 0)) return Kind::outside;   // outside, or touching only from outside
      if (out > 0) crossing.push_back(j);
      touching |= on > 0;
    }
    return crossing.empty() && !touching ? Kind::inside : Kind::boundary;
  }
  /*! Keep a grid tetrahedron inside the zone. These are kept before any other vertex
  is named, so a grid point not yet seen is a new vertex without comparing it with
  others, and, strictly inside the zone, it is on no boundary plane. */
  void keep(const tetrahedron & t) {
    std::array<size_t, 4> ids{};
    for (int k = 0; k < 4; ++k) {
      const auto & y = t[k];
      if (const auto it = grid_ids_.find(y); it != grid_ids_.end()) { ids[k] = it->second; continue; }
      const auto p = grid_point(y);
      const auto x = geometry().coordinates(p);
      const size_t id = vertices_.size();
      vertices_.push_back(x);
      named_.emplace_back(p);
      on_planes_.emplace_back();
      std::array<long long, 3> key{};
      for (int i = 0; i < 3; ++i) key[i] = std::llround(x[i] * 1e9);
      buckets_[key].push_back(id);
      names_.emplace(name_of(p), id);
      grid_ids_.emplace(y, id);
      ids[k] = id;
    }
    tetrahedra_.push_back(ids);
  }

  //! The part of grid tetrahedron t inside the zone, cut by the boundary planes `crossing`
  [[nodiscard]] std::optional<Piece> clip_piece(const tetrahedron & t, const std::vector<size_t> & crossing) const {
    const auto & planes = boundary_.planes();
    std::array<Plane, 4> fp;
    for (int k = 0; k < 4; ++k) {
      std::array<point, 3> f;
      int j{0};
      for (int l = 0; l < 4; ++l) if (l != k) f[j++] = t[l];
      fp[k] = face_plane(f[0], f[1], f[2], t[k]);
    }
    // corners named as grid points, which is how the vertices of whole tetrahedra are named
    std::array<Point, 4> corners;
    for (int k = 0; k < 4; ++k) corners[k] = grid_point(t[k]);
    auto P = Polytope::tetrahedron(geometry(), fp, corners);
    for (const auto j: crossing) P.cut(planes[j]);
    if (P.flat()) return std::nullopt;
    auto piece = prepare(P);
    piece.clipped = !crossing.empty();
    return piece;
  }

  //! whether x lies on the polytope (inside or on every plane)
  [[nodiscard]] bool within(const Polytope & P, const Point & x) const {
    return std::all_of(P.planes().begin(), P.planes().end(), [&](const Plane & p) { return geometry().side(x, p) <= 0; });
  }

  //! The polygon with the candidate points lying strictly inside its edges inserted
  [[nodiscard]] PieceFace with_points(const Polytope & P, const Polygon & poly, const Plane & face, const std::vector<Point> & candidates) const {
    PieceFace out;
    out.edge_planes = poly.edges;
    const size_t n = poly.vertices.size();
    // the candidates on the face, found once for all its edges
    std::vector<std::pair<const Point *, std::array<double, 3>>> on_face;
    for (const auto & x: candidates) if (geometry().on(x, face)) on_face.emplace_back(&x, geometry().coordinates(x));
    for (size_t k = 0; k < n; ++k) {
      const auto & a = poly.vertices[k];
      const auto & b = poly.vertices[(k + 1) % n];
      out.points.push_back(a);
      out.segment_edge.push_back(k);
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
        if (!geometry().same(x, out.points.back())) {
          out.points.push_back(x);
          out.segment_edge.push_back(k);
        }
    }
    return out;
  }

  //! Fan a polygon (vertex ids) from its centroid; a triangle is kept
  void fan(const std::vector<size_t> & poly, std::vector<std::array<size_t, 3>> & triangles) {
    if (poly.size() == 3) { triangles.push_back({poly[0], poly[1], poly[2]}); return; }
    const size_t c = centroid(poly);
    for (size_t k = 0; k < poly.size(); ++k) triangles.push_back({poly[k], poly[(k + 1) % poly.size()], c});
  }

  //! Everything a piece needs before its vertices are named: this changes nothing shared
  [[nodiscard]] Piece prepare(const Polytope & P) const {
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
    struct Part { Plane face; Polygon polygon; size_t original; int plane; bool boundary; };
    std::vector<Part> parts;
    bool changed = P.vertices().size() != 4;
    for (const auto & [f, bi]: faces) {
      const auto & face = P.planes()[static_cast<size_t>(f)];
      const auto poly = P.polygon(f);
      if (bi < 0) { parts.push_back({face, poly, poly.vertices.size(), f, false}); continue; }
      size_t count{0};
      for (const auto & cell: cells[static_cast<size_t>(bi)]) {
        Polygon q = poly;
        for (const auto & ce: cell.edges) {
          q = exact::clip(geometry(), q, face, ce);
          if (q.vertices.empty()) break;
        }
        if (q.vertices.empty()) continue;
        for (const auto & v: q.vertices) candidates.push_back(v);
        parts.push_back({face, q, poly.vertices.size(), f, true});
        ++count;
      }
      if (count > 1) changed = true;
    }
    Piece piece;
    for (const auto & part: parts) {
      auto face = with_points(P, part.polygon, part.face, candidates);
      if (face.points.size() != part.original) changed = true;
      face.plane = part.plane;
      face.boundary = part.boundary;
      piece.faces.push_back(std::move(face));
    }
    for (const auto & v: P.vertices()) { piece.corners.push_back(v.point); piece.incident.push_back(v.incident); }
    piece.changed = changed;
    return piece;
  }

  //! Whether no chord from point m runs along the face polygon's side: neither of m's
  //! neighbours is followed by a point in line with m (decided exactly)
  [[nodiscard]] bool fan_eligible(const PieceFace & face, const size_t m) const {
    const size_t n = face.points.size();
    if (n <= 3) return true;
    auto in_line = [&](const size_t segment) {   // m on the line of the segment from `segment` to the next point
      return geometry().on(face.points[m], face.edge_planes[face.segment_edge[segment]]);
    };
    return !in_line((m + 1) % n) && !in_line((m + n - 2) % n);
  }
  /*! The point to fan a face inside the zone from: the smallest-id eligible point.
  None if no point qualifies; the face is then fanned from its centroid. Both pieces
  sharing the face choose alike. */
  [[nodiscard]] std::optional<size_t> fan_vertex(const PieceFace & face, const std::vector<size_t> & ids) const {
    std::optional<size_t> best;
    for (size_t m = 0; m < ids.size(); ++m)
      if ((!best || ids[m] < ids[*best]) && fan_eligible(face, m)) best = m;
    return best;
  }
  /*! A label for a point that equivalent points share: the smallest, over the point
  group, of its image modulo the lattice, on a fine grid */
  [[nodiscard]] std::array<long long, 3> point_key(const std::array<double, 3> & x) const {
    constexpr double fine{1e7};
    const auto q = static_cast<long long>(fine);
    std::optional<std::array<long long, 3>> best;
    for (const auto & g: ops_) {
      std::array<long long, 3> k{};
      for (int i = 0; i < 3; ++i) {
        double y{0};
        for (int j = 0; j < 3; ++j) y += static_cast<double>(g[3 * i + j]) * x[j];
        k[i] = ((std::llround(y * fine) % q) + q) % q;
      }
      if (!best || k < *best) best = k;
    }
    return *best;
  }
  /*! The point to fan a face on the zone boundary from: the eligible point of
  smallest `point_key`. The face paired with it maps points to points with equal
  keys, so it chooses the image point and the two fans match. None if no point
  qualifies, or if two eligible points share the smallest key (a face paired with
  itself); the face is then fanned from its centroid. */
  [[nodiscard]] std::optional<size_t> boundary_fan_vertex(const PieceFace & face, const std::vector<size_t> & ids) const {
    std::optional<size_t> best;
    std::array<long long, 3> best_key{};
    bool tied{false};
    for (size_t m = 0; m < ids.size(); ++m) {
      if (!fan_eligible(face, m)) continue;
      const auto key = point_key(vertices_[ids[m]]);
      if (!best || key < best_key) { best = m; best_key = key; tied = false; }
      else if (key == best_key) tied = true;
    }
    if (tied) return std::nullopt;
    return best;
  }

  /*! Fan a face from its point `m`, leaving out the segments on a line through m.
  Whether a segment is on such a line is decided exactly, not from which polygon
  edge m was found on: a point where the polygon runs straight may be a polygon
  vertex in one piece and a point inserted on an edge in its neighbour, and both
  must fan the shared face alike. */
  void vertex_fan(const PieceFace & face, const std::vector<size_t> & ids, const size_t m, std::vector<std::array<size_t, 3>> & triangles) const {
    const size_t n = ids.size();
    for (size_t j = 0; j < n; ++j) {
      const size_t k = (j + 1) % n;
      if (j == m || k == m) continue;
      if (geometry().on(face.points[m], face.edge_planes[face.segment_edge[j]])) continue;
      triangles.push_back({ids[m], ids[j], ids[k]});
    }
  }

  /*! Name a piece's vertices and split it into tetrahedra.

  A face on the zone boundary is fanned from its `boundary_fan_vertex`, which its
  partner face chooses too; any other face from its `fan_vertex`, which the piece on
  its other side chooses too; either from its centroid if it has no fan vertex. The
  piece is a cone from a corner, over the faces not holding that corner, if some
  corner is the fan vertex of every face it is on (the cone fans those faces from
  it); otherwise it is a cone from its centroid. */
  void emit(const Piece & piece) {
    std::vector<size_t> corner_ids;
    for (const auto & p: piece.corners) corner_ids.push_back(named(p));
    if (!piece.changed) {
      tetrahedra_.push_back({corner_ids[0], corner_ids[1], corner_ids[2], corner_ids[3]});
      return;
    }
    std::vector<std::vector<size_t>> face_ids;
    for (const auto & face: piece.faces) {
      std::vector<size_t> ids;
      for (const auto & p: face.points) ids.push_back(named(p));
      face_ids.push_back(std::move(ids));
    }
    std::vector<std::optional<size_t>> fan_from;
    for (size_t f = 0; f < piece.faces.size(); ++f)
      fan_from.push_back(piece.faces[f].boundary ? boundary_fan_vertex(piece.faces[f], face_ids[f]) : fan_vertex(piece.faces[f], face_ids[f]));
    // the apex: the smallest-id corner that qualifies
    std::optional<size_t> apex;
    for (size_t c = 0; c < corner_ids.size(); ++c) {
      bool ok{true};
      for (size_t f = 0; f < piece.faces.size() && ok; ++f) {
        if (!piece.incident[c].count(piece.faces[f].plane)) continue;
        ok = fan_from[f] && face_ids[f][*fan_from[f]] == corner_ids[c];
      }
      if (ok && (!apex || corner_ids[c] < corner_ids[*apex])) apex = c;
    }

    std::vector<std::array<size_t, 3>> triangles;
    for (size_t f = 0; f < piece.faces.size(); ++f) {
      const auto & face = piece.faces[f];
      if (apex && piece.incident[*apex].count(face.plane)) continue;
      const auto & ids = face_ids[f];
      if (!fan_from[f]) { fan(ids, triangles); continue; }
      vertex_fan(face, ids, *fan_from[f], triangles);
    }
    const size_t c = apex ? corner_ids[*apex] : centroid(corner_ids);
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
