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
#include <cstdint>
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
  std::vector<std::optional<Point>> named_;          // while constructing: each vertex's name, if it has one
  std::vector<std::set<int>> on_planes_;              // boundary planes each vertex lies on
  std::vector<std::array<size_t, 4>> tetrahedra_;
  // refinement
  using edge_t = std::pair<size_t, size_t>;
  using tet_t = std::array<size_t, 4>;
  // Refinement's indexes: the tetrahedra (a slot each; freed slots are reused) and,
  // for each edge, the slots of the tetrahedra around it. Kept lean: a refined mesh
  // holds them for as long as it may be refined again.
  std::vector<tet_t> slots_;
  std::vector<std::uint32_t> free_slots_;
  std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> edge_tets_;
  bool indexed_{false};
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
    const auto shapes = cell_shapes(grid);
    int3 reach_lo{0, 0, 0}, reach_hi{0, 0, 0};
    for (const auto & shape: shapes)
      for (const auto & p: shape.vertices)
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
                for (size_t k = 0; k < shapes.size(); ++k)
                  switch (classify(shapes[k], {a, b, c}, crossing)) {
                    case Kind::inside:
                      for (const auto & t: shapes[k].tets) {
                        tetrahedron s{};
                        for (int v = 0; v < 4; ++v) s[v] = at(shapes[k].vertices[t[v]], {a, b, c});
                        inside[i].push_back(s);
                      }
                      break;
                    case Kind::boundary: queued[i].emplace_back(CellAt{k, {a, b, c}}, crossing); break;
                    case Kind::outside: break;
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
    // The cells at the boundary are clipped and prepared in parallel, touching
    // nothing shared, then named and split in order, so the result doesn't depend on
    // the thread count.
    std::vector<std::optional<Piece>> pieces(pending.size());
    {
      std::atomic<size_t> next{0};
      for (size_t w = 0; w < workers; ++w)
        pool->enqueue([&]() {
          for (size_t i = next++; i < pending.size(); i = next++)
            pieces[i] = clip_cell(shapes[pending[i].first.shape], pending[i].first.cell, pending[i].second);
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
    // Only construction looks points up by name or position, or tests them exactly:
    // for 16,000 vertices the lookups hold 12 MiB and the cached point data up to
    // 180 MiB (12 threads), several times what refinement needs.
    decltype(names_)().swap(names_);
    decltype(buckets_)().swap(buckets_);
    decltype(grid_ids_)().swap(grid_ids_);
    decltype(centroids_)().swap(centroids_);
    decltype(named_)().swap(named_);
    geometry().clear_caches();
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
    if (!indexed_) {
      for (size_t i = 0; i < vertices_.size(); ++i) index_position(i);
      for (auto t: tetrahedra_) add(t);
      indexed_ = true;
    }
    min_edge2_ = min_edge * min_edge;
    std::vector<tet_t> todo;
    for (auto t: marked) { std::sort(t.begin(), t.end()); if (has_tet(t) && refinable(t)) todo.push_back(t); }
    std::sort(todo.begin(), todo.end());
    for (const auto & t: todo) if (has_tet(t)) split_edge(chosen(t));
    // in sorted order, as before these indexes were lean
    tetrahedra_.clear();
    std::vector<bool> freed(slots_.size(), false);
    for (const auto f: free_slots_) freed[f] = true;
    for (size_t i = 0; i < slots_.size(); ++i) if (!freed[i]) tetrahedra_.push_back(slots_[i]);
    std::sort(tetrahedra_.begin(), tetrahedra_.end());
    geometry().clear_caches();
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
    std::vector<std::array<Point, 3>> grid_triangles;   //!< if not empty, the face keeps these (the grid's)
    std::optional<Point> preferred;            //!< fan from this point, if one of `points` and eligible
    std::optional<Point> centre;               //!< the grid triangles' shared centre: a possible apex
    int plane{-1};                             //!< the piece's plane holding the face
    bool boundary{false};                      //!< whether that plane is a zone boundary plane
  };
  //! A piece of a grid tetrahedron, ready to triangulate
  struct Piece {
    std::vector<Point> corners;
    std::vector<std::set<int>> incident;       //!< the piece planes each corner lies on
    std::vector<PieceFace> faces;
    bool changed{false};                       //!< not simply the grid cell
    std::vector<std::array<Point, 4>> tets;    //!< the grid cell's tetrahedra, if unchanged
    bool clipped{false};
  };

  /*! A cell of the grid pattern: a convex polytope split into grid tetrahedra, with
  its face planes (integer, in grid coordinates) and the grid's triangles on each */
  struct CellShape {
    std::vector<point> vertices;                   //!< every grid point of its tetrahedra
    std::vector<std::array<size_t, 4>> tets;        //!< its tetrahedra, by vertex
    std::vector<size_t> corners;                    //!< the vertices on three or more faces
    struct Face {
      int3 normal{};                                //!< N·y <= offset inside, y in grid coordinates
      long long offset{0};
      std::vector<size_t> corners;                  //!< its corners
      std::vector<std::array<size_t, 3>> triangles; //!< the grid's triangles on it
      std::optional<size_t> preferred;              //!< the lower end of its diagonal (a Kuhn cell's parallelogram)
      std::optional<size_t> centre;                 //!< the vertex all its grid triangles share (a split face's centre)
    };
    std::vector<Face> faces;
  };
  //! The plane N·y = k through grid points a, b, c, with `inside` on the side N·y < k
  static std::pair<int3, long long> plane_through(const point & a, const point & b, const point & c, const point & inside) {
    const int3 u{b[0] - a[0], b[1] - a[1], b[2] - a[2]}, v{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    int3 N{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
    const long long g = std::gcd(std::gcd(std::llabs(N[0]), std::llabs(N[1])), std::llabs(N[2]));
    for (auto & x: N) x /= g;
    long long k = N[0] * a[0] + N[1] * a[1] + N[2] * a[2];
    if (N[0] * inside[0] + N[1] * inside[1] + N[2] * inside[2] > k) { for (auto & x: N) x = -x; k = -k; }
    return {N, k};
  }
  /*! The cells the grid is clipped by: the larger they are, the fewer grid edges cross
  the zone boundary and the fewer boundary points the mesh has, but they must be as
  symmetric as the grid, or the boundary would not match across paired faces.
  - A degenerate grid's cells are its Delaunay cells, which are.
  - A Kuhn grid's tetrahedra make up the parallelepipeds spanned by three of the four
    superbase vectors; those without v_k are symmetric if every operation takes v_k to
    ±v_k (inversion maps each to a lattice translate of itself).
  - Otherwise (e.g. body- and face-centred cubic grids) each tetrahedron is a cell. */
  static std::vector<std::vector<tetrahedron>> cell_tetrahedra(const Grid & grid) {
    std::vector<std::vector<tetrahedron>> cells;
    if (grid.degenerate()) {
      for (const auto & members: grid.cells()) {
        std::vector<tetrahedron> tets;
        for (const auto m: members) tets.push_back(grid.pattern()[m]);
        cells.push_back(std::move(tets));
      }
      return cells;
    }
    const auto & v = grid.superbase();
    for (size_t k = 0; k < 4; ++k) {
      const bool fixed = std::all_of(grid.operations().begin(), grid.operations().end(), [&](const mat3i & g) {
        const auto w = Grid::act(g, v[k]);
        return w == v[k] || w == point{-v[k][0], -v[k][1], -v[k][2]};
      });
      if (!fixed) continue;
      std::array<size_t, 3> others{};
      size_t j{0};
      for (size_t i = 0; i < 4; ++i) if (i != k) others[j++] = i;
      std::vector<tetrahedron> tets;
      std::array<size_t, 3> perm{0, 1, 2};
      do {
        const auto & a = v[others[perm[0]]];
        const auto & b = v[others[perm[1]]];
        const auto & c = v[others[perm[2]]];
        tetrahedron t{};
        for (int i = 0; i < 3; ++i) {
          t[1][i] = scale * a[i];
          t[2][i] = scale * (a[i] + b[i]);
          t[3][i] = scale * (a[i] + b[i] + c[i]);
        }
        tets.push_back(t);
      } while (std::next_permutation(perm.begin(), perm.end()));
      return {tets};
    }
    for (const auto & t: grid.pattern()) cells.push_back({t});
    return cells;
  }
  static std::vector<CellShape> cell_shapes(const Grid & grid) {
    std::vector<CellShape> out;
    for (const auto & cell: cell_tetrahedra(grid)) {
      CellShape shape;
      std::map<point, size_t> index;
      auto vertex = [&](const point & p) {
        const auto [it, added] = index.emplace(p, shape.vertices.size());
        if (added) shape.vertices.push_back(p);
        return it->second;
      };
      for (const auto & t: cell) shape.tets.push_back({vertex(t[0]), vertex(t[1]), vertex(t[2]), vertex(t[3])});
      // the cell's faces: triangles of its tetrahedra that only one of them has, by plane
      std::map<std::array<size_t, 3>, int> count;
      for (const auto & t: shape.tets)
        for (int skip = 0; skip < 4; ++skip) {
          std::array<size_t, 3> f{};
          int j{0};
          for (int k = 0; k < 4; ++k) if (k != skip) f[j++] = t[k];
          std::sort(f.begin(), f.end());
          ++count[f];
        }
      std::map<std::pair<int3, long long>, size_t> by_plane;
      for (const auto & [f, k]: count) {
        if (k != 1) continue;
        const auto & a = shape.vertices[f[0]], & b = shape.vertices[f[1]], & c = shape.vertices[f[2]];
        // any cell vertex off the triangle's plane is inside
        size_t off{0};
        for (size_t i = 0; i < shape.vertices.size(); ++i) {
          const auto [N, kk] = plane_through(a, b, c, shape.vertices[i]);
          const auto & y = shape.vertices[i];
          if (N[0] * y[0] + N[1] * y[1] + N[2] * y[2] != kk) { off = i; break; }
        }
        const auto key = plane_through(a, b, c, shape.vertices[off]);
        const auto [it, added] = by_plane.emplace(key, shape.faces.size());
        if (added) shape.faces.push_back({key.first, key.second, {}, {}, std::nullopt, std::nullopt});
        shape.faces[it->second].triangles.push_back(f);
      }
      // corners: vertices on three or more face planes
      for (size_t i = 0; i < shape.vertices.size(); ++i) {
        const auto & y = shape.vertices[i];
        size_t on{0};
        for (auto & face: shape.faces)
          if (face.normal[0] * y[0] + face.normal[1] * y[1] + face.normal[2] * y[2] == face.offset) ++on;
        if (on >= 3) shape.corners.push_back(i);
      }
      for (auto & face: shape.faces)
        for (const auto i: shape.corners) {
          const auto & y = shape.vertices[i];
          if (face.normal[0] * y[0] + face.normal[1] * y[1] + face.normal[2] * y[2] == face.offset) face.corners.push_back(i);
        }
      // A Kuhn parallelepiped's six tetrahedra run from its origin to its far corner, and
      // each face holds one of the two, where its diagonal ends. Its preferred point is the
      // diagonal's other end from the far corner: the origin of this cell, and the same
      // point for the cell on the face's other side, which has it at its far-corner end.
      // A face split into three or more triangles from one vertex: that vertex, its centre.
      std::optional<size_t> far;
      if (shape.tets.size() == 6 && std::all_of(shape.tets.begin(), shape.tets.end(), [&](const auto & t) { return t[3] == shape.tets[0][3] && t[0] == shape.tets[0][0]; }))
        far = shape.tets[0][3];
      for (auto & face: shape.faces) {
        if (face.triangles.size() == 2 && far) {
          std::vector<size_t> shared;
          const auto & t0 = face.triangles[0], & t1 = face.triangles[1];
          for (const auto i: t0) if (std::find(t1.begin(), t1.end(), i) != t1.end()) shared.push_back(i);
          if (shared.size() == 2) face.preferred = shared[0] == *far ? shared[1] : shared[0];
        } else if (face.triangles.size() >= 3) {
          for (const auto i: face.triangles[0])
            if (std::all_of(face.triangles.begin(), face.triangles.end(), [&](const auto & t) { return std::find(t.begin(), t.end(), i) != t.end(); }))
              face.centre = i;
        }
      }
      out.push_back(std::move(shape));
    }
    return out;
  }
  //! A cell of the pattern, moved to the grid cell at `cell`
  struct CellAt {
    size_t shape{0};
    int3 cell{};
  };
  [[nodiscard]] static point at(const point & y, const int3 & cell) { return Grid::shifted(y, cell); }

  using pending_t = std::vector<std::pair<CellAt, std::vector<size_t>>>;
  enum class Kind { outside, inside, boundary };
  /*! Whether a grid cell is outside the zone, inside it, or crosses or touches its
  boundary; then `crossing` holds the boundary planes it crosses */
  [[nodiscard]] Kind classify(const CellShape & shape, const int3 & cell, std::vector<size_t> & crossing) const {
    std::vector<point> ys;
    std::vector<Point> corners;
    for (const auto i: shape.corners) { ys.push_back(at(shape.vertices[i], cell)); corners.push_back(grid_point(ys.back())); }
    const auto & planes = boundary_.planes();
    crossing.clear();
    bool touching{false};
    const auto count = static_cast<int>(ys.size());
    for (size_t j = 0; j < planes.size(); ++j) {
      int out{0}, on{0};
      for (int k = 0; k < count; ++k) { const int s = grid_side(ys[static_cast<size_t>(k)], corners[static_cast<size_t>(k)], j); out += s > 0; on += s == 0; }
      if (out == count || (out + on == count && out > 0)) return Kind::outside;   // outside, or touching only from outside
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

  //! What a cut cell's piece needs of the whole cell
  struct CellInfo {
    std::vector<std::array<Point, 4>> tets;                        //!< its grid tetrahedra
    std::vector<std::vector<std::array<Point, 3>>> face_triangles; //!< the grid's triangles on face f (piece plane f)
    std::vector<size_t> face_corners;                              //!< how many corners face f has
    std::vector<std::optional<Point>> face_preferred;              //!< where to fan face f from, if a polygon point
    std::vector<std::optional<Point>> face_centre;                 //!< the centre face f's grid triangles share
    size_t corners{0};
  };
  //! Whether a point is named as a grid point (by integer planes D x_i = y_i)
  [[nodiscard]] bool is_grid_point(const Point & p) const {
    for (int i = 0; i < 3; ++i) {
      const auto & q = p.planes[static_cast<size_t>(i)];
      if (q.kind != Plane::Kind::integer || q.orientation != 1) return false;
      for (int k = 0; k < 3; ++k) if (q.a[static_cast<size_t>(k)] != (k == i ? denominator_ : 0)) return false;
    }
    return true;
  }

  //! The part of a grid cell inside the zone, cut by the boundary planes `crossing`
  [[nodiscard]] std::optional<Piece> clip_cell(const CellShape & shape, const int3 & cell, const std::vector<size_t> & crossing) const {
    const auto & planes = boundary_.planes();
    const int3 shift{cell[0] * scale, cell[1] * scale, cell[2] * scale};
    CellInfo info;
    info.corners = shape.corners.size();
    std::vector<Plane> fp;
    for (const auto & face: shape.faces) {
      const long long k = face.offset + face.normal[0] * shift[0] + face.normal[1] * shift[1] + face.normal[2] * shift[2];
      fp.push_back(Plane::integer_plane({face.normal[0] * denominator_, face.normal[1] * denominator_, face.normal[2] * denominator_}, k));
      std::vector<std::array<Point, 3>> triangles;
      for (const auto & t: face.triangles)
        triangles.push_back({grid_point(at(shape.vertices[t[0]], cell)), grid_point(at(shape.vertices[t[1]], cell)), grid_point(at(shape.vertices[t[2]], cell))});
      info.face_triangles.push_back(std::move(triangles));
      info.face_corners.push_back(face.corners.size());
      info.face_preferred.push_back(face.preferred ? std::optional<Point>(grid_point(at(shape.vertices[*face.preferred], cell))) : std::nullopt);
      info.face_centre.push_back(face.centre ? std::optional<Point>(grid_point(at(shape.vertices[*face.centre], cell))) : std::nullopt);
    }
    for (const auto & t: shape.tets)
      info.tets.push_back({grid_point(at(shape.vertices[t[0]], cell)), grid_point(at(shape.vertices[t[1]], cell)),
                           grid_point(at(shape.vertices[t[2]], cell)), grid_point(at(shape.vertices[t[3]], cell))});
    // the cell as a polytope: its corners, named as grid points, on their face planes
    std::vector<Polytope::Vertex> vertices;
    for (const auto i: shape.corners) {
      std::set<int> incident;
      for (size_t f = 0; f < shape.faces.size(); ++f)
        if (std::find(shape.faces[f].corners.begin(), shape.faces[f].corners.end(), i) != shape.faces[f].corners.end()) incident.insert(static_cast<int>(f));
      vertices.push_back({grid_point(at(shape.vertices[i], cell)), incident});
    }
    Polytope P(geometry(), fp, vertices);
    for (const auto j: crossing) P.cut(planes[j]);
    if (P.flat()) return std::nullopt;
    auto piece = prepare(P, info);
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
  [[nodiscard]] Piece prepare(const Polytope & P, const CellInfo & info) const {
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
    // whether a face is split into boundary cells or has points inserted on its edges
    bool changed{false};
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
      // A whole face of the cell, uncut and with nothing inserted, keeps the grid's
      // triangles: a neighbouring cell that is kept whole has them too
      const auto f = static_cast<size_t>(part.plane);
      if (!part.boundary && f < info.face_triangles.size()) {
        if (face.points.size() == info.face_corners[f]
            && std::all_of(face.points.begin(), face.points.end(), [&](const Point & p) { return is_grid_point(p); })) {
          face.grid_triangles = info.face_triangles[f];
          face.centre = info.face_centre[f];
        }
        face.preferred = info.face_preferred[f];
      }
      piece.faces.push_back(std::move(face));
    }
    for (const auto & v: P.vertices()) { piece.corners.push_back(v.point); piece.incident.push_back(v.incident); }
    // A piece whose faces are neither split nor given more points is a single
    // tetrahedron, or else the whole cell (its corners, all grid points), which keeps the
    // grid's tetrahedra; anything else is split from a vertex or centroid (emit)
    const bool whole_cell = P.vertices().size() == info.corners
                            && std::all_of(P.vertices().begin(), P.vertices().end(), [&](const auto & v) { return is_grid_point(v.point); });
    piece.changed = changed || !(whole_cell || P.vertices().size() == 4);
    if (!piece.changed) {
      if (whole_cell) piece.tets = info.tets;
      else piece.tets = {{P.vertices()[0].point, P.vertices()[1].point, P.vertices()[2].point, P.vertices()[3].point}};
    }
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
  /*! The point to fan a face inside the zone from: its preferred point if that is
  eligible, else the smallest-id eligible point. None if no point qualifies; the face
  is then fanned from its centroid. Both pieces sharing the face choose alike. */
  [[nodiscard]] std::optional<size_t> fan_vertex(const PieceFace & face, const std::vector<size_t> & ids) const {
    // A grid face's preferred point, which both cells on it know, comes first: a cell's
    // lowest corner is then the fan vertex of all its faces, so it can be the apex
    if (face.preferred) {
      const auto & q = *face.preferred;
      for (size_t m = 0; m < face.points.size(); ++m)
        if (name_of(face.points[m]) == name_of(q) && fan_eligible(face, m)) return m;
    }
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
    if (!piece.changed) {
      for (const auto & t: piece.tets) tetrahedra_.push_back({named(t[0]), named(t[1]), named(t[2]), named(t[3])});
      return;
    }
    std::vector<size_t> corner_ids;
    for (const auto & p: piece.corners) corner_ids.push_back(named(p));
    std::vector<std::vector<size_t>> face_ids;
    std::vector<std::vector<std::array<size_t, 3>>> grid_triangles;   // sorted, for the faces that keep the grid's
    for (const auto & face: piece.faces) {
      std::vector<size_t> ids;
      for (const auto & p: face.points) ids.push_back(named(p));
      face_ids.push_back(std::move(ids));
      std::vector<std::array<size_t, 3>> tris;
      for (const auto & t: face.grid_triangles) {
        std::array<size_t, 3> tri{named(t[0]), named(t[1]), named(t[2])};
        std::sort(tri.begin(), tri.end());
        tris.push_back(tri);
      }
      std::sort(tris.begin(), tris.end());
      grid_triangles.push_back(std::move(tris));
    }
    std::vector<std::optional<size_t>> fan_from;
    for (size_t f = 0; f < piece.faces.size(); ++f)
      fan_from.push_back(!grid_triangles[f].empty() ? std::nullopt
                         : piece.faces[f].boundary ? boundary_fan_vertex(piece.faces[f], face_ids[f]) : fan_vertex(piece.faces[f], face_ids[f]));
    // whether fanning face f from corner id c gives its triangles
    auto fans_from = [&](const size_t f, const size_t c) {
      const auto & ids = face_ids[f];
      if (grid_triangles[f].empty()) return fan_from[f] && ids[*fan_from[f]] == c;
      const auto m = std::find(ids.begin(), ids.end(), c);
      if (m == ids.end()) return false;
      std::vector<std::array<size_t, 3>> tris;
      vertex_fan(piece.faces[f], ids, static_cast<size_t>(m - ids.begin()), tris);
      for (auto & t: tris) std::sort(t.begin(), t.end());
      std::sort(tris.begin(), tris.end());
      return tris == grid_triangles[f];
    };
    // the apex: the smallest-id corner that qualifies
    std::optional<size_t> apex;
    for (size_t c = 0; c < corner_ids.size(); ++c) {
      if (apex && corner_ids[c] >= corner_ids[*apex]) continue;
      bool ok{true};
      for (size_t f = 0; f < piece.faces.size() && ok; ++f)
        if (piece.incident[c].count(piece.faces[f].plane)) ok = fans_from(f, corner_ids[c]);
      if (ok) apex = c;
    }

    std::vector<std::array<size_t, 3>> triangles;
    for (size_t f = 0; f < piece.faces.size(); ++f) {
      const auto & face = piece.faces[f];
      if (apex && piece.incident[*apex].count(face.plane)) continue;
      const auto & ids = face_ids[f];
      if (!grid_triangles[f].empty()) { triangles.insert(triangles.end(), grid_triangles[f].begin(), grid_triangles[f].end()); continue; }
      if (!fan_from[f]) { fan(ids, triangles); continue; }
      vertex_fan(face, ids, *fan_from[f], triangles);
    }
    // Without a corner to cone from, the centre of a face kept whole may do: it is on
    // that face alone, whose grid triangles are its fan, and is a vertex already
    std::optional<size_t> centre_face;
    if (!apex)
      for (size_t f = 0; f < piece.faces.size() && !centre_face; ++f)
        if (piece.faces[f].centre) centre_face = f;
    if (centre_face) {
      triangles.clear();
      for (size_t f = 0; f < piece.faces.size(); ++f) {
        if (f == *centre_face) continue;
        const auto & face = piece.faces[f];
        const auto & ids = face_ids[f];
        if (!grid_triangles[f].empty()) { triangles.insert(triangles.end(), grid_triangles[f].begin(), grid_triangles[f].end()); continue; }
        if (!fan_from[f]) { fan(ids, triangles); continue; }
        vertex_fan(face, ids, *fan_from[f], triangles);
      }
    }
    const size_t c = apex ? corner_ids[*apex] : centre_face ? named(*piece.faces[*centre_face].centre) : centroid(corner_ids);
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
  static std::uint64_t edge_key(const edge_t & e) { return (static_cast<std::uint64_t>(e.first) << 32) | static_cast<std::uint64_t>(e.second); }
  [[nodiscard]] bool has_edge(const edge_t & e) const { return edge_tets_.count(edge_key(e)) > 0; }
  //! The tetrahedra around an edge, sorted
  [[nodiscard]] std::vector<tet_t> tets_around(const edge_t & e) const {
    std::vector<tet_t> out;
    if (const auto it = edge_tets_.find(edge_key(e)); it != edge_tets_.end()) for (const auto i: it->second) out.push_back(slots_[i]);
    std::sort(out.begin(), out.end());
    return out;
  }
  //! The slot of tetrahedron t (sorted), if it is in the mesh
  [[nodiscard]] std::optional<std::uint32_t> slot_of(const tet_t & t) const {
    const auto it = edge_tets_.find(edge_key(edge(t[0], t[1])));
    if (it == edge_tets_.end()) return std::nullopt;
    for (const auto i: it->second) if (slots_[i] == t) return i;
    return std::nullopt;
  }
  [[nodiscard]] bool has_tet(const tet_t & t) const { return slot_of(t).has_value(); }
  void add(tet_t t) {
    std::sort(t.begin(), t.end());
    std::uint32_t i;
    if (free_slots_.empty()) { i = static_cast<std::uint32_t>(slots_.size()); slots_.push_back(t); }
    else { i = free_slots_.back(); free_slots_.pop_back(); slots_[i] = t; }
    for (const auto & e: edges(t)) edge_tets_[edge_key(e)].push_back(i);
  }
  void remove(const tet_t & t) {
    const auto i = slot_of(t);
    if (!i) return;
    for (const auto & e: edges(t)) {
      auto it = edge_tets_.find(edge_key(e));
      auto & v = it->second;
      v.erase(std::find(v.begin(), v.end(), *i));
      if (v.empty()) edge_tets_.erase(it);
    }
    free_slots_.push_back(*i);
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
    const auto it = edge_tets_.find(edge_key(edge(f[0], f[1])));
    if (it == edge_tets_.end()) return false;
    int n{0};
    for (const auto i: it->second) if (std::count(slots_[i].begin(), slots_[i].end(), f[2])) ++n;
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
            if (has_edge(y) && boundary_edge(y)) seen.insert(y);
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
    if (!has_edge(e)) return;
    const bool on_boundary = boundary_edge(e);
    // e must be the refinement edge of its boundary triangles and of every tetrahedron
    // around it. Splitting another edge can change either, so recheck both each time.
    while (has_edge(e)) {
      std::optional<edge_t> other;
      if (on_boundary)
        for (const auto & t: tets_around(e)) {
          for (const auto & f: boundary_triangles(t))
            if (std::count(f.begin(), f.end(), e.first) && std::count(f.begin(), f.end(), e.second)) {
              const auto c = chosen_2d(f);
              if (c != e) { other = c; break; }
            }
          if (other) break;
        }
      if (!other)
        for (const auto & t: tets_around(e)) { const auto c = chosen(t); if (c != e) { other = c; break; } }
      if (!other) break;
      split_edge(*other, depth + 1);
    }
    if (!has_edge(e)) return;
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
    on_planes_.push_back(common_planes({e.first, e.second}));
    index_position(m);
    const auto around = tets_around(e);
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
