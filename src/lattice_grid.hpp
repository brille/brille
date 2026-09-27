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
#ifndef BRILLE_LATTICE_GRID_HPP_
#define BRILLE_LATTICE_GRID_HPP_
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

namespace brille::latticetri {
using int3 = std::array<long long, 3>;
using mat3i = std::array<long long, 9>;   // row-major, acts on coordinate columns
using vec3 = std::array<double, 3>;

/*! Grid points are integer coordinates in the grid lattice's basis, multiplied by
`scale`, so that the centroids used to split degenerate Delaunay cells (of up to
14 vertices) are exact. */
constexpr long long scale{840};
using point = int3;
using tetrahedron = std::array<point, 4>;

//! floor(a / b) for b > 0
constexpr long long floor_div(const long long a, const long long b) {
  return a >= 0 ? a / b : -((-a + b - 1) / b);
}

/*! \brief A periodic triangulation of a lattice, invariant under its symmetry

The grid for the structured mesh of the irreducible Brillouin zone (see the design
note). Its vertices are the points of the grid lattice L, and its tetrahedra are
L's Delaunay triangulation: the Kuhn tetrahedra {0, v_a, v_a+v_b, v_a+v_b+v_c} of
an obtuse superbase from Selling reduction. The Delaunay triangulation is unique,
so invariant under L's symmetries, unless a Selling parameter is zero (e.g.
hexagonal and monoclinic lattices). Then the Delaunay cells are larger polytopes,
which are split symmetrically: each face from its centroid, the cell from its own.

Point location takes constant time: scale x to lattice coordinates, take floor,
and test the tetrahedra of the few cells that can contain it.
*/
class Grid {
  std::array<double, 9> basis_;        // rows: grid lattice vectors (Cartesian)
  std::array<double, 9> inverse_;      // coordinates y = x · inverse_
  std::vector<mat3i> ops_;
  std::array<int3, 4> superbase_{};    // integer coordinates, sum zero
  std::array<double, 6> selling_{};    // -v_i·v_j, i < j
  bool degenerate_{false};
  std::vector<tetrahedron> pattern_;   // tetrahedra of the cells anchored in [0,1)^3
  std::vector<std::vector<size_t>> cells_;   // the pattern tetrahedra of each cell (convex)
  int3 lo_{}, hi_{};                   // extent of the pattern, in whole cells
  double size_{0};

public:
  Grid(const std::array<double, 9> & basis_rows, std::vector<mat3i> ops, const double tie_tolerance = 1e-9)
      : basis_(basis_rows), ops_(std::move(ops)) {
    invert();
    reduce(tie_tolerance);
    if (degenerate_) split_degenerate(tie_tolerance); else kuhn();
    extent();
  }
  [[nodiscard]] bool degenerate() const { return degenerate_; }
  [[nodiscard]] const std::vector<tetrahedron> & pattern() const { return pattern_; }
  //! The obtuse superbase (integer coordinates, summing to zero) the grid is built from
  [[nodiscard]] const std::array<int3, 4> & superbase() const { return superbase_; }
  //! The cells of the pattern: each a convex polytope, split into the listed pattern tetrahedra
  [[nodiscard]] const std::vector<std::vector<size_t>> & cells() const { return cells_; }
  [[nodiscard]] const std::array<double, 6> & selling_parameters() const { return selling_; }
  [[nodiscard]] const std::vector<mat3i> & operations() const { return ops_; }

  [[nodiscard]] vec3 cartesian(const point & p) const {
    vec3 x{0, 0, 0};
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) x[k] += static_cast<double>(p[i]) / scale * basis_[3 * i + k];
    return x;
  }
  [[nodiscard]] vec3 coordinates(const vec3 & x) const {
    vec3 y{0, 0, 0};
    for (int k = 0; k < 3; ++k) for (int i = 0; i < 3; ++i) y[i] += x[k] * inverse_[3 * k + i];
    return y;
  }
  [[nodiscard]] static point act(const mat3i & o, const point & p) {
    point q{0, 0, 0};
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) q[i] += o[3 * i + k] * p[k];
    return q;
  }
  [[nodiscard]] static point shifted(const point & p, const int3 & cell) {
    return {p[0] + cell[0] * scale, p[1] + cell[1] * scale, p[2] + cell[2] * scale};
  }

  //! All tetrahedra of the grid whose centroid is within `radius` of the origin
  [[nodiscard]] std::vector<tetrahedron> patch(const double radius) const {
    std::vector<tetrahedron> out;
    const auto span = static_cast<long long>(std::ceil(radius * max_inverse_row())) + 2;
    for (long long a = -span - hi_[0]; a <= span - lo_[0]; ++a)
      for (long long b = -span - hi_[1]; b <= span - lo_[1]; ++b)
        for (long long c = -span - hi_[2]; c <= span - lo_[2]; ++c)
          for (const auto & t: pattern_) {
            tetrahedron s{};
            for (int k = 0; k < 4; ++k) s[k] = shifted(t[k], {a, b, c});
            const auto x = centroid_cartesian(s);
            if (x[0] * x[0] + x[1] * x[1] + x[2] * x[2] <= radius * radius) out.push_back(s);
          }
    return out;
  }

  /*! \brief The grid tetrahedron containing x (Cartesian) and x's barycentric weights

  Weights may be as small as -tolerance, for points on a face.
  \return false if no tetrahedron contains x, which can't happen for finite x
  */
  bool locate(const vec3 & x, tetrahedron & found, std::array<double, 4> & weights, const double tolerance = 1e-12) const {
    const auto y = coordinates(x);
    int3 f{};
    for (int i = 0; i < 3; ++i) f[i] = static_cast<long long>(std::floor(y[i]));
    for (long long a = f[0] - hi_[0]; a <= f[0] - lo_[0]; ++a)
      for (long long b = f[1] - hi_[1]; b <= f[1] - lo_[1]; ++b)
        for (long long c = f[2] - hi_[2]; c <= f[2] - lo_[2]; ++c)
          for (const auto & t: pattern_) {
            tetrahedron s{};
            for (int k = 0; k < 4; ++k) s[k] = shifted(t[k], {a, b, c});
            if (barycentric(s, x, weights) && *std::min_element(weights.begin(), weights.end()) >= -tolerance) {
              found = s;
              return true;
            }
          }
    return false;
  }

  //! Whether every tetrahedron within `radius` maps onto a grid tetrahedron under every operation
  [[nodiscard]] bool invariant(const double radius) const {
    auto inner = patch(radius);
    auto outer = patch(radius + 2 * size_);
    std::set<std::array<point, 4>> all;
    for (auto t: outer) { std::sort(t.begin(), t.end()); all.insert(t); }
    for (const auto & o: ops_)
      for (const auto & t: inner) {
        tetrahedron s{};
        for (int k = 0; k < 4; ++k) s[k] = act(o, t[k]);
        std::sort(s.begin(), s.end());
        if (!all.count(s)) return false;
      }
    return true;
  }

private:
  void invert() {
    const auto & m = basis_;
    const double det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
    if (det == 0) throw std::invalid_argument("the grid lattice basis is singular");
    inverse_ = {(m[4] * m[8] - m[5] * m[7]) / det, (m[2] * m[7] - m[1] * m[8]) / det, (m[1] * m[5] - m[2] * m[4]) / det,
                (m[5] * m[6] - m[3] * m[8]) / det, (m[0] * m[8] - m[2] * m[6]) / det, (m[2] * m[3] - m[0] * m[5]) / det,
                (m[3] * m[7] - m[4] * m[6]) / det, (m[1] * m[6] - m[0] * m[7]) / det, (m[0] * m[4] - m[1] * m[3]) / det};
    for (int i = 0; i < 3; ++i) size_ = std::max(size_, std::sqrt(m[3 * i] * m[3 * i] + m[3 * i + 1] * m[3 * i + 1] + m[3 * i + 2] * m[3 * i + 2]));
  }
  [[nodiscard]] double max_inverse_row() const {
    // |y_i| <= |x| * |column i of inverse_|, a bound on the cells a ball can reach
    double out{0};
    for (int i = 0; i < 3; ++i) out = std::max(out, std::sqrt(inverse_[i] * inverse_[i] + inverse_[3 + i] * inverse_[3 + i] + inverse_[6 + i] * inverse_[6 + i]));
    return out;
  }
  [[nodiscard]] vec3 cart(const int3 & n) const {
    vec3 x{0, 0, 0};
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) x[k] += static_cast<double>(n[i]) * basis_[3 * i + k];
    return x;
  }
  static double dot(const vec3 & a, const vec3 & b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

  //! Selling reduction to an obtuse superbase (v_i·v_j <= 0)
  void reduce(const double tolerance) {
    std::array<int3, 4> n{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, -1, -1}}};
    auto v = [&](int i) { return cart(n[i]); };
    double top{0};
    for (int i = 0; i < 4; ++i) top = std::max(top, dot(v(i), v(i)));
    for (int iteration = 0; iteration < 10000; ++iteration) {
      bool changed{false};
      for (int i = 0; i < 4 && !changed; ++i)
        for (int j = i + 1; j < 4 && !changed; ++j)
          if (dot(v(i), v(j)) > 1e-13 * top) {
            int k{-1}, l{-1};
            for (int m = 0; m < 4; ++m) if (m != i && m != j) (k < 0 ? k : l) = m;
            for (int c = 0; c < 3; ++c) { n[k][c] += n[i][c]; n[l][c] += n[i][c]; n[i][c] = -n[i][c]; }
            changed = true;
          }
      if (!changed) {
        superbase_ = n;
        int p{0};
        double biggest{0};
        for (int i = 0; i < 4; ++i) biggest = std::max(biggest, dot(v(i), v(i)));
        for (int i = 0; i < 4; ++i) for (int j = i + 1; j < 4; ++j) selling_[p++] = -dot(v(i), v(j));
        degenerate_ = *std::min_element(selling_.begin(), selling_.end()) <= tolerance * biggest;
        return;
      }
    }
    throw std::runtime_error("Selling reduction did not converge");
  }

  [[nodiscard]] std::vector<tetrahedron> kuhn_cell() const {
    std::vector<tetrahedron> out;
    std::array<int, 3> perm{0, 1, 2};
    do {
      const auto & a = superbase_[perm[0]];
      const auto & b = superbase_[perm[1]];
      const auto & c = superbase_[perm[2]];
      tetrahedron t;
      t[0] = {0, 0, 0};
      for (int i = 0; i < 3; ++i) {
        t[1][i] = scale * a[i];
        t[2][i] = scale * (a[i] + b[i]);
        t[3][i] = scale * (a[i] + b[i] + c[i]);
      }
      out.push_back(t);
    } while (std::next_permutation(perm.begin(), perm.end()));
    return out;
  }
  void kuhn() {
    pattern_ = kuhn_cell();
    cells_ = {{0, 1, 2, 3, 4, 5}};
  }

  [[nodiscard]] vec3 centroid_cartesian(const tetrahedron & t) const {
    vec3 x{0, 0, 0};
    for (const auto & p: t) { const auto c = cartesian(p); for (int k = 0; k < 3; ++k) x[k] += c[k] / 4; }
    return x;
  }
  [[nodiscard]] vec3 circumcentre(const tetrahedron & t) const {
    std::array<vec3, 4> p;
    for (int k = 0; k < 4; ++k) p[k] = cartesian(t[k]);
    std::array<double, 9> a{};
    vec3 rhs{};
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) a[3 * r + c] = 2 * (p[r + 1][c] - p[0][c]);
      rhs[r] = dot(p[r + 1], p[r + 1]) - dot(p[0], p[0]);
    }
    const double det = a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6]) + a[2] * (a[3] * a[7] - a[4] * a[6]);
    vec3 x{};
    for (int c = 0; c < 3; ++c) {
      auto m = a;
      for (int r = 0; r < 3; ++r) m[3 * r + c] = rhs[r];
      x[c] = (m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6])) / det;
    }
    return x;
  }
  [[nodiscard]] static point centroid_exact(const std::vector<point> & pts) {
    point s{0, 0, 0};
    for (const auto & p: pts) for (int i = 0; i < 3; ++i) s[i] += p[i];
    const auto m = static_cast<long long>(pts.size());
    for (int i = 0; i < 3; ++i) {
      if (s[i] % m != 0) throw std::runtime_error("a Delaunay cell centroid is not representable on the grid");
      s[i] /= m;
    }
    return s;
  }

  //! Merge the Kuhn tetrahedra into Delaunay cells and split each symmetrically
  void split_degenerate(const double) {
    // Kuhn tetrahedra of a block of cells, grouped by circumcentre
    std::map<std::array<long long, 3>, std::vector<tetrahedron>> cells;
    const auto base = kuhn_cell();
    const double key = size_ * 1e-7;
    for (long long a = -3; a <= 3; ++a)
      for (long long b = -3; b <= 3; ++b)
        for (long long c = -3; c <= 3; ++c)
          for (const auto & t: base) {
            tetrahedron s{};
            for (int k = 0; k < 4; ++k) s[k] = shifted(t[k], {a, b, c});
            const auto x = circumcentre(s);
            cells[{std::llround(x[0] / key), std::llround(x[1] / key), std::llround(x[2] / key)}].push_back(s);
          }
    pattern_.clear();
    cells_.clear();
    for (const auto & [k, members]: cells) {
      // keep the cells anchored in [0,1)^3, by their exact vertex centroid
      std::set<point> vs;
      for (const auto & t: members) vs.insert(t.begin(), t.end());
      point sum{0, 0, 0};
      for (const auto & p: vs) for (int i = 0; i < 3; ++i) sum[i] += p[i];
      const auto m = static_cast<long long>(vs.size());
      bool anchored{true};
      for (int i = 0; i < 3; ++i) anchored &= floor_div(sum[i], m * scale) == 0;
      if (!anchored) continue;
      const size_t first = pattern_.size();
      if (members.size() == 1) pattern_.push_back(members[0]);
      else split_cell(members);
      std::vector<size_t> cell(pattern_.size() - first);
      std::iota(cell.begin(), cell.end(), first);
      cells_.push_back(std::move(cell));
    }
  }

  void split_cell(const std::vector<tetrahedron> & members) {
    // faces of the cell: member triangles that occur once
    std::map<std::array<point, 3>, int> count;
    std::vector<point> vertices;
    for (const auto & t: members) {
      for (int skip = 0; skip < 4; ++skip) {
        std::array<point, 3> f;
        int j{0};
        for (int k = 0; k < 4; ++k) if (k != skip) f[j++] = t[k];
        std::sort(f.begin(), f.end());
        count[f] += 1;
      }
      vertices.insert(vertices.end(), t.begin(), t.end());
    }
    std::sort(vertices.begin(), vertices.end());
    vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
    const auto c = centroid_exact(vertices);
    // merge coplanar boundary triangles into polygons
    struct polygon { vec3 n; double d; std::vector<std::array<point, 3>> triangles; };
    std::vector<polygon> polys;
    for (const auto & [f, k]: count) {
      if (k != 1) continue;
      const auto a = cartesian(f[0]), b = cartesian(f[1]), cc = cartesian(f[2]);
      vec3 u{b[0] - a[0], b[1] - a[1], b[2] - a[2]}, v{cc[0] - a[0], cc[1] - a[1], cc[2] - a[2]};
      vec3 n{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
      const double len = std::sqrt(dot(n, n));
      for (auto & x: n) x /= len;
      const double d = dot(n, a);
      bool merged{false};
      for (auto & p: polys)
        if (std::abs(std::abs(dot(p.n, n)) - 1) < 1e-9 && std::abs(dot(p.n, a) - p.d) < 1e-9 * (1 + std::abs(d))) {
          p.triangles.push_back(f);
          merged = true;
          break;
        }
      if (!merged) polys.push_back({n, d, {f}});
    }
    for (const auto & p: polys) {
      if (p.triangles.size() == 1) {
        const auto & f = p.triangles[0];
        pattern_.push_back({f[0], f[1], f[2], c});
        continue;
      }
      std::set<point> pv;
      std::map<std::array<point, 2>, int> edges;
      for (const auto & f: p.triangles) {
        pv.insert(f.begin(), f.end());
        for (int i = 0; i < 3; ++i) for (int j = i + 1; j < 3; ++j) edges[{std::min(f[i], f[j]), std::max(f[i], f[j])}] += 1;
      }
      const auto fc = centroid_exact(std::vector<point>(pv.begin(), pv.end()));
      for (const auto & [e, k]: edges) if (k == 1) pattern_.push_back({e[0], e[1], fc, c});
    }
  }

  void extent() {
    lo_ = {0, 0, 0};
    hi_ = {0, 0, 0};
    for (const auto & t: pattern_)
      for (const auto & p: t)
        for (int i = 0; i < 3; ++i) {
          lo_[i] = std::min(lo_[i], floor_div(p[i], scale));
          hi_[i] = std::max(hi_[i], -floor_div(-p[i], scale));   // ceiling
        }
  }

  bool barycentric(const tetrahedron & t, const vec3 & x, std::array<double, 4> & w) const {
    std::array<vec3, 4> p;
    for (int k = 0; k < 4; ++k) p[k] = cartesian(t[k]);
    auto det3 = [](const vec3 & a, const vec3 & b, const vec3 & c) {
      return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
    };
    auto sub = [](const vec3 & a, const vec3 & b) { return vec3{a[0] - b[0], a[1] - b[1], a[2] - b[2]}; };
    const double vol = det3(sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]));
    if (vol == 0) return false;
    w[1] = det3(sub(x, p[0]), sub(p[2], p[0]), sub(p[3], p[0])) / vol;
    w[2] = det3(sub(p[1], p[0]), sub(x, p[0]), sub(p[3], p[0])) / vol;
    w[3] = det3(sub(p[1], p[0]), sub(p[2], p[0]), sub(x, p[0])) / vol;
    w[0] = 1 - w[1] - w[2] - w[3];
    return true;
  }
};
}
#endif
