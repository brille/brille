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
#ifndef BRILLE_EXACT_GEOMETRY_HPP_
#define BRILLE_EXACT_GEOMETRY_HPP_
/*! \file
\brief Exact planes and points for the structured mesh's zone boundary

Every plane the boundary construction meets is one of two kinds, named by
integers (coordinates are in the primitive reciprocal lattice basis):

- a metric plane  x·Gσ = uᵀGσ  with σ ∈ ℤ³ and u ∈ ½ℤ³ (stored as 2u): the zone
  planes x·Gτ = τᵀGτ/2 and their images;
- an integer plane  c·x = k  with c ∈ ℤ³, k ∈ ℤ: the wedge planes (integer
  Dirichlet-cone normals) and their images.

A point is the intersection of three planes. Mapping by x → g x + t (g an
integer unimodular matrix of the point group, t ∈ ℤ³) takes names to names with
integer arithmetic only, because G is invariant under g: σ → gσ, 2u → g(2u) + 2t;
c → g⁻ᵀc, k → k + (g⁻ᵀc)·t.

Predicates (which side of a plane a point is on) are signs of polynomials in the
entries of G with integer coefficients. They are evaluated exactly with floating
point expansions (Shewchuk), so there are no tolerances, and no big integers.
*/
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace brille::exact {

// --- expansion arithmetic -----------------------------------------------------
/*! A number held exactly as the sum of doubles, in increasing magnitude with no
overlap (Shewchuk's expansions). Sums and products of doubles and expansions are
exact; the sign is the sign of the largest component. */
class Expansion {
  std::vector<double> c_;
  static void two_sum(const double a, const double b, double & x, double & y) {
    x = a + b;
    const double bv = x - a;
    y = (a - (x - bv)) + (b - bv);
  }
  static void two_product(const double a, const double b, double & x, double & y) {
    x = a * b;
    y = std::fma(a, b, -x);
  }
  //! add the double b (Shewchuk's grow-expansion); keeps the non-overlapping order
  void grow(const double b) {
    std::vector<double> out;
    out.reserve(c_.size() + 1);
    double q = b;
    for (const double e: c_) {
      double s, h;
      two_sum(q, e, s, h);
      if (h != 0) out.push_back(h);
      q = s;
    }
    if (q != 0 || out.empty()) out.push_back(q);
    c_.swap(out);
  }
public:
  Expansion() = default;
  explicit Expansion(const double x) { if (x != 0) c_.push_back(x); }
  [[nodiscard]] int sign() const {
    for (auto it = c_.rbegin(); it != c_.rend(); ++it) if (*it != 0) return *it > 0 ? 1 : -1;
    return 0;
  }
  [[nodiscard]] double estimate() const { double s{0}; for (const double x: c_) s += x; return s; }
  Expansion & operator+=(const Expansion & o) { for (const double x: o.c_) grow(x); return *this; }
  Expansion & operator-=(const Expansion & o) { for (const double x: o.c_) grow(-x); return *this; }
  friend Expansion operator+(Expansion a, const Expansion & b) { a += b; return a; }
  friend Expansion operator-(Expansion a, const Expansion & b) { a -= b; return a; }
  friend Expansion operator-(Expansion a) { for (auto & x: a.c_) x = -x; return a; }
  //! multiply by a double exactly
  [[nodiscard]] Expansion scaled(const double b) const {
    Expansion out;
    for (const double x: c_) {
      double p, e;
      two_product(x, b, p, e);
      out.grow(e);
      out.grow(p);
    }
    return out;
  }
  friend Expansion operator*(const Expansion & a, const Expansion & b) {
    Expansion out;
    for (const double x: b.c_) out += a.scaled(x);
    return out;
  }
};

using int3 = std::array<long long, 3>;
using mat3i = std::array<long long, 9>;   // row-major, acts on coordinate columns

//! the metric G (reciprocal lattice coordinates), with exact matrix-vector products
class Metric {
  std::array<double, 9> g_;
public:
  explicit Metric(const std::array<double, 9> & g) : g_(g) {}
  [[nodiscard]] double operator()(const int i, const int j) const { return g_[3 * i + j]; }
  [[nodiscard]] const std::array<double, 9> & values() const { return g_; }
  //! (G v)_i for integer v, exactly
  [[nodiscard]] Expansion row_dot(const int i, const int3 & v) const {
    Expansion out;
    for (int j = 0; j < 3; ++j) if (v[j]) out += Expansion(g_[3 * i + j]).scaled(static_cast<double>(v[j]));
    return out;
  }
  //! aᵀ G b for integer a, b, exactly
  [[nodiscard]] Expansion form(const int3 & a, const int3 & b) const {
    Expansion out;
    for (int i = 0; i < 3; ++i) if (a[i]) out += row_dot(i, b).scaled(static_cast<double>(a[i]));
    return out;
  }
  //! whether gᵀ G g == G exactly, for the point group operation g
  [[nodiscard]] bool invariant_under(const mat3i & g) const {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        const int3 a{g[i], g[3 + i], g[6 + i]};   // column i of g
        const int3 b{g[j], g[3 + j], g[6 + j]};
        if ((form(a, b) - Expansion(g_[3 * i + j])).sign() != 0) return false;
      }
    return true;
  }
};

// --- planes ---------------------------------------------------------------------
/*! A named plane n·x = d, with n and d exact. Oriented: the half-space n·x <= d is
"inside" when a plane bounds a region. */
struct Plane {
  enum class Kind { metric, integer };
  Kind kind{Kind::integer};
  int3 a{0, 0, 0};      // metric: σ; integer: c
  int3 b{0, 0, 0};      // metric: 2u; integer: (k, 0, 0)
  int orientation{1};   // multiplies n and d: -1 flips the inside

  static Plane metric_plane(const int3 & sigma, const int3 & twice_u) { return {Kind::metric, sigma, twice_u, 1}; }
  static Plane integer_plane(const int3 & c, const long long k) { return {Kind::integer, c, {k, 0, 0}, 1}; }
  [[nodiscard]] Plane flipped() const { Plane p = *this; p.orientation = -p.orientation; return p; }

  //! the coefficients (n, d), exactly
  [[nodiscard]] std::array<Expansion, 4> coefficients(const Metric & G) const {
    std::array<Expansion, 4> out;
    const auto s = static_cast<double>(orientation);
    if (kind == Kind::metric) {
      for (int i = 0; i < 3; ++i) out[i] = G.row_dot(i, a).scaled(s);
      out[3] = G.form(b, a).scaled(0.5 * s);
    } else {
      for (int i = 0; i < 3; ++i) out[i] = Expansion(static_cast<double>(a[i]) * s);
      out[3] = Expansion(static_cast<double>(b[0]) * s);
    }
    return out;
  }
  //! the image under x -> g x + t; g must leave G invariant
  [[nodiscard]] Plane mapped(const mat3i & g, const int3 & t) const {
    Plane p = *this;
    if (kind == Kind::metric) {
      p.a = apply(g, a);
      const auto gb = apply(g, b);
      for (int i = 0; i < 3; ++i) p.b[i] = gb[i] + 2 * t[i];
    } else {
      p.a = apply(inverse_transpose(g), a);
      p.b = {b[0] + p.a[0] * t[0] + p.a[1] * t[1] + p.a[2] * t[2], 0, 0};
    }
    return p;
  }
  static int3 apply(const mat3i & g, const int3 & v) {
    return {g[0] * v[0] + g[1] * v[1] + g[2] * v[2], g[3] * v[0] + g[4] * v[1] + g[5] * v[2], g[6] * v[0] + g[7] * v[1] + g[8] * v[2]};
  }
  //! g⁻ᵀ for a unimodular integer matrix
  static mat3i inverse_transpose(const mat3i & m) {
    const long long det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
    if (det != 1 && det != -1) throw std::invalid_argument("point group operations must be unimodular");
    // cofactor matrix / det is the inverse transpose
    mat3i c{m[4] * m[8] - m[5] * m[7], m[5] * m[6] - m[3] * m[8], m[3] * m[7] - m[4] * m[6],
            m[2] * m[7] - m[1] * m[8], m[0] * m[8] - m[2] * m[6], m[1] * m[6] - m[0] * m[7],
            m[1] * m[5] - m[2] * m[4], m[2] * m[3] - m[0] * m[5], m[0] * m[4] - m[1] * m[3]};
    for (auto & x: c) x *= det;
    return c;
  }
};

//! 3×3 determinant of expansions
inline Expansion det3(const std::array<std::array<Expansion, 3>, 3> & m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
       - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
       + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

// --- points ---------------------------------------------------------------------
//! The point where three planes meet
struct Point {
  std::array<Plane, 3> planes;
};

/*! \brief Exact geometry on named planes and points for one metric */
class Geometry {
  Metric G_;
public:
  explicit Geometry(const Metric & G) : G_(G) {}
  [[nodiscard]] const Metric & metric() const { return G_; }

  //! sign of det of the three planes' normals: zero if they don't meet in a point
  [[nodiscard]] int independent(const Point & p) const {
    std::array<std::array<Expansion, 3>, 3> n;
    for (int r = 0; r < 3; ++r) { const auto c = p.planes[r].coefficients(G_); for (int k = 0; k < 3; ++k) n[r][k] = c[k]; }
    return det3(n).sign();
  }
  /*! sign of n·x - d for the plane q at the point p: +1 outside, 0 on, -1 inside */
  [[nodiscard]] int side(const Point & p, const Plane & q) const {
    std::array<std::array<Expansion, 4>, 3> rows;
    for (int r = 0; r < 3; ++r) rows[r] = p.planes[r].coefficients(G_);
    const auto qc = q.coefficients(G_);
    // x = N⁻¹ d; n_q·x - d_q = (n_q · adj(N) d - d_q det N) / det N
    std::array<std::array<Expansion, 3>, 3> N;
    for (int r = 0; r < 3; ++r) for (int k = 0; k < 3; ++k) N[r][k] = rows[r][k];
    const Expansion det = det3(N);
    const int ds = det.sign();
    if (ds == 0) throw std::runtime_error("a point's three planes do not meet in a point");
    // n_q · adj(N) d = det of N with row... use Cramer: x_k det = det(N with column k replaced by d)
    Expansion numerator = -(qc[3] * det);
    for (int k = 0; k < 3; ++k) {
      auto M = N;
      for (int r = 0; r < 3; ++r) M[r][k] = rows[r][3];
      numerator += qc[k] * det3(M);
    }
    return numerator.sign() * ds;
  }
  [[nodiscard]] bool on(const Point & p, const Plane & q) const { return side(p, q) == 0; }
  //! whether two named points are the same point
  [[nodiscard]] bool same(const Point & p, const Point & q) const {
    return on(p, q.planes[0]) && on(p, q.planes[1]) && on(p, q.planes[2]);
  }
  //! whether two planes are the same plane (orientation ignored)
  [[nodiscard]] bool same(const Plane & a, const Plane & b) const {
    const auto x = a.coefficients(G_), y = b.coefficients(G_);
    // parallel normals: n_a × n_b = 0; same offset: d_a n_b = d_b n_a
    for (int i = 0; i < 3; ++i) {
      const int j = (i + 1) % 3, k = (i + 2) % 3;
      if ((x[j] * y[k] - x[k] * y[j]).sign() != 0) return false;
    }
    for (int i = 0; i < 3; ++i) if ((x[3] * y[i] - y[3] * x[i]).sign() != 0) return false;
    return true;
  }
  //! approximate coordinates, for output
  [[nodiscard]] std::array<double, 3> coordinates(const Point & p) const {
    std::array<std::array<double, 4>, 3> r;
    for (int i = 0; i < 3; ++i) { const auto c = p.planes[i].coefficients(G_); for (int k = 0; k < 4; ++k) r[i][k] = c[k].estimate(); }
    auto d3 = [](double a, double b, double c, double d, double e, double f, double g, double h, double i) {
      return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    };
    const double det = d3(r[0][0], r[0][1], r[0][2], r[1][0], r[1][1], r[1][2], r[2][0], r[2][1], r[2][2]);
    return {d3(r[0][3], r[0][1], r[0][2], r[1][3], r[1][1], r[1][2], r[2][3], r[2][1], r[2][2]) / det,
            d3(r[0][0], r[0][3], r[0][2], r[1][0], r[1][3], r[1][2], r[2][0], r[2][3], r[2][2]) / det,
            d3(r[0][0], r[0][1], r[0][3], r[1][0], r[1][1], r[1][3], r[2][0], r[2][1], r[2][3]) / det};
  }
  [[nodiscard]] Point mapped(const Point & p, const mat3i & g, const int3 & t) const {
    return {{p.planes[0].mapped(g, t), p.planes[1].mapped(g, t), p.planes[2].mapped(g, t)}};
  }
};
}
#endif
