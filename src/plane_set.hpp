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
#ifndef BRILLE_PLANE_SET_HPP_
#define BRILLE_PLANE_SET_HPP_
#include <algorithm>
#include <array>
#include <vector>
#include "array_l_.hpp"
#include "approx_float.hpp"
#include "math.hpp"
#include "utilities.hpp"

namespace brille {
/*! \brief Planes, each through three points in lattice coordinates, as plain numbers

`point_inside_all_planes(a, b, c, x)` on lattice vectors builds lattice-aware
temporaries (views, differences, cross and dot products) for every plane and
point. Profiling showed that reference counting and allocation, not arithmetic,
took most of its time, and that sharing those counters stopped it from scaling
with threads. PlaneSet holds the same planes as plain numbers and repeats the
same arithmetic, in the same order and with the same helpers, so its results
are identical, without temporaries.
*/
class PlaneSet {
public:
  using v3 = std::array<double, 3>;
private:
  std::vector<v3> a_, b_, c_;
  std::array<double, 9> own_metric_{};   // metric of the points' LengthUnit, used by dot
  std::array<double, 9> other_metric_{}; // metric of the other LengthUnit, used by star()
  double cross_scale_{0};                // lattice volume / 2π, as in cross(LVec, LVec)
  double tol_{0};
  int n_{1};
public:
  //! Planes through (a[j], b[j], c[j]); tolerances as for point_inside_all_planes
  template<class T, template<class> class L>
  PlaneSet(const L<T>& a, const L<T>& b, const L<T>& c, const double tol, const int n): tol_(tol), n_(n) {
    const auto& lat = a.lattice();
    const auto type = a.type();
    const auto other = type == LengthUnit::angstrom ? LengthUnit::inverse_angstrom : LengthUnit::angstrom;
    own_metric_ = lat.metric(type);
    other_metric_ = lat.metric(other);
    cross_scale_ = static_cast<double>(lat.volume(type) / math::two_pi);
    for (ind_t j=0; j<a.size(0); ++j){
      a_.push_back({static_cast<double>(a.val(j,0)), static_cast<double>(a.val(j,1)), static_cast<double>(a.val(j,2))});
      b_.push_back({static_cast<double>(b.val(j,0)), static_cast<double>(b.val(j,1)), static_cast<double>(b.val(j,2))});
      c_.push_back({static_cast<double>(c.val(j,0)), static_cast<double>(c.val(j,1)), static_cast<double>(c.val(j,2))});
    }
  }
  [[nodiscard]] size_t size() const { return a_.size(); }
  //! pseudo_orient3d(a, b, c, q) for plane j: dot(a-q, cross(b-q, c-q)) in lattice coordinates
  [[nodiscard]] double orient(const size_t j, const v3& q) const {
    v3 u, v, w, cr, s, mu;
    for (int k=0; k<3; ++k){ u[k] = a_[j][k] - q[k]; v[k] = b_[j][k] - q[k]; w[k] = c_[j][k] - q[k]; }
    utils::vector_cross<double,double,double,3>(cr.data(), v.data(), w.data());
    for (auto& x: cr) x *= cross_scale_;                                          // cross_star *= V/2π
    utils::multiply_matrix_vector(s.data(), other_metric_.data(), cr.data());      // .star()
    for (auto& x: s) x /= math::two_pi;
    utils::mul_mat_vec(mu.data(), 3u, own_metric_.data(), u.data());              // same_lattice_dot
    double out{0};
    for (int k=0; k<3; ++k) out += mu[k] * s[k];
    return out;
  }
  //! point_inside_all_planes(a, b, c, q, tol, n)
  [[nodiscard]] bool inside(const v3& q) const {
    double v = orient(0, q);
    for (size_t j=1; j<a_.size(); ++j) v = std::min(v, orient(j, q));
    return v > 0 || approx_float::scalar(v, 0., tol_, tol_, n_);
  }
};
}
#endif
