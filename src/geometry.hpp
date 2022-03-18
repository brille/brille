/* This file is part of brille

Copyright © 2019-2022 Greg Tucker <gregory.tucker@ess.eu>

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

#ifndef BRILLE_GEOMETRY_H_
#define BRILLE_GEOMETRY_H_
/*! \file
    \author Greg Tucker
    \brief A class representing convex polyhedra in three dimensions
*/
#include <random>
#include <utility>
#include <csignal>
#include "tetgen.h"
#include "array_.hpp" // defines bArray
#include "array_l_.hpp"
#include "utilities.hpp"
namespace brille {


    /*! \brief Identify the relative orientation of three points
     *
     * Using the same principle as Schewchuk's `orient2dfast`, this algorithm is
     * able to determine whether a test point is on the left or right of a line
     * defined by two points, when viewed from above.
     *
     * \param a The first point of the line(s)
     * \param b The second point of the line(s)
     * \param c The test point(s)
     *
     * \note This function broadcasts singleton lines or singleton points to
     *       match the input points or lines, respectively. If multiple lines
     *       and multiple points are provided no broadcasting is performed and
     *       their multiplicities must match.
    */
    template<class T, template<class> class A>
    std::enable_if_t<isBareArray<T,A>, std::vector<T>>
    pseudo_orient2d(const A<T>& a, const A<T>& b, const A<T>& c){
      auto predicate = [](const T& ax, const T& ay, const T& bx, const T& by, const T& cx, const T& cy){
        // possibly replace this with the floating-point predicate?
        return (ax - cx) * (by - cy) - (ay - cy) * (bx - cx);
      };
      assert(a.ndim() == 2u && b.ndim() == 2u && c.ndim() == 2u);
      assert(a.size(1u) == 2u && b.size(1u) == 2u && c.size(1u) == 2u);
      assert(b.size(0u) == a.size(0u));
      // the BroadcastIt2 will raise an error if a & c are incompatible
      std::array<ind_t,2> ash{{a.size(0u), 1u}}, csh{{c.size(0u), 1u}};
      auto broadcaster = BroadcastIt2(ash, csh);
      std::vector<T> out;
      out.reserve(broadcaster.size());
      for (auto [outerSub, aSub, cSub]: broadcaster){
        const auto & ax{a[aSub]};
        const auto & bx{b[aSub]};
        const auto & cx{c[cSub]};
        aSub[1] = 1u;
        cSub[1] = 1u;
        out.push_back(predicate(ax, a[aSub], bx, b[aSub], cx, c[aSub]));
      }
      return out;
    }
    template<class T, template<class> class A>
    std::enable_if_t<isLatVec<T,A>, std::vector<T>>
    pseudo_orient2d(const A<T>& a, const A<T>& b, const A<T>& c){
      assert(a.same_lattice(b) && a.same_lattice(c));
      return pseudo_orient2d(a.xyz(), b.xyz(), c.xyz());
    }

    template<class T>
    T predicate3d(const T& ax, const T& ay, const T& az,
                   const T& bx, const T& by, const T& bz,
                   const T& cx, const T& cy, const T& cz,
                   const T& dx, const T& dy, const T& dz){
    // possibly replace this with the floating-point predicate?
    auto adx = ax - dx;
    auto bdx = bx - dx;
    auto cdx = cx - dx;
    auto ady = ay - dy;
    auto bdy = by - dy;
    auto cdy = cy - dy;
    auto adz = az - dz;
    auto bdz = bz - dz;
    auto cdz = cz - dz;
    return adx * (bdy * cdz - bdz * cdy) + bdx * (cdy * adz - cdz * ady) + cdx * (ady * bdz - adz * bdy);
  };

    /*! \brief Identify the relative orientation of four points
     *
     * Using the same principle as Schewchuk's `orient3dfast`, this algorithm is
     * able to determine whether a test point is above or below a plane
     * defined by three points. The returned value is positive if the test point
     * is below the plane, zero if it is on the plane, and negative if it is
     * above. The points of the plane are ordered counterclockwise when viewed
     * from above.
     *
     * \param a The first point on the plane(s)
     * \param b The second point on the plane(s)
     * \param c The third point on the plane(s)
     * \param d The test points(s)
     *
     * \note This function broadcasts singleton planes or singleton points to
     *       match the input points or planes, respectively. If multiple planes
     *       and multiple points are provided no broadcasting is performed and
     *       their multiplicities must match.
     */
    template<class T, template<class> class A>
    std::enable_if_t<isBareArray<T,A>, std::vector<T>>
    pseudo_orient3d(const A<T>& a, const A<T>& b, const A<T>& c, const A<T>& d){
      assert(a.ndim() == 2u && b.ndim() == 2u && c.ndim() == 2u && d.ndim() == 2u);
      assert(a.size(1u) == 3u && b.size(1u) == 3u && c.size(1u) == 3u && d.size(1u) == 3u);
      assert(b.size(0u) == a.size(0u) && c.size(0u) == a.size(0u));
      // the BroadcastIt2 will raise an error if a & d are incompatible
      std::array<ind_t,2> ash{{a.size(0u), 1u}}, dsh{{d.size(0u), 1u}};
      auto broadcaster = BroadcastIt2(ash, dsh);
      std::vector<T> out;
      out.reserve(broadcaster.size());
      for (auto [outerSub, aSub, dSub]: broadcaster){
//        const auto & ax{a[aSub]};
//        const auto & bx{b[aSub]};
//        const auto & cx{c[aSub]};
//        const auto & dx{d[dSub]};
//        aSub[1] = 1u;
//        dSub[1] = 1u;
//        const auto & ay{a[aSub]};
//        const auto & by{b[aSub]};
//        const auto & cy{c[aSub]};
//        const auto & dy{d[dSub]};
//        aSub[1] = 2u;
//        dSub[2] = 2u;
//        const auto & az{a[aSub]};
//        const auto & bz{b[aSub]};
//        const auto & cz{c[aSub]};
//        const auto & dz{d[dSub]};
//        auto pred = predicate3d(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);

        const auto ar{a.view(aSub[0])};
        const auto br{b.view(aSub[0])};
        const auto cr{c.view(aSub[0])};
        const auto dr{d.view(dSub[0])};
        auto tp = dot(ar-dr, cross(br-dr, cr-dr)).sum();
        out.push_back(tp);
      }
      return out;
    }
    template<class T, template<class> class A>
    std::enable_if_t<isLatVec<T,A>, std::vector<T>>
    pseudo_orient3d(const A<T>& a, const A<T>& b, const A<T>& c, const A<T>& d){
      assert(a.same_lattice(b) && a.same_lattice(c) && d.same_lattice(d));
      return pseudo_orient3d(a.xyz(), b.xyz(), c.xyz(), d.xyz());
    }

//  template<class T, template<class> class A>
//  std::enable_if_t<isBareArray<T,A>, std::tuple<A<T>, A<T>>>
//  plane_points_from_normal(const A<T> & n, const A<T> & p) {
//    A<T> a(n.shape(), T(1)), b;
//    auto n_norm = norm(n); // this *should* be all 1, but maybe it's not
//    auto nn = n / n_norm;
//    auto a_dot_n = dot(a, nn) * nn;
//    auto is_n = a_dot_n.row_is(brille::cmp::eq, nn * norm(a));
//    if (is_n.any()) {
//      for (ind_t i = 0; i < a.size(0); ++i)
//        if (is_n.val(i)) {
//          // set to (2, 0, 1) so that subtracting (1, 1, 1) gives (1, -1, 0)
//          a[{i, 0}] = T(2);
//          a[{i, 1}] = T(0);
//          for (ind_t j = 2; j < a.size(1); ++j) a[{i, j}] = T(1);
//        }
//    }
//    a -= a_dot_n;
//    a /= norm(a);
//    b = cross(nn, a);  // so that the points (p, a, b) define a plane with normal n
//    // the returned a and b should have the property (x - p) . n == 0
//    return std::make_tuple(a + p, b + p);
//  }
//
//  template<class T, template<class> class A>
//  std::enable_if_t<isLatVec<T,A>, std::tuple<A<T>, A<T>>>
//  plane_points_from_normal(const A<T> & n, const A<T> & p) {
//    assert(n.same_lattice(p));
//    auto [b, c] = plane_points_from_normal(n.xyz(), p.xyz());
//    return std::make_tuple(A<T>::from_invA(p.lattice(), b), A<T>::from_invA(p.lattice(), c));
//  }

  template<class T, template<class> class A>
  std::enable_if_t<isArray<T,A>, std::tuple<A<T>, A<T>, A<T>>>
  plane_points_from_normal(const A<T> & n, const A<T> & p) {
    verbose_update("Find planes from, normals\n", n.to_string(), "points\n", p.to_string());
    auto a = T(1) + T(0) * n;
    auto b = T(0) * n;
    auto n_norm = norm(n); // this *should* be all 1, but maybe it's not
    auto nn = n / n_norm;
    auto a_dot_n = dot(a, nn) * nn;
    auto is_n = a_dot_n.row_is(brille::cmp::eq, a);
    if (is_n.any()) {
      for (ind_t i = 0; i < a.size(0); ++i)
        if (is_n.val(i)) {
          // set to (2, 0, 1) so that subtracting (1, 1, 1) gives (1, -1, 0)
          a[{i, 0}] = T(2);
          a[{i, 1}] = T(0);
          for (ind_t j = 2; j < a.size(1); ++j) a[{i, j}] = T(1);
        }
    }
    a -= a_dot_n;
    a /= norm(a);
    b = cross(nn, a);  // so that the points (p, a, b) define a plane with normal n
    // the returned a and b should have the property (x - p) . n == 0
    return std::make_tuple(p, a + p, b + p);
  }

//  template<class T, template<class> class A>
//  std::enable_if_t<isBareArray<T,A>, std::tuple<A<T>, A<T>, A<T>>>
//  plane_points_from_normal(const A<T>& n, const A<T>& p){
//    A<T> xhat({1u, 3u}, T(0)), yhat({1u, 3u}, T(0)), zhat({1u, 3u}, T(0));
//    xhat[{0u, 0u}] = T(1);
//    yhat[{0u, 1u}] = T(1);
//    zhat[{0u, 2u}] = T(1);
//    auto a = T(0) * p;
//    auto b = T(0) * p;
//    auto c = T(0) * p;
//    for (ind_t i=0; i < n.size(0); ++i){
//      auto nx = dot(xhat, n.view(i)).sum();
//      auto ny = dot(yhat, n.view(i)).sum();
//      auto nz = dot(zhat, n.view(i)).sum();
//      auto np = dot(p.view(i), n.view(i)).sum();
//      a[{i, 0u}] = nx ? np / nx : T(1);
//      b[{i, 1u}] = ny ? np / ny : T(1);
//      c[{i, 2u}] = nz ? np / nz : T(1);
//      if (!nx && !ny && !nz) throw std::logic_error("Zero-length normal vectors not allowed");
//      if (!nx && ny && nz) a.set(i, a.view(i) + b.view(i));
//      if (nx && !ny && nz) b.set(i, b.view(i) + c.view(i));
//      if (nx && ny && !nz) c.set(i, c.view(i) + a.view(i));
//      if (!nx && !ny && nz) {
//        a.set(i, a.view(i) + c.view(i));
//        b.set(i, b.view(i) + c.view(i));
//      }
//      if (!nx && ny && !nz) {
//        a.set(i, a.view(i) + b.view(i));
//        c.set(i, c.view(i) + b.view(i));
//      }
//      if (nx && !ny && !nz) {
//        b.set(i, b.view(i) + a.view(i));
//        c.set(i, c.view(i) + a.view(i));
//      }
//      if (dot(three_point_normal(a.view(i), b.view(i), c.view(i)), n.view(i)).sum() < 0){
//        auto t = T(1) * b.view(i);
//        b.set(i, c.view(i));
//        c.set(i, t.view(0));
//      }
//    }
//    return std::make_tuple(a, b, c);
//  }

//  template<class T, template<class> class A>
//  std::enable_if_t<isLatVec<T,A>, std::tuple<A<T>, A<T>, A<T>>>
//  plane_points_from_normal(const A<T>& n, const A<T>& p){
//    auto nxyz = n.xyz();
//    auto pxyz = p.xyz();
//    auto [a, b, c] = plane_points_from_normal(nxyz, pxyz);
//    auto al = A<T>::from_invA(n.lattice(), a);
//    auto bl = A<T>::from_invA(n.lattice(), b);
//    auto cl = A<T>::from_invA(n.lattice(), c);
//    return std::make_tuple(al, bl, cl);
//  }


  template<class T, template<class> class A> std::enable_if_t<isArray<T,A>, A<T>>
  three_point_normal(const A<T> &a, const A<T> &b, const A<T> &c) {
    auto n = cross(b - a, c - b);
    return n / norm(n);
  }

  template<class T, template<class> class A, class I> std::enable_if_t<isArray<T,A>, A<T>>
  three_point_normal(const A<T> &p, const I a, const I b, const I c) {
    return three_point_normal(p.view(a), p.view(b), p.view(c));
  }

  template<class T, template<class> class A, class I> std::enable_if_t<isArray<T,A>, A<T>>
  three_point_normal(const A<T> &p, const std::vector<I> &f) {
    return three_point_normal(p.view(f[0]), p.view(f[1]), p.view(f[2]));
  }

  template<class T, class I, template<class> class A, template<class> class B>
  std::enable_if_t<isArray<T,A> && isBareArray<I,B>, A<T>>
  three_point_normal(const A<T> &p, const B<I> &t) {
    A<T> out(t.size(0), 3u);
    for (ind_t i = 0; i < t.size(0); ++i)
      out.set(i, three_point_normal(p, t[{i, 0}], t[{i, 1}], t[{i, 2}]));
    return out;
  }


//! Determine if a facet of a Polyhedron is part of a convex volume
  template<typename T, typename I>
  static bool is_not_dangling(const std::vector<I> &counts, const std::vector<T> &face) {
    return std::all_of(face.begin(), face.end(), [counts](const T &x) { return counts[x] > 2u; });
  }

  template<class T, template<class> class A> static std::enable_if_t<isArray<T,A>, void>
  on_plane_vector(const A<T>& x, const A<T>& y, const A<T>& v3, T* v2){
    v2[0] = dot(v3, x).sum();
    v2[1] = dot(v3, y).sum();
  }

//  //! Determine the winding angles for the vertices of a Polyhedron facet
//  template<class T>
//  std::vector<T> bare_winding_angles(const bArray<T> &vecs, const ind_t i, const bArray<T> &n) {
//    if (vecs.ndim() != 2 || vecs.size(1) != 3)
//      throw std::runtime_error("Finding a winding angle requires the cross product, which is only defined in 3D");
//    // vecs should be normalized already
//    std::vector<T> angles(vecs.size(0), 0.);
//    T dotij, y_len, angij;
//    bArray<T> x(1u, 3u), y(1u, 3u); // ensure all have memory allocated
//    T crsij[3]; // to hold the cross product
//    for (ind_t j = 0; j < vecs.size(0); ++j)
//      if (i != j) {
//        dotij = vecs.dot(i, j);
//        brille::utils::vector_cross(crsij, vecs.ptr(i), vecs.ptr(j));
//        x = dotij * vecs.view(i);
//        y = vecs.view(j) - x;
//        y_len = y.norm(0) * (std::signbit(brille::utils::vector_dot(crsij, n.ptr(0))) ? -1 : 1);
//        angij = std::atan2(y_len, dotij);
//        angles[j] = angij < 0 ? angij + 2 * brille::math::pi : angij;
//      }
//    return angles;
//  }

  template<class T, class R, template<class> class A, template<class> class B>
  [[nodiscard]]
  std::enable_if_t<isArray<T,A> && isArray<R,B>, std::vector<T>>
  bare_winding_angles(const A<T>& vectors, const ind_t i, const B<R>& n){
    T b[2];
    // double check that we're dealing with vectors in a plane:
    auto nn = norm(n);
    auto pv = vectors - dot(vectors, n / nn) * (n / nn);
//    verbose_update("vectors\n", vectors.to_string(), "-> pv\n", pv.to_string());
    auto x = pv.view(i) / norm(pv.view(i));
    auto y = cross(n / nn, x);
    std::vector<T> angles(vectors.size(0), 0.);
    for (ind_t j=0; j< vectors.size(0); ++j) if(i != j) {
        on_plane_vector(x, y, pv.view(j), b);
//        verbose_update("j:", j, "b: [", b[0], ", ", b[1], "]");
//        auto ang_ij = std::atan2(b[1], b[0]);
//        angles[j] = ang_ij < 0 ? ang_ij + 2 * brille::math::pi : ang_ij;
        // Use fractions of pi instead of radian
        auto ang_ij = std::atan2(b[1], b[0]) / std::atan2(0., -1.);
        angles[j] = ang_ij < 0 ? ang_ij + 2 : ang_ij;
      }
//    verbose_update("Relative to index ", i," winding angles", angles);
    return angles;
  }


  //! Check whether the points defining a Polyhedron facet encompass a finite area
  template<class T, template<class> class A>
  std::enable_if_t<isArray<T,A>, bool>
  face_has_area(const A<T> &points, const bool strict = false) {
    // first verify that all points are coplanar
    // pick the first three points to define a plane, then ensure all points are in it
    if (points.size(0) < 3) return -2; // can't be a face
    // move to the 2-D face coordinate system so that we can use orient2d
    auto centre = points.sum(0) / static_cast<T>(points.size(0));
    auto facet = points - centre;
    // pick the first on-face vector to be our x-axis
    auto x = facet.view(0).decouple(); // decouple since we're going to normalise this
    auto z = three_point_normal(facet, 0, 1, 2);
    auto y = cross(z, x);
    x /= norm(x);
    y /= norm(y);
    T a[2] = {1, 0}, b[2], c[2];
    T s{0};
    for (ind_t i = 1; i < facet.size(0) - 1; ++i) {
      on_plane_vector(x, y, facet.view(i), b);
      on_plane_vector(x, y, facet.view(i+1), c);
      auto piece = orient2d(a, b, c);
      if (!strict) piece = std::abs(piece);
      s += piece;
    }
    return (s > T(0));
  }

  /*! \brief Check if a point is 'behind' a plane

  \param a the first point defining the plane
  \param b the second point defining the plane
  \param c the third point defining the plane
  \param x the point or points to check
  \returns whether `x` is on or behind the plane

  \note Uses the geometry predicates orient3d to perform the check within machine precision
  */
  template<class T, template<class> class A>
  std::enable_if_t<isBareArray<T,A>, std::vector<bool>>
  point_inside_plane(const A<T>& a, const A<T>& b, const A<T>& c, const A<T>& x) {
    assert(a.numel() == 3 && b.numel() == 3 && c.numel() == 3);
    assert(a.is_contiguous() && b.is_contiguous() && c.is_contiguous());
    assert(x.is_contiguous() && x.is_row_ordered() && x.size(1u) == 3);
    std::vector<bool> pip;
    pip.reserve(x.size(0u));
    const auto a_ptr{a.ptr(0)};
    const auto b_ptr{b.ptr(0)};
    const auto c_ptr{c.ptr(0)};
    for (ind_t i = 0; i < x.size(0u); ++i) {
      auto o3d = orient3d(a_ptr, b_ptr, c_ptr, x.ptr(i));
      pip.push_back(o3d >= 0.);
    }
    return pip;
  }

  template<class T, template<class> class A>
  std::enable_if_t<isBareArray<T,A>, std::vector<bool>>
  point_inside_planes(const A<T>& a, const A<T>& b, const A<T>& c, const A<T>& x){
    auto n_planes = a.size(0);
    assert(n_planes == b.size(0u) && n_planes == c.size(0u));
    assert(a.size(1u) == 3u && b.size(1u) == 3u && c.size(1u) == 3u);
    assert(a.is_contiguous() && b.is_contiguous() && c.is_contiguous());
    assert(a.is_row_ordered() && b.is_row_ordered() && c.is_row_ordered());
    assert(x.numel() == 3 && x.is_contiguous());
    std::vector<bool> pip;
    pip.reserve(n_planes);
    const auto x_ptr{x.ptr(0)};
    for (ind_t i=0; i<n_planes; ++i){
      auto o3d = orient3d(a.ptr(i), b.ptr(i), c.ptr(i), x_ptr);
      pip.push_back(o3d >= 0.);
    }
    return pip;
  }

  template<class T, template<class> class A>
  std::enable_if_t<isBareArray<T,A>, bool>
  point_inside_all_planes(const A<T>& a, const A<T>& b, const A<T>& c, const A<T>& x){
      auto o3d = pseudo_orient3d(a, b, c, x);
      auto v = std::min_element(o3d.begin(), o3d.end());
//      if (*v < 0 && *v > -1e-10){
//        std::vector<T> below;
//        std::copy_if(o3d.begin(), o3d.end(), std::back_inserter(below), [](const T& z){return z < 0;});
//        std::transform(below.begin(), below.end(), below.begin(), [](const T& z){return z * 1e10;});
//        info_update("The point ", x.to_string(), "has negative parallelepiped volumes ", below, "e10");
//      }
      return *v > 0 || approx::scalar(*v, 0.);
  }

  template<class T, class R, template<class> class A, template<class> class B>
  std::enable_if_t<isLatVec<T,A> && isLatVec<R,B>, bool>
  point_in_plane_lattice_check(const A<T>& a, const A<T>& b, const A<T>& c, const B<R>& x){
    bool abc = a.same_lattice(b) && a.same_lattice(c);
    bool ax = x.same_lattice(a) || x.star_lattice(a);
    verbose_update_if(!abc, "a, b, and c should have the same lattice!");
    verbose_update_if(!ax, "a and x should have the same or dual lattices");
    return ax && abc;
  }

  #define POINT_IN_PLANE(FUNCTION, OUT_TYPE) \
  template<class T, class R, template<class> class A, template<class> class B>\
  std::enable_if_t<isLatVec<T,A> && isLatVec<R,B>, OUT_TYPE>\
  FUNCTION(const A<T>& a, const A<T>& b, const A<T>& c, const B<R>& x){\
    assert(point_in_plane_lattice_check(a, b, c, x));\
    return FUNCTION(a.xyz(), b.xyz(), c.xyz(), x.xyz());\
  }
  POINT_IN_PLANE(point_inside_plane, std::vector<bool>)
  POINT_IN_PLANE(point_inside_planes, std::vector<bool>)
  POINT_IN_PLANE(point_inside_all_planes, bool)

  #undef POINT_IN_PLANE

  template<class I>
  std::vector<I> two_in_one_indexes(const std::vector<std::vector<I>> &one, const std::vector<I> &two) {
    std::vector<I> indexes;
    indexes.reserve(one.size());
    auto has = [one](const size_t &i, const I &x) {
      return std::find(one[i].begin(), one[i].end(), x) != one[i].end();
    };
    for (size_t index = 0; index < one.size(); ++index) {
      if (std::any_of(two.begin(), two.end(), [index, has](const I &x) { return has(index, x); }))
        indexes.push_back(static_cast<I>(index));
    }
    return indexes;
  }

  template<class I>
  std::vector<I> keep_to_cut_list(const std::vector<bool> &keep, const std::vector<std::vector<I>> &faces) {
    std::vector<I> to_remove;
    for (size_t j = 0; j < keep.size(); ++j) if (!keep[j]) to_remove.push_back(static_cast<I>(j));
    return two_in_one_indexes(faces, to_remove);
  }

  /*! \brief Return edge vertex pairs plus containing face indexes for cut edges
   *
   * \returns std::vector<std::tuple<first_face, second_face, edge_vertex_pair>>
   * */
  template<class I>
  std::vector<std::tuple<size_t, size_t, std::pair<I,I>>> keep_to_cut_edge_list(const std::vector<bool> & keep, const std::vector<std::vector<std::pair<I,I>>> & edges){
    std::vector<std::tuple<size_t, size_t, std::pair<I,I>>> out;
    for (size_t j=0; j < keep.size(); ++j) {
      if (!keep[j]) { // if this vertex is removed, check all faces for an edge with it as the *first* vertex
        for (size_t first = 0; first < edges.size(); ++first) {
          for (const auto &edge: edges[first]) {
            // if the polyhedron is closed then every edge has an opposite pair on another face
            // therefore we only check if the *first* edge vertex is j.
            if (edge.first == j && keep[edge.second]) {
              // now look for the opposite edge containing face
              for (size_t second = 0; second < edges.size(); ++second) {
                if (second != first) {
                  for (const auto &oppo: edges[second]) {
                    if (oppo.first == edge.second && oppo.second == j) {
                      out.emplace_back(first, second, edge);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    return out;
  }


  /*! \brief Find the common-edge vertex indices between two facets

  \param  one the first facet vertex indexes
  \param  two the second facet vertex indexes
  \param  strict whether the two vertex indexes must be neighbours and in opposite order on the two facets
  \returns A tuple containing a success value and the two indexes; if there is no common edge the success value is false.
  */
  template<class I>
  std::tuple<bool, I, I> common_edge(const std::vector<I> &one, const std::vector<I> &two, const bool strict = false) {
    bool ok{false};
    I a{0}, b{0};
    if (strict) {
      for (size_t i = 0; !ok && i < one.size(); ++i) {
        a = one[i];
        b = one[(i + 1) % one.size()];
        for (size_t j = 0; !ok && j < two.size(); ++j) {
          ok = b == two[j] && a == two[(j + 1) % two.size()];
        }
      }
    } else {
      for (size_t ia = 0; !ok && ia < one.size(); ++ia) {
        a = one[ia];
        for (size_t ib = ia + 1; !ok && ib < one.size() + 1; ++ib) {
          b = one[ib % one.size()];
          for (size_t j = 0; !ok && j < two.size(); ++j) {
            const auto &ta{two[j]};
            for (size_t k = j + 1; !ok && k < two.size() + 1; ++k) {
              const auto &tb{two[k % two.size()]};
              ok = (ta == a && tb == b) || (tb == a && ta == b);
            }
          }
        }
      }
    }
    return std::make_tuple(ok, a, b);
  }

  /*! \brief Find the intersection of the edge between two facets and a plane defined by three points
  \param v      the vertices of the facets
  \param vpf    the vertex indices for each known facet
  \param j      the first facet index (accesses vpf)
  \param k      the second facet index (accesses vpf)
  \param a      the first point on the plane
  \param b      the second point on the plane
  \param c      the third point on the plane
  \param strict whether the facet vertex indices are sorted correctly
  \returns A bArray containing the intersection point, or an empty bArray if there is no intersection
  */
  template<class T, template<class> class A, class I>
  std::enable_if_t<isBareArray<T,A>, A<T>> edge_plane_intersection(
    const A<T> &v, const std::vector<I> &one, const std::vector<I> &two,
    const A<T> &a, const A<T> &b, const A<T> &c, const int tol = 1, const bool strict = false) {
    // find the correct pair of vertices which form the edge:
    auto[ok, j, k] = common_edge(one, two, strict);
    if (ok) {
      // check whether they're on opposite sides of the plane (a, b, c)
//      auto o3dj = orient3d(a.ptr(0), b.ptr(0), c.ptr(0), v.ptr(j));
//      auto o3dk = orient3d(a.ptr(0), b.ptr(0), c.ptr(0), v.ptr(k));
      auto o3dj = pseudo_orient3d(a, b, c, v.view(j))[0];
      auto o3dk = pseudo_orient3d(a, b, c, v.view(k))[0];
      // both on one side, both on the other side, or both *in* the plane all preclude a single intersection
      if ((o3dj < 0 && o3dk < 0) || (o3dj > 0 && o3dk > 0) || (o3dj == 0 && o3dk == 0)) ok = false;
      if (brille::approx::scalar(o3dj, 0., tol) && brille::approx::scalar(o3dk, 0., tol)) ok = false;
    }
    A<T> at;
    if (ok) {
      auto plane_n = three_point_normal(a, b, c);
      auto line_u = v.view(k) - v.view(j);
      auto denominator = dot(plane_n, line_u);
      auto numerator = dot(plane_n, a - v.view(j));
      assert(denominator.abs().sum() > 0);
      at = v.view(j) + (numerator * line_u) / denominator;
//      auto scalar = (numerator / denominator).sum();
//      at = v.view(j) + scalar * line_u;
//      debug_update_if(orient3d(a.ptr(0), b.ptr(0), c.ptr(0), at.ptr(0)) != 0.,
//                      "The found intersection point ", at.to_string(0),
//                      " is off the plane, proportional to ", orient3d(a.ptr(0), b.ptr(0), c.ptr(0), at.ptr(0)));
    }
    return at;
  }

  template<class T, class R, template<class> class A, template<class> class B, class I>
  std::enable_if_t<isLatVec<T,A> && isLatVec<R,B>, A<T>> edge_plane_intersection(
      const A<T> &v, const std::vector<I> &one, const std::vector<I> &two,
      const B<R> &a, const B<R> &b, const B<R> &c, const int tol=1, const bool strict = false) {
    assert(a.same_lattice(b) && a.same_lattice(c));
    if (v.same_lattice(a) || v.star_lattice(a)) {
      auto at = edge_plane_intersection(v.xyz(), one, two, a.xyz(), b.xyz(), c.xyz(), tol, strict);
      // replaced ::from_invA by from_xyz_like
      return at.size(0) > 0 ? from_xyz_like(v, at) : A<T>(v.type(), v.lattice(), 0u);
    }
    throw std::runtime_error("");
  }

  template<class T, template<class> class A, class I>
  std::enable_if_t<isArray<T,A>, A<T>> edge_plane_intersection(const A<T>& v, const std::pair<I,I> &edge, const A<T>& a, const A<T>& b, const A<T>& c){
    // we *know* that edge.first is on the wrong side of the plane  and edge.second is not.
    // we *know* that there must be an intersection between the two which might be either endpoint?
    auto plane_n = three_point_normal(a, b, c);
    auto line_u = v.view(edge.second) - v.view(edge.first);
    auto denominator = dot(plane_n, line_u);
    if (denominator.abs().sum() == 0){
      // this edge is *in* the plane ... it doesn't intersect, but one is expected. return an end?
      return v.view(edge.first);
    }
    return v.view(edge.first) + (dot(plane_n, a - v.view(edge.first)) * line_u) / denominator;
  }

  template<class T, template<class> class A, class I>
  std::enable_if_t<isArray<T,A>, std::vector<std::tuple<size_t, size_t, std::pair<I,I>, A<T>>>>
  valid_edge_plane_intersections(const A<T>& a, const A<T>& b, const A<T>& c, const A<T>& v, const std::vector<std::tuple<size_t, size_t, std::pair<I,I>>>& edges){
    auto plane_n = three_point_normal(a, b, c);
    std::vector<std::tuple<size_t, size_t, std::pair<I,I>, A<T>>> out;
    for (const auto & [first, second, edge]: edges){
      auto line_u = v.view(edge.second) - v.view(edge.first);
      auto denominator =  dot(plane_n, line_u);
      if (denominator.abs().sum() == 0) continue;
      auto factor = (dot(plane_n, a - v.view(edge.first)) / denominator).sum();
      if (factor < T(0) && approx::scalar(factor, T(0))){
        out.emplace_back(first, second, edge, v.view(edge.first));
        continue;
      }
      if (factor > T(1) && approx::scalar(factor, T(1))){
        out.emplace_back(first, second, edge, v.view(edge.second));
        continue;
      }
      if (factor < T(0) || factor > T(1)) {
        info_update("The factor ", factor, " is outside the valid range for edge ", my_to_string(edge));
        continue;
      }
      auto at = v.view(edge.first) + factor * line_u;
      out.emplace_back(first, second, edge, at);
    }
    return out;
  }

  template<class T, template<class> class A>
  std::enable_if_t<isBareArray<T, A>, std::tuple<Array2<T>, Array2<T>, Array2<T>>>
  find_convex_hull_planes(const A<T>& points){
    if (points.size(0) < 4)
      throw std::runtime_error("Not enough points to form a Convex Hull");
    auto bc = utils::binomial_coefficient(points.size(0), 3u);
    if (bc > static_cast<unsigned long long>(std::numeric_limits<ind_t>::max()))
      throw std::runtime_error("Too many points to count all possible bounding planes with an `ind_t` integer");
    auto ibc = static_cast<ind_t>(bc);
    A<T> n(ibc, 3u), a(ibc, 3u), b(ibc, 3u), c(ibc, 3u);
    ind_t count{0}, no{points.size(0)};
    for (ind_t i=0; i < no-2; ++i) for (ind_t j=i+1; j < no-1; ++j) for (ind_t k=j+1; k<no; ++k) {
      bool positive{true}, negative{true};
      for (ind_t r=0; (positive || negative) && r < no; ++r){
        if (r == i || r == j || r == k) continue;
        auto d = orient3d(points.ptr(i), points.ptr(j), points.ptr(k), points.ptr(r));
        positive &= d >= 0;
        negative &= d <= 0;
      }
      if (positive ^ negative) {
        // (T,T) -> co-linear (i,j,k) => all points are co-planar
        // (F,F) -> plane (i,j,k) does not divide space into point-full and point-less
        a.set(count, points.view(positive ? i : j));
        b.set(count, points.view(positive ? j : i));
        c.set(count, points.view(k));
        n.set(count, three_point_normal(a.view(count), b.view(count), c.view(count)));
        ++count;
      }
    }
    if (ibc < count)
      throw std::logic_error("Too many planes found");
    n.resize(count);
    a.resize(count);
    b.resize(count);
    c.resize(count);
    auto u = n.is_unique();
    return std::make_tuple(a.extract(u), b.extract(u), c.extract(u));
  }

  template<class T, template<class> class A>
  std::enable_if_t<isLatVec<T,A>, std::tuple<A<T>, A<T>, A<T>>>
  find_convex_hull_planes(const A<T>& points){
    auto [a, b, c] = find_convex_hull_planes(points.xyz());
    return std::make_tuple(from_xyz_like(points, a), from_xyz_like(points, b), from_xyz_like(points, c));
  }

  template<class T, template<class> class A>
  std::enable_if_t<isBareArray<T, A>, std::vector<std::vector<ind_t>>>
  find_planes_containing_point(const A<T>& a, const A<T>& b, const A<T>& c, const A<T>& points){
    std::vector<std::vector<ind_t>> result;
    result.reserve(points.size(0));
    for (ind_t i=0; i < points.size(0); ++i){
      const auto d{points.ptr(i)};
      std::vector<ind_t> plane;
      plane.reserve(a.size(0));
      for (ind_t j=0; j < a.size(0); ++j){
        if (orient3d(a.ptr(j), b.ptr(j), c.ptr(j), d) == 0) plane.push_back(j);
      }
      result.push_back(plane);
    }
    return result;
  }

  template<class T, template<class> class A>
  std::enable_if_t<isLatVec<T,A>, std::vector<std::vector<ind_t>>>
  find_planes_containing_point(const A<T>& a, const A<T>& b, const A<T>& c, const A<T>& points){
    // if they're all in the same lattice, can we act directly on (hkl)?
    // No. But maybe we can only worry about the co-angles and ignore the basis vector lengths?
    return find_planes_containing_point(a.xyz(), b.xyz(), c.xyz(), points.xyz());
  }

  template<class T, template<class> class A>
  std::enable_if_t<isBareArray<T,A>, std::vector<ind_t>>
  sort_convex_polygon_face(const A<T>& a, const A<T>& b, const A<T>& c, const std::vector<ind_t>& face, const A<T>& points){
    auto vectors = points.extract(face);
    auto centre = vectors.sum(0) / static_cast<double>(vectors.size(0));
    auto good = vectors.row_is(cmp::neq, centre);
    std::vector<ind_t> pruned;
    pruned.reserve(face.size());
    if (!good.all()){
      debug_update("Not all vectors are non-centre vectors?\n\tvectors\n", vectors.to_string(),"centre ", centre.to_string(0));
      // do something about the bad one(s)
      vectors = vectors.extract(good);
      for (ind_t i=0; i<face.size(); ++i) if (good.val(i)) pruned.push_back(face[i]);
    } else {
      std::copy(face.begin(), face.end(), std::back_inserter(pruned));
    }
    // this is why we removed the bad point(s)
    vectors -= centre;
    vectors /= norm(vectors);

    auto normal = three_point_normal(a, b, c);

//    debug_update("Input face ", face, " gives pruned ", pruned, " centre ", centre.to_string(0));
//    debug_update("plus normal ", normal.to_string(0), " vectors\n", vectors.to_string());

    std::vector<ind_t> perm(pruned.size(), 0), used(1, 0);
    used.reserve(pruned.size());
    auto angles = bare_winding_angles(vectors, 0, normal);
    angles[0] = (std::numeric_limits<T>::max)();
    for (ind_t i=1; i < pruned.size(); ++i){
      auto itr = std::min_element(angles.begin(), angles.end());
      perm[i] = std::distance(angles.begin(), itr);
      *itr = angles[0];
    }

//    for (ind_t i=1; i < pruned.size(); ++i){
//      auto angles = bare_winding_angles(vectors, perm[i-1], normal);
//      double min_angle = 1e3;
//      auto min_idx = static_cast<ind_t>(pruned.size());
//      for (ind_t k=0; k < pruned.size(); ++k){
//        auto unused = std::find(used.begin(), used.end(), k) == used.end();
//        auto not_zero = !approx::scalar(angles[k], 0.);
//        if (unused && not_zero && angles[k] < min_angle){
//          min_idx = k;
//          min_angle = angles[k];
//        }
//      }
//      if (min_idx >= pruned.size()){
//        std::string msg = "Error finding minimum winding angle vertex\n";
//        for (ind_t i0=0; i0<vectors.size(0)-1; ++i0) for (ind_t i1=i0+1; i1 < vectors.size(0); ++i1){
//          auto val = (vectors.view(i0) - vectors.view(i1)).abs().sum();
//          if (val < 1e-8) {
//            msg += "(" + std::to_string(i0) + "," + std::to_string(i1) + ") -> ";
//            msg += std::to_string(val *1e15) + "x10^-15\n";
//          }
//        }
//        msg += "points\n";
//        for (const auto & v: pruned) msg += points.to_string(v) + "\n";
//        msg += "vectors\n" + vectors.to_string();
//        msg += "winding angles " + list_to_string(angles) + "\n";
//        throw std::runtime_error(msg);
//      }
//      perm[i] = min_idx;
//      used.push_back(min_idx);
//    }

    std::vector<ind_t> sorted;
    sorted.reserve(pruned.size());
    for (const auto & p: perm) sorted.push_back(pruned[p]);
//    verbose_update("The permutation", perm, " sorts the vertices to ", sorted);

    // finally, go through sorted to ensure there are no extraneous points hanging about
    for (ind_t i=0, j; sorted.size() > 3 && i < sorted.size();){
      j = (sorted.size() + i - 1) % sorted.size(); // why not just (i-1) % sorted.size()?
      auto prev = points.view(sorted[i]) - points.view(sorted[j]);
      j = (sorted.size() + i + 1) % sorted.size();
      auto next = points.view(sorted[j]) - points.view(sorted[i]);
      if (dot(normal, cross(prev, next)).all(cmp::gt, 0.)){
        // left turn; we keep this point
        ++i;
      } else{
        // right turn, we remove this point
        sorted.erase(sorted.begin() + i);
      }
    }
    return sorted;
  }

  template<class T, template<class> class A>
  std::enable_if_t<isLatVec<T,A>, std::vector<ind_t>>
  sort_convex_polygon_face(const A<T>& a, const A<T>& b, const A<T>& c, const std::vector<ind_t>& face, const A<T>& points){
    assert(a.same_lattice(b) && a.same_lattice(c) && a.same_lattice(points));
    return sort_convex_polygon_face(a.xyz(), b.xyz(), c.xyz(), face, points.xyz());
  }

  template<class T, class I, template<class> class A>
  std::enable_if_t<isArray<T,A>, std::vector<ind_t>>
  sort_convex_polygon_face(const std::vector<I>& face, const A<T>& vertices){
    auto a = vertices.view(face[0]);
    auto b = vertices.view(face[1]);
    auto c = vertices.view(face[2]);
    return sort_convex_polygon_face(a, b, c, face, vertices);
  }

  template<class T, template<class> class A, class I>
  std::enable_if_t<isArray<T, A>, std::vector<std::vector<ind_t>>>
  polygon_faces(const std::vector<std::vector<I>>& faces, const A<T>& points){
    std::vector<std::vector<I>> reduced_faces;
    reduced_faces.reserve(faces.size());
    for (const auto & face: faces) if (face_has_area(points.extract(face)))
      reduced_faces.push_back(sort_convex_polygon_face(face, points));
    return reduced_faces;
  }

  template<class T, template<class> class A, class I>
  std::enable_if_t<isArray<T,A>, std::vector<std::vector<I>>>
  polygon_faces(const A<T>& a, const A<T>& b, const A<T>& c, const std::vector<std::vector<I>>& faces_per_point, const A<T>& points){
    auto faces = utils::invert_lists(faces_per_point);

    auto no_faces = faces.size();
    if (a.size(0) != no_faces || b.size(0) != no_faces || c.size(0) != no_faces)
      throw std::runtime_error("Wrong number of plane points for the number of faces");

    std::vector<bool> is_polygon;
    is_polygon.reserve(no_faces);
    std::transform(faces.begin(), faces.end(), std::back_inserter(is_polygon),
                   [points](const auto & face) { return face_has_area(points.extract(face)); });

    std::vector<std::vector<I>> reduced_faces;
    reduced_faces.reserve(no_faces);
    for (ind_t i=0; i<no_faces; ++i) if (is_polygon[i]) {
      reduced_faces.push_back(sort_convex_polygon_face(a.view(i), b.view(i), c.view(i), faces[i], points));
    }

    return reduced_faces;
  }


  template<class T, template<class> class A>
  std::enable_if_t<isArray<T,A>, bool>
  polygon_face_vertex_purge(A<T>& points, std::vector<std::vector<ind_t>>& faces){
    std::vector<bool> keep(points.size(0), false);
    for (ind_t i=0; i<points.size(0); ++i) for(const auto & face: faces) {
      if (std::find(face.begin(), face.end(), i) != face.end()) {
        keep[i] = true;
        break;
      }
    }
    auto total = std::count(keep.begin(), keep.end(), true);
    if (total < points.size(0)){
      ind_t count{0};
      std::vector<ind_t> nim;
      nim.reserve(total);
      for (const auto & b: keep) nim.push_back(b ? count++ : total);
      for (auto & face: faces) for (auto & v: face) if (nim[v] < total) v = nim[v];
      points = points.extract(keep);
      return true;
    }
    return false;
  }

  template<class T, template<class> class A>
  std::enable_if_t<isArray<T,A>, std::tuple<A<T>, std::vector<std::vector<ind_t>>>>
  remove_duplicate_points_and_update_face_indexing(const A<T>& points, const std::vector<std::vector<ind_t>> faces){
    auto are_unique = points.is_unique();
    if(std::find(are_unique.begin(), are_unique.end(), false) != are_unique.end()){
      std::vector<ind_t> index;
      ind_t count{0}, no{points.size(0)};
      index.reserve(no);
      for (const auto & b: are_unique) index.push_back(b ? count++ : no);
      for (ind_t i=0; i < no; ++i) {
        if (index[i] >= no) {
          for (ind_t j = 0; j < no; ++j) {
            if (i != j && index[j] < no && points.match(i, j)) index[i] = index[j];
          }
        }
      }
      std::vector<std::vector<ind_t>> new_faces;
      new_faces.reserve(faces.size());
      for (auto & face: faces){
        std::vector<ind_t> one;
        one.reserve(face.size());
        for (auto & f: face) one.push_back(index[f]);
        new_faces.push_back(one);
      }
      return std::make_tuple(points.extract(are_unique), new_faces);
    }
    return std::make_tuple(points, faces);
  }

  template<class T, template<class> class A>
  std::enable_if_t<isArray<T,A>, std::vector<ind_t>>
  remove_middle_colinear_points(const A<T>& points, const std::vector<ind_t>& face, const int tol){
    std::vector<ind_t> updated;
    updated.reserve(face.size());
    updated.push_back(face[0]);
    const auto s{face.size()};
    for (ind_t i=1; i < s; ++i){
      auto a = points.view(face[(i-1) % s]);
      auto b = points.view(face[i]);
      auto c = points.view(face[(i+1) % s]);
      if (!cross(b-a, c-b).all(cmp::eq, 0., tol)) updated.push_back(face[i]);
    }
    return updated;
  }
  template<class T, template<class> class A>
  std::enable_if_t<isArray<T,A>, std::vector<std::vector<ind_t>>>
  remove_middle_colinear_points(const A<T>& points, const std::vector<std::vector<ind_t>>& faces, const int tol){
    std::vector<std::vector<ind_t>> updated;
    updated.reserve(faces.size());
    for (const auto & face : faces){
      auto one = remove_middle_colinear_points(points, face, tol);
      if (one.size() > 2) updated.push_back(one);
    }
//    verbose_update("After removing colinear face points\n", updated);
    return updated;
  }

  template<class T, template<class> class A>
  std::enable_if_t<isArray<T,A>, std::tuple<A<T>, std::vector<std::vector<ind_t>>>>
  remove_points_and_update_face_indexing(const std::vector<bool> & keep, const A<T>& points, const std::vector<std::vector<ind_t>> & faces, const bool sort=true){
    if (std::find(keep.begin(), keep.end(), false) != keep.end()) {
      verbose_update("Prune points; keep=", keep);
      std::vector<ind_t> map(points.size(0), points.size(0));
      ind_t count{0};
      for (ind_t i=0; i < points.size(0); ++i) if (keep[i]) map[i] = count++;
      verbose_update("Old to new mapping: ", map);
      auto kept = points.extract(keep);

      std::vector<std::vector<ind_t>> updated;
      for (auto & face: faces) {
        std::vector<ind_t> one;
        for (const auto & x: face) if (keep[x]) one.push_back(map[x]);
        if (one.size() < 3) continue;
        if (sort) {
          if (!(keep[face[0]] && keep[face[1]] && keep[face[2]])) {
            auto old_normal = three_point_normal(points, face);
            auto new_normal = three_point_normal(kept, one);
            if (dot(old_normal, new_normal).all(cmp::le, 0.)) {
              std::swap(one[1], one[2]);
            }
          }
          one = sort_convex_polygon_face(one, kept);
        }
        updated.push_back(one);
      }
      verbose_update("Old faces\n", faces, "becomes\n", updated);
      return std::make_tuple(kept, updated);
    }
    return std::make_tuple(points, faces);
  }

  template<class I>
  std::vector<std::vector<I>>
  remove_dangling_faces(const I no, const std::vector<std::vector<I>>& faces){
    // does a face have the index
    auto has = [](const auto & face, const I & index){
      return std::find(face.begin(), face.end(), index) != face.end();
    };
    // copy the faces which are ok to a new vector of vectors
    auto loop = [no, has](const auto & source){
      if (source.size() < 1) return std::make_tuple(false, source);
      verbose_update("De-dangling loop for ", source);
      // count how many times each vertex appears in the faces lists
      std::vector<I> co(no);
      std::iota(co.begin(), co.end(), 0u);
      verbose_update("vertices ", co);
      std::transform(co.begin(), co.end(), co.begin(), [source, has](const I i){
        I c{0};
        for (const auto & f : source) if (has(f, i)) ++c;
        return c;
      });
      verbose_update("  counts ", co);
      // for each face, remove any points which appear on less than three faces
      // and then remove faces with less than three vertices
      auto checked_face = [co](const auto & f){
        std::vector<I> face;
        face.reserve(f.size());
        for (const auto & x: f) if (co[x] > 2u) face.push_back(x);
        return face;
      };
      std::vector<std::vector<I>> checked;
      checked.reserve(source.size());
      std::transform(source.begin(), source.end(), std::back_inserter(checked), checked_face);
      std::vector<std::vector<I>> valid;
      valid.reserve(checked.size());
      std::copy_if(checked.begin(), checked.end(), std::back_inserter(valid), [](const auto & f){return f.size() > 2u;});
      //
      return std::make_tuple(source.size() > valid.size(), valid);

      /*
       The following version works if all faces are fully formed without
       extraneous edge points. If there are extraneous edge points they may
       appear on two faces which should be retained after the extraneous points
       are removed.

      // copy a face if all of its vertices appear at least twice
      auto check = [co](const auto & f){
        return std::all_of(f.begin(), f.end(), [co](const I &x) {return co[x] > 2u;});
      };
      std::vector<std::vector<I>> sink;
      sink.reserve(source.size());
      std::copy_if(source.begin(), source.end(), std::back_inserter(sink), check);
      return std::make_tuple(source.size() > sink.size(), sink);
      */
    };
    // actually run the good-face reduction until all faces are good (or gone)
    auto [again, out] = loop(faces);
    while (again) std::tie(again, out) = loop(out);
    verbose_update("De-dangled faces:\n", out);
    return out;
  }

  template<class T, template<class> class A, class I>
  std::enable_if_t<isArray<T,A>, std::tuple<A<T>, std::vector<std::vector<I>>>>
  remove_faceless_points(const A<T>& points, const std::vector<std::vector<I>> & faces){
    std::vector<bool> keep(points.size(0), false);
    for (const auto & face: faces) for (const auto & index: face) keep[index] = true;
    if (std::find(keep.begin(), keep.end(), false) != keep.end()) {
      I count{0}, no{points.size(0)};
      std::vector<I> map;
      map.reserve(no);
      for (const auto & b: keep) map.push_back(b ? count++ : no);
      std::vector<std::vector<I>> updated;
      updated.reserve(faces.size());
      for (const auto & face: faces) {
        std::vector<I> one;
        std::transform(face.begin(), face.end(), std::back_inserter(one), [map](const auto & i){return map[i];});
        updated.push_back(one);
      }
      verbose_update("Keep ", keep, " face-full points\nnp.array(", get_xyz(points), "),\n", updated);
      return std::make_tuple(points.extract(keep), updated);
    }
    return std::make_tuple(points, faces);
  }

  template<class I> bool face_not_in_faces(const std::vector<I>& candidate, const std::vector<std::vector<I>>& faces){
    // check for all-matching vertices
    if (utils::unordered_list_in_lists(candidate, faces)) return false;
    // check for disallowed situations caused by rounding errors:
    for (const auto & face: faces){
      size_t count{0};
      // should we worry about faces that visit the same point more than once?
      for (const auto & v: candidate) if (std::find(face.begin(), face.end(), v) != face.end()) ++count;
      // any faces that share more than two (unique, non-co-linear) points must be coplanar
      verbose_update_if(count > 2, "The candidate face ", candidate, " matches existing face ", face);
      if (count > 2) return false;
    }
    return true;
  }

}
#endif