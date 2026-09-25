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
#ifndef BRILLE_TETRAHEDRON_OVERLAP_HPP_
#define BRILLE_TETRAHEDRON_OVERLAP_HPP_
#include <array>
#include "tetgen.h"

namespace brille {
//! Exact geometric tests on tetrahedra given as four vertex pointers
namespace overlap {
  inline int sign(const double x) { return (x > 0) - (x < 0); }
  using tet_t = std::array<const double *, 4>;
  // faces as vertex indices, and the vertex opposite each
  constexpr std::array<std::array<int, 3>, 4> faces{{{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}}};
  constexpr std::array<std::array<int, 2>, 6> edges{{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};

  //! orient3d of face f of t with p, signed so that t's own opposite vertex is positive
  inline int side(const tet_t &t, const int f, const double *p, const int orientation) {
    const auto &v = faces[f];
    return orientation * sign(orient3d(t[v[0]], t[v[1]], t[v[2]], p));
  }
  //! p is inside or on the closed tetrahedron t
  inline bool contains(const tet_t &t, const int orientation, const double *p) {
    for (int f = 0; f < 4; ++f) if (side(t, f, p, orientation) < 0) return false;
    return true;
  }
  //! Segment pq meets closed triangle abc, or is coplanar with it
  inline bool segment_meets_triangle(const double *p, const double *q, const double *a, const double *b, const double *c) {
    const int sp = sign(orient3d(a, b, c, p));
    const int sq = sign(orient3d(a, b, c, q));
    if (sp == 0 && sq == 0) return true; // coplanar: report as touching
    if (sp == sq) return false;          // strictly on one side
    const int s0 = sign(orient3d(p, q, a, b));
    const int s1 = sign(orient3d(p, q, b, c));
    const int s2 = sign(orient3d(p, q, c, a));
    const bool positive = s0 > 0 || s1 > 0 || s2 > 0;
    const bool negative = s0 < 0 || s1 < 0 || s2 < 0;
    return !(positive && negative);
  }
  inline bool edges_meet_faces(const tet_t &e, const tet_t &f) {
    for (const auto &edge: edges)
      for (const auto &face: faces)
        if (segment_meets_triangle(e[edge[0]], e[edge[1]], f[face[0]], f[face[1]], f[face[2]])) return true;
    return false;
  }
}

/*! \brief Whether two closed tetrahedra share any point

Every decision is the sign of Shewchuk's adaptive `orient3d`, which is exact for
the given coordinates, so no new points are computed and no tolerance is needed.
Coplanar cases that would need a two-dimensional test are reported as touching:
callers use the answer to find candidate tetrahedra, where an extra candidate
costs a little time but a missing one loses points.

\param a the four vertices of one tetrahedron
\param b the four vertices of the other
*/
inline bool tetrahedra_overlap(const overlap::tet_t &a, const overlap::tet_t &b) {
  using namespace overlap;
  // orient3d(face, opposite vertex) has the same sign for every face of a tetrahedron
  const int oa = sign(orient3d(a[faces[0][0]], a[faces[0][1]], a[faces[0][2]], a[0]));
  const int ob = sign(orient3d(b[faces[0][0]], b[faces[0][1]], b[faces[0][2]], b[0]));
  if (oa == 0 || ob == 0) return true; // degenerate: cannot separate, so report touching
  // a face plane of one with all of the other strictly outside separates them
  for (int f = 0; f < 4; ++f) {
    bool all_outside{true};
    for (int k = 0; k < 4 && all_outside; ++k) all_outside = side(a, f, b[k], oa) < 0;
    if (all_outside) return false;
    all_outside = true;
    for (int k = 0; k < 4 && all_outside; ++k) all_outside = side(b, f, a[k], ob) < 0;
    if (all_outside) return false;
  }
  for (int k = 0; k < 4; ++k) if (contains(a, oa, b[k]) || contains(b, ob, a[k])) return true;
  // otherwise their boundaries cross, so an edge of one meets a face of the other
  return edges_meet_faces(a, b) || edges_meet_faces(b, a);
}
}
#endif
