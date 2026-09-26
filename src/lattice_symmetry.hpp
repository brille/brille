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
#ifndef BRILLE_LATTICE_SYMMETRY_HPP_
#define BRILLE_LATTICE_SYMMETRY_HPP_
#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace brille {
/*! \brief The number of symmetry operations of a lattice, to a tolerance

Counts the integer matrices R with det R = ±1 and RᵀGR = G, within `tolerance`
times the largest element of G, for G the metric of a primitive basis. These are
the lattice's own symmetries (its holohedry), whatever the crystal's symmetry:
2, 4, 8, 12, 16, 24 or 48.

Counting at a tight and at a loose tolerance tells whether a lattice is close
to a more symmetric one. Near such a lattice the Brillouin zone has faces or
edges far smaller than itself, where they would appear or vanish, and meshes
of it are slow and poorly shaped.

\param metric the metric tensor of a primitive basis, row-major
\param tolerance the largest change allowed, relative to the metric
*/
inline size_t lattice_symmetry_count(std::array<double, 9> metric, const double tolerance) {
  auto g = [&metric](const int i, const int j) -> double & { return metric[3 * i + j]; };
  // reduce the basis pairwise (bᵢ ← bᵢ - n bⱼ) so that its symmetries have small
  // entries; only the metric of the reduced basis is needed
  for (bool changed = true; changed;) {
    changed = false;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        if (i == j) continue;
        const auto n = std::round(g(i, j) / g(j, j));
        if (n == 0.) continue;
        const double before = g(i, i);
        // column and row operations: G ← TᵀGT with T = I - n eⱼeᵢᵀ
        for (int k = 0; k < 3; ++k) g(k, i) -= n * g(k, j);
        for (int k = 0; k < 3; ++k) g(i, k) -= n * g(j, k);
        changed |= g(i, i) < before * (1 - 1e-12);
      }
  }
  double scale{0};
  for (const auto &x: metric) scale = std::max(scale, std::abs(x));
  const double limit = tolerance * scale;
  auto product = [&metric](const std::array<int, 3> &u, const std::array<int, 3> &v) {
    double x{0};
    for (int k = 0; k < 3; ++k)
      for (int l = 0; l < 3; ++l) x += u[k] * metric[3 * k + l] * v[l];
    return x;
  };
  // lattice vectors that could be images of each basis vector
  std::vector<std::array<int, 3>> vectors;
  for (int a = -2; a <= 2; ++a)
    for (int b = -2; b <= 2; ++b)
      for (int c = -2; c <= 2; ++c)
        if (a || b || c) vectors.push_back({a, b, c});
  std::array<std::vector<std::array<int, 3>>, 3> candidates;
  for (int i = 0; i < 3; ++i)
    for (const auto &v: vectors)
      if (std::abs(product(v, v) - metric[4 * i]) <= limit) candidates[i].push_back(v);
  size_t count{0};
  for (const auto &u: candidates[0])
    for (const auto &v: candidates[1]) {
      if (std::abs(product(u, v) - metric[1]) > limit) continue;
      for (const auto &w: candidates[2]) {
        if (std::abs(product(u, w) - metric[2]) > limit || std::abs(product(v, w) - metric[5]) > limit) continue;
        const int det = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0]);
        if (det == 1 || det == -1) ++count;
      }
    }
  return count;
}

//! The name of the holohedry with this many operations
inline std::string holohedry_name(const size_t operations) {
  switch (operations) {
    case 2: return "-1 (triclinic)";
    case 4: return "2/m (monoclinic)";
    case 8: return "mmm (orthorhombic)";
    case 12: return "-3m (rhombohedral)";
    case 16: return "4/mmm (tetragonal)";
    case 24: return "6/mmm (hexagonal)";
    case 48: return "m-3m (cubic)";
    default: return std::to_string(operations) + " operations";
  }
}
}
#endif
