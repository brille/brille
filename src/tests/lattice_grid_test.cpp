#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <random>
#include "lattice_dual.hpp"
#include "lattice_grid.hpp"

using namespace brille;
using namespace brille::latticetri;

namespace {
//! A grid of a lattice's (conventional) reciprocal lattice, with its point group
Grid reciprocal_grid(const std::array<double, 3> & lengths, const std::array<double, 3> & angles, const std::string & symmetry) {
  auto lat = lattice::Direct<double>(lengths, angles, symmetry);
  const auto b = lat.reciprocal_basis_vectors();     // columns
  std::array<double, 9> rows{};
  for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) rows[3 * i + k] = b[3 * k + i];
  // reciprocal lattice coordinates transform by Wᵀ
  std::vector<mat3i> ops;
  const auto ps = lat.pointgroup_symmetry();
  for (size_t j = 0; j < ps.size(); ++j) {
    const auto w = ps.get(j);
    mat3i r{};
    for (int a = 0; a < 3; ++a) for (int c = 0; c < 3; ++c) r[3 * a + c] = w[3 * c + a];
    ops.push_back(r);
  }
  return Grid(rows, ops);
}
}

TEST_CASE("Lattice grids are invariant under the point group", "[lattice_grid]") {
  struct Case { std::array<double, 3> l, a; std::string s; bool degenerate; };
  for (const auto & c: {Case{{4., 4., 4.}, {90., 90., 90.}, "-P 4 2 3", true},
                        Case{{3.3, 3.3, 3.3}, {90., 90., 90.}, "-I 4 2 3", true},   // the conventional reciprocal lattice is cubic P
                        Case{{3., 3., 5.}, {90., 90., 120.}, "-P 6 2", true},
                        Case{{4.02, 4.90, 3.29}, {98.97, 86.24, 88.47}, "-P 1", false}}) {
    auto g = reciprocal_grid(c.l, c.a, c.s);
    CAPTURE(c.s);
    REQUIRE(g.degenerate() == c.degenerate);
    REQUIRE(g.invariant(3.0));
  }
}

TEST_CASE("Lattice grid point location", "[lattice_grid]") {
  for (const auto & s: {std::string("-P 6 2"), std::string("-P 1")}) {
    auto g = s == "-P 1" ? reciprocal_grid({4.02, 4.90, 3.29}, {98.97, 86.24, 88.47}, s) : reciprocal_grid({3., 3., 5.}, {90., 90., 120.}, s);
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(-4, 4);
    for (int i = 0; i < 2000; ++i) {
      vec3 x{u(rng), u(rng), u(rng)};
      tetrahedron t{};
      std::array<double, 4> w{};
      REQUIRE(g.locate(x, t, w));
      vec3 y{0, 0, 0};
      for (int k = 0; k < 4; ++k) { const auto p = g.cartesian(t[k]); for (int j = 0; j < 3; ++j) y[j] += w[k] * p[j]; }
      for (int j = 0; j < 3; ++j) REQUIRE_THAT(y[j], Catch::Matchers::WithinAbs(x[j], 1e-10));
    }
  }
}
