#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "exact_geometry.hpp"

using namespace brille::exact;

TEST_CASE("Expansions are exact", "[exact]") {
  // (1e16 + 1) - 1e16 is 1, but in doubles it is 0 or 2
  const Expansion big(1e16), one(1.0);
  REQUIRE(((big + one) - big).sign() == 1);
  REQUIRE(((big + one) - big - one).sign() == 0);
  // (a + b)² - a² - 2ab - b² == 0 exactly, for awkward a, b
  const Expansion a(0.1), b(3e-17);
  const Expansion s = a + b;
  REQUIRE((s * s - a * a - (a * b).scaled(2.0) - b * b).sign() == 0);
  REQUIRE((s * s - a * a - (a * b).scaled(2.0)).sign() == 1);
}

TEST_CASE("Expansions longer than their inline storage stay exact", "[exact]") {
  // 18 non-overlapping components, more than are stored in place: 54 bits apart, so
  // no two fit in one double, and spanning few enough that squares neither overflow
  // nor underflow
  Expansion sum;
  std::vector<double> parts;
  for (int k = 0; k < 18; ++k) parts.push_back(std::ldexp(1.0, 480 - 54 * k));
  for (const double x: parts) sum += Expansion(x);
  REQUIRE(sum.size() == 18u);
  REQUIRE(sum.sign() == 1);
  // removing all but the smallest leaves exactly the smallest
  Expansion rest = sum;
  for (size_t k = 0; k + 1 < parts.size(); ++k) rest -= Expansion(parts[k]);
  REQUIRE((rest - Expansion(parts.back())).sign() == 0);
  REQUIRE(rest.sign() == 1);
  // and a product of long expansions is exact too: (s)(s) - s² computed termwise
  const Expansion square = sum * sum;
  Expansion termwise;
  for (const double x: parts) termwise += sum.scaled(x);
  REQUIRE((square - termwise).sign() == 0);
  REQUIRE(((sum - sum)).sign() == 0);
}

namespace {
// a hexagonal reciprocal metric, invariant under the 6-fold rotation below
Metric hexagonal() {
  const double a = 1.3, c = 0.7;
  return Metric({a, a / 2, 0, a / 2, a, 0, 0, 0, c});
}
const mat3i sixfold{0, -1, 0, 1, 1, 0, 0, 0, 1};   // on reciprocal coordinates
}

TEST_CASE("Metric invariance is checked exactly", "[exact]") {
  REQUIRE(hexagonal().invariant_under(sixfold));
  auto g = hexagonal().values();
  g[1] = std::nextafter(g[1], 1.0);
  g[3] = g[1];
  REQUIRE_FALSE(Metric(g).invariant_under(sixfold));
}

TEST_CASE("Mapping named points agrees with mapping coordinates", "[exact]") {
  const Geometry geom(hexagonal());
  REQUIRE(geom.metric().invariant_under(sixfold));
  // a zone vertex: three zone planes x·Gτ = τᵀGτ/2
  const Point p{{Plane::metric_plane({1, 0, 0}, {1, 0, 0}), Plane::metric_plane({0, 1, 0}, {0, 1, 0}), Plane::metric_plane({0, 0, 1}, {0, 0, 1})}};
  REQUIRE(geom.independent(p) != 0);
  const auto x = geom.coordinates(p);
  const int3 t{1, -2, 1};
  const auto q = geom.mapped(p, sixfold, t);
  const auto y = geom.coordinates(q);
  for (int i = 0; i < 3; ++i) {
    double expected = t[i];
    for (int k = 0; k < 3; ++k) expected += static_cast<double>(sixfold[3 * i + k]) * x[k];
    REQUIRE_THAT(y[i], Catch::Matchers::WithinAbs(expected, 1e-12));
  }
  // every plane through p maps to a plane through q
  for (const auto & pl: p.planes) REQUIRE(geom.on(q, pl.mapped(sixfold, t)));
  // and an integer plane through the origin stays exact under the mapping
  const Point o{{Plane::integer_plane({1, 0, 0}, 0), Plane::integer_plane({0, 1, 0}, 0), Plane::integer_plane({0, 0, 1}, 0)}};
  const auto o2 = geom.mapped(o, sixfold, t);
  const auto z = geom.coordinates(o2);
  for (int i = 0; i < 3; ++i) REQUIRE(z[i] == static_cast<double>(t[i]));
}

TEST_CASE("Points are cached per geometry, not shared between metrics", "[exact]") {
  // the same named point under two metrics has different coordinates and sides
  const Point p{{Plane::metric_plane({1, 0, 0}, {1, 0, 0}), Plane::metric_plane({0, 1, 0}, {0, 1, 0}), Plane::metric_plane({0, 0, 1}, {0, 0, 1})}};
  const Geometry cubic(Metric({1, 0, 0, 0, 1, 0, 0, 0, 1}));
  const Geometry hexagonal_one(hexagonal());
  const auto x = cubic.coordinates(p);
  const auto y = hexagonal_one.coordinates(p);
  REQUIRE_THAT(x[0], Catch::Matchers::WithinAbs(0.5, 1e-15));
  REQUIRE_THAT(y[0], Catch::Matchers::WithinAbs(1.0 / 3.0, 1e-15));   // x + y/2 = 1/2 and x/2 + y = 1/2
  // a plane through one point's position but not the other's
  const Plane half = Plane::integer_plane({2, 0, 0}, 1);   // x = 1/2
  REQUIRE(cubic.on(p, half));
  REQUIRE_FALSE(hexagonal_one.on(p, half));
  // and again, now that both are cached
  REQUIRE(cubic.on(p, half));
  REQUIRE(hexagonal_one.side(p, half) < 0);
  // a new geometry at the same address gets its own entries
  for (int k = 0; k < 3; ++k) {
    const Geometry g(k % 2 ? hexagonal() : Metric({1, 0, 0, 0, 1, 0, 0, 0, 1}));
    REQUIRE(g.on(p, half) == (k % 2 == 0));
  }
}
