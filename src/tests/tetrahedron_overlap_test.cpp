#include <catch2/catch_test_macros.hpp>
#include <array>
#include "tetrahedron_overlap.hpp"

using namespace brille;
using p3 = std::array<double, 3>;
using t4 = std::array<p3, 4>;

static bool overlaps(const t4 &a, const t4 &b) {
  overlap::tet_t pa{a[0].data(), a[1].data(), a[2].data(), a[3].data()};
  overlap::tet_t pb{b[0].data(), b[1].data(), b[2].data(), b[3].data()};
  const auto ab = tetrahedra_overlap(pa, pb);
  REQUIRE(ab == tetrahedra_overlap(pb, pa)); // symmetric
  return ab;
}

TEST_CASE("Tetrahedra overlap","[tetrahedron]"){
  const t4 unit{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
  SECTION("identical"){ REQUIRE(overlaps(unit, unit)); }
  SECTION("one inside the other"){
    const t4 small{{{.1, .1, .1}, {.2, .1, .1}, {.1, .2, .1}, {.1, .1, .2}}};
    REQUIRE(overlaps(unit, small));
  }
  SECTION("touching counts as overlapping"){
    const t4 face{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}}};
    const t4 vertex{{{1, 0, 0}, {2, 0, 0}, {1, 1, 0}, {1, 0, 1}}};
    REQUIRE(overlaps(unit, face));
    REQUIRE(overlaps(unit, vertex));
  }
  SECTION("a round-off gap is a gap"){
    const t4 apart{{{1 + 1e-12, 0, 0}, {2, 0, 0}, {1, 1, 0}, {1, 0, 1}}};
    REQUIRE_FALSE(overlaps(unit, apart));
  }
  SECTION("crossing edges, with no vertex inside the other tetrahedron"){
    const t4 a{{{-2, 0, 0}, {2, 0, 0}, {0, .2, 1}, {0, -.2, 1}}};
    const t4 b{{{0, -2, .5}, {0, 2, .5}, {.2, 0, -.5}, {-.2, 0, -.5}}};
    REQUIRE(overlaps(a, b));
  }
  SECTION("either vertex order"){
    const t4 flipped{{{0, 0, 0}, {0, 1, 0}, {1, 0, 0}, {0, 0, 1}}};
    const t4 far{{{3, 3, 3}, {4, 3, 3}, {3, 4, 3}, {3, 3, 4}}};
    REQUIRE(overlaps(flipped, unit));
    REQUIRE_FALSE(overlaps(flipped, far));
  }
}
