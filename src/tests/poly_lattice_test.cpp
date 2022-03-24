#include <catch2/catch.hpp>
#include <filesystem>

#include "array_.hpp" // defines bArray<T> as Array2<T> or Array<T>
#include "polyhedron_flex.hpp"
#include "bz.hpp"

using namespace brille;
using namespace brille::polyhedron;
using namespace brille::lattice;
using namespace brille::math;

TEST_CASE("Lattice Polyhedron","[polyhedron][lattice]") {
  double a_star{4.}, c_star{5.};
  double volume{half_root_three * a_star * a_star * c_star};
  auto lat = Reciprocal<double>({a_star, a_star, c_star}, {90., 90., 60.}, "P 1");
  auto bz = BrillouinZone(lat);
  auto fbz = bz.get_polyhedron();
  auto fbz_xyz = Poly(fbz.vertices().xyz(), fbz.faces());

  debug_update(fbz.python_string());

  REQUIRE(fbz_xyz.volume() == Approx(volume));
  REQUIRE(fbz.volume() == Approx(volume));
}