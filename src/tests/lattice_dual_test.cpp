#include <catch2/catch.hpp>

#include "lattice_dual.hpp"

using namespace brille::lattice;
using namespace brille::math;

TEST_CASE("LatticeDual convenience constructors","[latticedual]"){
  std::array<double, 3> len {two_pi, two_pi, two_pi}, ang {half_pi, half_pi, half_pi}, rlen {1, 1, 1};
  auto d = Direct(len, ang, "P 1");
  auto r = Reciprocal(rlen, ang, "P 1");
  REQUIRE(d == r);
}

/*
Add a test case to expose the bug associated with
      github.com/brille/brille/issues/28
which raises an exception due to the valid Hall symbol 'P -2zb' not having
an assigned Hall 'number' and `Lattice::inner_star` not written to handle
this valid possibility.
*/
TEST_CASE("LatticeDual inner_star", "[latticedual]"){
  double h = half_pi;
  std::string hall_symbol = "P -2zb";
  auto dlat = Direct<double>({4.0, 4.0, 4.0}, {h, h, h}, hall_symbol);
  auto rlat = Reciprocal<double>({h, h, h}, {h, h, h}, hall_symbol);
  REQUIRE(rlat == dlat);
}

TEST_CASE("LatticeDual with Seitz symbol","[latticedual]"){
  double h = half_pi;
  std::string seitz_symbol = "x, y+1/2, -z";
  auto dlat = Direct<double>({4.0, 4.0, 4.0}, {h, h, h}, seitz_symbol);
  auto rlat = Reciprocal<double>({h, h, h}, {h, h, h}, seitz_symbol);
  REQUIRE(rlat == dlat);
}
