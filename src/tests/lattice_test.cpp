#include <catch2/catch.hpp>

#include "lattice.hpp"

TEST_CASE("Lattice and its subclasses are tested","[lattice]"){
  Lattice l(2*PI,2*PI,2*PI,PI/2,PI/2,PI/2);
  Direct d(l);
  Reciprocal r(1,1,1,PI/2,PI/2,PI/2);
  REQUIRE(d.issame(r.star()));
}

/*
Add a test case to expose the bug associated with
      github.com/g5t/brille/issues/28
which raises an exception due to the valid Hall symbol 'P -2zb' not having
an assigned Hall 'number' and `Lattice::inner_star` not written to handle
this valid possibility.
*/
TEST_CASE("Lattice inner_star", "[lattice]"){
  double h = PI/2;
  std::string hall_symbol = "P -2zb";
  Direct dlat(4.0, 4.0, 4.0, h, h, h, hall_symbol);
  Reciprocal rlat(h, h, h, h, h, h, hall_symbol);
  REQUIRE(rlat.issame(dlat.star()));
}
