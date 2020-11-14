#include <catch2/catch.hpp>

#include "lattice.hpp"

using namespace brille;

TEST_CASE("Lattice and its subclasses are tested","[lattice]"){
  Lattice l(2*brille::pi,2*brille::pi,2*brille::pi,brille::halfpi,brille::halfpi,brille::halfpi);
  Direct d(l);
  Reciprocal r(1,1,1,brille::halfpi,brille::halfpi,brille::halfpi);
  REQUIRE(d.issame(r.star()));
}

/*
Add a test case to expose the bug associated with
      github.com/brille/brille/issues/28
which raises an exception due to the valid Hall symbol 'P -2zb' not having
an assigned Hall 'number' and `Lattice::inner_star` not written to handle
this valid possibility.
*/
TEST_CASE("Lattice inner_star", "[lattice]"){
  double h = brille::halfpi;
  std::string hall_symbol = "P -2zb";
  Direct dlat(4.0, 4.0, 4.0, h, h, h, hall_symbol);
  Reciprocal rlat(h, h, h, h, h, h, hall_symbol);
  REQUIRE(rlat.issame(dlat.star()));
}

TEST_CASE("Lattice with Seitz symbol","[lattice]"){
  double h = brille::halfpi;
  std::string seitz_symbol = "x, y+1/2, -z";
  Direct dlat(4., 4., 4., h, h, h, seitz_symbol);
  Reciprocal rlat(h,h,h, h,h,h, seitz_symbol);
  REQUIRE(rlat.issame(dlat.star()));
}
