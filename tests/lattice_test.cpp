// #include "catch.hpp"
#include <catch.hpp>

#include "lattice.h"

TEST_CASE("Lattice and its subclasses are tested","[lattice]"){
  Lattice l(2*PI,2*PI,2*PI,PI/2,PI/2,PI/2);
  Direct d(l);
  Reciprocal r(1,1,1,PI/2,PI/2,PI/2);
  REQUIRE(d.issame(r.star()));
}
