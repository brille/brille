#include <catch.hpp>

#include "grid.h"

TEST_CASE("Testing MapGrid3 instantiation"){
  MapGrid3<int> a;
  REQUIRE(a.numel() == 0);
  size_t n[3] = {3,4,5};
  MapGrid3<double> b(n);
  for (int i=0; i<3; i++) REQUIRE( b.size(i) == n[i] );
  MapGrid3<double> c(b);
  for (int i=0; i<3; i++) REQUIRE( c.size(i) == n[i] );
}
