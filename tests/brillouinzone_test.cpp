#include <catch.hpp>
#include "bz.h"

TEST_CASE("Testing BrillouinZone instantiation"){
  Direct d(3.,3.,3.,PI/2,PI/2,PI/2*3);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
}
