#include <catch.hpp>

#include "grid.hpp"
#include "bz_grid.hpp"
#include "interpolation.hpp"

TEST_CASE("MapGrid3 instantiation","[grid]"){
  MapGrid3<int> a;
  REQUIRE(a.numel() == 0);
  size_t n[3] = {3,4,5};
  MapGrid3<double> b(n);
  for (int i=0; i<3; i++) REQUIRE( b.size(i) == n[i] );
  MapGrid3<double> c(b);
  for (int i=0; i<3; i++) REQUIRE( c.size(i) == n[i] );
}
TEST_CASE("BrillouinZoneGrid3 instantiation","[grid]"){
  Direct d(3.,3.,3.,PI/2,PI/2,PI/2*3);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  size_t half[3]={2,2,2};
  BrillouinZoneGrid3<double> bzg(bz,half);
}

TEST_CASE("BrillouinZoneGrid4 instantiation","[grid]"){
  Direct d(3.,3.,3.,PI/2,PI/2,PI/2*3);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double spec[3]={0.,1.,10.};
  size_t half[3]={2,2,2};
  BrillouinZoneGrid4<double> bzg(bz,spec,half);

}
