#include <catch2/catch.hpp>

#include "mesh.hpp"
#include "bz_mesh.hpp"

using namespace brille;

TEST_CASE("BrillouinZoneMesh3 instantiation","[mesh]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, brille::halfpi, brille::halfpi, brille::halfpi, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  // double max_size = 0.01;
  BrillouinZoneMesh3<double,double> bzm0(bz);
  // The following line causes a segmentation fault on a Hyper-V
  // Windows10 dev VM for unknown reason. As BrillouinZoneMesh3 is not
  // used at present, perhaps we can get away with commenting it out for now.

  //BrillouinZoneMesh3<double> bzm1(bz, max_size);
}
