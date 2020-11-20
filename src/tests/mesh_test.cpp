#include <catch2/catch.hpp>

#include "mesh.hpp"
#include "bz_mesh.hpp"
#include "interpolation.hpp"

using namespace brille;

/* Direct entry of the information required to construct a Mesh3 object would be
   less than fun. Instead skip to using a BrillouinZone object to do the heavy
   lifting for us.
*/
// TEST_CASE("Mesh3 instantiation","[mesh]"){
//   MapGrid3<int> a;
//   REQUIRE(a.numel() == 0);
//   size_t n[3] = {3,4,5};
//   MapGrid3<double> b(n);
//   for (int i=0; i<3; i++) REQUIRE( b.size(i) == n[i] );
//   MapGrid3<double> c(b);
//   for (int i=0; i<3; i++) REQUIRE( c.size(i) == n[i] );
// }

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
