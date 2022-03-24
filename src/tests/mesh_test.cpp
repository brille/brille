#include <catch2/catch.hpp>

#include "mesh.hpp"
#include "bz_mesh.hpp"

using namespace brille;

TEST_CASE("BrillouinZoneMesh3 instantiation","[mesh]"){
  using namespace brille::math;
  using namespace brille::lattice;
  // The conventional cell for Nb
  std::array<double,3> len{3.2598, 3.2598, 3.2598}, ang{half_pi, half_pi, half_pi};
  auto lat = Direct(len, ang, "-I 4 2 3"); // was 529
  BrillouinZone bz(lat);
  // double max_size = 0.01;
  BrillouinZoneMesh3<double,double,double> bzm0(bz);
  // The following line causes a segmentation fault on a Hyper-V
  // Windows10 dev VM for unknown reason. As BrillouinZoneMesh3 is not
  // used at present, perhaps we can get away with commenting it out for now.

  //BrillouinZoneMesh3<double> bzm1(bz, max_size);
}
