#include <catch.hpp>

#include "bz_trellis.h"

TEST_CASE("BrillouinZoneTrellis3 instantiation","[trellis]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, PI/2, PI/2, PI/2, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double node_volume_fraction = 0.01;
  BrillouinZoneTrellis3<double> bzt0(bz);
  BrillouinZoneTrellis3<double> bzt1(bz, node_volume_fraction);
}
TEST_CASE("BrillouinZoneTrellis3 vertex accessors","[trellis]"){
  Direct d(10.75, 10.75, 10.75, PI/2, PI/2, PI/2, 525);
  BrillouinZone bz(d.star());
  double node_volume_fraction = 0.002;
  BrillouinZoneTrellis3<double> bzt(bz, node_volume_fraction);

  SECTION("get_xyz"){auto verts = bzt.get_xyz(); REQUIRE(verts.size() > 0u);}
  SECTION("get_hkl"){auto verts = bzt.get_hkl(); REQUIRE(verts.size() > 0u);}
  SECTION("get_outer_xyz"){auto verts = bzt.get_outer_xyz(); REQUIRE(verts.size() > 0u);}
  SECTION("get_outer_hkl"){auto verts = bzt.get_outer_hkl(); REQUIRE(verts.size() > 0u);}
  SECTION("get_inner_xyz"){auto verts = bzt.get_inner_xyz(); REQUIRE(verts.size() > 0u);}
  SECTION("get_inner_hkl"){auto verts = bzt.get_inner_hkl(); REQUIRE(verts.size() > 0u);}
}

// TEST_CASE("Simple BrillouinZoneMesh3 interpolation","[mesh]"){
//   // The conventional cell for Nb
//   Direct d(3.2598, 3.2598, 3.2598, PI/2, PI/2, PI/2, 529);
//   Reciprocal r = d.star();
//   BrillouinZone bz(r);
//   double max_size = 0.01;
//   BrillouinZoneMesh3<double> bzm(bz, max_size);
//
//   ArrayVector<double> Qmap = bzm.get_mesh_hkl();
//   bzm.replace_data( Qmap.get_xyz() );
//
//   // In order to have easily-interpretable results we need to ensure we only
//   // interpolate at points within the irreducible meshed volume.
//   // So let's stick to points that are random linear interpolations between
//   // neighbouring mesh vertices
//   std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
//   std::uniform_real_distribution<double> distribution(0.,1.);
//
//   size_t nQmap = Qmap.size(), nQ = 10;//10000;
//   LQVec<double> Q(r,nQ);
//   double rli;
//   for (size_t i=0; i<nQ; ++i){
//     rli = distribution(generator);
//     Q.set(i, rli*Qmap.extract(i%nQmap) + (1-rli)*Qmap.extract((i+1)%nQmap) );
//   }
//
//   ArrayVector<double> intres = bzm.interpolate_at(Q);
//   ArrayVector<double> antres = Q.get_xyz();
//
//   ArrayVector<double> diff = intres - antres;
//   // printf("\nInterpolation results:\n");
//   // intres.print();
//   // printf("\nExpected results:\n");
//   // antres.print();
//   // printf("\nRounded difference:\n");
//   // diff.round().print();
//
//   REQUIRE( diff.round().all_zero() ); // this is not a great test :(
//   for (size_t i=0; i<diff.size(); ++i)
//   for (size_t j=0; j<diff.numel(); ++j)
//   REQUIRE( abs(diff.getvalue(i,j))< 2E-14 );
// }
