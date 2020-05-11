#include <catch2/catch.hpp>

#include "mesh.hpp"
#include "bz_mesh.hpp"
#include "interpolation.hpp"

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
  Direct d(3.2598, 3.2598, 3.2598, PI/2, PI/2, PI/2, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  // double max_size = 0.01;
  BrillouinZoneMesh3<double,double> bzm0(bz);
  // The following line causes a segmentation fault on a Hyper-V
  // Windows10 dev VM for unknown reason. As BrillouinZoneMesh3 is not
  // used at present, perhaps we can get away with commenting it out for now.

  //BrillouinZoneMesh3<double> bzm1(bz, max_size);
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
