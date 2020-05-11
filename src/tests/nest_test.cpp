#include <catch2/catch.hpp>
#include <tuple>
#include "bz_nest.hpp"

TEST_CASE("BrillouinZoneNest3 instantiation","[nest]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, PI/2, PI/2, PI/2, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.01;
  size_t number_rho = 100;
  size_t max_branchings = 5;
  BrillouinZoneNest3<double,double> bzn0(bz, number_rho, max_branchings);
  BrillouinZoneNest3<double,double> bzn1(bz, max_volume, max_branchings);
}
TEST_CASE("BrillouinZoneNest3 vertex accessors","[nest]"){
  Direct d(10.75, 10.75, 10.75, PI/2, PI/2, PI/2, 525);
  BrillouinZone bz(d.star());
  size_t number_rho = 1000;
  size_t max_branchings = 5;
  BrillouinZoneNest3<double,double> bzn(bz, number_rho, max_branchings);

  SECTION("get_xyz"){auto verts = bzn.get_xyz(); REQUIRE(verts.size() > 0u);}
  SECTION("get_hkl"){auto verts = bzn.get_hkl(); REQUIRE(verts.size() > 0u);}
  SECTION("get_all_xyz"){auto verts = bzn.get_all_xyz(); REQUIRE(verts.size() > 0u);}
  SECTION("get_all_hkl"){auto verts = bzn.get_all_hkl(); REQUIRE(verts.size() > 0u);}
}

TEST_CASE("Simple BrillouinZoneNest3 interpolation","[nest]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, PI/2, PI/2, PI/2, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.01;
  size_t max_branchings = 5;
  BrillouinZoneNest3<double,double> bzn(bz, max_volume, max_branchings);

  ArrayVector<double> Qmap = bzn.get_hkl();
  std::vector<size_t> shape{Qmap.size(), 1, 3};
  std::array<size_t,3> elements{0,3,0};
  RotatesLike rt = RotatesLike::Reciprocal;
  bzn.replace_value_data( bzn.get_xyz(), shape, elements, rt );

  // In order to have easily-interpretable results we need to ensure we only
  // interpolate at points within the irreducible meshed volume.
  // So let's stick to points that are random linear interpolations between
  // neighbouring mesh vertices
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(0.,1.);

  size_t nQmap = Qmap.size(), nQ = 10;//10000;
  LQVec<double> Q(r,nQ);
  double rli;
  for (size_t i=0; i<nQ; ++i){
    rli = distribution(generator);
    Q.set(i, rli*Qmap.extract(i%nQmap) + (1-rli)*Qmap.extract((i+1)%nQmap) );
  }

  ArrayVector<double> intres, dummy;
  std::tie(intres, dummy) = bzn.ir_interpolate_at(Q,1);
  ArrayVector<double> antres = Q.get_xyz();

  ArrayVector<double> diff = intres - antres;
  // printf("\nInterpolation results:\n");
  // intres.print();
  // printf("\nExpected results:\n");
  // antres.print();
  // printf("\nRounded difference:\n");
  // diff.round().print();

  REQUIRE( diff.round().all_zero() ); // this is not a great test :(
  for (size_t i=0; i<diff.size(); ++i)
  for (size_t j=0; j<diff.numel(); ++j)
  REQUIRE( abs(diff.getvalue(i,j))< 2E-14 );
}
