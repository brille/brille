#include <catch2/catch.hpp>
#include <tuple>
#include "bz_nest.hpp"

TEST_CASE("BrillouinZoneNest3 instantiation","[nest]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, brille::halfpi, brille::halfpi, brille::halfpi, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.01;
  size_t number_rho = 100;
  size_t max_branchings = 5;
  BrillouinZoneNest3<double,double> bzn0(bz, number_rho, max_branchings);
  BrillouinZoneNest3<double,double> bzn1(bz, max_volume, max_branchings);
}
TEST_CASE("BrillouinZoneNest3 vertex accessors","[nest]"){
  Direct d(10.75, 10.75, 10.75, brille::halfpi, brille::halfpi, brille::halfpi, 525);
  BrillouinZone bz(d.star());
  size_t number_rho = 1000;
  size_t max_branchings = 5;
  BrillouinZoneNest3<double,double> bzn(bz, number_rho, max_branchings);

  SECTION("get_xyz"){auto verts = bzn.get_xyz(); REQUIRE(verts.size(0) > 0u);}
  SECTION("get_hkl"){auto verts = bzn.get_hkl(); REQUIRE(verts.size(0) > 0u);}
  SECTION("get_all_xyz"){auto verts = bzn.get_all_xyz(); REQUIRE(verts.size(0) > 0u);}
  SECTION("get_all_hkl"){auto verts = bzn.get_all_hkl(); REQUIRE(verts.size(0) > 0u);}
}

TEST_CASE("Simple BrillouinZoneNest3 interpolation","[nest]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, brille::halfpi, brille::halfpi, brille::halfpi, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.01;
  size_t max_branchings = 5;
  BrillouinZoneNest3<double,double> bzn(bz, max_volume, max_branchings);

  auto Qmap = bzn.get_hkl();
  auto tostore = bzn.get_xyz().reshape({Qmap.size(0),1u,Qmap.size(1)});
  std::array<unsigned,3> elements{0,3,0};
  RotatesLike rt = RotatesLike::Reciprocal;
  bzn.replace_value_data( tostore, elements, rt );

  // If branches ≠ 1 then the vector is not treated as a vector!
  REQUIRE(bzn.data().branches() == 1u);

  // In order to have easily-interpretable results we need to ensure we only
  // interpolate at points within the irreducible meshed volume.
  // So let's stick to points that are random linear interpolations between
  // neighbouring mesh vertices
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(0.,1.);

  brille::ind_t nQmap = Qmap.size(0), nQ = 10;//10000;
  LQVec<double> Q(r,nQ);
  double rli;
  for (size_t i=0; i<nQ; ++i){
    rli = distribution(generator);
    Q.set(i, rli*Qmap.view(i%nQmap) + (1-rli)*Qmap.view((i+1)%nQmap) );
  }

  auto [intres, dummy] =  bzn.ir_interpolate_at(Q,1);
  brille::shape_t antshp{nQ, 1u, 3u};
  brille::Array<double> QinvA = Q.get_xyz();
  brille::Array<double> antres = QinvA.reshape(antshp);

  brille::Array<double> diff = intres - antres;

  if (!diff.round().all(brille::cmp::eq, 0.)) for (size_t i = 0; i < nQ; ++i) {
      info_update_if(!diff.view(i).round().all(0.,0),
        "The interpolation point Q = ", Q.to_string(i),
        "            returned result ", intres.to_string(i),
        "                 instead of ", antres.to_string(i));
  }
  REQUIRE( diff.round().all(brille::cmp::eq, 0.) ); // this is not a great test :(
  for (auto i: brille::ArrayIt(diff))
    REQUIRE(std::abs(i) < 2E-14);
}

TEST_CASE("Random BrillouinZoneNest3 interpolation","[nest]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, brille::halfpi, brille::halfpi, brille::halfpi, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.01;
  BrillouinZoneNest3<double,double> bzn(bz, max_volume);

  auto Qmap = bzn.get_hkl();
  auto tostore = bzn.get_xyz().reshape({Qmap.size(0),1u,Qmap.size(1)});
  std::array<unsigned,3> elements{0,3,0};
  RotatesLike rt = RotatesLike::Reciprocal;
  bzn.replace_value_data( tostore, elements, rt );

  // If branches ≠ 1 then the vector is not treated as a vector!
  REQUIRE(bzn.data().branches() == 1u);

  brille::ind_t nQ = 10;
  // In order to have easily-interpretable results we need to ensure we only
  // interpolate at points within the irreducible meshed volume.
  // Use points distributed randomly in the Irreducible Polyhedron
  Polyhedron irp = bz.get_ir_polyhedron();
  auto Qxyz = irp.rand_rejection(nQ);
  auto Q = LQVec<double>::from_invA(r, Qxyz);
  // Q are now random points in the irreducible Brillouin zone polyhedron

  // We may run into problems if any of the points are too close to the irBz
  // boundary where the components expressed in the primitive lattice are 0.5
  // Since different "std" library **versions** round 0.5 differently (I'm looking at you MSVC).
  // If 0.5 is rounded to 1 this simple test *will* fail since bz.moveinto rounds
  // the components of Q to try and find an equivalent q and tau.


  brille::Array<double> antres=Q.get_xyz().reshape({nQ,1u,3u});
  auto [intres, dummy] = bzn.ir_interpolate_at(Q, 1 /*thread*/);

  brille::Array<double> diff = intres - antres;
  // info_update("\nInterpolation results Expected results:\n",antres.to_string(intres));
  // info_update("\nRounded difference:\n",diff.to_string());

  if (!diff.round().all(0.,0)) for (size_t i = 0; i < nQ; ++i) {
      info_update_if(!diff.view(i).round().all(0.,0),
        "\nThe interpolation point Q = ", Q.to_string(i),
        "\n            returned result ", intres.to_string(i),
        "\n                 instead of ", antres.to_string(i), "\n");
  }
  REQUIRE( diff.round().all(0.,0.) ); // this is not a great test :(
  for (auto i: brille::ArrayIt(diff))
    REQUIRE(std::abs(i) < 2E-10);
}
