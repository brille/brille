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
  std::vector<size_t> shape{Qmap.size(), 3}; // was {Qmap.size(), 1, 3}
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

  if (!diff.round().all_zero()) for (size_t i = 0; i < nQ; ++i) {
      info_update_if(!diff.extract(i).round().all_zero(),
        "\nThe interpolation point Q = ", Q.to_string(i),
        "\n            returned result ", intres.to_string(i),
        "\n                 instead of ", antres.to_string(i), "\n");
  }
  REQUIRE( diff.round().all_zero() ); // this is not a great test :(
  for (size_t i=0; i<diff.size(); ++i)
  for (size_t j=0; j<diff.numel(); ++j)
  REQUIRE( abs(diff.getvalue(i,j))< 2E-14 );
}

TEST_CASE("Random BrillouinZoneNest3 interpolation","[nest]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, PI/2, PI/2, PI/2, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.01;
  BrillouinZoneNest3<double,double> bzn(bz, max_volume);

  ArrayVector<double> Qmap = bzn.get_hkl();
  std::vector<size_t> shape{Qmap.size(), 3};
  std::array<unsigned long,3> elements{0,3,0};
  RotatesLike rl = RotatesLike::Reciprocal;
  bzn.replace_value_data( bzn.get_xyz(), shape, elements, rl);

  size_t nQ = 10;
  // In order to have easily-interpretable results we need to ensure we only
  // interpolate at points within the irreducible meshed volume.
  // Use points distributed randomly in the Irreducible Polyhedron
  Polyhedron irp = bz.get_ir_polyhedron();
  ArrayVector<double> Qxyz = irp.rand_rejection(nQ);
  LQVec<double> Q(r, nQ);
  double fromxyz[9];
  r.get_inverse_xyz_transform(fromxyz);
  for (size_t i=0; i<nQ; ++i)
    multiply_matrix_vector<double,double,double,3>(Q.data(i), fromxyz, Qxyz.data(i));
  // Q are now random points in the irreducible Brillouin zone polyhedron

  // We may run into problems if any of the points are too close to the irBz
  // boundary where the components expressed in the primitive lattice are 0.5
  // Since different "std" library **versions** round 0.5 differently (I'm looking at you MSVC).
  // If 0.5 is rounded to 1 this simple test *will* fail since bz.moveinto rounds
  // the components of Q to try and find an equivalent q and tau.


  ArrayVector<double> intres, dummy, antres=Q.get_xyz();
  std::tie(intres, dummy) = bzn.ir_interpolate_at(Q, 1 /*thread*/);

  ArrayVector<double> diff = intres - antres;
  // info_update("\nInterpolation results Expected results:\n",antres.to_string(intres));
  // info_update("\nRounded difference:\n",diff.to_string());

  if (!diff.round().all_zero()) for (size_t i = 0; i < nQ; ++i) {
      info_update_if(!diff.extract(i).round().all_zero(),
        "\nThe interpolation point Q = ", Q.to_string(i),
        "\n            returned result ", intres.to_string(i),
        "\n                 instead of ", antres.to_string(i), "\n");
  }
  REQUIRE( diff.round().all_zero() ); // this is not a great test :(
  for (size_t i=0; i<diff.size(); ++i)
  for (size_t j=0; j<diff.numel(); ++j)
  REQUIRE( abs(diff.getvalue(i,j))< 2E-10 );
}
