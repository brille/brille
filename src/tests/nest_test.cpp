#include <catch2/catch_test_macros.hpp>

#include <tuple>
#include <catch2/matchers/catch_matchers_quantifiers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "bz_nest.hpp"

using namespace brille;
using namespace brille::math;
using namespace brille::lattice;

TEST_CASE("BrillouinZoneNest3 instantiation","[nest]"){
  // The conventional cell for Nb
  std::array<double,3> len{3.2598, 3.2598, 3.2598}, ang{half_pi, half_pi, half_pi};
  auto lat = Direct(len, ang, "-I 4 2 3"); // was 529
  BrillouinZone bz(lat);
  double max_volume = 0.01;
  ind_t number_rho = 100;
  ind_t max_branchings = 5;
  BrillouinZoneNest3<double,double,double> bzn0(bz, number_rho, max_branchings);
  BrillouinZoneNest3<double,double,double> bzn1(bz, max_volume, max_branchings);
}
TEST_CASE("BrillouinZoneNest3 vertex accessors","[nest]"){
  auto lat = Direct<double>({10.75, 10.75, 10.75}, {half_pi, half_pi, half_pi}, "F 4d 2 3 -1d"); // was 525
  BrillouinZone bz(lat);
  ind_t number_rho = 1000;
  ind_t max_branchings = 5;
  BrillouinZoneNest3<double,double,double> bzn(bz, number_rho, max_branchings);

  SECTION("get_xyz"){auto verts = bzn.get_xyz(); REQUIRE(verts.size(0) > 0u);}
  SECTION("get_hkl"){auto verts = bzn.get_hkl(); REQUIRE(verts.size(0) > 0u);}
  SECTION("get_all_xyz"){auto verts = bzn.get_all_xyz(); REQUIRE(verts.size(0) > 0u);}
  SECTION("get_all_hkl"){auto verts = bzn.get_all_hkl(); REQUIRE(verts.size(0) > 0u);}
}


/* This test is not very good since it attempts to verify that interpolating a vector-valued field with cubic symmetry
 * works. The problem being that no real vector can 'rotate like' a reciprocal lattice vector *and* follow the symmetry
 * of a cube. The operations need to map (x,y,z) to (-x, -y, z), (x, -y, -z), (-x, y, -z), (y, x, z), (y, z, x), etc.
 * But a finite vector field can not obey the vast majority of the symmetries.
 * This test can only ever work if the interpolation points stay within the irreducible Brillouin zone, which rounding
 * errors can easily disrupt; otherwise a non-identity rotation matrix due to non-zero lattice translation will cause
 * mixing of the components of the interpolated vector which will not work.
 * */
TEST_CASE("Simple BrillouinZoneNest3 interpolation","[nest][macos-arm]"){
  std::array<double,3> len{3.2598, 3.2598, 3.2598}, ang{half_pi, half_pi, half_pi};
  auto lat = Direct(len, ang, "-I 4 2 3"); // was 529
  BrillouinZone bz(lat);
  double max_volume = 0.01;
  ind_t max_branchings = 5;
  BrillouinZoneNest3<double,double,double> bzn(bz, max_volume, max_branchings);

  auto Qmap = bzn.get_hkl();
  auto Qxyz = bzn.get_xyz();
  // At present converting from Array2 to Array is not possible, we must copy data by hand
  brille::shape_t tostoreshape{Qmap.size(0), 1u, Qmap.size(1)};
  brille::Array<double> tostore(tostoreshape);
  for (auto i: tostore.subItr()) tostore[i] = Qxyz.val(i[0], i[2]);
  std::array<unsigned,3> elements{0,3,0};
  RotatesLike rt = RotatesLike::vector;
  LengthUnit lu = LengthUnit::reciprocal_lattice;
  bzn.replace_value_data( tostore, elements, rt, lu);

  // If branches ≠ 1 then the vector is not treated as a vector!
  REQUIRE(bzn.data().branches() == 1u);

  std::default_random_engine generator(static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count()));
  std::uniform_real_distribution<double> distribution(0.,1.);

  auto Q = LQVec<double>(lat);
  brille::ind_t nQ{0};

//  SECTION("Random points between mapped Q points") {
//    std::cout << "Random points between mapped Q points\n";
//    // construct random Q points that are linear combinations of the mapped nest Q points
//    brille::ind_t nQmap = Qmap.size(0);
//    nQ = 10;//10000;
//    Q.resize(nQ, 3u);
//    //Q = LQVec<double>(lat, nQ);
//    double rli;
//    for (ind_t i = 0; i < nQ; ++i) {
//      rli = distribution(generator);
//      Q.set(i, rli * Qmap.view(i % nQmap) + (1 - rli) * Qmap.view((i + 1) % nQmap));
//    }
//  }
//  SECTION("Check point(s) that failed on macOS-14 (arm64) before") {
//    // like all generated points, this one is on the surface of the ir_bz polyhedron
//    // but this one is not inside the ir_polyhedron according to ir_moveinto.
//    // _really_ strangely Q = inv(transpose(R)) q_ir + tau is such that Q == q_ir to within tolerance
//    nQ = 1u;
//    Q.resize(nQ, 3u);
//    // auto bad_Q = std::array<double, 3>{0.6362324202966432,  0.3637675797033568, -0.0000000000000011};
//    Q.set(0u, std::array<double, 3>{0.6362324202966432, 0.3637675797033568, -0.0000000000000011});
//  }
  SECTION("Random points inside the irreducible Brillouin zone polyhedron") {
    nQ = 10;
    Q = bz.get_ir_polyhedron().rand_rejection(nQ);
  }

  // verify that those points are all inside of the Brillouin zone
  REQUIRE_THAT(bz.isinside(Q), Catch::Matchers::AllTrue());
  // and in the irreducible wedge
  REQUIRE_THAT(bz.isinside_wedge(Q), Catch::Matchers::AllTrue());

  REQUIRE_THAT(bz.isinside_wedge_outer(Q, true), Catch::Matchers::AllTrue());

  auto ir_poly = bz.get_ir_polyhedron();
  auto in_ir_poly = ir_poly.contains(Q);
  REQUIRE_THAT(in_ir_poly, Catch::Matchers::AllTrue());

  auto first_poly = bz.get_polyhedron();
  auto in_first_poly = first_poly.contains(Q);
  REQUIRE_THAT(in_first_poly, Catch::Matchers::AllTrue());

  // Extra check that the points do not get moved by ir_moveinto (that they're not so close that rounding is a problem)
  auto q_ir = 0 * Q;
  auto tau = LQVec<int>(lat);
  auto rot = std::vector<size_t>(nQ, 0u);
  auto inv_rot = std::vector<size_t>(nQ, 0u);
  bz.ir_moveinto(Q, q_ir, tau, rot, inv_rot);
  REQUIRE(Q == q_ir); // because of how we constructed Q, this can be true even if tau != 0 which is surprising
  // if any tau is not 0, the rest of this test _will_ fail
  REQUIRE(tau == 0 * tau);
  // this should be an equivalent test to tau == 0
  REQUIRE_THAT(rot, Catch::Matchers::AllMatch(Catch::Matchers::WithinAbs(bz.get_pointgroup_symmetry().find_identity_index(), 0)));

  auto Q_invA = Q.xyz();
  brille::Array<double> expected(brille::shape_t{nQ, 1u, 3u});
  for (auto i: expected.subItr()) expected[i] = Q_invA.val(i[0], i[2]);

  auto check_interpolate = [&](const auto & result){
    auto diff = result - expected;
    if (!(diff.round().all(brille::cmp::eq, 0))) for (ind_t i = 0; i < nQ; ++i) {
        info_update_if(!(diff.view(i).round().all(0,0)),
                       "The interpolation point Q = ", Q.to_string(i), "\n",
                       "           returned result ", result.to_string(i), "\n",
                       "                instead of ", expected.to_string(i), "\n");
      }
    REQUIRE( diff.round().all(brille::cmp::eq, 0) ); // this is not a great test :(
    for (auto i: diff.valItr()) REQUIRE(std::abs(i) < 2E-14);
  };

  auto [no_move, unused] = bzn.ir_interpolate_at<true>(Q, 1);
  check_interpolate(no_move);

  //Now the real check, allow ir_interpolate_at to move Q into the irreducible wedge
  auto [move, also_unused] = bzn.ir_interpolate_at(Q, 1);
  check_interpolate(move);
}

TEST_CASE("Random BrillouinZoneNest3 interpolation","[nest][nb][random]"){
  // The conventional cell for Nb
  std::array<double,3> len{3.2598, 3.2598, 3.2598}, ang{half_pi, half_pi, half_pi};
  auto lat = Direct(len, ang, "-I 4 2 3"); // was 529

  BrillouinZone bz(lat);
  double max_volume = 0.01;
  BrillouinZoneNest3<double,double,double> bzn(bz, max_volume);

  auto Qmap = bzn.get_hkl();
  auto Qxyz = bzn.get_xyz();
  // At present converting from Array2 to Array is not possible, so the
  // original reshape will not work. Instead we must copy data by hand
  brille::shape_t tostoreshape{Qmap.size(0), 1u, Qmap.size(1)};
  brille::Array<double> tostore(tostoreshape);
  for (auto i: tostore.subItr()) tostore[i] = Qxyz.val(i[0], i[2]);
  std::array<unsigned,3> elements{0,3,0};
  RotatesLike rt = RotatesLike::vector;
  LengthUnit lu = LengthUnit::reciprocal_lattice;
  bzn.replace_value_data( tostore, elements, rt, lu);

  // If branches ≠ 1 then the vector is not treated as a vector!
  REQUIRE(bzn.data().branches() == 1u);

  brille::ind_t nQ = 100;
  // In order to have easily-interpretable results we need to ensure we only
  // interpolate at points within the irreducible meshed volume.
  // Use points distributed randomly in the Irreducible Polyhedron
  auto irp = bz.get_ir_polyhedron();

  auto interpolate_at_random_points = [&](unsigned int seed){
    auto Q = irp.rand_rejection(nQ, seed); // already reciprocal lattice vectors
    // Q are now random points in the irreducible Brillouin zone polyhedron

    auto Q_in_irp = irp.contains(Q);
    REQUIRE(std::find(Q_in_irp.begin(), Q_in_irp.end(), false) == Q_in_irp.end());

    // We may run into problems if any of the points are too close to the irBz
    // boundary where the components expressed in the primitive lattice are 0.5
    // Since different "std" library **versions** round 0.5 differently (I'm looking at you MSVC).
    // If 0.5 is rounded to 1 this simple test *will* fail since bz.moveinto rounds
    // the components of Q to try and find an equivalent q and tau.

    auto [intres, dummy] = bzn.ir_interpolate_at(Q, 1 /*thread*/);

    auto QinvA = Q.xyz();
    // QinvA is (at present) a brille::Array2<double> and so can not be reshaped
    // to 3D (maybe introducing conversion routines is a good idea?)
    // Instead make a new brille::Array and copy the Q values by hand:
    brille::shape_t antshp{nQ, 1u, 3u};
    brille::Array<double> antres(antshp);
    for (auto i: antres.subItr()) antres[i] = QinvA.val(i[0],i[2]);

    auto diff = intres - antres;
    // info_update("\nInterpolation results Expected results:\n",antres.to_string(intres));
    // info_update("\nRounded difference:\n",diff.to_string());

    if (!diff.round().all(0,0)) for (ind_t i = 0; i < nQ; ++i) {
        info_update_if(!diff.view(i).round().all(0,0),
                       "\nThe interpolation point Q = ", Q.to_string(i),
                       "\n            returned result ", intres.to_string(i),
                       "\n                 instead of ", antres.to_string(i), "\n");
      }
    REQUIRE( diff.round().all(0,0) ); // this is not a great test :(
    for (auto i: diff.valItr()) REQUIRE(std::abs(i) < 2E-10);
  };
  interpolate_at_random_points(20220718u); // seed for reproducible testing
  interpolate_at_random_points(0u); // 0 -> use time since epoch as seed
}



TEST_CASE("BrillouinZoneNest3 interpolation along face","[nest][nb][gb976690]"){
  // The conventional cell for Nb
  std::array<double,3> len{3.2598, 3.2598, 3.2598}, ang{half_pi, half_pi, half_pi};
  auto lat = Direct(len, ang, "-I 4 2 3"); // was 529

  BrillouinZone bz(lat);
  double max_volume = 0.01;
  BrillouinZoneNest3<double,double,double> bzn(bz, max_volume);

  auto Qmap = bzn.get_hkl();
  auto Qxyz = bzn.get_xyz();
  // At present converting from Array2 to Array is not possible, so the
  // original reshape will not work. Instead we must copy data by hand
  brille::shape_t tostoreshape{Qmap.size(0), 1u, Qmap.size(1)};
  brille::Array<double> tostore(tostoreshape);
  for (auto i: tostore.subItr()) tostore[i] = Qxyz.val(i[0], i[2]);
  std::array<unsigned,3> elements{0,3,0};
  RotatesLike rt = RotatesLike::vector;
  LengthUnit lu = LengthUnit::reciprocal_lattice;
  bzn.replace_value_data( tostore, elements, rt, lu);

  // If branches ≠ 1 then the vector is not treated as a vector!
  REQUIRE(bzn.data().branches() == 1u);

  brille::ind_t nQ = 100;
  // Replicate errant 'random' point selection algorithm:
  auto irp = bz.get_ir_polyhedron();

  auto v_min = irp.vertices().min(0);
  auto v_del = irp.vertices().max(0) - v_min;
  auto Q =  0 * irp.vertices().view(0);
  double step{0.64/(static_cast<double>(nQ-1))};
  // skip the first point for Apple clang++ to be happy :(
  auto actual_nQ = nQ - 1u;
  Q.resize(actual_nQ);
  for (brille::ind_t i=0; i<actual_nQ; ++i){
    Q.set(i, v_min + v_del * static_cast<double>(i+1) * step);
  }
  // Q are now points along one face of the irreducible Brillouin zone

  auto Q_in_irp = irp.contains(Q);
  REQUIRE(std::find(Q_in_irp.begin(), Q_in_irp.end(), false) == Q_in_irp.end());

  // We may run into problems if any of the points are too close to the irBz
  // boundary where the components expressed in the primitive lattice are 0.5
  // Since different "std" library **versions** round 0.5 differently (I'm looking at you MSVC).
  // If 0.5 is rounded to 1 this simple test *will* fail since bz.moveinto rounds
  // the components of Q to try and find an equivalent q and tau.

  auto [intres, dummy] = bzn.ir_interpolate_at(Q, 1 /*thread*/);

  auto QinvA = Q.xyz();
  // QinvA is (at present) a brille::Array2<double> and so can not be reshaped
  // to 3D (maybe introducing conversion routines is a good idea?)
  // Instead make a new brille::Array and copy the Q values by hand:
  brille::shape_t antshp{actual_nQ, 1u, 3u};
  brille::Array<double> antres(antshp);
  for (auto i: antres.subItr()) antres[i] = QinvA.val(i[0],i[2]);

  auto diff = intres - antres;
  // info_update("\nInterpolation results Expected results:\n",antres.to_string(intres));
  // info_update("\nRounded difference:\n",diff.to_string());

  if (!diff.round().all(0,0)) for (ind_t i = 0; i < nQ; ++i) {
      info_update_if(!diff.view(i).round().all(0,0),
                     "\nThe interpolation point Q = ", Q.to_string(i),
                     "\n            returned result ", intres.to_string(i),
                     "\n                 instead of ", antres.to_string(i), "\n");
    }
  REQUIRE( diff.round().all(0,0) ); // this is not a great test :(
  for (auto i: diff.valItr()) REQUIRE(std::abs(i) < 2E-10);
}
