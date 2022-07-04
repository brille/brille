#include <catch2/catch.hpp>

#include "lattice_dual.hpp"

using namespace brille::lattice;
using namespace brille::math;
using namespace brille::linear_algebra;
using namespace brille;

TEST_CASE("LatticeDual convenience constructors","[latticedual]"){
  std::array<double, 3> len {two_pi, two_pi, two_pi}, ang {half_pi, half_pi, half_pi}, rlen {1, 1, 1};
  auto d = Direct(len, ang, "P 1");
  auto r = Reciprocal(rlen, ang, "P 1");
  REQUIRE(d == r);
}

/*
Add a test case to expose the bug associated with
      github.com/brille/brille/issues/28
which raises an exception due to the valid Hall symbol 'P -2zb' not having
an assigned Hall 'number' and `Lattice::inner_star` not written to handle
this valid possibility.
*/
TEST_CASE("LatticeDual inner_star", "[latticedual]"){
  double h = half_pi;
  std::string hall_symbol = "P -2zb";
  auto dlat = Direct<double>({4.0, 4.0, 4.0}, {h, h, h}, hall_symbol);
  auto rlat = Reciprocal<double>({h, h, h}, {h, h, h}, hall_symbol);
  REQUIRE(rlat == dlat);
}

TEST_CASE("LatticeDual with Seitz symbol","[latticedual]"){
  double h = half_pi;
  std::string seitz_symbol = "x, y+1/2, -z";
  auto dlat = Direct<double>({4.0, 4.0, 4.0}, {h, h, h}, seitz_symbol);
  auto rlat = Reciprocal<double>({h, h, h}, {h, h, h}, seitz_symbol);
  REQUIRE(rlat == dlat);
}


TEST_CASE("Verify Lattice property correctness","[latticedual]"){
  auto perform_checks = [] (Lattice<double>& lat, bool B_is_triangular){
    //we have constructed A and B *both* as column vector matrices,
    // so they are related by B = 2π (A⁻¹)ᵀ   !!Note the transpose!!
    auto A = lat.real_basis_vectors();
    auto B = lat.reciprocal_basis_vectors();
    Lattice<double>::matrix_t I{1,0,0, 0,1,0, 0,0,1};
    auto AB = mul_mat_mat(A,transpose(B));
    for (size_t i=0; i < I.size(); ++i) REQUIRE(I[i]*two_pi == Approx(AB[i]));

    if (B_is_triangular){
//      info_update("A and B should be triangular!");
//      info_update("A: ", A);
//      info_update("B: ", B);
      // we constructed B internally such that it is an Upper Triangular matrix
      // it's elements below the diagonal *must* be zero!
      REQUIRE(B[3] == 0.);
      REQUIRE(B[6] == 0.);
      REQUIRE(B[7] == 0.);
      // And as the transpose of the inverse, A *must* be Lower Triangular
      REQUIRE(A[1] == 0.);
      REQUIRE(A[2] == 0.);
      REQUIRE(A[5] == 0.);
    }

    // to_xyz(angstrom) == A && to_xyz(inverse_angstrom) == B
    auto A1 = lat.to_xyz(LengthUnit::angstrom);
    auto B1 = lat.to_xyz(LengthUnit::inverse_angstrom);
    REQUIRE(approx_float::equal(A, A1));
    REQUIRE(approx_float::equal(B, B1));

    // verify that we provide mutually inverse to and from cartesian coordinate matrices
    auto to_from_check = [&](LengthUnit lu){
      auto t = lat.to_xyz(lu);
      auto f = lat.from_xyz(lu);
      auto tf = mul_mat_mat(t, f);
      for (size_t i=0; i<I.size(); ++i) REQUIRE(I[i] == Approx(tf[i]));
    };
    to_from_check(LengthUnit::angstrom);
    to_from_check(LengthUnit::inverse_angstrom);

    // the real space and reciprocal space metrics should be related
    auto G = lat.metric(LengthUnit::angstrom);
    auto invG = lat.metric(LengthUnit::inverse_angstrom);
    auto GinvG = mul_mat_mat(G, invG);
    for (auto & x: GinvG) x /= two_pi * two_pi;
    REQUIRE(approx_float::equal(I, GinvG)); // switch to brille::approx to overcome inaccuracies from inversion?
  };
  //
  SECTION("Rhombohedral"){
    std::string hall_symbol = "-R 3";
    // the other version below has slightly-off b vector length
    std::array<double,3> lens{10.699417459999999, 10.699417459999999, 4.528988750000000};
    std::array<double,3> angs{90,90,120};
    //
    auto lat = Direct<double>(lens, angs, hall_symbol);
    perform_checks(lat, true);
  }
  SECTION("Primitive Hexagonal"){
    auto lat = Direct<double>({3., 3., 3.}, {90., 90., 120.}, "P 1");
    perform_checks(lat, true);
  }
}
