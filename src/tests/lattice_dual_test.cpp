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
  std::string hall_symbol = "-R 3";
  // the other version below has slightly-off b vector length
  std::array<double,3> lens{10.699417459999999, 10.699417459999999, 4.528988750000000};
  std::array<double,3> angs{90,90,120};
  //
  auto lat = Direct<double>(lens, angs, hall_symbol);

  // we have constructed A and B *both* as column vector matrices,
  // so they are related by B = 2π (A⁻¹)ᵀ   !!Note the transpose!!
  auto A = lat.real_basis_vectors();
  auto B = lat.reciprocal_basis_vectors();
  Lattice<double>::matrix_t I{1,0,0, 0,1,0, 0,0,1};
  auto AB = mul_mat_mat(A,transpose(B));
  for (size_t i=0; i < I.size(); ++i) REQUIRE(I[i]*two_pi == Approx(AB[i]));

  // verify that we provide mutually inverse to and from cartesian coordinate matrices
  auto to_from_check = [&](LengthUnit lu){
    auto t = lat.to_xyz(lu);
    auto f = lat.from_xyz(lu);
    auto tf = mul_mat_mat(t, f);
    for (size_t i=0; i<I.size(); ++i) REQUIRE(I[i] == Approx(tf[i]));
  };
  to_from_check(LengthUnit::angstrom);
  to_from_check(LengthUnit::inverse_angstrom);

}