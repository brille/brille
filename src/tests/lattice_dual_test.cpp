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



TEST_CASE("La2Zr2O7 construction off-symmetry basis vector input","[lattice][la2zr2o7][64]"){
  // Basis vectors from https://github.com/pace-neutrons/Euphonic/blob/aa3cc28786797bb3052f898dd63d4928d6f27ee2/tests_and_analysis/test/data/force_constants/LZO_force_constants.json#L186
  std::array<double,9> latmat {7.583912824349999, 1.8412792137035698e-32, 0.,
                               3.791956412170034, 3.791956412170034, 5.362636186024768,
                               3.791956412170034,-3.791956412170034, 5.362636186024768};
  // Symmetry information from CASTEP file via brilleu
  // the generators are: 4-fold [1 -1 1], 2-fold [-1 1 1], 3-fold [1 1 -3], -𝟙
  // row-ordered generator matrices
  //  std::vector<std::array<int,9>> W {
  //      {{ 0,-1, 0,  0, 0,-1,  1, 1, 1}},
  //      {{-1,-1,-1,  0, 0, 1,  0, 1, 0}},
  //      {{-1,-1,-1,  1, 0, 0,  0, 0, 1}},
  //      {{-1, 0, 0,  0,-1, 0,  0, 0,-1}}
  //  };
  //  std::vector<std::array<double,3>> w{
  //      {{0.0, 0.0, 0.5}},
  //      {{0.5, 0.0, 0.0}},
  //      {{0.5, 0.0, 0.0}},
  //      {{0.0, 0.0, 0.0}}
  //  };

  std::vector<std::array<int,9>> W {
      {{ 0,-1, 0,  1, 1, 1, -1, 0, 0}},
      {{ 0, 0, 1, -1,-1,-1,  1, 0, 0}},
      {{ 0, 0, 1,  1, 0, 0,  0, 1, 0}},
      {{-1, 0, 0,  0,-1, 0,  0, 0,-1}}
  };
  std::vector<std::array<double,3>> w{
      {{0.0, 0.5, 0.0}},
      {{0.0, 0.5, 0.0}},
      {{0.0, 0.0, 0.0}},
      {{0.0, 0.0, 0.0}}
  };

  Symmetry::Motions mots;
  mots.reserve(W.size());
  for (size_t i=0; i<W.size(); ++i) mots.push_back(Motion<int,double>(W[i], w[i]));
  Symmetry sym(mots);
  //
  LengthUnit lu { LengthUnit::angstrom };
  AngleUnit au {AngleUnit::degree };
  //
  auto wrong_lat = Direct<double>(latmat, MatrixVectors::row, sym, Basis(), false);
  auto wrong_lengths = wrong_lat.lengths(lu);
  auto wrong_angles = wrong_lat.angles(lu, au);

  auto right_lat = Direct<double>(latmat, MatrixVectors::row, sym, Basis(), true);
  auto right_lengths = right_lat.lengths(lu);
  auto right_angles = right_lat.angles(lu, au);

  auto std_dev = [](const auto & v){
    double m{0};
    for (const auto & x: v) m += x;
    m /= static_cast<double>(v.size());
    double s{0};
    for (const auto & x: v) s += (x-m)*(x-m);
    s = std::sqrt(s) / static_cast<double>(v.size() - 1);
    return s;
  };
  // snapping should improve the standard deviation in the parameters
  REQUIRE(std_dev(right_lengths) <= std_dev(wrong_lengths));
  REQUIRE(std_dev(right_angles) <= std_dev(wrong_angles));
	//
	auto rlu = LengthUnit::inverse_angstrom;
  auto wrong_mat = wrong_lat.to_xyz(rlu);
  auto right_mat = right_lat.to_xyz(rlu);
  for (size_t i=0; i<9u; ++i){
    // snapping to symmetry should preserve the basis vector orientation
		REQUIRE(std::abs(wrong_mat[i] - right_mat[i]) < 1e-12);
  }
}

TEST_CASE("Basis vectors corrected","[limits][snap_to_symmetry]"){
  auto test = [](double eps){
    double a{3.95}, c{12.9};
    double two_pi_over_volume{math::two_pi / (a * a * c)};
    std::array<double, 3> exact{a * c, a * c, a * a}, len{a+eps, a-eps, c}, ang{90.,90.,90.};
    for (auto & x: exact) x *= two_pi_over_volume;
    // info_update("len=",len);
    auto lat = Direct<double>(len, ang, "I4/mmm", Basis(), true);
    // info_update(lat.to_string(LengthUnit::angstrom));
    REQUIRE(approx_float::vector(exact.data(), lat.lengths(LengthUnit::inverse_angstrom).data()));
    REQUIRE(approx_float::vector(ang.data(), lat.angles(LengthUnit::inverse_angstrom, AngleUnit::degree).data()));
  };

  for (int p=-17; p<-9; ++p){
    for (int x=1; x<10; x+=2){
      auto eps = static_cast<double>(x) * std::pow(10.0, static_cast<double>(p));
      // std::cout << "p = " << p << ", eps = " << eps << std::endl;
      test(eps);
    }
  }


}

TEST_CASE("snap_to_symmetry should not change parameters unnecessarily", "[lattice][snap_to_symmetry]"){
    std::array<double,3> len{3.2598, 3.2598, 3.2598}, ang{half_pi, half_pi, half_pi};
    auto n = Direct(len, ang, "-I 4 2 3", Basis(), false);
    auto w = Direct(len, ang, "-I 4 2 3", Basis(), true);

    auto diff = [](const std::array<double, 9>& a, const std::array<double, 9>& b){
        std::array<double,9> d{};
        for (size_t i=0; i<9u; ++i) d[i] = a[i] - b[i];
        return d;
    };
    auto sum = [](const std::array<double, 9> & d){
        double s{0};
        for (const auto & i: d) s += std::abs(i);
        return s;
    };

    REQUIRE(n == w);
    REQUIRE(sum(diff(n.metric(LengthUnit::angstrom), w.metric(LengthUnit::angstrom))) == 0.);
}