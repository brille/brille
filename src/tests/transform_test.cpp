#include <random>
#include <chrono>
#include <catch2/catch.hpp>
#include <array>

#include "spg_database.hpp"
#include "utilities.hpp"
#include "lattice_dual.hpp"
#include "array_.hpp"
#include "transform.hpp"
#include "approx.hpp"

using namespace brille;
using namespace brille::lattice;
using namespace brille::math;

TEST_CASE("primitive transforms","[transform]"){
  Bravais c{Bravais::_};
  SECTION("Body centring"){         c = Bravais::I; }
  SECTION("All-face centring"){     c = Bravais::F; }
  SECTION("A-face centring"){       c = Bravais::A; }
  SECTION("B-face centring"){       c = Bravais::B; }
  SECTION("C-face centring"){       c = Bravais::C; }
  SECTION("Rhombohedral centring"){ c = Bravais::R; }
  PrimitiveTransform PT(c);
  REQUIRE( !PT.does_nothing() );
  std::array<PrimitiveTraits::sixP,9> sixP = PT.get_6P();
  std::array<PrimitiveTraits::invP,9> invP = PT.get_invP();
  int I[9]={6,0,0, 0,6,0, 0,0,6};
  int res[9];
  brille::utils::multiply_matrix_matrix<int,PrimitiveTraits::invP,PrimitiveTraits::sixP,3>(res,invP.data(),sixP.data());
  REQUIRE( brille::approx::matrix<int, int, 3>(res,I) );
}

TEST_CASE("primitive vector transforms","[transform]"){
  std::string spgr;
  SECTION("Primitive spacegroup"){    spgr = "P 1";   }
  SECTION("Body-centred spacegroup"){ spgr = "Im-3m"; }
  SECTION("Face-centred spacegroup"){ spgr = "Fmm2";  }
  auto lat = Direct<double>({1,1,1}, {half_pi, half_pi, half_pi}, spgr);

  std::default_random_engine gen(static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count()));
  std::uniform_real_distribution<double> dst(-5.0,5.0);

  int nQ = 33;
  std::vector<std::array<double,3>> rawQ;
  for (int i=0; i<nQ; ++i) rawQ.push_back({dst(gen), dst(gen), dst(gen)});
  auto V = LDVec<double>(lat, bArray<double>::from_std(rawQ));
  auto Q = LQVec<double>(lat, bArray<double>::from_std(rawQ));

  auto Vp = transform_to_primitive(lat, V);

  // Test 1: make sure that |Vᵢ|==|Vpᵢ| and |Qᵢ|==|Qpᵢ|
  for (int i=0; i<nQ; ++i)
    REQUIRE( Vp.norm(i) == Approx(V.norm(i)) );
  // Test 2: Check compoments, expressed in inverse Angstrom:
  auto Vxyz = V.xyz();
  auto Vpxyz = Vp.xyz();
  for (int i=0; i<nQ; ++i)
    REQUIRE( Vxyz.norm(i) == Approx(Vpxyz.norm(i)) );
  // // THE FOLLOWING WILL FAIL, since the xyz coordinate system is arbitrarily
  // // aligned with x along a, and the direction of a and a' are not guaranteed
  // // to be the same!
  // for (auto i: SubIt(V.shape())) REQUIRE(Vpxyz[i] == Approx(Vxyz[i]));

  // Test 3: check transfrom_from_primitive(transform_to_primitive(X)) == X
  auto pVp = transform_from_primitive(lat,Vp);
  for (auto i: V.subItr()) REQUIRE(pVp[i] == Approx(V[i]));

  auto Qp = transform_to_primitive(lat, Q);
  // Test 1: make sure that |Qᵢ|==|Qpᵢ|
  for (int i=0; i<nQ; ++i)
    REQUIRE( Qp.norm(i) == Approx(Q.norm(i)) );
  // Test 2: Check compoments, expressed in inverse Angstrom:
  bArray<double> Qxyz = Q.xyz();
  bArray<double> Qpxyz = Qp.xyz();
  for (int i=0; i<nQ; ++i)
    REQUIRE( Qxyz.norm(i) == Approx(Qpxyz.norm(i)) );
  // // THE FOLLOWING WILL FAIL, since the xyz coordinate system is arbitrarily
  // // aligned with x along a, and the direction of a and a' are not guaranteed
  // // to be the same!
  // for (auto i: SubIt(Q.shape())) REQUIRE(Qpxyz[i] == Approx(Qxyz[i]));

  // Test 3: check transfrom_from_primitive(transform_to_primitive(X)) == X
  auto pQp = transform_from_primitive(lat, Qp);
  for (auto i: Q.subItr()) REQUIRE(pQp[i] == Approx(Q[i]));
}
