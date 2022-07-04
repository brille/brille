#include <random>
#include <chrono>
#include <catch2/catch.hpp>
#include <array>

#include "spg_database.hpp"
#include "utilities.hpp"
#include "lattice_dual.hpp"
#include "array_.hpp"
#include "transform.hpp"
#include "approx_float.hpp"

using namespace brille;
using namespace brille::lattice;
using namespace brille::math;
using namespace brille::approx_float;

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
  REQUIRE( matrix<int, int, 3>(res,I) );
}

TEST_CASE("primitive vector transforms","[transform]"){
  auto run_tests = [] (const std::string& spgr, const LengthUnit lu){
    auto lat = Direct<double>({1,1,1}, {half_pi, half_pi, half_pi}, spgr);

    std::default_random_engine gen(static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<double> dst(-5.0,5.0);

    int nQ = 500;
    std::vector<std::array<double,3>> rawQ;
    rawQ.reserve(nQ);
    for (int i=0; i<nQ; ++i) rawQ.push_back({dst(gen), dst(gen), dst(gen)});
//    for (int i=0; i<nQ; ++i) rawQ.push_back({static_cast<double>(i), static_cast<double>(i*i), static_cast<double>(i*i*i)});

    auto V = LVec<double>(lu, lat, bArray<double>::from_std(rawQ));
    auto Vp = transform_to_primitive(lat, V);
    // Test 1: make sure that |Vᵢ|==|Vpᵢ|
    for (int i=0; i<nQ; ++i)
      REQUIRE( Vp.norm(i) == Approx(V.norm(i)) );
    // Test 2: Check components, expressed in absolute LengthUnit units
    auto Vxyz = V.xyz();
    auto Vpxyz = Vp.xyz();
    for (int i=0; i<nQ; ++i) REQUIRE( Vxyz.norm(i) == Approx(Vpxyz.norm(i)) );
    // Test 2a: Verify that the individual components are the same in the cartesian coordinate system
    for (auto i: V.subItr()) REQUIRE(Vpxyz[i] == Approx(Vxyz[i]));

    // Test 3: check transfrom_from_primitive(transform_to_primitive(X)) == X
    auto pVp = transform_from_primitive(lat,Vp);
    for (auto i: V.subItr()) REQUIRE(pVp[i] == Approx(V[i]));
  };
  SECTION("Primitive, direct"){    run_tests("P 1", LengthUnit::angstrom);   }
  SECTION("Body-centred, direct"){ run_tests("Im-3m", LengthUnit::angstrom); }
  SECTION("Face-centred, direct"){ run_tests("Fmm2", LengthUnit::angstrom);  }
  SECTION("Primitive, reciprocal"){    run_tests("P 1", LengthUnit::inverse_angstrom);   }
  SECTION("Body-centred, reciprocal"){ run_tests("Im-3m", LengthUnit::inverse_angstrom); }
  SECTION("Face-centred, reciprocal"){ run_tests("Fmm2", LengthUnit::inverse_angstrom);  }
}
