#include <catch2/catch.hpp>

#include "lattice_dual.hpp"
#include "array_l_.hpp"

using namespace brille;
using namespace brille::lattice;

TEST_CASE("Lattice Vector tests","[latvec]"){
  auto lattice = Direct<double>({1,1,1}, {90,90,90}, "P 1");

  std::vector<std::array<double,3>> values{{1,0,0},{0,1,0},{0,0,1}};
  auto q = LQVec<double>(lattice, bArray<double>::from_std(values));

  SECTION("length via .norm(integer)"){
    for(int i=0; i<3; i++)
    REQUIRE( q.norm(i) == Approx(brille::math::two_pi) );
  }
  SECTION("length via norm()"){
    auto normq = norm(q);
    REQUIRE( normq.size(0) == q.size(0) );
    REQUIRE( normq.size(1) == 1u );
    for(size_t i=0; i<q.size(0); ++i)
      REQUIRE( q.norm(i) == normq[i] );
  }
  SECTION("Single-element access via [integer] and cross product .cross(integer,integer)"){
    auto z = q.view(2);
    auto xy = q.cross(0,1)*(norm(z)/norm(q.view(1))/norm(q.view(0)));
    REQUIRE( z.is(xy) );
  }
}

////////////////////////////////////////////////////////////////////////////////
//  TODO: Implement tests for:                                                //
//============================================================================//
// Class methods:                                                             //
//----------------------------------------------------------------------------//
// instantiation via L[DQ]Vec(::lattice,::Array)                        //
// copy via L[DQ]Vec(::L[DQ]Vec)                                              //
// assignment via =                                                           //
// .get_lattice()                                                             //
// starlattice()                                                              //
// samelattice()                                                              //
// .get(integer)                                                              //
// .get_hkl()                                                                 //
// .get_xyz()                                                                 //
// .star()                                                                    //
// .dot(integer,integer)                                                      //
// [+-]=(L[DQ]Vec)                                                            //
// [+-/*]=(scalar)                                                            //
//============================================================================//
// Functions operating on class objects:                                      //
//----------------------------------------------------------------------------//
// cross(L[DQ]Vec, L[DQ]Vec)                                                  //
// dot(L[DQ]Vec, L[DQ]Vec)                                                    //
// dot(L[DQ]Vec, L[QD]Vec)                                                    //
// star(L[DQ]Vec)                                                             //
// [+-/*](L[DQ]Vec,L[DQ]Vec)                                                  //
// [+-/*](L[DQ]Vec,Array)                                               //
// [+-/*](Array,L[DQ]Vec)                                               //
////////////////////////////////////////////////////////////////////////////////
