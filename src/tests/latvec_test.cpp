#include <catch2/catch.hpp>

#include "lattice.hpp"
#include "latvec.hpp"

TEST_CASE("Lattice Vector tests","[latvec]"){
  Direct d(1.,1.,1.);
  Reciprocal r = d.star();

  double values[] = { 1,0,0, 0,1,0, 0,0,1 };
  LQVec<double> q(r,3,values);

  SECTION("length via .norm(integer)"){
    for(size_t i=0; i<3; i++)
    REQUIRE( q.norm(i) == Approx(2*PI) );
  }
  SECTION("length via norm()"){
    ArrayVector<double> normq = norm(q);
    REQUIRE( normq.size() == q.size() );
    REQUIRE( normq.numel() == 1u );
    for(size_t i=0; i<q.size(); ++i)
      REQUIRE( q.norm(i) == normq.getvalue(i) );
  }
  SECTION("Single-element access via [integer] and cross product .cross(integer,integer)"){
    LQVec<double> z = q[2];
    LQVec<double> xy = q.cross(0,1)*(norm(q[2])/norm(q[1])/norm(q[0]));
    REQUIRE( z.isapprox(xy) );
  }
}

////////////////////////////////////////////////////////////////////////////////
//  TODO: Implement tests for:                                                //
//============================================================================//
// Class methods:                                                             //
//----------------------------------------------------------------------------//
// instantiation via L[DQ]Vec(::lattice,::ArrayVector)                        //
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
// [+-/*](L[DQ]Vec,ArrayVector)                                               //
// [+-/*](ArrayVector,L[DQ]Vec)                                               //
////////////////////////////////////////////////////////////////////////////////
