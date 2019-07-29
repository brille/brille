#include <random>
#include <chrono>
#include <catch.hpp>
#include "bz.h"

TEST_CASE("BrillouinZone instantiation","[brillouinzone]"){
  Direct d(3.,3.,3.,PI/2,PI/2,PI/2*3);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
}

TEST_CASE("BrillouinZone moveinto","[brillouinzone]"){
  std::string spgr;
  SECTION("Primitive spacegroup"){
    spgr = "P 1";
  }
  SECTION("Body-centred spacegroup"){
    spgr = "Im-3m";
  }
  SECTION("Face-centred spacegroup"){
    spgr = "Fmm2";
  }
  Direct d(2.87,2.87,2.87,PI/2,PI/2,PI/2,spgr);
  Reciprocal r = d.star();
  BrillouinZone bz(r);

  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(-5.0,5.0);

  int nQ = 3333;
  double* rawQ = new double[nQ*3]();
  for (int i=0; i<3*nQ; ++i) rawQ[i] = distribution(generator);
  LQVec<double> Q(r,nQ,rawQ);
  delete[] rawQ;

  LQVec<double> q(r,nQ);
  LQVec<int> tau(r,nQ);

  REQUIRE( bz.moveinto(Q,q,tau) ); // success indicated by return of true
  REQUIRE( bz.isinside(q).all_true() );
  LQVec<double> Qmq = Q-q;
  LQVec<double> Qmqmtau = Q-(q+tau);
  for (size_t i=0; i<Q.size(); ++i)
  for (size_t j=0; j<Q.numel(); ++j)
  REQUIRE( Q.getvalue(i,j) == Approx( q.getvalue(i,j) + tau.getvalue(i,j) ) );
}
