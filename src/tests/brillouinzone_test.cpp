#include <random>
#include <chrono>
#include <catch2/catch.hpp>
#include "bz.hpp"

TEST_CASE("BrillouinZone instantiation","[brillouinzone]"){
  Direct d(3.,3.,3.,90.,90.,120.);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
}

TEST_CASE("BrillouinZone moveinto","[brillouinzone]"){
  std::string spgr;
  SECTION("Primitive"){
    spgr = "P 1";
  }
  SECTION("Body-centred"){
    spgr = "Im-3m";
  }
  SECTION("Face-centred"){
    spgr = "Fd-3c"; // the previous choice, Fmm2 is orthorhombic. Trying to use it with a=b=c prevents the BrillouinZone algorithm from working
  }
  Direct d(2.87,2.87,2.87,90.,90.,90.,spgr);
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

TEST_CASE("BrillouinZone moveinto hexagonal","[brillouinzone][moveinto]"){
  Direct d(3.,3.,9.,90.,90.,120.,1);
  Reciprocal r(d.star());
  BrillouinZone bz(r);
  std::vector<std::array<double,3>> rawQ;
  for (int i=0; i<21; i++){
    double x = static_cast<double>(i)/20.0;
    rawQ.push_back({x,x,x});
  }
  LQVec<double> Q(r,rawQ);
  LQVec<double> q(r,Q.size());
  LQVec<int> tau(r,Q.size());
  REQUIRE(bz.moveinto(Q,q,tau));
}
