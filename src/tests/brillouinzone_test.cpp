#include <random>
#include <chrono>
#include <catch2/catch.hpp>
#include "bz.hpp"

TEST_CASE("Primitive Cubic BrillouinZone instantiation","[brillouinzone]"){
  Direct d(2*brille::pi,2*brille::pi,2*brille::pi,90.,90.,90.);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
}

TEST_CASE("Primitive Hexagonal BrillouinZone instantiation","[brillouinzone]"){
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
  std::vector<std::array<double,3>> rawQ(nQ);
  for (auto& a: rawQ) for (auto& q: a) q = distribution(generator);
  auto data = brille::Array<double>::from_std(rawQ); // performs a copy into contiguous memory
  LQVec<double> Q(r,data);

  LQVec<double> q(r,nQ);
  LQVec<int> tau(r,nQ);

  REQUIRE( bz.moveinto(Q,q,tau) ); // success indicated by return of true
  auto inside = bz.isinside(q);
  REQUIRE(std::count(inside.begin(), inside.end(), false) == 0);
  LQVec<double> Qmq = Q-q;
  LQVec<double> Qmqmtau = Q-(q+tau);
  for (auto i: SubIt(Q.shape())){
    REQUIRE(Q[i] == Approx(q[i] + tau[i]));
    REQUIRE(brille::approx::scalar(Qmq[i], static_cast<double>(tau[i])));
    REQUIRE(brille::approx::scalar(std::abs(Qmqmtau[i]), 0.));
  }
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
  LQVec<double> Q(r, brille::Array<double>::from_std(rawQ));
  LQVec<double> q(r,Q.size(0));
  LQVec<int> tau(r,Q.size(0));
  REQUIRE(bz.moveinto(Q,q,tau));
}
