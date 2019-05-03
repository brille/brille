#include <random>
#include <chrono>
#include <catch.hpp>

#include "grid.h"
#include "bz_grid.h"

ArrayVector<double> f_of_Q(const ArrayVector<double>& Q){
  ArrayVector<double> out( 9u, Q.size() );
  for (size_t i=0; i<Q.size(); ++i){
    for (size_t j=0; j<9u; ++j) out.insert(0.0, i,j);
    for (size_t j=0; j<3u; ++j) out.insert( Q.getvalue(i,j), i, j*4u);
  }
  return out;
}

TEST_CASE("Testing BrillouinZoneGrid3 Interpolation"){
  Reciprocal r(1.,1.,1., PI/2, PI/2, PI/2);
  BrillouinZone bz(r);
  size_t halfN[3] = {10,10,10};
  BrillouinZoneGrid3<double> bzg(bz,halfN);

  bzg.replace_data( f_of_Q( bzg.get_mapped_xyz() ) ); // maybe mapped_hkl instead?

  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(-0.5,0.5);

  int nQ = 100;
  double* rawQ = new double[nQ*3]();
  for (int i=0; i<nQ; ++i) rawQ[i] = distribution(generator);
  LQVec<double> Q(r,nQ,rawQ);
  delete[] rawQ;

  ArrayVector<double> intres = bzg.linear_interpolate_at(Q);
  ArrayVector<double> antres = f_of_Q( Q.get_xyz() );

  ArrayVector<double> diff = intres - antres;
  REQUIRE( diff.round().areallzero() ); // this is not a great test :(
  for (size_t i=0; i<diff.size(); ++i)
  for (size_t j=0; j<diff.numel(); ++j)
  REQUIRE( abs(diff.getvalue(i,j))< 2E-14 );


  intres = bzg.parallel_linear_interpolate_at(Q,0);
  diff = intres - antres;
  REQUIRE( diff.round().areallzero() ); // this is not a great test :(
  for (size_t i=0; i<diff.size(); ++i)
  for (size_t j=0; j<diff.numel(); ++j)
  REQUIRE( abs(diff.getvalue(i,j))< 2E-14 );
}
