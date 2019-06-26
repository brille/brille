#include <random>
#include <chrono>
#include <catch.hpp>

#include "grid.h"
#include "bz_grid.h"
#include "interpolation.h"

ArrayVector<double> f_of_Q(const ArrayVector<double>& Q){
  size_t i,j,n = 40u;
  ArrayVector<double> out( n, Q.size() );
  for (i=0; i<Q.size(); ++i){
    for (j=0; j<n; ++j) out.insert(0.0, i,j);
    // mode 0
    out.insert( Q.getvalue(i,0)+Q.getvalue(i,1)+Q.getvalue(i,2), i, 0); // first scalar
    for (j=0; j<3u; ++j) out.insert( Q.getvalue(i,j), i, 1u+j*4u); // first matrix
    // mode 1
    out.insert( Q.getvalue(i,0)+Q.getvalue(i,1)-Q.getvalue(i,2), i, 10u);
    for (j=0; j<3u; ++j) out.insert( -1.0*Q.getvalue(i,j), i, 11u+j*4u);
    // mode 2
    out.insert( Q.getvalue(i,0)-Q.getvalue(i,1)-Q.getvalue(i,2), i, 20u);
    for (j=0; j<3u; ++j) out.insert( 2.0*Q.getvalue(i,j), i, 21u+j*4u);
    // mode 3
    out.insert( -Q.getvalue(i,0)-Q.getvalue(i,1)-Q.getvalue(i,2), i, 30u);
    for (j=0; j<3u; ++j) out.insert( -2.0*Q.getvalue(i,j), i, 31u+j*4u);
  }
  return out;
}
ArrayVector<double> f_of_Q_mats(const ArrayVector<double>& Q){
  size_t i,j,n = 36u;
  ArrayVector<double> out( n, Q.size() );
  for (i=0; i<Q.size(); ++i){
    for (j=0; j<n; ++j) out.insert(0.0, i,j);
    // mode 0
    for (j=0; j<3u; ++j) out.insert( Q.getvalue(i,j), i, j*4u); // first matrix
    // mode 1
    for (j=0; j<3u; ++j) out.insert(-1.0*Q.getvalue(i,j), i,  9u+j*4u);
    // mode 2
    for (j=0; j<3u; ++j) out.insert( 2.0*Q.getvalue(i,j), i, 18u+j*4u);
    // mode 3
    for (j=0; j<3u; ++j) out.insert(-2.0*Q.getvalue(i,j), i, 27u+j*4u);
  }
  return out;
}
ArrayVector<double> f_of_Q_eigs(const ArrayVector<double>& Q){
  size_t i,j,n = 4u;
  ArrayVector<double> out( n, Q.size() );
  for (i=0; i<Q.size(); ++i){
    for (j=0; j<n; ++j) out.insert(0.0, i,j);
    // mode 0
    out.insert( Q.getvalue(i,0)+Q.getvalue(i,1)+Q.getvalue(i,2), i, 0); // first scalar
    // mode 1
    out.insert( Q.getvalue(i,0)+Q.getvalue(i,1)-Q.getvalue(i,2), i, 1u);
    // mode 2
    out.insert( Q.getvalue(i,0)-Q.getvalue(i,1)-Q.getvalue(i,2), i, 2u);
    // mode 3
    out.insert(-Q.getvalue(i,0)-Q.getvalue(i,1)-Q.getvalue(i,2), i, 3u);
  }
  return out;
}


TEST_CASE("BrillouinZoneGrid3 Interpolation","[interpolation]"){
  Reciprocal r(1.,1.,1., PI/2, PI/2, PI/2);
  BrillouinZone bz(r);
  size_t halfN[3] = {10,10,10};
  BrillouinZoneGrid3<double> bzg(bz,halfN);

  ArrayVector<double> Qmap = bzg.get_mapped_xyz();
  bzg.replace_data( f_of_Q(Qmap) );

  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(-0.5,0.5);

  int nQ = 10;//10000;
  double* rawQ = new double[nQ*3]();
  for (int i=0; i<3*nQ; ++i) rawQ[i] = distribution(generator);
  LQVec<double> Q(r,nQ,rawQ);
  delete[] rawQ;

  ArrayVector<double> intres = bzg.linear_interpolate_at(Q);
  ArrayVector<double> antres = f_of_Q( Q.get_xyz() );

  ArrayVector<double> diff = intres - antres;
  // printf("\nInterpolation results:\n");
  // intres.print();
  // printf("\nExpected results:\n");
  // antres.print();
  // printf("\nRounded difference:\n");
  // diff.round().print();

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

TEST_CASE("BrillouinZoneGrid3 Sorting","[munkres]"){
  Reciprocal r(1.,1.,1., PI/2, PI/2, PI/2);
  BrillouinZone bz(r);
  size_t halfN[3] = {10,10,10};
  BrillouinZoneGrid3<double> bzg(bz,halfN);

  ArrayVector<size_t> newshape(1u,4u); // (nQ, 3,3, nModes)
  ArrayVector<double> Qmap = bzg.get_mapped_xyz();
  newshape.insert( Qmap.size(), 0u );
  newshape.insert( 3u, 1u );
  newshape.insert( 3u, 2u );
  newshape.insert( 4u, 3u );
  bzg.replace_data( f_of_Q_mats( Qmap ), newshape ); // maybe mapped_hkl instead?

  ArrayVector<size_t> sortperm = bzg.centre_sort_perm();
}

ArrayVector<double> fe_dispersion(const ArrayVector<double>& Q){
  ArrayVector<double> out( 1u, Q.size() );
  double J = 16, delta = 0.1;
  for (size_t i=0; i<Q.size(); ++i)
    out.insert( delta+8*J*(1-std::cos(PI*Q.getvalue(i,0))*std::cos(PI*Q.getvalue(i,1))*std::cos(PI*Q.getvalue(i,2))), i,0);
  return out;
}
ArrayVector<std::complex<double>> complex_fe_dispersion(const ArrayVector<double>& Q){
  ArrayVector<std::complex<double>> out( 1u, Q.size() );
  double J = -16, delta = -0.1;
  for (size_t i=0; i<Q.size(); ++i)
    out.insert( delta+8*J*(1-std::cos(PI*Q.getvalue(i,0))*std::cos(PI*Q.getvalue(i,1))*std::cos(PI*Q.getvalue(i,2))), i,0);
  return out;
}

TEST_CASE("primitive BrillouinZoneGrid3 Interpolation","[interpolation]"){
  std::string spgr = "Im-3m";
  Direct d(2.87,2.87,2.87,PI/2,PI/2,PI/2,spgr);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  size_t halfN[3] = {30,30,30};
  BrillouinZoneGrid3<std::complex<double>> bzg(bz,halfN);

  bzg.replace_data( fe_dispersion(bzg.get_mapped_xyz()) );

  LQVec<double> Q=LQVec<double>(r,bzg.get_mapped_hkl());

  ArrayVector<std::complex<double>> intres;
  SECTION("Single-thread interpolation"){
    intres = bzg.linear_interpolate_at(Q);
  }
  SECTION("Multi-thread interpolation"){
    intres = bzg.parallel_linear_interpolate_at(Q,10);
  }
  ArrayVector<std::complex<double>> fedisp = fe_dispersion(Q.get_xyz());
  ArrayVector<std::complex<double>> diff = intres - fedisp;
  ArrayVector<std::complex<double>> sum = intres+fedisp;
  for (size_t i=0; i<diff.size(); ++i)
  for (size_t j=0; j<diff.numel(); ++j)
  REQUIRE( Approx(abs(diff.getvalue(i,j))) <= 1e-10*abs(sum.getvalue(i,j)) );
}

// The following test is not a great general test-case since we're checking
// if a linear interpolation has done a good job approximating a
// decidedly-not-linear function.
// The tag [.] indicates to Catch that this test should be hidden and only run
// if its second tag is specified.
TEST_CASE("primitive BrillouinZoneGrid3 Interpolation, random","[.][random]"){
  std::string spgr = "Im-3m";
  Direct d(2.87,2.87,2.87,PI/2,PI/2,PI/2,spgr);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  size_t halfN[3] = {30,30,30};
  BrillouinZoneGrid3<double> bzg(bz,halfN);

  bzg.replace_data( fe_dispersion(bzg.get_mapped_xyz()) );

  ArrayVector<double> intres;

  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(-5,5);

  int nQ = 333;
  double* rawQ = new double[nQ*3]();
  for (int i=0; i<3*nQ; ++i) rawQ[i] = distribution(generator);
  LQVec<double> Q(r,nQ,rawQ);
  delete[] rawQ;
  LQVec<double> q(r,nQ);
  LQVec<int> tau(r,nQ);
  REQUIRE( bz.moveinto(Q,q,tau) );
  // printf("Q = \n"); Q.get_xyz().print();
  // printf("q = \n"); q.get_xyz().print();

  SECTION("Single-thread interpolation"){
    intres = bzg.linear_interpolate_at(q.get_xyz());
  }
  SECTION("Multi-thread interpolation"){
    intres = bzg.parallel_linear_interpolate_at(q.get_xyz(),1);
  }

  ArrayVector<double> fedisp = fe_dispersion(q.get_xyz());
  ArrayVector<double> diff = intres - fedisp;
  ArrayVector<double> percent = 100*diff/fedisp;


  for (size_t i=0; i<diff.size(); ++i)
  for (size_t j=0; j<diff.numel(); ++j)
  REQUIRE( abs(percent.getvalue(i,j)) <= Approx(35) );
}
