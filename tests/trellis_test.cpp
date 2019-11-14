#include <catch.hpp>
#include <omp.h>
#include <complex>
#include "debug.hpp"
#include "bz_trellis.hpp"

TEST_CASE("BrillouinZoneTrellis3 instantiation","[trellis]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, PI/2, PI/2, PI/2, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.01;
  // BrillouinZoneTrellis3<double> bzt0(bz); !! No default maximum volume
  BrillouinZoneTrellis3<double> bzt1(bz, max_volume);
}
TEST_CASE("BrillouinZoneTrellis3 vertex accessors","[trellis]"){
  Direct d(10.75, 10.75, 10.75, PI/2, PI/2, PI/2, 525);
  BrillouinZone bz(d.star());
  double max_volume = 0.002;
  BrillouinZoneTrellis3<double> bzt(bz, max_volume);

  SECTION("get_xyz"){auto verts = bzt.get_xyz(); REQUIRE(verts.size() > 0u);}
  SECTION("get_hkl"){auto verts = bzt.get_hkl(); REQUIRE(verts.size() > 0u);}
  SECTION("get_outer_xyz"){auto verts = bzt.get_outer_xyz(); REQUIRE(verts.size() > 0u);}
  SECTION("get_outer_hkl"){auto verts = bzt.get_outer_hkl(); REQUIRE(verts.size() > 0u);}
  SECTION("get_inner_xyz"){auto verts = bzt.get_inner_xyz(); REQUIRE(verts.size() > 0u);}
  SECTION("get_inner_hkl"){auto verts = bzt.get_inner_hkl(); REQUIRE(verts.size() > 0u);}
}

TEST_CASE("Simple BrillouinZoneTrellis3 interpolation","[trellis]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, PI/2, PI/2, PI/2, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.0001;
  BrillouinZoneTrellis3<double> bzt(bz, max_volume);

  ArrayVector<double> Qmap = bzt.get_hkl();
  std::vector<size_t> shape{Qmap.size(), 3};
  std::array<unsigned long,4> elements{0,0,3,0};
  bzt.replace_data( bzt.get_xyz(), shape, elements);

  // In order to have easily-interpretable results we need to ensure we only
  // interpolate at points within the irreducible meshed volume.
  // So let's stick to points that are random linear interpolations between
  // neighbouring mesh vertices
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(0.,1.);

  size_t nQmap = Qmap.size(), nQ = 10000;//10000;
  LQVec<double> Q(r,nQ);
  double rli;
  for (size_t i=0; i<nQ; ++i){
    rli = distribution(generator);
    Q.set(i, rli*Qmap.extract(i%nQmap) + (1-rli)*Qmap.extract((i+1)%nQmap) );
  }

  ArrayVector<double> intres, antres=Q.get_xyz();
  intres = bzt.interpolate_at(Q, 1 /*thread*/);

  ArrayVector<double> diff = intres - antres;
  // printf("\nInterpolation results:\n");
  // intres.print();
  // printf("\nExpected results:\n");
  // antres.print();
  // printf("\nRounded difference:\n");
  // diff.round().print();

  REQUIRE( diff.round().all_zero() ); // this is not a great test :(
  for (size_t i=0; i<diff.size(); ++i)
  for (size_t j=0; j<diff.numel(); ++j)
  REQUIRE( abs(diff.getvalue(i,j))< 2E-10 );
}

TEST_CASE("BrillouinZoneTrellis3 interpolation timing","[.][trellis][timing]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, PI/2, PI/2, PI/2, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.0001;

  // In order to have easily-interpretable results we need to ensure we only
  // interpolate at points within the irreducible meshed volume.
  // So let's stick to points that are random linear interpolations between
  // neighbouring mesh vertices
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(-5.,5.);

  BrillouinZoneTrellis3<std::complex<double>> bzt(bz, max_volume);
  ArrayVector<double> Qmap = bzt.get_hkl();

  ArrayVector<std::complex<double>> data((1+3)*3, Qmap.size());
  for (size_t i=0; i<data.size(); ++i) for (size_t j=0; j<data.numel(); ++j)
    data.insert( std::complex<double>(distribution(generator), distribution(generator)), i, j);
  std::vector<size_t> shape{data.size(), 3, 4};
  std::array<unsigned long,4> elements{1,3,0,0};

  bzt.replace_data(data, shape, elements);

  size_t nQ = 10000;//10000;
  LQVec<double> Q(r,nQ);
  for (size_t i=0; i<nQ; ++i) for (size_t j=0; j<3; ++j)
    Q.insert(distribution(generator), i,j);

  ArrayVector<std::complex<double>> intres;
  Stopwatch timer = Stopwatch<>();
  for (int threads=1; threads<10; ++threads){
    bool again = true;
    timer.tic();
    while (again && timer.elapsed()<10000){
      intres = bzt.interpolate_at(Q, threads);
      timer.toc();
      again = timer.jitter()/timer.average() > 0.002;
    }
    info_update("Interpolation of ",nQ," points performed by ",threads, " threads in ",timer.average(),"+/-",timer.jitter()," msec");
  }
}
