#include <catch2/catch.hpp>
#include <tuple>
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
  BrillouinZoneTrellis3<double,double> bzt1(bz, max_volume);
}
TEST_CASE("BrillouinZoneTrellis3 vertex accessors","[trellis]"){
  Direct d(10.75, 10.75, 10.75, PI/2, PI/2, PI/2, 525);
  BrillouinZone bz(d.star());
  double max_volume = 0.002;
  BrillouinZoneTrellis3<double,double> bzt(bz, max_volume);

  SECTION("get_xyz"){auto verts = bzt.get_xyz(); REQUIRE(verts.size() > 0u);}
  SECTION("get_hkl"){auto verts = bzt.get_hkl(); REQUIRE(verts.size() > 0u);}
  SECTION("get_outer_xyz"){auto verts = bzt.get_outer_xyz(); REQUIRE(verts.size() > 0u);}
  SECTION("get_outer_hkl"){auto verts = bzt.get_outer_hkl(); REQUIRE(verts.size() > 0u);}
  SECTION("get_inner_xyz"){auto verts = bzt.get_inner_xyz(); REQUIRE(verts.size() > 0u);}
  SECTION("get_inner_hkl"){auto verts = bzt.get_inner_hkl(); REQUIRE(verts.size() > 0u);}
}

TEST_CASE("Simple BrillouinZoneTrellis3 interpolation","[trellis][debugging]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, PI/2, PI/2, PI/2, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.0001;
  BrillouinZoneTrellis3<double,double> bzt(bz, max_volume);

  ArrayVector<double> Qmap = bzt.get_hkl();
  std::vector<size_t> shape{Qmap.size(), 3};
  std::array<unsigned long,3> elements{0,3,0};
  RotatesLike rl = RotatesLike::Reciprocal;
  bzt.replace_value_data( bzt.get_xyz(), shape, elements, rl);

  size_t nQ = 10;
  // In order to have easily-interpretable results we need to ensure we only
  // interpolate at points within the irreducible meshed volume.
  // Use points distributed randomly in the Irreducible Polyhedron
  Polyhedron irp = bz.get_ir_polyhedron();
  ArrayVector<double> Qxyz = irp.rand_rejection(nQ);
  LQVec<double> Q(r, nQ);
  double fromxyz[9];
  r.get_inverse_xyz_transform(fromxyz);
  for (size_t i=0; i<nQ; ++i)
    multiply_matrix_vector<double,double,double,3>(Q.data(i), fromxyz, Qxyz.data(i));
  // Q are now random points in the irreducible Brillouin zone polyhedron

  // We may run into problems if any of the points are too close to the irBz
  // boundary where the components expressed in the primitive lattice are 0.5
  // Since different "std" library **versions** round 0.5 differently (I'm looking at you MSVC).
  // If 0.5 is rounded to 1 this simple test *will* fail since bz.moveinto rounds
  // the components of Q to try and find an equivalent q and tau.


  ArrayVector<double> intres, dummy, antres=Q.get_xyz();
  std::tie(intres, dummy) = bzt.ir_interpolate_at(Q, 1 /*thread*/);

  ArrayVector<double> diff = intres - antres;
  // info_update("\nInterpolation results Expected results:\n",antres.to_string(intres));
  // info_update("\nRounded difference:\n",diff.to_string());

  if (!diff.round().all_zero()) for (size_t i = 0; i < nQ; ++i) {
      info_update_if(!diff.extract(i).round().all_zero(),
        "\nThe interpolation point Q = ", Q.to_string(i),
        "\n            returned result ", intres.to_string(i),
        "\n                 instead of ", antres.to_string(i), "\n");
  }
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

  BrillouinZoneTrellis3<double,std::complex<double>> bzt(bz, max_volume);
  ArrayVector<double> Qmap = bzt.get_hkl();

  size_t n_modes{3u};
  ArrayVector<double> eigenvalues(1*n_modes, Qmap.size());
  ArrayVector<std::complex<double>> eigenvectors(3*n_modes, Qmap.size());
  for (size_t i=0; i<Qmap.size(); ++i) for (size_t b=0; b<n_modes; ++b){
    eigenvalues.insert( distribution(generator), i, b);
    for (size_t j=0; j<3u; ++j)
      eigenvectors.insert( std::complex<double>(distribution(generator), distribution(generator)), i, b*3u+j);
  }
  std::vector<size_t> vals_sh{Qmap.size(), n_modes, 1}, vecs_sh{Qmap.size(), n_modes, 3};
  std::array<unsigned long,3> vals_el{{1,0,0}}, vecs_el{{0,3,0}};
  RotatesLike vecs_rt{RotatesLike::Reciprocal};

  bzt.replace_value_data(eigenvalues, vals_sh, vals_el); // scalars do not rotate, so any RotatesLike value is fine
  bzt.replace_vector_data(eigenvectors, vecs_sh, vecs_el, vecs_rt);

  size_t nQ = 10000;//10000;
  LQVec<double> Q(r,nQ);
  for (size_t i=0; i<nQ; ++i) for (size_t j=0; j<3; ++j)
    Q.insert(distribution(generator), i,j);

  ArrayVector<double> intvals;
  ArrayVector<std::complex<double>> intvecs;
  auto timer = Stopwatch<>();
  for (int threads=1; threads<7; ++threads){
    bool again = true;
    timer.tic();
    while (again && timer.elapsed()<10000){
      std::tie(intvals, intvecs) = bzt.ir_interpolate_at(Q, threads);
      timer.toc();
      again = timer.jitter()/timer.average() > 0.002;
    }
    info_update("Interpolation of ",nQ," points performed by ",threads, " threads in ",timer.average(),"+/-",timer.jitter()," msec");
  }
}

TEST_CASE("BrillouinZoneTrellis3 interpolation profiling","[.][trellis][profiling]"){
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

  BrillouinZoneTrellis3<double,std::complex<double>> bzt(bz, max_volume);
  ArrayVector<double> Qmap = bzt.get_hkl();

  size_t n_modes{3u};
  ArrayVector<double> eigenvalues(1*n_modes, Qmap.size());
  ArrayVector<std::complex<double>> eigenvectors(3*n_modes, Qmap.size());
  for (size_t i=0; i<Qmap.size(); ++i) for (size_t b=0; b<n_modes; ++b){
    eigenvalues.insert( distribution(generator), i, b);
    for (size_t j=0; j<3u; ++j)
      eigenvectors.insert( std::complex<double>(distribution(generator), distribution(generator)), i, b*3u+j);
  }
  std::vector<size_t> vals_sh{Qmap.size(), n_modes, 1}, vecs_sh{Qmap.size(), n_modes, 3};
  std::array<unsigned long,3> vals_el{{1,0,0}}, vecs_el{{0,3,0}};
  RotatesLike vecs_rt{RotatesLike::Reciprocal};

  bzt.replace_value_data(eigenvalues, vals_sh, vals_el); // scalars do not rotate, so any RotatesLike value is fine
  bzt.replace_vector_data(eigenvectors, vecs_sh, vecs_el, vecs_rt);

  size_t nQ = 100000;//10000;
  LQVec<double> Q(r,nQ);
  for (size_t i=0; i<nQ; ++i) for (size_t j=0; j<3; ++j)
    Q.insert(distribution(generator), i,j);

  ArrayVector<double> vals_intres;
  ArrayVector<std::complex<double>> vecs_intres;
  int threads = omp_get_max_threads();
  auto timer = Stopwatch<>();
  timer.tic();
  std::tie(vals_intres, vecs_intres) = bzt.ir_interpolate_at(Q, threads);
  timer.toc();
  info_update("Interpolation of ",nQ," points performed by ",threads, " threads in ",timer.average(),"+/-",timer.jitter()," msec");
}

TEST_CASE("BrillouinZoneTrellis3 creation time","[.][trellis][creation_profiling]"){
  Direct d(3.,3.,3., PI/2, PI/2, PI/2, 1);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.0001;
  auto timer = Stopwatch<>();
  timer.tic();
  BrillouinZoneTrellis3<double,std::complex<double>> bzt(bz, max_volume);
  timer.toc();
  info_update("Creation of Cubic BrilluoinZoneTrellis3 object with max_volume ",max_volume," in ",timer.average()," msec");

  Direct quartz_d(4.85235, 4.85235, 5.350305, PI/2, PI/2, 2*PI/3, 443);
  BrillouinZone quartz_bz(quartz_d.star());
  timer.tic();
  BrillouinZoneTrellis3<double,std::complex<double>> quartz_bzt(quartz_bz, max_volume);
  timer.toc();
  info_update("Creation of Quartz BrilluoinZoneTrellis3 object with max_volume ",max_volume," in ",timer.average()," msec");
}


TEST_CASE("BrillouinZoneTrellis3 contains Gamma","[trellis][gamma]"){
  Direct d(3.,3.,3., PI/2, PI/2, PI/2, 1);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.001;
  BrillouinZoneTrellis3<double,std::complex<double>> bzt(bz, max_volume);

  auto diff = find(norm(bzt.vertices()).is_approx(Comp::eq,0.));
  REQUIRE(diff.size() == 1u);
}
