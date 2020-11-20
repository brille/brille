#include <catch2/catch.hpp>
#include <tuple>
#include <omp.h>
#include <complex>
#include "debug.hpp"
#include "bz_trellis.hpp"

using namespace brille;

TEST_CASE("BrillouinZoneTrellis3 instantiation","[trellis]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, brille::halfpi, brille::halfpi, brille::halfpi, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.01;
  // BrillouinZoneTrellis3<double> bzt0(bz); !! No default maximum volume
  BrillouinZoneTrellis3<double,double> bzt1(bz, max_volume);
}
TEST_CASE("BrillouinZoneTrellis3 vertex accessors","[trellis][accessors]"){
  Direct d(10.75, 10.75, 10.75, brille::halfpi, brille::halfpi, brille::halfpi, 525);
  BrillouinZone bz(d.star());
  double max_volume = 0.002;
  BrillouinZoneTrellis3<double,double> bzt(bz, max_volume);

  SECTION("get_xyz"){auto verts = bzt.get_xyz(); REQUIRE(verts.size(0) > 0u);}
  SECTION("get_hkl"){auto verts = bzt.get_hkl(); REQUIRE(verts.size(0) > 0u);}
  SECTION("get_outer_xyz"){auto verts = bzt.get_outer_xyz(); REQUIRE(verts.size(0) > 0u);}
  SECTION("get_outer_hkl"){auto verts = bzt.get_outer_hkl(); REQUIRE(verts.size(0) > 0u);}
  SECTION("get_inner_xyz"){auto verts = bzt.get_inner_xyz(); REQUIRE(verts.size(0) > 0u);}
  SECTION("get_inner_hkl"){auto verts = bzt.get_inner_hkl(); REQUIRE(verts.size(0) > 0u);}
}

TEST_CASE("Simple BrillouinZoneTrellis3 interpolation","[trellis][debugging]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, brille::halfpi, brille::halfpi, brille::halfpi, 529);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.0001;
  BrillouinZoneTrellis3<double,double> bzt(bz, max_volume);

  auto Qmap = bzt.get_hkl();
  auto Qxyz = bzt.get_xyz();
  // At present converting from Array2 to Array is not possible, so the
  // original reshape will not work. Instead we must copy data by hand
  brille::shape_t tostoreshape{Qmap.size(0), 1u, Qmap.size(1)};
  brille::Array<double> tostore(tostoreshape);
  for (auto i: tostore.subItr()) tostore[i] = Qxyz.val(i[0], i[2]);

  std::array<unsigned,3> elements{0,3,0};
  RotatesLike rl = RotatesLike::Reciprocal;
  // make sure we store an (nQ, 1, 3) array to have one mode per Q
  bzt.replace_value_data(tostore , elements, rl);

  REQUIRE(bzt.data().branches() == 1u);

  brille::ind_t nQ = 10;
  // In order to have easily-interpretable results we need to ensure we only
  // interpolate at points within the irreducible meshed volume.
  // Use points distributed randomly in the Irreducible Polyhedron
  Polyhedron irp = bz.get_ir_polyhedron();
  LQVec<double> Q = LQVec<double>::from_invA(r, irp.rand_rejection(nQ));
  // Q are now random points in the irreducible Brillouin zone polyhedron

  // We may run into problems if any of the points are too close to the irBz
  // boundary where the components expressed in the primitive lattice are 0.5
  // Since different "std" library **versions** round 0.5 differently (I'm looking at you MSVC).
  // If 0.5 is rounded to 1 this simple test *will* fail since bz.moveinto rounds
  // the components of Q to try and find an equivalent q and tau.


  auto QinvA = Q.get_xyz();
  // QinvA is (at present) a brille::Array2<double> and so can not be reshaped
  // to 3D (maybe introducing conversion routines is a good idea?)
  // Instead make a new brille::Array and copy the Q values by hand:
  brille::shape_t antshp{nQ, 1u, 3u};
  brille::Array<double> antres(antshp);
  for (auto i: antres.subItr()) antres[i] = QinvA.val(i[0],i[2]);

  auto [intres, dummy] = bzt.ir_interpolate_at(Q, 1 /*thread*/);

  auto diff = intres - antres;
  // info_update("\nInterpolation results Expected results:\n",antres.to_string(intres));
  // info_update("\nRounded difference:\n",diff.to_string());

  if (!diff.round().all(0.,0)) for (size_t i = 0; i < nQ; ++i) {
      info_update_if(!diff.view(i).round().all(0.,0),
        "\nThe interpolation point Q = ", Q.to_string(i),
        "\n            returned result ", intres.to_string(i),
        "\n                 instead of ", antres.to_string(i), "\n");
  }
  REQUIRE( diff.round().all(0.,0) ); // this is not a great test :(
  for (auto i: diff.valItr()) REQUIRE(std::abs(i) < 2E-10);
}

TEST_CASE("BrillouinZoneTrellis3 interpolation timing","[.][trellis][timing]"){
  // The conventional cell for Nb
  Direct d(3.2598, 3.2598, 3.2598, brille::halfpi, brille::halfpi, brille::halfpi, 529);
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
  auto Qmap = bzt.get_hkl();

  brille::ind_t n_modes{3u}, n_pt{Qmap.size(0)};
  brille::shape_t valsh{n_pt, n_modes, 1u}, vecsh{n_pt, n_modes, 3u};
  brille::Array<double> eigenvalues(valsh);
  brille::Array<std::complex<double>> eigenvectors(vecsh);
  for (auto& i: eigenvalues.valItr()) i = distribution(generator);
  for (auto& i: eigenvectors.valItr()) i = std::complex<double>(distribution(generator), distribution(generator));
  std::array<unsigned,3> vals_el{{1,0,0}}, vecs_el{{0,3,0}};
  RotatesLike vecs_rt{RotatesLike::Reciprocal};

  bzt.replace_value_data(eigenvalues, vals_el); // scalars do not rotate, so any RotatesLike value is fine
  bzt.replace_vector_data(eigenvectors, vecs_el, vecs_rt);

  REQUIRE(bzt.data().branches() == n_modes);

  size_t nQ = 10000;//10000;
  LQVec<double> Q(r,nQ);
  for (auto& i: Q.valItr()) i = distribution(generator);

  brille::Array<double> intvals;
  brille::Array<std::complex<double>> intvecs;
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
  Direct d(3.2598, 3.2598, 3.2598, brille::halfpi, brille::halfpi, brille::halfpi, 529);
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
  auto Qmap = bzt.get_hkl();

  brille::ind_t n_modes{3u}, n_pt{Qmap.size(0)};
  brille::shape_t valsh{n_pt, n_modes, 1u}, vecsh{n_pt, n_modes, 3u};
  brille::Array<double> eigenvalues(valsh);
  brille::Array<std::complex<double>> eigenvectors(vecsh);
  for (auto& i: eigenvalues.valItr()) i = distribution(generator);
  for (auto& i: eigenvectors.valItr()) i = std::complex<double>(distribution(generator), distribution(generator));
  std::array<unsigned,3> vals_el{{1,0,0}}, vecs_el{{0,3,0}};
  RotatesLike vecs_rt{RotatesLike::Reciprocal};

  bzt.replace_value_data(eigenvalues, vals_el); // scalars do not rotate, so any RotatesLike value is fine
  bzt.replace_vector_data(eigenvectors, vecs_el, vecs_rt);

  REQUIRE(bzt.data().branches() == n_modes);

  size_t nQ = 100000;//10000;
  LQVec<double> Q(r,nQ);
  for (auto& i: Q.valItr()) i = distribution(generator);

  int threads = omp_get_max_threads();
  auto timer = Stopwatch<>();
  timer.tic();
  auto [vals_intres, vecs_intres] = bzt.ir_interpolate_at(Q, threads);
  timer.toc();
  info_update("Interpolation of ",nQ," points performed by ",threads, " threads in ",timer.average(),"+/-",timer.jitter()," msec");
}

TEST_CASE("BrillouinZoneTrellis3 creation time","[.][trellis][creation_profiling]"){
  Direct d(3.,3.,3., brille::halfpi, brille::halfpi, brille::halfpi, 1);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.0001;
  auto timer = Stopwatch<>();
  timer.tic();
  BrillouinZoneTrellis3<double,std::complex<double>> bzt(bz, max_volume);
  timer.toc();
  info_update("Creation of Cubic BrilluoinZoneTrellis3 object with max_volume ",max_volume," in ",timer.average()," msec");

  Direct quartz_d(4.85235, 4.85235, 5.350305, brille::halfpi, brille::halfpi, 2*brille::pi/3, 443);
  BrillouinZone quartz_bz(quartz_d.star());
  timer.tic();
  BrillouinZoneTrellis3<double,std::complex<double>> quartz_bzt(quartz_bz, max_volume);
  timer.toc();
  info_update("Creation of Quartz BrilluoinZoneTrellis3 object with max_volume ",max_volume," in ",timer.average()," msec");
}


TEST_CASE("BrillouinZoneTrellis3 contains Gamma","[trellis][gamma]"){
  Direct d(3.,3.,3., brille::halfpi, brille::halfpi, brille::halfpi, 1);
  Reciprocal r = d.star();
  BrillouinZone bz(r);
  double max_volume = 0.1;
  BrillouinZoneTrellis3<double,std::complex<double>> bzt(bz, max_volume);

  auto diff = norm(bzt.vertices()).find(brille::cmp::eq,0.);
  REQUIRE(diff.size() == 1u);
}

TEST_CASE("PolyhedronTrellis creation","[polyhedron][trellis]"){
  std::vector<std::array<double,3>> va_verts{
    {0.16911452700846194,    -0.11958202884495578,     0.62136644889376591},
    {5.9764222286586582e-18, -1.6903874740503940e-17, -9.2979218529888708e-34},
    {2.1958566835792056e-13, -2.4632350748338648e-13,  0.62136644889376336},
    {0.33822905401652814,    -0.059791014422781175,    0.51780537407842375},
    {0.16911452700828469,     0.23916405768938562,     0.41424429926267342},
    {0.33822905401644060,    -0.23916405768994600,     0.41424429926267337}
  };
  auto verts = bArray<double>::from_std(va_verts);
  std::vector<std::vector<int>> verts_per_face{
    {2,0,3,4},
    {0,5,3},
    {2,4,1},
    {3,5,1,4},
    {0,2,1,5}
  };
  Polyhedron poly(verts, verts_per_face);

  for (auto x: {0.1, 0.01, 0.001, 0.0001})
    REQUIRE_NOTHROW(PolyhedronTrellis<double,std::complex<double>>(poly, x));
}


TEST_CASE("PolyhedronTrellis construction with 'sharp' polyhedron input","[trellis][polyhedron][intersection][la2zr2o7]"){
  // lattice information and generators of spacegroup via brilleu/CASTEP
  std::vector<double> latmat {7.583912824349999, 1.8412792137035698e-32, 0.,
                              3.791956412170034, 3.791956412170034, 5.362636186024768,
                              3.791956412170034,-3.791956412170034, 5.362636186024768};
  //
  std::vector<int> strides{3*sizeof(double), sizeof(double)};
  Direct dlat(latmat.data(), strides, "P_1"); // not P₁ but it *is* primitive
  // the generators are: 4-fold [1 -1 1], 2-fold [-1 1 1], 3-fold [1 1 -3], -𝟙
  // row-ordered generator matrices
  std::vector<std::array<int,9>> W {
    {{ 0,-1, 0,  0, 0,-1,  1, 1, 1}},
    {{-1,-1,-1,  0, 0, 1,  0, 1, 0}},
    {{-1,-1,-1,  1, 0, 0,  0, 0, 1}},
    {{-1, 0, 0,  0,-1, 0,  0, 0,-1}}
  };
  std::vector<std::array<double,3>> w{
    {{0.0, 0.0, 0.5}},
    {{0.5, 0.0, 0.0}},
    {{0.5, 0.0, 0.0}},
    {{0.0, 0.0, 0.0}}
  };
  Symmetry::Motions mots;
  mots.reserve(W.size());
  for (size_t i=0; i<W.size(); ++i) mots.push_back(Motion<int,double>(W[i], w[i]));
  Symmetry sym(mots);
  //
  dlat.set_spacegroup_symmetry(sym.generate());
  //
  BrillouinZone bz(dlat.star());
  // If the tolerance used in determine_tols is too small the following will
  // fail due to one or more missing vertices in a triangulated cube
  REQUIRE_NOTHROW(BrillouinZoneTrellis3<double,std::complex<double>>(bz, 0.001));
  BrillouinZoneTrellis3<double,std::complex<double>> bzt(bz, 0.001);

  // The total volume of all of the nodes should be the same as the irreducible
  // Brillouin zone polyhedron volume if the nodes fully fill the polyhedron.
  // When a polyhedron is 'sharp', as this one is, it may intersect with a node
  // cube *without* having any of its vertices in the cube or any of the cube
  // vertices in the polyhedron.
  // The PolyhedronTrellis constructor would happily construct a trellis with
  // imporant (overlapping) nodes inserted as NullNode objects instead of
  // PolyNode objects like they should be.
  REQUIRE(bzt.total_node_volume() == Approx(bz.get_ir_polyhedron().get_volume()));
}


TEST_CASE("PolyhedronTrellis construction from 'P1' hexagonal system must contain Gamma","[trellis][brillem]"){
  // this test case was exposed as problematic by the tests in brillem
  Reciprocal rlat(1.154701, 1.154701, 1, 90, 90, 60);
  BrillouinZone bz(rlat);
  REQUIRE_NOTHROW(BrillouinZoneTrellis3<double,double>(bz, 0.002));
}
