#include <catch2/catch.hpp>
#include <tuple>
#include <omp.h>
#include <complex>
#include <filesystem>
#include "debug.hpp"
#include "bz_trellis.hpp"
#include "lattice_dual.hpp"

using namespace brille;
using namespace brille::lattice;
using namespace brille::math;

#ifdef USE_HIGHFIVE
template<class T, class R, class S>
bool write_read_test(const BrillouinZoneTrellis3<T,R,S>& source, const std::string& name){
  namespace fs = std::filesystem;
  auto temp_dir = fs::temp_directory_path();
  fs::path filepath = temp_dir;
  filepath /= fs::path("brille"+std::to_string(processid())+".h5");
  auto filename = filepath.string();

  // write the BrillouinZoneTrellis3 to disk:
  if(!source.to_hdf(filename, name))
    throw std::runtime_error("Problem writing to HDF file?");

  auto sink = BrillouinZoneTrellis3<T,R,S>::from_hdf(filename, name);

  return (source == sink);
}
#else
template<class T, class R, class S>
bool write_read_test(const BrillouinZoneTrellis3<T,R,S>&, const std::string&){
  return true;
}
#endif

TEST_CASE("BrillouinZoneTrellis3 instantiation","[trellis][simple]"){
  // The conventional cell for Nb
  std::string spg("-I 4 2 3"); // == 529
  auto lat = Direct<double>({3.2598, 3.2598, 3.2598},
                            {half_pi, half_pi, half_pi}, spg);
  BrillouinZone bz(lat);
  double max_volume = 0.01;
  // BrillouinZoneTrellis3<double> bzt0(bz); !! No default maximum volume
  BrillouinZoneTrellis3<double,double,double> bzt1(bz, max_volume);
  REQUIRE(write_read_test(bzt1, "first_bzt"));
}
TEST_CASE("BrillouinZoneTrellis3 vertex accessors","[trellis][accessors]"){
  std::string spg("F 4d 2 3 -1d"); // == 525
  auto lat = Direct<double>({10.75, 10.75, 10.75},
                            {half_pi, half_pi, half_pi}, spg);
  BrillouinZone bz(lat);
  double max_volume{bz.get_ir_polyhedron().volume()/50}; // previous maximum, 0.002, gave ~8 intersecting nodes
  BrillouinZoneTrellis3<double,double,double> bzt(bz, max_volume);

  // there should be trellis vertices of some sort:
  auto all_hkl = bzt.get_hkl();
  auto all_xyz = bzt.get_xyz();
  REQUIRE(all_hkl.size(0) == all_xyz.size(0));
  REQUIRE(all_hkl.size(0) > 0);

  // there should be triangulated polygon vertices:
  auto out_hkl = bzt.get_outer_hkl();
  auto out_xyz = bzt.get_outer_xyz();
  REQUIRE(out_hkl.size(0) == out_xyz.size(0));
  REQUIRE(out_hkl.size(0) > 0);

  // there *probably* should be non-triangulated cubic vertices
  // there's no definite relationship between the sizes of the three types of vectors
  // since one vertex can be used for both triangulated and non-triangulated nodes
  auto in_hkl = bzt.get_inner_hkl();
  auto in_xyz = bzt.get_inner_xyz();
  REQUIRE(in_hkl.size(0) == in_xyz.size(0));
  REQUIRE(in_hkl.size(0) > 0);
}

TEST_CASE("Simple BrillouinZoneTrellis3 interpolation","[trellis][debugging]"){
  // The conventional cell for Nb
  std::string spg("-I 4 2 3"); // == 529
  auto lat = Direct<double>({3.2598, 3.2598, 3.2598},  {half_pi, half_pi, half_pi}, spg);
  BrillouinZone bz(lat);
  double max_volume = 0.0001;
  BrillouinZoneTrellis3<double,double,double> bzt(bz, max_volume);

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
  auto irp = bz.get_ir_polyhedron();
  auto Q = irp.rand_rejection(nQ);
  // Q are now random points in the irreducible Brillouin zone polyhedron

  // We may run into problems if any of the points are too close to the irBz
  // boundary where the components expressed in the primitive lattice are 0.5
  // Since different "std" library **versions** round 0.5 differently (I'm looking at you MSVC).
  // If 0.5 is rounded to 1 this simple test *will* fail since bz.moveinto rounds
  // the components of Q to try and find an equivalent q and tau.


  auto QinvA = Q.xyz();
  // QinvA is (at present) a brille::Array2<double> and so can not be reshaped
  // to 3D (maybe introducing conversion routines is a good idea?)
  // Instead make a new brille::Array and copy the Q values by hand:
  brille::shape_t antshp{nQ, 1u, 3u};
  brille::Array<double> antres(antshp);
  for (auto i: antres.subItr()) antres[i] = QinvA.val(i[0],i[2]);

  auto [intres, dummy] = bzt.ir_interpolate_at(Q, false, 1 /*thread*/);

  auto diff = intres - antres;
  // info_update("\nInterpolation results Expected results:\n",antres.to_string(intres));
  // info_update("\nRounded difference:\n",diff.to_string());

  if (!diff.round().all(0,0)) for (ind_t i = 0; i < nQ; ++i) {
      info_update_if(!diff.view(i).round().all(0,0),
        "\nThe interpolation point Q = ", Q.to_string(i),
        "\n            returned result ", intres.to_string(i),
        "\n                 instead of ", antres.to_string(i), "\n");
  }
  REQUIRE( diff.round().all(0,0) ); // this is not a great test :(
  for (auto i: diff.valItr()) REQUIRE(std::abs(i) < 2E-10);
}

TEST_CASE("BrillouinZoneTrellis3 interpolation timing","[.][trellis][timing]"){
  // The conventional cell for Nb
  std::string spg("-I 4 2 3"); // == 529
  auto lat = Direct<double>({3.2598, 3.2598, 3.2598},  {half_pi, half_pi, half_pi}, spg);
  BrillouinZone bz(lat);
  double max_volume = 0.0001;

  // In order to have easily-interpretable results we need to ensure we only
  // interpolate at points within the irreducible meshed volume.
  // So let's stick to points that are random linear interpolations between
  // neighbouring mesh vertices
  std::default_random_engine generator(static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count()));
  std::uniform_real_distribution<double> distribution(-5.,5.);

  BrillouinZoneTrellis3<double,std::complex<double>,double> bzt(bz, max_volume);
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
  auto Q = LVec<double>(LengthUnit::inverse_angstrom, lat, nQ);
  for (auto& i: Q.valItr()) i = distribution(generator);

  brille::Array<double> intvals;
  brille::Array<std::complex<double>> intvecs;
  auto timer = Stopwatch<>();
  for (int threads=1; threads<7; ++threads){
    bool again = true;
    timer.tic();
    while (again && timer.elapsed()<10000){
      std::tie(intvals, intvecs) = bzt.ir_interpolate_at(Q, false, threads);
      timer.toc();
      again = timer.jitter()/timer.average() > 0.002;
    }
    info_update("Interpolation of ",nQ," points performed by ",threads, " threads in ",timer.average(),"+/-",timer.jitter()," msec");
  }
}

TEST_CASE("BrillouinZoneTrellis3 interpolation profiling","[.][trellis][profiling]"){
  // The conventional cell for Nb
  std::string spg("-I 4 2 3"); // == 529
  auto lat = Direct<double>({3.2598, 3.2598, 3.2598},  {half_pi, half_pi, half_pi}, spg);
  BrillouinZone bz(lat);
  double max_volume = 0.0001;

  // In order to have easily-interpretable results we need to ensure we only
  // interpolate at points within the irreducible meshed volume.
  // So let's stick to points that are random linear interpolations between
  // neighbouring mesh vertices
  std::default_random_engine generator(static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count()));
  std::uniform_real_distribution<double> distribution(-5.,5.);

  BrillouinZoneTrellis3<double,std::complex<double>,double> bzt(bz, max_volume);
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
  auto Q = LQVec<double>(lat, nQ);
  for (auto& i: Q.valItr()) i = distribution(generator);

  int threads = omp_get_max_threads();
  auto timer = Stopwatch<>();
  timer.tic();
  auto [vals_intres, vecs_intres] = bzt.ir_interpolate_at(Q, false, threads);
  timer.toc();
  info_update("Interpolation of ",nQ," points performed by ",threads, " threads in ",timer.average(),"+/-",timer.jitter()," msec");
}

TEST_CASE("BrillouinZoneTrellis3 creation time","[.][trellis][creation_profiling]"){
  auto lat = Direct<double>({3,3,3}, {half_pi, half_pi, half_pi}, "P 1");
  BrillouinZone bz(lat);
  double max_volume = 0.0001;
  auto timer = Stopwatch<>();
  timer.tic();
  BrillouinZoneTrellis3<double,std::complex<double>,double> bzt(bz, max_volume);
  timer.toc();
  info_update("Creation of Cubic BrilluoinZoneTrellis3 object with max_volume ",max_volume," in ",timer.average()," msec");

  auto quartz = Direct<double>({4.85235, 4.85235, 5.350305}, {half_pi, half_pi, two_thirds_pi}, "P 32 2\"");
  BrillouinZone quartz_bz(quartz);
  timer.tic();
  BrillouinZoneTrellis3<double,std::complex<double>,double> quartz_bzt(quartz_bz, max_volume);
  timer.toc();
  info_update("Creation of Quartz BrilluoinZoneTrellis3 object with max_volume ",max_volume," in ",timer.average()," msec");
}


TEST_CASE("BrillouinZoneTrellis3 contains Gamma","[trellis][gamma]"){
  auto lat = Direct<double>({3,3,3}, {half_pi, half_pi, half_pi}, "P 1");
  BrillouinZone bz(lat);
  double max_volume = 0.1;
  BrillouinZoneTrellis3<double,std::complex<double>,double> bzt(bz, max_volume);

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
  polyhedron::Faces face(verts_per_face);
  polyhedron::Poly<double, bArray> poly(verts, face);

  for (auto x: {0.1, 0.01, 0.001, 0.0001})
    REQUIRE_NOTHROW(polytrellis::PolyTrellis<double,std::complex<double>,double,bArray>(poly, x));
}


TEST_CASE("PolyhedronTrellis construction with 'sharp' polyhedron input","[trellis][polyhedron][intersection][la2zr2o7]"){
  // lattice information and generators of spacegroup via brilleu/CASTEP
  std::array<double,9> latmat {7.583912824349999, 1.8412792137035698e-32, 0.,
                              3.791956412170034, 3.791956412170034, 5.362636186024768,
                              3.791956412170034,-3.791956412170034, 5.362636186024768};
  //
  auto lat = Direct<double>(latmat, MatrixVectors::row, "P 1", ""); // not P₁ but it *is* primitive
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
  lat.spacegroup_symmetry(sym.generate());
  //
  BrillouinZone bz(lat);
  // If the tolerance used in determine_tols is too small the following will
  // fail due to one or more missing vertices in a triangulated cube
  REQUIRE_NOTHROW(BrillouinZoneTrellis3<double,std::complex<double>,double>(bz, 0.001));
  BrillouinZoneTrellis3<double,std::complex<double>,double> bzt(bz, 0.001);

  auto fbz = bz.get_polyhedron();
  auto irp = bz.get_ir_polyhedron();
  auto verts = bzt.get_xyz();
  info_update("First Bz\n", fbz.python_string(), "\nIrreducible Bz\n", irp.python_string());
  info_update("Triangulated points\nnp.array(", verts.to_string(), ")");
  info_update("Triangulated tetrahedra\n", bzt.get_vertices_per_tetrahedron());


  // The total volume of all of the nodes should be the same as the irreducible
  // Brillouin zone polyhedron volume if the nodes fully fill the polyhedron.
  // When a polyhedron is 'sharp', as this one is, it may intersect with a node
  // cube *without* having any of its vertices in the cube or any of the cube
  // vertices in the polyhedron.
  // The PolyhedronTrellis constructor would happily construct a trellis with
  // imporant (overlapping) nodes inserted as NullNode objects instead of
  // PolyNode objects like they should be.
  REQUIRE(bzt.total_node_volume() == Approx(bz.get_ir_polyhedron().volume()));
}


TEST_CASE("PolyhedronTrellis construction from 'P1' hexagonal system must contain Gamma","[trellis][brillem]"){
  // this test case was exposed as problematic by the tests in brillem
  auto lat = Reciprocal<double>({1.154701, 1.154701, 1}, {90, 90, 60}, "P 1");
  BrillouinZone bz(lat);
  REQUIRE_NOTHROW(BrillouinZoneTrellis3<double,double,double>(bz, 0.002));
}

TEST_CASE("PolyNode inclusion rounding error","[trellis][quartz][polynode][61]"){
  auto quartz = Direct<double>({4.85235, 4.85235, 5.350305}, {half_pi, half_pi, two_thirds_pi}, "P 32 2\"");
  BrillouinZone quartz_bz(quartz);
  auto max_volume = quartz_bz.get_ir_polyhedron().volume()/2000.;
  BrillouinZoneTrellis3<double,double,double> quartz_bzt(quartz_bz, max_volume);

  Array<double> zeros(quartz_bzt.get_hkl().size(0), 1u);
  std::array<ind_t,3> elements{{1, 0, 0}};
  RotatesLike rl{RotatesLike::Reciprocal};
  Interpolator<double> val(zeros, elements, rl);
  quartz_bzt.replace_data(val, val);

  // error identified for delta=1e-9
  // Check for more powers of 10 for extra assurances
  for (int i=-15; i<0; ++i){
    auto delta = std::pow(10., i);
    std::vector<std::array<double,3>> values{{-0.1+delta, -0.1, 0.}};
    auto q = LQVec<double>(quartz, bArray<double>::from_std(values));
    REQUIRE_NOTHROW(quartz_bzt.ir_interpolate_at(q, false, 1));
  }
}

TEST_CASE("BrillouinZoneTrellis3 inclusion data race error","[trellis][la2zr2o7][omp][60]"){
  // lattice information and generators of spacegroup via brilleu/CASTEP
//  std::array<double,9> latmat {7.583912824349999, 1.8412792137035698e-32, 0.,
//                              3.791956412170034, 3.791956412170034, 5.362636186024768,
//                              3.791956412170034,-3.791956412170034, 5.362636186024768};
  //
  std::vector<int> strides{3*sizeof(double), sizeof(double)};
  std::array<double,3> len{7.583912824346341, 7.583912824346341, 7.583912824346341};
  std::array<double,3> ang{60,60,60};
//  Direct dlat(latmat.data(), strides, "P_1"); // not P₁ but it *is* primitive
  auto lat = Direct<double>(len, ang, "P 1");
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
  lat.spacegroup_symmetry(sym.generate());

  BrillouinZone bz(lat);
  info_update("First Brillouin zone\n", bz.get_polyhedron().python_string());
  info_update("Irreducible Brillouin zone\n", bz.get_ir_polyhedron().python_string());

  auto max_volume = bz.get_ir_polyhedron().volume()/20.;
  BrillouinZoneTrellis3<double, double, double> bzt(bz, max_volume);

  Array<double> zeros(bzt.get_hkl().size(0), 1u);
  std::array<ind_t,3> elements{{1, 0, 0}};
  RotatesLike rl{RotatesLike::Reciprocal};
  Interpolator<double> val(zeros, elements, rl);
  bzt.replace_data(val, val);

  int n{500};
  std::vector<std::array<double,3>> values;
  values.reserve(n);
  auto frac = [n](double from, double to, int j){
    auto step = (to - from) / static_cast<double>(n);
    return from + j*step;
  };
  for (int i=0; i<n; ++i){
    values.push_back({1.5, frac(-1.5, 2.0, i), frac(4.0, -3.0, i)});
  }
  auto Q = LQVec<double>(lat, bArray<double>::from_std(values));
  auto q = LQVec<double>(lat, Q.shape());
  auto tau = LQVec<int>(lat, Q.shape());
  std::vector<size_t> r, invr;

  for (int i=0; i<4; ++i){
    auto nthread = std::pow(2, i);
    verbose_update("ir_moveinto ",nthread," threads");
    REQUIRE_NOTHROW(bz.ir_moveinto(Q, q, tau, r, invr, nthread));
    verbose_update("ir_interpolate_at ",nthread," threads");
    REQUIRE_NOTHROW(bzt.ir_interpolate_at(Q, false, nthread));
  }
}


TEST_CASE("Equivalent atom error for CaHgO2","[trellis][interpolation][63]"){
  std::array<double,9> row_vectors{
    1.7797736800000001, 1.027552813333333, 6.219738366666666,
   -1.7797736800000001, 1.027552813333333, 6.219738366666666,
    0.0,               -2.055105626666667, 6.219738366666666
  };
  std::vector<std::array<double,3>> atom_positions {
      // Symmetries below require all atoms have positions (x, x, x)
      // but Oxygen and Mercury atoms differ in 1e-8 digits for y compared
      // to x and z. We need a tolerance worse than twice the difference to
      // match since we will need to match y with (1-y) which has a difference
      // of twice the difference
      {0.89401258, 0.89401259, 0.89401258}, // O
      {0.10598742, 0.10598741, 0.10598742}, // O
      {0.5, 0.5, 0.5}, // Ca
      {0.0, 0.99999999, 0.0} // Hg
  };
  std::vector<ind_t> atom_types{2,2,0,1};

  std::vector<std::array<int,9>> rotations {
    { 1, 0, 0,   0, 1, 0,   0, 0, 1},
    {-1, 0, 0,   0,-1, 0,   0, 0,-1},
    { 0, 0, 1,   1, 0, 0,   0, 1, 0},
    { 0, 0,-1,  -1, 0, 0,   0,-1, 0},
    { 0, 1, 0,   0, 0, 1,   1, 0, 0},
    { 0,-1, 0,   0, 0,-1,  -1, 0, 0},
    { 0, 0,-1,   0,-1, 0,  -1, 0, 0},
    { 0, 0, 1,   0, 1, 0,   1, 0, 0},
    { 0,-1, 0,  -1, 0, 0,   0, 0,-1},
    { 0, 1, 0,   1, 0, 0,   0, 0, 1},
    {-1, 0, 0,   0, 0,-1,   0,-1, 0},
    { 1, 0, 0,   0, 0, 1,   0, 1, 0}
  };
  std::vector<std::array<double,3>> translations {
    {0.,0.,0.},
    {0.,0.,0.},
    {0.,0.,0.},
    {0.,0.,0.},
    {0.,0.,0.},
    {0.,0.,0.},
    {0.,0.,0.},
    {0.,0.,0.},
    {0.,0.,0.},
    {0.,0.,0.},
    {0.,0.,0.},
    {0.,0.,0.}
  };
  std::vector<Motion<int,double>> motions;
  for (size_t i=0; i<rotations.size(); ++i) motions.emplace_back(rotations[i], translations[i]);
  auto sym = Symmetry(motions);
  auto bas = Basis(atom_positions, atom_types);
  auto lat = Direct<double>(row_vectors, MatrixVectors::row, sym, bas);

  auto bz = BrillouinZone(lat);
  auto goal_node_volume = bz.get_ir_polyhedron().volume() / static_cast<double>(1000);

  // create a local configuration
  auto cfg = approx_float::Config().direct(2.1e-8).reciprocal(5e-10);
  auto bzt = BrillouinZoneTrellis3<double,std::complex<double>,double>(bz, goal_node_volume, false, cfg);

//  // or modify the global one, which annoyingly has to be done in two steps
//  approx_float::config.direct(2.1e-8);
//  approx_float::config.reciprocal(5e-10);
//  auto bzt = BrillouinZoneTrellis3<double,std::complex<double>,double>(bz, goal_node_volume, false);

  auto n_modes = static_cast<ind_t>(3u * atom_types.size());
  std::vector<ind_t> val_shape{bzt.get_hkl().size(0), n_modes, 1u};
  Array<double> values(val_shape, 1.0);
  std::vector<ind_t> vec_shape{bzt.get_hkl().size(0), n_modes, n_modes};
  Array<std::complex<double>> vectors(vec_shape, std::complex<double>(1.0));

  std::array<ind_t,3> val_elements{{1, 0, 0}};
  Interpolator<double> val(values, val_elements, RotatesLike::Real);
  std::array<ind_t,3> vec_elements{{0, n_modes, 0}};
  Interpolator<std::complex<double>> vec(vectors, vec_elements, RotatesLike::Gamma);
  bzt.replace_data(val, vec);

  std::vector<std::array<double,3>> std_qpts {{0.05, 0.05, 0.05}};
  auto qpts = LQVec<double>(lat, bArray<double>::from_std(std_qpts));

  REQUIRE_NOTHROW(bzt.ir_interpolate_at(qpts, false, 1));
}