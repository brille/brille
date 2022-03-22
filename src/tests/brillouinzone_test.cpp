#include <random>
#include <chrono>
#include <catch2/catch.hpp>
#include "bz.hpp"
#include <filesystem>

using namespace brille;
using namespace brille::lattice;

#ifdef USE_HIGHFIVE
bool write_read_test(const BrillouinZone& source, const std::string& name){
    namespace fs = std::filesystem;
    auto temp_dir = fs::temp_directory_path();
    fs::path filepath = temp_dir;
    filepath /= fs::path("brille"+std::to_string(processid())+".h5");

    // write the BrillouinZone to disk:
    auto wrote_ok = source.to_hdf(filepath.string(), name);
    if (!wrote_ok) return false;

    auto sink = BrillouinZone::from_hdf(filepath.string(), name);

    return (source == sink);
}
#else
bool write_read_test(const BrillouinZone&, const std::string&){
	return true;
}
#endif

TEST_CASE("Primitive Cubic BrillouinZone instantiation","[bz_]"){
  auto lat = Direct<double>(
      {brille::math::two_pi, brille::math::two_pi, brille::math::two_pi},
      {90., 90., 90.} ,"P 1");
  BrillouinZone bz(lat);
  REQUIRE(write_read_test(bz, "primitive_cubic"));
}

TEST_CASE("P4 Cubic BrillouinZone instantiation","[bz_]"){
  auto lat = Direct<double>(
      {brille::math::two_pi, brille::math::two_pi, brille::math::two_pi},
      {90., 90., 90.} ,"P 2 2");
  BrillouinZone bz(lat);
  REQUIRE(write_read_test(bz, "primitive_cubic"));
  info_update("First Brillouin zone\n", bz.get_polyhedron().python_string());
  info_update("Irreducible Brillouin zone\n", bz.get_ir_polyhedron().python_string());
}


TEST_CASE("Primitive Hexagonal BrillouinZone instantiation","[bz_]"){
  auto lat = Direct<double>({3., 3., 3.}, {90., 90., 120.}, "P 1");
  BrillouinZone bz(lat);
  info_update("First Brillouin zone\n", bz.get_polyhedron().python_string());
  info_update("Irreducible Brillouin zone\n", bz.get_ir_polyhedron().python_string());
  REQUIRE(write_read_test(bz, "primitive_hexagonal"));
}

TEST_CASE("Rhombohedral Brillouin zone","[bz_]"){
  auto lat = Direct<double>({math::two_pi, math::two_pi, math::pi}, {90,90,120}, "-R 3");
  verbose_update(lat.to_verbose_string());
  auto ac = approx_float::Config().digit(10000); // 2e-8?
  auto cfg = BrillouinZoneConfig();
  cfg.divide_primitive(false); // can we leave this as true?
  BrillouinZone bz(lat, cfg, ac);
  REQUIRE(bz.check_ir_polyhedron());
  auto fbz = bz.get_polyhedron();
  auto irp = bz.get_ir_polyhedron();
  REQUIRE(irp.volume() == Approx(fbz.volume()/6));
  REQUIRE(write_read_test(bz, "mp-147"));
  info_update("First Brillouin zone\n", fbz.python_string());
  info_update("Irreducible Brillouin zone\n", irp.python_string());
}

//
//TEST_CASE("BrillouinZone moveinto","[bz_]"){
//  std::string spgr;
//  SECTION("Primitive"){
//    spgr = "P 1";
//  }
//  SECTION("Body-centred"){
//    spgr = "Im-3m";
//  }
//  SECTION("Face-centred"){
//    spgr = "Fd-3c"; // the previous choice, Fmm2 is orthorhombic. Trying to use it with a=b=c prevents the BrillouinZone algorithm from working
//  }
//  Direct d(2.87,2.87,2.87,90.,90.,90.,spgr);
//  Reciprocal r = d.star();
//  BrillouinZone bz(r);
//
//  std::default_random_engine generator(static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count()));
//  std::uniform_real_distribution<double> distribution(-5.0,5.0);
//
//  int nQ = 3333;
//  LQVec<double> Q(r, nQ); // holds an (nQ, 3) shaped Array or Array2
//  for (auto& v: Q.valItr()) v = distribution(generator);
//
//  LQVec<double> q(r,nQ);
//  LQVec<int> tau(r,nQ);
//
//  REQUIRE( bz.moveinto(Q,q,tau) ); // success indicated by return of true
//  auto inside = bz.isinside(q);
//  REQUIRE(std::count(inside.begin(), inside.end(), false) == 0);
//  LQVec<double> Qmq = Q-q;
//  LQVec<double> Qmqmtau = Q-(q+tau);
//  for (auto i: Q.subItr()){
//    REQUIRE(Q[i] == Approx(q[i] + tau[i]));
//    REQUIRE(brille::approx::scalar(Qmq[i], static_cast<double>(tau[i])));
//    REQUIRE(brille::approx::scalar(std::abs(Qmqmtau[i]), 0.));
//  }
//  REQUIRE(write_read_test(bz, spgr));
//}


TEST_CASE("BrillouinZone moveinto","[bz_]"){
  std::string spgr;
  auto run_tests = [](const std::string & spacegroup){
    auto lat =Direct<double>({2.87, 2.87, 2.87}, {90., 90., 90.},spacegroup);
    BrillouinZone bz(lat);

    std::default_random_engine generator(static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<double> distribution(-5.0,5.0);

    int nQ = 50;
    auto Q = LQVec<double>(lat, nQ); // holds an (nQ, 3) shaped Array or Array2
    for (auto& v: Q.valItr()) v = distribution(generator);

    auto q = LQVec<double>(lat, nQ);
    auto tau = LQVec<int>(lat, nQ);

    REQUIRE( bz.moveinto(Q,q,tau) ); // success indicated by return of true
    auto inside = bz.isinside(q);
    REQUIRE(std::count(inside.begin(), inside.end(), false) == 0);
    auto Qmq = Q-q;
    auto Qmqmtau = Q-(q+tau);
    for (auto i: Q.subItr()){
      REQUIRE(Q[i] == Approx(q[i] + tau[i]));
      REQUIRE(brille::approx::scalar(Qmq[i], static_cast<double>(tau[i])));
      REQUIRE(brille::approx::scalar(std::abs(Qmqmtau[i]), 0.));
    }
    REQUIRE(write_read_test(bz, spacegroup));
  };
  SECTION("Primitive"){
    spgr = "P 1";
    run_tests(spgr);
  }
  SECTION("Body-centred"){
    spgr = "Im-3m";
    run_tests(spgr);
  }
  SECTION("Face-centred"){
    spgr = "Fd-3c"; // the previous choice, Fmm2 is orthorhombic. Trying to use it with a=b=c prevents the BrillouinZone algorithm from working
    run_tests(spgr);
  }

}

TEST_CASE("BrillouinZone moveinto hexagonal","[bz_][moveinto]"){
  auto lat = Direct<double>({3,3,9}, {90,90,120}, "P 1");
  BrillouinZone bz(lat);
  std::vector<std::array<double,3>> rawQ;
  for (int i=0; i<21; i++){
    double x = static_cast<double>(i)/20.0;
    rawQ.push_back({x,x,x});
  }
  auto Q = LQVec<double>(lat, bArray<double>::from_std(rawQ));
  auto q = LQVec<double>(lat, Q.size(0));
  auto tau = LQVec<int>(lat,  Q.size(0));
  REQUIRE(bz.moveinto(Q,q,tau));
}

TEST_CASE("BrillouinZone moveinto hexagonal extended","[bz_][moveinto][.timing]"){
  auto lat = Direct<double>({3, 3, 9}, {90, 90, 120}, "P 1");
  BrillouinZone bz(lat);
  std::vector<std::array<double,3>> rawQ;

  std::default_random_engine g(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> d(-4.0, 4.0);
  auto random_array = [&]() -> std::array<double,3> { return {d(g), d(g), d(g)}; };

  size_t count{0}, total_points{10000};
  while (++count < total_points)
    rawQ.push_back(random_array());

  auto Q = LQVec<double>(lat, bArray<double>::from_std(rawQ));
  auto q = LQVec<double>(lat, Q.size(0));
  auto tau = LQVec<int>(lat, Q.size(0));
  REQUIRE_NOTHROW(bz.moveinto(Q,q,tau));
  
  // construct a timing object, defined in debug.hpp
  // has methods:
  //    tic       -- start/reset the timer
  //    toc       -- emit a split time and increment the number of splits counter
  //    average   -- divide the total elapsed time by the number of splits
  //    jitter    -- estimate the timing uncertainty assuming counting statistics in both timer and splits
  auto timer = Stopwatch<>();
  // setup explicit-number-of-threads region:
  omp_set_dynamic(0); // disables dynamic teams
  auto max_threads = omp_get_max_threads();
  std::vector<double> times;
  for (int threads=1; threads<max_threads; ++threads){
    omp_set_num_threads(threads);
    timer.tic();
    bz.moveinto(Q,Q,tau);
    times.push_back(timer.toc()); // keep the split-time (in msec)
    info_update("bz_::moveinto of ",total_points," points performed by ",threads, " threads in ",timer.average(),"+/-",timer.jitter()," msec");
  }
  for (size_t i=1u; i<times.size(); ++i)
    REQUIRE(times[i] <= times[0]);
}

TEST_CASE("Irreducible Brillouin zone for mp-147 alt","[bz_][materialsproject]"){
  // The spacegroup for elemental Se, from https://www.materialsproject.org/materials/mp-147/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-147.html
  // If brille::approx::scalar() is too strict this will fail.
  //
  std::string hall_symbol = "-R 3";
  // the other version below has slightly-off b vector length
  std::array<double,3> lens{10.699417459999999, 10.699417459999999, 4.528988750000000};
  std::array<double,3> angs{90,90,120};
  //
  auto lat = Direct<double>(lens, angs, hall_symbol);
  verbose_update(lat.to_verbose_string());
  auto ac = approx_float::Config().digit(10000);
  BrillouinZone bz(lat, ac);
  REQUIRE(bz.check_ir_polyhedron());
  auto fbz = bz.get_polyhedron();
  auto irp = bz.get_ir_polyhedron();
  REQUIRE(irp.volume() == Approx(fbz.volume()/6));
  REQUIRE(write_read_test(bz, "mp-147"));
  info_update("First Brillouin zone\n", fbz.python_string());
  info_update("Irreducible Brillouin zone\n", irp.python_string());
}

TEST_CASE("Irreducible Brillouin zone for mp-147 imprecise failure","[bz_][materialsproject]"){
  // The spacegroup for elemental Se, from https://www.materialsproject.org/materials/mp-147/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-147.html
  // If brille::approx::scalar() is too strict this will fail.
  std::array<double, 9> lattice_vectors {
    10.699417459999999, 0.000000050000000, 0.000000000000000,
    -5.349708760000000, 9.265967300000000, 0.000000000000000,
    0.000000000000000, 0.000000000000000, 4.528988750000000};
  //
  std::string hall_symbol = "-R 3";
  //
  auto lat = Direct<double>(lattice_vectors, MatrixVectors::row, hall_symbol);
  auto ac = approx_float::Config().digit(10000);
  REQUIRE_THROWS_AS(BrillouinZone(lat, ac), std::runtime_error);
}

TEST_CASE("Irreducible Brillouin zone for mp-147","[bz_][materialsproject]"){
  // The spacegroup for elemental Se, from https://www.materialsproject.org/materials/mp-147/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-147.html
  // If brille::approx::scalar() is too strict this will fail.
  std::array<double, 9> lattice_vectors{
      10.699417459999999, 0.000000000000000, 0.000000000000000,
      -5.349708729999997, 9.265967326054772, 0.000000000000000,
      0.000000000000000, 0.000000000000000, 4.528988750000000};
  //
  std::string hall_symbol = "-R 3";
  //
  auto lat = Direct<double>(lattice_vectors, MatrixVectors::row, hall_symbol);
  auto ac = approx_float::Config().digit(10000); // increase the equivalency sloppiness from 1000 to 10000 * epsilon
  BrillouinZone bz(lat, ac);
  REQUIRE(bz.check_ir_polyhedron());
  auto fbz = bz.get_polyhedron();
  auto irp = bz.get_ir_polyhedron();
  REQUIRE(irp.volume() == Approx(fbz.volume()/6));
  REQUIRE(write_read_test(bz, "mp-147"));
}

TEST_CASE("Irreducible Brillouin zone for mp-306","[bz_][materialsproject]"){
  // The spacegroup for B₂O₃, from https://www.materialsproject.org/materials/mp-306/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-306.html
  // If brille::approx::scalar() is too strict this will fail.
  std::array<double, 9> lattice_vectors{
     4.350620350000000, 0.000000000000000, 0.000000000000000,
    -2.175310180000000, 3.767747740000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 8.339814740000000};
  //
  std::string hall_symbol = "P 31 2\"";
  //
  auto lat = Direct<double>(lattice_vectors, MatrixVectors::row, hall_symbol);
  BrillouinZone bz(lat);
  REQUIRE(bz.check_ir_polyhedron());
  auto fbz = bz.get_polyhedron();
  auto irp = bz.get_ir_polyhedron();
  REQUIRE(irp.volume() == Approx(fbz.volume()/6));
  REQUIRE(write_read_test(bz, "mp-306"));
}

TEST_CASE("Irreducible Brillouin zone for mp-661","[bz_][materialsproject]"){
  // The spacegroup for AlN, from https://www.materialsproject.org/materials/mp-661/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-661.html
  // If brille::approx::scalar() is too strict this will fail.
  std::array<double, 9> lattice_vectors{
     3.113692560000000, 0.000000000000000, 0.000000000000000,
    -1.556846270000000, 2.696536850000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 4.983221820000000};
  //
  std::string hall_symbol = "P 6c -2c";
  //
  auto lat = Direct<double>(lattice_vectors, MatrixVectors::row, hall_symbol);
  BrillouinZone bz(lat);
  REQUIRE(bz.check_ir_polyhedron());
  auto fbz = bz.get_polyhedron();
  auto irp = bz.get_ir_polyhedron();
  REQUIRE(irp.volume() == Approx(fbz.volume()/12));
  REQUIRE(write_read_test(bz, "mp-661"));
}

TEST_CASE("Irreducible Brillouin zone for mp-7041","[bz_][materialsproject]"){
  // The spacegroup for CaHgO₂, from https://www.materialsproject.org/materials/mp-7041/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-7041.html
  // If brille::approx::scalar() is too strict this will fail.
  std::array<double, 9> lattice_vectors{
      3.559547360, 0.0,         0.0,
     -1.779773680, 3.082658440, 0.0,
      0.0,         0.0,        18.659215100000001};
  //
  std::string hall_symbol = "-R 3 2\"";
  //
  auto lat = Direct<double>(lattice_vectors, MatrixVectors::row, hall_symbol);
  BrillouinZone bz(lat);
  REQUIRE(bz.check_ir_polyhedron());
  auto fbz = bz.get_polyhedron();
  auto irp = bz.get_ir_polyhedron();
  REQUIRE(irp.volume() == Approx(fbz.volume()/12));
  REQUIRE(write_read_test(bz, "mp-7041"));
}

TEST_CASE("Irreducible Brillouin zone for mp-917 atl","[bz_][materialsproject]"){
  // The spacegroup for CaC₂, from https://www.materialsproject.org/materials/mp-917/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-917.html
  // If brille::approx::scalar() is too strict this will fail.
  std::array<double,3> len{7.16797000, 3.84256200, 7.46845574}, ang{90, 106.87385147, 90};
  std::string hall_symbol = "-C 2y";
  //
  auto lat = Direct<double>(len, ang, hall_symbol);
  auto ac = approx_float::Config().digit(10000);
  BrillouinZone bz(lat, ac);
  REQUIRE(bz.check_ir_polyhedron());
  auto fbz = bz.get_polyhedron();
  auto irp = bz.get_ir_polyhedron();
  REQUIRE(irp.volume() == Approx(fbz.volume()/4));
  REQUIRE(write_read_test(bz, "mp-917"));
  info_update("First Brillouin zone\n", fbz.python_string());
  info_update("Irreducible Brillouin zone\n", irp.python_string());
}

TEST_CASE("Irreducible Brillouin zone for mp-917","[bz_][materialsproject]"){
  // The spacegroup for CaC₂, from https://www.materialsproject.org/materials/mp-917/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-917.html
  // If brille::approx::scalar() is too strict this will fail.
  std::array<double, 9> lattice_vectors{
      7.068488010, 0.0,         0.002640550,
      0.0,         3.774933910, 0.0,
     -2.13512680,  0.0,         7.03226650};
  //
  std::string hall_symbol = "-C 2y";
  //
  auto lat = Direct<double>(lattice_vectors, MatrixVectors::row, hall_symbol);
  BrillouinZone bz(lat);
  REQUIRE(bz.check_ir_polyhedron());
  auto fbz = bz.get_polyhedron();
  auto irp = bz.get_ir_polyhedron();
  REQUIRE(irp.volume() == Approx(fbz.volume()/4));
  REQUIRE(write_read_test(bz, "mp-917"));
}

TEST_CASE("No irreducible Brillouin zone for inconsistent parameters and symmetry","[bz_]"){
  double a{3.5}, c{12.9}, alpha{90}, gamma{120};
  std::string spacegroup = "P 4";
  auto lat = Direct<double>({a,a,c}, {alpha,alpha,gamma}, spacegroup);
  REQUIRE_THROWS_WITH( BrillouinZone(lat), "Failed to find an irreducible Brillouin zone.");
}

TEST_CASE("Nb irreducible Brillouin Zone", "[bz_]"){
  // The conventional cell for Nb.
  std::string spg("-I 4 2 3");
  auto lat = Direct<double>({3.2598, 3.2598, 3.2598},
      {brille::math::half_pi, brille::math::half_pi, brille::math::half_pi}, spg);
  BrillouinZone bz(lat);
  auto fbz = bz.get_polyhedron();
  auto irbz = bz.get_ir_polyhedron();
  auto mult = bz.get_pointgroup_symmetry().size();
  REQUIRE(irbz.volume() == Approx(fbz.volume() / mult));
}