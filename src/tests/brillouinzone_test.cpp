#include <random>
#include <chrono>
#include <catch2/catch.hpp>
#include "bz.hpp"
#include <filesystem>
#include "process_id.hpp"

using namespace brille;
using namespace brille::lattice;

#ifdef USE_HIGHFIVE
bool write_read_test(const BrillouinZone& source, const std::string& name){
    namespace fs = std::filesystem;
    auto temp_dir = fs::temp_directory_path();
    fs::path filepath = temp_dir;
    filepath /= fs::path(pid_filename("brille", ".h5"));

    // write the BrillouinZone to disk:
    auto wrote_ok = source.to_hdf(filepath.string(), name);
    if (!wrote_ok) {
      fs::remove(filepath);
      return false;
    }
    auto sink = BrillouinZone::from_hdf(filepath.string(), name);
    fs::remove(filepath);
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
  debug_update("First Brillouin zone\n", bz.get_polyhedron().python_string());
  debug_update("Irreducible Brillouin zone\n", bz.get_ir_polyhedron().python_string());
}


TEST_CASE("Primitive Hexagonal BrillouinZone instantiation","[bz_]"){
  auto lat = Direct<double>({3., 3., 3.}, {90., 90., 120.}, "P 1");
  BrillouinZone bz(lat);
  debug_update("First Brillouin zone\n", bz.get_polyhedron().python_string());
  debug_update("Irreducible Brillouin zone\n", bz.get_ir_polyhedron().python_string());
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
  debug_update("First Brillouin zone\n", fbz.python_string());
  debug_update("Irreducible Brillouin zone\n", irp.python_string());
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
//    REQUIRE(brille::approx_float::scalar(Qmq[i], static_cast<double>(tau[i])));
//    REQUIRE(brille::approx_float::scalar(std::abs(Qmqmtau[i]), 0.));
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
      REQUIRE(brille::approx_float::scalar(Qmq[i], static_cast<double>(tau[i])));
      REQUIRE(brille::approx_float::scalar(std::abs(Qmqmtau[i]), 0.));
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
  // moveinto is slower when multithreaded for some reason ...
  for (size_t i=1u; i<times.size(); ++i)
    REQUIRE(times[i] <= times[0]);
}

TEST_CASE("Irreducible Brillouin zone for mp-147 alt","[bz_][materialsproject]"){
  // The spacegroup for elemental Se, from https://www.materialsproject.org/materials/mp-147/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-147.html
  // If brille::approx_float::scalar() is too strict this will fail.
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
  debug_update("First Brillouin zone\n", fbz.python_string());
  debug_update("Irreducible Brillouin zone\n", irp.python_string());
}

TEST_CASE("Irreducible Brillouin zone for mp-147","[bz_][materialsproject][mp-147]"){
  auto run_tests = [&](BrillouinZone bz){
    REQUIRE(bz.check_ir_polyhedron());
    auto fbz = bz.get_polyhedron();
    auto irp = bz.get_ir_polyhedron();
    REQUIRE(irp.volume() == Approx(fbz.volume()/6));
    return bz;
  };
  // The spacegroup for elemental Se, from https://www.materialsproject.org/materials/mp-147/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-147.html
  // If brille::approx_float::scalar() is too strict this will fail.
  std::array<double, 9> lattice_vectors{
      10.699417459999999, 0.000000000000000, 0.000000000000000,
      -5.349708729999997, 9.265967326054772, 0.000000000000000,
      0.000000000000000, 0.000000000000000, 4.528988750000000};
  std::string hall_symbol = "-R 3";
  auto lat = Direct<double>(lattice_vectors, MatrixVectors::row, hall_symbol);
  SECTION("Use default floating point tolerance"){
    auto bz = run_tests(BrillouinZone(lat));
    // we only need to verify writing & reading for one version of the BZ
    REQUIRE(write_read_test(bz, "mp-147"));
  }
  SECTION("Use relaxed floating point tolerance"){
    auto ac = approx_float::Config().digit(10000); // increase the equivalency sloppiness from 1000 to 10000 * epsilon
    run_tests(BrillouinZone(lat, ac));
  }
}

TEST_CASE("Irreducible Brillouin zone for mp-306","[bz_][materialsproject]"){
  // The spacegroup for B₂O₃, from https://www.materialsproject.org/materials/mp-306/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-306.html
  // If brille::approx_float::scalar() is too strict this will fail.
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
  // If brille::approx_float::scalar() is too strict this will fail.
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
  // If brille::approx_float::scalar() is too strict this will fail.
  std::array<double, 9> lattice_vectors{
      3.559547360, 0.0,         0.0,
     -1.779773680, 3.082658440, 0.0,
      0.0,         0.0,        18.659215100000001};
  //
  std::string hall_symbol = "-R 3 2\"";
  //
  auto lat = Direct<double>(lattice_vectors, MatrixVectors::row, hall_symbol);
  auto ac = approx_float::Config().reciprocal(1e-10);
  BrillouinZone bz(lat, ac);
  REQUIRE(bz.check_ir_polyhedron());
  auto fbz = bz.get_polyhedron();
  auto irp = bz.get_ir_polyhedron();
  REQUIRE(irp.volume() == Approx(fbz.volume()/12));
  REQUIRE(write_read_test(bz, "mp-7041"));
}

TEST_CASE("Irreducible Brillouin zone for mp-917 atl","[bz_][materialsproject]"){
  // The spacegroup for CaC₂, from https://www.materialsproject.org/materials/mp-917/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-917.html
  // If brille::approx_float::scalar() is too strict this will fail.
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
  debug_update("First Brillouin zone\n", fbz.python_string());
  debug_update("Irreducible Brillouin zone\n", irp.python_string());
}

TEST_CASE("Irreducible Brillouin zone for mp-917","[bz_][materialsproject]"){
  // The spacegroup for CaC₂, from https://www.materialsproject.org/materials/mp-917/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-917.html
  // If brille::approx_float::scalar() is too strict this will fail.
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
  REQUIRE_THROWS( BrillouinZone(lat));
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


TEST_CASE("From python tests of all Hall symbols","[bz_]"){
  auto run_test = [](const std::string & hall) {
    auto lat = Direct<double>({5, 10, 15}, {90, 90, 90}, hall);
    BrillouinZone bz(lat);
    auto fbz = bz.get_polyhedron();
    auto irbz = bz.get_ir_polyhedron();
    auto mult = bz.get_pointgroup_symmetry().size();
    REQUIRE(irbz.volume() == Approx(fbz.volume() / mult));
  };
  // too-precise tolerances cause Convex Hull errors when the IR wedge is doubled
  // due to two sets of three on-face points giving normals (x, y, +ε) and (x, y, -ε)
  //
  SECTION("I 2 2"){run_test("I 2 2");}
  // lack of tolerance control in Poly::operator+ means we can identify co-planar faces
  // but can not do anything about it. This *could* be improved at some point?
  SECTION("C 2c 2"){run_test("C 2c 2");}
  SECTION("C 2 2"){run_test("C 2 2");}
  SECTION("F 2 2"){run_test("F 2 2");}
}

TEST_CASE("Primitive trigonal from python test of all Hall symbols", "[bz_]"){
  auto run_test = [](const std::string & hall) {
    auto lat = Direct<double>({5, 5, 5}, {60, 60, 60}, hall);
    BrillouinZone bz(lat);
    auto fbz = bz.get_polyhedron();
    auto irbz = bz.get_ir_polyhedron();
    auto mult = bz.get_pointgroup_symmetry().size();
    REQUIRE(irbz.volume() == Approx(fbz.volume() / mult));
  };
  // Failed previously due to attempting to use to points co-linear with the stationary axis to define
  // the wedge planes, which should have given 0/0 errors but did not due to rounding errors.
  // Now a check is made that the points are not parallel to the stationary axis
  SECTION("Hall -P 3*"){run_test("-P 3*");}
  // All wedge-finding attempts used the 2-fold axes first which produced bad results due to [UNKNOWN REASON]
  // compared to before. (Sorting instabilities caused by changes in the compiled binary? :-/)
  // Two new attempts are made after the other 1-6 have failed which turn off the special handling of 2-fold axes
  // *and* sort the operations so that 3+-fold operators are handled first.
  SECTION("Hall P 3* 2"){run_test("P 3* 2");}
}

TEST_CASE("Monoclinic from python test of all Hall symbols", "[bz_]"){
  auto run_test = [](const std::string & hall) {
    auto lat = Direct<double>({5, 10.0000000000000003, 5}, {90, 106.163894819260035, 90}, hall);
    BrillouinZone bz(lat);
    auto fbz = bz.get_polyhedron();
    auto irbz = bz.get_ir_polyhedron();
    auto mult = bz.get_pointgroup_symmetry().size();
    REQUIRE(irbz.volume() == Approx(fbz.volume() / mult));
  };
  SECTION("Hall -C 2yc"){run_test("-C 2yc");}

}


TEST_CASE("Rhombohedral from python test of all Hall symbols", "[bz_]"){
  auto run_test = [](const std::string & hall) {
    auto lat = Direct<double>({5, 5, 10}, {90, 90, 120}, hall);
    BrillouinZone bz(lat);
    auto fbz = bz.get_polyhedron();
    auto irbz = bz.get_ir_polyhedron();
    auto mult = bz.get_pointgroup_symmetry().size();
    REQUIRE(irbz.volume() == Approx(fbz.volume() / mult));
  };
  SECTION("Hall -R 3"){run_test("-R 3");}
}


TEST_CASE("Aflow lattices from python test", "[bz_][aflow]"){
  auto run_test = [](const std::array<double,3>& len,
                     const std::array<double,3>& ang,
                     const std::string & hall
                     ){
    auto lat = Direct(len, ang, hall);
//    info_update(lat.to_verbose_string());
//    info_update(lat.real_basis_vectors());
//    info_update(lat.reciprocal_basis_vectors());
    BrillouinZone bz(lat);
    auto fbz = bz.get_polyhedron();
//    info_update(fbz.python_string());
    auto irb = bz.get_ir_polyhedron();
    auto mlt = bz.get_pointgroup_symmetry().size();
    REQUIRE(irb.volume() == Approx(fbz.volume() / mlt));
  };
  // 579 of 582 passed the python tests, these three did not
  /* This R 3 (rhombohedrally-centred hexagonal) lattice has a primitive cell which is very nearly cubic
   * As a result the tolerances must be very small to capture all 24 vertices of the first Brillouin zone
   * */
  // One can replace c by sqrt(3)/sqrt(2) * a to make the primitive cell Cubic
  SECTION("Hall R 3"){run_test({6.885854708699999, 6.885854708699999, 8.4334670046}, {90, 90, 120}, "R 3");}
  /* These body-centered cubic lattices require more-thorough Convex hull creation */
  SECTION("Hall I -4 2 3, 1"){run_test({8.866400000000002, 8.866400000000002, 8.866400000000002}, {90, 90, 90}, "I -4 2 3");}
  SECTION("Hall I -4 2 3, 2"){run_test({8.911, 8.911, 8.911}, {90, 90, 90}, "I -4 2 3");}
  /* These face-centred tetragonal/orthorhombic systems exposed more problems with Convex Hull creation
   * including colinear points being used to define a hull bounding plane.*/
  SECTION("Hall F 2 2 2, 1"){run_test({6.479306892400001, 9.056209633700002, 10.0259106652}, {90, 90, 90}, "F 2 2 2");}
  SECTION("Hall F 2 2 2, 2"){run_test({5.540029163200001, 5.487028884200001, 5.195027347099999}, {90, 90, 90}, "F 2 2 2");}
}

TEST_CASE("Find limiting tolerance", "[.][bz_][aflow]"){
  auto lat = Direct<double>({6.885854708699999, 6.885854708699999, 8.4334670046}, {90, 90, 120}, "R 3");
  auto run_test = [&](const double tol=1e-10){
    auto ac = approx_float::Config(1000, 1e-10, tol);
    ind_t n{0};
    try {
      BrillouinZone bz(lat, ac);
      auto fbz = bz.get_polyhedron();
      n = fbz.vertices().size(0);
    } catch (const std::runtime_error &) {
      n = 0;
    } catch (...) {
      n = 1;
    }
    return n;
  };
  std::vector<double> powers{-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1};
  std::vector<ind_t> counts;
  for (auto p: powers){
    p = std::pow(10., p);
    counts.push_back(run_test(p));
  }
  for (ind_t i=0; i<powers.size(); ++i){
    info_update("10^", static_cast<int>(powers[i]), " : ", counts[i]);
  }

}



TEST_CASE("La2Zr2O7 BrillouinZone construction from off-symmetry basis vector input","[bz_][la2zr2o7][64][!mayfail][!nonportable]"){
  // Basis vectors from https://github.com/pace-neutrons/Euphonic/blob/aa3cc28786797bb3052f898dd63d4928d6f27ee2/tests_and_analysis/test/data/force_constants/LZO_force_constants.json#L186
  std::array<double,9> latmat {7.583912824349999, 1.8412792137035698e-32, 0.,
                               3.791956412170034, 3.791956412170034, 5.362636186024768,
                               3.791956412170034,-3.791956412170034, 5.362636186024768};
  // Symmetry information from CASTEP file via brilleu
  // the generators are: 4-fold [1 -1 1], 2-fold [-1 1 1], 3-fold [1 1 -3], -𝟙
  std::vector<std::array<int,9>> W {
      {{ 0,-1, 0,  1, 1, 1, -1, 0, 0}},
      {{ 0, 0, 1, -1,-1,-1,  1, 0, 0}},
      {{ 0, 0, 1,  1, 0, 0,  0, 1, 0}},
      {{-1, 0, 0,  0,-1, 0,  0, 0,-1}}
  };
  std::vector<std::array<double,3>> w{
      {{0.0, 0.5, 0.0}},
      {{0.0, 0.5, 0.0}},
      {{0.0, 0.0, 0.0}},
      {{0.0, 0.0, 0.0}}
  };

  Symmetry::Motions mots;
  mots.reserve(W.size());
  for (size_t i=0; i<W.size(); ++i) mots.push_back(Motion<int,double>(W[i], w[i]));
  Symmetry sym(mots);
  // Without 'snap_to_symmetry' the basis vectors are off on the order of 2e-12
  auto wrong = Direct<double>(latmat, MatrixVectors::row, sym, Basis(), false);

  SECTION("Default tolerance, wrong basis vectors (no exception thrown on macOS)"){
  // Without specifying a larger-than-normal tolerance the BrillouinZone
  // One of two things happens on different OS/machines:
  // linux & windows: construction fails and a Runtime Error is thrown
  // macOS: construction succeeds, no error is thrown
  // For this reason we CHECK instead of REQUIRE:
  CHECK_THROWS_AS(BrillouinZone(wrong), std::runtime_error);
  }
  SECTION("Relaxed tolerance, wrong basis vectors"){
  // With a non-standard tolerance BrillouinZone construction succeeds
  auto ac = approx_float::Config().reciprocal(2e-12);
  REQUIRE_NOTHROW(BrillouinZone(wrong, ac));
  }
  SECTION("Default tolerance, corrected basis vectors"){
  // Or turning-on 'snap_to_symmetry' corrects the basis vectors to follow
  // the provided symmetry operations
  auto lat = Direct<double>(latmat, MatrixVectors::row, sym, Basis(), true);
  // Such that the BrillouinZone construction succeeds with standard tolerances
  REQUIRE_NOTHROW(BrillouinZone(lat));
  }
}
