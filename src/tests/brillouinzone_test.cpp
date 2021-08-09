#include <random>
#include <chrono>
#include <catch2/catch.hpp>
#include "bz.hpp"

using namespace brille;

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
  LQVec<double> Q(r, nQ); // holds an (nQ, 3) shaped Array or Array2
  for (auto& v: Q.valItr()) v = distribution(generator);

  LQVec<double> q(r,nQ);
  LQVec<int> tau(r,nQ);

  REQUIRE( bz.moveinto(Q,q,tau) ); // success indicated by return of true
  auto inside = bz.isinside(q);
  REQUIRE(std::count(inside.begin(), inside.end(), false) == 0);
  LQVec<double> Qmq = Q-q;
  LQVec<double> Qmqmtau = Q-(q+tau);
  for (auto i: Q.subItr()){
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
  LQVec<double> Q(r, bArray<double>::from_std(rawQ));
  LQVec<double> q(r,Q.size(0));
  LQVec<int> tau(r,Q.size(0));
  REQUIRE(bz.moveinto(Q,q,tau));
}

TEST_CASE("BrillouinZone moveinto hexagonal extended","[brillouinzone][moveinto][.timing]"){
  Direct dlat(3.,3.,9.,90.,90.,120.,1);
  Reciprocal rlat(dlat.star());
  BrillouinZone bz(rlat);
  std::vector<std::array<double,3>> rawQ;

  std::default_random_engine g(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> d(-4.0, 4.0);
  auto random_array = [&]() -> std::array<double,3> { return {d(g), d(g), d(g)}; };

  size_t count{0}, total_points{10000};
  while (++count < total_points)
    rawQ.push_back(random_array());

  LQVec<double> Q(rlat, bArray<double>::from_std(rawQ));
  LQVec<double> q(rlat, Q.size(0));
  LQVec<int> tau(rlat, Q.size(0));
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
    info_update("brillouinzone::moveinto of ",total_points," points performed by ",threads, " threads in ",timer.average(),"+/-",timer.jitter()," msec");
  }
  for (size_t i=1u; i<times.size(); ++i)
    REQUIRE(times[i] <= times[0]);
}

TEST_CASE("Irreducible Brillouin zone for mp-147","[brillouinzone][materialsproject]"){
  // The spacegroup for elemental Se, from https://www.materialsproject.org/materials/mp-147/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-147.html
  // If brille::approx::scalar() is too strict this will fail.
  double lattice_vectors[9] = {
    10.699417459999999, 0.000000050000000, 0.000000000000000,
    -5.349708760000000, 9.265967300000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 4.528988750000000};
  //
  std::string hall_symbol = "-R 3";
  //
  Direct dlat(lattice_vectors, hall_symbol);
  Reciprocal rlat = dlat.star();
  BrillouinZone bz(rlat);
  REQUIRE(bz.check_ir_polyhedron());
  Polyhedron fbz = bz.get_polyhedron();
  Polyhedron irp = bz.get_ir_polyhedron();
  REQUIRE(irp.get_volume() == Approx(fbz.get_volume()/6));
}

TEST_CASE("Irreducible Brillouin zone for mp-306","[brillouinzone][materialsproject]"){
  // The spacegroup for B₂O₃, from https://www.materialsproject.org/materials/mp-306/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-306.html
  // If brille::approx::scalar() is too strict this will fail.
  double lattice_vectors[9] = {
     4.350620350000000, 0.000000000000000, 0.000000000000000,
    -2.175310180000000, 3.767747740000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 8.339814740000000};
  //
  std::string hall_symbol = "P 31 2\"";
  //
  Direct dlat(lattice_vectors, hall_symbol);
  Reciprocal rlat = dlat.star();
  BrillouinZone bz(rlat);
  REQUIRE(bz.check_ir_polyhedron());
  Polyhedron fbz = bz.get_polyhedron();
  Polyhedron irp = bz.get_ir_polyhedron();
  REQUIRE(irp.get_volume() == Approx(fbz.get_volume()/6));
}

TEST_CASE("Irreducible Brillouin zone for mp-661","[brillouinzone][materialsproject]"){
  // The spacegroup for AlN, from https://www.materialsproject.org/materials/mp-661/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-661.html
  // If brille::approx::scalar() is too strict this will fail.
  double lattice_vectors[9] = {
     3.113692560000000, 0.000000000000000, 0.000000000000000,
    -1.556846270000000, 2.696536850000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 4.983221820000000};
  //
  std::string hall_symbol = "P 6c -2c";
  //
  Direct dlat(lattice_vectors, hall_symbol);
  Reciprocal rlat = dlat.star();
  BrillouinZone bz(rlat);
  REQUIRE(bz.check_ir_polyhedron());
  Polyhedron fbz = bz.get_polyhedron();
  Polyhedron irp = bz.get_ir_polyhedron();
  REQUIRE(irp.get_volume() == Approx(fbz.get_volume()/12));
}

TEST_CASE("Irreducible Brillouin zone for mp-7041","[brillouinzone][materialsproject]"){
  // The spacegroup for CaHgO₂, from https://www.materialsproject.org/materials/mp-7041/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-7041.html
  // If brille::approx::scalar() is too strict this will fail.
  double lattice_vectors[9] = {
      3.559547360, 0.0,         0.0,
     -1.779773680, 3.082658440, 0.0,
      0.0,         0.0,        18.659215100000001};
  //
  std::string hall_symbol = "-R 3 2\"";
  //
  Direct dlat(lattice_vectors, hall_symbol);
  Reciprocal rlat = dlat.star();
  BrillouinZone bz(rlat);
  REQUIRE(bz.check_ir_polyhedron());
  Polyhedron fbz = bz.get_polyhedron();
  Polyhedron irp = bz.get_ir_polyhedron();
  REQUIRE(irp.get_volume() == Approx(fbz.get_volume()/12));
}

TEST_CASE("Irreducible Brillouin zone for mp-917","[brillouinzone][materialsproject]"){
  // The spacegroup for CaC₂, from https://www.materialsproject.org/materials/mp-917/
  // via http://phonondb.mtl.kyoto-u.ac.jp/ph20180417/d000/mp-917.html
  // If brille::approx::scalar() is too strict this will fail.
  double lattice_vectors[9] = {
      7.068488010, 0.0,         0.002640550,
      0.0,         3.774933910, 0.0,
     -2.13512680,  0.0,         7.03226650};
  //
  std::string hall_symbol = "-C 2y";
  //
  Direct dlat(lattice_vectors, hall_symbol);
  Reciprocal rlat = dlat.star();
  BrillouinZone bz(rlat);
  REQUIRE(bz.check_ir_polyhedron());
  Polyhedron fbz = bz.get_polyhedron();
  Polyhedron irp = bz.get_ir_polyhedron();
  REQUIRE(irp.get_volume() == Approx(fbz.get_volume()/4));
}

TEST_CASE("No irreducible Brillouin zone for inconsistent parameters and symmetry","[brillouinzone]"){
  double a{3.5}, c{12.9}, alpha{90}, gamma{120};
  std::string spacegroup = "P 4";
  Direct dlat(a,a,c, alpha,alpha,gamma, spacegroup);
  Reciprocal rlat = dlat.star();
  REQUIRE_THROWS_WITH( BrillouinZone(rlat), "Failed to find an irreducible Brillouin zone."); 
}
