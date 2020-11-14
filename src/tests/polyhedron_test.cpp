#include <catch2/catch.hpp>

#include "array_latvec.hpp" // defines bArray<T> as Array2<T> or Array<T>
#include "polyhedron.hpp"

using namespace brille;

TEST_CASE("Polyhedron instantiation","[polyhedron]"){
  std::vector<std::array<double,3>> va_verts{{1,1,0},{2,0,0},{1,1,1},{0,0,0}};
  auto verts = bArray<double>::from_std(va_verts);
  std::vector<std::vector<int>> vpf{{0,1,3},{0,2,1},{0,3,2},{1,2,3}};
  Polyhedron poly;
  SECTION("Convex Hull creation"){
    poly = Polyhedron(verts);
  }
  SECTION("Vertices and Faces creation"){
    poly = Polyhedron(verts, vpf);
  }
  // verify that the anticipated faces were found:
  size_t num_faces_matching{0};
  for (auto cf: poly.get_vertices_per_face()){
    for (auto f: vpf){
      if (cf.size()==f.size() &&
         ((cf[0]==f[0] && cf[1]==f[1] && cf[2]==f[2]) ||
          (cf[0]==f[1] && cf[1]==f[2] && cf[2]==f[0]) ||
          (cf[0]==f[2] && cf[1]==f[0] && cf[2]==f[1])   )
       ) ++num_faces_matching;
    }
  }
  REQUIRE( num_faces_matching == vpf.size() );
  REQUIRE( poly.get_volume() == Approx(2./6.)); // [(200)×(110)]⋅(111)/6
}

TEST_CASE("Polyhedron intersection","[polyhedron]"){
  double a = 0.96373785;
  std::vector<std::array<double,3>> verts{{a,a,0},{2*a,0,0},{a,a,a},{0,0,0}};
  Polyhedron poly = Polyhedron(bArray<double>::from_std(verts));
  REQUIRE(poly.get_volume() == Approx(2.*a*a*a/6.));
  double x = 0.143963;
  std::array<double,3> boxmin{2*x,0,0}, boxmax{3*x,x,x};
  Polyhedron box = polyhedron_box(boxmin,boxmax);
  Polyhedron poly_box = poly.intersection(box);
  Polyhedron box_poly = box.intersection(poly);
  REQUIRE(poly_box.get_volume() == Approx(box.get_volume()/2));
  REQUIRE(box_poly.get_volume() == Approx(box.get_volume()/2));
}

TEST_CASE("Polyhedron bisect","[polyhedron]"){
  double a{0.457225}, c{0.753723};
  std::vector<std::array<double,3>> v{{a,0,0},{0,0,c},{a,a,c},{a,0,c},{a,a,0},{0,0,0}};
  auto poly = Polyhedron(bArray<double>::from_std(v));
  auto doubled_poly = poly + poly.mirror();
  bArray<double> n({1,3}, 0.), p({1,3}, 0.);
  SECTION("Hanging line"){
    n.val(0,0) = -1;
  }
  SECTION("Hanging triangular face"){
    n.val(0,2) = -1;
  }
  SECTION("Hanging rectangular face"){
    n.val(0,0) = -1;
    n.val(0,1) = 1;
  }
  auto cut = Polyhedron::bisect(doubled_poly, n, p);
  REQUIRE(cut.get_volume() == Approx(poly.get_volume()));
  REQUIRE(cut.get_vertices().size(0) == 6u);
}

TEST_CASE("Polyhedron random point distribution","[polyhedron]"){
  std::vector<std::array<double,3>> va_verts{{1,1,0},{2,0,0},{1,1,1},{0,0,0}};
  auto verts = bArray<double>::from_std(va_verts);
  Polyhedron poly(verts);

  size_t seed=10, npoints=1000;
  auto r0 = poly.rand_rejection(npoints, seed);
  // check that the generated points are inside the polyhedron
  auto in = poly.contains(r0);
  size_t nin = std::count(in.begin(), in.end(), true);
  REQUIRE(nin == npoints);
  // that using the same seed reproduces the results
  auto r1 = poly.rand_rejection(npoints, seed);
  REQUIRE((r1-r0).all(0.,0));

  // and that specifying a clock-based seed does not give the same results
  r0 = poly.rand_rejection(npoints); // no seed specified == 0
  r1 = poly.rand_rejection(npoints);
  REQUIRE(!(r1-r0).all(0.,0));
}

TEST_CASE("Polyhedron intersection La2Zr2O7","[polyhedron][intersection]"){
  // The La2Zr2O7 irreducible Brillouin zone:
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
  // And a cube near its corner:
  std::array<double,3> min {{5.9764222286586582e-18, -0.047832811538213366, 0.53259981333751372}};
  std::array<double,3> max {{0.084557263504132035, 0.047832811537652953, 0.62136644889376602}};
  Polyhedron cube = polyhedron_box(min, max);
  // Must have an intersection no larger in volume than the cube
  Polyhedron cbp = cube.intersection(poly);
  REQUIRE(cbp.get_volume() == Approx(0.000366739));
  REQUIRE(cbp.get_volume() <= cube.get_volume());
  REQUIRE(cbp.get_vertices().size(0) == 10u);
}

TEST_CASE("Small face polyhedron convex hull","[!shouldfail][polyhedron][convexhull]"){
  std::vector<std::array<double,3>> va_verts
  {{0.0846, 0.0217, 0.1775},
   {0.0846, 0.0652, 0.1775},
   {0.0846, 0.0217, 0.1506},
   {0.101 , 0.0217, 0.1775},
   {0.0846, 0.0652, 0.1757},
   {0.0857, 0.0652, 0.1775}};
  auto verts = bArray<double>::from_std(va_verts);
  Polyhedron cv_poly(verts);

  std::vector<std::vector<int>> vpf
  {{0,1,4,2},
   {0,3,5,1},
   {0,2,3},
   {2,4,5,3},
   {4,1,5}};
  Polyhedron fct_poly(verts, vpf);

  //REQUIRE(cv_poly.get_vertices_per_face().size() == fct_poly.get_vertices_per_face().size());
  REQUIRE(cv_poly.get_volume() == Approx(fct_poly.get_volume()));
}
