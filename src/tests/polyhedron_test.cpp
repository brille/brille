#include <catch2/catch.hpp>

#include "arrayvector.hpp"
#include "polyhedron.hpp"

TEST_CASE("Polyhedron instantiation","[polyhedron]"){
  std::vector<std::array<double,3>> va_verts{{1,1,0},{2,0,0},{1,1,1},{0,0,0}};
  ArrayVector<double> verts(va_verts);
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
  Polyhedron poly = Polyhedron(ArrayVector<double>(verts));
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
    auto poly = Polyhedron(ArrayVector<double>(v));
    auto doubled_poly = poly + poly.mirror();
    ArrayVector<double> n(3u,1u,0.), p(3u, 1u, 0.);
    SECTION("Hanging line"){
        n.insert(-1.0,0,0);
    }
    SECTION("Hanging triangular face"){
        n.insert(-1,0,2);
    }
    SECTION("Hanging rectangular face"){
        n.insert(-1.,0,0);
        n.insert(1.,0,1);
    }
    auto cut = Polyhedron::bisect(doubled_poly, n, p);
    REQUIRE(cut.get_volume() == Approx(poly.get_volume()));
    REQUIRE(cut.get_vertices().size() == 6u);
}

TEST_CASE("Polyhedron random point distribution","[polyhedron]"){
  std::vector<std::array<double,3>> va_verts{{1,1,0},{2,0,0},{1,1,1},{0,0,0}};
  ArrayVector<double> verts(va_verts);
  Polyhedron poly(verts);

  size_t seed=10, npoints=1000;
  ArrayVector<double> r0 = poly.rand_rejection(npoints, seed);
  // check that the generated points are inside the polyhedron
  REQUIRE(poly.contains(r0).count_true() == npoints);
  // that using the same seed reproduces the results
  ArrayVector<double> r1 = poly.rand_rejection(npoints, seed);
  REQUIRE((r1-r0).all_approx(0.));

  // and that specifying a clock-based seed does not give the same results
  r0 = poly.rand_rejection(npoints); // no seed specified == 0
  r1 = poly.rand_rejection(npoints);
  REQUIRE(!(r1-r0).all_approx(0.));
}
