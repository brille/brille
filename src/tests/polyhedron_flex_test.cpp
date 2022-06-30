#include <catch2/catch.hpp>
#include <filesystem>

#include "array_.hpp" // defines bArray<T> as Array2<T> or Array<T>
#include "polyhedron_flex.hpp"
#include "process_id.hpp"

using namespace brille;
using namespace brille::polyhedron;

template<class T, class R, class S=std::common_type_t<T,R>>
static bool lists_match(const std::vector<std::vector<T>>& a, const std::vector<std::vector<R>>& b){
  if (a.size() != b.size()) return false;
  auto s = a.size();
  for (size_t i=0; i<s; ++i){
    const auto & ai{a[i]};
    const auto & bi{b[i]};
    if (ai.size() != bi.size()) return false;
    for (size_t j=0; j<ai.size(); ++j)
      if (static_cast<S>(ai[j]) != static_cast<S>(bi[j])) return false;
  }
  return true;
}

template<class T, class R, class S=std::common_type_t<T,R>>
static bool equivalent_permutations(const std::vector<T>& a, const std::vector<R>& b){
    if (a.size() != b.size()) return false;
    for (size_t roll=0; roll < a.size(); ++roll){
        bool ok{true};
        for (size_t i=0; i < a.size(); ++i) ok &= static_cast<S>(a[i]) == static_cast<S>(b[(i + roll) % a.size()]);
        if (ok) return ok;
    }
    return false;
}

TEST_CASE("Polyhedron instantiation","[polyhedron]"){
  std::vector<std::array<double,3>> va_verts{{1,1,0},{2,0,0},{1,1,1},{0,0,0}};
  auto verts = bArray<double>::from_std(va_verts);
  std::vector<std::vector<int>> vpf{{0,1,3},{0,2,1},{0,3,2},{1,2,3}};
  std::vector<std::array<double,3>> va_norms{{0,0,-1},{1,1,0},{-1,1,0},{0,-1,1}};
  auto norms = bArray<double>::from_std(va_norms);
  norms /= norm(norms);  // ensure the normal vectors have length unity
  auto checks = [verts, norms, vpf](const Poly<double,bArray>& p){
    // verify that the vertices are the same
    auto pv = p.vertices();
    for (const auto & i: verts.subItr()){
      REQUIRE(verts[i] == Approx(pv[i]));
    }
    // verify that the anticipated faces were found:
    size_t num_faces_matching{0};
    auto pf = p.faces();
    auto pff = pf.faces();
    auto pn = p.face_normals();
    for (size_t i=0; i < pff.size(); ++i){
      for (size_t j=0; j < vpf.size(); ++j){
        if (equivalent_permutations(pff[i], vpf[j])){
          ++num_faces_matching;
          REQUIRE(pn.view(i) == norms.view(j));
        }
      }
    }
    REQUIRE( num_faces_matching == vpf.size() );
    REQUIRE( p.volume() == Approx(2./6.)); // [(200)×(110)]⋅(111)/6
  };
  SECTION("Convex Hull creation"){
    auto poly = Poly<double,bArray>(verts);
    checks(poly);
  }
  SECTION("Vertices and Faces creation"){
    auto faces = Faces(vpf);
    auto poly = Poly<double,bArray>(verts, faces);
    checks(poly);
  }

}

TEST_CASE("Polyhedron intersection","[polyhedron]"){
  double a = 0.96373785;
  std::vector<std::array<double,3>> verts{{a,a,0},{2*a,0,0},{a,a,a},{0,0,0}};
  auto poly = Poly(bArray<double>::from_std(verts));
  REQUIRE(poly.volume() == Approx(2.*a*a*a/6.));
  double x = 0.143963;
  std::vector<std::array<double,3>> body_diagonal{{2*x, 0, 0}, {3*x, x, x}};
  auto box = bounding_box(bArray<double>::from_std(body_diagonal));
  REQUIRE(box.volume() == Approx(x*x*x));
  auto poly_box = poly.intersection(box);
  auto box_poly = box.intersection(poly);
  REQUIRE(poly_box.volume() == Approx(box.volume()/2));
  REQUIRE(box_poly.volume() == Approx(box.volume()/2));
}

TEST_CASE("Polyhedron bisect","[polyhedron]"){
  double a{0.457225}, c{0.753723};
  std::vector<std::array<double,3>> v{{a,0,0},{0,0,c},{a,a,c},{a,0,c},{a,a,0},{0,0,0}};
  auto poly = Poly(bArray<double>::from_std(v));
  auto doubled_poly = poly + poly.mirror();
  bArray<double> pa({1, 3}, 0.), pb({1, 3}, 0.), pc({1, 3}, 0.);
  pa.val(0, 1) = 1.0;
  SECTION("Hanging line"){
    // n = (-1 0 0) --> (0, 1, 0), (0, 0, 0), (0, 0, 1)
    pc.val(0, 2) = 1;
  }
  SECTION("Hanging triangular face"){
    // n = (0 0 -1) -> (0, 1, 0), (0, 0, 0), (-1, 0, 0)
    pc.val(0, 0) = -1;
  }
  SECTION("Hanging rectangular face"){
    // n = (-1, 1, 0) -> (1, 1, 0), (0, 0, 0), (0, 0, 1)
    pa.val(0, 0) = 1;
    pc.val(0, 2) = 1;
  }
  auto cut = doubled_poly.cut(pa, pb, pc);
  REQUIRE(cut.volume() == Approx(poly.volume()));
  REQUIRE(cut.vertices().size(0) == 6u);
}

TEST_CASE("Polyhedron random point distribution","[polyhedron]"){
  std::vector<std::array<double,3>> va_verts{{1,1,0},{2,0,0},{1,1,1},{0,0,0}};
  auto verts = bArray<double>::from_std(va_verts);
  auto poly = Poly(verts);

  ind_t seed=10, npoints=1000;
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
  auto poly = Poly(verts, Faces(verts_per_face));
  // And a cube near its corner:
  std::vector<std::array<double,3>> body = {{5.9764222286586582e-18, -0.047832811538213366, 0.53259981333751372},
                                            {0.084557263504132035, 0.047832811537652953, 0.62136644889376602}};
  auto cube = bounding_box(bArray<double>::from_std(body));
  // Must have an intersection no larger in volume than the cube
  auto cbp = cube.intersection(poly);
  REQUIRE(cbp.volume() == Approx(0.000366739));
  REQUIRE(cbp.volume() <= cube.volume());
  REQUIRE(cbp.vertices().size(0) == 10u);
}

TEST_CASE("Small face polyhedron convex hull","[polyhedron][convexhull]"){
  std::vector<std::array<double,3>> va_verts
  {{0.0846, 0.0217, 0.1775},
   {0.0846, 0.0652, 0.1775},
   {0.0846, 0.0217, 0.1506},
   {0.101 , 0.0217, 0.1775},
   {0.0846, 0.0652, 0.1757},
   {0.0857, 0.0652, 0.1775}};
  auto verts = bArray<double>::from_std(va_verts);
  auto cv_poly = Poly(verts);

  std::vector<std::vector<int>> vpf
  {{0,1,4,2},
   {0,3,5,1},
   {0,2,3},
   {2,4,5,3},
   {4,1,5}};
  auto fct_poly = Poly(verts, Faces(vpf));

  //REQUIRE(cv_poly.get_vertices_per_face().size() == fct_poly.get_vertices_per_face().size());
  REQUIRE(cv_poly.volume() == Approx(fct_poly.volume()));
}

TEST_CASE("Polyhedron IO","[polyhedron][io]"){
    std::vector<std::array<double,3>> va_verts{{1,1,0},{2,0,0},{1,1,1},{0,0,0}};
    auto verts = bArray<double>::from_std(va_verts);
    std::vector<std::vector<int>> vpf{{0,1,3},{0,2,1},{0,3,2},{1,2,3}};
    auto poly = Poly(verts, Faces(vpf));

#ifdef USE_HIGHFIVE
    // write the Polyhedron to the file:
    namespace fs=std::filesystem;
    auto tdir = fs::temp_directory_path();
    fs::path filepath = tdir;
    filepath /= fs::path(pid_filename("brille",".h5"));
    auto filename = filepath.string();
    std::string dataset="/polyhedron";
    REQUIRE(poly.to_hdf(filename, dataset));

    // read-back the file's Polyhedron
    auto read = Poly<double, bArray>::from_hdf(filename, dataset);

    // ensure that the two Polyhedra are identical
    auto rv = read.vertices();
    auto rp = read.face_points();
    auto rn = read.face_normals();
    auto rf = read.faces();
    auto rfpv = rf.faces_per_vertex();
    auto rvpf = rf.faces(); // vertex indexes for each face

    auto p = poly.face_points();
    auto n = poly.face_normals();
    auto fpv = poly.faces().faces_per_vertex();

    for (unsigned int i=0; i<2u; ++i){
      REQUIRE(verts.size(i) == rv.size(i));
      REQUIRE(p.size(i) == rp.size(i));
      REQUIRE(n.size(i) == rn.size(i));
    }
    for (const auto & s: verts.subItr()) REQUIRE(verts[s] == rv[s]);
    for (const auto & s: p.subItr()) REQUIRE(p[s] == rp[s]);
    for (const auto & s: n.subItr()) REQUIRE(n[s] == rn[s]);

    REQUIRE(lists_match(fpv, rfpv));
    REQUIRE(lists_match(vpf, rvpf));

    // (over)write to the file again, just to ensure it doesn't raise an error
    REQUIRE(poly.to_hdf(filename, dataset));

    fs::remove(filepath);
#endif //USE_HIGHFIVE
}

TEST_CASE("Plane convention conversion","[plane]"){
  std::vector<std::array<double, 3>> raw_v{
    {1, 0, 0}, {0, 1, 0}, {0, 1, 1}, {1, 1, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 1}, {1, 1, -2}, {1, -1, 0}, {2, 1, 0}
  };
  auto v = bArray<double>::from_std(raw_v);
  auto [a, b, c] = plane_points_from_normal(v, v);
  auto n = three_point_normal(a, b, c);
  for (ind_t i=0; i < v.size(0); ++i){
    REQUIRE( dot(v.view(i), n.view(i)).sum() == Approx((norm(v.view(i)) * norm(n.view(i))).sum()));
  }
}

TEST_CASE("Doubled face","[polyhedron]"){
  std::vector<std::array<double,3>> poly_vertices {
    { 0.00000000000000000000000000000000000000000000000000,  0.66666693318353020814015508221928030252456665039062,  0.49999999999999994448884876874217297881841659545898},
    { 0.00000000000000000000000000000000000000000000000000,  0.66666693318353020814015508221928030252456665039062, -0.49999999999999972244424384371086489409208297729492},
    { 0.57735050000000009973177839128766208887100219726562, -0.33333346659176532611468246614094823598861694335938,  0.49999999999999988897769753748434595763683319091797},
    { 0.57735050000000009973177839128766208887100219726562,  0.33333346659176521509238000362529419362545013427734,  0.49999999999999994448884876874217297881841659545898},
    { 0.57735050000000009973177839128766208887100219726562, -0.33333346659176521509238000362529419362545013427734, -0.49999999999999972244424384371086489409208297729492},
    { 0.57735050000000009973177839128766208887100219726562,  0.33333346659176477100317015356267802417278289794922, -0.49999999999999966693309261245303787291049957275391},
    {-0.57735050000000009973177839128766208887100219726562,  0.33333346659176521509238000362529419362545013427734,  0.49999999999999983346654630622651893645524978637695},
    {-0.57735050000000009973177839128766208887100219726562,  0.33333346659176510407007754110964015126228332519531, -0.49999999999999983346654630622651893645524978637695},
    {-0.00000000000000011102230246251565404236316680908203, -0.66666693318353020814015508221928030252456665039062,  0.49999999999999977795539507496869191527366638183594},
    {-0.57735050000000009973177839128766208887100219726562, -0.33333346659176499304777507859398610889911651611328,  0.49999999999999977795539507496869191527366638183594},
    {-0.00000000000000022204460492503130808472633361816406, -0.66666693318353020814015508221928030252456665039062, -0.49999999999999983346654630622651893645524978637695},
    {-0.57735050000000009973177839128766208887100219726562, -0.33333346659176499304777507859398610889911651611328, -0.49999999999999988897769753748434595763683319091797}
  };
  std::vector<std::vector<int>> poly_faces {
    {2, 3, 0, 6, 9, 8},
    { 1, 5, 4, 10, 11, 7 },
    { 4, 2, 8, 10 },
    { 0, 1, 7, 6 },
    { 0, 3, 5, 1 },
    { 3, 2, 4, 5 },
    { 6, 7, 11, 9 },
    { 8, 9, 11, 10 }
  };

  std::vector<std::array<double,3>> cube_vertices {
    {-0.23094020000000003989271135651506483554840087890625, -0.66666693318353020814015508221928030252456665039062, -0.49999999999999988897769753748434595763683319091797},
    {-0.23094020000000003989271135651506483554840087890625, -0.54545476351379740265201689908280968666076660156250, -0.49999999999999988897769753748434595763683319091797},
    {-0.23094020000000003989271135651506483554840087890625, -0.54545476351379740265201689908280968666076660156250, -0.37499999999999988897769753748434595763683319091797},
    {-0.23094020000000003989271135651506483554840087890625, -0.66666693318353020814015508221928030252456665039062, -0.37499999999999988897769753748434595763683319091797},
    {-0.11547010000000001994635567825753241777420043945312, -0.66666693318353020814015508221928030252456665039062, -0.49999999999999988897769753748434595763683319091797},
    {-0.11547010000000001994635567825753241777420043945312, -0.54545476351379740265201689908280968666076660156250, -0.49999999999999988897769753748434595763683319091797},
    {-0.11547010000000001994635567825753241777420043945312, -0.54545476351379740265201689908280968666076660156250, -0.37499999999999988897769753748434595763683319091797},
    {-0.11547010000000001994635567825753241777420043945312, -0.66666693318353020814015508221928030252456665039062, -0.37499999999999988897769753748434595763683319091797}
  };
  std::vector<std::vector<int>> cube_faces {
    {3, 0, 4, 7},
    { 3, 2, 1, 0 },
    { 0, 1, 5, 4 },
    { 3, 7, 6, 2 },
    { 7, 4, 5, 6 },
    { 2, 6, 5, 1 }
  };

  std::vector<std::array<double, 3>> node_vertices {
    {-0.11547010000000001994635567825753241777420043945312, -0.54545476351379740265201689908280968666076660156250, -0.49999999999999988897769753748434595763683319091797},
    {-0.11547010000000001994635567825753241777420043945312, -0.54545476351379740265201689908280968666076660156250, -0.37499999999999988897769753748434595763683319091797},
    {-0.11547010000000001994635567825753241777420043945312, -0.60000023986517725393952105150674469769001007080078, -0.49999999999999988897769753748434595763683319091797},
    {-0.20994563636363644532067951331555377691984176635742, -0.54545476351379740265201689908280968666076660156250, -0.49999999999999988897769753748434595763683319091797},
    {-0.11547010000000001994635567825753241777420043945312, -0.60000023986517725393952105150674469769001007080078, -0.37499999999999988897769753748434595763683319091797},
    {-0.20994563636363644532067951331555377691984176635742, -0.54545476351379740265201689908280968666076660156250, -0.37499999999999988897769753748434595763683319091797}
  };
  std::vector<std::vector<int>> node_faces {
    {0, 2, 3},
    { 1, 5, 4 },
    { 0, 1, 4, 2 },
    { 1, 0, 3, 5 },
    { 3, 2, 4, 5 }
  };

  auto poly = Poly(bArray<double>::from_std(poly_vertices), Faces(poly_faces));
  auto cube = Poly(bArray<double>::from_std(cube_vertices), Faces(cube_faces));
  auto node = Poly(bArray<double>::from_std(node_vertices), Faces(node_faces));
  auto res = cube.intersection(poly);
//  info_update("Poly\n", poly.python_string(), "\ncube\n", cube.python_string());
//  info_update("res\n", res.python_string(), "\nnode\n", node.python_string());
  REQUIRE(res == node);
  REQUIRE(node == res);
}

/*  The following test requires HDF5 support *and* an appropriately generated HDF5 file with test objects.
 *  Sorting-out file-loading from a known location is 'easy' enough, but when we only know, e.g., that the file
 *  is located somewhere relative in the git repository the problem seems insurmountable. Still, testing objects without
 *  relying on text-based serialisation seems like a worthwhile functionality to preserve. The following serves as a
 *  pseudo-template for how such testing *could* be achieved.
 */

//TEST_CASE("Doubled face from file","[polyhedron]"){
//  std::vector<std::array<double, 3>> node_vertices {
//    {-0.11547010000000001994635567825753241777420043945312, -0.54545476351379740265201689908280968666076660156250, -0.49999999999999988897769753748434595763683319091797},
//    {-0.11547010000000001994635567825753241777420043945312, -0.54545476351379740265201689908280968666076660156250, -0.37499999999999988897769753748434595763683319091797},
//    {-0.11547010000000001994635567825753241777420043945312, -0.60000023986517725393952105150674469769001007080078, -0.49999999999999988897769753748434595763683319091797},
//    {-0.20994563636363644532067951331555377691984176635742, -0.54545476351379740265201689908280968666076660156250, -0.49999999999999988897769753748434595763683319091797},
//    {-0.11547010000000001994635567825753241777420043945312, -0.60000023986517725393952105150674469769001007080078, -0.37499999999999988897769753748434595763683319091797},
//    {-0.20994563636363644532067951331555377691984176635742, -0.54545476351379740265201689908280968666076660156250, -0.37499999999999988897769753748434595763683319091797}
//  };
//  std::vector<std::vector<int>> node_faces {
//    {0, 2, 3},
//    { 1, 5, 4 },
//    { 0, 1, 4, 2 },
//    { 1, 0, 3, 5 },
//    { 3, 2, 4, 5 }
//  };
//  std::string filepath{"/home/.../file.h5"};
//
//  auto poly = Poly<double,bArray>::from_hdf(filepath, "/poly");
//  auto cube = Poly<double,bArray>::from_hdf(filepath, "/cube");
//  auto node = Poly(bArray<double>::from_std(node_vertices), node_faces);
//  auto res = cube.intersection(poly);
//  REQUIRE(res == node);
//  REQUIRE(node == res);
//}
