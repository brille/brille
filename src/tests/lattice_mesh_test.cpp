#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "lattice_mesh.hpp"

using namespace brille;

namespace {
bArray<double> at(const double x, const double y, const double z) {
  return bArray<double>::from_std(std::vector<std::array<double, 3>>{{x, y, z}});
}
double total(const std::vector<std::pair<ind_t, double>> & vw) {
  double t{0};
  for (const auto & [v, w]: vw) t += w;
  return t;
}
// the unit cube as the six Kuhn tetrahedra along its main diagonal
LatticeMesh cube() {
  std::vector<std::array<double, 3>> corners{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
  std::vector<std::array<ind_t, 4>> tets{{0, 1, 3, 7}, {0, 1, 5, 7}, {0, 2, 3, 7}, {0, 2, 6, 7}, {0, 4, 5, 7}, {0, 4, 6, 7}};
  return {bArray<double>::from_std(corners), bArray<ind_t>::from_std(tets)};
}
}

TEST_CASE("The lattice mesh locates points in constant time", "[lattice_mesh]") {
  const auto mesh = cube();
  REQUIRE_THAT(mesh.volume(), Catch::Matchers::WithinAbs(1.0, 1e-14));
  SECTION("inside, with weights that reproduce the point") {
    const std::array<double, 3> p{0.3, 0.6, 0.2};
    const auto vw = mesh.locate(at(p[0], p[1], p[2]));
    REQUIRE_FALSE(vw.empty());
    REQUIRE_THAT(total(vw), Catch::Matchers::WithinAbs(1.0, 1e-12));
    std::array<double, 3> q{0, 0, 0};
    for (const auto & [v, w]: vw) for (int i = 0; i < 3; ++i) q[i] += w * mesh.get_vertex_positions().val(v, i);
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(q[i], Catch::Matchers::WithinAbs(p[i], 1e-12));
  }
  SECTION("on a vertex, an edge and a shared face") {
    REQUIRE(mesh.locate(at(1, 1, 1)).size() == 1u);
    REQUIRE(mesh.locate(at(0.5, 0.5, 0.5)).size() == 2u);
    REQUIRE_THAT(total(mesh.locate(at(0.4, 0.4, 0.1))), Catch::Matchers::WithinAbs(1.0, 1e-12));
  }
  SECTION("outside by round-off") {
    for (const auto & p: {at(1 + 1e-13, 0.37, 0.61), at(0.42, -2e-14, 0.5), at(1 + 1e-14, 1 + 1e-14, 0.5)}) {
      const auto vw = mesh.locate(p);
      REQUIRE_FALSE(vw.empty());
      REQUIRE_THAT(total(vw), Catch::Matchers::WithinAbs(1.0, 1e-12));
      for (const auto & [v, w]: vw) REQUIRE(w > 0.);
    }
  }
  SECTION("outside") {
    REQUIRE(mesh.locate(at(1.001, 0.37, 0.61)).empty());
    REQUIRE(mesh.locate(at(-0.5, 0.5, 0.5)).empty());
  }
}

TEST_CASE("The lattice mesh fills a zone's irreducible part", "[lattice_mesh]") {
  using namespace brille::lattice;
  using namespace brille::math;
  // body-centred cubic Nb, given by its conventional cell
  std::array<double, 3> len{3.2598, 3.2598, 3.2598}, ang{half_pi, half_pi, half_pi};
  const BrillouinZone bz(Direct(len, ang, "-I 4 2 3"));
  const auto ir = bz.get_ir_polyhedron();
  const auto mesh = LatticeMesh::from_zone(bz, ir.volume() / 200);
  REQUIRE_THAT(mesh.volume(), Catch::Matchers::WithinRel(ir.volume(), 1e-10));
  REQUIRE_FALSE(mesh.refinement_limited());
  // every corner of the irreducible polyhedron is a mesh vertex (to round-off:
  // the zone and the mesh compute the corners differently)
  const auto corners = ir.vertices().xyz();
  for (ind_t i = 0; i < corners.size(0); ++i) {
    const auto vw = mesh.locate(corners.view(i));
    REQUIRE_FALSE(vw.empty());
    const auto top = *std::max_element(vw.begin(), vw.end(), [](const auto & a, const auto & b) { return a.second < b.second; });
    REQUIRE_THAT(top.second, Catch::Matchers::WithinAbs(1.0, 1e-9));
    for (int k = 0; k < 3; ++k)
      REQUIRE_THAT(mesh.get_vertex_positions().val(top.first, k), Catch::Matchers::WithinAbs(corners.val(i, k), 1e-12));
  }
}

#include "bz_mesh.hpp"
TEST_CASE("Refining a filled mesh appends the new vertices' data", "[lattice_mesh]") {
  using namespace brille::lattice;
  using namespace brille::math;
  std::array<double, 3> len{3.0, 3.0, 5.0}, ang{half_pi, half_pi, 2 * pi / 3};
  const BrillouinZone bz(Direct(len, ang, "-P 6"));
  BrillouinZoneMesh3<double, std::complex<double>, double> mesh(bz, bz.get_ir_polyhedron().volume() / 200);
  const auto nv = mesh.size();
  // one scalar per vertex: its index
  brille::Array<double> values(brille::shape_t{nv, 1u});
  for (ind_t i = 0; i < nv; ++i) values.val(brille::shape_t{i, 0u}) = static_cast<double>(i);
  brille::Array<std::complex<double>> vectors(brille::shape_t{nv, 1u}, std::complex<double>(0, 0));
  brille::Interpolator<double> iv(values, {1, 0, 0}, RotatesLike::vector, LengthUnit::real_lattice);
  brille::Interpolator<std::complex<double>> ie(vectors, {1, 0, 0}, RotatesLike::vector, LengthUnit::real_lattice);
  mesh.replace_data(iv, ie);
  std::vector<ind_t> where;
  for (ind_t t = 0; t < mesh.get_mesh_tetrehedra().size(0); t += 7) where.push_back(t);
  const auto planned = mesh.refinement_points(where, 0.0);
  const auto added = planned.size(0);
  REQUIRE(added > 0);
  brille::Array<double> more(brille::shape_t{added, 1u}, -1.0);
  brille::Array<std::complex<double>> more_vectors(brille::shape_t{added, 1u}, std::complex<double>(0, 0));
  const auto points = mesh.refine(where, 0.0, more, more_vectors);
  REQUIRE(points.size(0) == added);
  REQUIRE(mesh.size() == nv + added);
  REQUIRE(mesh.data().values().data().size(0) == nv + added);
  REQUIRE(mesh.data().values().data().val(0u, 0u) == 0.0);
  REQUIRE(mesh.data().values().data().val(nv, 0u) == -1.0);
}
