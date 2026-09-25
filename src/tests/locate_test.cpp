#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mesh.hpp" // triangulation_layers.hpp needs its includes

using namespace brille;

TEST_CASE("Points round-off outside the mesh are located","[locate]"){
  // a unit cube meshed in three layers
  std::vector<std::array<double,3>> corners{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
  auto vertices = bArray<double>::from_std(corners);
  std::vector<std::vector<int>> faces{{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}};
  auto mesh = triangulate(vertices, faces, 0.01, 3);
  auto at = [](const double x, const double y, const double z){
    return bArray<double>::from_std(std::vector<std::array<double,3>>{{x, y, z}});
  };
  auto total = [](const auto & vw){
    double t{0};
    for (const auto & [v, w]: vw) t += w;
    return t;
  };
  SECTION("inside"){
    auto vw = mesh.locate(at(0.3, 0.6, 0.2));
    REQUIRE_FALSE(vw.empty());
    REQUIRE_THAT(total(vw), Catch::Matchers::WithinAbs(1.0, 1e-12));
  }
  SECTION("outside by round-off"){
    for (const auto & p: {at(1 + 1e-13, 0.37, 0.61), at(0.42, -2e-14, 0.5), at(1 + 1e-14, 1 + 1e-14, 0.5)}){
      auto vw = mesh.locate(p);
      REQUIRE_FALSE(vw.empty());
      REQUIRE_THAT(total(vw), Catch::Matchers::WithinAbs(1.0, 1e-12));
      for (const auto & [v, w]: vw) REQUIRE(w > 0.);
    }
  }
  SECTION("outside"){
    REQUIRE(mesh.locate(at(1.001, 0.37, 0.61)).empty());
    REQUIRE(mesh.locate(at(-0.5, 0.5, 0.5)).empty());
  }
}
