#include <catch2/catch.hpp>
#include <string>
#include <array>
#include <vector>
#include <utility>

#include "symmetry.hpp"
#include "hall_symbol.hpp"

using namespace brille;

TEST_CASE("Symmetry CIF xyz format","[Symmetry]"){
  std::array<int,9> matrix{{1,0,0, 0,-1,0, 0,0,1}};
  std::array<double,3> vector{{0, 0.5, 0.5}};
  std::string xyz;

  SECTION("No trailing motion delimeter"){ xyz = "x,1/2-y,1/2+z";}
  SECTION("Trailing ; motion delimeter"){ xyz = "x,1/2-y,1/2+z;";}

  Motion<int,double> from_matrix_vector(matrix, vector);
  std::vector<Motion<int,double>> motions{from_matrix_vector};
  Symmetry from_motion(motions);

  Symmetry from_xyz(xyz);

  REQUIRE(from_xyz == from_motion);
}

TEST_CASE("Spacegroup CIF xyz format", "[Symmetry]"){
  // A random -P 6c (IT: P 6_3/m) CIF provides:
  // Fluorapatite pulled from https://aoterodelaroza.github.io/devnotes/space-group-collection/
  std::string xyz = "x,y,z;  -x,-y,-z;  x-y,x,-1/2+z;  -x+y,-x,-1/2-z;  -y,x-y,z;  y,-x+y,-z;  -x,-y,-1/2+z;  x,y,-1/2-z;  -x+y,-x,z;  x-y,x,-z;  y,-x+y,-1/2+z;  -y,x-y,-1/2-z";
  Symmetry from_xyz(xyz);
  // which *should* match the spacegroup we get from the HallSymbol:
  std::string hall_symbol = "-P 6c";
  HallSymbol hs(hall_symbol);
  Symmetry from_hall = hs.get_generators().generate();

  REQUIRE(from_xyz == from_hall);
}

