#include <catch2/catch.hpp>

#include "lattice_dual.hpp"
#include "basis.hpp"
#include "symmetry.hpp"

using namespace brille::lattice;
using namespace brille::math;

TEST_CASE("CaHgO2 Lattice","[lattice][basis][issue63]"){
  Lattice<double>::matrix_t vectors {{
      1.7797736800000001, 1.027552813333333, 6.219738366666666,
     -1.7797736800000001, 1.027552813333333, 6.219738366666666,
      0.0,               -2.055105626666667, 6.219738366666666
    }};
  std::vector<std::array<double,3>> wrong_positions {
      {{0.89401258, 0.89401259, 0.89401258}},
      {{0.10598742, 0.10598741, 0.10598742}},
      {{0.5, 0.5, 0.5}},
      {{0., 0.99999999, 0.}}
  };
  std::vector<brille::ind_t> types {0, 0, 1, 2};

  brille::Basis basis(wrong_positions, types);

  std::vector<brille::Motion<int, double>> gen_M {
      {{-1,  0,  0,  0, -1,  0,  0,  0, -1}, {0., 0., 0.}},
      {{ 0,  0,  1,  1,  0,  0,  0,  1,  0}, {0., 0., 0.}},
      {{ 0,  0, -1,  0, -1,  0, -1,  0,  0}, {0., 0., 0.}}
  };
  brille::Symmetry gens(gen_M);
  auto symmetry = gens.generate();

  auto bmv = brille::MatrixVectors::row;

  // If false provided to snap_to_symmetry the Lattice basis matches the input
  auto bad_lat = Direct<double>(vectors, bmv,  symmetry, basis, false);
  auto bad_basis = bad_lat.basis();
  for (size_t i=0; i<4; ++i){
    for (size_t j=0; j<3; ++j){
      REQUIRE(wrong_positions[i][j] == bad_basis.position(i)[j]);
    }
  }

  // If true is provided for `snap_to_symmetry` the basis is modified
  std::vector<std::array<double,3>> right_positions {
      {{0.8940125833333333, 0.8940125833333333, 0.8940125833333333}},
      {{0.1059874166666667, 0.1059874166666666, 0.1059874166666667}},
      {{0.5, 0.5, 0.5}},
      {{0., 1., 0.}}
  };

  auto good_lat = Direct<double>(vectors, bmv, symmetry, basis, true);
  for (size_t i=0; i<4; ++i){
    for (size_t j=0; j<3; ++j){
      REQUIRE(brille::approx_float::scalar(right_positions[i][j], good_lat.basis().position(i)[j]));
    }
  }
}
