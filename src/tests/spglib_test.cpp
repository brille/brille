#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "hall_symbol.hpp"
#include "symmetry.hpp"
#include "spglib_dat.hpp"

using namespace brille;

namespace Catch {
  template<>
  struct StringMaker<brille::Symmetry> {
    static std::string convert( brille::Symmetry const& value ) {
      return value.to_string();
    }
  };
}

TEST_CASE("Check for consistency between Spglib sourced and HallSymbol decoded Symmetry objects","[spglib]"){
  // check every defined Hall group in the Spglib database
  auto hall_no = GENERATE(Catch::Generators::range(1, 531));
  // construct a Spacegroup object to have access to the Hall symbol associated with this number`
  Spacegroup spg(hall_no);
  // construct a HallSymbol object from the Spacegroup's Hall symbol
  HallSymbol hall_symbol(spg.get_hall_symbol());
  // convert the HallSymbol to a set of generators
  Symmetry hall_symbol_generators = hall_symbol.get_generators();
  // and use these to build a full spacegroup Symmetry
  Symmetry hall_symbol_sym = hall_symbol_generators.generate();
  // pull the spacegroup symmetries from the database
  Symmetry spglib_sym = make_spacegroup_symmetry_object(hall_no);
  // Now we can compare the Spglib sourced and decoded Symmetry objects
  if (hall_symbol_sym != spglib_sym){
    std::cout << "Failure for " << hall_no << ": " << spg.get_hall_symbol() << "\n";
  }
  REQUIRE(hall_symbol_sym == spglib_sym); // does not check Motion order
}
