#include <catch2/catch.hpp>
#include <string>
#include <vector>
#include <utility>

#include "hall_symbol.hpp"
#include "spg_database.hpp"

using namespace brille;

TEST_CASE("Ensure HallSymbol and Symmetry(xyz string) concur","[HallSymbol]"){
  std::vector<std::pair<std::string, std::string>> hall_xyz;
  hall_xyz.push_back(std::make_pair("P -2zb","x, y+1/2, -z"));
  for (auto hs: hall_xyz){
    HallSymbol h(hs.first);
    Symmetry s(hs.second);
    // Verify equivalent motions exist in both generated Symmetry groups
    // Order of the Motions is not important.
    REQUIRE( h.get_generators().generate() == s.generate() );
  }
}

TEST_CASE("HallSymbols validate works","[HallSymbol]"){
  std::vector<std::string> strs{"x, y+1/2, -z","P -2zb"};
  for (auto s: strs){
  HallSymbol hs(s);
    if (hs.validate())
      REQUIRE_NOTHROW(hs.get_generators());
    else
      REQUIRE_THROWS(hs.get_generators());
  }
}
