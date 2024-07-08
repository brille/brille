#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>


TEMPLATE_TEST_CASE("catch_test","[catch][macos-arm]", float, double){
  /* *
   * This is a simple test to ensure that the catch2 testing framework is working
   * as expected on the macOS ARM platform, for some small subset of suspicious cases.
   */
  TestType a{1}, b{2};
  REQUIRE(a+b == Approx(TestType(3)));
  REQUIRE(a-b == Approx(TestType(-1)));
  REQUIRE(a*b == Approx(TestType(2)));
  REQUIRE(a/b == Approx(TestType(0.5)));
  REQUIRE(0*a == Approx(TestType(0)));
  REQUIRE(0*a == Approx(TestType(-0.0)));
  REQUIRE(-0.0 == Approx(TestType(0.0)));
  REQUIRE(-0.0 == Approx(TestType(-0.0)));

}