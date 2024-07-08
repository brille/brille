#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEMPLATE_TEST_CASE("catch_test","[catch][macos-arm]", float, double){
  /* *
   * This is a simple test to ensure that the catch2 testing framework is working
   * as expected on the macOS ARM platform, for some small subset of suspicious cases.
   */
  TestType a{1}, b{2}, tol{0.0001};
  REQUIRE_THAT(a+b, Catch::Matchers::WithinRel(TestType(3), tol));
  REQUIRE_THAT(a-b, Catch::Matchers::WithinRel(TestType(-1), tol));
  REQUIRE_THAT(a*b, Catch::Matchers::WithinRel(TestType(2), tol));
  REQUIRE_THAT(a/b, Catch::Matchers::WithinRel(TestType(0.5), tol));
  REQUIRE_THAT(0*a - tol, Catch::Matchers::WithinAbs(TestType(0), tol));
  REQUIRE_THAT(0*a + tol, Catch::Matchers::WithinAbs(TestType(-0.0), tol));
  REQUIRE_THAT(-0.0 - tol, Catch::Matchers::WithinAbs(TestType(0.0), tol));
  REQUIRE_THAT(-0.0 + tol, Catch::Matchers::WithinAbs(TestType(-0.0), tol));
}