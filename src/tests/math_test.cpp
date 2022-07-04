#include <catch2/catch.hpp>
#include "math.hpp"
#include "debug.hpp"

using namespace brille;
using namespace brille::math;

TEMPLATE_TEST_CASE("sine","[math]", float, double){
  TestType two{2}, three{3};
  TestType s22{std::sqrt(two)/two}, s32{std::sqrt(three)/two};
  auto eps = std::numeric_limits<TestType>::epsilon();
  std::vector<TestType> degrees{-360, -330, -315, -300, -270, -240, -225, -210,
                                -180, -150, -135, -120,  -90,  -60,  -45,  -30,
                                   0,   30,   45,   60,   90,  120,  135,  150,
                                 180,  210,  225,  240,  270,  300,  315,  330};
  std::vector<TestType>   sines{  0.,  0.5,  s22,  s32,   1.,  s32,  s22,  0.5,
                                  0., -0.5, -s22, -s32,  -1., -s32, -s22, -0.5,
                                  0.,  0.5,  s22,  s32,   1.,  s32,  s22,  0.5,
                                  0., -0.5, -s22, -s32,  -1., -s32, -s22, -0.5};
  for (size_t i=0; i<degrees.size(); ++i){
    auto calc = sind(degrees[i]);
    auto comparison = std::abs(calc + sines[i]) * eps;
    auto diff = std::abs(calc - sines[i]);
    info_update_if(diff > comparison, degrees[i], " => ",
                   my_to_string(calc), " != ", my_to_string(sines[i]),
                   " due to ",
                   my_to_string(diff), " > ", my_to_string(comparison));
    REQUIRE(diff <= comparison);
  }
}

TEMPLATE_TEST_CASE("cosine","[math]", float, double){
  TestType two{2}, three{3};
  TestType s22{std::sqrt(two)/two}, s32{std::sqrt(three)/two};
  auto eps = std::numeric_limits<TestType>::epsilon();
  std::vector<TestType> degrees{-360, -330, -315, -300, -270, -240, -225, -210,
                                -180, -150, -135, -120,  -90,  -60,  -45,  -30,
                                   0,   30,   45,   60,   90,  120,  135,  150,
                                 180,  210,  225,  240,  270,  300,  315,  330};
  std::vector<TestType> cosines{  1.,  s32,  s22,  0.5,   0., -0.5, -s22, -s32,
                                 -1., -s32, -s22, -0.5,   0.,  0.5,  s22,  s32,
                                  1.,  s32,  s22,  0.5,   0., -0.5, -s22, -s32,
                                 -1., -s32, -s22, -0.5,   0.,  0.5,  s22,  s32};
  for (size_t i=0; i<degrees.size(); ++i){
    auto calc = cosd(degrees[i]);
    auto comparison = std::abs(calc + cosines[i]) * eps;
    auto diff = std::abs(calc - cosines[i]);
    info_update_if(diff > comparison, degrees[i], " => ",
                   my_to_string(calc), " != ", my_to_string(cosines[i]),
                   " due to ",
                   my_to_string(diff), " > ", my_to_string(comparison));
    REQUIRE(diff <= comparison);
  }
}