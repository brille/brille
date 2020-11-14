#include<memory>
#include<array>
#include <catch2/catch.hpp>

#include "munkres.hpp"

TEST_CASE("Munkres' Assignment algorithm","[munkres]"){
  // double cost[9] = {1.,2.,3., 2.,4.,6., 3.,6.,9.};
  std::vector<double> cost = {1.,2.,3., 2.,4.,6., 3.,6.,9.};

  size_t expected_assignment[3] = {2,1,0};
  size_t found_assignment[3];
  brille::assignment::Munkres<double> m_one(3,cost);
  // check that the assignment was made, and get the result
  REQUIRE( m_one.get_assignment(found_assignment));
  for (size_t i=0; i<3u; ++i)
    REQUIRE( found_assignment[i] == expected_assignment[i] );

  brille::assignment::Munkres<double> m_two(3);
  for (size_t i=0; i<9u; ++i) m_two.get_cost()[i] = cost[i];
  for (size_t i=0; i<3u; ++i) found_assignment[i] = 0u;

  REQUIRE( m_two.get_assignment(found_assignment) );
  for (size_t i=0; i<3u; ++i)
    REQUIRE( found_assignment[i] == expected_assignment[i] );
}
