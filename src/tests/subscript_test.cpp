#include <catch2/catch.hpp>
#include "subscript.hpp"
#include "utilities.hpp"

using namespace brille;

TEST_CASE("Full array iterator","[subscript]"){
  std::vector<std::vector<int>> expected{
    {0,0,0},{0,0,1},{0,0,2},
    {0,1,0},{0,1,1},{0,1,2},
    {0,2,0},{0,2,1},{0,2,2},
    {1,0,0},{1,0,1},{1,0,2},
    {1,1,0},{1,1,1},{1,1,2},
    {1,2,0},{1,2,1},{1,2,2},
  };
  SubIt<int> sit({2,3,3});
  size_t n{0};
  for (auto x : sit){
    REQUIRE(n < expected.size());
    REQUIRE(x.size() == expected[n].size());
    REQUIRE(brille::approx::vector(x.size(), x.data(), expected[n++].data()));
  }
  REQUIRE( expected.size() == n);
}

TEST_CASE("Fixed-index array iterator", "[subscript]"){
  std::vector<std::vector<int>> expected{
    {0,4,0},{0,4,1},{0,4,2},
    {1,4,0},{1,4,1},{1,4,2},
    {2,4,0},{2,4,1},{2,4,2},
  };
  SubIt<int> sit({3,5,3},{3,4,3});
  size_t n{0};
  for (auto x : sit){
    REQUIRE(n < expected.size());
    REQUIRE(x.size() == expected[n].size());
    REQUIRE(brille::approx::vector(x.size(), x.data(), expected[n++].data()));
  }
  REQUIRE( expected.size() == n);
}

TEST_CASE("Auto broadcasting","[subscript]"){
  using namespace brille;
  std::vector<int> a{3,4}, b;
  std::vector<std::vector<int>>
    o_expected{{0,0},{0,1},{0,2},{0,3},{1,0},{1,1},{1,2},{1,3},{2,0},{2,1},{2,2},{2,3}},
    b_expected;
  SECTION("row vector"){
    b={1,4};
    b_expected = {{0,0},{0,1},{0,2},{0,3},{0,0},{0,1},{0,2},{0,3},{0,0},{0,1},{0,2},{0,3}};
  }
  SECTION("column vector"){
    b={3,1};
    b_expected = {{0,0},{0,0},{0,0},{0,0},{1,0},{1,0},{1,0},{1,0},{2,0},{2,0},{2,0},{2,0}};
  }
  BroadcastIt<int> ab_itr(a, b);
  REQUIRE(ab_itr.shape().size() == a.size());
  REQUIRE(approx::vector(a.size(), a.data(), ab_itr.shape().data()));
  size_t n{0};
  for (auto [outer, ax, bx]: ab_itr){
    // std::cout << "{ ";
    // for(auto i: outer) std::cout << i << " ";
    // std::cout << "} <- { ";
    // for(auto i: ax) std::cout << i << " ";
    // std::cout << "} & { ";
    // for(auto i: bx) std::cout << i << " ";
    // std::cout << "}" << std::endl;
    //
    REQUIRE(outer.size() == ax.size());
    REQUIRE(outer.size() == bx.size());
    REQUIRE(approx::vector(outer.size(), outer.data(), ax.data())); // specific to this case;
    //
    REQUIRE(n < o_expected.size());
    REQUIRE(outer.size() == o_expected[n].size());
    REQUIRE(approx::vector(outer.size(), outer.data(), o_expected[n].data()));
    //
    REQUIRE(bx.size() == b_expected[n].size());
    REQUIRE(approx::vector(bx.size(), bx.data(), b_expected[n].data()));
    //
    ++n;
  }
  REQUIRE(o_expected.size() == n);
}

TEST_CASE("Equal-shape 'broadcasting'","[subscript]"){
  using namespace brille;
  std::vector<unsigned> a{3,3}, b{3,3};
  std::vector<std::vector<unsigned>>
  expected{
    {0,0},{0,1},{0,2}, {1,0},{1,1},{1,2}, {2,0},{2,1},{2,2}
  };
  BroadcastIt<unsigned> ab_itr(a, b);
  REQUIRE(ab_itr.shape().size() == a.size());
  REQUIRE(approx::vector(a.size(), a.data(), ab_itr.shape().data()));
  size_t n{0};
  for (auto [ox, ax, bx]: ab_itr){
    // std::cout << "{ ";
    // for(auto i: ox) std::cout << i << " ";
    // std::cout << "} <- { ";
    // for(auto i: ax) std::cout << i << " ";
    // std::cout << "} & { ";
    // for(auto i: bx) std::cout << i << " ";
    // std::cout << "}" << std::endl;
    //
    REQUIRE(ox.size() == ax.size());
    REQUIRE(ox.size() == bx.size());
    //
    REQUIRE(n < expected.size());
    REQUIRE(ox.size() == expected[n].size());
    REQUIRE(approx::vector(ox.size(), ox.data(), expected[n].data()));
    //
    REQUIRE(ax.size() == expected[n].size());
    REQUIRE(approx::vector(ax.size(), ax.data(), expected[n].data()));
    //
    REQUIRE(bx.size() == expected[n].size());
    REQUIRE(approx::vector(bx.size(), bx.data(), expected[n].data()));
    //
    ++n;
  }
  REQUIRE(expected.size() == n);
}

TEST_CASE("Higher dimension broadcasting","[subscript]"){
  using namespace brille;
  std::vector<int> o{3,2,4}, a{1,2,4}, b{3,2,1};
  std::vector<std::vector<int>>
  o_expected{
    {0,0,0},{0,0,1},{0,0,2},{0,0,3}, {0,1,0},{0,1,1},{0,1,2},{0,1,3},
    {1,0,0},{1,0,1},{1,0,2},{1,0,3}, {1,1,0},{1,1,1},{1,1,2},{1,1,3},
    {2,0,0},{2,0,1},{2,0,2},{2,0,3}, {2,1,0},{2,1,1},{2,1,2},{2,1,3},
  },
  a_expected{
    {0,0,0},{0,0,1},{0,0,2},{0,0,3}, {0,1,0},{0,1,1},{0,1,2},{0,1,3},
    {0,0,0},{0,0,1},{0,0,2},{0,0,3}, {0,1,0},{0,1,1},{0,1,2},{0,1,3},
    {0,0,0},{0,0,1},{0,0,2},{0,0,3}, {0,1,0},{0,1,1},{0,1,2},{0,1,3},
  },
  b_expected{
    {0,0,0},{0,0,0},{0,0,0},{0,0,0}, {0,1,0},{0,1,0},{0,1,0},{0,1,0},
    {1,0,0},{1,0,0},{1,0,0},{1,0,0}, {1,1,0},{1,1,0},{1,1,0},{1,1,0},
    {2,0,0},{2,0,0},{2,0,0},{2,0,0}, {2,1,0},{2,1,0},{2,1,0},{2,1,0},
  };
  BroadcastIt<int> ab_itr(a, b);
  REQUIRE(ab_itr.shape().size() == o.size());
  REQUIRE(approx::vector(o.size(), o.data(), ab_itr.shape().data()));
  size_t n{0};
  for (auto [ox, ax, bx]: ab_itr){
    // std::cout << "{ ";
    // for(auto i: ox) std::cout << i << " ";
    // std::cout << "} <- { ";
    // for(auto i: ax) std::cout << i << " ";
    // std::cout << "} & {";
    // for(auto i: bx) std::cout << i << " ";
    // std::cout << "}" << std::endl;
    //
    REQUIRE(ox.size() == ax.size());
    REQUIRE(ox.size() == bx.size());
    //
    REQUIRE(n < o_expected.size());
    REQUIRE(ox.size() == o_expected[n].size());
    REQUIRE(approx::vector(ox.size(), ox.data(), o_expected[n].data()));
    //
    REQUIRE(ax.size() == a_expected[n].size());
    REQUIRE(approx::vector(ax.size(), ax.data(), a_expected[n].data()));
    //
    REQUIRE(bx.size() == b_expected[n].size());
    REQUIRE(approx::vector(bx.size(), bx.data(), b_expected[n].data()));
    //
    ++n;
  }
  REQUIRE(o_expected.size() == n);
}


TEST_CASE("Full Array2 iterator","[subscript2]"){
  std::vector<std::vector<int>> expected{
    {0,0},{0,1},{0,2},{0,3},
    {1,0},{1,1},{1,2},{1,3},
    {2,0},{2,1},{2,2},{2,3},
  };
  SubIt2<int> sit({3,4});
  size_t n{0};
  for (auto x : sit){
    REQUIRE(n < expected.size());
    REQUIRE(x.size() == expected[n].size());
    REQUIRE(brille::approx::vector(x.size(), x.data(), expected[n++].data()));
  }
  REQUIRE( expected.size() == n);
}

TEST_CASE("Fixed-index Array2 iterator", "[subscript2]"){
  std::vector<std::vector<int>> expected{ {4,0},{4,1},{4,2} };
  SubIt<int> sit({5,3},{4,3});
  size_t n{0};
  for (auto x : sit){
    REQUIRE(n < expected.size());
    REQUIRE(x.size() == expected[n].size());
    REQUIRE(brille::approx::vector(x.size(), x.data(), expected[n++].data()));
  }
  REQUIRE( expected.size() == n);
}

TEST_CASE("Auto Array2 broadcasting","[subscript2]"){
  using namespace brille;
  std::array<int,2> a({3,4}), b;
  std::vector<std::vector<int>>
    o_expected{{0,0},{0,1},{0,2},{0,3},{1,0},{1,1},{1,2},{1,3},{2,0},{2,1},{2,2},{2,3}},
    b_expected;
  SECTION("row vector"){
    b={{1,4}};
    b_expected = {{0,0},{0,1},{0,2},{0,3},{0,0},{0,1},{0,2},{0,3},{0,0},{0,1},{0,2},{0,3}};
  }
  SECTION("column vector"){
    b={{3,1}};
    b_expected = {{0,0},{0,0},{0,0},{0,0},{1,0},{1,0},{1,0},{1,0},{2,0},{2,0},{2,0},{2,0}};
  }
  BroadcastIt2<int> ab_itr(a, b);
  REQUIRE(ab_itr.shape().size() == a.size());
  REQUIRE(approx::vector(a.size(), a.data(), ab_itr.shape().data()));
  size_t n{0};
  for (auto [outer, ax, bx]: ab_itr){
    // std::cout << "{ ";
    // for(auto i: outer) std::cout << i << " ";
    // std::cout << "} <- { ";
    // for(auto i: ax) std::cout << i << " ";
    // std::cout << "} & { ";
    // for(auto i: bx) std::cout << i << " ";
    // std::cout << "}" << std::endl;
    //
    REQUIRE(outer.size() == ax.size());
    REQUIRE(outer.size() == bx.size());
    REQUIRE(approx::vector(outer.size(), outer.data(), ax.data())); // specific to this case;
    //
    REQUIRE(n < o_expected.size());
    REQUIRE(outer.size() == o_expected[n].size());
    REQUIRE(approx::vector(outer.size(), outer.data(), o_expected[n].data()));
    //
    REQUIRE(bx.size() == b_expected[n].size());
    REQUIRE(approx::vector(bx.size(), bx.data(), b_expected[n].data()));
    //
    ++n;
  }
  REQUIRE(o_expected.size() == n);
}

TEST_CASE("Equal-shape Array2 'broadcasting'","[subscript2]"){
  using namespace brille;
  std::array<unsigned,2> a({3,3}), b({3,3});
  std::vector<std::vector<unsigned>>
  expected{
    {0,0},{0,1},{0,2}, {1,0},{1,1},{1,2}, {2,0},{2,1},{2,2}
  };
  BroadcastIt2<unsigned> ab_itr(a, b);
  REQUIRE(ab_itr.shape().size() == a.size());
  REQUIRE(approx::vector(a.size(), a.data(), ab_itr.shape().data()));
  size_t n{0};
  for (auto [ox, ax, bx]: ab_itr){
    // std::cout << "{ ";
    // for(auto i: ox) std::cout << i << " ";
    // std::cout << "} <- { ";
    // for(auto i: ax) std::cout << i << " ";
    // std::cout << "} & { ";
    // for(auto i: bx) std::cout << i << " ";
    // std::cout << "}" << std::endl;
    //
    REQUIRE(ox.size() == ax.size());
    REQUIRE(ox.size() == bx.size());
    //
    REQUIRE(n < expected.size());
    REQUIRE(ox.size() == expected[n].size());
    REQUIRE(approx::vector(ox.size(), ox.data(), expected[n].data()));
    //
    REQUIRE(ax.size() == expected[n].size());
    REQUIRE(approx::vector(ax.size(), ax.data(), expected[n].data()));
    //
    REQUIRE(bx.size() == expected[n].size());
    REQUIRE(approx::vector(bx.size(), bx.data(), expected[n].data()));
    //
    ++n;
  }
  REQUIRE(expected.size() == n);
}
