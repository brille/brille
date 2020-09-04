#include <random>
#include <chrono>

#include <catch2/catch.hpp>

#include "array.hpp"
#include "latvec.hpp" // otherwise we don't pick-up the common operators

TEST_CASE("Array creation","[array]"){
  brille::Array<double> novalues(3,3);
  REQUIRE(novalues.ndim() == 2);
  REQUIRE(novalues.numel() == 9);
  REQUIRE(novalues.size(0) == 3);
  REQUIRE(novalues.size(1) == 3);
  for (auto x: brille::ArrayIt(novalues)) REQUIRE(x == 0);

  double tmp[9] = {1,2,3,4,5,6,7,8,9};
  brille::shape_t sp{3,3}, st{3,1};
  bool brille_owns_data = false;
  brille::Array<double> withvalue(tmp, 9, brille_owns_data, sp, st);
  REQUIRE(withvalue.ndim() == 2);
  REQUIRE(withvalue.numel() == 9);
  REQUIRE(withvalue.size(0) == 3);
  REQUIRE(withvalue.size(1) == 3);
  size_t idx{0};
  for (auto x: brille::ArrayIt(withvalue)) REQUIRE(x == tmp[idx++]);

  brille::Array<double> copy_initialized(withvalue);
  REQUIRE(copy_initialized.ndim() == 2);
  REQUIRE(copy_initialized.numel() == 9);
  REQUIRE(copy_initialized.size(0) == 3);
  REQUIRE(copy_initialized.size(1) == 3);
  idx=0;
  for(auto x: brille::ArrayIt(copy_initialized)) REQUIRE(x == tmp[idx++]);

  brille::Array<double> to_be_assigned;
  REQUIRE(to_be_assigned.numel()==0);
  to_be_assigned = withvalue;
  REQUIRE(to_be_assigned.ndim() == 2);
  REQUIRE(to_be_assigned.numel() == 9);
  REQUIRE(to_be_assigned.size(0) == 3);
  REQUIRE(to_be_assigned.size(1) == 3);
  idx=0;
  for(auto x: brille::ArrayIt(to_be_assigned)) REQUIRE(x == tmp[idx++]);

  for (size_t i = 0; i<to_be_assigned.size(0); ++i){
    auto just_one = withvalue.view(i);
    REQUIRE(just_one.ndim() == 2);
    REQUIRE(just_one.numel() == 3);
    REQUIRE(just_one.size(0) == 1);
    REQUIRE(just_one.size(1) == 3);
    idx = 3*i;
    for (auto x: brille::ArrayIt(just_one)) REQUIRE(x == tmp[idx++]);
  }
}

TEST_CASE("Array scalar math operations","[array]"){
  double values[9] = {1,2,3,4,5,6,7,8,9};
  brille::Array<double> x(values, 9, /*owns*/ false, {3,3}, {3,1});
  auto orig = brille::Array<double>(x).decouple(); // make a distinct copy

  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(1.0,100.0);
  double rval = distribution(generator);
  auto itr = SubIt(x.shape());
  auto ait = brille::ArrayIt(x);
  size_t idx{0};
  SECTION("In-place adding a scalar effects all elements"){
    x += rval;
    for (auto i: itr) REQUIRE(x[i] == orig[i]+rval);
    for (auto v: ait) REQUIRE(v == values[idx++]);
  }
  SECTION("In-place subtracting a scalar effects all elements"){
    x -= rval;
    for (auto i: itr) REQUIRE(x[i] == orig[i]-rval);
    for (auto v: ait) REQUIRE(v == values[idx++]);
  }
  SECTION("In-place multiplying a scalar effects all elements"){
    x *= rval;
    for (auto i: itr) REQUIRE(x[i] == orig[i]*rval);
    for (auto v: ait) REQUIRE(v == values[idx++]);
  }
  SECTION("In-place dividing a scalar effects all elements"){
    x /= rval;
    for (auto i: itr) REQUIRE(x[i] == orig[i]/rval);
    for (auto v: ait) REQUIRE(v == values[idx++]);
  }
}

TEST_CASE("Array Array math operations","[array]"){
  double values[9] = {1,2,3,4,5,6,7,8,9};
  brille::Array<double> x(values, 9, /*owns*/ false, {3,3}, {3,1});
  auto orig = brille::Array<double>(x).decouple(); // make a distinct copy

  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(1.0,100.0);

  brille::shape_t sp{3,3};
  brille::Array<double> y(sp);
  // fill y with random values
  for (auto& v: brille::ArrayIt(y)) v = distribution(generator);

  auto itr = SubIt(x.shape());
  SECTION("In-place adding a second Array is element-wise addition"){
    x += y;
    for (auto i: itr) REQUIRE(x[i] == orig[i] + y[i]);
  }
  SECTION("In-place subtracting a second Array is element-wise subtraction"){
    x -= y;
    for (auto i: itr) REQUIRE(x[i] == orig[i] - y[i]);
  }
  SECTION("In-place multiplying a second Array is element-wise multiplication"){
    x *= y;
    for (auto i: itr) REQUIRE(x[i] == orig[i] * y[i]);
  }
  SECTION("In-place dividing a second Array is element-wise division"){
    x /= y;
    for (auto i: itr) REQUIRE(x[i] == orig[i] / y[i]);
  }
  SECTION("Array + Array is element-wise and scales-up single vectors"){
    brille::Array<double> z = x+y;
    for (auto i: itr) REQUIRE(z[i] == x[i]+y[i]);
    for (size_t j=0; j<y.size(0); ++j){
      auto yj = y.view(j);
      REQUIRE(yj.ndim() == 2);
      REQUIRE(yj.size(0) == 1);
      REQUIRE(yj.size(1) == 3);
      z = x+yj;
      REQUIRE(z.ndim() == 2);
      REQUIRE(z.size(0) == 3);
      REQUIRE(z.size(1) == 3);
      for (auto [iz,ix,iy]: BroadcastIt(x.shape(), yj.shape()))
        REQUIRE(z[iz] == x[ix] + yj[iy]);
    }
  }
  SECTION("Array - Array is element-wise and scales-up single vectors"){
    brille::Array<double> z = x-y;
    for (auto i: itr) REQUIRE(z[i] == x[i]-y[i]);
    for (size_t j=0; j<y.size(0); ++j){
      auto yj = y.view(j);
      REQUIRE(yj.ndim() == 2);
      REQUIRE(yj.size(0) == 1);
      REQUIRE(yj.size(1) == 3);
      z = x-yj;
      REQUIRE(z.ndim() == 2);
      REQUIRE(z.size(0) == 3);
      REQUIRE(z.size(1) == 3);
      for (auto [iz,ix,iy]: BroadcastIt(x.shape(), yj.shape()))
        REQUIRE(z[iz] == x[ix] - yj[iy]);
    }
  }
  SECTION("Array * Array is element-wise and scales-up single vectors"){
    brille::Array<double> z = x*y;
    for (auto i: itr) REQUIRE(z[i] == x[i]*y[i]);
    for (size_t j=0; j<y.size(0); ++j){
      auto yj = y.view(j);
      REQUIRE(yj.ndim() == 2);
      REQUIRE(yj.size(0) == 1);
      REQUIRE(yj.size(1) == 3);
      z = x*yj;
      REQUIRE(z.ndim() == 2);
      REQUIRE(z.size(0) == 3);
      REQUIRE(z.size(1) == 3);
      for (auto [iz,ix,iy]: BroadcastIt(x.shape(), yj.shape()))
        REQUIRE(z[iz] == x[ix] * yj[iy]);
    }
  }
  SECTION("Array / Array is element-wise and scales-up single vectors"){
    brille::Array<double> z = x/y;
    for (auto i: itr) REQUIRE(z[i] == x[i]/y[i]);
    for (size_t j=0; j<y.size(0); ++j){
      auto yj = y.view(j);
      REQUIRE(yj.ndim() == 2);
      REQUIRE(yj.size(0) == 1);
      REQUIRE(yj.size(1) == 3);
      z = x/yj;
      REQUIRE(z.ndim() == 2);
      REQUIRE(z.size(0) == 3);
      REQUIRE(z.size(1) == 3);
      for (auto [iz,ix,iy]: BroadcastIt(x.shape(), yj.shape()))
        REQUIRE(z[iz] == x[ix] / yj[iy]);
    }
  }
}

TEST_CASE("Append Array(s)","[array]"){
  double values[9] = {0,1,2,4,5,6,8,9,10};
  brille::Array<double> x(values, 9, /*owns*/ false, {3,3}, {3,1});
  double extra[3] = {3,7,11};
  brille::Array<double> y(extra, 3, /*owns*/ false, {3,1}, {1,0});
  auto z = cat(1,x,y);
  for (brille::ind_t i=0; i<12; ++i)
    REQUIRE(z[i] == Approx(static_cast<double>(i)));
}
