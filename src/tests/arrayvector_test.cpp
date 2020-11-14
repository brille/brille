#include <random>
#include <chrono>

#include <catch2/catch.hpp>
// for now we want to be able to switch between Array and Array2
// array_latvec.hpp defines a bArray<T> as one or the other:
#include "array_latvec.hpp"

using namespace brille;

TEST_CASE("Array creation","[array]"){
  bArray<double> novalues(3,3);
  REQUIRE(novalues.ndim() == 2);
  REQUIRE(novalues.numel() == 9);
  REQUIRE(novalues.size(0) == 3);
  REQUIRE(novalues.size(1) == 3);
  for (auto x: novalues.valItr()) REQUIRE(x == 0);

  double tmp[9] = {1,2,3,4,5,6,7,8,9};
  bArray<double>::shape_t sp({3,3}), st({3,1});
  bool brille_owns_data = false;
  bArray<double> withvalue(tmp, 9, brille_owns_data, sp, st);
  REQUIRE(withvalue.ndim() == 2);
  REQUIRE(withvalue.numel() == 9);
  REQUIRE(withvalue.size(0) == 3);
  REQUIRE(withvalue.size(1) == 3);
  size_t idx{0};
  for (auto x: withvalue.valItr()) REQUIRE(x == tmp[idx++]);

  bArray<double> copy_initialized(withvalue);
  REQUIRE(copy_initialized.ndim() == 2);
  REQUIRE(copy_initialized.numel() == 9);
  REQUIRE(copy_initialized.size(0) == 3);
  REQUIRE(copy_initialized.size(1) == 3);
  idx=0;
  for(auto x: copy_initialized.valItr()) REQUIRE(x == tmp[idx++]);

  bArray<double> to_be_assigned;
  REQUIRE(to_be_assigned.numel()==0);
  to_be_assigned = withvalue;
  REQUIRE(to_be_assigned.ndim() == 2);
  REQUIRE(to_be_assigned.numel() == 9);
  REQUIRE(to_be_assigned.size(0) == 3);
  REQUIRE(to_be_assigned.size(1) == 3);
  idx=0;
  for(auto x: to_be_assigned.valItr()) REQUIRE(x == tmp[idx++]);

  for (size_t i = 0; i<to_be_assigned.size(0); ++i){
    auto just_one = withvalue.view(i);
    REQUIRE(just_one.ndim() == 2);
    REQUIRE(just_one.numel() == 3);
    REQUIRE(just_one.size(0) == 1);
    REQUIRE(just_one.size(1) == 3);
    idx = 3*i;
    for (auto x: just_one.valItr()) REQUIRE(x == tmp[idx++]);
  }
}

TEST_CASE("Array scalar math operations","[array]"){
  double values[9] = {1,2,3,4,5,6,7,8,9};
  bArray<double> x(values, 9, /*owns*/ false, {3,3}, {3,1});
  auto orig = bArray<double>(x).decouple(); // make a distinct copy

  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(1.0,100.0);
  double rval = distribution(generator);
  size_t idx{0};
  SECTION("In-place adding a scalar effects all elements"){
    x += rval;
    for (auto i: x.subItr()) REQUIRE(x[i] == orig[i]+rval);
    for (auto v: x.valItr()) REQUIRE(v == values[idx++]);
  }
  SECTION("In-place subtracting a scalar effects all elements"){
    x -= rval;
    for (auto i: x.subItr()) REQUIRE(x[i] == orig[i]-rval);
    for (auto v: x.valItr()) REQUIRE(v == values[idx++]);
  }
  SECTION("In-place multiplying a scalar effects all elements"){
    x *= rval;
    for (auto i: x.subItr()) REQUIRE(x[i] == orig[i]*rval);
    for (auto v: x.valItr()) REQUIRE(v == values[idx++]);
  }
  SECTION("In-place dividing a scalar effects all elements"){
    x /= rval;
    for (auto i: x.subItr()) REQUIRE(x[i] == orig[i]/rval);
    for (auto v: x.valItr()) REQUIRE(v == values[idx++]);
  }
}

TEST_CASE("Array Array math operations","[array]"){
  double values[9] = {1,2,3,4,5,6,7,8,9};
  bArray<double> x(values, 9, /*owns*/ false, {3,3}, {3,1});
  auto orig = bArray<double>(x).decouple(); // make a distinct copy

  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(1.0,100.0);

  bArray<double>::shape_t sp({3,3});
  bArray<double> y(sp);
  // fill y with random values
  for (auto& v: y.valItr()) v = distribution(generator);

  SECTION("In-place adding a second Array is element-wise addition"){
    x += y;
    for (auto i: x.subItr()) REQUIRE(x[i] == orig[i] + y[i]);
  }
  SECTION("In-place subtracting a second Array is element-wise subtraction"){
    x -= y;
    for (auto i: x.subItr()) REQUIRE(x[i] == orig[i] - y[i]);
  }
  SECTION("In-place multiplying a second Array is element-wise multiplication"){
    x *= y;
    for (auto i: x.subItr()) REQUIRE(x[i] == orig[i] * y[i]);
  }
  SECTION("In-place dividing a second Array is element-wise division"){
    x /= y;
    for (auto i: x.subItr()) REQUIRE(x[i] == orig[i] / y[i]);
  }
  SECTION("Array + Array is element-wise and scales-up single vectors"){
    bArray<double> z = x+y;
    for (auto i: x.subItr()) REQUIRE(z[i] == x[i]+y[i]);
    for (size_t j=0; j<y.size(0); ++j){
      auto yj = y.view(j);
      REQUIRE(yj.ndim() == 2);
      REQUIRE(yj.size(0) == 1);
      REQUIRE(yj.size(1) == 3);
      z = x+yj;
      REQUIRE(z.ndim() == 2);
      REQUIRE(z.size(0) == 3);
      REQUIRE(z.size(1) == 3);
      for (auto [iz,ix,iy]: x.broadcastItr(yj))
        REQUIRE(z[iz] == x[ix] + yj[iy]);
    }
  }
  SECTION("Array - Array is element-wise and scales-up single vectors"){
    bArray<double> z = x-y;
    for (auto i: x.subItr()) REQUIRE(z[i] == x[i]-y[i]);
    for (size_t j=0; j<y.size(0); ++j){
      auto yj = y.view(j);
      REQUIRE(yj.ndim() == 2);
      REQUIRE(yj.size(0) == 1);
      REQUIRE(yj.size(1) == 3);
      z = x-yj;
      REQUIRE(z.ndim() == 2);
      REQUIRE(z.size(0) == 3);
      REQUIRE(z.size(1) == 3);
      for (auto [iz,ix,iy]: x.broadcastItr(yj))
        REQUIRE(z[iz] == x[ix] - yj[iy]);
    }
  }
  SECTION("Array * Array is element-wise and scales-up single vectors"){
    bArray<double> z = x*y;
    for (auto i: x.subItr()) REQUIRE(z[i] == x[i]*y[i]);
    for (size_t j=0; j<y.size(0); ++j){
      auto yj = y.view(j);
      REQUIRE(yj.ndim() == 2);
      REQUIRE(yj.size(0) == 1);
      REQUIRE(yj.size(1) == 3);
      z = x*yj;
      REQUIRE(z.ndim() == 2);
      REQUIRE(z.size(0) == 3);
      REQUIRE(z.size(1) == 3);
      for (auto [iz,ix,iy]: x.broadcastItr(yj))
        REQUIRE(z[iz] == x[ix] * yj[iy]);
    }
  }
  SECTION("Array / Array is element-wise and scales-up single vectors"){
    bArray<double> z = x/y;
    for (auto i: x.subItr()) REQUIRE(z[i] == x[i]/y[i]);
    for (size_t j=0; j<y.size(0); ++j){
      auto yj = y.view(j);
      REQUIRE(yj.ndim() == 2);
      REQUIRE(yj.size(0) == 1);
      REQUIRE(yj.size(1) == 3);
      z = x/yj;
      REQUIRE(z.ndim() == 2);
      REQUIRE(z.size(0) == 3);
      REQUIRE(z.size(1) == 3);
      for (auto [iz,ix,iy]: x.broadcastItr(yj))
        REQUIRE(z[iz] == x[ix] / yj[iy]);
    }
  }
}

TEST_CASE("Append Array(s)","[array]"){
  double values[9] = {0,1,2,4,5,6,8,9,10};
  bArray<double> x(values, 9, /*owns*/ false, {3,3}, {3,1});
  double extra[3] = {3,7,11};
  bArray<double> y(extra, 3, /*owns*/ false, {3,1}, {1,0});
  // std::cout << "x = " << std::endl << x.to_string() << "y = " << std::endl << y.to_string() << std::endl;
  auto z = cat(1,x,y);
  // std::cout << "cat(1,x,y)" << std::endl << z.to_string() << std::endl;
  bArray<double>::shape_t expected({3,4});
  REQUIRE(z.ndim() == 2);
  REQUIRE(z.size(0) == expected[0]);
  REQUIRE(z.size(1) == expected[1]);
  for (brille::ind_t i=0; i<12; ++i)
    REQUIRE(z[i] == Approx(static_cast<double>(i)));
}
