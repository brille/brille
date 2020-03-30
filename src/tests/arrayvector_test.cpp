#include <random>
#include <chrono>

#include <catch2/catch.hpp>

#include "arrayvector.hpp"

TEST_CASE("ArrayVector creation","[arrayvector]"){
  ArrayVector<double> novalues(3,3);
  REQUIRE( novalues.numel() == 3);
  REQUIRE( novalues.size() == 3);
  for (size_t i=0; i<novalues.size(); ++i)
    for (size_t j=0; j<novalues.numel(); ++j)
      REQUIRE( novalues.getvalue(i,j) == 0);

  double tmp[9] = {1,2,3,4,5,6,7,8,9};
  ArrayVector<double> withvalue(3,3,tmp);
  REQUIRE( withvalue.numel()==3);
  REQUIRE( withvalue.size() ==3);
  for (size_t i=0; i<withvalue.size(); ++i)
    for (size_t j=0; j<withvalue.numel(); ++j)
      REQUIRE( withvalue.getvalue(i,j) == tmp[i*withvalue.numel()+j]);

  ArrayVector<double> copy_initialized = withvalue;
  REQUIRE( copy_initialized.numel()==3);
  REQUIRE( copy_initialized.size() ==3);
  for (size_t i=0; i<copy_initialized.size(); ++i)
    for (size_t j=0; j<copy_initialized.numel(); ++j)
      REQUIRE( copy_initialized.getvalue(i,j) == tmp[i*copy_initialized.numel()+j]);

  ArrayVector<double> to_be_assigned;
  REQUIRE(to_be_assigned.numel()==0);
  REQUIRE(to_be_assigned.size()==0);
  to_be_assigned = withvalue;
  REQUIRE( to_be_assigned.numel()==3);
  REQUIRE( to_be_assigned.size() ==3);
  for (size_t i=0; i<to_be_assigned.size(); ++i)
    for (size_t j=0; j<to_be_assigned.numel(); ++j)
      REQUIRE( to_be_assigned.getvalue(i,j) == tmp[i*to_be_assigned.numel()+j]);

  ArrayVector<double> just_one;
  for (size_t i = 0; i<to_be_assigned.size(); ++i){
      just_one = withvalue.extract(i);
      REQUIRE(just_one.size() == 1);
      REQUIRE(just_one.numel() == withvalue.numel() );
      for (size_t j=0; j<just_one.numel(); ++j)
        REQUIRE( just_one.getvalue(0,j) == withvalue.getvalue(i,j) );
  }

}

TEST_CASE("ArrayVector scalar math operations","[arrayvector]"){
    double values[9] = {1,2,3,4,5,6,7,8,9};
    ArrayVector<double> x(3,3,values);

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> distribution(1.0,100.0);
    double rval = distribution(generator);
    SECTION("In-place adding a scalar effects all elements"){
        x += rval;
        for (size_t i=0; i<x.size(); ++i)
            for (size_t j=0; j<x.numel(); ++j)
                REQUIRE(x.getvalue(i,j) == values[i*x.numel()+j] + rval);
    }
    SECTION("In-place subtracting a scalar effects all elements"){
        x -= rval;
        for (size_t i=0; i<x.size(); ++i)
            for (size_t j=0; j<x.numel(); ++j)
                REQUIRE(x.getvalue(i,j) == values[i*x.numel()+j] - rval);
    }
    SECTION("In-place multiplying a scalar effects all elements"){
        x *= rval;
        for (size_t i=0; i<x.size(); ++i)
            for (size_t j=0; j<x.numel(); ++j)
                REQUIRE(x.getvalue(i,j) == values[i*x.numel()+j] * rval);
    }
    SECTION("In-place dividing a scalar effects all elements"){
        x /= rval;
        for (size_t i=0; i<x.size(); ++i)
            for (size_t j=0; j<x.numel(); ++j)
                REQUIRE(x.getvalue(i,j) == values[i*x.numel()+j] / rval);
    }
}

TEST_CASE("ArrayVector ArrayVector math operations","[arrayvector]"){
    double values[9] = {1,2,3,4,5,6,7,8,9};
    ArrayVector<double> x(3,3,values);

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> distribution(1.0,100.0);
    // double rval = distribution(generator);

    ArrayVector<double> y(3,3);
    for (size_t i=0; i<y.size(); ++i)
        for (size_t j=0; j<y.numel(); ++j)
            y.insert(distribution(generator),i,j); // set y to random values

    SECTION("In-place adding a second ArrayVector is element-wise addition"){
        x += y;
        for (size_t i=0; i<x.size(); ++i)
            for (size_t j=0; j<x.numel(); ++j)
                REQUIRE( x.getvalue(i,j) == values[i*x.numel()+j] + y.getvalue(i,j));
    }
    SECTION("In-place subtracting a second ArrayVector is element-wise subtraction"){
        x -= y;
        for (size_t i=0; i<x.size(); ++i)
            for (size_t j=0; j<x.numel(); ++j)
                REQUIRE( x.getvalue(i,j) == values[i*x.numel()+j] - y.getvalue(i,j));
    }
    SECTION("In-place multiplying a second ArrayVector is element-wise multiplication"){
        x *= y;
        for (size_t i=0; i<x.size(); ++i)
            for (size_t j=0; j<x.numel(); ++j)
                REQUIRE( x.getvalue(i,j) == values[i*x.numel()+j] * y.getvalue(i,j));
    }
    SECTION("In-place dividing a second ArrayVector is element-wise division"){
        x /= y;
        for (size_t i=0; i<x.size(); ++i)
            for (size_t j=0; j<x.numel(); ++j)
                REQUIRE( x.getvalue(i,j) == values[i*x.numel()+j] / y.getvalue(i,j));
    }
    SECTION("ArrayVector + ArrayVector is element-wise and scales-up single vectors"){
        ArrayVector<double> z = x+y;
        for (size_t i=0; i<x.size(); ++i)
            for (size_t j=0; j<x.numel(); ++j)
                REQUIRE( z.getvalue(i,j) == x.getvalue(i,j)+y.getvalue(i,j));
        for (size_t k=0; k<y.size(); ++k){
            z = x + y.extract(k);
            REQUIRE( z.size() == x.size() );
            REQUIRE( z.numel() == x.numel() );
            REQUIRE( y.size() == x.size() );
            REQUIRE( y.numel() == x.numel() );
            for (size_t i=0; i<x.size(); ++i)
                for (size_t j=0; j<x.numel(); ++j)
                    REQUIRE( z.getvalue(i,j) == x.getvalue(i,j)+y.getvalue(k,j));
        }
    }
    SECTION("ArrayVector - ArrayVector is element-wise and scales-up single vectors"){
        ArrayVector<double> z = x-y;
        for (size_t i=0; i<x.size(); ++i)
            for (size_t j=0; j<x.numel(); ++j)
                REQUIRE( z.getvalue(i,j) == x.getvalue(i,j)-y.getvalue(i,j));
        for (size_t k=0; k<y.size(); ++k){
            z = x - y.extract(k);
            REQUIRE( z.size() == x.size() );
            REQUIRE( z.numel() == x.numel() );
            REQUIRE( y.size() == x.size() );
            REQUIRE( y.numel() == x.numel() );
            for (size_t i=0; i<x.size(); ++i)
                for (size_t j=0; j<x.numel(); ++j)
                    REQUIRE( z.getvalue(i,j) == x.getvalue(i,j)-y.getvalue(k,j));
        }
    }
    SECTION("ArrayVector * ArrayVector is element-wise and scales-up single vectors"){
        ArrayVector<double> z = x*y;
        for (size_t i=0; i<x.size(); ++i)
            for (size_t j=0; j<x.numel(); ++j)
                REQUIRE( z.getvalue(i,j) == x.getvalue(i,j)*y.getvalue(i,j));
        for (size_t k=0; k<y.size(); ++k){
            z = x * y.extract(k);
            REQUIRE( z.size() == x.size() );
            REQUIRE( z.numel() == x.numel() );
            REQUIRE( y.size() == x.size() );
            REQUIRE( y.numel() == x.numel() );
            for (size_t i=0; i<x.size(); ++i)
                for (size_t j=0; j<x.numel(); ++j)
                    REQUIRE( z.getvalue(i,j) == x.getvalue(i,j)*y.getvalue(k,j));
        }
    }
    SECTION("ArrayVector / ArrayVector is element-wise and scales-up single vectors"){
        ArrayVector<double> z = x/y;
        for (size_t i=0; i<x.size(); ++i)
            for (size_t j=0; j<x.numel(); ++j)
                REQUIRE( z.getvalue(i,j) == x.getvalue(i,j)/y.getvalue(i,j));
        for (size_t k=0; k<y.size(); ++k){
            z = x / y.extract(k);
            REQUIRE( z.size() == x.size() );
            REQUIRE( z.numel() == x.numel() );
            REQUIRE( y.size() == x.size() );
            REQUIRE( y.numel() == x.numel() );
            for (size_t i=0; i<x.size(); ++i)
                for (size_t j=0; j<x.numel(); ++j)
                    REQUIRE( z.getvalue(i,j) == x.getvalue(i,j)/y.getvalue(k,j));
        }
    }

}
