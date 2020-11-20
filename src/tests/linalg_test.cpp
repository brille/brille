#include <catch2/catch.hpp>
#include "utilities.hpp"
#include "approx.hpp"


using namespace brille;
using namespace brille::utils;

TEST_CASE("trace","[linalg]"){
  int mat[9] = {1,2,3, 4,5,6, 7,8,9};
  REQUIRE(trace<int,3>(mat) == 1+5+9);
}
TEST_CASE("array copying and equivalency","[linalg]"){
  int source[16] = {1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16};
  int dest[16] = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0};
  SECTION("copy_array"){
    copy_array<int,8,2>(dest,source);
    REQUIRE( equal_array<int,8,2>(dest,source) );
  }
  SECTION("copy_matrix"){
    copy_matrix<int,4>(dest,source);
    REQUIRE( equal_matrix<int,4>(dest,source) );
  }
  SECTION("copy_vector"){
    copy_vector<int,16>(dest,source);
    REQUIRE( equal_vector<int,16>(dest,source) );
  }
}
TEST_CASE("approx","[approx]"){
  // for some reason this compared float(1e-7) with double(1e-15) which has a difference
  // *of* 1e-7 and was only deemed approximately the same because a check was made whether
  // the result was less than 0.002384, which is not particularly useful.
  // Switching this 1e-7 to 1e-9 checks whether (1e-9 - 1e-15) < 3e-9, the absolute tolerance
  float  f[16] = {1e-9,1/2,-1/3,1/4,1/5,1/6,-1/7,1/8,1/9,-1/10,1/11,1/12,-1/13,1/14,1/15,1/16};
  double d[16] ={1e-15,1/2,-1/3,1/4,1/5,1/6,-1/7,1/8,1/9,-1/10,1/11,1/12,-1/13,1/14,1/15,1/16};
  SECTION("brille::approx::scalar"){ for (int i=0; i<16; i++) REQUIRE( brille::approx::scalar(f[i],d[i]) ); }
  SECTION("brille::approx::array"){  REQUIRE( brille::approx::array<float,double,2,8>(f,d) ); }
  SECTION("brille::approx::matrix"){ REQUIRE( brille::approx::matrix<float,double,4>(f,d) ); }
  SECTION("brille::approx::vector"){ REQUIRE( brille::approx::vector<float,double,16>(f,d) ); }
}

TEST_CASE("array multiplication","[linalg]"){
  int A[9] = {1,2,3,4,5,6,7,8,9};
  int B[9] = {10,-2,14,5,-3,8,-18,0,4};
  int W[12] = {1,2,3,4,5,6,7,8,9,10,11,12};
  int v[3] = {1,2,3};
  SECTION("multiply_arrays"){
    int R[12];
    SECTION("(3,3)*(3,4)"){
      int expected[12]={126,148,170,192, 62,72,82,92, 18,4,-10,-24};
      multiply_arrays<int,int,int,3,3,4>(R,B,W);
      REQUIRE( equal_array<int,3,4>(R,expected) );
    }
    SECTION("(4,3)*(3,3)"){
      int expected[12]={-34,-8,42, -43,-23,120, -52,-38,198, -61,-53,276};
      multiply_arrays<int,int,int,4,3,3>(R,W,B);
      REQUIRE( equal_array<int,4,3>(R,expected) );
    }
  }
  SECTION("multiply_matrix_matrix"){
    int C[9], expected[9]={-34,-8,42,-43,-23,120,-52,-38,198};
    multiply_matrix_matrix<int,int,int,3>(C,A,B);
    REQUIRE( equal_matrix<int,3>(C,expected) );
  }
  SECTION("multiply_matrix_vector"){
    int x[3], expected[3]={48,23,-6};
    multiply_matrix_vector<int,int,int,3>(x,B,v);
    REQUIRE( equal_vector<int,3>(x,expected) );
  }
  SECTION("multiply_vector_matrix"){
    int x[3], expected[3]={-34,-8,42};
    multiply_vector_matrix<int,int,int,3>(x,v,B);
    REQUIRE( equal_vector<int,3>(x,expected) );
  }
}

TEST_CASE("array addition","[linalg]"){
  int A[16]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
  double B[16]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.01,0.11,0.21,0.31,0.41,0.51,0.61};
  double C[16], expected[16]={1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9,10.01,11.11,12.21,13.31,14.41,15.51,16.61};
  SECTION("add_arrays"){
    add_arrays<double,int,double,2,8>(C,A,B);
    REQUIRE( equal_array<double,2,8>(C,expected) );
  }
  SECTION("add_matrix"){
    add_matrix<double,int,double,4>(C,A,B);
    REQUIRE( equal_matrix<double,4>(C,expected) );
  }
}

TEST_CASE("specialized (rounding for double->int) casting","[linalg]"){
  double cD[9], D[9] = {-1.9, -1.1, -0.3, 0.6, 2.2, 1.4, 0.3, -0.7, -2.1};
  int    cI[9], I[9] = {-2,   -1,    0,   1,   2,   1,   0,   -1,   -2  };
  SECTION("my_cast"){
    for(int i=0; i<9; i++){
      REQUIRE( my_cast<double>(D[i]) == D[i] );
      REQUIRE( my_cast<int>(D[i]) == I[i] );
    }
  }
  SECTION("cast_matrix"){
    cast_matrix<double,double,3>(cD,D);
    REQUIRE( equal_matrix<double,3>(cD,D) );
    cast_matrix<int,double,3>(cI,D);
    REQUIRE( equal_matrix(cI,I) );
  }
  SECTION("cast_vector"){
    cast_vector<double,double,9>(cD,D);
    REQUIRE( equal_vector<double,9>(cD,D) );
    cast_vector<int,double,9>(cI,D);
    REQUIRE( equal_vector<int,9>(cI,I) );
  }
}

TEST_CASE("determinant and inverse of matrices","[linalg]"){
  double M[9] = {1,0,1, 0,1,0, 0,1,1};
  double invM[9], expected_invM[9] = {1,1,-1, 0,1,0, 0,-1,1};
  REQUIRE( matrix_determinant<double>(M) == Approx(1) );
  REQUIRE( matrix_inverse<double>(invM, M) ); // require that the inverse exists
  REQUIRE( brille::approx::matrix<double,double,3>(invM, expected_invM) ); // and that the result is correct
}

TEST_CASE("similar matrix","[linalg]"){
  double A[9] = {1,2,3,4,5,6,7,8,9};
  double M[9] = {1,0,1, 0,1,0, 0,1,1};
  double S[9], expected[9] = {-2,-1,-2,4,11,10,3,6,6};
  REQUIRE( similar_matrix(S,A,M) );
  REQUIRE( brille::approx::matrix<double,double,3>(S,expected) );
}

TEST_CASE("array transpose","[linalg]"){
  int M[12] = {1,2,3,4, 5,6,7,8, 9,10,11,12};
  int Mt[12];
  SECTION("(3,4)->(4,3)"){
    int expected[12] = {1,5,9, 2,6,10, 3,7,11, 4,8,12};
    array_transpose<int,3,4>(Mt,M);
    REQUIRE( equal_array<int,4,3>(Mt,expected) );
  }
  SECTION("(4,3)->(3,4)"){
    int expected[12] = {1,4,7,10, 2,5,8,11, 3,6,9,12};
    array_transpose<int,4,3>(Mt,M);
    REQUIRE( equal_array<int,3,4>(Mt,expected) );
  }
}
TEST_CASE("matrix transpose","[linalg]"){
  int M[9] = {1,2,3,4,5,6,7,8,9};
  int expected[9] = {1,4,7, 2,5,8, 3,6,9};
  SECTION("store transpose in second matrix"){
    int Mt[9];
    matrix_transpose<int,3>(Mt,M);
    REQUIRE( equal_matrix<int,3>(Mt,expected) );
  }
  SECTION("transpose in-place"){
    matrix_transpose(M);
    REQUIRE( equal_matrix<int,3>(M,expected) );
  }
}

TEST_CASE("matrix metric","[linalg]"){
  double B[9] = {10,-2,14,5,-3,8,-18,0,4};
  double Bm[9], expected[9]={449,-35,108,-35,13,-52,108,-52,276};
  matrix_metric<double,3>(Bm,B);
  REQUIRE( brille::approx::matrix<double,double,3>(Bm,expected) );
}
TEST_CASE("vector norm squared","[linalg]"){
  double v[3] = {1,-2,3};
  REQUIRE( vector_norm_squared<double,3>(v) == Approx(14.0) );
}
TEST_CASE("vector cross products","[linalg]"){
  double v1[3] = {1,0,0}, v2[3] = {0,1,0}, v3[3] = {0,0,1};
  double cross[3];
  SECTION("(100)×(010)"){
    vector_cross<double,double,double,3>(cross,v1,v2);
    REQUIRE( brille::approx::vector<double,double,3>(cross,v3) );
  }
  SECTION("(010)×(001)"){
    vector_cross<double,double,double,3>(cross,v2,v3);
    REQUIRE( brille::approx::vector<double,double,3>(cross,v1) );
  }
  SECTION("(001)×(100)"){
    vector_cross<double,double,double,3>(cross,v3,v1);
    REQUIRE( brille::approx::vector<double,double,3>(cross,v2) );
  }
}
TEST_CASE("vector dot products","[linalg]"){
  double v1[3] = {1,1,0}, v2[3] = {0,1,1}, v3[3] = {1,0,1};
  // double dot[3];
  SECTION("(110)⋅(011)"){ REQUIRE( vector_dot<double,3>(v1,v2) == Approx(1.0) );}
  SECTION("(011)⋅(101)"){ REQUIRE( vector_dot<double,3>(v2,v3) == Approx(1.0) );}
  SECTION("(101)⋅(110)"){ REQUIRE( vector_dot<double,3>(v3,v1) == Approx(1.0) );}
}
