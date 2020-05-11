/* Copyright 2019 Greg Tucker
//
// This file is part of brille.
//
// brille is free software: you can redistribute it and/or modify it under the
// terms of the GNU Affero General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// brille is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with brille. If not, see <https://www.gnu.org/licenses/>.            */

#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <cmath>
#include <type_traits>
#include <limits>
#include <stdexcept>  // for std::overflow_error
#include <string>
#include <math.h>
// #include <complex> // included by debug
#include <numeric>
#include "debug.hpp"

// A mulitplier for the approximate-comparison tolerance
// 10000 is too big for Monoclinic system (5.7224, 5.70957, 4.13651),(90,90.498,90),'C -2y'
// but setting it any lower (9000 tried) causes other test lattices, namely,
// (7.189, 4.407, 5.069) (90,90.04,90) '-C 2y' to throw a runtime error
const int TOL_MULT=10000;

// We include numeric for std::gcd in C++17 (are we using it for anything else?)
/* Pre-C++17 gcd was not part of the standard but may be present in an
// experimental header -- which of course doesn't exist for MSVC
*/
#if defined(__cplusplus) && __cplusplus > 201700L
  #define STD_GCD
#endif

const double PI = std::atan(1.0)*4;
const double PICUBED = PI*PI*PI;
const double PIOVERTWO = PI/2.0;

template<typename T>
using IsInteger = std::enable_if_t<std::is_integral<T>::value,T>;
template<typename T>
using IsNotInteger = std::enable_if_t<!std::is_integral<T>::value,T>;

//! trace of a square matrix
template<typename T, int N=3> T trace(const T *M);

//! copying array from source (S) to destination (D)
template<typename T, int N, int M> void copy_array(T *D, const T *S);
template<typename T, int N=3> void copy_matrix(T *M, const T *A);
template<typename T, int N=3> void copy_vector(T *V, const T *A);

//! element-wise checking of approximate equality
template<typename T, int N, int M> bool equal_array(const T *A, const T *B, const T tol=0);
template<typename T, int N=3> bool equal_matrix(const T *A, const T *B, const T tol=0);
template<typename T, int N=3> bool equal_vector(const T *A, const T *B, const T tol=0);

/*
    This terrible class-based function specialization is all to get MSVC to stop complaining about a constant conditional expression.
    The four void classes tacked onto the end of the template are all to allow for specializations where either or both of T and R
    are floating point or integer numbers. If enable_if_t provides a type which is notvoid, then the specializations are:

       | T | R | TFloat  |  TInt   | RFloat  |  Rint   |
       |---|---|---------|---------|---------|---------|
       | 1.| 1.| notvoid |  void   | notvoid |  void   |
       | 1 | 1.|  void   | notvoid | notvoid |  void   |
       | 1.| 1 | notvoid |  void   |  void   | notvoid |
       | 1 | 1 |  void   | notvoid |  void   | notvoid |

*/
/*isfpT | isfpR | which? | why?
  ------|-------|--------|-----
    0       0     either    both Ttol and Rtol are 0
    1       1      Ttol     R is convertible to T
    0       1      Rtol     Ttol is 0, so use Rtol
    1       0      Ttol     Rtol is 0, so use Ttol
*/
/*template<class T, class R, class TFloat = void, class TInt = void, class RFloat = void, class RInt = void>
class determine_tols_impl {
public:
    static std::tuple<bool, bool, T, R> determine_tols(const int i);
};
template<class T, class R>
class determine_tols_impl<T,R,std::enable_if_t<std::is_floating_point<T>::value>, void, std::enable_if_t<std::is_floating_point<R>::value>, void> {
public:
    static std::tuple<bool, bool, T, R> determine_tols(const int tol) {
        T Ttol = static_cast<T>(tol)* static_cast<T>(TOL_MULT)* std::numeric_limits<T>::epsilon();
        R Rtol = static_cast<R>(tol)* static_cast<R>(TOL_MULT)* std::numeric_limits<R>::epsilon();
        return std::make_tuple(std::is_convertible<T, R>::value, true, Ttol, Rtol);
    }
};
template<class T, class R>
class determine_tols_impl <T,R,std::enable_if_t<std::is_floating_point<T>::value>, void, void, std::enable_if_t<std::is_integral<R>::value>> {
public:
    static std::tuple<bool, bool, T, R> determine_tols(const int tol) {
        T Ttol = static_cast<T>(tol)* static_cast<T>(TOL_MULT)* std::numeric_limits<T>::epsilon();
        return std::make_tuple(true, true, Ttol, R(0));
    }
};
template<class T, class R>
class determine_tols_impl<T,R, void, std::enable_if_t<std::is_integral<T>::value>, std::enable_if_t<std::is_floating_point<R>::value>, void> {
public:
    static std::tuple<bool, bool, T, R> determine_tols(const int tol) {
        R Rtol = static_cast<R>(tol)* static_cast<R>(TOL_MULT)* std::numeric_limits<R>::epsilon();
        return std::make_tuple(true, false, T(0), Rtol); // R can't generally be converted to integer T, but this was the behaviour before
    }
};
template<class T, class R>
class determine_tols_impl<T,R, void, std::enable_if_t<std::is_integral<T>::value>, void, std::enable_if_t<std::is_integral<R>::value>> {
public:
    static std::tuple<bool, bool, T, R> determine_tols(const int) {
        return std::make_tuple(true, false, T(0), R(0));
    }
};

template<class T, class R> std::tuple<bool, bool, T, R> determine_tols(const int i=1) {
    return determine_tols_impl<T, R>::determine_tols(i);
}
*/

template<typename T, typename R> std::tuple<bool,bool,T,R> determine_tols(const int=1);
template<typename T, typename R> bool approx_scalar(const T a, const R b, const int tol=1);
template<typename T, typename R> bool approx_array(const size_t N, const size_t M,const T *A, const R *B, const int tol=1);
template<typename T, typename R> bool approx_matrix(const size_t N, const T *A, const R *B, const int tol=1);
template<typename T, typename R> bool approx_vector(const size_t N, const T *A, const R *B, const int tol=1);


template<typename T, typename R, size_t N, size_t M> bool approx_array(const T *A, const R *B, const int tol=1);
template<typename T, typename R, size_t N=3> bool approx_matrix(const T *A, const R *B, const int tol=1);
template<typename T, typename R, size_t N=3> bool approx_vector(const T *A, const R *B, const int tol=1);


//! array multiplication C = A * B -- where C is (N,M), A is (N,I) and B is (I,M)
template<typename T, typename R, typename S, size_t N, size_t I, size_t M> void multiply_arrays(T *C, const R *A, const S *B);
template<typename T, typename R, typename S, size_t N=3> void multiply_matrix_matrix(T *C, const R *A, const S *B);
template<typename T, typename R, typename S, size_t N=3> void multiply_matrix_vector(T *C, const R *A, const S *b);
template<typename T, typename R, typename S, size_t N=3> void multiply_vector_matrix(T *C, const R *a, const S *B);

// array multiplication with specializations for complex array type(s)
template<class T, class R, class S> void mul_arrays(T* C, const size_t n, const size_t l, const size_t m, const R* A, const S* B);
template<class T, class R, class S> void mul_arrays(std::complex<T>* C, const size_t n, const size_t l, const size_t m, const R* A, const std::complex<S>* B);
template<class T, class R, class S> void mul_arrays(std::complex<T>* C, const size_t n, const size_t l, const size_t m, const std::complex<R>* A, const S* B);
template<class T, class R, class S> void mul_arrays(std::complex<T>* C, const size_t n, const size_t l, const size_t m, const std::complex<R>* A, const std::complex<S>* B);
// specializations for matrix*matrix, matrix*vector, and vector*matrix
template<class T, class R, class S> void mul_mat_mat(T* C, const size_t n, const R* A, const S* B){ mul_arrays(C, n, n, n, A, B);}
template<class T, class R, class S> void mul_mat_vec(T* C, const size_t n, const R* A, const S* B){ mul_arrays(C, n, n, 1, A, B);}
template<class T, class R, class S> void mul_vec_mat(T* C, const size_t n, const R* A, const S* B){ mul_arrays(C, 1, n, n, A, B);}

//! array element-wise addition
template<typename T, typename R, typename S, int N, int M> void add_arrays(T *C, const R *A, const S *B);
template<typename T, typename R, typename S, int N=3> void add_matrix(T *C, const R *A, const S *B);
//! array element-wise subtraction
template<typename T, typename R, typename S, int N, int M> void subtract_arrays(T *C, const R *A, const S *B);
template<typename T, typename R, typename S, int N> void subtract_matrix(T *C, const R *A, const S *B);

/*! special casting of floating point values to integer

 effectively rounds doubles to their nearest integer **not** by truncating.
 So 0.6=>1 not 0.
*/
template<typename T,typename R> T my_cast(const R a);
//! element-wise re-(special)-casting of arrays
template<typename T, typename R, int N, int M> void cast_array(T *A, const R *B);
template<typename T, typename R, int N> void cast_matrix(T *A, const R *B);
template<typename T, typename R, int N> void cast_vector(T *A, const R *B);

//! The cofactor array C_ij of the N by M array A is a M-1 by N-1 array of transpose(A with row i and column j) missing:
template<typename R> void array_cofactor(R *C, const R *A, const int i, const int j, const int N=3, const int M=3);
template<typename R> void matrix_cofactor(R *C, const R *A, const int i, const int j, const int N=3);

//! the determinant is only defined for square matrices
template<typename R> R matrix_determinant(const R *M, const int N=3);
//! the adjoint is also only defined for a square matrix
template<typename R> void matrix_adjoint(R *A, const R *M, const int N=3);
//! the inverse is, yet again, only defined for square matrices (with non-zero determinants)
template<typename R> R matrix_determinant_and_inverse(R *invM, const R *M, const R tol=0, const int N=3);
template<typename R> bool matrix_inverse(R *invM, const R *M, R tol=0, const int N=3);

template<typename T, int N=3> bool similar_matrix(T *M, const T *A, const T *B, const T tol=0);

template<typename R, int N, int M> void array_transpose(R *D, const R *S);
template<typename R, int N=3> void matrix_transpose(R *D, const R *S);
template<typename R, int N=3> void matrix_transpose(R *B);

template<typename R, int N=3> void matrix_metric(R *M, const R *L);

template<typename R, int N=3> R vector_norm_squared(const R * v);

template<typename T, typename R, typename S, int N> class vector_cross_impl {
public:
    static void vector_cross(T*, const R*, const S*) {
        throw std::runtime_error("The cross product is only defined for 3-vectors");
    }
};
template<typename T, typename R, typename S> class vector_cross_impl<T,R,S,3>{
public:
    static void vector_cross(T * c, const R * a, const S * b) {
        c[0] = static_cast<T>(a[1])* static_cast<T>(b[2]) - static_cast<T>(a[2])* static_cast<T>(b[1]);
        c[1] = static_cast<T>(a[2])* static_cast<T>(b[0]) - static_cast<T>(a[0])* static_cast<T>(b[2]);
        c[2] = static_cast<T>(a[0])* static_cast<T>(b[1]) - static_cast<T>(a[1])* static_cast<T>(b[0]);
    }
};
template<typename T, typename R, typename S, int N = 3>
void vector_cross(T* c, const R* a, const S* b) {
    vector_cross_impl<T, R, S, N>::vector_cross(c, a, b);
}
template<typename R, int N=3> R vector_dot(const R *a, const R *b);

template<typename T> T mod1(const T a);

template<typename T, typename R, int N=3> bool is_int_matrix(const T * A, const R tol);
template<typename R> bool is_int_matrix(const int *, const R);

/*! \brief The "distance" between two matrices using the Frobenius norm

The [Frobenius norm](http://mathworld.wolfram.com/FrobeniusNorm.html)
is the matrix norm of an 𝑚×𝑛 matrix 𝑎 given by

    |A|ᶠ = √∑ᵢᵐ∑ⱼⁿ|𝑎ᵢⱼ|²

or, for a square matrix,

    |A|ᶠ = √tr(𝑎𝑎ᴴ)

where 𝑎ᴴ is the conjugate transpose of 𝑎.
This function calculates the Frobenius norm of the matrix A-B.
@param A A n×n square matrix
@param B A n×n square matrix
@returns |A|ᶠ
*/
template<typename T> T frobenius_distance(const size_t n, const T* A, const T* B);
template<typename T> T frobenius_distance(const size_t n, const std::complex<T>* A, const std::complex<T>* B);

/*! \brief The general n-dimensional angle between two real-valued vectors

Calculate and return the angle θ between two vectors given by

    cos(θ) = <A,B>/|A||B|

@param n The dimensionality of the vector space
@param A A pointer to the first vector
@param B A pointer to the second vector
@returns θ
*/
template<typename T> T vector_angle(const size_t n, const T* A, const T* B);
//! A convenience function calling euclidean_angle
template<typename T> T vector_angle(const size_t n, const std::complex<T>* A, const std::complex<T>* B);
/*! \brief The general n-dimensional Euclidean angle between two complex-valued vectors

Calculate and return the Euclidean angle between to vectors in a complex vector
space Vᶜ(≃Cₙ, n∈N, n≥2) given by

    cos(θᵣ) = <A,B>ᵣ/|A||B|

where the real inner product of two complex vectors is performed in the
real vector space Vʳ (≃R₂ₙ) isometric to Vᶜ.

@param n The dimensionality of the complex vector space Vᶜ
@param A A pointer to the first complex vector
@param B A pointer to the second complex vector
@returns θᵣ
*/
template<typename T> T euclidean_angle(const size_t n, const std::complex<T>* A, const std::complex<T>* B);
/*! \brief The general n-dimensional Hermitian angle between two complex-valued vectors

Calculate and return the Hermitian angle between to vectors given by

    cos(θₕ) = |<A,B>|/|A||B|

where the inner product of two complex-valued vectors is
[given by](https://www.mathphysicsbook.com/mathematics/abstract-algebra/generalizing-vectors/norms-of-vectors/)
<A,B> = ¼(|A+B|²-|A-B|²+𝑖|A-𝑖B|²-𝑖|A+𝑖B|²)

@param n The dimensionality of the vector space
@param A A pointer to the first complex vector
@param B A pointer to the second complex vector
@returns θₕ
*/
template<typename T> T hermitian_angle(const size_t n, const std::complex<T>* A, const std::complex<T>* B);
template<typename T> T hermitian_angle(const size_t n, const T* A, const T* B);

template<typename T> T hermitian_product(const size_t n, const T* a, const T* b);
template<typename T> std::complex<T> hermitian_product(const size_t n, const T* a, const std::complex<T>* b);
template<typename T> std::complex<T> hermitian_product(const size_t n, const std::complex<T>* a, const T* b);
template<typename T> std::complex<T> hermitian_product(const size_t n, const std::complex<T>* a, const std::complex<T>* b);

template<typename T> T vector_distance(const size_t n, const T* a, const T* b);
template<typename T> T vector_distance(const size_t n, const std::complex<T>* a, const std::complex<T>* b);

template<typename T> T vector_product(const size_t n, const T* a, const T* b);
template<typename T> T vector_product(const size_t n, const T* a, const std::complex<T>* b);
template<typename T> T vector_product(const size_t n, const std::complex<T>* a, const T* b);
template<typename T> T vector_product(const size_t n, const std::complex<T>* a, const std::complex<T>* b);

template<typename T> T inner_product(const size_t n, const T* a, const T* b);
template<typename T> T inner_product(const size_t n, const std::complex<T>* a, const std::complex<T>* b);

template<class T> T squared_distance(const T&A, const T& B);
template<class T> T squared_distance(const std::complex<T>& A, const std::complex<T>& B);
template<class T> T magnitude(const T a);
template<class T> T magnitude(const std::complex<T> a);

template<class T> size_t encode_array_signs(const size_t n, const T* a);
template<class T> size_t encode_array_signs(const size_t n, const std::complex<T>* a);

template<class T> int make_eigenvectors_equivalent(const size_t n, const T* v0, T* v1);
template<class T> int make_eigenvectors_equivalent(const size_t n, const std::complex<T>* v0, std::complex<T> v1);

template<class T, class R>
enable_if_t<std::is_unsigned<T>::value&&std::is_unsigned<R>::value, unsigned long long>
binomial_coefficient(const T n, const R k);

template<typename S,typename U> S unsigned_to_signed(const U u);
template<typename U,typename S> U signed_to_unsigned(const S s);

#include "utilities.tpp"


#endif
