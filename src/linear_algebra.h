/*! \file */
#ifndef __linear_algebra_H__
#define __linear_algebra_H__

#ifndef SPGCONST
#define SPGCONST
#endif

#include <stdio.h>
#include <type_traits>
#include <limits>
#include <math.h>
#include <complex>

const double PI = std::atan(1.0)*4;
const double PICUBED = PI*PI*PI;
const double PIOVERTWO = PI/2.0;

//namespace rsm{

//! trace of a square matrix
template<typename T, int N=3> T trace(const T *M);

//! copying array from source (S) to destination (D)
template<typename T, int N, int M> void copy_array(T *D, const T *S);
template<typename T, int N=3> void copy_matrix(T *M, const T *A);
template<typename T, int N=3> void copy_vector(T *V, const T *A);

//! absolute value of a scalar
template<typename T> T my_abs(const T a);

//! element-wise checking of approximate equality
template<typename T, int N, int M> bool equal_array(const T *A, const T *B, const T tol=0);
template<typename T, int N=3> bool equal_matrix(const T *A, const T *B, const T tol=0);
template<typename T, int N=3> bool equal_vector(const T *A, const T *B, const T tol=0);

template<typename T, typename R> bool approx_scalar(const T a, const R b);
template<typename T, typename R, int N, int M> bool approx_array(const T *A, const R *B);
template<typename T, typename R, int N=3> bool approx_matrix(const T *A, const R *B);
template<typename T, typename R, int N=3> bool approx_vector(const T *A, const R *B);


//! array multiplication C = A * B -- where C is (N,M), A is (N,I) and B is (I,M)
template<typename T, typename R, typename S, int N, int I, int M> void multiply_arrays(T *C, const R *A, const S *B);
template<typename T, typename R, typename S, int N=3> void multiply_matrix_matrix(T *C, const R *A, const S *B);
template<typename T, typename R, typename S, int N=3> void multiply_matrix_vector(T *C, const R *A, const S *b);
template<typename T, typename R, typename S, int N=3> void multiply_vector_matrix(T *C, const R *a, const S *B);

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
template<typename T, typename R, typename S, int N=3> void vector_cross(T *c, const R *a, const S *b);
template<typename R, int N=3> R vector_dot(const R *a, const R *b);

template<typename T> T mod1(const T a);

template<typename T, typename R, int N=3> bool is_int_matrix(const T * A, const R tol);
template<typename R> bool is_int_matrix(const int *, const R);

template<typename T> T frobenius_distance(const T* A, const T* B, const size_t n);
template<typename T> T frobenius_distance(const std::complex<T>* A, const std::complex<T>* B, const size_t n);

#include "linear_algebra.hpp"

//} // namespace

#endif
