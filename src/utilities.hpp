/* This file is part of brille.

Copyright © 2019,2020 Greg Tucker <greg.tucker@stfc.ac.uk>

brille is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

brille is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with brille. If not, see <https://www.gnu.org/licenses/>.            */

#ifndef BRILLE_UTILITIES_HPP_
#define BRILLE_UTILITIES_HPP_
/*! \file
    \author Greg Tucker
    \brief Various utility functions, mostly linear algebra
*/
// #include <cmath>
// #include <type_traits>
// #include <limits>
// #include <stdexcept>  // for std::overflow_error
// #include <string>
// #include <math.h>
// #include <numeric>
// #include <iostream>
// #include "debug.hpp"
#include "approx.hpp"


namespace brille{
  const double             pi = 3.14159265358979323846;
  const double         halfpi = 1.57079632679489661923;
  // const double      quarterpi = 0.785398163397448309616;
  // const double      inversepi = 0.318309886183790671538;
  // const double      twooverpi = 0.636619772367581343076;
  // const double  twooversqrtpi = 1.12837916709551257390;
  // const double        sqrttwo = 1.41421356237309504880;
  // const double sqrttwoovertwo = 0.707106781186547524401;
  // const double              e = 2.71828182845904523536;
  // const double          log2e = 1.44269504088896340736;
  // const double         log10e = 0.434294481903251827651;
  // const double            ln2 = 0.693147180559945309417;
  // const double           ln10 = 2.30258509299404568402;

  //! Utility functions for `brille`, mostly implementing linear algebra
  namespace utils{

//! trace of a square matrix
template<typename T, int N=3> T trace(const T *M);

//! copy an `N`×`M` array from source to destination
template<typename T, int N, int M> void copy_array(T *D /*!< destination */, const T *S /*!< source */);
//! copy a `N`×`N` matrix from source to destination
template<typename T, int N=3> void copy_matrix(T *M /*!< destination */, const T *A /*!< source */);
//! copy a `N` vector from source to destination
template<typename T, int N=3> void copy_vector(T *V /*!< desitnation */, const T *A /*!< source */);

//! Exact or approximate equality check for `N`×`M` arrays
template<typename T, int N, int M> bool equal_array(const T *A /*!< a*/, const T *B/*!< b*/, const T tol=0/*!< tolerance multiplier*/);
//! Exact or approximate equality check for `N`×`N` matrices
template<typename T, int N=3> bool equal_matrix(const T *A /*!< a*/, const T *B/*!< b*/, const T tol=0/*!< tolerance multiplier*/);
//! Exact or approximate equality check for `N` vectors
template<typename T, int N=3> bool equal_vector(const T *A /*!< a*/, const T *B/*!< b*/, const T tol=0/*!< tolerance multiplier*/);


/*! \brief Array multiplation

\f$ C = A B \f$

\param[out] C pointer to the first element of a `N`×`M` array
\param[in]  A pointer to the first element of a `N`×`I` array
\param[in]  B pointer to the first element of a `I`×`M` array
*/
template<typename T, typename R, typename S, size_t N, size_t I, size_t M> void multiply_arrays(T *C, const R *A, const S *B);
/*! \brief Matrix multiplation

\f$ C = A B \f$

\param[out] C pointer to the first element of a `N`×`N` array
\param[in]  A pointer to the first element of a `N`×`N` array
\param[in]  B pointer to the first element of a `N`×`M` array
*/
template<typename T, typename R, typename S, size_t N=3> void multiply_matrix_matrix(T *C, const R *A, const S *B);
/*! \brief Matrix-vector multiplation

\f$ C = A v \f$

\param[out] C pointer to the first element of a `N` vector
\param[in]  A pointer to the first element of a `N`×`N` array
\param[in]  v pointer to the first element of a `N` vector
*/
template<typename T, typename R, typename S, size_t N=3> void multiply_matrix_vector(T *C, const R *A, const S *v);
/*! \brief Vector-matrix multiplation

\f$ C = v^T B \f$

\param[out] C pointer to the first element of a `N` vector
\param[in]  v pointer to the first element of a `N` vector
\param[in]  B pointer to the first element of a `N`×`N` array
*/
template<typename T, typename R, typename S, size_t N=3> void multiply_vector_matrix(T *C, const R *v, const S *B);

/*! \brief Array multiplation with runtime specified dimensions

\f$ C = A B \f$

\param[out] C pointer to the first element of a `n`×`m` array
\param      n The first dimension size of the `C` and `A` arrays
\param      l The second dimension size of `A` and first dimension size of `B`
\param      m The second dimension size of the `C` and `B` arrays
\param[in]  A pointer to the first element of a `n`×`l` array
\param[in]  B pointer to the first element of a `l`×`m` array

\note This function is overloaded for complex-valued `A` or `B` (or both) in
      which case the output `C` must also be complex.
*/
template<class T, class R, class S, class I> void mul_arrays(T* C, const I n, const I l, const I m, const R* A, const S* B);
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T, class R, class S, class I> void mul_arrays(std::complex<T>* C, const I n, const I l, const I m, const R* A, const std::complex<S>* B);
template<class T, class R, class S, class I> void mul_arrays(std::complex<T>* C, const I n, const I l, const I m, const std::complex<R>* A, const S* B);
template<class T, class R, class S, class I> void mul_arrays(std::complex<T>* C, const I n, const I l, const I m, const std::complex<R>* A, const std::complex<S>* B);
#endif
/*! \brief Matrix multiplation with runtime specified dimension

\f$ C = A B \f$

\param[out] C pointer to the first element of a `n`×`n` matrix
\param      n The dimension of the system
\param[in]  A pointer to the first element of a `n`×`n` matrix
\param[in]  B pointer to the first element of a `n`×`n` matrix

\see mul_arrays
*/
template<class T, class R, class S, class I> void mul_mat_mat(T* C, const I n, const R* A, const S* B){ mul_arrays(C, n, n, n, A, B);}
/*! \brief Matrix-vector multiplation with runtime specified dimension

\f$ C = A v \f$

\param[out] C pointer to the first element of a `n` vector
\param      n The dimension of the system
\param[in]  A pointer to the first element of a `n`×`n` matrix
\param[in]  B pointer to the first element of a `n` vector

\see mul_arrays
*/
template<class T, class R, class S, class I> void mul_mat_vec(T* C, const I n, const R* A, const S* B){ mul_arrays(C, n, n, I(1), A, B);}
/*! \brief Vector-matrix multiplation with runtime specified dimension

\f$ C = v B \f$

\param[out] C pointer to the first element of a `n` vector
\param      n The dimension of the system
\param[in]  A pointer to the first element of a `n` vector
\param[in]  B pointer to the first element of a `n`×`n` matrix

\see mul_arrays
*/
template<class T, class R, class S, class I> void mul_vec_mat(T* C, const I n, const R* A, const S* B){ mul_arrays(C, I(1), n, n, A, B);}

// template<class T,class R,class I,class S=std::common_type_t<T,R>>
// std::enable_if_t<std::is_unsigned_v<I>&&std::is_same_v<S,R>>
// mul_mat_vec_inplace(const I n, const T* A, R* B){
//   S* C = new S[n]();
//   mul_arrays<S,T,R,I>(C, n, n, I(1), A, B);
//   for (I i=0; i<n; ++i) B[n] = C[n];
//   delete[] C;
// }

/*! \brief In-place Matrix-vector multiplication with runtime specified dimension

Performs the assignment \f$ v = A v \f$.

\param         n the dimension of the system
\param[in]     A pointer to the first element of a `n`×`n` matrix
\param[in,out] v pointer to the first element of a `n` vector

\note A temporary array is created and destroyed within this function
*/
template<class T,class R,class I,class S=std::common_type_t<T,R>>
std::enable_if_t<std::is_unsigned_v<I>>
mul_mat_vec_inplace(const I n, const T* A, R* v){
  S* c = new S[n]();
  mul_arrays(c, n, n, I(1), A, v);
  std::copy(c,c+n,v);
  delete[] c;
}

/*! \brief Addition of two arrays

\f$ C = A + B \f$

\param[out] C pointer to the first element of a `N`×`M` array
\param[in]  A pointer to the first element of a `N`×`M` array
\param[in]  B pointer to the first element of a `N`×`M` array
*/
template<typename T, typename R, typename S, int N, int M> void add_arrays(T *C, const R *A, const S *B);
/*! \brief Addition of two matrices

\f$ C = A + B \f$

\param[out] C pointer to the first element of a `N`×`N` matrix
\param[in]  A pointer to the first element of a `N`×`N` matrix
\param[in]  B pointer to the first element of a `N`×`N` matrix
*/
template<typename T, typename R, typename S, int N=3> void add_matrix(T *C, const R *A, const S *B);
/*! \brief Subtraction of two arrays

\f$ C = A - B \f$

\param[out] C pointer to the first element of a `N`×`M` array
\param[in]  A pointer to the first element of a `N`×`M` array
\param[in]  B pointer to the first element of a `N`×`M` array
*/
template<typename T, typename R, typename S, int N, int M> void subtract_arrays(T *C, const R *A, const S *B);
/*! \brief Subtraction of two matrices

\f$ C = A - B \f$

\param[out] C pointer to the first element of a `N`×`N` matrix
\param[in]  A pointer to the first element of a `N`×`N` matrix
\param[in]  B pointer to the first element of a `N`×`N` matrix
*/
template<typename T, typename R, typename S, int N> void subtract_matrix(T *C, const R *A, const S *B);

/*! \brief Special casting of floating point values to integer

 Effectively rounds doubles to their nearest integer **not** by truncating.
 So 0.6=>1 not 0.

 \param a the value to be cast to type `T`
 \returns The closest representation of `a` in type `T`.
*/
template<typename T,typename R> T my_cast(const R a);
/*! \brief Special casting of the elements of an array

\param[out] A pointer to the first element of a `N`×`M` array
\param[in]  B pointer to the first element of a `N`×`M` array
\see my_cast
*/
template<typename T, typename R, int N, int M> void cast_array(T *A, const R *B);
/*! \brief Special casting of the elements of a matrix

\param[out] A pointer to the first element of a `N`×`N` matrix
\param[in]  B pointer to the first element of a `N`×`N` matrix
\see my_cast
*/
template<typename T, typename R, int N> void cast_matrix(T *A, const R *B);
/*! \brief Special casting of the elements of an vector

\param[out] A pointer to the first element of a `N` vector
\param[in]  B pointer to the first element of a `N` vector
\see my_cast
*/
template<typename T, typename R, int N> void cast_vector(T *A, const R *B);

/*! \brief Find the cofactor of an array

The cofactor array \f$C_{ij}\f$ is a `M-1`×`N-1` array which is the transpose of
a `N-1`×`M-1` array \f$A_{ij}\f$, which is itself the elements of a `N`×`M`
array \f$A\f$ with row \f$i\f$ and column \f$j\f$ j missing.

\param[out] C pointer to the first element of a `M-1`×`N-1` array
\param[in]  A poitner to the first element of a `N`×`M` array
\param      i the row of `A` to omit, `i` ∈ (0,`N`]
\param      j the column of `A` to omit, `j` ∈ (0,`M`]
\param      N the first dimension of `A`
\param      M the second dimension of `A`
*/
template<typename R> void array_cofactor(R *C, const R *A, const int i, const int j, const int N=3, const int M=3);
/*! \brief Find the cofactor of a matrix

The cofactor matrix \f$C_{ij}\f$ is a `N-1`×`N-1` matrix which is the transpose
of a `N-1`×`N-1` matrix \f$A_{ij}\f$, which is itself the elements of a `N`×`N`
matrix \f$A\f$ with row \f$i\f$ and column \f$j\f$ j missing.

\param[out] C pointer to the first element of a `N-1`×`N-1` matrix
\param[in]  A poitner to the first element of a `N`×`N` matrix
\param      i the row of `A` to omit, `i` ∈ (0,`N`]
\param      j the column of `A` to omit, `j` ∈ (0,`N`]
\param      N the size of `A`
*/
template<typename R> void matrix_cofactor(R *C, const R *A, const int i, const int j, const int N=3);

/*! \brief Find the determinant of a matrix

\param[in] M pointer to the first element of a `N`×`N` matrix
\param     N the size of the matrix `M`
\returns the determinant, \f$|M|\f$, of `M`
*/
template<typename R> R matrix_determinant(const R *M, const int N=3);
/*! \brief Find the adjoint of a matrix

\param[out] A pointer to the first element of a `N`×`N` matrix
\param[in]  M pointer to the first element of a `N`×`N` matrix
\param      N the size of the matrix `M`
*/
template<typename R> void matrix_adjoint(R *A, const R *M, const int N=3);
/*! \brief Find the inverse of a matrix and return the determinant of the matrix

\param[out] invM pointer to the first element of a `N`×`N` matrix
\param[in]  M    pointer to the first element of a `N`×`N` matrix
\param      tol  minimum magnitude of the determinant of `M`,
                 above which `M` is invertable
\param      N    the size of the matrix `M`
\returns The determinant of `M`

\note The inverse is only calculated if the absolute value of the determinant of
      `M` is greater than `tol` to avoid dividing by zero. If the inverse is not
      calculated `invM` is unchanged by this function.
*/
template<typename R> R matrix_determinant_and_inverse(R *invM, const R *M, const R tol=0, const int N=3);
/*! \brief Find the inverse of a matrix

\param[out] invM pointer to the first element of a `N`×`N` matrix
\param[in]  M    pointer to the first element of a `N`×`N` matrix
\param      tol  minimum magnitude of the determinant of `M`,
                 above which `M` is invertable
\param      N    the size of the matrix `M`
\returns Whether the matrix was inverted

\note The inverse is only calculated if the absolute value of the determinant of
      `M` is greater than `tol` to avoid dividing by zero. If the inverse is not
      calculated `invM` is unchanged by this function.
*/
template<typename R> bool matrix_inverse(R *invM, const R *M, R tol=0, const int N=3);
/*! \brief Construct the similar matrix

\f$ M = B^{-1} A B \f$

\param[out] M   pointer to the first element of a `N`×`N` matrix
\param[in]  A   pointer to the first element of a `N`×`N` matrix
\param[in]  B   pointer to the first element of a `N`×`N` matrix
\param      tol minimum magnitude of the determinant of `B` for it to be invertible
\returns Whether a similar matrix was calculated
*/
template<typename T, int N=3> bool similar_matrix(T *M, const T *A, const T *B, const T tol=0);

/*! \brief Exchange the rows and columns of an array

\f$ D = S^T \f$

\param[out] D pointer to a `M`×`N` array
\param[in]  S pointer to a `N`×`M` array
*/
template<typename R, int N, int M> void array_transpose(R *D, const R *S);
/*! \brief Exchange the rows and columns of a matrix

\f$ D = S^T \f$

\param[out] D pointer to a `N`×`N` array
\param[in]  S pointer to a `N`×`N` array
*/
template<typename R, int N=3> void matrix_transpose(R *D, const R *S);
/*! \brief Exchange the rows and columns of a matrix in place

Performs \f$\frac{N(N-1)}{2}\f$ swaps using a single temporary scalar.

\param[in,out] B pointer to a `N`×`N` array
*/
template<typename R, int N=3> void matrix_transpose(R *B);
/*! \brief Calculate the metric of a matrix

\f$ M = L^T L \f$

\param[out] M pointer to a `N`×`N` matrix
\param[in]  L pointer to a `N`×`N` matrix
\note Allocates and deallocates a temporary `N`×`N` matrix
*/
template<typename R, int N=3> void matrix_metric(R *M, const R *L);
/*! \brief Calculate the dot product of a vector with itself

\param[in] v poitner to the first element of a `N` vector
\returns the value of `v`⋅`v` ≡ ∑ᵢvᵢvᵢ
*/
template<typename R, int N=3> R vector_norm_squared(const R * v);

template<typename T, typename R, typename S, int N> class vector_cross_impl {
public:
    static void vector_cross(T*, const R*, const S*) {
        throw std::runtime_error("The cross product is only defined for 3-vectors");
    }
};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<typename T, typename R, typename S> class vector_cross_impl<T,R,S,3>{
public:
    static void vector_cross(T * c, const R * a, const S * b) {
        c[0] = static_cast<T>(a[1])* static_cast<T>(b[2]) - static_cast<T>(a[2])* static_cast<T>(b[1]);
        c[1] = static_cast<T>(a[2])* static_cast<T>(b[0]) - static_cast<T>(a[0])* static_cast<T>(b[2]);
        c[2] = static_cast<T>(a[0])* static_cast<T>(b[1]) - static_cast<T>(a[1])* static_cast<T>(b[0]);
    }
};
#endif

/*! \brief Take the vector cross product between two 3-vectors

Perform the cross prduct betwen two vectors, \f$\vec{a} = \left(a_x, a_y, a_z\left) \f$
and \f$\vec{b}=\left(b_x,b_y,b_z\right)\f$:
\f[
\vec{a} \times \vec{b} \equiv \left(a_y b_z - b_y a_z, a_z b_x - a_x b_z, a_x b_y - b_x a_y \right)
\f]

\param[out] c The pointer to the first element of the result vector
\param[in]  a The pointer to the first element of the first vector, \f$a_x\f$
\param[in]  b The pointer to the first element of the second vecotr, \f$b_x\f$

\note The dimensionality of the vectors is defined by the template parameter `N`.
      Since the cross product is only defined in 3-dimensions this function will
      throw a `std::runtime_error` for `N` ≠ 3.
*/
template<typename T, typename R, typename S, int N = 3>
void vector_cross(T* c, const R* a, const S* b) {
    vector_cross_impl<T, R, S, N>::vector_cross(c, a, b);
}
/*! \brief Calculate the dot product of two vectors whose size is only known at runtime

\param     n the size of the two vectors
\param[in] a pointer to the first element of a `n` vector
\param[in] b pointer to the first element of a `n` vector
\returns the value of `a`⋅`b` ≡ ∑ᵢⁿaᵢbᵢ
*/
template<typename S,typename T, typename R> S vector_dot(const size_t n, const T* a, const R* b);
/*! \brief Calculate the dot product of two vectors whose size is known at compilation

\param[in] a pointer to the first element of a `N` vector
\param[in] b pointer to the first element of a `N` vector
\returns the value of `a`⋅`b` ≡ ∑ᵢⁿaᵢbᵢ
*/
template<typename R, int N=3> R vector_dot(const R *a, const R *b);
//! Return the remainder of a/1, ensuring the result is within the range (0,1]
template<typename T> T mod1(const T a);

/*! \brief Determine if a floating point matrix holds only integer values

\param[in] A   pointer to the first element of a `N`×`N` matrix
\param     tol tolerance within which a value can deviate from its closest exact integer
\returns whether all elements of `A` are no more than `tol` from their closest
         exact integer value.
*/
template<typename T, typename R, int N=3> bool is_int_matrix(const T * A, const R tol);
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<typename R> bool is_int_matrix(const int *, const R);
#endif

/*! \brief The "distance" between two matrices using the Frobenius norm

The [Frobenius norm](http://mathworld.wolfram.com/FrobeniusNorm.html)
is the matrix norm of an 𝑚×𝑛 matrix 𝑎 given by

    |A|ᶠ = √∑ᵢᵐ∑ⱼⁿ|𝑎ᵢⱼ|²

or, for a square matrix,

    |A|ᶠ = √tr(𝑎𝑎ᴴ)

where 𝑎ᴴ is the conjugate transpose of 𝑎.
This function calculates the Frobenius norm of the matrix A-B.
@param n The size of the matrices
@param A A n×n square matrix
@param B A n×n square matrix
@returns |A|ᶠ
*/
template<class I, class T> T frobenius_distance(const I n, const T* A, const T* B);
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class I, class T> T frobenius_distance(const I n, const std::complex<T>* A, const std::complex<T>* B);
#endif

/*! \brief The general n-dimensional angle between two real-valued vectors

Calculate and return the angle θ between two vectors given by

    cos(θ) = <A,B>/|A||B|

@param n The dimensionality of the vector space
@param A A pointer to the first vector
@param B A pointer to the second vector
@returns θ
*/
template<class I, class T> T vector_angle(const I n, const T* A, const T* B);
/*! \brief The general 2n-dimensional angle between two complex-valued n-vectors

@param n The dimensionality of the vector space
@param A A pointer to the first vector
@param B A pointer to the second vector
@see euclidian_angle
*/
template<class I, class T> T vector_angle(const I n, const std::complex<T>* A, const std::complex<T>* B);
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
template<class I, class T> T euclidean_angle(const I n, const std::complex<T>* A, const std::complex<T>* B);
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
template<class I, class T> T hermitian_angle(const I n, const std::complex<T>* A, const std::complex<T>* B);
/*! \brief The general n-dimensional Hermitian angle between two real-valued vectors

@param n The dimensionality of the vector space
@param A A pointer to the first vector
@param B A pointer to the second vector
@returns θₕ
@see vector_angle
*/
template<class I, class T> T hermitian_angle(const I n, const T* A, const T* B);
/*! \brief The general n-dimensional Hermitian product betwen two real-valued vectors

\param n the dimensionality of the vector space
\param a poitner to the a real-valued `n` vector
\param b poitner to the a real-valued `n` vector
\returns `a`⋅`b`
*/
template<class I, class T> T hermitian_product(const I n, const T* a, const T* b);
/*! \brief The general n-dimensional Hermitian product betwen a real- and complex-valued vector

\param n the dimensionality of the vector space
\param a poitner to the a real-valued `n` vector
\param b poitner to the a complex-valued `n` vector
\returns `a`⋅`b`
*/
template<class I, class T> std::complex<T> hermitian_product(const I n, const T* a, const std::complex<T>* b);
/*! \brief The general n-dimensional Hermitian product betwen a complex- and real-valued vector

\param n the dimensionality of the vector space
\param a poitner to the a complex-valued `n` vector
\param b poitner to the a real-valued `n` vector
\returns `std::conj(a)`⋅`b`
*/
template<class I, class T> std::complex<T> hermitian_product(const I n, const std::complex<T>* a, const T* b);
/*! \brief The general n-dimensional Hermitian product betwen two complex--valued vectors

\param n the dimensionality of the vector space
\param a poitner to the a complex-valued `n` vector
\param b poitner to the a complex-valued `n` vector
\returns `std::conj(a)`⋅\`b`
*/
template<class I, class T> std::complex<T> hermitian_product(const I n, const std::complex<T>* a, const std::complex<T>* b);

/*! \brief Calculate the length of the vector connecting two real points

\param n the dimensionality of the vector space
\param a pointer to a real-valued `n` vector
\param b pointer to a real-valued `n` vector
\returns |`a`-`b`|
*/
template<class I, class T> T vector_distance(const I n, const T* a, const T* b);
/*! \brief Calculate the length of the vector connecting two complex points

\param n the dimensionality of the vector space
\param a pointer to a complex-valued `n` vector
\param b pointer to a complex-valued `n` vector
\returns |`a`-`b`|
*/
template<class I, class T> T vector_distance(const I n, const std::complex<T>* a, const std::complex<T>* b);
/*! \brief Calculate the vector product of two real-valued vectors

\param n the dimensionality of the vector space
\param a pointer to a real-valued `n` vector
\param b pointer to a real-valued `n` vector
\returns `a`⋅`b`
*/
template<class I, class T> T vector_product(const I n, const T* a, const T* b);
/*! \brief Calculate the vector product of a real- and complex-valued vector

\param n the dimensionality of the vector space
\param a pointer to a real-valued `n` vector
\param b pointer to a complex-valued `n` vector
\returns the real part of `a`⋅`b`
*/
template<class I, class T> T vector_product(const I n, const T* a, const std::complex<T>* b);
/*! \brief Calculate the vector product of a complex- and real-valued vector

\param n the dimensionality of the vector space
\param a pointer to a complex-valued `n` vector
\param b pointer to a real-valued `n` vector
\returns the real part of `std::conj(a)`⋅`b`
*/
template<class I, class T> T vector_product(const I n, const std::complex<T>* a, const T* b);
/*! \brief Calculate the vector product of a complex-valued vectors

\param n the dimensionality of the vector space
\param a pointer to a complex-valued `n` vector
\param b pointer to a complex-valued `n` vector
\returns the real part of `std::conj(a)`⋅`b`
*/
template<class I, class T> T vector_product(const I n, const std::complex<T>* a, const std::complex<T>* b);

/*! \brief Find the squared difference between two values in the complex plane

The distance between two real values in the complex plane is the difference
of their values.
\param A the first value
\param B the second value
\returns |`A`-`B`|²
*/
template<class T> T squared_distance(const T&A, const T& B);
/*! \brief Find the squared difference between two values in the complex plane

\param A the first complex value
\param B the second complex value
\returns |`A`-`B`|²
*/
template<class T> T squared_distance(const std::complex<T>& A, const std::complex<T>& B);
/*! \brief Find the distance of a value from the origin of the complex plane

\param a the point on the real-axis to find the magnitude of
\returns |`a`|
*/
template<class T> T magnitude(const T a);
/*! \brief Find the distance of a value from the origin of the complex plane

\param a the point in the complex plane to find the magnitude of
\returns |`a`| ≡ `std::sqrt(std::real(std::conj(a)*a))`
*/
template<class T> T magnitude(const std::complex<T> a);


//! Convert unsigned integers to signed integers
template<typename S,typename U> S u2s(const U u);
//! Convert signed integers to unsigned integers
template<typename U,typename S> U s2u(const S s);

#include "utilities.tpp"

  } // namespace utils
} // namespace brille


#endif
