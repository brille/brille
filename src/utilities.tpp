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

// trace of a square matrix
template<typename T, int N> T trace(const T *M){
  T out = T(0);
  for (int i=0; i<N; i++) out += M[i*(1+N)];
  return out;
}

// copying array from source (S) to destination (D)
template<typename T, int N, int M> void copy_array(T *D, const T *S){ for (int i=0; i<N*M; i++) D[i]=S[i]; }
template<typename T, int N> void copy_matrix(T *M, const T *A){ copy_array<T,N,N>(M,A); }
template<typename T, int N> void copy_vector(T *V, const T *A){ copy_array<T,N,1>(V,A); }

// element-wise checking of approximate equality
template<typename T, int N, int M> bool equal_array(const T *A, const T *B, const T tol){
  for (int i=0; i<N*M; i++) if ( std::abs(A[i]-B[i]) > tol ) return false;
  return true;
}
template<typename T, int N> bool equal_matrix(const T *A, const T *B, const T tol){ return equal_array<T,N,N>(A,B,tol); }
template<typename T, int N> bool equal_vector(const T *A, const T *B, const T tol){ return equal_array<T,N,1>(A,B,tol); }

// array multiplication C = A * B -- where C is (N,M), A is (N,I) and B is (I,M)
template<typename T, typename R, typename S, size_t N, size_t I, size_t M> void multiply_arrays(T *C, const R *A, const S *B){
  for (size_t i=0;i<N*M;i++) C[i]=T(0);
  for (size_t i=0;i<N;i++) for (size_t j=0;j<M;j++) for (size_t k=0;k<I;k++) C[i*M+j] += T(A[i*I+k]*B[k*M+j]);
}
template<typename T, typename R, typename S, size_t N> void multiply_matrix_matrix(T *C, const R *A, const S *B){ multiply_arrays<T,R,S,N,N,N>(C,A,B); }
template<typename T, typename R, typename S, size_t N> void multiply_matrix_vector(T *C, const R *A, const S *b){ multiply_arrays<T,R,S,N,N,1>(C,A,b); }
template<typename T, typename R, typename S, size_t N> void multiply_vector_matrix(T *C, const R *a, const S *B){ multiply_arrays<T,R,S,1,N,N>(C,a,B); }


// array multiplication specialization for non-complex * complex arrays.
template<class T, class R, class S, class I>
// std::enable_if_t<std::is_unsigned_v<I>>
void
mul_arrays(T* C, const I n, const I l, const I m, const R* A, const S* B){
  for (I i=0;i<n*m;i++) C[i]=T(0);
  for (I i=0;i<n;i++) for (I j=0;j<m;j++) for (I k=0;k<l;k++) C[i*m+j] += static_cast<T>(A[i*l+k]*B[k*m+j]);
}
template<class T, class R, class S, class I>
// std::enable_if_t<std::is_unsigned_v<I>>
void
mul_arrays(std::complex<T>* C, const I n, const I l, const I m, const R* A, const std::complex<S>* B){
  for (I i=0;i<n*m;i++) C[i]=std::complex<T>(0);
  for (I i=0;i<n;i++) for (I j=0;j<m;j++) for (I k=0;k<l;k++) C[i*m+j] += static_cast<S>(A[i*l+k])*B[k*m+j];
}
template<class T, class R, class S, class I>
// std::enable_if_t<std::is_unsigned_v<I>>
void
mul_arrays(std::complex<T>* C, const I n, const I l, const I m, const std::complex<R>* A, const S* B){
  for (I i=0;i<n*m;i++) C[i]=std::complex<T>(0);
  for (I i=0;i<n;i++) for (I j=0;j<m;j++) for (I k=0;k<l;k++) C[i*m+j] += A[i*l+k]*static_cast<R>(B[k*m+j]);
}
template<class T, class R, class S, class I>
// std::enable_if_t<std::is_unsigned_v<I>>
void
mul_arrays(std::complex<T>* C, const I n, const I l, const I m, const std::complex<R>* A, const std::complex<S>* B){
  for (I i=0;i<n*m;i++) C[i]=std::complex<T>(0);
  for (I i=0;i<n;i++) for (I j=0;j<m;j++) for (I k=0;k<l;k++) C[i*m+j] += A[i*l+k]*B[k*m+j];
}


// array element-wise addition
template<typename T, typename R, typename S, int N, int M> void add_arrays(T *C, const R *A, const S *B){ for (int i=0; i<N*M; i++) C[i] = T(A[i]+B[i]); }
template<typename T, typename R, typename S, int N> void add_matrix(T *C, const R *A, const S *B){ add_arrays<T,R,S,N,N>(C,A,B); }
// array element-wise subtraction
template<typename T, typename R, typename S, int N, int M> void subtract_arrays(T *C, const R *A, const S *B){ for (int i=0; i<N*M; i++) C[i] = T(A[i]-B[i]); }
template<typename T, typename R, typename S, int N> void subtract_matrix(T *C, const R *A, const S *B){ subtract_arrays<T,R,S,N,N>(C,A,B); }

// specialized casting of doubles to ints
template<typename T,typename R> T my_cast(const R a){ return T(a); }
// inlined to avoid redefining this specialization each time the header is called
template<> inline int my_cast<int,double>(const double a){ return int( a<0 ? a-0.5 : a+0.5); }
// element-wise re-(special)-casting of arrays
template<typename T, typename R, int N, int M> void cast_array(T *A, const R *B){ for (int i=0; i<N*M; i++) A[i] = my_cast<T,R>(B[i]); }

template<typename T, typename R, int N> void cast_matrix(T *A, const R *B){ cast_array<T,R,N,N>(A,B); }
template<typename T, typename R, int N> void cast_vector(T *A, const R *B){ cast_array<T,R,N,1>(A,B); }


// The cofactor array C_ij of the N by M array A is a M-1 by N-1 array of transpose(A with row i and column j) missing:
template<typename R> void array_cofactor(R *C, const R *A, const int i, const int j, const int N, const int M){
  int k = 0;
  for (int m=0; m<M; m++) for (int n=0; n<N; n++) if (i!=n && j!=m) C[k++] = A[n+N*m];
}
template<typename R> void matrix_cofactor(R *C, const R *A, const int i, const int j, const int N){ array_cofactor<R>(C,A,i,j,N,N); }

// the determinant is only defined for square matrices
template<typename R> R matrix_determinant(const R *M, const int N){
  if (1==N) return M[0];
  // the determinant:
  R d{0};
  // temporary cofactor storage
  R *cof = nullptr;
  cof = new R[(N-1)*(N-1)]();

  // loop over one row/column (let's go with row, for fun)
  R pm = 1; // +/-1 for alternating elements
  for (int i=0; i<N; i++){
    matrix_cofactor(cof,M,0,i,N);
    d += pm * M[i*N] * matrix_determinant(cof,N-1);
    pm *= -1;
  }
  delete[] cof;
  return d;
}

// the adjoint is also only defined for a square matrix
template<typename R> void matrix_adjoint(R *A, const R *M, const int N){
  if (1==N){
    A[0]=R(1);
    return;
  }
  R pm=1, *cof =nullptr;
  cof = new R[(N-1)*(N-1)]();
  for (int i=0; i<N; i++) for (int j=0; j<N; j++){
      // the cofactor matrix for M[i,j]
      matrix_cofactor(cof,M,i,j,N);
      // if i+j is even, multiply the determinant by +1, otherwise -1
      pm = ( (i+j)%2==0 ) ? R(1) : R(-1);
      // We want to save the transpose, so A[j+i*n]
      A[j+i*N] = pm * matrix_determinant(cof,N-1);
  }
  delete[] cof;
}

// the inverse is, yet again, only defined for square matrices (with non-zero determinants)
template<typename R> R matrix_determinant_and_inverse(R *invM, const R *M, const R tol, const int N){
  R d = matrix_determinant(M,N);
  if (std::abs(d) > tol ){
    // the matrix is *not* singular and has an inverse
    matrix_adjoint(invM,M,N); // at this point invM is the Adjoint(M)
    // The inv(M) = adjoint(M)/det(M)
    int nn = N*N;
    for (int i=0; i<nn; i++) invM[i] /= d;
    // now invM is inv(M)!
  }
  return d;
}
template<typename R> bool matrix_inverse(R *invM, const R *M, R tol, const int N){ return ( std::abs( matrix_determinant_and_inverse(invM, M, tol, N) ) > tol) ; }

template<typename T, int N> bool similar_matrix(T *M, const T *A, const T *B, const T tol){
  T *C = nullptr;
  C = new T[N*N]();
  bool ok = matrix_inverse(C,B,tol,N);
  if ( ok ){
    T *P = nullptr;
    P = new T[N*N]();
    multiply_matrix_matrix<T,T,T,N>(P,A,B);
    multiply_matrix_matrix<T,T,T,N>(M,C,P);
    delete[] P;
  } else {
    printf("spglib: No similar matrix due to 0 determinant.\n");
  }
  delete[] C;
  return ok;
}

template<typename R, int N, int M> void array_transpose(R *D, const R *S){
  for (int i=0; i<N; ++i) for (int j=0; j<M; ++j)  D[i+j*N] = S[j+i*M];
}
template<typename R, int N> void matrix_transpose(R *D, const R *S){ array_transpose<R,N,N>(D,S); }
//
template<typename R, int N, int M> void complex_conj_array_transpose(std::complex<R> *D, const std::complex<R> *S){
  for (int i=0; i<N; ++i) for (int j=0; j<M; ++j)  D[i+j*N] = std::conj(S[j+i*M]);
}
template<typename R, int N> void complex_conj_matrix_transpose(std::complex<R> *D, const std::complex<R> *S){ complex_conj_array_transpose<R,N,N>(D,S); }
// in place transpose:
template<typename R, int N> void matrix_transpose(R *B){
  R t;
  for (int i=0; i<N; i++) for(int j=i+1; j<N; j++) {
    t = B[i+j*3];
    B[i+j*3] = B[j+i*3];
    B[j+i*3] = t;
  }
}

template<typename R, int N> void matrix_metric(R *M, const R *L){
  R *Lt = nullptr;
  Lt = new R[N*N]();
  matrix_transpose<R,N>(Lt,L);
  multiply_matrix_matrix<R,R,R,N>(M,Lt,L);
  delete[] Lt;
}

template<typename R, int N> R vector_norm_squared(const R *v){
  R vv = 0;
  for (int i=0; i<N; i++) vv += v[i]*v[i];
  return vv;
}
//template<typename T, typename R, typename S, int N> void vector_cross(T *c, const R *a, const S *b){
//  if (3!=N)
//    throw std::domain_error("The cross product is only defined for 3-vectors");
//  c[0] = static_cast<T>(a[1])*static_cast<T>(b[2]) - static_cast<T>(a[2])*static_cast<T>(b[1]);
//  c[1] = static_cast<T>(a[2])*static_cast<T>(b[0]) - static_cast<T>(a[0])*static_cast<T>(b[2]);
//  c[2] = static_cast<T>(a[0])*static_cast<T>(b[1]) - static_cast<T>(a[1])*static_cast<T>(b[0]);
//}
template<typename S, typename T, typename R>
S vector_dot(const size_t n, const T* a, const R* b){
  S out = 0;
  for (size_t i=0; i<n; ++i) out += static_cast<S>(a[i])*static_cast<S>(b[i]);
  return out;
}

template<typename R, int N> R vector_dot(const R *a, const R *b){
  R out = 0;
  for (int i=0; i<N; i++) out += a[i]*b[i];
  return out;
}

template<typename T> T mod1(const T a){
  T b = a - my_cast<int,T>(a);
  return ( (b < T(0) - 1e-10 ) ? b + T(1) : b );
}

template<typename T, typename R, int N> bool is_int_matrix(const T * A, const R tol){
  for (int i=0; i<N*N; i++) if ( std::abs(my_cast<int,T>(A[i]) - A[i]) > tol ) return false;
  return true;
}
template<typename R> bool is_int_matrix(const int *, const R){ return true; }


template<class I, class T>
T
frobenius_distance(const I n, const T* A, const T* B) {
  // sqrt(tr([A-B]*[A-B]')) == sqrt(A00^2 + A01^2 + A02^2 ... Ann^2)
    return vector_distance(n*n, A, B);
}

template<class I, class T>
T
frobenius_distance(const I n, const std::complex<T>* A, const std::complex<T>* B) {
    // sqrt( tr( (A-B)x(A-B)') ) is identically the same as the vector distance of 9-vectors
    return vector_distance(n*n, A, B);
}


template<class I, class T> T vector_angle(const I n, const T* A, const T* B){
  T AA=0, BB=0, AB=0;
  for (I i=0; i<n; ++i){
    AA += A[i]*A[i];
    BB += B[i]*B[i];
    AB += A[i]*B[i];
  }
  T nA, nB, c_t;
  nA = std::sqrt(AA);
  nB = std::sqrt(BB);
  if (nA && nB){
    c_t = AB/(nA*nB);
  } else {
    // if both are zero their angle is strictly undefined, but for our
    // purposes it would be better to set the angle to 0 [and cos(0)=1]
    c_t = (nA || nB) ? T(0) : T(1);
  }
  T act = std::abs(c_t);
  if (brille::approx::scalar(act, 1.0) && act>1){
    c_t /= act;
    act = std::abs(c_t);
  }
  if (act>1){
    std::string msg = "vector angle cos(theta)="
    + my_to_string(c_t) + " = " + my_to_string(AB) + "/("
    + my_to_string(nA) + "×" + my_to_string(nB) + ")";
    throw std::runtime_error(msg);
  }
  return std::acos(c_t);
}

template<class I, class T> T vector_angle(const I n, const std::complex<T>* A, const std::complex<T>* B){
  // return hermitian_angle(n,A,B);
  return euclidean_angle(n,A,B);
}

template<class I, class T> T euclidean_angle(const I n, const std::complex<T>* A, const std::complex<T>* B){
  T AB{0}, nA{0}, nB{0}, c_t;
  // Compute the products of complex n-vectors as if they were real 2n-vectors
  for (I i=0; i<n; ++i){
    AB += std::real(A[i])*std::real(B[i]) + std::imag(A[i])*std::imag(B[i]);
    nA += std::real(A[i])*std::real(A[i]) + std::imag(A[i])*std::imag(A[i]);
    nB += std::real(B[i])*std::real(B[i]) + std::imag(B[i])*std::imag(B[i]);
  }
  nA = std::sqrt(nA);
  nB = std::sqrt(nB);
  if (nA && nB){
    c_t = AB/(nA*nB);
  } else {
    // if both are zero their angle is strictly undefined, but for our
    // purposes it would be better to set the angle to 0 [and cos(0)=1]
    c_t = (nA || nB) ? T(0) : T(1);
  }
  T act = std::abs(c_t);
  if (brille::approx::scalar(act, 1.0) && act>1){
    c_t /= act;
    act = std::abs(c_t);
  }
  if (act>1){
    std::string msg = "Complex Euclidean angle cos(theta)="
    + my_to_string(c_t) + " = " + my_to_string(AB) + "/("
    + my_to_string(nA) + "×" + my_to_string(nB) + ")";
    throw std::runtime_error(msg);
  }
  return std::acos(c_t);
}

template<class I, class T> T hermitian_product(const I n, const T* a, const T* b){
  T h_dot{0};
  for (I i=0; i<n; ++i) h_dot += a[i]*b[i];
  return h_dot;
}
template<class I, class T> std::complex<T> hermitian_product(const I n, const T* a, const std::complex<T>* b){
  std::complex<T> h_dot{0,0};
  for (I i=0; i<n; ++i) h_dot += a[i]*b[i];
  return h_dot;
}
template<typename T> std::complex<T> complex_prod(const std::complex<T>& x, T y){
  return std::conj(x)*y;
}
template<typename T> std::complex<T> complex_prod(const std::complex<T>& x, const std::complex<T>& y){
  T xr{x.real()}, xi{x.imag()}, yr{y.real()}, yi{y.imag()};
  // a*×b = (a'+ia")*×(b'+ib") = (a'-ia")×(b'+ib") = a'b' + a"b" + i(a'b" -a"b')
  return std::complex<T>(xr*yr+xi*yi, xr*yi-xi*yr);
}
template<class I, class T> std::complex<T> hermitian_product(const I n, const std::complex<T>* a, const T* b){
  std::complex<T> h_dot{0,0};
  // for (size_t i=0; i<n; ++i) h_dot += std::conj(a[i])*b[i];
  for (I i=0; i<n; ++i) h_dot += complex_prod(a[i], b[i]);
  return h_dot;
}
template<class I, class T> std::complex<T> hermitian_product(const I n, const std::complex<T>* a, const std::complex<T>* b){
  T hr{0}, hi{0};
  long long sn = u2s<long long>(n);
  #pragma omp parallel for shared(a,b) reduction(+: hr,hi)
  for (long long si=0; si<sn; ++si){
    std::complex<T> h = std::conj(a[si])*b[si];
    hr += h.real();
    hi += h.imag();
  }
  return std::complex<T>(hr,hi);
}



template<class I, class T> T hermitian_angle(const I n, const T* A, const T* B){
  return vector_angle(n,A,B);
}
template<class I, class T> T hermitian_angle(const I n, const std::complex<T>* A, const std::complex<T>* B){
  T nAB, nA, nB, c_t;
  nAB = std::sqrt(vector_product(n,A,B));
  nA  = std::sqrt(std::real(hermitian_product(n,A,A)));
  nB  = std::sqrt(std::real(hermitian_product(n,B,B)));
  if (nA && nB){
    c_t = nAB/(nA*nB);
  } else {
    // if both are zero their angle is strictly undefined, but for our
    // purposes it would be better to set the angle to 0 [and cos(0)=1]
    c_t = (nA || nB) ? T(0) : T(1);
  }
  T act = std::abs(c_t);
  if (brille::approx::scalar(act,1.0) && act>1){
    c_t/=act; // force close-to-one values to one, maintaining the sign
    act = std::abs(c_t);
  }
  if (act>1){
    std::string msg = "Complex Hermitian angle cos(theta)="
    + my_to_string(c_t) + " = " + my_to_string(nAB) + "/("
    + my_to_string(nA) + "×" + my_to_string(nB) + ")";
    throw std::runtime_error(msg);
  }
  return std::acos(c_t);
}

template<class I, class T> T vector_distance(const I n, const T* a, const T* b){
  T d, s=0;
  for (I i=0; i<n; ++i){
    d = a[i]-b[i];
    s += d*d;
  }
  return std::sqrt(s);
}
template<class I, class T> T vector_distance(const I n, const std::complex<T>* a, const std::complex<T>* b){
  std::complex<T> d;
  T s=0;
  for (I i=0; i<n; ++i){
    d = a[i]-b[i];
    s += std::real(d*std::conj(d));
  }
  return std::sqrt(s);
}

template<class I, class T> T vector_product(const I n, const T* a, const T* b){
  T h_dot{0};
  for (size_t i=0; i<n; ++i) h_dot += a[i]*b[i];
  return h_dot;
}
template<class I, class T> T vector_product(const I n, const T* a, const std::complex<T>* b){
  std::complex<T> h_dot = hermitian_product(n,a,b);
  return std::real(h_dot*std::conj(h_dot));
}
template<class I, class T> T vector_product(const I n, const std::complex<T>* a, const T* b){
  std::complex<T> h_dot = hermitian_product(n,a,b);
  return std::real(h_dot*std::conj(h_dot));
}
template<class I, class T> T vector_product(const I n, const std::complex<T>* a, const std::complex<T>* b){
  std::complex<T> h_dot = hermitian_product(n,a,b);
  return std::real(h_dot*std::conj(h_dot));
}

template<class I, class T> T inner_product(const I n, const T* a, const T* b){
  T h_dot{0};
  for (I i=0; i<n; ++i) h_dot += a[i]*b[i];
  return h_dot;
}
template<class I, class T> T inner_product(const I n, const std::complex<T>* a, const std::complex<T>* b){
  std::complex<T> h_dot = hermitian_product(n,a,b);
  return std::real(h_dot);
}

template<class T> T squared_distance(const T&A, const T& B){
  return (A-B)*(A-B);
}
template<class T> T squared_distance(const std::complex<T>& A, const std::complex<T>& B){
  T r = std::real(A)-std::real(B);
  T i = std::imag(A)-std::imag(B);
  return r*r + i*i;
}

template<class T> T magnitude(const T a){ return std::abs(a); }
template<class T> T magnitude(const std::complex<T> a){
  return std::sqrt(std::real(a*std::conj(a)));
}


/*! Determine an unsigned integer that encodes the signs of a complex vector.

Each element of a complex valued vector has two signs associated with it, one
for the real part and one for the imaginary part. Some properties of the vector,
notably the squared modulus of its dot product with a real vector, depend on
the relative signs of each element. This function encodes the four possible
combinations into a quaternary based system, (++,+-,-+ ,--)→(0,1,2,3)₄ and
combines the `n` values into a single unsigned integer using 2`n` bits.

As an example, a vector (+a-𝑖b, -c+𝑖d, -e-𝑖f) has signs (++,-+,--) which give
the quaternary numbers (0,2,3)₄ ⇒ 001011₂ = 9.
*/
template<class T> size_t encode_array_signs(const size_t n, const std::complex<T>* a){
  size_t signs=0;
  for (size_t i=0; i<n; ++i){
    if (std::signbit(std::imag(a[i]))) signs += 1 << 2*(n-1-i);
    if (std::signbit(std::real(a[i]))) signs += 2 << 2*(n-1-i);
  }
  return signs;
}
/*! Determine an unsigned integer that encodes the signs of a real vector

Less useful than it's complex-valued pair. Only really defined in case some
system produces strictly-real eigenvectors.
*/
template<class T> size_t encode_array_signs(const size_t n, const T* a){
  size_t signs=0;
  for (size_t i=0; i<n; ++i){
    if (std::signbit(a[i])) signs += 2 << 2*(n-1-i);
  }
  return signs;
}


template<class T> int make_eigenvectors_equivalent(const size_t n, const T* v0, T* v1){
  size_t s0, s1;
  s0 = encode_array_signs(n,v0);
  s1 = encode_array_signs(n,v1);
  if (s0 == s1) return 0;
  // the only valid permutation for real values is by 2.
  size_t onesign, s1mod{0};
  for (size_t j=0; j<n; ++j){
    // extract a single sign quaternary
    onesign = (s1 >> 2*(n-1-j)) - ((s1 >> 2*(n-j)) << 2*(n-j));
    // permute it by 2
    onesign = (onesign+2)%4;
    // and stick it in s1mod
    s1mod += onesign << 2*(n-1-j);
  }
  if (s0 == s1mod){
    T m1{-1};
    for (size_t j=0; j<n; ++j) v1[j] *= m1;
    return 0;
  }
  return 1;
}
template<class T> int make_eigenvectors_equivalent(const size_t n, const std::complex<T>* v0, std::complex<T> v1){
  size_t s0, s1;
  s0 = encode_array_signs(n,v0);
  s1 = encode_array_signs(n,v1);
  if (s0 == s1) return 0;
  // The signs of each vector are not the same, so check for equivalence:
  size_t onesign;
  // we can permute each sign quaternary number up to three times
  for (size_t i=1; i<4; ++i){
    size_t s1mod{0};
    for (size_t j=0; j<n; ++j){
      // extract a single sign quaternary
      onesign = (s1 >> 2*(n-1-j)) - ((s1 >> 2*(n-j)) << 2*(n-j));
      // permute it
      onesign = (onesign+i)%4;
      // and stick it in s1mod
      s1mod += onesign << 2*(n-1-j);
    }
    if (s0 == s1mod){
      // The two are equivalent for our purposes, but the signs of v1 need to be changed!
      T m1{-1};
      switch (i){
        case 1: // exchange the imaginary sign
        for (size_t j=0; j<n; ++j) v1[j] = std::complex<T>(std::real(v1[j]), m1*std::imag(v1[j]));
        break;
        case 2: // exchange both real and imaginary signs
        for (size_t j=0; j<n; ++j) v1[j] = std::complex<T>(m1*std::real(v1[j]), m1*std::imag(v1[j]));
        break;
        case 3: // exchange the real sign
        for (size_t j=0; j<n; ++j) v1[j] = std::complex<T>(m1*std::real(v1[j]), std::imag(v1[j]));
        break;
      }
      return 0;
    }
  }
}


template<class T> T antiphase(const T){
  //return std::signbit(z) ? T(-1) : T(1);
  return T(1);
}
template<class T> std::complex<T> antiphase(const std::complex<T> z){
  return std::polar(T(1),-std::arg(z));
}

template<class I, class T>
T
antiphase(const I, const T*, const T*){
  return T(1);
}
template<class I, class T>
std::complex<T> antiphase(const I n, const std::complex<T>* a, const std::complex<T>* b){
  T real_dot{0}, imag_dot{0};
  for (I i=0; i<n; ++i){
    T areal{a[i].real()}, aimag{a[i].imag()}, breal{b[i].real()}, bimag{b[i].imag()};
    real_dot += areal * breal + aimag * bimag;
    imag_dot += areal * bimag - aimag * breal;
  }
  return std::polar(T(1), T(-1)*std::atan2(imag_dot, real_dot));
}

template<class I, class T>
void
inplace_antiphase(const I, const T*, const T*, T*)
{}
template<class I, class T>
void
inplace_antiphase(const I n, const std::complex<T>* a, const std::complex<T>* b, std::complex<T>* phased)
{
  std::complex<T> eith = antiphase(n, a, b);
  for (I i=0; i<n; ++i) phased[i] = eith*b[i];
}

template<class T, template<class> class A>
T
unsafe_antiphase(const A<T>&, const A<T>&){
  return T(1);
}

template<class T, template<class> class A>
std::complex<T>
unsafe_antiphase(const A<std::complex<T>>& a, const A<std::complex<T>>& b){
  T real_dot{0}, imag_dot{0};
  // a safe version of this would first ensure that a and b are the same shape
  for (auto [ox, ax, bx]: a.broadcastItr(b)){
    std::complex<T> av{a[ax]}, bv{b[bx]};
    T ar{av.real()}, ai{av.imag()}, br{bv.real()}, bi{bv.imag()};
    real_dot += ar*br + ai*bi;
    imag_dot += ar*bi - ai*br;
  }
  return std::polar(T(1), T(-1)*std::atan2(imag_dot, real_dot));
}

template<class T, class R, class S = typename std::common_type<T,R>::type>
S coth_over_en(const T en, const R beta){
  S Sen = static_cast<S>(en);
  S Sbeta = static_cast<S>(beta);
  S den = std::tanh(Sen*Sbeta/S(2))*Sen;
  return S(1)/den;
}
template<class T, class R, class S = typename std::common_type<T,R>::type>
S coth_over_en(const std::complex<T> en, const R beta){
  S den = static_cast<S>(std::real(std::tanh(en*beta*0.5)*en));
  return S(1)/den;
}

template<class T>
std::enable_if_t<std::is_unsigned<T>::value, T>
gcd(const T a, const T b){
  // We include numeric for std::gcd in C++17
  /* Pre-C++17 gcd was not part of the standard but may be present in an
  // experimental header -- which of course doesn't exist for MSVC        */
  #if defined(__cplusplus) && __cplusplus > 201700L
    return std::gcd(a,b);
  #else
    if (b==0) return a;
    return gcd(b, a % b);
  #endif
}

template<class T, class R>
std::enable_if_t<std::is_unsigned<T>::value&&std::is_unsigned<R>::value, unsigned long long>
binomial_coefficient(const T n, const R k){
  unsigned long long ans{1}, num{1}, den{1}, comdiv;
  if (k>n){
    std::string msg = "The binomial coefficient requires more choices than selections.";
    msg += "\n binomial_coefficient(" + std::to_string(n) + ":" + std::to_string(k) + ") is invalid";
    throw std::domain_error(msg);
  }
  // the Binomial coefficient is symmetric due to the denominator k!(n-k)!
  R m = (n-k < k) ? static_cast<R>(n-k) : k;
  bool overflow = false;
  for (R i=0; i<m; ++i){
    unsigned long long lastnum{num}, lastden{den};
    num *= static_cast<unsigned long long>(n-i);
    den *= static_cast<unsigned long long>(i)+1;
    if (lastnum > num || lastden > den){
      comdiv = gcd(lastnum, lastden);
      if (comdiv > 1){
        num = lastnum/comdiv;
        den = lastden/comdiv;
        --i;
      } else {
        overflow = true;
      }
    }
    if (overflow) break;
  }

  if (overflow){
    long double dans{1};
    for (T i=0; i<m; ++i)
      dans *= static_cast<long double>(n-i)/(static_cast<long double>(i)+1);
    ans = std::llround(dans);
  } else {
    ans = num/den;
  }
  return ans;
}


template<typename S,typename U>
S u2s(const U u){
  S smax = (std::numeric_limits<S>::max)();
  U usmax = static_cast<U>(smax);
  if (u > usmax){
    std::string msg = "unsigned_to_signed:: Value " + std::to_string(u)
                    + " can not be stored in requested signed type"
                    + " (max value "+ std::to_string(smax) + ")";
    throw std::overflow_error(msg);
  }
  return static_cast<S>(u);
}

template<typename U,typename S>
U s2u(const S s){
  if (s < 0
    //|| static_cast<U>((std::numeric_limits<S>::max)()) > (std::numeric_limits<U>::max)() // the largest element of type S is too big for U
    //|| s > static_cast<S>((std::numeric_limits<U>::max)()) // s, specifically is too big to be expressed in type U
  ){
    std::string msg = "signed_to_unsiged:: Value " + std::to_string(s)
                    + " can not be stored in requested unsiged type"
                    + " (max value "+ std::to_string((std::numeric_limits<U>::max)()) + ")";
    throw std::overflow_error(msg);
  }
  return static_cast<U>(s);
}


// from https://en.cppreference.com/w/cpp/language/sizeof...
template<typename... Ts>
constexpr auto make_array(Ts&&... ts)
    -> std::array<std::common_type_t<Ts...>,sizeof...(Ts)>
{
    return { std::forward<Ts>(ts)... };
}

template<typename... Ts>
constexpr auto make_vector(Ts&&... ts)
    -> std::vector<std::common_type_t<Ts...>>
{
    return { std::forward<Ts>(ts)... };
}

template<class Head, class... Tail>
using are_same = std::conjunction<std::is_same<Head, Tail>...>;
