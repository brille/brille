#define ZERO_PREC 1e-10

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

// absolute value of a scalar
template<typename T> T my_abs(const T a){ return (a<0 ? T(-1)*a : a); }

// element-wise checking of approximate equality
template<typename T, int N, int M> bool equal_array(const T *A, const T *B, const T tol){
  for (int i=0; i<N*M; i++) if ( my_abs(A[i]-B[i]) > tol ) return false;
  return true;
}
template<typename T, int N> bool equal_matrix(const T *A, const T *B, const T tol){ return equal_array<T,N,N>(A,B,tol); }
template<typename T, int N> bool equal_vector(const T *A, const T *B, const T tol){ return equal_array<T,N,1>(A,B,tol); }

template<typename T, typename R> bool approx_scalar(const T a, const R b){
  bool isfpT, isfpR;
  isfpT = std::is_floating_point<T>::value;
  isfpR = std::is_floating_point<R>::value;
  T Ttol = std::numeric_limits<T>::epsilon(); // zero for integer-type T
  R Rtol = std::numeric_limits<R>::epsilon(); // zero for integer-type R
  bool useTtol = false;
  if ( isfpT || isfpR ){
    if (isfpT && isfpR){
      if (std::is_convertible<R,T>::value) useTtol=true;
      else if (!std::is_convertible<T,R>::value) return false; // they can't be equal in this case
    } else if ( isfpT ) useTtol=true;
  }
  // if both a and b are close to epsilon for its type, our comparison of |a-b| to |a+b| might fail
  bool answer;
  if ( my_abs(a) <= 100*Ttol && my_abs(b) <= 100*Rtol )
    answer = my_abs(a-b) < 100*(useTtol ? Ttol :Rtol);
  else
    answer = my_abs(a-b) < 100*(useTtol ? Ttol :Rtol)*my_abs(a+b);
  return answer;
  // return ( my_abs(a-b) > 100*(useTtol ? Ttol :Rtol)*my_abs(a+b) ) ? false : true;
}
template<typename T, typename R, int N, int M> bool approx_array(const T *A, const R *B){
  bool isfpT, isfpR;
  isfpT = std::is_floating_point<T>::value;
  isfpR = std::is_floating_point<R>::value;
  T Ttol = std::numeric_limits<T>::epsilon(); // zero for integer-type T
  R Rtol = std::numeric_limits<R>::epsilon(); // zero for integer-type R
  bool useTtol = false;
  if ( isfpT || isfpR ){
    if (isfpT && isfpR){
      if (std::is_convertible<R,T>::value) useTtol=true;
      else if (!std::is_convertible<T,R>::value) return false; // they can't be equal in this case
    } else if ( isfpT ) useTtol=true;
  }
  // tol is defined for any combinations of types (might be zero).
  // since tol is epsilon, we need to make sure the value we check is scaled
  // by the sum of the items we're looking at the difference of
  bool answer=true;
  for (int i=0; i<N*M; i++){
    // if both a and b are close to epsilon for its type, our comparison of |a-b| to |a+b| might fail
    if ( my_abs(A[i]) <= 100*Ttol && my_abs(B[i]) <= 100*Rtol )
      answer = my_abs(A[i]-B[i]) < 100*(useTtol ? Ttol :Rtol);
    else
      answer = my_abs(A[i]-B[i]) < 100*(useTtol ? Ttol :Rtol)*my_abs(A[i]+B[i]);
  }
  return answer;
}
template<typename T, typename R, int N> bool approx_matrix(const T *A, const R *B){return approx_array<T,R,N,N>(A,B);}
template<typename T, typename R, int N> bool approx_vector(const T *A, const R *B){return approx_array<T,R,N,1>(A,B);}


// array multiplication C = A * B -- where C is (N,M), A is (N,I) and B is (I,M)
template<typename T, typename R, typename S, int N, int I, int M> void multiply_arrays(T *C, const R *A, const S *B){
  for (int i=0;i<N*M;i++) C[i]=T(0);
  for (int i=0;i<N;i++) for (int j=0;j<M;j++) for (int k=0;k<I;k++) C[i*M+j] += T(A[i*I+k]*B[k*M+j]);
}
template<typename T, typename R, typename S, int N> void multiply_matrix_matrix(T *C, const R *A, const S *B){ multiply_arrays<T,R,S,N,N,N>(C,A,B); }
template<typename T, typename R, typename S, int N> void multiply_matrix_vector(T *C, const R *A, const S *b){ multiply_arrays<T,R,S,N,N,1>(C,A,b); }
template<typename T, typename R, typename S, int N> void multiply_vector_matrix(T *C, const R *a, const S *B){ multiply_arrays<T,R,S,1,N,N>(C,a,B); }

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
  R *cof = safealloc<R>( (N-1)*(N-1) );

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
  R pm=1, *cof = safealloc<R>((N-1)*(N-1));
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
  if (my_abs(d) > tol ){
    // the matrix is *not* singular and has an inverse
    matrix_adjoint(invM,M,N); // at this point invM is the Adjoint(M)
    // The inv(M) = adjoint(M)/det(M)
    int nn = N*N;
    for (int i=0; i<nn; i++) invM[i] /= d;
    // now invM is inv(M)!
  }
  return d;
}
template<typename R> bool matrix_inverse(R *invM, const R *M, R tol, const int N){ return ( my_abs( matrix_determinant_and_inverse(invM, M, tol, N) ) > tol) ; }

template<typename T, int N> bool similar_matrix(T *M, const T *A, const T *B, const T tol){
  T *C = safealloc<T>(N*N);
  bool ok = matrix_inverse(C,B,tol,N);
  if ( ok ){
    T *P = safealloc<T>(N*N);
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
  R *Lt = safealloc<R>(N*N);
  matrix_transpose<R,N>(Lt,L);
  multiply_matrix_matrix<R,R,R,N>(M,Lt,L);
  delete[] Lt;
}

template<typename R, int N> R vector_norm_squared(const R *v){
  R vv = 0;
  for (int i=0; i<N; i++) vv += v[i]*v[i];
  return vv;
}
template<typename T, typename R, typename S, int N> void vector_cross(T *c, const R *a, const S *b){
  if (3==N){
    c[0] = (R)(a[1]*b[2] - a[2]*b[1]);
    c[1] = (R)(a[2]*b[0] - a[0]*b[2]);
    c[2] = (R)(a[0]*b[1] - a[1]*b[0]);
  }
}
template<typename R, int N> R vector_dot(const R *a, const R *b){
  R out = 0;
  for (int i=0; i<N; i++) out += a[i]*b[i];
  return out;
}

template<typename T> T mod1(const T a){
  T b = a - my_cast<int,T>(a);
  return ( (b < T(0) - ZERO_PREC ) ? b + T(1) : b );
}

template<typename T, typename R, int N> bool is_int_matrix(const T * A, const R tol){
  for (int i=0; i<N*N; i++) if ( my_abs(my_cast<int,T>(A[i]) - A[i]) > tol ) return false;
  return true;
}
template<typename R> bool is_int_matrix(const int *, const R){ return true; }


template<typename T> T frobenius_distance(const size_t n, const T* A, const T* B){
  // A-B
  T* AmB = safealloc<T>(n*n);
  for (size_t i=0; i<n*n; ++i) AmB[i] = A[i]-B[i];
  // (A-B)'
  T* AmBt = safealloc<T>(n*n);
  for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) AmBt[i*n+j] = AmB[i+j*n];
  // (A-B)x(A-B)'
  T* mult = safealloc<T>(n*n);
  for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j){
    mult[i*n+j] = T(0);
    for (size_t k=0; k<n; ++k) mult[i*n+j] += AmB[i*n+k]*AmBt[k*n+j];
  }
  delete[] AmB; delete[] AmBt;
  // tr( (A-B)x(A-B)')
  T tr = T(0);
  for (int i=0; i<n; i++) tr += mult[i*(1+n)];
  delete[] mult;
  // sqrt( tr( (A-B)x(A-B)') )
  return std::sqrt(tr);
}

template<typename T> T frobenius_distance(const size_t n, const std::complex<T>* A, const std::complex<T>* B){
  // A-B
  std::complex<T>* AmB = safealloc<std::complex<T>>(n*n);
  for (size_t i=0; i<n*n; ++i) AmB[i] = A[i]-B[i];
  // (A-B)'
  std::complex<T>* AmBt = safealloc<std::complex<T>>(n*n);
  for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) AmBt[i*n+j] = std::conj(AmB[i+j*n]);
  // (A-B)x(A-B)'
  std::complex<T>* mult = safealloc<std::complex<T>>(n*n);
  for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j){
    mult[i*n+j] = std::complex<T>(0);
    for (size_t k=0; k<n; ++k) mult[i*n+j] += AmB[i*n+k]*AmBt[k*n+j];
  }
  delete[] AmB; delete[] AmBt;
  // tr( (A-B)x(A-B)'), which is guaranteed to be real
  T tr{0};
  for (int i=0; i<n; ++i) tr += std::real(mult[i*(1+n)]);
  delete[] mult;
  // sqrt( tr( (A-B)x(A-B)') )
  return std::sqrt(tr);
}


template<typename T> T vector_angle(const size_t n, const T* A, const T* B){
  T AA=0, BB=0, AB=0;
  for (size_t i=0; i<n; ++i){
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
  if (approx_scalar(act, 1.0) && act>1){
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

template<typename T> T vector_angle(const size_t n, const std::complex<T>* A, const std::complex<T>* B){
  // return hermitian_angle(n,A,B);
  return euclidean_angle(n,A,B);
}

template<typename T> T euclidean_angle(const size_t n, const std::complex<T>* A, const std::complex<T>* B){
  T AB{0}, nA{0}, nB{0}, c_t;
  // Compute the products of complex n-vectors as if they were real 2n-vectors
  for (size_t i=0; i<n; ++i){
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
  if (approx_scalar(act, 1.0) && act>1){
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

template<typename T> T hermitian_product(const size_t n, const T* a, const T* b){
  T h_dot{0};
  for (size_t i=0; i<n; ++i) h_dot += a[i]*b[i];
  return h_dot;
}
template<typename T> std::complex<T> hermitian_product(const size_t n, const std::complex<T>* a, const std::complex<T>* b){
  std::complex<T> h_dot{0,0};
  for (size_t i=0; i<n; ++i) h_dot += std::conj(a[i])*b[i];
  return h_dot;
}
template<typename T> T hermitian_angle(const size_t n, const T* A, const T* B){
  return vector_angle(n,A,B);
}
template<typename T> T hermitian_angle(const size_t n, const std::complex<T>* A, const std::complex<T>* B){
  std::complex<T> AA=hermitian_product(n,A,A);
  std::complex<T> BB=hermitian_product(n,B,B);
  std::complex<T> AB=hermitian_product(n,A,B);

  T nAB, nA, nB, c_t;
  nAB = std::sqrt(std::real(AB*std::conj(AB)));
  nA = std::sqrt(std::real(AA));
  nB = std::sqrt(std::real(BB));
  if (nA && nB){
    c_t = nAB/(nA*nB);
  } else {
    // if both are zero their angle is strictly undefined, but for our
    // purposes it would be better to set the angle to 0 [and cos(0)=1]
    c_t = (nA || nB) ? T(0) : T(1);
  }
  T act = std::abs(c_t);
  if (approx_scalar(act,1.0) && act>1){
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


template<typename T> const std::string my_to_string(const T x){
  return (x>0?"+":"") + std::to_string(x);
}
template<typename T> const std::string my_to_string(const std::complex<T> x){
  return (std::real(x)>0?"+":"") + std::to_string(std::real(x)) + (std::imag(x)>0?"+i":"-i") + std::to_string(std::abs(std::imag(x)));
}

template<typename T> T vector_distance(const size_t n, const T* a, const T* b){
  T d, sum=0;
  for (size_t i=0; i<n; ++i){
    d = a[i]-b[i];
    sum += d*d;
  }
  return std::sqrt(sum);
}
template<typename T> T vector_distance(const size_t n, const std::complex<T>* a, const std::complex<T>* b){
  std::complex<T> d;
  T sum=0;
  for (size_t i=0; i<n; ++i){
    d = a[i]-b[i];
    sum += std::real(d*std::conj(d));
  }
  return std::sqrt(sum);
}

template<typename T> T vector_product(const size_t n, const T* a, const T* b){
  T h_dot{0};
  for (size_t i=0; i<n; ++i) h_dot += a[i]*b[i];
  return h_dot;
}
template<typename T> T vector_product(const size_t n, const std::complex<T>* a, const std::complex<T>* b){
  std::complex<T> h_dot = hermitian_product(n,a,b);
  return std::real(h_dot*std::conj(h_dot));
}

template<typename T> T inner_product(const size_t n, const T* a, const T* b){
  T h_dot{0};
  for (size_t i=0; i<n; ++i) h_dot += a[i]*b[i];
  return h_dot;
}
template<typename T> T inner_product(const size_t n, const std::complex<T>* a, const std::complex<T>* b){
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
  size_t onesign, s1mod;
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
  size_t onesign, s1mod;
  // we can permute each sign quaternary number up to three times
  for (size_t i=1; i<4; ++i){
    s1mod = 0;
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


template<class T> T antiphase(const T z){
  return std::signbit(z) ? T(-1) : T(1);
}
template<class T> std::complex<T> antiphase(const std::complex<T> z){
  return std::polar(T(1),-std::arg(z));
}
