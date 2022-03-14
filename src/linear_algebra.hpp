#ifndef _BRILLE_LINEAR_ALGEBRA_HPP_
#define _BRILLE_LINEAR_ALGEBRA_HPP_
#include <vector>
#include <array>
#include "math.hpp"
namespace brille::linear_algebra {
  template<class S, class I1, class I2>
    S dot(I1 a, const I1& end, I2 b){
    auto d = (*a) * (*b);
    while (++a != end) d += (*a) * (*(++b));
    return d;
  }
  template<class T,class R,class S=std::common_type_t<T,R>>
  S dot(const std::vector<T>& a, const std::vector<R>& b){
    assert(a.size() == b.size());
    return dot<S>(a.begin(), a.end(), b.begin());
  }
  template<class T, class R, size_t N, class S=std::common_type_t<T,R>>
  S dot(const std::array<T,N>& a, const std::array<R,N>& b){
    return dot<S>(a.begin(), a.end(), b.begin());
  }
  template<class T> T norm(const std::vector<T>& a){
    return std::sqrt(dot<T>(a.begin(), a.end(), a.begin()));
  }
  template<class T, size_t N> T norm(const std::array<T,N>& a){
    return std::sqrt(dot<T>(a.begin(), a.end(), a.begin()));
  }

  template<class T, size_t N, size_t M> std::array<T,M> mul_mat_vec(const std::array<T,N>& A, const std::array<T,M>& x){
    // multiply a matrix (stored as a row-ordered std::array) and a vector: Matrix * Vector
    assert(N == M * M);
    std::array<T,M> b;
    auto row = A.begin();
    for (size_t i=0; i<M; ++i) {
      b[i] = dot<T>(x.begin(), x.end(), row); // row takes 0--(M-1), then M--(2M-1), then ...
    }
    return b;
  }
  template<class T, size_t N, size_t M> std::array<T,M> mul_vec_mat(const std::array<T,N>& A, const std::array<T,M>& x){
    // multiply a vector and a matrix (stored as a row-ordered std::array): (row)-Vector * Matrix
    assert(N == M * M);
    std::array<T,M> b;
    auto At = transpose(A);
    auto row = At.begin();
    for (size_t i=0; i<M; ++i) {
      b[i] = dot<T>(x.begin(), x.end(), row); // row takes 0--(M-1), then M--(2M-1), then ...
    }
    return b;
  }

  template<class T, class R, size_t N, class S=std::common_type_t<T,R>>
    std::array<S,N> mul_mat_mat(const std::array<T,N> & A, const std::array<R,N>& B){
    // multiply two square matrices (stored as row-ordered std::arrays)
    auto M = static_cast<size_t>(std::sqrt(N));
    assert(N == M*M);
    std::array<S,N> C;
    C.fill(0);
    for (size_t i=0;i<M;i++) for (size_t j=0;j<M;j++) for (size_t k=0;k<M;k++) C[i*M+j] += S(A[i*M+k]*B[k*M+j]);
    return C;
  }

  template<class T, size_t N> std::tuple<bool, T, std::array<T,N>> LU_decomposition(const std::array<T,N> & A){
    auto M = static_cast<size_t>(std::sqrt(N));
    assert(N == M*M);
    std::array<T,N> LU;
    std::copy(A.begin(), A.end(), LU.begin());
    std::vector<size_t> refs(M);
    std::iota(refs.begin(), refs.end(), 0u);
    T det{1};
    auto s2i = [&](const size_t & i, const size_t & j){return refs[i] * N + j;};

    // LU factorisation:
    for (size_t p=0; p < N-1; ++p){
      // locate pivot element
      for (size_t i=p; i < N; ++i){
        if (std::abs(LU[s2i(i, p)]) > std::abs(LU[s2i(p, i)])) {
          // switch the index for the p -1 pivot row
          std::swap(refs[p], refs[i]);
          det = -det;
        }
      }
      if (LU[s2i(p, p)] == 0) {
        return std::make_tuple(false, T(0), LU); // singular matrix
      }
      // multiply diagonal
      det *= LU[s2i(p, p)];
      // form multiplier
      for (size_t i=p+1; i < N; ++i){
        LU[s2i(i, p)] /= LU[s2i(p, p)];
        // eliminate [p-1]
        for (size_t j=p+1; j < N; ++j){
          LU[s2i(i, j)] -= LU[s2i(i, p)] * LU[s2i(p, j)];
        }
      }
    }
    det = det * LU[s2i(N-1, N-1)];
    return std::make_tuple(det != 0, det, LU);
  }

  template<class T, size_t N> std::array<T, N> mat_inverse(const std::array<T,N>& A){
    auto M = static_cast<size_t>(std::sqrt(N));
    assert(N == M * M);
    std::array<T,N> inv; inv.fill(0);
    utils::matrix_inverse(inv.data(), A.data());
    return inv;
  }
}
#endif