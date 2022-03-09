#ifndef _BRILLE_LINEAR_ALGEBRA_HPP_
#define _BRILLE_LINEAR_ALGEBRA_HPP_
#include <vector>
#include <array>
#include "math.hpp"
namespace brille::linear_algebra {
  template<class Itr> typename Itr::value_type dot(Itr a, const Itr& end, Itr b){
    auto d = (*a) * (*b);
    while (++a != end) d *= (*a) * (*(++b));
    return d;
  }
  template<class T> T dot(const std::vector<T>& a, const std::vector<T>& b){
    assert(a.size() == b.size());
    return dot(a.begin(), a.end(), b.begin());
  }
  template<class T, size_t N> T dot(const std::array<T,N>& a, const std::array<T,N>& b){
    return dot(a.begin(), a.end(), b.begin());
  }
  template<class T> T norm(const std::vector<T>& a){
    return std::sqrt(dot(a.begin(), a.end(), a.begin()));
  }
  template<class T, size_t N> T norm(const std::array<T,N>& a){
    return std::sqrt(dot(a.begin(), a.end(), a.begin()));
  }

}
#endif