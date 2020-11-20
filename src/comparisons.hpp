/* This file is part of brille.

Copyright © 2020 Greg Tucker <greg.tucker@stfc.ac.uk>

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
#ifndef BRILLE_COMPARISONS_HPP_
#define BRILLE_COMPARISONS_HPP_
#include <functional>
#include "approx.hpp"
namespace brille {

enum class cmp {
  lt,    //< less than
  gt,    //< greater than
  le,    //< less than or equal
  ge,    //< greater than or equal
  eq,    //< equal
  nle,   //< not less than or equal
  nge,   //< not greather than or equal
  neq,   //< not equal
  le_ge, //< (less than or equal) OR (greater than or equal)
};

template<class T, class R>
class Comparer{
private:
  bool useT;
  T relT;
  R relR;
  T absT;
  R absR;
  std::function<bool(const T&,const R&)> scalar;
  std::function<bool(const size_t&,const T*,const size_t&,const R*,const size_t&)> vector;
public:
  Comparer(const cmp op){
    // predetermine tolerances and which we should use:
    bool c;
    std::tie(c, this->useT, this->relT, this->relR, this->absT, this->absR) = brille::approx::tols<T,R>();
    // set the comparison function
    switch(op){
      case cmp::lt:
      scalar = [&](const T& a, const R& b){
        return !brille::approx::_scalar(a,b,useT,relT,relR,absT,absR) && a<b;
      };
      break;
      case cmp::gt:
      scalar = [&](const T& a, const R& b){
        return !brille::approx::_scalar(a,b,useT,relT,relR,absT,absR) && a>b;
      };
      break;
      case cmp::le:
      scalar = [&](const T& a, const R& b){
        return brille::approx::_scalar(a,b,useT,relT,relR,absT,absR) || a<b;
      };
      break;
      case cmp::ge:
      scalar = [&](const T& a, const R& b){
        return brille::approx::_scalar(a,b,useT,relT,relR,absT,absR) || a>b;
      };
      break;
      case cmp::eq:
      scalar = [&](const T& a, const R& b){
        return brille::approx::_scalar(a,b,useT,relT,relR,absT,absR);
      };
      break;
      case cmp::nle:
      scalar = [&](const T& a, const R& b){
        return !brille::approx::_scalar(a,b,useT,relT,relR,absT,absR) && a>b;
      };
      break;
      case cmp::nge:
      scalar = [&](const T& a, const R& b){
        return !brille::approx::_scalar(a,b,useT,relT,relR,absT,absR) && a<b;
      };
      break;
      case cmp::neq:
      scalar = [&](const T& a, const R& b){
        return !brille::approx::_scalar(a,b,useT,relT,relR,absT,absR);
      };
      break;
      default:
      throw std::runtime_error("Unhandled comparison");
    }
    switch (op){
      // A vector is NOT X if any of the elements are NOT X
      case cmp::neq:
      case cmp::nge:
      case cmp::nle:
      vector = [&](const size_t n, const T* a, const size_t sa, const R* b, const size_t sb){
        bool ret{false};
        for (size_t i=0; i<n; ++i) ret |= this->scalar(a[i*sa], b[i*sb]);
        return ret;
      };
      break;
      // A vector IS X if all of the elements ARE X
      default:
      vector = [&](const size_t n, const T* a, const size_t sa, const R* b, const size_t sb){
        bool ret{true};
        for (size_t i=0; i<n; ++i) ret &= this->scalar(a[i*sa], b[i*sb]);
        return ret;
      };
    }
  }
  // call the comparison function for two scalars
  bool operator()(const T a, const R b) const {return this->scalar(a,b);}
  // call the comparison function for two raw vectors
  bool operator()(const size_t n, const T* a, const R* b) const {
    return this->vector(n, a, 1u, b, 1u);
  }
  // call the comparison function for two strided raw vectors
  bool operator()(const size_t n, const T* a, const size_t sa, const R* b, const size_t sb) const {
    return this->vector(n, a, sa, b, sb);
  }
};


enum class ops {
  plus,  //< plus
  minus, //< minus
  times, //< times
  rdiv,  //< right divide
  ldiv,  //< left divide
};
template<class T>
class RawBinaryOperator{
  ops operation;
  std::function<T(T,T)> binop;
public:
  RawBinaryOperator(ops op): operation(op) {
    switch (operation){
      case ops::plus:  binop = [](const T& a, const T& b){return a+b;}; break;
      case ops::minus: binop = [](const T& a, const T& b){return a-b;}; break;
      case ops::times: binop = [](const T& a, const T& b){return a*b;}; break;
      case ops::rdiv:  binop = [](const T& a, const T& b){return a/b;}; break;
      case ops::ldiv:  binop = [](const T& a, const T& b){return b/a;}; break;
      default: throw std::runtime_error("Unhandled operation!");
    }
  }
  T operator()(const T a, const T b) const {return this->binop(a,b);}
  void operator()(const size_t n, T* c, const T* a, const T* b) const {
    for (size_t i=0; i<n; ++i) c[i] = this->binop(a[i],b[i]);
  }
  void operator()(const size_t n, T* c, const size_t sc, const T* a, const size_t sa, const T* b, const size_t sb) const {
    for (size_t i=0; i<n; ++i) c[i*sc] = this->binop(a[i*sa],b[i*sb]);
  }
};


template<class T>
std::string to_string(const T& x){
  return std::to_string(x);
}
template<>
std::string to_string<cmp>(const cmp& c);
template<>
std::string to_string<ops>(const ops& o);

} // namespace brille



#endif
