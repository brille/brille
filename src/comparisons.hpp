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
/*! \file
    \author Greg Tucker
    \brief Defines a comparisons for `brille`
*/
#include <functional>
#include "approx.hpp"
namespace brille {

/** \brief Binary comparison operators

For runtime-selectable binary comparisons between scalars or arrays these
enumerated values are used in `brille::Comparer`.
*/
enum class cmp {
  lt,    /*!< less than */
  gt,    /*!< greater than */
  le,    /*!< less than or equal */
  ge,    /*!< greater than or equal */
  eq,    /*!< equal */
  nle,   /*!< not less than or equal */
  nge,   /*!< not greather than or equal */
  neq,   /*!< not equal */
  le_ge, /*!< (less than or equal) OR (greater than or equal) */
};

/*! \brief An object to compare scalars or vectors

The comparison-type is set at construction time and two comparison functions
are defined to handle scalar or vector comparisons utilizing approximate
floating point equivalency.
*/
template<class T, class R>
class Comparer{
private:
  bool useT; /*!< from `brille::approx::tols<T,R>` */
  T relT; /*!< from `brille::approx::tols<T,R>` */
  R relR; /*!< from `brille::approx::tols<T,R>` */
  T absT; /*!< from `brille::approx::tols<T,R>` */
  R absR; /*!< from `brille::approx::tols<T,R>` */
  std::function<bool(const T&,const R&)> scalar; /*!< comparison function for scalars */
  std::function<bool(const size_t&,const T*,const size_t&,const R*,const size_t&)> vector; /*!< comparison function for strided vectors */
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
  /*! \brief Compare two scalars

  \param a the first operand for the comparison
  \param b the second operand for the comparison
  \returns A value depending on the comparison for the scalars
  */
  bool operator()(const T a, const R b) const {return this->scalar(a,b);}
  /*! \brief Compare two contiguous vectors

  \param n the number of vector elements
  \param a the first operand to the comparison
  \param b the second operand to the comparison
  \return A value depending on the comparison between the vectors
  */
  bool operator()(const size_t n, const T* a, const R* b) const {
    return this->vector(n, a, 1u, b, 1u);
  }
  /*! \brief Compare two strided vectors

  \param n the number of vector elements
  \param a the first operand to the comparison
  \param sa the stride for `a`
  \param b the second operand to the comparison
  \param sb the stride for `b`
  \return A value depending on the comparison between the vectors
  \note All stride values should be expressed in units of of `sizeof(T)`.
  */
  bool operator()(const size_t n, const T* a, const size_t sa, const R* b, const size_t sb) const {
    return this->vector(n, a, sa, b, sb);
  }
};

/*! \brief Binary operators

For runtime-selectable binary operations between scalars or arrays these
enumerated values are used in `brille::RawBinaryOperator`.
*/
enum class ops {
  plus,  /*!< plus */
  minus, /*!< minus */
  times, /*!< times */
  rdiv,  /*!< right divide */
  ldiv,  /*!< left divide */
};
/*! \brief An object to compute a binary operation on two scalars or vectors

A binary operation function is set at construction-time based on a provided
enumerated binary operation value. The object then provides operator() access
to this function overloaded to operate on scalars, contiguous vectors, and
strided vectors
*/
template<class T>
class RawBinaryOperator{
  ops operation; /*!< The binary operation that this object performs */
  std::function<T(T,T)> binop; /*!< The binary operation function */
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
  /*! \brief Apply the binary operator to two scalars

  \param a the first operand for the binary operation
  \param b the second operand for the binary operation
  \return the result of the binary operation
  */
  T operator()(const T a, const T b) const {return this->binop(a,b);}
  /*! \brief Apply the binary operator to two contigous vectors

  \param n the number of vector elements
  \param[out] c a contigous vector where the binary operation results are stored
  \param a the first operand to the binary operation
  \param b the second operand to the binary operation
  */
  void operator()(const size_t n, T* c, const T* a, const T* b) const {
    for (size_t i=0; i<n; ++i) c[i] = this->binop(a[i],b[i]);
  }
  /*! \brief Apply the binary operator to two strided vectors

  \param n the number of vector elements
  \param[out] c a contigous vector where the binary operation results are stored
  \param sc the stride for `c`
  \param a the first operand to the binary operation
  \param sa the stride for `a`
  \param b the second operand to the binary operation
  \param sb the stride for `b`
  \note All stride values should be expressed in units of of `sizeof(T)`.
  */
  void operator()(const size_t n, T* c, const size_t sc, const T* a, const size_t sa, const T* b, const size_t sb) const {
    for (size_t i=0; i<n; ++i) c[i*sc] = this->binop(a[i*sa],b[i*sb]);
  }
};

/*! \brief Overload `std::to_string` within the `brille` namespace

\param x Any object for which std::to_string is defined
\return std::to_string(x)
*/
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
