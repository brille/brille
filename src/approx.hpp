/* Copyright 2019-2020 Greg Tucker
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

#ifndef _APPROX_HPP_
#define _APPROX_HPP_

#include <tuple>
#include <cmath>
#include <type_traits>
#include <limits>
#include <stdexcept>  // for std::overflow_error
#include <string>
#include <math.h>
// #include <complex> // included by debug
#include <numeric>
#include <iostream>
#include "debug.hpp"


namespace brille{
  template<typename T>
  std::enable_if_t<!std::is_unsigned_v<T>, T>
  abs(const T x){ return std::abs(x); }

  template<typename T>
  std::enable_if_t<std::is_unsigned_v<T>, T>
  abs(const T x){ return x; }

  namespace approx{
    // A mulitplier for the approximate-comparison tolerance
    // 10000 is too big for Monoclinic system (5.7224, 5.70957, 4.13651),(90,90.498,90),'C -2y'
    // but setting it any lower (9000 tried) causes other test lattices, namely,
    // (7.189, 4.407, 5.069) (90,90.04,90) '-C 2y' to throw a runtime error
    #if defined(_MSC_VER) || defined(__MINW32__)
    const int TOL_MULT=20000;
    #else
    const int TOL_MULT=10000;
    #endif

    /* isfpT | isfpR | which? | why?
       ------|-------|--------|-----
         0       0     either    both Ttol and Rtol are 0
         1       1      Ttol     R is convertible to T
         0       1      Rtol     Ttol is 0, so use Rtol
         1       0      Ttol     Rtol is 0, so use Ttol
    */

    /*! \brief Returns tuple of tolerance information for approximate comparisons for two datatypes, T and R

    The tuple contains four elements, the first is true if either T or R is an integer or if R can be converted to T.
    The second is true if T is a floating point datatype.
    The third is proportional to epsilon of the datatype T.
    The fourth is proportional to epsilon of the datatype R.
    */
    template<class T, class R>
    std::tuple<bool,bool,T,R> tols(const int tol=1){
      T Ttol = std::numeric_limits<T>::epsilon(); // zero for integer-type T
      R Rtol = std::numeric_limits<R>::epsilon(); // zero for integer-type R
      bool TorRisInteger = Ttol*Rtol==0 || std::is_convertible<T,R>::value;
      bool TisFloatingPt = Ttol > 0;
      Ttol *= static_cast<T>(tol)*static_cast<T>(TOL_MULT);
      Rtol *= static_cast<R>(tol)*static_cast<R>(TOL_MULT);
      return std::make_tuple(TorRisInteger, TisFloatingPt, Ttol, Rtol);
    }



    /* \NOTE For future Greg:
    The following _* templates *CAN NOT BE PREDECLARED* in the hpp file.
    Doing so may not upset the compiler -- after all, it finds a suitable
    definition of _* for every call -- but will anger the linker since the
    actual definitions get overlooked entirely and the symbols never get built.

    Leave these templates alone if you can. If you can't try not to be clever.
    */
    // Catchall case to capture unsigned T and unsigned R
    template<typename T, typename R>
    std::enable_if_t<std::is_integral_v<T> && std::is_integral_v<R>, bool>
    _scalar(const T a, const R b, const bool, const T, const R){
      return a==b;
    }
    // Integral T and Floating Point R
    template<typename T, typename R>
    std::enable_if_t<std::is_integral_v<T> && !std::is_integral_v<R>, bool>
    _scalar(const T a, const R b, const bool, const T, const R Rtol){
      if ( a == T(0) && brille::abs(b) <= Rtol )
        return true;
      else
        return brille::abs(a-b) <= Rtol*brille::abs(a+b);
    }
    // Floating Point T and Integral R
    template<typename T, typename R>
    std::enable_if_t<!std::is_integral_v<T> && std::is_integral_v<R>, bool>
    _scalar(const T a, const R b, const bool, const T Ttol, const R){
      if ( brille::abs(a) <= Ttol && b == R(0) )
        return true;
      else
        return brille::abs(a-b) <= Ttol*brille::abs(a+b);
    }
    // Floating Point T and R
    template<typename T, typename R>
    std::enable_if_t<!std::is_integral_v<T> && !std::is_integral_v<R>, bool>
    _scalar(const T a, const R b, const bool useTtol, const T Ttol, const R Rtol){
      // if both a and b are close to epsilon for its type, our comparison of |a-b| to |a+b| might fail
      if ( brille::abs(a) <= Ttol && brille::abs(b) <= Rtol )
        return brille::abs(a-b) <= (useTtol ? Ttol :Rtol);
      else
        return brille::abs(a-b) <= (useTtol ? Ttol :Rtol)*brille::abs(a+b);
    }

    template<class T, class R>
    bool scalar(const T a, const R b, const int tol=1){
      auto [convertible, useT, Ttol, Rtol] = tols<T,R>(tol);
      return convertible && _scalar(a, b, useT, Ttol, Rtol);
    }

    //! The array comparitor with meta information predetermined
    template<class T, class R>
    bool _array(const size_t NM, const T* a, const R* b, const bool useT, const T Ttol, const R Rtol){
      bool answer=true;
      // we need <= in case T and R are integer, otherwise this is *always* false since 0 !< 0
      if (useT){
        for (size_t i=0; i<NM; ++i){
          // if both a and b are close to epsilon for its type, our comparison of |a-b| to |a+b| might fail
          if ( brille::abs(a[i]) <= Ttol && brille::abs(b[i]) <= Rtol )
          answer &= brille::abs(a[i]-b[i]) <= Ttol;
          else
          answer &= brille::abs(a[i]-b[i]) <= Ttol*brille::abs(a[i]+b[i]);
        }
      } else {
        for (size_t i=0; i<NM; ++i){
          // if both a and b are close to epsilon for its type, our comparison of |a-b| to |a+b| might fail
          if ( brille::abs(a[i]) <= Ttol && brille::abs(b[i]) <= Rtol )
          answer &= brille::abs(a[i]-b[i]) <= Rtol;
          else
          answer &= brille::abs(a[i]-b[i]) <= Rtol*brille::abs(a[i]+b[i]);
        }
      }
      return answer;
    }

    //! Compare N by M arrays for approximate equivalency
    template<class T, class R>
    bool array(const size_t N, const size_t M,const T *a, const R *b, const int tol=1){
      auto [convertible, useT, Ttol, Rtol] = tols<T,R>(tol);
      return convertible && _array(N*M, a, b, useT, Ttol, Rtol);
    }

    //! Compare N by N matrices for approximate equivalency
    template<class T, class R>
    bool matrix(const size_t N, const T *a, const R *b, const int tol=1){
      auto [convertible, useT, Ttol, Rtol] = tols<T,R>(tol);
      return convertible && _array(N*N, a, b, useT, Ttol, Rtol);
    }

    //! Compare length N vectors for approximate equivalency
    template<class T, class R>
    bool vector(const size_t N, const T *a, const R *b, const int tol=1){
      auto [convertible, useT, Ttol, Rtol] = tols<T,R>(tol);
      return convertible && _array(N, a, b, useT, Ttol, Rtol);
    }


    //! Compare N by M arrays when the size is known at compilation-time
    template<class T, class R, size_t N, size_t M>
    bool array(const T *a, const R *b, const int tol=1){
      return brille::approx::array(N,M,a,b,tol);
    }
    //! Compare N by N matrices when the size is known at compilation-time, defaulting to N=3
    template<class T, class R, size_t N=3>
    bool matrix(const T *a, const R *b, const int tol=1){
      return brille::approx::matrix(N,a,b,tol);
    }
    //! Compare length N vectors when the size is known at compilation-time, defaulting to N=3
    template<class T, class R, size_t N=3>
    bool vector(const T *a, const R *b, const int tol=1){
      return brille::approx::vector(N,a,b,tol);
    }

  } // approx::
} // brille::

#endif // _APPROX_HPP_
