/* This file is part of brille.

Copyright © 2019-2022 Greg Tucker <gregory.tucker@ess.eu>

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
/*! \file
    \author Greg Tucker
    \brief Approximate floating point matching methods
*/

#ifndef BRILLE_APPROX_FLOAT_HPP_
#define BRILLE_APPROX_FLOAT_HPP_

#include <tuple>
#include <cmath>
#include <type_traits>
#include <limits>
#include <stdexcept>  // for std::overflow_error
//#include <math.h>
#include <numeric>
#include "debug.hpp"


namespace brille::approx_float{
  // // courtesy of https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
  // template<class T>
  // std::enable_if_t<!std::numeric_limits<T>::is_integer, bool>
  // almost_equal(T x, T y, int ulp){
  // // the machine epsilon has to be scaled to the magnitude of the values used
  // // and multiplied by the desired precision in ULPs (units in the last place)
  // return std::fabs(x-y) <= std::numeric_limits<T>::epsilon() * std::fabs(x+y) * ulp
  //     // unless the result is subnormal
  //     || std::fabs(x-y) < std::numeric_limits<T>::min();
  // }

  #if !defined(DOXYGEN_SHOULD_SKIP_THIS)
  /*! \brief A mulitplier for the approximate-comparison tolerance */
  // 10000 is too big for Monoclinic system (5.7224, 5.70957, 4.13651),(90,90.498,90),'C -2y'
  // but setting it any lower (9000 tried) causes other test lattices, namely,
  // (7.189, 4.407, 5.069) (90,90.04,90) '-C 2y' to throw a runtime error
    #if defined(_MSC_VER) || defined(__MINW32__)
    const int TOL_MULT=20000;
    #else
    const int TOL_MULT=10000;
    #endif
  #endif

  /* isfpT | isfpR | which? | why?
     ------|-------|--------|-----
       0       0     either    both Trel and Rrel are 0
       1       1      Trel     R is convertible to T
       0       1      Rrel     Trel is 0, so use Rrel
       1       0      Trel     Rrel is 0, so use Trel
  */

  /*! \brief Returns tuple of tolerance information for approximate comparisons for two datatypes, `T` and `R`

  \param tol An optional tolerance multiplier defaulting to 1
  \return A tuple containing
  -# `boolean`, true if either `T` or `R` is an integer or if `R` can be converted to `T`
  -# `boolean`, true if `T` is a floating point datatype
  -# `Trel`, proportional to epsilon of the datatype `T`
  -# `Rrel`, proportional to epsilon of the datatype `R`.
  -# `Tabs`, a small value `T(5e-12)`
  -# `Rabs`, a small value `R(5e-12)`
  */
  template<class T, class R>
  std::tuple<bool,bool,T,R,T,R> tols(const T Ttol, const R Rtol, const int tol=1){
    T Trel = std::numeric_limits<T>::epsilon(); // zero for integer-type T
    R Rrel = std::numeric_limits<R>::epsilon(); // zero for integer-type R
    /* 10¹⁵ is fine for a double (but not float ... )     *
     * clang reports that integer 10¹⁵ becomes 999999986991104 if not wrapped
     * with the type-declaration T or R
     * We stick with 5×10⁻¹⁵ to avoid having to chase down where 1×10⁻¹⁵ was
     * too small.
     * */
    T Tabs = T(5)/T(1000000000000000);
    R Rabs = R(5)/R(1000000000000000);
    bool TorRisInteger = Trel*Rrel==0 || std::is_convertible<T,R>::value;
    bool TisFloatingPt = Trel > 0;
    Trel *= static_cast<T>(tol)*static_cast<T>(TOL_MULT);
    Rrel *= static_cast<R>(tol)*static_cast<R>(TOL_MULT);
    if (Ttol > Trel) Trel = Ttol;
    if (Rtol > Rrel) Rrel = Rtol;
    if (Ttol > Tabs) Tabs = Ttol;
    if (Rtol > Rabs) Rabs = Rtol;
//    info_update_if(Ttol > 0 || Rtol > 0, "Ttol=", Ttol, " Rtol=", Rtol," tol=", tol, " gives Trel=", Trel, " Rrel=", Rrel, " Tabs=", Tabs, " Rabs=", Rabs);
    return std::make_tuple(TorRisInteger, TisFloatingPt, Trel, Rrel, Tabs, Rabs);
  }

  #if defined(DOXYGEN_SHOULD_SKIP_THIS)
    /*! \brief Approximate equivalency of two scalars `a ≈ b`

    For relative tolerance \f$r\f$ and absolute tolerance \f$\delta\f$
    determine whether \f$ |a-b| <= \delta + r |a+b| \f$.

    Overloaded to handle comparisons between
    - integer, integer
    - integer, floating point
    - floating point, integer
    - floating point, floating point

    \param a    the first value to compare, datatype `T`
    \param b    the second value to compare, datatype `R`
    \param useT whether `T` tolerances should be used
    \param Trel A relative difference tolerance in datatype `T`
    \param Rrel A relative difference tolerance in datatype `R`
    \param Tabs An absolute difference tolerance in datatype `T`
    \param Rabs An absolute difference tolerance in datatype `R`
    \return boolean

    You must determine `useT`, `Trel`, `Rrel`, `Tabs` and `Rabs` before you
    can use this overloaded function. You may note that those are the exact
    outputs of `brille::approx_float::tols`.
    If you intend to compare many values of the same type, e.g., the elements
    of two arrays, it may be advantageous to use this function directly;
    otherwise you can save yourself a bit of headache by using the wrapper
    `brille::approx_float::scalar`.

    \see scalar, tols
    */
    template<typename T, typename R>
    bool
    _scalar(const T a, const R b, const bool useT, const T Trel, const R Rrel, const T Tabs, const R Rabs);

  #else
    /*************************************************************************
    *                          Future developers:                            *
    **************************************************************************
    * The templates comprising this function *CAN NOT BE PREDECLARED* in the *
    * header file. Doing so may not upset the compiler -- after all, it      *
    * finds a suitable definition for every call -- but will anger the       *
    * linker since the actual definitions get overlooked entirely and the    *
    * symbols never get built.                                               *
    * Leave the templates alone if you can. Otherwise try not to be clever.  *
    *************************************************************************/
    // Catchall case to capture unsigned T and unsigned R
    template<typename T, typename R>
    std::enable_if_t<std::is_integral_v<T> && std::is_integral_v<R>, bool>
    _scalar(const T a, const R b, const bool, const T, const R, const T, const R){
      return a==b;
    }
    // Integral T and Floating Point R
    template<typename T, typename R>
    std::enable_if_t<std::is_integral_v<T> && !std::is_integral_v<R>, bool>
    _scalar(const T a, const R b, const bool, const T, const R Rrel, const T, const R Rabs){
      // if ( a == T(0) && brille::abs(b) <= Rrel ) return true;
      auto x = brille::abs(a-b);
      return x <= Rabs + Rrel*brille::abs(a+b) || x < std::numeric_limits<R>::min();
    }
    // Floating Point T and Integral R
    template<typename T, typename R>
    std::enable_if_t<!std::is_integral_v<T> && std::is_integral_v<R>, bool>
    _scalar(const T a, const R b, const bool, const T Trel, const R, const T Tabs, const R){
      // if ( brille::abs(a) <= Trel && b == R(0) ) return true;
      auto x = brille::abs(a-b);
      return x <= Tabs + Trel*brille::abs(a+b) || x < std::numeric_limits<T>::min();
    }
    // same-type Floating Point T and R
    template<typename T, typename R>
    std::enable_if_t<!std::is_integral_v<T>, bool>
    _scalar(const T a, const T b, const bool, const T Trel, const T, const T Tabs, const T){
      // if both a and b are close to epsilon for its type, our comparison of |a-b| to |a+b| might fail
      auto x = brille::abs(a-b);
      return x <= Tabs + Trel*brille::abs(a+b) || x < std::numeric_limits<T>::min();
    }
    // Floating Point T and R
    template<typename T, typename R>
    std::enable_if_t<!std::is_integral_v<T> && !std::is_integral_v<R>, bool>
    _scalar(const T a, const R b, const bool useT, const T Trel, const R Rrel, const T Tabs, const R Rabs){
      // if both a and b are close to epsilon for its type, our comparison of |a-b| to |a+b| might fail
      auto x = brille::abs(a-b);
      if (useT)
        return x <= Tabs + Trel*brille::abs(a+b) || x < std::numeric_limits<T>::min();
      else
        return x <= Rabs + Rrel*brille::abs(a+b) || x < std::numeric_limits<R>::min();
    }
  #endif

  /*! \brief Gateway function to determine approximate equivalency of two scalars `a ≈ b`

  Determine appropriate absolute and relative tolerances for the datatypes `T`
  and `R`, then compare whether `a` and `b` are approximately the same.

  \param a the first value to compare, datatype `T`
  \param b the second value to compare, datatype `R`
  \param tol an optional tolerance multiplier
  \return boolean

  \see tols, _scalar
  */
  template<class T, class R>
  bool scalar(const T a, const R b, const T Ttol=T(0), const R Rtol=R(0), const int tol=1){
    auto [convertible, useT, Trel, Rrel, Tabs, Rabs] = tols<T,R>(Ttol, Rtol, tol);
    return convertible && _scalar(a, b, useT, Trel, Rrel, Tabs, Rabs);
  }

  /*! \brief Approximate equivalency of two arrays

  \param NM   the number of array elements to compare
  \param a    a pointer to the first value to compare, datatype `T*`
  \param b    a pointer to the second value to compare, datatype `R*`
  \param useT whether `T` tolerances should be used
  \param Trel A relative difference tolerance in datatype `T`
  \param Rrel A relative difference tolerance in datatype `R`
  \param Tabs An absolute difference tolerance in datatype `T`
  \param Rabs An absolute difference tolerance in datatype `R`
  \return boolean
  */
  template<class T, class R>
  bool _array(const size_t NM, const T* a, const R* b, const bool useT, const T Trel, const R Rrel, const T Tabs, const R Rabs){
    bool answer=true;
    // we need <= in case T and R are integer, otherwise this is *always* false since 0 !< 0
    if (useT){
      T Tmin = std::numeric_limits<T>::min();
      for (size_t i=0; i<NM; ++i){
        auto x = brille::abs(a[i]-b[i]);
        // if both a and b are close to epsilon for its type, our comparison of |a-b| to |a+b| might fail
        // if ( brille::abs(a[i]) <= Trel && brille::abs(b[i]) <= Rrel ) answer &= x <= Trel;
        answer &= x <= Tabs + Trel*brille::abs(a[i]+b[i]) || x < Tmin;
      }
    } else {
      R Rmin = std::numeric_limits<R>::min();
      for (size_t i=0; i<NM; ++i){
        auto x = brille::abs(a[i]-b[i]);
        // if both a and b are close to epsilon for its type, our comparison of |a-b| to |a+b| might fail
        // if ( brille::abs(a[i]) <= Trel && brille::abs(b[i]) <= Rrel ) answer &= x <= Rrel;
        answer &= x <= Rabs + Rrel*brille::abs(a[i]+b[i]) || x < Rmin;
      }
    }
    return answer;
  }

  /*! \brief Gateway function to determine approximate equivalency of two arrays

  Determine appropriate absolute and relative tolerances for the datatypes `T`
  and `R`, then compare whether `a` and `b` are approximately the same.

  \param N the first array dimension size
  \param M the second array dimension size
  \param a a pointer to the first value to compare, datatype `T*`
  \param b a pointer to the second value to compare, datatype `R*`
  \param tol an optional tolerance multiplier
  \return boolean

  \see tols, _array
  */
  template<class T, class R>
  bool array(const size_t N, const size_t M,const T *a, const R *b, const T Ttol=T(0), const R Rtol=R(0), const int tol=1){
    auto [convertible, useT, Trel, Rrel, Tabs, Rabs] = tols<T,R>(Ttol, Rtol, tol);
    return convertible && _array(N*M, a, b, useT, Trel, Rrel, Tabs, Rabs);
  }

  /*! \brief Gateway function to determine approximate equivalency of two matrices

  Determine appropriate absolute and relative tolerances for the datatypes `T`
  and `R`, then compare whether `a` and `b` are approximately the same.

  \param N the first and second array dimension size
  \param a a pointer to the first value to compare, datatype `T*`
  \param b a pointer to the second value to compare, datatype `R*`
  \param tol an optional tolerance multiplier
  \return boolean

  \see tols, _array
  */
  template<class T, class R>
  bool matrix(const size_t N, const T *a, const R *b, const T Ttol=T(0), const R Rtol=R(0), const int tol=1){
    auto [convertible, useT, Trel, Rrel, Tabs, Rabs] = tols<T,R>(Ttol, Rtol, tol);
    return convertible && _array(N*N, a, b, useT, Trel, Rrel, Tabs, Rabs);
  }

  /*! \brief Gateway function to determine approximate equivalency of two vectors

  Determine appropriate absolute and relative tolerances for the datatypes `T`
  and `R`, then compare whether `a` and `b` are approximately the same.

  \param N the vector size
  \param a a pointer to the first value to compare, datatype `T*`
  \param b a pointer to the second value to compare, datatype `R*`
  \param tol an optional tolerance multiplier
  \return boolean

  \see tols, _array
  */
  template<class T, class R>
  bool vector(const size_t N, const T *a, const R *b, const T Ttol=T(0), const R Rtol=R(0), const int tol=1){
    auto [convertible, useT, Trel, Rrel, Tabs, Rabs] = tols<T,R>(Ttol, Rtol, tol);
    return convertible && _array(N, a, b, useT, Trel, Rrel, Tabs, Rabs);
  }


  /*! \brief Gateway function to determine approximate equivalency of two arrays

  Determine appropriate absolute and relative tolerances for the datatypes `T`
  and `R`, then compare whether `a` and `b` are approximately the same.

  \param a a pointer to the first value to compare, datatype `T*`
  \param b a pointer to the second value to compare, datatype `R*`
  \param tol an optional tolerance multiplier
  \return boolean

  \see tols, _array, array
  \deprecated As no benefit is provided by knowing the array size at
              compile-time, this function should probably not be used.
  */
  template<class T, class R, size_t N, size_t M>
  bool array(const T *a, const R *b, const T Ttol=T(0), const R Rtol=R(0), const int tol=1){
    return brille::approx_float::array(N,M,a,b,Ttol,Rtol,tol);
  }
  /*! \brief Gateway function to determine approximate equivalency of two matrices

  Determine appropriate absolute and relative tolerances for the datatypes `T`
  and `R`, then compare whether `a` and `b` are approximately the same.

  \param a a pointer to the first value to compare, datatype `T*`
  \param b a pointer to the second value to compare, datatype `R*`
  \param tol an optional tolerance multiplier
  \return boolean

  \see tols, _array, matrix
  \deprecated As no benefit is provided by knowing the matrix size at
              compile-time, this function should probably not be used unless
              if using the default size of \f$ 3 \times 3 \f$.
  */
  template<class T, class R, size_t N=3>
  bool matrix(const T *a, const R *b, const T Ttol=T(0), const R Rtol=R(0), const int tol=1){
    return brille::approx_float::matrix(N,a,b,Ttol,Rtol,tol);
  }
  /*! \brief Gateway function to determine approximate equivalency of two vectors

  Determine appropriate absolute and relative tolerances for the datatypes `T`
  and `R`, then compare whether `a` and `b` are approximately the same.

  \param a a pointer to the first value to compare, datatype `T*`
  \param b a pointer to the second value to compare, datatype `R*`
  \param tol an optional tolerance multiplier
  \return boolean

  \see tols, _array, vector
  \deprecated As no benefit is provided by knowing the vector size at
              compile-time, this function should probably not be used unless
              if using the default size of \f$ 3 \f$.
  */
  template<class T, class R, size_t N=3>
  bool vector(const T *a, const R *b, const T Ttol=T(0), const R Rtol=R(0), const int tol=1){
    return brille::approx_float::vector(N,a,b,Ttol,Rtol,tol);
  }

  template<class T, class R, size_t N>
  bool equal(const std::array<T,N>& a, const std::array<R,N>& b, const T Ttol=T(0), const R Rtol=R(0), const int tol=1){
    return brille::approx_float::vector(N, a.data(), b.data(), Ttol, Rtol, tol);
  }
  template<class T, class R>
  bool equal(const std::vector<T>& a, const std::vector<R>& b, const T Ttol=T(0), const R Rtol=R(0), const int tol=1){
    size_t N{a.size()};
    return b.size() == N && brille::approx_float::vector(N, a.data(), b.data(), Ttol, Rtol, tol);
  }

}

#endif