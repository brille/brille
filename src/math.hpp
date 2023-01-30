#ifndef BRILLE_MATH_HPP
#define BRILLE_MATH_HPP

#include <cmath>

namespace brille::math{
//const double             pi = 3.14159265358979323846;
//const double         halfpi = 1.57079632679489661923;
// const double      quarterpi = 0.785398163397448309616;
// const double      inversepi = 0.318309886183790671538;
// const double      twooverpi = 0.636619772367581343076;
// const double  twooversqrtpi = 1.12837916709551257390;
// const double        sqrttwo = 1.41421356237309504880;
// const double sqrttwoovertwo = 0.707106781186547524401;
// const double sqrtthreeovertwo = 0.8660254037844386;
// const double              e = 2.71828182845904523536;
// const double          log2e = 1.44269504088896340736;
// const double         log10e = 0.434294481903251827651;
// const double            ln2 = 0.693147180559945309417;
// const double           ln10 = 2.30258509299404568402;

  static const long double long_pi = std::acos(-1.0L);
  static const long double long_two_pi = 2.0L * long_pi;
  static const long double long_one_thirds_pi = long_pi / 3.0L;
  static const long double long_half_pi = std::acos(0.0L);
  static const long double long_two_thirds_pi = long_two_pi / 3.0L;
  static const long double long_half_root_two = std::sqrt(2.0L) / 2.0L;
  static const long double long_half_root_three = std::sqrt(3.0L) / 2.0L;


  static const double pi = static_cast<double>(long_pi);
  static const double two_pi = static_cast<double>(long_two_pi);
  static const double one_thirds_pi = static_cast<double>(long_one_thirds_pi);
  static const double half_pi = static_cast<double>(long_half_pi);
  static const double two_thirds_pi = static_cast<double>(long_two_thirds_pi);
  static const double half_root_two = static_cast<double>(long_half_root_two);
  static const double half_root_three = static_cast<double>(long_half_root_three);

template<class T>
static T degree_to_radian(T degree){
  return static_cast<T>((degree / T(180)) * pi);
}

template<class T>
T cosd(T degree){
  if (!std::isfinite(degree)) return std::cos(degree);
  int quotient;
  T degree90 = std::remquo(std::fabs(degree), T(90), &quotient);
  T radian = degree_to_radian(degree90);
  switch (quotient % 4) {
  case 0:
    return std::cos(radian);
  case 1:
    return std::sin(-radian + T(0)); // add 0 to avoid -0.
  case 2:
    return -std::cos(radian);
  case 3:
    return std::sin(radian + T(0));
  }
  return T(0);
}

template<class T>
T sind(T degree){
  if (!std::isfinite(degree)) return std::sin(degree);
  if (degree < T(0)) return -sind(-degree);
  int quotient;
  T degree90 = std::remquo(std::fabs(degree), T(90), &quotient);
  T radian = degree_to_radian(degree90);
  switch (quotient % 4) {
  case 0:
    return std::sin(radian * T(1));
  case 1:
    return std::cos(radian);
  case 2:
    return std::sin(-radian * T(1));
  case 3:
    return -std::cos(radian);
  }
  return T(0);
}

template<class T> std::tuple<T,T> check_cos_sin_consistency(T c, T s){
//  if (c*c + s*s != T(1)){
//    auto new_c = std::sqrt(T(1) - s*s);
//    if (c < T(0)) new_c *= -1;
//    return std::make_tuple(new_c, s);
//  }
  return std::make_tuple(c, s);
}

template<class T> std::tuple<T, T> cos_and_sin(T r){
  T zero{0}, one{1}, half{T(1)/T(2)};
  if (T(            0) == r) return std::make_tuple(  one,            zero);
  if (T(one_thirds_pi) == r) return std::make_tuple( half, half_root_three);
  if (T(      half_pi) == r) return std::make_tuple( zero,             one);
  if (T(two_thirds_pi) == r) return std::make_tuple(-half, half_root_three);
  if (T(           pi) == r) return std::make_tuple( -one,            zero);
  return check_cos_sin_consistency(std::cos(r), std::sin(r));
}

template<class T> std::tuple<T, T> cos_and_sin_d(T d){
  T zero{0}, one{1}, half{T(1)/T(2)};
  if (T(  0) == d) return std::make_tuple(             one,            zero);
  if (T( 30) == d) return std::make_tuple( half_root_three,            half);
  if (T( 45) == d) return std::make_tuple(   half_root_two,   half_root_two);
  if (T( 60) == d) return std::make_tuple(            half, half_root_three);
  if (T( 90) == d) return std::make_tuple(            zero,             one);
  if (T(120) == d) return std::make_tuple(           -half, half_root_three);
  if (T(135) == d) return std::make_tuple(  -half_root_two,   half_root_two);
  if (T(150) == d) return std::make_tuple(-half_root_three,            half);
  if (T(180) == d) return std::make_tuple(            -one,            zero);
  return check_cos_sin_consistency(cosd(d), sind(d));
}

}

#endif