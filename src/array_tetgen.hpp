#ifndef BRILLE_ARRAY_TETGET_HPP_
#define BRILLE_ARRAY_TETGET_HPP_

//#include <Eigen/Core>
//#include <Eigen/Cholesky>

#include "tetgen.h"
#include "array_.hpp"

namespace brille {

  template<class T, template<class> class A>
  inline std::enable_if_t<isBareArray<T,A>, std::array<T,4>>
  four_point_sphere(const A<T>& a, const A<T>& b, const A<T>& c, const A<T>& d){
    std::array<T,4> z;
    tetgenmesh tgm;
    tgm.circumsphere(a.ptr(0), b.ptr(0), c.ptr(0), d.ptr(0), z.data(), z.data()+3);
    return z;
  }

  template<class T, template<class> class A>
  inline std::enable_if_t<isLatVec<T,A>, std::array<T,4>>
  four_point_sphere(const A<T>& a, const A<T>& b, const A<T>& c, const A<T>& d){
    return four_point_sphere(a.xyz(), b.xyz(), c.xyz(), d.xyz());
  }
  template<class T, class I, template<class> class A, template<class> class B>
  inline std::enable_if_t<(isArray<T,A> && isBareArray<I,B> && std::is_integral_v<I>), std::array<T,4>>
  four_point_sphere(const A<T>& v, const B<I>& i){
    return four_point_sphere(v.view(i[0]), v.view(i[1]), v.view(i[2]), v.view(i[3]));
  }

}
#endif