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

// In order to overload operations between bArrays and LatVecs unambiguously
// these definitions can only be made after all types have been defined:

// 'pure' bArray<T> [+-*/] 'pure' bArray<R>
// (or derived classes with only extra static properties or methods)

#define ARRAY_ARRAY_OP(X) template<class T, class R, template<class> class A, class S = std::common_type_t<T,R>>\
std::enable_if_t<bareArrays<T,A,R,A>, A<S>>\
operator X (const A<T>& a, const A<R>& b){\
  auto itr = a.broadcastItr(b);\
  A<S> out(itr.shape());\
  for (auto [ox, ax, bx]: itr) out[ox] = a[ax] X b[bx];\
  return out;\
}
ARRAY_ARRAY_OP(+)
ARRAY_ARRAY_OP(-)
ARRAY_ARRAY_OP(*)
ARRAY_ARRAY_OP(/)
#undef ARRAY_ARRAY_OP

// LatVec [+-*/] LatVec
// (or derived classes with only extra static properties or methods)
#define ARRAY_ARRAY_OP(X) template<class T, class R, template<class> class A, class S = std::common_type_t<T,R>>\
std::enable_if_t<bothLatVecs<T,A,R,A>, A<S>>\
operator X (const A<T>& a, const A<R>& b){\
  assert(a.samelattice(b));\
  auto itr = a.broadcastItr(b);\
  A<S> out(a.get_lattice(), itr.shape());\
  for (auto [ox, ax, bx]: itr) out[ox] = a[ax] X b[bx];\
  return out;\
}
ARRAY_ARRAY_OP(+)
ARRAY_ARRAY_OP(-)
//ARRAY_ARRAY_OP(*)
//ARRAY_ARRAY_OP(/)
#undef ARRAY_ARRAY_OP

// The old calling syntax (which only partly worked?)
// A<S> out(a, /*mutable*/ true, /*decouple*/ true);
// replaced by
// A<S> out = A<S>(a).decouple();
// below

// any derived bArray<T> class with a copy constructor [+-*/] scalar:
#define ARRAY_SCALAR_OP(X,Y) template<class T, class R, template<class> class A, class S = std::common_type_t<T,R>>\
std::enable_if_t<isArray<T,A>, A<S>>\
operator X (const A<T>& a, const R b){\
  A<S> out = A<S>(a).decouple(); \
  out Y b;\
  return out;\
}
ARRAY_SCALAR_OP(+,+=)
ARRAY_SCALAR_OP(-,-=)
ARRAY_SCALAR_OP(*,*=)
ARRAY_SCALAR_OP(/,/=)
#undef ARRAY_SCALAR_OP

// scalar [+-*] any derived bArray<T> class with a copy constructor:
#define SCALAR_ARRAY_OP(X,Y) template<class T, class R, template<class> class A, class S = std::common_type_t<T,R>>\
std::enable_if_t<isArray<T,A>, A<S>>\
operator X (const R b, const A<T>& a){\
  A<S> out = A<S>(a).decouple();\
  out Y b;\
  return out;\
}
SCALAR_ARRAY_OP(+,+=)
SCALAR_ARRAY_OP(-,-=)
SCALAR_ARRAY_OP(*,*=)
#undef SCALAR_ARRAY_OP

// any derived bArray<T> with a copy constructor [+-*/] std::array<R,N>
// broadcasts the std::array but only if N matches the size of the last dimension of the bArray
#define ARRAY_STDA_OP(X,Y) template<class T, template<class> class A, class R, class S = std::common_type_t<T,R>, size_t Nel>\
std::enable_if_t<isArray<T,A>, A<S>>\
operator X (const A<T>& a, const std::array<R,Nel>& b){\
  assert(a.shape().back() == Nel);\
  A<S> out = A<S>(a).decouple();\
  auto sh = a.shape();\
  sh.back()=0;/*fix the last dimension of the iterator*/\
  for (auto x: out.subItr(sh)) for (size_t i=0; i<Nel; ++i){\
    x.back() = static_cast<brille::ind_t>(i);\
    out[x] Y b[i];\
  }\
  return out;\
}
ARRAY_STDA_OP(+,+=)
ARRAY_STDA_OP(-,-=)
ARRAY_STDA_OP(*,*=)
ARRAY_STDA_OP(/,/=)
#undef ARRAY_STDA_OP


// LatVec {+,-,*,/} 'pure' bArray
#define LATVEC_ARRAY_BINARY_OP(X) \
template<class T, class R, template<class> class L, template<class> class A, class S = std::common_type_t<T,R>>\
std::enable_if_t<(isLatVec<T,L>&&isBareArray<R,A>), L<S>>\
operator X (const L<T>& a, const A<R>& b){\
  auto itr = a.broadcastItr(b);\
  L<S> out(a.get_lattice(), itr.shape());\
  for (auto [ox, ax, bx]: itr) out[ox] = a[ax] X b[bx];\
  return out;\
}
LATVEC_ARRAY_BINARY_OP(+)
LATVEC_ARRAY_BINARY_OP(-)
LATVEC_ARRAY_BINARY_OP(*)
LATVEC_ARRAY_BINARY_OP(/)
#undef LATVEC_ARRAY_BINARY_OP

// bArray {+,-,*,/} LatVec
#define ARRAY_LATVEC_BINARY_OP(X) \
template<class T, class R, template<class> class L, template<class> class A, class S = std::common_type_t<T,R>>\
std::enable_if_t<(isLatVec<T,L>&&isBareArray<R,A>), L<S>>\
operator X (const A<R>& b, const L<T>& a){\
  auto itr = b.broadcastItr(a);\
  L<S> out(a.get_lattice(), itr.shape());\
  for (auto [ox, bx, ax]: itr) out[ox] = b[bx] X a[ax];\
  return out;\
}
ARRAY_LATVEC_BINARY_OP(+)
ARRAY_LATVEC_BINARY_OP(-)
ARRAY_LATVEC_BINARY_OP(*)
ARRAY_LATVEC_BINARY_OP(/)
#undef ARRAY_LATVEC_BINARY_OP

// // cross (LatVec × LatVec)
// template<class T, class R, template<class> class L>
// std::enable_if_t<bothLatVecs<T,L,R,L>, L<double>>
// cross(const L<T>& a, const L<R>& b) {
//   assert( a.samelattice(b) && a.ndim() == b.ndim());
//   assert( a.size(a.ndim()-1) == 3 && b.size(b.ndim()-1)==3 );
//   assert( a.is_row_ordered() && b.is_row_ordered() && a.is_contiguous() && b.is_contiguous() );
//   auto ashape = a.shape(); // (N,M,…,3)
//   auto bshape = b.shape(); // (O,P,…,3) (N≡O unless N≡1 or O≡1, M≡P unless M≡1 or P≡1, etc.)
//   // verify that broadcasting is possible and find the outer array shape
//   auto itr = a.broadcastItr(b);
//   auto oshape = itr.shape();
//   // to use brille::utils::vector_cross we need to iterate over all but the last dimension:
//   ashape.pop_back();
//   bshape.pop_back();
//   auto realitr = L<T>::bItr(ashape, bshape);
//   // perform the cross product storing the result into a temporary array
//   bArray<double> oarray(oshape); // row-ordered contiguous
//   for (auto [ox, ax, bx]: realitr)
//     brille::utils::vector_cross<double,T,R,3>(oarray.ptr(ox), a.ptr(ax), b.ptr(bx));
//   // setup the output array
//   auto lat = a.get_lattice();
//   typename LatVecTraits<L<T>, double>::star cross_star(lat.star(), oarray);
//   cross_star *= lat.get_volume()/2.0/brille::pi;
//   return cross_star.star();
// }
// // cross (Array × Array)
// template<class T, class R, template<class> class L>
// std::enable_if_t<bareArrays<T,L,R,L>, L<double>>
// cross(const L<T>& a, const L<R>& b) {
//   assert( a.ndim() == b.ndim() );
//   assert( a.size(a.ndim()-1) == 3 && b.size(b.ndim()-1)==3 );
//   assert( a.is_row_ordered() && b.is_row_ordered() && a.is_contiguous() && b.is_contiguous() );
//   auto ashape = a.shape(); // (N,M,…,3)
//   auto bshape = b.shape(); // (O,P,…,3) (N≡O unless N≡1 or O≡1, M≡P unless M≡1 or P≡1, etc.)
//   // verify that broadcasting is possible and find the outer array shape
//   auto itr = L<T>::bItr(ashape, bshape);
//   auto oshape = itr.shape();
//   // to use brille::utils::vector_cross we need to iterate over all but the last dimension:
//   ashape.pop_back();
//   bshape.pop_back();
//   auto realitr = L<T>::bItr(ashape, bshape);
//   // perform the cross product storing the result into a temporary array
//   bArray<double> oarray(oshape); // row-ordered contiguous
//   for (auto [ox, ax, bx]: realitr)
//     brille::utils::vector_cross<double,T,R,3>(oarray.ptr(ox), a.ptr(ax), b.ptr(bx));
//   return oarray;
// }
// // dot (Array ⋅ Array)
// // this version requires equal last dimension size
// template<class T, class R, template<class> class A>
// std::enable_if_t<bareArrays<T,A,R,A>, A<double>>
// dot(const A<T>& a, const A<R>& b) {
//   assert( a.ndim() == b.ndim() );
//   assert( a.size(a.ndim()-1) == b.size(b.ndim()-1) );
//   assert( a.is_row_ordered() && b.is_row_ordered() && a.is_contiguous() && b.is_contiguous() );
//   auto ashape = a.shape(); // (N,M,…,3)
//   auto bshape = b.shape(); // (O,P,…,3) (N≡O unless N≡1 or O≡1, M≡P unless M≡1 or P≡1, etc.)
//   auto ndot = ashape.back();
//   ashape.back()=1;
//   bshape.back()=1;
//   // verify that broadcasting is possible and find the outer array shape
//   auto itr = A<T>::bItr(ashape, bshape);
//   auto oshape = itr.shape();
//   // perform the dot product storing the result into a temporary array
//   A<double> oarray(oshape); // row-ordered contiguous
//   for (auto [ox, ax, bx]: itr)
//     oarray[ox] = brille::utils::vector_dot<double,T,R>(ndot, a.ptr(ax), b.ptr(bx));
//   return oarray;
// }
// // // [bArray] dot
// // // this version broadcasts over the last dimension, which seems like a bad idea
// // template<class T, class R, template<class> class A, class S = std::common_type_t<T,R>>
// // std::enable_if_t<bareArrays<T,A,R,A>, A<S>>
// // dot(const A<T> &a, const A<R> &b){
// //   auto itr = A<T>::bItr(a.shape(), b.shape());
// //   auto outshape = itr.shape();
// //   size_t last = outshape.size()-1;
// //   outshape[last] = 1u;
// //   A<S> out(outshape, S(0));
// //   for (auto [ix, ax, bx]: itr){
// //     ix[last] = 0u;
// //     out[last] += a[ax]*b[bx];
// //   }
// //   return out;
// // }
// template<class T, class R, template<class> class A, class S>
// std::enable_if_t<bothArrays<T,A,R,A>, double>
// same_lattice_dot(const A<R>& x, const A<T>& y, const std::vector<S>& len, const std::vector<S>& ang){
//   brille::shape_t a{0,0}, b{0,1}, c{0,2};
//   S x0{static_cast<S>(x[a])}, x1{static_cast<S>(x[b])}, x2{static_cast<S>(x[c])};
//   S y0{static_cast<S>(y[a])}, y1{static_cast<S>(y[b])}, y2{static_cast<S>(y[c])};
//   S out = x0*y0*len[0]*len[0] + x1*y1*len[1]*len[1] + x2*y2*len[2]*len[2]
//         + (x0*y1+x1*y0)*len[0]*len[1]*cos(ang[2])
//         + (x1*y2+x2*y1)*len[1]*len[2]*cos(ang[0])
//         + (x2*y0+x0*y2)*len[2]*len[0]*cos(ang[1]);
//   return out;
// }
// dot [LatVec ⋅ LatVec]
// template<class T, class R, template<class> class L1, template<class> class L2>
// std::enable_if_t<bothLatVecs<T,L1,R,L2>, bArray<double>>
// dot(const L1<T> &a, const L2<R> &b){
//   bool issame = a.samelattice(b);
//   if (!( issame || a.starlattice(b) )){
//     debug_update("Incompatible lattices\n",a.get_lattice().string_repr(),"\n",b.get_lattice().string_repr());
//     throw std::runtime_error("the dot product between Lattice Vectors requires same or starred lattices");
//   }
//   if (!issame){
//     auto out = dot(a.get_hkl(), b.get_hkl());
//     out *= 2*brille::pi; // to avoid making a copy of the array, compared to 2π*dot(...)
//     return out;
//   }
//   assert( a.ndim() == b.ndim());
//   assert( a.size(a.ndim()-1) == 3 && b.size(b.ndim()-1)==3 );
//   assert( a.is_row_ordered() && b.is_row_ordered() && a.is_contiguous() && b.is_contiguous() );
//   auto ashape = a.shape(); // (N,M,…,3)
//   auto bshape = b.shape(); // (O,P,…,3) (N≡O unless N≡1 or O≡1, M≡P unless M≡1 or P≡1, etc.)
//   // to use same_lattice_dot we need to iterate over all but the last dimension:
//   ashape.back()=1;
//   bshape.back()=1;
//   auto itr = L1<T>::bItr(ashape, bshape);
//   auto oshape = itr.shape();
//   // perform the dot product storing the result into the output
//   bArray<double> oarray(oshape); // row-ordered contiguous
//   auto lat = a.get_lattice();
//   std::vector<double> len{lat.get_a(), lat.get_b(), lat.get_c()};
//   std::vector<double> ang{lat.get_alpha(), lat.get_beta(), lat.get_gamma()};
//   for (auto [ox, ax, bx]: itr)
//     oarray[ox] = same_lattice_dot(a.view(ax), b.view(bx), len, ang);
//   return oarray;
// }

/*==============================================================================
  The broadcasting of cross is intimately dependent on the implementation of
  the Array class. If Array2 is used then pop() is not a method of shape_t!
  So this variant must be used instead:
==============================================================================*/
// cross (Array × Array)
template<class T, class R, template<class> class L>
std::enable_if_t<bareArrays<T,L,R,L>, L<double>>
cross(const L<T>& a, const L<R>& b) {
  using namespace brille::utils;
  assert( a.size(1) == 3 && b.size(1)==3 );
  assert( a.is_row_ordered() && b.is_row_ordered() && a.is_contiguous() && b.is_contiguous() );
  brille::ind_t aN=a.size(0), bN=b.size(0);
  assert( 1u==aN || 1u==bN || aN==bN );
  brille::ind_t oO = (1u == aN) ? bN : aN;
  bArray<double> oarray(oO, 3u); // row-ordered contiguous
  if (1u == aN || 1u == bN){
    if (1u == aN){
      for (brille::ind_t j=0; j<bN; ++j) vector_cross<double,T,R,3>(oarray.ptr(j), a.ptr(0), b.ptr(j));
    } else {
      for (brille::ind_t j=0; j<aN; ++j) vector_cross<double,T,R,3>(oarray.ptr(j), a.ptr(j), b.ptr(0));
    }
  } else {
    for (brille::ind_t j=0; j<aN; ++j) vector_cross<double,T,R,3>(oarray.ptr(j), a.ptr(j), b.ptr(j));
  }
  return oarray;
}
// cross (LatVec × LatVec)
template<class T, class R, template<class> class L>
std::enable_if_t<bothLatVecs<T,L,R,L>, L<double>>
cross(const L<T>& a, const L<R>& b) {
  assert( a.samelattice(b) );
  auto oarray = cross(a.get_hkl(), b.get_hkl());
  // setup the output array
  auto lat = a.get_lattice();
  typename LatVecTraits<L<T>, double>::star cross_star(lat.star(), oarray);
  cross_star *= lat.get_volume()/2.0/brille::pi;
  return cross_star.star();
}


// dot (Array2 ⋅ Array2)
template<class T, class R, template<class> class A>
std::enable_if_t<bareArrays<T,A,R,A>, A<double>>
dot(const A<T>& a, const A<R>& b) {
  using namespace brille::utils;
  assert( a.size(1) == b.size(1) );
  assert( a.is_row_ordered() && b.is_row_ordered() && a.is_contiguous() && b.is_contiguous() );
  brille::ind_t aN=a.size(0), bN=b.size(0), d=a.size(1);
  assert( 1u==aN || 1u==bN || aN==bN );
  brille::ind_t oO = (1u == aN) ? bN : aN;
  bArray<double> oarray(oO, 1u);
  if (1u==aN || 1u==bN) {
    if (1u==aN){
      for (brille::ind_t i=0; i<bN; ++i) oarray.val(i,0) = vector_dot<double,T,R>(d, a.ptr(0), b.ptr(i));
    } else {
      for (brille::ind_t i=0; i<aN; ++i) oarray.val(i,0) = vector_dot<double,T,R>(d, a.ptr(i), b.ptr(0));
    }
  } else {
    for (brille::ind_t i=0; i<aN; ++i) oarray.val(i,0) = vector_dot<double,T,R>(d, a.ptr(i), b.ptr(i));
  }
  return oarray;
}

template<class T, class R, template<class> class A, class S>
std::enable_if_t<bothArrays<T,A,R,A>, double>
same_lattice_dot(const A<R>& x, const A<T>& y, const std::vector<S>& len, const std::vector<S>& ang){
  S x0{static_cast<S>(x.val(0,0))}, x1{static_cast<S>(x.val(0,1))}, x2{static_cast<S>(x.val(0,2))};
  S y0{static_cast<S>(y.val(0,0))}, y1{static_cast<S>(y.val(0,1))}, y2{static_cast<S>(y.val(0,2))};
  S out = x0*y0*len[0]*len[0] + x1*y1*len[1]*len[1] + x2*y2*len[2]*len[2]
        + (x0*y1+x1*y0)*len[0]*len[1]*cos(ang[2])
        + (x1*y2+x2*y1)*len[1]*len[2]*cos(ang[0])
        + (x2*y0+x0*y2)*len[2]*len[0]*cos(ang[1]);
  return out;
}

// dot [LatVec ⋅ LatVec]
template<class T, class R, template<class> class L1, template<class> class L2>
std::enable_if_t<bothLatVecs<T,L1,R,L2>, bArray<double>>
dot(const L1<T> &a, const L2<R> &b){
  bool issame = a.samelattice(b);
  if (!( issame || a.starlattice(b) )){
    debug_update("Incompatible lattices\n",a.get_lattice().string_repr(),"\n",b.get_lattice().string_repr());
    throw std::runtime_error("the dot product between Lattice Vectors requires same or starred lattices");
  }
  if (!issame){
    auto out = dot(a.get_hkl(), b.get_hkl());
    out *= 2*brille::pi; // to avoid making a copy of the array, compared to 2π*dot(...)
    return out;
  }
  assert( a.size(1) == 3 && b.size(1)==3 );
  assert( a.is_row_ordered() && b.is_row_ordered() && a.is_contiguous() && b.is_contiguous() );
  brille::ind_t aN=a.size(0), bN=b.size(0);
  assert( 1u==aN || 1u==bN || aN==bN );
  brille::ind_t oO = (1u == aN) ? bN : aN;
  bArray<double> oarray(oO, 1u);
  auto lat = a.get_lattice();
  std::vector<double> len{lat.get_a(), lat.get_b(), lat.get_c()};
  std::vector<double> ang{lat.get_alpha(), lat.get_beta(), lat.get_gamma()};
  if (1u==aN || 1u==bN){
    if (1u==aN) {
      for (brille::ind_t i=0; i<bN; ++i) oarray.val(i,0) = same_lattice_dot(a, b.view(i), len, ang);
    } else {
      for (brille::ind_t i=0; i<aN; ++i) oarray.val(i,0) = same_lattice_dot(a.view(i), b, len, ang);
    }
  } else {
    for (brille::ind_t i=0; i<aN; ++i) oarray.val(i,0) = same_lattice_dot(a.view(i), b.view(i), len, ang);
  }
  return oarray;
}

// [bArray] norm
template<class T, template<class> class L>
std::enable_if_t<isBareArray<T,L>, L<double>>
norm(const L<T> &a){
  L<double> out = dot(a,a);
  for (auto& x : out.valItr()) x = std::sqrt(x);
  return out;
}
// [LatVec] norm
template<class T, template<class> class L>
std::enable_if_t<isLatVec<T,L>, bArray<double>>
norm(const L<T> &a){
  bArray<double> out = dot(a,a);
  for (auto& x: out.valItr()) x = std::sqrt(x);
  return out;
}

// [LatVec],[Array] cat
template<class T, class R, template<class> class LV, template<class> class BA, class S = std::common_type_t<T,R>>
std::enable_if_t<(isLatVec<T,LV>&&isBareArray<R,BA>), BA<S>>
cat(const brille::ind_t dim, const LV<T>& a, const BA<R>& b){
  // BA<S> to ensure type conversion happens *before* the append
  return BA<S>(a.get_hkl()).append(dim, b);
}
// [Array],[LatVec] cat
template<class T, class R, template<class> class LV, template<class> class BA, class S = std::common_type_t<T,R>>
std::enable_if_t<(isLatVec<T,LV>&&isBareArray<R,BA>), BA<S>>
cat(const brille::ind_t dim, const BA<R>& a, const LV<T>& b){
  return BA<S>(a.decouple()).append(dim, b.get_hkl());
}
// [bArray],[bArray] cat
template<class T, class R, template<class> class L, class S = std::common_type_t<T,R>>
std::enable_if_t<bareArrays<T,L,R,L>, L<S>>
cat(const brille::ind_t dim, const L<T>& a, const L<R>& b){
  return L<S>(a.decouple()).append(dim, b);
}
// [LatVec] cat
template<class T, class R, template<class> class L, class S = std::common_type_t<T,R>>
std::enable_if_t<bothLatVecs<T,L,R,L>, L<S>>
cat(const brille::ind_t dim, const L<T>& a, const L<R>& b){
  if (!a.samelattice(b))
    throw std::runtime_error("LatVec cat requires vectors in the same lattice.");
  // use cat([Array],[Array]) to handle possible type-conversion
  return L<S>(a.get_lattice(), cat(dim, a.get_hkl(), b.get_hkl()));
}
// poor-man's variadic Array concatenation
template<class T, class R, template<class> class A, class S=std::common_type_t<T,R>, class... Args>
std::enable_if_t<bothArrays<T,A,R,A>, A<S>>
cat(const brille::ind_t dim, const A<T>& a, const A<R>& b, Args... args){
  return cat(dim,cat(dim,a,b),args...);
}
// star
template<class T, template<class> class L>
std::enable_if_t<isLatVec<T,L>, LatVecStar<L<T>,double> >
star(const L<T>& v){
  std::vector<double> cvmt(9);
  v.get_lattice().get_covariant_metric_tensor(cvmt.data());
  LatVecStar<L<T>,double> vstar( v.get_lattice().star(), v.shape() );
  auto vsh = v.shape();
  vsh.back() = 0;
  for (auto x: v.subItr(vsh)) // iterate over all but the last dimension
    brille::utils::multiply_matrix_vector(vstar.ptr(x), cvmt.data(), v.ptr(x));
  return vstar;
}


// Matrix * LatVec
template<class T, class R, template<class> class L, class S = std::common_type_t<T,R>>
std::enable_if_t<(isLatVec<T,L> && std::is_floating_point<S>::value), L<S>>
operator*(const std::array<R,9>& m, const L<T>& a){
  auto ashape = a.shape();
  assert(a.is_row_ordered() && a.is_contiguous() && ashape.back() == 3);
  L<S> out(a.get_lattice(), ashape);
  ashape.back() = 0;
  for (auto x : a.subItr(ashape))
    brille::utils::multiply_matrix_vector<S,R,T,3>(out.ptr(x), m.data(), a.ptr(x));
  return out;
}


// -LatVec
template<class T, template<class> class L>
std::enable_if_t<isLatVec<T,L>, L<T>>
operator-(const L<T>& a){
  L<T> out(a.get_lattice(), -(a.get_hkl()));
  return out;
}
// -('pure' brille:Array)
template<class T, template<class> class L>
std::enable_if_t<(isBareArray<T,L>&&!std::is_unsigned_v<T>), L<T>>
operator-(const L<T>& a){
  return -1*a;
}

// sum(bArray)
template<class T, template<class> class A>
std::enable_if_t<isBareArray<T,A>, A<T>>
sum(const A<T>& a, typename A<T>::ind_t dim){
  return a.sum(dim);
}
// sum(LatVec)
template<class T, template<class> class A>
std::enable_if_t<isLatVec<T,A>, A<T>>
sum(const A<T>& a, typename A<T>::ind_t dim){
  return A<T>(a.get_lattice(), a.sum(dim));
}


// abs(bArray)
template<class T, template<class> class L>
std::enable_if_t<isBareArray<T,L>, L<T>>
abs(const L<T>& a){
  L<T> out(a.shape());
  for (auto x : a.subItr()) out[x] = brille::utils::magnitude(a[x]);
  return out;
}
// abs(LatVec)
template<class T, template<class> class L>
std::enable_if_t<isLatVec<T,L>, L<T>>
abs(const L<T>& a){
  L<T> out(a.get_lattice(), a.shape());
  for (auto x : a.subItr()) out[x] = brille::utils::magnitude(a[x]);
  return out;
}
