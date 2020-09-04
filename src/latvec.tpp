/* Copyright 2019 Greg Tucker
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

// [brille::Array] dot
template<class T, class R, template<class> class L1, template<class> class L2,
         typename=typename std::enable_if<std::is_base_of<brille::Array<T>,L1<T>>::value&&!std::is_base_of<LatVec,L1<T>>::value>::type,
         typename=typename std::enable_if<std::is_base_of<brille::Array<R>,L2<R>>::value&&!std::is_base_of<LatVec,L2<R>>::value>::type,
         class S = typename std::common_type<T,R>::type>
brille::Array<S> dot(const L1<T> &a, const L2<R> &b){
  auto itr = BroadcastIt(a.shape(), b.shape());
  auto outshape = itr.shape();
  size_t last = outshape.size()-1;
  outshape[last] = 1u;
  brille::Array<S> out(outshape, S(0));
  for (auto [ix, ax, bx]: itr){
    ix[last] = 0u;
    out[last] += a[ax]*b[bx];
  }
  return out;
}
// [brille::Array] norm
template<class T, template<class> class L,
         typename=typename std::enable_if<std::is_base_of<brille::Array<T>,L<T>>::value>::type
        >
brille::Array<double> norm(const L<T> &a){
  brille::Array<double> out = dot(a,a);
  for (auto x : brille::ArrayIt(out)) out[x] = std::sqrt(out[x]);
  return out;
}

// cross (LatVec × LatVec)
template<class T, class R, template<class> class L>
std::enable_if_t<(std::is_base_of_v<LatVec,L<T>>&&std::is_base_of_v<LatVec,L<R>>), L<double> >
cross(const L<T>& a, const L<R>& b) {
  assert( a.samelattice(b) && a.ndim() == b.ndim());
  assert( a.size(a.ndim()-1) == 3 && b.size(b.ndim()-1)==3 );
  assert( a.is_row_ordered() && b.is_row_ordered() && a.is_contiguous() && b.is_contiguous() );
  auto ashape = a.shape(); // (N,M,…,3)
  auto bshape = b.shape(); // (O,P,…,3) (N≡O unless N≡1 or O≡1, M≡P unless M≡1 or P≡1, etc.)
  // verify that broadcasting is possible and find the outer array shape
  auto itr = BroadcastIt(ashape, bshape);
  auto oshape = itr.shape();
  // to use brille::utils::vector_cross we need to iterate over all but the last dimension:
  ashape.pop_back();
  bshape.pop_back();
  auto realitr = BroadcastIt(ashape, bshape);
  // perform the cross product storing the result into a temporary array
  brille::Array<double> oarray(oshape); // row-ordered contiguous
  for (auto [ox, ax, bx]: realitr)
    brille::utils::vector_cross<double,T,R,3>(oarray.ptr(ox), a.ptr(ax), b.ptr(bx));
  // setup the output array
  auto lat = a.get_lattice();
  typename LatVecTraits<L<T>, double>::star cross_star(lat.star(), oarray);
  cross_star *= lat.get_volume()/2.0/brille::pi;
  return cross_star.star();
}

// cross (Array × Array)
template<class T, class R, template<class> class L>
std::enable_if_t<!(std::is_base_of_v<LatVec,L<T>>||std::is_base_of_v<LatVec,L<R>>), L<double> >
cross(const L<T>& a, const L<R>& b) {
  assert( a.ndim() == b.ndim() );
  assert( a.size(a.ndim()-1) == 3 && b.size(b.ndim()-1)==3 );
  assert( a.is_row_ordered() && b.is_row_ordered() && a.is_contiguous() && b.is_contiguous() );
  auto ashape = a.shape(); // (N,M,…,3)
  auto bshape = b.shape(); // (O,P,…,3) (N≡O unless N≡1 or O≡1, M≡P unless M≡1 or P≡1, etc.)
  // verify that broadcasting is possible and find the outer array shape
  auto itr = BroadcastIt(ashape, bshape);
  auto oshape = itr.shape();
  // to use brille::utils::vector_cross we need to iterate over all but the last dimension:
  ashape.pop_back();
  bshape.pop_back();
  auto realitr = BroadcastIt(ashape, bshape);
  // perform the cross product storing the result into a temporary array
  brille::Array<double> oarray(oshape); // row-ordered contiguous
  for (auto [ox, ax, bx]: realitr)
    brille::utils::vector_cross<double,T,R,3>(oarray.ptr(ox), a.ptr(ax), b.ptr(bx));
  return oarray;
}


template<typename T,typename R> double same_lattice_dot(const R* x, const T* y, const double* len, const double* ang){
  double out = double(x[0])*double(y[0])*len[0]*len[0]
             + double(x[1])*double(y[1])*len[1]*len[1]
             + double(x[2])*double(y[2])*len[2]*len[2]
             + (double(x[0])*double(y[1])+double(x[1])*double(y[0]))*len[0]*len[1]*cos(ang[2])
             + (double(x[1])*double(y[2])+double(x[2])*double(y[1]))*len[1]*len[2]*cos(ang[0])
             + (double(x[2])*double(y[0])+double(x[0])*double(y[2]))*len[2]*len[0]*cos(ang[1]);
  return out;
}

// dot [LatVec ⋅ LatVec]
template<class T, class R, template<class> class L1, template<class> class L2>
std::enable_if_t<(std::is_base_of_v<LatVec,L1<T>>&&std::is_base_of_v<LatVec,L2<R>>), brille::Array<double> >
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
  assert( a.ndim() == b.ndim());
  assert( a.size(a.ndim()-1) == 3 && b.size(b.ndim()-1)==3 );
  assert( a.is_row_ordered() && b.is_row_ordered() && a.is_contiguous() && b.is_contiguous() );
  auto ashape = a.shape(); // (N,M,…,3)
  auto bshape = b.shape(); // (O,P,…,3) (N≡O unless N≡1 or O≡1, M≡P unless M≡1 or P≡1, etc.)
  // to use same_lattice_dot we need to iterate over all but the last dimension:
  ashape.back()=1;
  bshape.back()=1;
  auto itr = BroadcastIt(ashape, bshape);
  auto oshape = itr.shape();
  // perform the dot product storing the result into the output
  brille::Array<double> oarray(oshape); // row-ordered contiguous
  auto lat = a.get_lattice();
  std::vector<double> len{lat.get_a(), lat.get_b(), lat.get_c()};
  std::vector<double> ang{lat.get_alpha(), lat.get_beta(), lat.get_gamma()};
  for (auto [ox, ax, bx]: itr)
    oarray[ox] = same_lattice_dot(a.ptr(ax), b.ptr(bx), len, ang);
  return oarray;
}

// dot (Array ⋅ Array)
template<class T, class R, template<class> class L>
std::enable_if_t<!(std::is_base_of_v<LatVec,L<T>>||std::is_base_of_v<LatVec,L<R>>), L<double> >
dot(const L<T>& a, const L<R>& b) {
  assert( a.ndim() == b.ndim() );
  assert( a.size(a.ndim()-1) == b.size(b.ndim()-1) );
  assert( a.is_row_ordered() && b.is_row_ordered() && a.is_contiguous() && b.is_contiguous() );
  auto ashape = a.shape(); // (N,M,…,3)
  auto bshape = b.shape(); // (O,P,…,3) (N≡O unless N≡1 or O≡1, M≡P unless M≡1 or P≡1, etc.)
  auto ndot = ashape.back();
  ashape.back()=1;
  bshape.back()=1;
  // verify that broadcasting is possible and find the outer array shape
  auto itr = BroadcastIt(ashape, bshape);
  auto oshape = itr.shape();
  // perform the dot product storing the result into a temporary array
  brille::Array<double> oarray(oshape); // row-ordered contiguous
  for (auto [ox, ax, bx]: itr)
    oarray[ox] = brille::utils::vector_dot<double,T,R>(ndot, a.ptr(ax), b.ptr(bx));
  return oarray;
}

// // [LatVec] norm
// template<class T, template<class> class L,
//          typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type
//         >
// brille::Array<double> norm(const L<T> &a){
//   brille::Array<double> out = dot(a,a);
//   for (size_t i=0; i<out.size(); ++i) out.insert(sqrt(out.getvalue(i)),i);
//   return out;
// }

// [LatVec] cat
template<class T, class R, template<class> class L,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type,
         class S = typename std::common_type<T,R>::type>
L<S> cat(const L<T>& a, const L<R>& b){
  if (!a.samelattice(b))
    throw std::runtime_error("LatVec cat requires vectors in the same lattice.");
  brille::Array<S> out(a.get_rlu());
  out.append(0, b.get_rlu()); // cat always along the first dimension
  return L<S>(a.get_lattice(), out);
}

// star
template<class T, template <class R> class L,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type>
typename LatVecTraits<L<T>,double>::star star(const L<T>& v){
  std::vector<double> cvmt(9);
  v.get_lattice().get_covariant_metric_tensor(cvmt.data());
  typename LatVecTraits<L<T>,double>::star vstar( v.get_lattice().star(), v.shape() );
  auto vsh = v.shape();
  vsh.back() = 0;
  for (auto x: SubIt(v.shape(), vsh)) // iterate over all but the last dimension
    brille::utils::multiply_matrix_vector(vstar.ptr(x), cvmt.data(), v.ptr(x));
  return vstar;
}


// Matrix * LatVec
template<class T, class R, template<class> class L,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type,
         class S = typename std::common_type<T,R>::type,
         typename=typename std::enable_if<std::is_floating_point<S>::value>::type>
L<S> operator*(const std::array<R,9>& m, const L<T>& a){
  auto ashape = a.shape();
  assert(a.is_row_ordered() && a.is_contiguous() && ashape.back() == 3);
  L<S> out(a.get_lattice(), ashape);
  ashape.back() = 0;
  for (auto x : SubIt(a.shape(), ashape))
    brille::utils::multiply_matrix_vector<S,R,T,3>(out.ptr(x), m.data(), a.ptr(x));
  return out;
}


// -LatVec
template<class T, template<class> class L, typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type>
L<T> operator-(const L<T>& a){
  L<T> out(a.get_lattice(), -(a.get_hkl()));
  return out;
}


// sum(brille::Array)
template<class T>
brille::Array<T> sum(const brille::Array<T>& a, typename brille::ind_t dim){
  return a.sum(dim);
}
template<class T>
LDVec<T> sum(const LDVec<T>& a, typename brille::ind_t dim){
  return LDVec<T>(a.get_lattice(), a.sum(dim));
}
template<class T>
LQVec<T> sum(const LQVec<T>& a, typename brille::ind_t dim){
  return LQVec<T>(a.get_lattice(), a.sum(dim));
}
// template<class T, template<class> class L,
//          typename=typename std::enable_if<
//          std::is_base_of<brille::Array<T>,L<T>>::value &&
//         !std::is_base_of<LatVec,L<T>>::value
//          >::type>
// L<T> sum(const L<T>& a, typename brille::ind_t dim=0){
//   return a.sum(dim);
// }
// // sum(LatVec)
// template<class T, template<class> class L, typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type>
// L<T> sum(const L<T>& a, typename brille::ind_t dim){
//   return L<T>(a.get_lattice(), a.sum(dim));
// }

// abs(brille::Array)
template<class T, template<class> class L,
         typename=typename std::enable_if<std::is_base_of<brille::Array<T>,L<T>>::value&&!std::is_base_of<LatVec,L<T>>::value>::type
        >
L<T> abs(const L<T>& a){
  L<T> out(a.shape());
  for (auto x : SubIt(a.shape())) out[x] = magnitude(a[x]);
  return out;
}
// // abs(LatVec)
// template<class T, template<class> class L, typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type>
// L<T> abs(const L<T>& a){
//   L<T> out(a.get_lattice(), a.size());
//   for (size_t i=0; i<a.size(); ++i) for (size_t j=0; j<a.numel(); ++j)
//   out.insert(std::abs(a.getvalue(i,j), i,j));
//   return out;
// }

// LatVec {+,-,*,/} brille::Array
#define LATVEC_ARRAY_BINARY_OP(X) \
template<class T, class R, template<class> class L, template<class> class A, class S = std::common_type_t<T,R>>\
std::enable_if_t<(std::is_base_of_v<LatVec,L<T>>&&(std::is_base_of_v<brille::Array<R>,A<R>>&&!std::is_base_of_v<LatVec,A<R>>)), L<S>>\
operator X (const L<T>& a, const A<R>& b){\
  auto itr = BroadcastIt(a.shape(), b.shape());\
  L<S> out(a.get_lattice(), itr.shape());\
  for (auto [ox, ax, bx]: itr) out[ox] = a[ax] X b[bx];\
  return out;\
}
LATVEC_ARRAY_BINARY_OP(+)
LATVEC_ARRAY_BINARY_OP(-)
LATVEC_ARRAY_BINARY_OP(*)
LATVEC_ARRAY_BINARY_OP(/)
#undef LATVEC_ARRAY_BINARY_OP

// brille::Array {+,-,*,/} LatVec
#define ARRAY_LATVEC_BINARY_OP(X) \
template<class T, class R, template<class> class L, template<class> class A, class S = std::common_type_t<T,R>>\
std::enable_if_t<(std::is_base_of_v<LatVec,L<T>>&&(std::is_base_of_v<brille::Array<R>,A<R>>&&!std::is_base_of_v<LatVec,A<R>>)), L<S>>\
operator X (const A<R>& b, const L<T>& a){\
  auto itr = BroadcastIt(b.shape(), a.shape());\
  L<S> out(a.get_lattice(), itr.shape());\
  for (auto [ox, bx, ax]: itr) out[ox] = b[bx] X a[ax];\
  return out;\
}
ARRAY_LATVEC_BINARY_OP(+)
ARRAY_LATVEC_BINARY_OP(-)
ARRAY_LATVEC_BINARY_OP(*)
ARRAY_LATVEC_BINARY_OP(/)
#undef ARRAY_LATVEC_BINARY_OP
