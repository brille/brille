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

// [ArrayVector] dot
template<class T, class R, template<class> class L1, template<class> class L2,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,L1<T>>::value&&!std::is_base_of<LatVec,L1<T>>::value>::type,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<R>,L2<R>>::value&&!std::is_base_of<LatVec,L2<R>>::value>::type,
         class S = typename std::common_type<T,R>::type>
ArrayVector<S> dot(const L1<T> &a, const L2<R> &b){
  AVSizeInfo si = a.consistency_check(b);
  if (si.scalara||si.scalarb) throw std::runtime_error("dot product requires two vectors");
  ArrayVector<double> out(1,si.n);
  double tmp=0;
  for (size_t i=0; i<si.n; i++) {
    for (size_t j=0; j<si.m; j++) tmp += a.getvalue(si.oneveca?0:i,j) * b.getvalue(si.onevecb?0:i,j);
    out.insert(tmp, i);
  }
  return out;
}
// [ArrayVector] norm
template<class T, template<class> class L,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,L<T>>::value>::type
        >
ArrayVector<double> norm(const L<T> &a){
  ArrayVector<double> out = dot(a,a);
  for (size_t i=0; i<out.size(); ++i) out.insert(sqrt(out.getvalue(i)),i);
  return out;
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
// template<typename T,typename R> double same_lattice_dot(const R& x, const T& y, const double* len, const double* ang){
//   double out = double(x[0])*double(y[0])*len[0]*len[0]
//              + double(x[1])*double(y[1])*len[1]*len[1]
//              + double(x[2])*double(y[2])*len[2]*len[2]
//              + (double(x[0])*double(y[1])+double(x[1])*double(y[0]))*len[0]*len[1]*cos(ang[2])
//              + (double(x[1])*double(y[2])+double(x[2])*double(y[1]))*len[1]*len[2]*cos(ang[0])
//              + (double(x[2])*double(y[0])+double(x[0])*double(y[2]))*len[2]*len[0]*cos(ang[1]);
//   return out;
// }

// cross (LatVec × LatVec)
template<class T, class R, template<class> class L,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type
        >
L<double> cross(const L<T>& a, const L<R>& b) {
  assert( a.samelattice(b) );
  AVSizeInfo si = a.consistency_check(b);
  if (si.m != 3u) throw std::runtime_error("cross product is only defined for three vectors");
  typename LatticeTraits<L<T>>::type lat = a.get_lattice();
  typename LatVecTraits<L<T>,double>::star tmp( lat.star(), si.n);
  for (size_t i=0; i<si.n; i++)
    vector_cross<double,T,R,3>(tmp.data(i), a.data(si.oneveca?0:i), b.data(si.onevecb?0:i));
  tmp *= lat.get_volume()/2.0/PI;
  return tmp.star();
}

// dot [LatVec ⋅ LatVec]
template<class T, class R, template<class> class L1, template<class> class L2,
         typename=typename std::enable_if<std::is_base_of<LatVec,L1<T>>::value>::type,
         typename=typename std::enable_if<std::is_base_of<LatVec,L2<R>>::value>::type
        >
ArrayVector<double> dot(const L1<T> &a, const L2<R> &b){
  bool issame = a.samelattice(b);
  if (!( issame || a.starlattice(b) )){
    debug_update("Incompatible lattices\n",a.get_lattice().string_repr(),"\n",b.get_lattice().string_repr());
    throw std::runtime_error("the dot product between Lattice Vectors requires same or starred lattices");
  }
  AVSizeInfo si = a.consistency_check(b);
  if (si.m!=3u) throw std::runtime_error("Lattice dot product is only defined for three vectors");
  if (si.scalara||si.scalarb) throw std::runtime_error("Lattice dot product requires two three-vectors");
  ArrayVector<double> out(1,si.n);
  if (issame){
    typename LatticeTraits<L1<T>>::type lat = a.get_lattice();
    double len[3] = {lat.get_a(),lat.get_b(),lat.get_c()};
    double ang[3] = {lat.get_alpha(),lat.get_beta(),lat.get_gamma()};
    for (size_t i=0; i<si.n; i++)
      out.insert( same_lattice_dot( a.data(si.oneveca?0:i), b.data(si.onevecb?0:i), len, ang), i);
  } else {
    double tmp=0;
    for (size_t i=0; i<si.n; i++) {
      for (size_t j=0; j<si.m; j++) tmp += a.getvalue(si.oneveca?0:i,j) * b.getvalue(si.onevecb?0:i,j);
      out.insert(2*PI*tmp, i);
    }
  }
  return out;
}

// // [LatVec] norm
// template<class T, template<class> class L,
//          typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type
//         >
// ArrayVector<double> norm(const L<T> &a){
//   ArrayVector<double> out = dot(a,a);
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
  L<S> out(a.get_lattice(), a.size()+b.size());
  for (size_t i=0; i<a.size(); ++i) for (size_t j=0; j<3u; ++j)
    out.insert( static_cast<S>(a.getvalue(i,j)), i, j);
  for (size_t i=0; i<b.size(); ++i) for (size_t j=0; j<3u; ++j)
    out.insert( static_cast<S>(b.getvalue(i,j)), a.size()+i, j);
  return out;
}

// star
template<class T, template <class R> class L,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type>
typename LatVecTraits<L<T>,double>::star star(const L<T>& v){
  double *cvmt = nullptr;
  cvmt = new double[9]();
  v.get_lattice().get_covariant_metric_tensor(cvmt);
  typename LatVecTraits<L<T>,double>::star vstar( v.get_lattice().star(), v.size() );
  for (size_t i=0; i<v.size(); ++i)
    multiply_matrix_vector(vstar.data(i), cvmt, v.data(i));
  delete[] cvmt;
  return vstar;
}

// LatVec + ArrayVec
template<class T, class R, template<class> class L, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<R>,A<R>>::value&&!std::is_base_of<LatVec,A<R>>::value>::type,
         class S = typename std::common_type<T,R>::type>
L<S> operator+(const L<T>& a, const A<R>& b){
  AVSizeInfo si = a.consistency_check(b);
  if (si.m != 3u) throw std::runtime_error("lattice vectors should always have numel()==3");
  L<S> out(a.get_lattice(),si.n);
  for (size_t i=0; i<si.n; ++i)
    for (size_t j=0; j<si.m; ++j)
      out.insert( a.getvalue(si.oneveca?0:i,j) + b.getvalue(si.onevecb?0:i,si.scalarb?0:j), i,j);
  return out;
}
// LatVec - ArrayVec
template<class T, class R, template<class> class L, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<R>,A<R>>::value&&!std::is_base_of<LatVec,A<R>>::value>::type,
         class S = typename std::common_type<T,R>::type>
L<S> operator-(const L<T>& a, const A<R>& b){
  AVSizeInfo si = a.consistency_check(b);
  if (si.m != 3u) throw std::runtime_error("lattice vectors should always have numel()==3");
  L<S> out(a.get_lattice(),si.n);
  for (size_t i=0; i<si.n; ++i)
    for (size_t j=0; j<si.m; ++j)
      out.insert( a.getvalue(si.oneveca?0:i,j) - b.getvalue(si.onevecb?0:i,si.scalarb?0:j), i,j);
  return out;
}
// LatVec * ArrayVec
template<class T, class R, template<class> class L, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<R>,A<R>>::value&&!std::is_base_of<LatVec,A<R>>::value>::type,
         class S = typename std::common_type<T,R>::type>
L<S> operator*(const L<T>& a, const A<R>& b){
  AVSizeInfo si = a.consistency_check(b);
  if (si.m != 3u) throw std::runtime_error("lattice vectors should always have numel()==3");
  L<S> out(a.get_lattice(),si.n);
  for (size_t i=0; i<si.n; ++i)
    for (size_t j=0; j<si.m; ++j)
      out.insert( a.getvalue(si.oneveca?0:i,j) * b.getvalue(si.onevecb?0:i,si.scalarb?0:j), i,j);
  return out;
}
// LatVec / ArrayVec
template<class T, class R, template<class> class L, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type,
         // typename=typename std::enable_if<std::is_base_of<ArrayVector<R>,A<R>>::value&&!std::is_base_of<LatVec,A<R>>::value>::type,
         typename=typename std::enable_if<!std::is_base_of<LatVec,A<R>>::value>::type,
         class S = typename std::common_type<T,R>::type,
         typename=typename std::enable_if<std::is_floating_point<S>::value>::type>
L<S> operator/(const L<T>& a, const A<R>& b){
  AVSizeInfo si = a.consistency_check(b);
  if (si.m != 3u) throw std::runtime_error("lattice vectors should always have numel()==3");
  L<S> out(a.get_lattice(),si.n);
  for (size_t i=0; i<si.n; ++i)
    for (size_t j=0; j<si.m; ++j)
      out.insert( a.getvalue(si.oneveca?0:i,j) / b.getvalue(si.onevecb?0:i,si.scalarb?0:j), i,j);
  return out;
}

// ArrayVec + LatVec
template<class T, class R, template<class> class L, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<R>,A<R>>::value&&!std::is_base_of<LatVec,A<R>>::value>::type,
         class S = typename std::common_type<T,R>::type>
L<S> operator+(const A<R>& b, const L<T>& a){
  AVSizeInfo si = a.consistency_check(b);
  if (si.m != 3u) throw std::runtime_error("lattice vectors should always have numel()==3");
  L<S> out(a.get_lattice(),si.n);
  for (size_t i=0; i<si.n; ++i)
    for (size_t j=0; j<si.m; ++j)
      out.insert( b.getvalue(si.onevecb?0:i,si.scalarb?0:j) + a.getvalue(si.oneveca?0:i,j), i,j);
  return out;
}
// ArrayVec - LatVec
template<class T, class R, template<class> class L, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<R>,A<R>>::value&&!std::is_base_of<LatVec,A<R>>::value>::type,
         class S = typename std::common_type<T,R>::type>
L<S> operator-(const A<R>& b, const L<T>& a){
  AVSizeInfo si = a.consistency_check(b);
  if (si.m != 3u) throw std::runtime_error("lattice vectors should always have numel()==3");
  L<S> out(a.get_lattice(),si.n);
  for (size_t i=0; i<si.n; ++i)
    for (size_t j=0; j<si.m; ++j)
      out.insert( b.getvalue(si.onevecb?0:i,si.scalarb?0:j) - a.getvalue(si.oneveca?0:i,j), i,j);
  return out;
}
// LatVec * ArrayVec
template<class T, class R, template<class> class L, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<R>,A<R>>::value&&!std::is_base_of<LatVec,A<R>>::value>::type,
         class S = typename std::common_type<T,R>::type>
L<S> operator*(const A<R>& b, const L<T>& a){
  AVSizeInfo si = a.consistency_check(b);
  if (si.m != 3u) throw std::runtime_error("lattice vectors should always have numel()==3");
  L<S> out(a.get_lattice(),si.n);
  for (size_t i=0; i<si.n; ++i)
    for (size_t j=0; j<si.m; ++j)
      out.insert( b.getvalue(si.onevecb?0:i,si.scalarb?0:j) * a.getvalue(si.oneveca?0:i,j), i,j);
  return out;
}
// LatVec / ArrayVec
template<class T, class R, template<class> class L, template<class> class A,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<R>,A<R>>::value&&!std::is_base_of<LatVec,A<R>>::value>::type,
         class S = typename std::common_type<T,R>::type,
         typename=typename std::enable_if<std::is_floating_point<S>::value>::type>
L<S> operator/(const A<R>& b, const L<T>& a){
  if (!(a.samelattice(b)))
    throw std::runtime_error("Adding lattice vectors requires they're represented in the same lattice");
  AVSizeInfo si = a.consistency_check(b);
  if (si.m != 3u) throw std::runtime_error("lattice vectors should always have numel()==3");
  L<S> out(a.get_lattice(),si.n);
  for (size_t i=0; i<si.n; ++i)
    for (size_t j=0; j<si.m; ++j)
      out.insert( b.getvalue(si.onevecb?0:i,si.scalarb?0:j) / a.getvalue(si.oneveca?0:i,j), i,j);
  return out;
}

// Matrix * LatVec
template<class T, class R, template<class> class L,
         typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type,
         class S = typename std::common_type<T,R>::type,
         typename=typename std::enable_if<std::is_floating_point<S>::value>::type>
L<S> operator*(const std::array<R,9>& m, const L<T>& a){
  L<S> out(a.get_lattice(), a.size());
  S tmp[3];
  for (size_t i=0; i<a.size(); ++i)
    multiply_matrix_vector<S,R,T,3>(out.data(i), m.data(), a.data(i) );
  return out;
}

// -LatVec
template<class T, template<class> class L, typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type>
L<T> operator-(const L<T>& a){
  L<T> out(a.get_lattice(),a.size());
  for (size_t i=0; i<a.size(); ++i)
    for (size_t j=0; j<a.numel(); ++j)
      out.insert( -a.getvalue(i,j), i,j);
  return out;
}


// sum(ArrayVector)
template<class T, template<class> class L,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,L<T>>::value&&!std::is_base_of<LatVec,L<T>>::value>::type
        >
L<T> sum(const L<T>& a, const int dim=0){
  T tmp;
  L<T> out;
  switch (dim){
    case 1:
      out.refresh(1u, a.size());
      for (size_t i=0; i<a.size(); ++i){
        tmp = T(0);
        for (size_t j=0; j<a.numel(); ++j) tmp += a.getvalue(i,j);
        out.insert(tmp, i, 0);
      }
      break;
    default:
      out.refresh(a.numel(), 1u);
      for (size_t j=0; j<a.numel(); ++j){
        tmp = T(0);
        for (size_t i=0; i<a.size(); ++i) tmp += a.getvalue(i, j);
        out.insert(tmp, 0, j);
      }
      break;
    }
return out;
}
// sum(LatVec)
template<class T, template<class> class L, typename=typename std::enable_if<std::is_base_of<LatVec,L<T>>::value>::type>
L<T> sum(const L<T>& a){
  T tmp;
  L<T> out(a.get_lattice(), 1u);
  for (size_t j=0; j<a.numel(); ++j){
    tmp = T(0);
    for (size_t i=0; i<a.size(); ++i) tmp += a.getvalue(i, j);
    out.insert(tmp, 0, j);
  }
  return out;
}

// abs(ArrayVector)
template<class T, template<class> class L,
         typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,L<T>>::value&&!std::is_base_of<LatVec,L<T>>::value>::type
        >
L<T> abs(const L<T>& a){
  L<T> out(a.numel(), a.size());
  for (size_t i=0; i<a.size(); ++i) for (size_t j=0; j<a.numel(); ++j)
    out.insert(magnitude(a.getvalue(i,j)), i,j);
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
