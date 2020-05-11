// // In Place arithmetic ArrayVector +-*/ ArrayVector
// template<typename T> ArrayVector<T>& ArrayVector<T>:: operator += (const ArrayVector<T>& av){
//   AVSizeInfo si = this->inplace_consistency_check(av);
//   unsigned flg = static_cast<unsigned>(si.onevecb) + (static_cast<unsigned>(si.singular)<<1u);
//   switch(flg){
//     case 3: /*singular && onvecb*/
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) + av.getvalue(0,0), i,j); break;
//     case 2: /*singular*/
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) + av.getvalue(i,0), i,j); break;
//     case 1: /*onevecb*/
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) + av.getvalue(0,j), i,j); break;
//     default:
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) + av.getvalue(i,j), i,j);
//   }
//   return *this;
// }
// template<typename T> ArrayVector<T>& ArrayVector<T>:: operator -= (const ArrayVector<T>& av){
//   AVSizeInfo si = this->inplace_consistency_check(av);
//   unsigned flg = static_cast<unsigned>(si.onevecb) + (static_cast<unsigned>(si.singular)<<1u);
//   switch(flg){
//     case 3: /*singular && onvecb*/
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) - av.getvalue(0,0), i,j); break;
//     case 2: /*singular*/
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) - av.getvalue(i,0), i,j); break;
//     case 1: /*onevecb*/
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) - av.getvalue(0,j), i,j); break;
//     default:
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) - av.getvalue(i,j), i,j);
//   }
//   return *this;
// }
// template<typename T> ArrayVector<T>& ArrayVector<T>:: operator *= (const ArrayVector<T>& av){
//   AVSizeInfo si = this->inplace_consistency_check(av);
//   unsigned flg = static_cast<unsigned>(si.onevecb) + (static_cast<unsigned>(si.singular)<<1u);
//   switch(flg){
//     case 3: /*singular && onvecb*/
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) * av.getvalue(0,0), i,j); break;
//     case 2: /*singular*/
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) * av.getvalue(i,0), i,j); break;
//     case 1: /*onevecb*/
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) * av.getvalue(0,j), i,j); break;
//     default:
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) * av.getvalue(i,j), i,j);
//   }
//   return *this;
// }
// template<typename T> ArrayVector<T>& ArrayVector<T>:: operator /= (const ArrayVector<T>& av){
//   AVSizeInfo si = this->inplace_consistency_check(av);
//   unsigned flg = static_cast<unsigned>(si.onevecb) + (static_cast<unsigned>(si.singular)<<1u);
//   switch(flg){
//     case 3: /*singular && onvecb*/
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) / av.getvalue(0,0), i,j); break;
//     case 2: /*singular*/
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) / av.getvalue(i,0), i,j); break;
//     case 1: /*onevecb*/
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) / av.getvalue(0,j), i,j); break;
//     default:
//     #pragma omp simd collapse(2)
//     for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//       this->insert(this->getvalue(i,j) / av.getvalue(i,j), i,j);
//   }
//   return *this;
// }
//
// // In-place binary operators with scalars
// template<typename T> ArrayVector<T>& ArrayVector<T>:: operator += (const T& av){
//   #pragma omp simd collapse(2)
//   for (size_t i=0; i<this->size(); ++i) for (size_t j=0; j<this->numel(); ++j)
//     this->insert(this->getvalue(i,j) + av, i,j);
//   return *this;
// }
// template<typename T> ArrayVector<T>& ArrayVector<T>:: operator -= (const T& av){
//   #pragma omp simd collapse(2)
//   for (size_t i=0; i<this->size(); ++i) for (size_t j=0; j<this->numel(); ++j)
//     this->insert(this->getvalue(i,j) - av, i,j);
//   return *this;
// }
// template<typename T> ArrayVector<T>& ArrayVector<T>:: operator *= (const T& av){
//   #pragma omp simd collapse(2)
//   for (size_t i=0; i<this->size(); ++i) for (size_t j=0; j<this->numel(); ++j)
//     this->insert(this->getvalue(i,j) * av, i,j);
//   return *this;
// }
// template<typename T> ArrayVector<T>& ArrayVector<T>:: operator /= (const T& av){
//   #pragma omp simd collapse(2)
//   for (size_t i=0; i<this->size(); ++i) for (size_t j=0; j<this->numel(); ++j)
//     this->insert(this->getvalue(i,j) / av, i,j);
//   return *this;
// }
//
// // binary ArrayVector ArrayVector operators:
// template<class T, class R, template<class> class A,typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,class S = typename std::common_type<T,R>::type>
// A<S> operator + (const A<T>& a, const A<R>& b){
//   AVSizeInfo si = a.consistency_check(b);
//   A<S> out = si.aorb ? A<S>(a) : A<S>(b);
//   out.refresh(si.m,si.n); /*in case a.size == b.size but one is singular, or a.numel == b.numel but one is scalar*/
//   unsigned scalarflg = static_cast<unsigned>(si.scalara) + (static_cast<unsigned>(si.scalarb)<<1u);
//   unsigned onevecflg = static_cast<unsigned>(si.oneveca) + (static_cast<unsigned>(si.onevecb)<<1u);
//   switch (scalarflg){
//     case 3u: /*a and b are scalars*/
//     switch (onevecflg){
//       case 3u:
//       out.insert(a.getvalue(0,0) + b.getvalue(0,0), 0,0);
//       break;
//       case 2u:
//       #pragma omp simd
//       for (size_t i=0; i<si.n; ++i)
//         out.insert(a.getvalue(i,0) + b.getvalue(0,0), i,0);
//       break;
//       case 1u:
//       #pragma omp simd
//       for (size_t i=0; i<si.n; ++i)
//         out.insert(a.getvalue(0,0) + b.getvalue(i,0), i,0);
//       break;
//       default:
//       #pragma omp simd
//       for (size_t i=0; i<si.n; ++i)
//         out.insert(a.getvalue(i,0) + b.getvalue(i,0), i,0);
//     }
//     break;
//     case 2u: /*b is a scalar*/
//     switch (onevecflg){
//       case 3u:
//       #pragma omp simd
//       for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) + b.getvalue(0,0), 0,j);
//       break;
//       case 2u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) + b.getvalue(0,0), i,j);
//       break;
//       case 1u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) + b.getvalue(i,0), i,j);
//       break;
//       default:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) + b.getvalue(i,0), i,j);
//     }
//     break;
//     case 1u: /*a is a scalar*/
//     switch (onevecflg){
//       case 3u:
//       #pragma omp simd
//       for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,0) + b.getvalue(0,j), 0,j);
//       break;
//       case 2u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,0) + b.getvalue(0,j), i,j);
//       break;
//       case 1u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,0) + b.getvalue(i,j), i,j);
//       break;
//       default:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,0) + b.getvalue(i,j), i,j);
//     }
//     break;
//     default: /*neither a nor b are scalars*/
//     switch (onevecflg){
//       case 3u:
//       #pragma omp simd
//       for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) + b.getvalue(0,j), 0,j);
//       break;
//       case 2u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) + b.getvalue(0,j), i,j);
//       break;
//       case 1u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) + b.getvalue(i,j), i,j);
//       break;
//       default:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) + b.getvalue(i,j), i,j);
//     }
//   }
//   return out;
// }
// template<class T, class R, template<class> class A,typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,class S = typename std::common_type<T,R>::type>
// A<S> operator - (const A<T>& a, const A<R>& b){
//   AVSizeInfo si = a.consistency_check(b);
//   A<S> out = si.aorb ? A<S>(a) : A<S>(b);
//   out.refresh(si.m,si.n); /*in case a.size == b.size but one is singular, or a.numel == b.numel but one is scalar*/
//   unsigned scalarflg = static_cast<unsigned>(si.scalara) + (static_cast<unsigned>(si.scalarb)<<1u);
//   unsigned onevecflg = static_cast<unsigned>(si.oneveca) + (static_cast<unsigned>(si.onevecb)<<1u);
//   switch (scalarflg){
//     case 3u: /*a and b are scalars*/
//     switch (onevecflg){
//       case 3u:
//       out.insert(a.getvalue(0,0) - b.getvalue(0,0), 0,0);
//       break;
//       case 2u:
//       #pragma omp simd
//       for (size_t i=0; i<si.n; ++i)
//         out.insert(a.getvalue(i,0) - b.getvalue(0,0), i,0);
//       break;
//       case 1u:
//       #pragma omp simd
//       for (size_t i=0; i<si.n; ++i)
//         out.insert(a.getvalue(0,0) - b.getvalue(i,0), i,0);
//       break;
//       default:
//       #pragma omp simd
//       for (size_t i=0; i<si.n; ++i)
//         out.insert(a.getvalue(i,0) - b.getvalue(i,0), i,0);
//     }
//     break;
//     case 2u: /*b is a scalar*/
//     switch (onevecflg){
//       case 3u:
//       #pragma omp simd
//       for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) - b.getvalue(0,0), 0,j);
//       break;
//       case 2u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) - b.getvalue(0,0), i,j);
//       break;
//       case 1u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) - b.getvalue(i,0), i,j);
//       break;
//       default:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) - b.getvalue(i,0), i,j);
//     }
//     break;
//     case 1u: /*a is a scalar*/
//     switch (onevecflg){
//       case 3u:
//       #pragma omp simd
//       for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,0) - b.getvalue(0,j), 0,j);
//       break;
//       case 2u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,0) - b.getvalue(0,j), i,j);
//       break;
//       case 1u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,0) - b.getvalue(i,j), i,j);
//       break;
//       default:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,0) - b.getvalue(i,j), i,j);
//     }
//     break;
//     default: /*neither a nor b are scalars*/
//     switch (onevecflg){
//       case 3u:
//       #pragma omp simd
//       for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) - b.getvalue(0,j), 0,j);
//       break;
//       case 2u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) - b.getvalue(0,j), i,j);
//       break;
//       case 1u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) - b.getvalue(i,j), i,j);
//       break;
//       default:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) - b.getvalue(i,j), i,j);
//     }
//   }
//   return out;
// }
// template<class T, class R, template<class> class A,typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,class S = typename std::common_type<T,R>::type>
// A<S> operator * (const A<T>& a, const A<R>& b){
//   AVSizeInfo si = a.consistency_check(b);
//   A<S> out = si.aorb ? A<S>(a) : A<S>(b);
//   out.refresh(si.m,si.n); /*in case a.size == b.size but one is singular, or a.numel == b.numel but one is scalar*/
//   unsigned scalarflg = static_cast<unsigned>(si.scalara) + (static_cast<unsigned>(si.scalarb)<<1u);
//   unsigned onevecflg = static_cast<unsigned>(si.oneveca) + (static_cast<unsigned>(si.onevecb)<<1u);
//   switch (scalarflg){
//     case 3u: /*a and b are scalars*/
//     switch (onevecflg){
//       case 3u:
//       out.insert(a.getvalue(0,0) * b.getvalue(0,0), 0,0);
//       break;
//       case 2u:
//       #pragma omp simd
//       for (size_t i=0; i<si.n; ++i)
//         out.insert(a.getvalue(i,0) * b.getvalue(0,0), i,0);
//       break;
//       case 1u:
//       #pragma omp simd
//       for (size_t i=0; i<si.n; ++i)
//         out.insert(a.getvalue(0,0) * b.getvalue(i,0), i,0);
//       break;
//       default:
//       #pragma omp simd
//       for (size_t i=0; i<si.n; ++i)
//         out.insert(a.getvalue(i,0) * b.getvalue(i,0), i,0);
//     }
//     break;
//     case 2u: /*b is a scalar*/
//     switch (onevecflg){
//       case 3u:
//       #pragma omp simd
//       for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) * b.getvalue(0,0), 0,j);
//       break;
//       case 2u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) * b.getvalue(0,0), i,j);
//       break;
//       case 1u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) * b.getvalue(i,0), i,j);
//       break;
//       default:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) * b.getvalue(i,0), i,j);
//     }
//     break;
//     case 1u: /*a is a scalar*/
//     switch (onevecflg){
//       case 3u:
//       #pragma omp simd
//       for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,0) * b.getvalue(0,j), 0,j);
//       break;
//       case 2u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,0) * b.getvalue(0,j), i,j);
//       break;
//       case 1u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,0) * b.getvalue(i,j), i,j);
//       break;
//       default:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,0) * b.getvalue(i,j), i,j);
//     }
//     break;
//     default: /*neither a nor b are scalars*/
//     switch (onevecflg){
//       case 3u:
//       #pragma omp simd
//       for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) * b.getvalue(0,j), 0,j);
//       break;
//       case 2u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) * b.getvalue(0,j), i,j);
//       break;
//       case 1u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) * b.getvalue(i,j), i,j);
//       break;
//       default:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) * b.getvalue(i,j), i,j);
//     }
//   }
//   return out;
// }
// template<class T, class R, template<class> class A,typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,class S = typename std::common_type<T,R>::type>
// A<S> operator / (const A<T>& a, const A<R>& b){
//   AVSizeInfo si = a.consistency_check(b);
//   A<S> out = si.aorb ? A<S>(a) : A<S>(b);
//   out.refresh(si.m,si.n); /*in case a.size == b.size but one is singular, or a.numel == b.numel but one is scalar*/
//   unsigned scalarflg = static_cast<unsigned>(si.scalara) + (static_cast<unsigned>(si.scalarb)<<1u);
//   unsigned onevecflg = static_cast<unsigned>(si.oneveca) + (static_cast<unsigned>(si.onevecb)<<1u);
//   switch (scalarflg){
//     case 3u: /*a and b are scalars*/
//     switch (onevecflg){
//       case 3u:
//       out.insert(a.getvalue(0,0) / b.getvalue(0,0), 0,0);
//       break;
//       case 2u:
//       #pragma omp simd
//       for (size_t i=0; i<si.n; ++i)
//         out.insert(a.getvalue(i,0) / b.getvalue(0,0), i,0);
//       break;
//       case 1u:
//       #pragma omp simd
//       for (size_t i=0; i<si.n; ++i)
//         out.insert(a.getvalue(0,0) / b.getvalue(i,0), i,0);
//       break;
//       default:
//       #pragma omp simd
//       for (size_t i=0; i<si.n; ++i)
//         out.insert(a.getvalue(i,0) / b.getvalue(i,0), i,0);
//     }
//     break;
//     case 2u: /*b is a scalar*/
//     switch (onevecflg){
//       case 3u:
//       #pragma omp simd
//       for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) / b.getvalue(0,0), 0,j);
//       break;
//       case 2u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) / b.getvalue(0,0), i,j);
//       break;
//       case 1u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) / b.getvalue(i,0), i,j);
//       break;
//       default:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) / b.getvalue(i,0), i,j);
//     }
//     break;
//     case 1u: /*a is a scalar*/
//     switch (onevecflg){
//       case 3u:
//       #pragma omp simd
//       for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,0) / b.getvalue(0,j), 0,j);
//       break;
//       case 2u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,0) / b.getvalue(0,j), i,j);
//       break;
//       case 1u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,0) / b.getvalue(i,j), i,j);
//       break;
//       default:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,0) / b.getvalue(i,j), i,j);
//     }
//     break;
//     default: /*neither a nor b are scalars*/
//     switch (onevecflg){
//       case 3u:
//       #pragma omp simd
//       for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) / b.getvalue(0,j), 0,j);
//       break;
//       case 2u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) / b.getvalue(0,j), i,j);
//       break;
//       case 1u:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(0,j) / b.getvalue(i,j), i,j);
//       break;
//       default:
//       #pragma omp simd collapse(2)
//       for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j)
//         out.insert(a.getvalue(i,j) / b.getvalue(i,j), i,j);
//     }
//   }
//   return out;
// }

/*-----------------------------------------------------------------------------|
|        Unused preprocessor macros to make the above more manageable          |
|------------------------------------------------------------------------------|
|    The macros can not be used along with preprocessor directives necessary   |
|    to implement simd loop vectorization, sadly.                              |
|-----------------------------------------------------------------------------*/

// In Place arithmetic ArrayVector +-*/ ArrayVector
#define INPLACE_ARRAYVECTOR_ARRAYVECTOR_OPERATOR(X,Y) template<typename T> ArrayVector<T>& ArrayVector<T>:: operator X (const ArrayVector<T>& av){\
  AVSizeInfo si = this->inplace_consistency_check(av);\
  unsigned flg = static_cast<unsigned>(si.onevecb) + (static_cast<unsigned>(si.singular)<<1u);\
  switch(flg){\
    case 3: /*singular && onvecb*/\
    for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) this->insert(this->getvalue(i,j) Y av.getvalue(0,0), i,j); break;\
    case 2: /*singular*/\
    for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) this->insert(this->getvalue(i,j) Y av.getvalue(i,0), i,j); break;\
    case 1: /*onevecb*/\
    for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) this->insert(this->getvalue(i,j) Y av.getvalue(0,j), i,j); break;\
    default:\
    for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) this->insert(this->getvalue(i,j) Y av.getvalue(i,j), i,j);\
  }\
  return *this;\
}
// In-place binary operators with scalars
#define INPLACE_ARRAYVECTOR_SCALAR_OPERATOR(X,Y) template<typename T> ArrayVector<T>& ArrayVector<T>:: operator X (const T& av){\
  for (size_t i=0; i<this->size(); ++i) for (size_t j=0; j<this->numel(); ++j) this->insert(this->getvalue(i,j) Y av, i,j);\
  return *this;\
}

#define BINARY_ARRAYVECTOR_ARRAYVECTOR_OPERATOR(X) template<class T, class R, template<class> class A,typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type,class S = typename std::common_type<T,R>::type>\
A<S> operator X (const A<T>& a, const A<R>& b){\
  AVSizeInfo si = a.consistency_check(b);\
  A<S> out = si.aorb ? A<S>(a) : A<S>(b);\
  out.refresh(si.m,si.n); /*in case a.size == b.size but one is singular, or a.numel == b.numel but one is scalar*/\
  unsigned scalarflg = static_cast<unsigned>(si.scalara) + (static_cast<unsigned>(si.scalarb)<<1u);\
  unsigned onevecflg = static_cast<unsigned>(si.oneveca) + (static_cast<unsigned>(si.onevecb)<<1u);\
  switch (scalarflg){\
    case 3u: /*a and b are scalars*/\
    switch (onevecflg){\
      case 3u: out.insert(a.getvalue(0,0) X b.getvalue(0,0), 0,0); break;\
      case 2u: for (size_t i=0; i<si.n; ++i) out.insert(a.getvalue(i,0) X b.getvalue(0,0), i,0); break;\
      case 1u: for (size_t i=0; i<si.n; ++i) out.insert(a.getvalue(0,0) X b.getvalue(i,0), i,0); break;\
      default: for (size_t i=0; i<si.n; ++i) out.insert(a.getvalue(i,0) X b.getvalue(i,0), i,0);\
    }\
    break;\
    case 2u: /*b is a scalar*/\
    switch (onevecflg){\
      case 3u: for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(0,j) X b.getvalue(0,0), 0,j); break;\
      case 2u: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(i,j) X b.getvalue(0,0), i,j); break;\
      case 1u: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(0,j) X b.getvalue(i,0), i,j); break;\
      default: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(i,j) X b.getvalue(i,0), i,j);\
    }\
    break;\
    case 1u: /*a is a scalar*/\
    switch (onevecflg){\
      case 3u: for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(0,0) X b.getvalue(0,j), 0,j); break;\
      case 2u: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(i,0) X b.getvalue(0,j), i,j); break;\
      case 1u: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(0,0) X b.getvalue(i,j), i,j); break;\
      default: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(i,0) X b.getvalue(i,j), i,j);\
    }\
    break;\
    default: /*neither a nor b are scalars*/\
    switch (onevecflg){\
      case 3u: for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(0,j) X b.getvalue(0,j), 0,j); break;\
      case 2u: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(i,j) X b.getvalue(0,j), i,j); break;\
      case 1u: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(0,j) X b.getvalue(i,j), i,j); break;\
      default: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(i,j) X b.getvalue(i,j), i,j);\
    }\
  }\
  return out;\
}

#define BINARY_ARRAYVECTOR_STDARRAY_OPERATOR(X) template<class T, class R, template<class> class A,typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,A<T>>::value>::type, class S = typename std::common_type<T,R>::type, size_t Nel>\
A<S> operator X (const A<T>& a, const std::array<R,Nel>& b){\
  AVSizeInfo si = a.consistency_check(b);\
  A<S> out = A<S>(a);\
  out.refresh(si.m,si.n); /*in case a.size == b.size but one is singular, or a.numel == b.numel but one is scalar*/\
  unsigned scalarflg = static_cast<unsigned>(si.scalara) + (static_cast<unsigned>(si.scalarb)<<1u);\
  unsigned onevecflg = static_cast<unsigned>(si.oneveca) + (static_cast<unsigned>(si.onevecb)<<1u);\
  std::string interpret_error_msg = "I do not know how to interpret this std::array as an ArrayVector for binary operations";\
  switch (scalarflg){\
    case 3u: /*a and b are scalars*/\
    switch (onevecflg){\
      case 3u: out.insert(a.getvalue(0,0) X b[0], 0,0); break;\
      case 2u: for (size_t i=0; i<si.n; ++i) out.insert(a.getvalue(i,0) X b[0], i,0); break;\
      case 1u: for (size_t i=0; i<si.n; ++i) out.insert(a.getvalue(0,0) X b[i], i,0); break;\
      default: for (size_t i=0; i<si.n; ++i) out.insert(a.getvalue(i,0) X b[i], i,0);\
    }\
    break;\
    case 2u: /*b is a scalar*/\
    switch (onevecflg){\
      case 3u:                               for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(0,j) X b[0], 0,j); break;\
      case 2u: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(i,j) X b[0], i,j); break;\
      case 1u: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(0,j) X b[i], i,j); break;\
      default: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(i,j) X b[i], i,j);\
    }\
    break;\
    case 1u: /*a is a scalar*/\
    switch (onevecflg){\
      case 3u:                               for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(0,0) X b[j], 0,j); break;\
      case 2u: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(i,0) X b[j], i,j); break;\
      default: throw std::runtime_error(interpret_error_msg);\
    }\
    break;\
    default: /*neither a nor b are scalars*/\
    switch (onevecflg){\
      case 3u:                               for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(0,j) X b[j], 0,j); break;\
      case 2u: for (size_t i=0; i<si.n; ++i) for (size_t j=0; j<si.m; ++j) out.insert(a.getvalue(i,j) X b[j], i,j); break;\
      default: throw std::runtime_error(interpret_error_msg);\
    }\
  }\
  return out;\
}

INPLACE_ARRAYVECTOR_ARRAYVECTOR_OPERATOR(+=,+)
INPLACE_ARRAYVECTOR_ARRAYVECTOR_OPERATOR(-=,-)
INPLACE_ARRAYVECTOR_ARRAYVECTOR_OPERATOR(*=,*)
INPLACE_ARRAYVECTOR_ARRAYVECTOR_OPERATOR(/=,/)

INPLACE_ARRAYVECTOR_SCALAR_OPERATOR(+=,+)
INPLACE_ARRAYVECTOR_SCALAR_OPERATOR(-=,-)
INPLACE_ARRAYVECTOR_SCALAR_OPERATOR(*=,*)
INPLACE_ARRAYVECTOR_SCALAR_OPERATOR(/=,/)

BINARY_ARRAYVECTOR_ARRAYVECTOR_OPERATOR(+)
BINARY_ARRAYVECTOR_ARRAYVECTOR_OPERATOR(-)
BINARY_ARRAYVECTOR_ARRAYVECTOR_OPERATOR(*)
BINARY_ARRAYVECTOR_ARRAYVECTOR_OPERATOR(/)

BINARY_ARRAYVECTOR_STDARRAY_OPERATOR(+)
BINARY_ARRAYVECTOR_STDARRAY_OPERATOR(-)
BINARY_ARRAYVECTOR_STDARRAY_OPERATOR(*)
BINARY_ARRAYVECTOR_STDARRAY_OPERATOR(/)
