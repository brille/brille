#ifndef BRILLE_ARRAY_ATTRIBUTES_HPP_
#define BRILLE_ARRAY_ATTRIBUTES_HPP_

namespace brille{

/*! \brief Templated struct to help differentiate between Array2 and its derived
           classes LQVec and LDVec

Since the derived classes add Lattice information some functions behave
differently for 'raw' Array2 versus LQVec or LDVec objects.
This struct is used in template arguments to cause substitution failure-based
differentiation of otherwise identical function signatures.
*/
  template<class... T> struct ArrayTraits{
  static constexpr bool array = false;  //!< is this an Array2, LQVec, or LDVec
  static constexpr bool latvec = false; //!< is this an LQVec or LDVec
};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T> struct ArrayTraits<bArray<T>>{
  static constexpr bool array = true;
  static constexpr bool latvec = false;
};
#endif

/*! Easier access to ArrayTraits for template substitution

true for Array2, LQVec, LDVec
*/
template<class T, template<class> class A>
inline constexpr bool isArray = ArrayTraits<A<T>>::array;
/*! Easier access to ArrayTraits for template substitution

true for LQVec, LDVec
*/
template<class T, template<class> class A>
inline constexpr bool isLatVec = ArrayTraits<A<T>>::latvec;


/*! Easier access to ArrayTraits for template substitution

true for Array2
*/
template<class T, template<class> class A>
inline constexpr bool isBareArray = isArray<T,A> && !isLatVec<T,A>;

/*! Easier access to ArrayTraits for double template substitution

true for any combination of two Array2, LQVec, LDVec objects
*/
template<class T, template<class> class A, class R, template<class> class B>
inline constexpr bool bothArrays = isArray<T,A> && isArray<R,B>;
/*! Easier access to ArrayTraits for double template substitution

true for two Array2 objects
*/
template<class T, template<class> class A, class R, template<class> class B>
inline constexpr bool bareArrays = isBareArray<T,A> && isBareArray<R,B>;
/*! Easier access to ArrayTraits for double template substitution

true for any combination of two LQVec, LDVec objects
*/
template<class T, template<class> class A, class R, template<class> class B>
inline constexpr bool bothLatVecs = isLatVec<T,A> && isLatVec<R,B>;

}

#endif