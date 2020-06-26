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

#ifndef _ARRAYVECTOR_H_
#define _ARRAYVECTOR_H_

// #include <iostream>
// #include <string>
// #include <cmath>
#include <functional>
// #include <vector>
#include <algorithm>
// #include <numeric>
#include <complex> // for +=  support in unsafe_interpolate_to
#include "utilities.hpp"
// #include "debug.h" // ensurses __PRETTY_FUNCTION__ is defined for MSVC, provides debug_update()

/* @enum Comp
   @brief for simplified comparisons of ArrayVector
*/
enum class Comp {
  lt,    //< <
  gt,    //< >
  le,    //< ≤
  ge,    //< ≥
  eq,    //< ==
  nle,   //< ⩽̸
  nge,   //< ⩾̸
  neq,   //< ≠
  le_ge, //< ⪓
  plus,  //< +
  minus, //< -
  times, //< ×
  rdiv,  //< /
  ldiv   //< \.
};

/*!  \brief A class to hold a vector of arrays in contiguous memory

  This class holds N arrays, each with M elements in a N*M contiguous block
  of memory.
*/
template<typename T> class ArrayVector;
class LatVec; // forward declare so that we can prevent operator overloads from applying to LQVec and LDVec

/*! \brief A struct to hold information about two ArrayVectors

  Binary operations between two ArrayVectors require that they have compatible
  dimensions. This struct is the return type of the function(s) which determine
  ArrayVector compatibility and contains information about their sizes.
*/
typedef struct{
  size_t n;     //!< the common number of arrays
  size_t m;     //!< the common number of elements per array
  bool oneveca; //!< is "a" a "SingleVector"
  bool onevecb; //!< is "b" a "SingleVector"
  bool scalara; //!< is "a" an "ArrayScalar"
  bool scalarb; //!< is "b" an "ArrayScalar"
  bool singular;//!< is only "b" an "ArrayScalar"
  bool aorb;    //!< does "a" hold more Vector(/Scalar) elements
} AVSizeInfo;

/**********************************************************
 * The ArrayVector class is intended to hold a continuous *
 * block of equal-length vector-like arrays. Within this  *
 * block, the stride is the size of each array.           *
 * So N arrays of length M has its first array at element *
 * 0, its second at M, its third at 2M, ..., and its last *
 * at (N-1)*M                                             *
 **********************************************************/
 // Replace _data by a std::vector< std::array<T,M> > ... then ArrayVector<T,M> and the standard libraries define most things
template<typename T> class ArrayVector{
protected:
  size_t M; //!< The number of elements within each array
  size_t N; //!< The number of arrays in the ArrayVector
  T* _data;  //!< A pointer to the first element of the contiguous memory _data block
public:
  /*! Standard ArrayVector constructor
      @param m the number of elements within each array
      @param n the number of arrays in the ArrayVector
      @param d [optional] a pointer to a n*m block of data which is copied into the ArrayVector
  */
  ArrayVector(size_t m=0, size_t n=0, const T* d=nullptr) : M(m), N(n), _data(nullptr){
      if (m && n) _data = new T[m*n]();
      if (d && m && n) for(size_t i=0; i<m*n; i++) _data[i] = d[i];
  };
  /*! Single-value initializer ArrayVector constructor
      @param m the number of elements within each array
      @param n the number of arrays in the ArrayVector
      @param init the value which every element of every array is initialized to
  */
  ArrayVector(size_t m, size_t n, const T init): M(m), N(n), _data(nullptr){
    if (m&&n){
      _data =new T[m*n]();
      for (size_t i=0; i<m*n; ++i) _data[i] = init;
    }
  }
  /*! ArrayVector contstructor taking shape and stride information for the case
      of a non-contiguous and/or non-row-ordered array at the provided pointer.
      @param d a pointer to the n*m block of data, to be copied
      @param shape a vector of sizes along each dimension of the input array
      @param strides a vector of strides along each dimension of the input array
      @note The strides vector must indicate the number of bytes necessary to
            move from one point to the next along a given dimension of the input
            data array. The ArrayVector is strictly two dimensional and
            therefore a convention must be adopted for converting the
            D-dimensional input array [D=shape.size()].
            This constructor will attempt to combine the 2ⁿᵈ through Dᵗʰ
            dimensions such that (n,m) = (shape[0],shape[1]*…*shape[D-1])
  */
  template<class Integer, typename=typename std::enable_if<std::is_integral<Integer>::value>::type>
  ArrayVector(const T* d, const std::vector<Integer>& shape, const std::vector<Integer>& strides): _data(nullptr){
    if (shape.size()>0 && shape.size()==strides.size()){
      this->N = static_cast<size_t>(shape[0]);
      Integer nel=1; for (size_t i=1; i<shape.size(); ++i) nel *= shape[i];
      this->M = static_cast<size_t>(nel);
      if (this->N && this->M) this->_data = new T[this->N*this->M]();
      if (d && this->N && this->M){
        // we want to copy-in the _data to a row-ordered array (and flatten it)
        // so calculate the span along each dimension of that row-ordered array
        std::vector<size_t> spans(shape.size());
        spans[shape.size()-1]=1u;
        for (int i=static_cast<int>(shape.size())-2; i>-1; --i)
          spans[i] = spans[i+1]*shape[i+1];
        // if the calculated spans and input strides are equivalent, we can
        // skip calculating indicies:
        bool roword = true;
        for (size_t i=0; i<strides.size(); ++i)
          roword &= strides[i]/sizeof(T) == spans[i];
        if (roword){
          for (size_t i=0; i<this->N*this->M; ++i) this->_data[i] = d[i];
        } else {
          // size_t tmp, lin, idx;
          // loop over all linear indicies
          // calculate the subscripted indices using our new row-ordered span
          // and then calculate the linear index into the input _data using
          // the provided strides -- remembering to convert from bytes to index
          for (size_t i=0; i<this->N*this->M; ++i){
            size_t tmp = i, lin = 0;
            for (size_t j=0; j<shape.size(); ++j){
              size_t idx = tmp/spans[j];
              tmp -= idx*spans[j];
              lin += idx*strides[j]/sizeof(T);
            }
            this->_data[i] = d[lin];
          }
        } // end if not row-ordered and contiguous
      } // end if input _data is not null and N*M>0
    } // end if shape and strides contain equal number of elements
  }
  /*! Type converting ArrayVector constructor
      @param m the number of elements within each array
      @param n the number of arrays in the ArrayVector
      @param d a pointer to a n*m block of _data which is type converted and copied into the ArrayVector
  */
  template<class R, typename=typename std::enable_if<std::is_convertible<R,T>::value>::type>
  ArrayVector(size_t m=0, size_t n=0, const R* d=nullptr): M(m), N(n), _data(nullptr){
    if (m && n) _data = new T[m*n]();
    if (d && m && n) for (size_t i=0; i<m*n; i++) _data[i] = T(d[i]);
  }
  /*! Copy constructor
      @param vec another ArrayVector which is copied into the new object
      @note this is used by objects which wrap the ArrayVector,
            e.g., the lattice vectors LQVec and LDVec
  */
  ArrayVector(const ArrayVector<T>& vec): M(vec.numel()), N(vec.size()), _data(nullptr){
    size_t m = vec.numel();
    size_t n = vec.size();
    if (m && n){
      T *d = vec.data(); // if m*n==0 data will throw a ValueError
      _data = new T[m*n]();
      if (d) for(size_t i=0; i<m*n; i++) _data[i] = d[i];
    }
  }
  //! Constructor from a standard vector of standard arrays
  template<class R, size_t Nel> ArrayVector(const std::vector<std::array<R,Nel>>& va):
  M(Nel), N(va.size()), _data(nullptr){
    if (M && N){
      _data = new T[M*N]();
      for (size_t i=0; i<N; ++i) for (size_t j=0; j<M; ++j) _data[i*M+j] = static_cast<T>(va[i][j]);
    }
  }
  //! Construct a single-array ArrayVector from a standard array
  template<class R, size_t Nel> ArrayVector(const std::array<R,Nel>& va):
  M(Nel), N(1u), _data(nullptr){
    if (M && N){
      _data = new T[M*N]();
      for (size_t j=0; j<M; ++j) _data[0*M+j] = static_cast<T>(va[j]);
    }
  }
  //! Type converting copy constructor
  template<class R, typename=typename std::enable_if<std::is_convertible<R,T>::value>::type>
  ArrayVector(const ArrayVector<R>& vec): M(vec.numel()), N(vec.size()), _data(nullptr){
    size_t m = vec.numel();
    size_t n = vec.size();
    if (m && n){
      R *d = vec.data(); // if m*n==0 data will throw a ValueError
      _data = new T[m*n]();
      if (d) for(size_t i=0; i<m*n; i++) _data[i] = static_cast<T>(d[i]);
    }
  }
  // Assignment operator
  ArrayVector<T>& operator=(const ArrayVector<T>& other){
    if ( this != &other ){ // avoid self-assignment
      size_t m = other.numel();
      size_t n = other.size();
      // reuse the _data-block if we can. otherwise, refresh/resize it
      if ( m !=this->numel() ) this->refresh(m,n);
      if ( n !=this->size()  ) this->resize(n);
      // copy-over the _data (if it exists)
      if (m && n)
        for(size_t i=0; i<m*n; i++)
          this->_data[i] = other._data[i];
    }
    return *this;
  }
  //! Accessor for single-arrays -- preforms a memory copy
  ArrayVector<T> operator[](const size_t i) const{
    bool isok = i < this->size();
    ArrayVector<T> out(this->numel(), isok ? 1u: 0u);
    if (isok){
      out.set(1, this->data(i));
    }
    return out;
  }
  //! Custom deconstructor to deallocate heap memory
  ~ArrayVector() { if (M && N) delete[] _data; };
  //! Return the number of arrays
  size_t size() const {return N;};
  //! Return the number of elements in each array
  size_t numel() const {return M;};
  //! Returns the pointer to the ith array's jth element -- with bounds checking
  T* checked_data(size_t i=0, size_t j=0) const;
  //! Return the pointer to the ith array's jth element
  T* data(size_t i=0, size_t j=0) const;
  //! Return the value of the ith array's jth element
  T getvalue(const size_t i=0, const size_t j=0) const;
  //! Return the ith single-array ArrayVector
  ArrayVector<T> extract(const size_t i=0) const ;
  //! Return the first `num` arrays of the ArrayVector
  ArrayVector<T> first(const size_t num) const;
  /*! Return a collection of arrays from the ArrayVector
    @param n the number of arrays to return
    @param i a pointer to the first index of the n arrays to return
    @returns An ArrayVector containing the n arrays
  */
  ArrayVector<T> extract(const size_t n, const size_t *i) const;
  /*! Return a collection of arrays from the ArrayVector
    @param idx a reference to an ArrayVector containing the to-be-returned indices
    @returns An ArrayVector containing the indicies indicated by idx
  */
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  ArrayVector<T> extract(const ArrayVector<I>& idx) const;
  /*! Return a collection of arrays from the ArrayVector
    @param idx a reference to a standard vector containing the to-be-returned indices
    @returns An ArrayVector containing the indicies indicated by idx
  */
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value> >
  ArrayVector<T> extract(const std::vector<I>& idx) const;
  /*! Return a collection of array from the ArrayVector
    @param idx a reference to a standard array containing the to-be-returend indices
    @returns an ArrayVector containing the indicies indicated by idx
  */
  template<typename I, typename=std::enable_if_t<std::is_integral<I>::value>, size_t Nel>
  ArrayVector<T> extract(const std::array<I,Nel>& idx) const;
  /*! Return a collection of arrays from the ArrayVector
    @param tfvec a reference to an ArrayVector<bool> with true for the to-be-returned indices
    @returns An ArrayVector containing the indicies indicated by tfvec
  */
  ArrayVector<T> extract(const ArrayVector<bool>& idx) const;
  /*! Return a collection of arrays from the ArrayVector
    @param tfvec a reference to an std::vector<bool> with true for the to-be-returned indices
    @returns An ArrayVector containing the indicies indicated by tfvec
  */
  ArrayVector<T> extract(const std::vector<bool>& idx) const;
  /*! Copy-out a single array
    @param i the index of the array to copy
    @param[out] out a pointer where the array will be copied
    @returns a flag to indicate success
  */
  bool get(const size_t i, T* out) const;
  /*! Copy-in a single array
    @param i the index of the array to copy into
    @param[in] in a pointer where the array will be copied from
    @returns a flag to indicate success
  */
  bool set(const size_t i, const T* in);
  /*! Copy-in a single array
    @param i the index of the array to copy into
    @param[in] in a pointer to an ArrayVector the array will be copied from
    @returns a flag to indicate success
  */
  bool set(const size_t i, const ArrayVector<T>* in);
  /*! Copy-in a single array
    @param i the index of the array to copy into
    @param[in] in a reference to an ArrayVector the array will be copied from
    @returns a flag to indicate success
  */
  bool set(const size_t i, const ArrayVector<T>& in);
  /*! Insert a new value at the specified index
    @param[in] in the new value to store
    @param i the index of the array to hold the new value
    @param j [optional] the index into the ith array where the value will be stored
  */
  bool set(const size_t i, const std::vector<T>& in);
  template<size_t Nel> bool set(const size_t i, const std::array<T,Nel>& in);
  bool insert(const T in, const size_t i, const size_t j=0u);
  /*! Print a subset of the arrays to console using the specified format string
    @param[in] fmt A char pointer to the format string used to format a single array element, e.g., "%g"
    @param first The index of the first array to print
    @param last The index of the last array to print
    @param[in] after [optional] A char pointer to what (if anything) should be printed after each array (default="\n")
  */
  void printformatted(const char * fmt, const size_t first, const size_t last, const char *after = "\n") const;
  /*! Print all arrays to console */
  void print() const;
  //! Print the ith array to console
  void print(const size_t i) const;
  /*! Print a subset of the arrays to console
    @param first The index of the first array to print
    @param last The index of the last array to print
    @param[in] after [optional] A char pointer to what (if anything) should be printed after each array (default="\n")
  */
  void print(const size_t first, const size_t last, const char *after="\n") const;
  void printheader(const char* name="ArrayVector") const;
  /*! Return a subset of the contents of the ArrayVector as a std::string
    @param first The index of the first array to convert
    @param last The index of the last array to convert
    @param after [optional] will be included in the string after each array (default="\n")
    @note This method performs no bounds checking. Use to_string to include bounds checking on first and last
  */
  std::string unsafe_to_string(const size_t first, const size_t last, const std::string &after="\n") const;
  std::string to_string() const;
  //! Return a std::string containing the ith array
  std::string to_string(const size_t i) const;
  //! Return a std::string containing all arrays with specified seperator
  std::string to_string(const std::string &) const;
  //! Return a std::string containing the ith array with specified following string
  std::string to_string(const size_t i, const std::string &) const;
  //! Return a std::string containing the ith array with the specified following character
  /* This version prevents to_string(i,'\n') from being converted to
     to_string(i,10u) implicitly, since '\n' is char(10). If i≠10 such a case
     would give unexpected results. One could ensure that to_string(i, "\n") was
     used instead, but it is probably unrealistic to do so.*/
  void to_string(const size_t, const char) const {
    throw std::runtime_error("Invalid use of ArrayVector::to_string. Check for, e.g., x.to_string(i,'\\n')");
  }
  /*! Return a subset of the contents of the ArrayVector as a std::string
    @param first The index of the first array to convert
    @param last The index of the last array to convert
    @param after [optional] will be included in the string after each array (default="\n")
  */
  std::string to_string(const size_t first, const size_t last, const std::string &after="\n") const;
  //! Print two ArrayVectors together:
  template<class R> std::string to_string(const ArrayVector<R>& other, const size_t num=0) const;
  /*! Modify the number of arrays that the ArrayVector can hold
    @param newsize the desired number of arrays to hold
    @note A new _data block is created even if the size does not change. As much
          of the old _data block as will fit is copied into the new block.
  */
  size_t resize(size_t newsize);
  /*! Modify the number of elements per array and the number of arrays
    @param newnumel The new number of elements per array
    @param newsize  The new number of arrays to hold
    @note A new _data block is created even if newnumel*newsize == numel*size.
          No _data is copied from the old to new block.
  */
  size_t refresh(size_t newnumel, size_t newsize=0u);
  //! Returns true if all elements evaluate to true.
  bool all_true(const size_t n=0) const;
  size_t count_true(const size_t n=0) const;
  size_t first_true(const size_t n=0) const;
  size_t last_true(const size_t n=0) const;
  //! Returns true if any elements evaluate to true.
  bool any_true(const size_t n=0) const;
  // //! Returns true if all elements are greater or equal to zero
  // bool all_positive(const size_t n=0) const;
  //! Returns true if all elements evaluate to false
  bool all_zero(const size_t n=0) const;
  /*! Returns true if all elements are approximately equal to the passed value
    @param val The comparison value, the type of which determines comparision precision
  */
  bool all_approx(const T val, const size_t n=0) const;
  bool none_approx(const T val, const size_t n=0) const;
  bool all_approx(const Comp expr, const T val, const size_t n=0) const;
  bool any_approx(const Comp expr, const T val, const size_t n=0) const;
  ArrayVector<bool> is_approx(const Comp expr, const T val, const size_t n=0) const;
  ArrayVector<bool> is_approx(const Comp expr, const std::vector<T>& val) const;
  bool vector_approx(size_t i, size_t j, Comp op=Comp::eq, T val=0.) const;
  template<class R, size_t Nel> bool rotate_approx(const size_t i, const size_t j, const std::array<R,Nel>&, const int order=1) const;
  //! Round all elements using std::round
  ArrayVector<int> round() const;
  //! Find the floor of all elements using std::floor
  ArrayVector<int> floor() const;
  //! Find the ceiling of all elements using std::ceil
  ArrayVector<int> ceil() const;
  /*! Find the sum over arrays or over elements
    @param dim Which dimension to sum over.
    @note dim==1 sums over the elements, returning a one-array ArrayVector
          dim!=1 sums over the arrays, returning a one-element ArrayVector
  */
  ArrayVector<T> sum( const int dim=0 ) const;
  ArrayVector<T> prod( const int dim=0 ) const;
  ArrayVector<T> min(const int dim=0) const;
  ArrayVector<T> max(const int dim=0) const;
  ArrayVector<bool> is_unique(void) const;
  ArrayVector<size_t> unique_idx(void) const;
  ArrayVector<T> unique(void) const;
  //! Ensure that a second ArrayVector object is consistent for binary operations
  template<typename R, typename=typename std::enable_if<std::is_convertible<T,R>::value||std::is_convertible<R,T>::value>::type>
  AVSizeInfo consistency_check(const ArrayVector<R>& b) const {
    const ArrayVector<T>& a = *this;
    AVSizeInfo si;
    si.oneveca = a.size() ==1u;
    si.scalara = a.numel()==1u;
    si.onevecb = b.size() ==1u;
    si.scalarb = b.numel()==1u;
    si.singular = si.scalarb && b.numel()!=a.numel(); // if both have numel==1 don't set the singular flag
    if (!(si.scalara^si.scalarb) && a.numel()!=b.numel()){
      throw std::runtime_error("binary operation(a,b) requires a.numel()==b.numel() or b.numel()==1");
    }
    si.n = si.oneveca ? b.size() : a.size();
    si.m = si.scalara? b.numel() : a.numel();
    if (si.oneveca^si.onevecb){
      si.aorb = !si.oneveca; // hopefully they're both scalars (or non scalars and the same)
    } else {
      si.aorb = !si.scalara; // in reality we need to make sure out gets *resized* to (m,n) no matter what
    }
    return si;
  }
  //! Ensure that a std::array object is consistent for binary operations
  template<typename R, size_t Nel, typename=typename std::enable_if<std::is_convertible<T,R>::value||std::is_convertible<R,T>::value>::type>
  AVSizeInfo consistency_check(const std::array<R,Nel>&) const {
    const ArrayVector<T>& a = *this;
    AVSizeInfo si;
    si.oneveca = a.size() ==1u;
    si.scalara = a.numel()==1u;
    if (Nel!=1u && (a.numel()!=Nel && a.size()!=Nel)){
      std::string err_msg = "std::array has size " + std::to_string(Nel);
      err_msg += " but needs to have " + std::to_string(a.size()) + " or ";
      err_msg += std::to_string(a.numel()) + " for binary operation with this";
      err_msg += " ArrayVector";
      throw std::runtime_error(err_msg);
    }
    si.onevecb = 1u==Nel || a.size()  != Nel; // if a.size() ==Nel then b is a single vector (and Nel==a.numel() or 1)
    si.scalarb = 1u==Nel || a.numel() != Nel; // if a.numel()==Nel then b is a scalar array (and Nel==a.size() or 1)
    // a std::array is either a scalar array or single vector. it is *never* multidimensional
    // if it happens that *this is square, the orientation of the std::array is arbitrary
    // typically it is a single vector, so default to that interpretation
    if (a.size() == a.numel() && a.size() == Nel) si.onevecb = true;
    size_t bsize =si.scalarb ? 1u : Nel;
    size_t bnumel=si.onevecb ? 1u : Nel;
    si.singular = si.scalarb && Nel!=a.numel(); // if both have numel==1 don't set the singular flag
    if (!(si.scalara^si.scalarb) && a.numel()!=Nel){
      throw std::runtime_error("binary operation(a,b) requires a.numel()==b.numel() or b.numel()==1");
    }
    si.n = si.oneveca ? bsize : a.size();
    si.m = si.scalara? bnumel : a.numel();
    if (si.oneveca^si.onevecb){
      si.aorb = !si.oneveca; // hopefully they're both scalars (or non scalars and the same)
    } else {
      si.aorb = !si.scalara; // in reality we need to make sure out gets *resized* to (m,n) no matter what
    }
    return si;
  }
  //! Ensure that a second ArrayVector object is consistent for in-place binary operations
  template<typename R, typename=typename std::enable_if<std::is_convertible<R,T>::value>::type>
  AVSizeInfo inplace_consistency_check(const ArrayVector<R> &b) const{
    const ArrayVector<T>& a = *this;
    AVSizeInfo si;
    si.oneveca = false;
    si.scalara = a.numel()==1u;
    si.onevecb = b.size()==1u;
    si.scalarb = b.numel()==1u;
    si.singular = si.scalarb && b.numel()!=a.numel(); // if both have numel==1 don't set the singular flag
    if (!si.scalarb && a.numel()!=b.numel()){
      throw std::runtime_error("binary operation(a,b) requires a.numel()==b.numel() or b.numel()==1");
    }
    si.n = a.size();
    si.m = a.numel();
    if (!si.onevecb && b.size()!=a.size()){
      throw std::runtime_error("equal sized or second-singular arrays required");
    }
    si.aorb = true; // this doesn't matter here but will later
    return si;
  }
  /*! Truncate the arrays by removing elements
    @param from the first element which will be removed from each array
    @param to the last element which will be removed from each array
    @note If to is greater than the number of elements in each array, the
          removal will stop at the array end.
    @returns 0 if the removal is successful,
             1 if from is greater than the number of starting elements, or
             2 if from=1 and to>=[number of initial elements]-1
    @note A non-zero return indicates that nothing was done to the _data block,
          otherwise a copy has been preformed.
  */
  int removeelements(const size_t from, const size_t to){
    if (from < this->M){ // otherwise there's nothing to do
       size_t last = to < this->M-1 ? to : this->M-1;
       size_t remaining_elements = from + (this->M-1 - last);
       if (remaining_elements != this->M){ // otherwise we have nothing to do
         T *newdata = nullptr;
         newdata = new T[remaining_elements*this->N]();
         size_t idx;
         for (size_t i=0; i<this->N; ++i)
           for (size_t j=0; j<this->M; ++j){
            idx = 0;
            if ( j < from || j > last)
              newdata[i*remaining_elements + idx++] = this->_data[i*this->N + j];
         }
         delete[] this->_data; //before we loose its pointer
         this->M = remaining_elements;
         this->_data = newdata;
         return 0;
       }
       return 2;
    }
    return 1;
  }
  /*! Add extra elements to each array
    @param ntoadd how many elements should be added
    @param valtoadd the value which each new element is set to (default=0)
  */
  int addelements(const size_t ntoadd, const T valtoadd=T(0)){
    size_t newM = this->M + ntoadd;
    size_t i, j;
    T* newdata = nullptr;
    newdata = new T[newM*this->N]();
    for (i=0; i<this->N; ++i){
      for (j=0; j<this->M; ++j) newdata[i*newM+j] = this->_data[i*this->M+j];
      for (j=this->M; j<newM; ++j)  newdata[i*newM+j] = valtoadd;
    }
    delete[] this->_data;
    this->M = newM;
    this->_data = newdata;
    return 0;
  }

  ArrayVector<T> operator -() const;

  ArrayVector<T>& operator +=(const ArrayVector<T>& plus);
  ArrayVector<T>& operator -=(const ArrayVector<T>& minus);
  ArrayVector<T>& operator *=(const ArrayVector<T>& times);
  ArrayVector<T>& operator /=(const ArrayVector<T>& divide);

  ArrayVector<T>& operator +=(const T& plus);
  ArrayVector<T>& operator -=(const T& minus);
  ArrayVector<T>& operator *=(const T& times);
  ArrayVector<T>& operator /=(const T& divide);

  //! Determine whether two ArrayVector objects contain approximately the same information
  template<typename R> bool isapprox(const ArrayVector<R> &that) const;
  //! Determine whether two arrays at indexes i and j contain approximately the same information
  bool isapprox(const size_t i, const size_t j) const;
  /*! If the ArrayVector has numel()==3, calculate the cross product of two arrays
    @param i the index of the first array
    @param j the index of the second array
    @param out a pointer where the resulting 3-elements are to be stored
  */
  void cross(const size_t i, const size_t j, T* out) const;
  //! Determine the dot product between two arrays at indexes i and j
  T dot(const size_t i, const size_t j) const;
  //! Determine the length (or its equivalent) of the array at index i
  T norm(const size_t i) const;

  void permute(const std::vector<size_t>& p);
  bool swap(const size_t i, const size_t j);
  bool swap(const size_t i, const size_t a, const size_t b);
  std::vector<T> to_std() const;

  template<class I> void permute_modes(const size_t, std::vector<I>&);
  template<class I> void inverse_permute_modes(const size_t, std::vector<I>&);
};

#include "arrayvector_operators.tpp"
#include "arrayvector.tpp"


#endif
