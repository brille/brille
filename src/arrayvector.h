/*! \file */
#ifndef _ARRAYVECTOR_H_
#define _ARRAYVECTOR_H_

#include<iostream>
#include<string>
#include<cmath>
#include<functional>
#include<vector>
#include "linear_algebra.h"
#include "debug.h" // ensurses __PRETTY_FUNCTION__ is defined for MSVC, provides status_update()

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
 // Replace data by a std::vector< std::array<T,M> > ... then ArrayVector<T,M> and the standard libraries define most things
template<typename T> class ArrayVector{
protected:
  size_t M; //!< The number of elements within each array
  size_t N; //!< The number of arrays in the ArrayVector
  T* data;  //!< A pointer to the first element of the contiguous memory data block
public:
  /*! Standard ArrayVector constructor
      @param m the number of elements within each array
      @param n the number of arrays in the ArrayVector
      @param d [optional] a pointer to a n*m block of data which is copied into the ArrayVector
  */
  ArrayVector(size_t m=0, size_t n=0, const T* d=nullptr) : M(m), N(n){
      if (m && n) data = safealloc<T>(m*n);
      if (d && m && n) for(size_t i=0; i<m*n; i++) data[i] = d[i];
  };
  /*! ArrayVector contstructor taking shape and stride information for the case
      of a non-contiguous and/or non-row-ordered array at the provided pointer.
      @param d a pointer to the n*m block of data, to be copied
      @param shape a vector of sizes along each dimension of the input array
      @param strides a vector of strides along each dimension of the input array
      @note The strides vector must indicate the number of bytes necessary to
            move from one point to the next along a given dimension of the input
            data array. The ArrayVector is strictly two dimensional and
            therefore a convention must be adopted for converting the
            D-dimensional input array [D=s.size()].
            This constructor will attempt to combine the 2ⁿᵈ through Dᵗʰ
            dimensions such that (n,m) = (s[0],s[1]*…*s[D-1])
  */
  template<class Integer, typename=typename std::enable_if<std::is_integral<Integer>::value>::type>
  ArrayVector(const T* d, const std::vector<Integer>& shape, const std::vector<Integer>& strides){
    if (shape.size()>0 && shape.size()==strides.size()){
      this->N = static_cast<size_t>(shape[0]);
      Integer nel=1; for (size_t i=1; i<shape.size(); ++i) nel *= shape[i];
      this->M = static_cast<size_t>(nel);
      if (this->N && this->M) this->data = safealloc<T>(this->N*this->M);
      if (d && this->N && this->M){
        // we want to copy-in the data to a row-ordered array (and flatten it)
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
          for (size_t i=0; i<this->N*this->M; ++i) this->data[i] = d[i];
        } else {
          size_t tmp, lin, idx;
          // loop over all linear indicies
          // calculate the subscripted indices using our new row-ordered span
          // and then calculate the linear index into the input data using
          // the provided strides -- remembering to convert from bytes to index
          for (size_t i=0; i<this->N*this->M; ++i){
            tmp = i;
            lin = 0;
            for (size_t j=0; j<shape.size(); ++j){
              idx  = tmp/spans[j];
              tmp -= idx*spans[j];
              lin += idx*strides[j]/sizeof(T);
            }
            this->data[i] = d[lin];
          }
        } // end if not row-ordered and contiguous
      } // end if input data is not null and N*M>0
    } // end if shape and strides contain equal number of elements
  };
  /*! Type converting ArrayVector constructor
      @param m the number of elements within each array
      @param n the number of arrays in the ArrayVector
      @param d a pointer to a n*m block of data which is type converted and copied into the ArrayVector
  */
  template<class R, typename=typename std::enable_if<std::is_convertible<R,T>::value>::type>
  ArrayVector(size_t m=0, size_t n=0, const R* d=nullptr): M(m), N(n){
    if (m && n) data = safealloc<T>(m*n);
    if (d && m && n) for (size_t i=0; i<m*n; i++) data[i] = T(d[i]);
  }
  /*! Copy constructor
      @param vec another ArrayVector which is copied into the new object
      @note this is used by objects which wrap the ArrayVector, e.g., the lattice vectors
  */
  ArrayVector(const ArrayVector<T>& vec): M(vec.numel()), N(vec.size()), data(nullptr){
    size_t m = vec.numel();
    size_t n = vec.size();
    if (m && n){
      T *d = vec.datapointer(); // if m*n==0 datapointer will throw a ValueError
      data = safealloc<T>(m*n);
      if (d) for(size_t i=0; i<m*n; i++) data[i] = d[i];
    }
  };
  //! Constructor from std::vector<std::array<T,N>>
  template<class R, size_t Nel> ArrayVector(const std::vector<std::array<R,Nel>>& va):
  M(Nel), N(va.size()), data(nullptr){
    if (M && N){
      data = safealloc<T>(M*N);
      for (size_t i=0; i<N; ++i) for (size_t j=0; j<M; ++j) data[i*M+j] = static_cast<T>(va[i][j]);
    }
  }
  //! Type converting copy constructor
  template<class R, typename=typename std::enable_if<std::is_convertible<R,T>::value>::type>
  ArrayVector(const ArrayVector<R>& vec): M(vec.numel()), N(vec.size()), data(nullptr){
    size_t m = vec.numel();
    size_t n = vec.size();
    if (m && n){
      R *d = vec.datapointer(); // if m*n==0 datapointer will throw a ValueError
      data = safealloc<T>(m*n);
      if (d) for(size_t i=0; i<m*n; i++) data[i] = static_cast<T>(d[i]);
    }
  };
  // Assignment operator
  ArrayVector<T>& operator=(const ArrayVector<T>& other){
    if ( this != &other ){ // avoid self-assignment
      size_t m = other.numel();
      size_t n = other.size();
      // reuse the data-block if we can. otherwise, refresh/resize it
      if ( m !=this->numel() ) this->refresh(m,n);
      if ( n !=this->size()  ) this->resize(n);
      // copy-over the data (if it exists)
      if (other.data && m && n)
        for(size_t i=0; i<m*n; i++)
          this->data[i] = other.data[i];
    }
    return *this;
  };
  ArrayVector<T> operator[](const size_t i) const{
    bool isok = i < this->size();
    ArrayVector<T> out(this->numel(), isok ? 1u: 0u);
    if (isok){
      out.set(1, this->datapointer(i));
    }
    return out;
  }
  // Destructor
  ~ArrayVector() { if (M && N) delete[] data; };
  //! Returns the number of arrays
  size_t size() const {return N;};
  //! Returns the number of elements in each array
  size_t numel() const {return M;};
  //! Returns the pointer to the ith array's jth element
  T* datapointer(size_t i=0, size_t j=0) const;
  //! Returns the value of the ith array's jth element
  T getvalue(const size_t i=0, const size_t j=0) const;
  //! Return the ith single-array ArrayVector
  ArrayVector<T> extract(const size_t i=0) const ;
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
  ArrayVector<T> extract(const ArrayVector<size_t>& idx) const;
  /*! Return a collection of arrays from the ArrayVector
    @param idx a reference to a standard vector containing the to-be-returned indices
    @returns An ArrayVector containing the indicies indicated by idx
  */
  ArrayVector<T> extract(const std::vector<size_t>& idx) const;
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
  /*! Return a subset of the contents of the ArrayVector as a std::string
    @param first The index of the first array to convert
    @param last The index of the last array to convert
    @param after [optional] will be included in the string after each array (default="\n")
  */
  std::string to_string(const size_t first, const size_t last, const std::string &after="\n") const;
  /*! Modify the number of arrays that the ArrayVector can hold
    @param newsize the desired number of arrays to hold
    @note A new data block is created even if the size does not change. As much
          of the old data block as will fit is copied into the new block.
  */
  size_t resize(size_t newsize);
  /*! Modify the number of elements per array and the number of arrays
    @param newnumel The new number of elements per array
    @param newsize  The new number of arrays to hold
    @note A new data block is created even if newnumel*newsize == numel*size.
          No data is copied from the old to new block.
  */
  size_t refresh(size_t newnumel, size_t newsize=0u);
  //! Returns true if all elements evaluate to true.
  bool all_true(const size_t n=0) const;
  size_t count_true(const size_t n=0) const;
  //! Returns true if any elements evaluate to true.
  bool any_true(const size_t n=0) const;
  //! Returns true if all elements are greater or equal to zero
  bool all_positive(const size_t n=0) const;
  //! Returns true if all elements evaluate to false
  bool all_zero(const size_t n=0) const;
  /*! Returns true if all elements are approximately equal to the passed value
    @param val The comparison value, the type of which determines comparision precision
  */
  bool all_approx(const T val, const size_t n=0) const;
  bool none_approx(const T val, const size_t n=0) const;
  bool all_approx(const std::string& expr, const T val, const size_t n=0) const;
  bool any_approx(const std::string& expr, const T val, const size_t n=0) const;
  ArrayVector<bool> is_approx(const std::string& expr, const T val, const size_t n=0) const;
  bool vector_approx(const size_t i, const size_t j, const std::string& op="", const T val=0.) const;
  template<class R, size_t Nel> bool rotate_approx(const size_t i, const size_t j, const std::array<R,Nel>&) const;
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
  };
  /*! Truncate the arrays by removing elements
    @param from the first element which will be removed from each array
    @param to the last element which will be removed from each array
    @note If to is greater than the number of elements in each array, the
          removal will stop at the array end.
    @returns 0 if the removal is successful,
             1 if from is greater than the number of starting elements, or
             2 if from=1 and to>=[number of initial elements]-1
    @note A non-zero return indicates that nothing was done to the data block,
          otherwise a copy has been preformed.
  */
  int removeelements(const size_t from, const size_t to){
    size_t last;
    size_t remaining_elements;
    if (from < this->M){ // otherwise there's nothing to do
       last = to < this->M-1 ? to : this->M-1;
       remaining_elements = from + (this->M-1 - last);
       if (remaining_elements != this->M){ // otherwise we have nothing to do
         T *newdata = new T[remaining_elements*this->N]();
         size_t idx;
         for (size_t i=0; i<this->N; ++i)
           for (size_t j=0; j<this->M; ++j){
            idx = 0;
            if ( j < from || j > last)
              newdata[i*remaining_elements + idx++] = this->data[i*this->N + j];
         }
         delete[] this->data; //before we loose its pointer
         this->M = remaining_elements;
         this->data = newdata;
         return 0;
       }
       return 2;
    }
    return 1;
  };
  /*! Add extra elements to each array
    @param ntoadd how many elements should be added
    @param valtoadd the value which each new element is set to (default=0)
  */
  int addelements(const size_t ntoadd, const T valtoadd=T(0)){
    size_t newM = this->M + ntoadd;
    size_t i, j;
    T* newdata = new T[newM*this->N]();
    for (i=0; i<this->N; ++i){
      for (j=0; j<this->M; ++j) newdata[i*newM+j] = this->data[i*this->M+j];
      for (j=this->M; j<newM; ++j)  newdata[i*newM+j] = valtoadd;
    }
    delete[] this->data;
    this->M = newM;
    this->data = newdata;
    return 0;
  };

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

};


#include "arrayvector.hpp"


#endif
