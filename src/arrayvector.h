#ifndef _ARRAYVECTOR_H_
#define _ARRAYVECTOR_H_

#include<iostream>
#include<string>
#include<math.h>
#include<cmath>
#include<functional>
#include "safealloc.h"

template<typename T> class ArrayVector;
class LatVec; // forward declare so that we can prevent operator overloads from applying to LQVec and LDVec

typedef struct{
	size_t n; // the common number of Vectors
	size_t m; // the common number of elements per Vector
	bool oneveca; // is "a" a "SingleVector"
	bool onevecb; // is "b" a "SingleVector"
	bool scalara; // is "a" an "ArrayScalar"
	bool scalarb; // is "b" an "ArrayScalar"
	bool singular;// is only "b" an "ArrayScalar"
	bool aorb; // does "a" hold more Vector(/Scalar) elements
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
	size_t M;
	size_t N;
	T* data;
public:
	// Constructor
	ArrayVector(size_t m=0, size_t n=0, const T* d=nullptr) : M(m), N(n){
	    if (m*n) data = safealloc<T>(m*n);
	    if (d && m*n) for(size_t i=0; i<m*n; i++) data[i] = d[i];
	};
	// type converting constructor
	template<class R, typename=typename std::enable_if<std::is_convertible<R,T>::value>::type>
	ArrayVector(size_t m=0, size_t n=0, const R* d=nullptr): M(m), N(n){
		if (m*n) data = safealloc<T>(m*n);
		if (d && m*n) for (size_t i=0; i<m*n; i++) data[i] = T(d[i]);
	}
	// Copy-constructor (used by LDVec and LQVec constructors)
	ArrayVector(const ArrayVector<T>& vec): M(vec.numel()), N(vec.size()), data(nullptr){
		size_t m = vec.numel();
		size_t n = vec.size();
		T *d = vec.datapointer();
		if (m*n) data = safealloc<T>(m*n);
		if (d && m*n) for(size_t i=0; i<m*n; i++) data[i] = d[i];
	};
	// pseudo-copy constructors
	template<class R, typename=typename std::enable_if<std::is_convertible<R,T>::value>::type>
	ArrayVector(const ArrayVector<R>& vec): M(vec.numel()), N(vec.size()), data(nullptr){
		size_t m = vec.numel();
		size_t n = vec.size();
		R *d = vec.datapointer();
		if (m*n) data = safealloc<T>(m*n);
		if (d && m*n) for(size_t i=0; i<m*n; i++) data[i] = d[i];
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
			if (other.data && m*n)
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
	// ArrayVector<T>& operator[](const size_t i){
	// 	return *(this+i); // of course, this doesn't make sense.
	// }

	// Destructor
	~ArrayVector() { if (M*N) delete[] data; };
	// the number of arrays
	size_t size() const {return N;};
	// Number of elements in each array
	size_t numel() const {return M;};
	// get the pointer to the ith array's jth element
	T* datapointer(size_t i=0, size_t j=0) const;
	// get the value to the ith array's jth element
	T getvalue(const size_t i=0, const size_t j=0) const;
	// return a single element (but still an ArrayVector)
	ArrayVector<T> extract(const size_t i=0) const ;
	ArrayVector<T> extract(const size_t n, const size_t *i) const;
	ArrayVector<T> extract(const ArrayVector<size_t>& idx) const;
	// copy-out the ith array
	bool get(const size_t i, T* out) const;
	// copy-in the ith array
	bool set(const size_t i, const T* in);
	bool set(const size_t i, const ArrayVector<T>* in);
	bool set(const size_t i, const ArrayVector<T>& in);
	bool insert(const T in, const size_t i, const size_t j=0u);
	// print the array value to the console
	void printformatted(const char * fmt, const size_t first, const size_t last, const char *after = "\n") const;
	void print() const;
	void print(const size_t i) const;
	void print(const size_t first, const size_t last, const char *after="\n") const;
	void printheader(const char* name="ArrayVector") const;
	// Allow for initializing without knowing how big it needs to be:
	size_t resize(size_t newsize);
	// also we sometimes want to create a pointer with known numel and size
	size_t refresh(size_t newnumel, size_t newsize=0u);
	// treat elements as logical values:
	bool arealltrue(void) const;
	bool areanytrue(void) const;
	bool areallpositive(void) const;
	bool areallzero(void) const;
	ArrayVector<int> round() const;
	ArrayVector<int> floor() const;
	ArrayVector<int> ceil() const;
	ArrayVector<T> sum( const int dim=0 ) const;
	//operator overloading:
	template<typename R, typename=typename std::enable_if<std::is_convertible<T,R>::value||std::is_convertible<R,T>::value>::type>
	AVSizeInfo consistency_check(const ArrayVector<R>& b) const {
		const ArrayVector<T>& a = *this;
		AVSizeInfo si;
		si.oneveca = a.size() ==1u;
		si.scalara = a.numel()==1u;
		si.onevecb = b.size() ==1u;
		si.scalarb = b.numel()==1u;
		si.singular = si.scalarb && b.numel()!=a.numel(); // if both have numel==1 don't set the singular flag
		if (!(si.scalara^si.scalarb) && a.numel()!=b.numel()) throw std::runtime_error("binary operation(a,b) requires a.numel()==b.numel() or b.numel()==1");
		si.n = si.oneveca ? b.size() : a.size();
		si.m = si.scalara? b.numel() : a.numel();
		if (si.oneveca^si.onevecb){
			si.aorb = !si.oneveca; // hopefully they're both scalars (or non scalars and the same)
		} else {
			si.aorb = !si.scalara; // in reality we need to make sure out gets *resized* to (m,n) no matter what
		}
		return si;
	}
	template<typename R, typename=typename std::enable_if<std::is_convertible<R,T>::value>::type>
	AVSizeInfo inplace_consistency_check(const ArrayVector<R> &b) const{
		const ArrayVector<T>& a = *this;
		AVSizeInfo si;
		si.oneveca = false;
		si.scalara = a.numel()==1u;
		si.onevecb = b.size()==1u;
		si.scalarb = b.numel()==1u;
		si.singular = si.scalarb && b.numel()!=a.numel(); // if both have numel==1 don't set the singular flag
		if (!si.scalarb && a.numel()!=b.numel()) throw std::runtime_error("binary operation(a,b) requires a.numel()==b.numel() or b.numel()==1");
		si.n = a.size();
		si.m = a.numel();
		if (!si.onevecb && b.size()!=a.size()) throw std::runtime_error("equal sized or second-singular arrays required");
		si.aorb = true; // this doesn't matter here but will later
		return si;
	};

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
			 }
		}
	};
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

	template<typename R> bool isapprox(const ArrayVector<R> &that) const;
	bool isapprox(const size_t i, const size_t j) const;
	void cross(const size_t i, const size_t j, T* out) const;
	T dot(const size_t i, const size_t j) const;
	T norm(const size_t i) const;

};

#include "linear_algebra.h"
#include "arrayvector.hpp"


#endif
