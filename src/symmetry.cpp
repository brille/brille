#include "symmetry.h"

/*****************************************************************************\
| Symmetry class Member functions:                                            |
|-----------------------------------------------------------------------------|
|   set, setrot, settran   copy a given matrix and/or vector into R and/or T  |
|                          for the motion i.                                  |
|   get, getrot, gettran   copy the rotation and/or translation of motion i   |
|                          into the provided matrix and/or vector.            |
|   getrot, gettran        return a pointer to the first element of the       |
|                          rotation or translation of motion i.               |
|   resize                 grow or shrink the number of motions that the      |
|                          object can hold -- this necessitates a memory copy |
|                          if the old and new sizes are non-zero.             |
|   size                   return the number of motions the object can store. |
\*****************************************************************************/
bool Symmetry::set(const unsigned int i, const int *r, const double *t){
	if ( i>=N ) return false;
	for(unsigned int j=0; j<9; j++) this->R[i*N+j] = r[j];
	for(unsigned int j=0; j<3; j++) this->T[i*N+j] = t[j];
	return true;
}
bool Symmetry::setrot(const unsigned int i, const int *r){
	if ( i>=N ) return false;
	for(unsigned int j=0; j<9; j++) this->R[i*N+j] = r[j];
	return true;
}
bool Symmetry::settran(const unsigned int i, const double *t){
	if ( i>=N ) return false;
	for(unsigned int j=0; j<3; j++) this->T[i*N+j] = t[j];
	return true;
}
bool Symmetry::get(const unsigned int i, int *r, double *t){
	if ( i>=N ) return false;
	for(unsigned int j=0; j<9; j++) r[j]=this->R[i*N+j];
	for(unsigned int j=0; j<3; j++) t[j]=this->T[i*N+j];
	return true;
}
bool Symmetry::getrot(const unsigned int i, int *r){
	if ( i>=N ) return false;
	for(unsigned int j=0; j<9; j++) r[j]=this->R[i*N+j];
	return true;
}
bool Symmetry::gettran(const unsigned int i, double *t){
	if ( i>=N ) return false;
	for(unsigned int j=0; j<3; j++) t[j]=this->T[i*N+j];
	return true;
}
int    * Symmetry::getrot (const unsigned int i){ return (i<N) ? this->R + i*9 : nullptr; }
double * Symmetry::gettran(const unsigned int i){ return (i<N) ? this->T + i*3 : nullptr; }
unsigned int Symmetry::resize(const unsigned int newsize){
	bool std = newsize>0;
	int    *newR = nullptr;
	double *newT = nullptr;
	// allocate a new block of memory
	if (std) {
		newR = new int [newsize*9]();
		newT = new double [newsize*3]();
	}
	if (this->N) { // copy-over data :(
		unsigned int smallerN = this->N < newsize ? this->N : newsize;
		for (unsigned int i=0; i<smallerN*9; i++) newR[i] = this->R[i];
		for (unsigned int i=0; i<smallerN*3; i++) newT[i] = this->T[i];
		// hand-back the chunk of memory which data points to
		delete[] this->R;
		delete[] this->T;
	}
	// and set data to the newdata pointer;
	this->N = newsize;
	if (std) {
		this->R = newR;
		this->T = newT;
	}
	return newsize;
}

/*****************************************************************************\
| PointSymmetry class Member functions:                                       |
|-----------------------------------------------------------------------------|
|   set      copy a given matrix into R at index i.                           |
|   get      copy the rotation at index i into the provided matrix.           |
|   resize   grow or shrink the number of rotations that the object can hold  |
|            -- this causes a memory copy if the old and new sizes are finite.|
|   size     return the number of rotations the object can/does store.        |
\*****************************************************************************/
bool  PointSymmetry::set(const unsigned int i, const int *r){
	if ( i>=N ) return false;
	for(unsigned int j=0; j<9; j++) this->R[i*N+j] = r[j];
	return true;
}
bool  PointSymmetry::get(const unsigned int i, int *r){
	if ( i>=N ) return false;
	for(unsigned int j=0; j<9; j++) r[j]=this->R[i*N+j];
	return true;
}
int * PointSymmetry::get(const unsigned int i){ return (i<N) ? this->R+i*9 : nullptr; }
unsigned int PointSymmetry::resize(const unsigned int newsize){
	bool std = newsize>0;
	int    *newR = nullptr;
	// allocate a new block of memory
	if (std) {
		newR = new int [newsize*9]();
	}
	if (this->N) { // copy-over data :(
		unsigned int smallerN = this->N < newsize ? this->N : newsize;
		for (unsigned int i=0; i<smallerN*9; i++) newR[i] = this->R[i];
		// hand-back the chunk of memory which data points to
		delete[] this->R;
	}
	// and set data to the newdata pointer;
	this->N = newsize;
	if (std) {
		this->R = newR;
	}
	return newsize;
}
