/*! \file */
#ifndef _SYMMETRY_H_
#define _SYMMETRY_H_

/*****************************************************************************\
| Symmetry class                                                              |
|-----------------------------------------------------------------------------|
|   Holds N 3x3 matrices R and N 3 vectors T which are the motions of a       |
|   space group. A motion is the combination of a rotation and a translation. |
|-----------------------------------------------------------------------------|
| Member functions:                                                           |
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
class Symmetry{
  unsigned int N;
  int         *R;
  double      *T;
public:
  Symmetry(unsigned int n=0): N(n), R(nullptr), T(nullptr) { if (N){ R = new int [N*9](); T = new double [N*3](); } };
  ~Symmetry(){ if (N){ delete[] R; delete[] T; } };
  unsigned int size() { return N; };
  bool     set(const unsigned int i, const int *r, const double *t);
  bool     setrot(const unsigned int i, const int *r);
  bool     settran(const unsigned int i, const double *t);
  bool     get(const unsigned int i, int *r, double *t);
  bool     getrot(const unsigned int i, int *r);
  bool     gettran(const unsigned int i, double *t);
  int    * getrot(const unsigned int i);
  double * gettran(const unsigned int i);
  unsigned int resize(const unsigned int i);
};


/*****************************************************************************\
| PointSymmetry class                                                         |
|-----------------------------------------------------------------------------|
|   Holds N 3x3 rotation matrices R which comprise a point group symmetry.    |
|-----------------------------------------------------------------------------|
| Member functions:                                                           |
|   set      copy a given matrix into R at index i.                           |
|   get      copy the rotation at index i into the provided matrix.           |
|   resize   grow or shrink the number of rotations that the object can hold  |
|            -- this causes a memory copy if the old and new sizes are finite.|
|   size     return the number of rotations the object can/does store.        |
\*****************************************************************************/
class PointSymmetry{
  unsigned int N;
  int         *R;
public:
  PointSymmetry(unsigned int n=0): N(n), R(nullptr) {if (N) R = new int[N*9](); };
  ~PointSymmetry(){ if (N) delete[] R; };
  unsigned int size(){ return N; };
  bool  get(const unsigned int i, int *r);
  bool  set(const unsigned int i, const int *r);
  int * get(const unsigned int i);
  unsigned int resize(const unsigned int newsize);
};




#endif
