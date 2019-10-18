/*! \file */
#ifndef _SYMMETRY_H_
#define _SYMMETRY_H_

// #include <vector>
// #include <array>
#include <algorithm>
// #include <numeric>
#include "linear_algebra.h"

template<class T> using Matrix = std::array<T,9>;
template<class T> using Vector = std::array<T,3>;
template<class T> using Matrices = std::vector<Matrix<T>>;
template<class T> using Vectors = std::vector<Vector<T>>;

// class Motion{
//   Matrix<int> W;
//   Vector<double> w;
// public:
//   Motion(): W({1,0,0, 0,1,0, 0,0,1}), w({0.,0.,0.}) {}
//   Motion(const Matrix<int>& X): W(X), w({0.,0.,0.}) {}
//   Motion(const Vector<double>& x): W({1,0,0, 0,1,0, 0,0,1}), w(x) {}
//   Motion(const Matrix<int>& X, const Vector<double>& x): W(X), w(x) {}
//   Motion operator*(const Motion& m) const;
//   template<class T, typename S/*=promotion stuff*/>
//   Vector<S> operator*(const Vector<T>& x) const;
//   Matrix<int> rotation(void) const;
//   Vector<double> translation(void) const;
//   const Matrix<int>& rotation(void) const;
//   const Vector<double>& translation(void) const;
// }
// using Motions = std::vector<Motion>;

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
  Matrices<int> R;
  Vectors<double> T;
public:
  Symmetry(size_t n=0): R(n), T(n) { R.resize(n); T.resize(n); }
  const Matrices<int>&   getallr(void)               const { return this->R;  }
  const Vectors<double>& getallt(void)               const { return this->T;  }
  size_t                 size(void)                  const { return R.size(); }
  size_t                 resize(const size_t i)                                ;
  size_t                 add(const int *r, const double *t)                    ;
  size_t                 add(const Matrix<int>&, const Vector<double>&)        ;
  bool                   set(const size_t i, const int *r, const double *t)    ;
  bool                   get(const size_t i, int *r, double *t)           const;
  int *                  rdata(const size_t i)                                 ;
  double *               tdata(const size_t i)                                 ;
  const int *            rdata(const size_t i)                            const;
  const double *         tdata(const size_t i)                            const;
  Matrix<int>            getr(const size_t i)                             const;
  Vector<double>         gett(const size_t i)                             const;
  // const Matrix<int>&     getr(const size_t i)                             const;
  // const Vector<double>&  gett(const size_t i)                             const;
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
  Matrices<int> R;
public:
  PointSymmetry(size_t n=0): R(n) { R.resize(n);}
  PointSymmetry(const Matrices<int>& rots): R(rots){ this->sort(); }
  const Matrices<int>& getall(void)                  const { return this->R;  }
  size_t               size(void)                    const { return R.size(); }
  size_t               resize(const size_t newsize)                            ;
  size_t               add(const int *r)                                       ;
  size_t               add(const Matrix<int>&)                                 ;
  bool                 get(const size_t i, int *r)                        const;
  bool                 set(const size_t i, const int *r)                       ;
  int *                data(const size_t i)                                    ;
  const int *          data(const size_t i)                               const;
  Matrix<int>          pop(const size_t i=0)                                   ;
  size_t               erase(const size_t i)                                   ;
  bool                 has(const Matrix<int>&)                            const;
  Matrix<int>          get(const size_t i)                                const;
  Matrix<int>          get_proper(const size_t i)                         const;
  Matrix<int>          get_inverse(const size_t i)                        const;
  // const Matrix<int>&   get(const size_t i)                                const;
  void                 sort(const int ad=0)                                    ;
  void                 permute(const std::vector<size_t>&)                     ;
  int                  order(const size_t i)                              const;
  std::vector<int>     orders(void)                                       const;
  int                  isometry(const size_t i)                           const;
  std::vector<int>     isometries(void)                                   const;
  Vector<int>          axis(const size_t i)                               const;
  Vectors<int>         axes(void)                                         const;
  bool                 has_space_inversion(void)                          const;
  void                 print(const size_t i)                              const;
  void                 print(void)                                        const;
  PointSymmetry        generate(void)                                     const;
  PointSymmetry        generators(void)                                   const;
  PointSymmetry        nfolds(const int min_order=0)                      const;
  Vector<int>          perpendicular_axis(const size_t i)                 const;
  Vectors<int>         perpendicular_axes(void)                           const;
  PointSymmetry        higher(const int min_order=0)                      const;
};




#endif
