/*! \file */
#ifndef _SYMMETRY_H_
#define _SYMMETRY_H_

#include <vector>
#include <array>

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
  std::vector<std::array<int,9>> R;
  std::vector<std::array<double,3>> T;
public:
  Symmetry(size_t n=0): R(n), T(n) { R.resize(n); T.resize(n);};
  size_t               size() const { return R.size(); };
  size_t               add(const int *r, const double *t)                      ;
  size_t               add(const std::array<int,9>&, const std::array<double,3>&);
  bool                 set(const size_t i, const int *r, const double *t)      ;
  bool                 setrot(const size_t i, const int *r)                    ;
  bool                 settran(const size_t i, const double *t)                ;
  bool                 get(const size_t i, int *r, double *t)             const;
  bool                 getrot(const size_t i, int *r)                     const;
  bool                 gettran(const size_t i, double *t)                 const;
  int *                getrot(const size_t i)                                  ;
  double *             gettran(const size_t i)                                 ;
  const int *          getrot(const size_t i)                             const;
  const double *       gettran(const size_t i)                            const;
  std::array<int,9>    getrotarray(const size_t i)                        const;
  std::array<double,3> gettranarray(const size_t i)                       const;
  size_t               resize(const size_t i)                                  ;
  const std::vector<std::array<int,9>>& getallrots(void) const { return this->R; };
  const std::vector<std::array<double,3>>& getalltrans(void) const { return this->T; };
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
  std::vector<std::array<int,9>> R;
public:
  PointSymmetry(size_t n=0): R(n) { R.resize(n);};
  PointSymmetry(std::vector<std::array<int,9>>& rots): R(rots){ this->sort(); };
  size_t            size() const { return R.size(); };
  size_t            add(const int *r)                                          ;
  size_t            add(const std::array<int,9>&)                              ;
  bool              get(const size_t i, int *r)                           const;
  bool              set(const size_t i, const int *r)                          ;
  int *             get(const size_t i)                                        ;
  const int *       get(const size_t i)                                   const;
  std::array<int,9> getarray(const size_t i)                              const;
  size_t            resize(const size_t newsize)                               ;
  void              sort(const int ad=0)                                       ;
  const std::vector<std::array<int,9>>& getall(void) const { return this->R; };
  std::vector<int> orders(void) const;
  std::vector<int> isometries(void) const;
  std::vector<std::array<int,3>> axes(void) const;
  bool has_space_inversion() const;

  void print(const size_t i) const;
  void print(void) const;
};




#endif
