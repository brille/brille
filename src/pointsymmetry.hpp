/* This file is part of brille.

Copyright © 2019,2020 Greg Tucker <greg.tucker@stfc.ac.uk>

brille is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

brille is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with brille. If not, see <https://www.gnu.org/licenses/>.            */

/*! \file */

#ifndef BRILLE_POINTSYMMETRY_H_
#define BRILLE_POINTSYMMETRY_H_
#include <algorithm>
#include "utilities.hpp"
#include "symmetry_common.hpp" // defines Matrix, Vector, Matrices, Vectors
namespace brille {

// template<class T> using Matrix = std::array<T,9>;
// template<class T> using Vector = std::array<T,3>;
// template<class T> using Matrices = std::vector<Matrix<T>>;
// template<class T> using Vectors = std::vector<Vector<T>>;

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
  size_t               add(const std::string&)                                 ;
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
  size_t               get_inverse_index(const size_t i)                  const;
  size_t               find_index(const Matrix<int>&)                     const;
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

} // namespace brille
#endif
