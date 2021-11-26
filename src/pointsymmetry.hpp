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
/*! \file
    \author Greg Tucker
    \brief Classes for a lattice pointgroup symmetry operations
*/
#include "symmetry_common.hpp" // defines Matrix, Vector, Matrices, Vectors
#include "hdf_interface.hpp"
namespace brille {
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
/*! \brief Holds N 3x3 rotation matrices R which comprise a point group symmetry. */
class PointSymmetry{
  Matrices<int> R;
public:
  explicit PointSymmetry(size_t n=0): R(n) { R.resize(n);}
  explicit PointSymmetry(const Matrices<int>& rots): R(rots){ this->sort(); }
  [[nodiscard]] const Matrices<int>& getall()                  const { return this->R;  }
  [[nodiscard]] size_t               size()                    const { return R.size(); }
  size_t               resize(size_t newsize)                            ;
  size_t               add(const int *r)                                       ;
  size_t               add(const Matrix<int>&)                                 ;
  size_t               add(const std::string&)                                 ;
  bool                 get(size_t i, int *r)                        const;
  bool                 set(size_t i, const int *r)                       ;
  int *                data(size_t i)                                    ;
  [[nodiscard]] const int *          data(size_t i)                               const;
  Matrix<int>          pop(size_t i=0)                                   ;
  size_t               erase(size_t i)                                   ;
  [[nodiscard]] bool                 has(const Matrix<int>&)                            const;
  [[nodiscard]] Matrix<int>          get(size_t i)                                const;
  [[nodiscard]] Matrix<int>          get_proper(size_t i)                         const;
  [[nodiscard]] Matrix<int>          get_inverse(size_t i)                        const;
  [[nodiscard]] size_t               get_inverse_index(size_t i)                  const;
  [[nodiscard]] size_t               find_index(const Matrix<int>&)               const;
  [[nodiscard]] size_t               find_identity_index()                        const;
  // const Matrix<int>&   get(const size_t i)                                const;
  void                 sort(int ad=0)                                    ;
  void                 permute(const std::vector<size_t>&)                     ;
  [[nodiscard]] int                  order(size_t i)                              const;
  [[nodiscard]] std::vector<int>     orders()                                       const;
  [[nodiscard]] int                  isometry(size_t i)                           const;
  [[nodiscard]] std::vector<int>     isometries()                                   const;
  [[nodiscard]] Vector<int>          axis(size_t i)                               const;
  [[nodiscard]] Vectors<int>         axes()                                         const;
  [[nodiscard]] bool                 has_space_inversion()                          const;
  void                 print(size_t i)                              const;
  void                 print()                                        const;
  [[nodiscard]] PointSymmetry        generate()                                     const;
  [[nodiscard]] PointSymmetry        generators()                                   const;
  [[nodiscard]] PointSymmetry        nfolds(int min_order=0)                      const;
  [[nodiscard]] Vector<int>          perpendicular_axis(size_t i)                 const;
  [[nodiscard]] Vectors<int>         perpendicular_axes()                           const;
  [[nodiscard]] PointSymmetry        higher(int min_order=0)                      const;
#ifdef USE_HIGHFIVE
    // Output to HDF5 file/object
    template<class HF>
    std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, bool>
    to_hdf(HF& obj, const std::string& entry) const{
        std::vector<HF_Matrix<int>> hfm;
        for (const auto & x: R) hfm.emplace_back(x);
        auto ds = overwrite_data(obj, entry, hfm);
        return true;
    }
    // Input from HDF5 file/object
    template<class HF>
    static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, PointSymmetry>
    from_hdf(HF& obj, const std::string& entry){
        std::vector<HF_Matrix<int>> hfm;
        obj.getDataSet(entry).read(hfm);
        Matrices<int> m;
        for (const auto& x: hfm) m.push_back(x.array());
        return PointSymmetry(m);
    }
#endif
};

} // namespace brille
#endif
