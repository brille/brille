#ifndef _MESH_H_
#define _MESH_H_

#include "latvec.h"
#include <vector>
#include <array>
#include <omp.h>
#include "interpolation.h"

#include "triangulation.h" // defined Delaunay
// typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>     Vertex_Base;
// typedef CGAL::Delaunay_triangulation_cell_base_3<K>                  Cell_Base;
// typedef CGAL::Triangulation_data_structure_3<Vertex_Base, Cell_Base> Tds;
// typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>  Delaunay;


template<class T> class Mesh3{
protected:
  TetrahedralTriangulation mesh;
  ArrayVector<T> data;         //!< The stored ArrayVector indexed by `map`
  ArrayVector<size_t> shape;   //!< A second ArrayVector to indicate a possible higher-dimensional shape of each `data` array
  std::array<unsigned,4> elements; //! The number of [scalar, normalized eigenvector, vector, matrix] elements per data array
  size_t branches;             //!< The number of branches contained per data array
public:
  //constructors
  Mesh3(const ArrayVector<double>& verts,
        const std::vector<std::vector<int>>& facets,
        const double max_volume,
        const double min_angle=-1.0,
        const double max_angle=-1.0,
        const double min_ratio=-1.0,
        const int max_points=-1
      ):
        data(0,0),
        shape(1,0),
        elements({0,0,0,0}),
        branches(0) {
    this->mesh = triangulate(verts, facets, max_volume, min_angle, max_angle, min_ratio, max_points);
  }
  Mesh3(const Mesh3<T>& other){
    this->mesh = other.mesh;
    this->data = other.data;
    this->shape= other.shape;
    this->elements = other.elements;
    this->branches = other.branches;
  }
  Mesh3<T>& operator=(const Mesh3<T>& other){
    this->mesh = other.mesh;
    this->data = other.data;
    this->shape= other.shape;
    this->elements = other.elements;
    this->branches = other.branches;
    return *this;
  }
  //! Return the number of mesh vertices
  size_t size() const { return this->mesh.number_of_vertices(); }
  //! Return the positions of all vertices in the mesh
  const ArrayVector<double>& get_mesh_xyz() const{ return this->mesh.get_vertex_positions(); }
  //! Return the tetrahedron indices of the mesh
  const ArrayVector<size_t>& get_mesh_tetrehedra() const{ return this->mesh.get_vertices_per_tetrahedron();}
  //! Return a constant reference to the data ArrayVector
  const ArrayVector<T>& get_data() const { return this->data; }
  //! Return the number of dimensions of each array in the `data` ArrayVector
  size_t data_ndim(void) const {return this->shape.size();}
  //! Return the number of arrays in the `data` ArrayVector
  size_t num_data(void) const {return this->data.size();}
  //! Return the `shape` ArrayVector containing information about the arrays in the `data` ArrayVector
  ArrayVector<size_t> data_shape(void) const {return this->shape;}
  /*! \brief Sum over the data array

  Either add together the elements of the array stored at each mapped point
  or add together all of the arrays.
  @param axis The axis along which to perform the summation -- 1 adds arrays,
  0 (or not 1, really) adds elements of each array.
  @returns An ArrayVector with `numel()==1` for `axis=0` or `size()==1` for
  `axis=1`
  */
  ArrayVector<T> sum_data(const int axis) const{return this->data.sum(axis);}

  /*! Replace the data stored in the object
  @param newdata the new ArrayVector of data to be stored
  @param newshape the shape information of each array in `newdata`
  @param new_elements The number of scalar, eigenvector elements, vector elements, and matrix elements contained each newdata array
  @note `newshape` provides information about the N dimensions and extent of
        each array in `newdata`, the following four unsigned integers provide
        information on how many scalar, eigenvector, vector, and matrix elements
        (in that order) are contained within the N-1 dimensions of each array.
        `new_elements[0]`+`new_elements[1]`+`new_elements[2]`+
        `new_elements[3]`×`new_elements[3]` must be equal to the
        product of newshape[1:N-1].
  */
  int replace_data(const ArrayVector<T>& newdata,
                   const ArrayVector<size_t>& newshape,
                   const std::array<unsigned,4>& new_elements={0,0,0,0});
  /*! Replace the data stored in the object
  @param newdata the new ArrayVector of data to be stored
  @param new_elements The number of scalar, eigenvector elements, vector elements, and matrix elements contained each newdata array
  @note This version of the method assumes each array in `newdata` is a vector
  */
  int replace_data(const ArrayVector<T>& newdata, const std::array<unsigned,4>& new_elements={0,0,0,0});

  //! A convenience function to find Q in Å⁻¹ before calling the other debye_waller_sum
  template<class R, class S = typename std::common_type<typename CostTraits<T>::type,R>::type>
  ArrayVector<S> debye_waller_sum(const LQVec<R>& Q, const R t_K) const;
  //! Perform the summation part of the Debye Waller factor
  template<class R, class S = typename std::common_type<typename CostTraits<T>::type,R>::type>
  ArrayVector<S> debye_waller_sum(const ArrayVector<R>& Q, const R t_K) const;
  /*! \brief Calculate the Debye-Waller factor for one or more Q points.

  @param Q An array of N 3-vectors expressed in reciprocal lattice units or
           inverse angstrom. If rlu are provided Q must be a LQVec such that Å⁻¹
           can be calculated.
  @param M The ion masses in meV⋅s²⋅Å⁻² // or should we take amu and convert? 1 amu == 1.03642688E-25 meV s² Å⁻²
  @param T The temperature in K.
  @returns The Debye-Waller factor at each of the N Q points, which is unitless.
  */
  template<class R, template<class> class A, class S = typename std::common_type<typename CostTraits<T>::type,R>::type>
  ArrayVector<S> debye_waller(const A<R>& Q, const std::vector<R>& M, const R t_K) const;

  //! Perform sanity checks before attempting to interpolate
  template<typename R> unsigned int check_before_interpolating(const ArrayVector<R>& x) const;
  //! Perform linear interpolation at the specified Reciprocal lattice points
  template<typename R> ArrayVector<T> interpolate_at(const LQVec<R>& x) const {return this->interpolate_at(x.get_xyz());}
  //! Perform linear interpolating at the specified points in the mesh's orthonormal frame
  template<typename R> ArrayVector<T> interpolate_at(const ArrayVector<R>& x) const;
private:
  //! Ensure that the provided elements make sense
  void check_elements(void);
};

#include "mesh.hpp"

#endif // _MESH_H_
