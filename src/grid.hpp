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
#ifndef _GRID_H_
#define _GRID_H_
typedef long slong; // ssize_t is only defined for gcc?

#include "latvec.hpp"
#include "neighbours.hpp"
#include "interpolation.hpp"
#include <tuple>
#include <queue>
#include <omp.h>
// #include <complex>
// #include <memory>
#include "interpolation_data.hpp"
#include "utilities.hpp"
#include "permutation.hpp"

// A grid is a 3 (or 4) dimensional object that for a given index, e.g.,
// [i][j][k], contains the (linear) index into a second ArrayVector object.

const size_t default_n[3] = { 0u,0u,0u };

/*! \brief A class holding a 3D grid which maps into an internal ArrayVector

The MapGrid3 holds a mapping grid where each entry is either a valid index into
a held ArrayVector or is `-1` to indicate that the grid point does not map.
*/
template<class T> class MapGrid3{
protected:
  size_t N[3];                 //!< The number of points along each axis of the grid
  size_t span[3];              //!< The span along each axis, allowing conversion from subscripted to linear indexing
  slong *map;                  //!< The mapping grid
  InterpolationData<T> data_;
public:
  // constructors
  MapGrid3(const size_t *n=default_n): map(nullptr)
    { this->set_size(n); }
  MapGrid3(const size_t *n, const ArrayVector<T>& av): map(nullptr)
    { this->set_size(n); this->replace_data(av); }
  MapGrid3(const size_t *n, const slong *inmap, const ArrayVector<T>& av): map(nullptr)
    { this->set_size(n); this->replace_data(av); this->set_map(inmap,n,3u); }
  // copy constructor
  MapGrid3(const MapGrid3<T>& other): map(nullptr) {
    this->resize(other.size(0),other.size(1),other.size(2)); // sets N, calculates span, frees/allocates map memory if necessary
    for (size_t i=0; i<other.numel(); i++) this->map[i] = other.map[i];
    this->data_ = other.data();
  }
  // destructor
  ~MapGrid3(){
    if ( numel()>0 && map!=nullptr) delete[] map;
  } // everything else is handled by itself
  // Copy constructor:
  // Assignment operator:
  MapGrid3<T>& operator=(const MapGrid3<T> &other){
    if (this != &other){
      this->resize(other.size(0),other.size(1),other.size(2)); // sets N, calculates span, frees/allocates map memory if necessary
      for (size_t i=0; i<other.numel(); i++) this->map[i] = other.map[i];
      this->data_ = other.data();
    }
    return *this;
  }
  //! Print the number of points along each axis to the console
  void print_N(const bool nl=false) const;
  //! Print the span along each axis to the console
  void print_span(const bool nl=false) const;
  //! Print the whole map to the console as stacked slices in the first two dimensions
  void print_map(void) const;
  //! Set the mapping index of each grid point to its linear index
  int set_map(void);
  /*! Copy-in a mapping grid
  @param inmap The existant mapping grid to copy into the object
  @param n The number of points along each axis in `inmap`
  @param d The number of dimensions in `n` and `inmap`
  */
  int set_map(const slong* inmap, const size_t* n, const size_t d);
  //! Copy-in a mapping grid without checking that it is the right size
  int unsafe_set_map(slong *inmap);
  //! Copy-out the mapping grid without verifying that the storage location can hold it
  size_t unsafe_get_map(slong *outmap) const;
  //
  /*! Determine the maximum mapping index of a map
  @param map2check The mapping array to search through
  @param num2check The number of entries in the mapping index
  */
  size_t maximum_mapping(const slong *map2check, const size_t num2check) const;
  /*! Determine the maximum mapping index of a map, assuming the size of `map2check` is the same as `map`
  @param map2check The mapping array to search through
  */
  size_t maximum_mapping(const slong *map2check) const;
  //! Determine the maximum mapping index of `map`
  size_t maximum_mapping(void) const;
  //! Return the number of valid mappings within `map`
  size_t valid_mapping_count(void) const;
  /*! Determine if the interal mapping is consistent with an ArrayVector
  @param data2check Will be compared against the maximum mapping in `map`
  */
  int check_map(const ArrayVector<T>& data2check) const;
  //! Determine if `map` is consistent with `data`
  int check_map(void) const;
  //! Replace the data stored in the object
  template<typename... A> int replace_data(A... args) { this->data_.replace_data(args...); return this->check_map();}
  //! Calculate the linear index of a point given its three subscripted indices
  size_t sub2lin(const size_t i, const size_t j, const size_t k) const;
  /*! Calculate the linear index of a point given an array of its three subscripted indices
  @param s the subscripted index to convert
  @param[out] l the calculated linear index
  @returns 0 if successful, 1 if unsuccessful
  */
  int sub2lin(const size_t *s, size_t *l) const;
  /*! Calculate the subscripted indices of a point given its linear index
  @param l the linear index
  @param[out] s a pointer to store the subscripted index
  @returns 0 if successful, 1 if unsuccessful.
  */
  int lin2sub(const size_t  l, size_t *s) const;
  /*! Find the mapping of a point given an array of its three subscripted indices
  @param s the subscripted index
  @param[out] m the mapped index
  @returns 0 if successful,
           1 if `s` is not a valid subscripted index,
          -1 if `s` is not a valid mapping index
  */
  int sub2map(const size_t* s, size_t& m) const;
  /*! Find the mapping of a point given an array of its three subscripted indices
  @param s the subscripted indices
  @returns the valid mapping index or one more than the maximum mapping if s is invalid.
  */
  size_t sub2map(const size_t *s) const;
  // /*! Find the mapping of a point given its linear index
  // @param l the linear index
  // @param[out] m a pointer to store the mapping index at
  // @returns 0 if successful,
  //          1 if `l` is not a valid linear index,
  //         -1 if `l` is not a valid mapping index
  // */
  // int lin2map(const size_t  l, size_t *m) const;
  /*! Find the mapping of a point given its linear index
  @param l the linear index
  @param[out] m a reference to store the mapping index
  @returns 0 if successful,
           1 if `l` is not a valid linear index,
          -1 if `l` is not a valid mapping index
  */
  int lin2map(const size_t  l, size_t& m) const;
  /*! Find the linear index of a point given its mapping, by searching through
  the mapping array -- this is almost certainly highly ineficcient and should
  be avoided.
  @param m the mapping index
  @param[out] l a reference to store the linear index
  @return 0 if `l` was found successfully,
          1 if `m` is not a valid mapping
         -1 if `m` is valid by `l` was not found.
  */
  int map2lin(const size_t m, size_t& l) const;
  //
  //! Return the total number of elements in the mapping grid
  size_t numel(void) const;
  //! Return the number of elements in the mapping grid along a given dimension
  size_t size(const size_t i) const;
  // //! Return the span of elements neccessary to move between neighbouring points
  // size_t span(const size_t i) const;
  //! Return the stride (in bytes) necessary to move between neighbouring points
  size_t stride(const size_t i) const;
  //! Change the size of the mapping grid given three new sizes
  size_t resize(const size_t n0, const size_t n1, const size_t n2);
  //! Change the size of the mapping grid given an array of the new sizes
  size_t resize(const size_t *n);
  // Get a constant reference to the stored data
  const InterpolationData<T>& data(void) const {return data_;}
  //! Return the three sizes of the mapping grid as an ArrayVector
  ArrayVector<size_t> get_N(void) const;
  /*! Determine the neighbouring grid points of a given grid linear index
  @param centre The linear index to a point in the mapping grid
  @returns Up to 26 neighbouring linear indices, skipping points that are out
           of the grid or do not represent a valid mapping.
  */
  ArrayVector<size_t> get_neighbours(const size_t centre) const;

  /*! \brief Determine which neighbours have been sorted and, of those,
             which can be used in partial-derivative sorting.

  Given the linear index of a grid point which is itself a valid mapped point,
  first find all valid mapped neighbouring points. For each of the valid mapped
  neighbouring points identify if there is a second mapped point along the same
  direction which can be used for partial-derivative sorting.
  @param map The MapGrid3 object containing data undergoing sorting.
  @param sorted A logical vector indicating which mapped points have been sorted
  @param clin The linear index for the central valid-mapped-point.
  @returns A std::vector which
           (a) is empty, in which case no sorted neighbours exist,
           (b) has a single map index of a sorted neighbour, in which case no
               sorted neighbours could be used in partial-derivative sorting
           (c) has two map indexes: a sorted neighbour *and* a second sorted
               point suitable for partial-derivative sorting, in that order.
           In the case of (b) or (c) there may be additional sorted neighbours.
  */
  std::vector<size_t> find_sorted_neighbours(const std::vector<bool>& sorted,
                                             const size_t clin) const;
  std::vector<size_t> find_unsorted_neighbours(const std::vector<bool>& sorted,
                                             const size_t centre_map_idx) const;
  /*!
  \brief Assign the elements at one point to the elements at a second point,
         storing the permutation array.

  When sorting the data stored in a grid, it is necessary to make a sorting
  assignment based soley on the difference in values at two points if a partial
  derivative can not be calculated.
  @param scaleS weight factor for scalar cost
  @param scaleE weight factor for eigenvector cost
  @param scaleV weight factor for vector cost
  @param scaleM weight factor for matrix cost
  @param span the span of the elements of one object, nS*nV*nM*nM
  @param nobj the number of object per grid point
  @param[out] perm the permutation per grid point
  @param cidx the index of the central grid point
  @param nidx the index of the neighbouring grid point, which has been sorted
  @param ecf which eigenvector cost function to use.
  @param vcf which vector cost function to use.
  */
  template<typename R>
  bool sort_difference(const R, const R, const R,
                       const size_t, const size_t,
                       ArrayVector<size_t>&,
                       const size_t, const size_t, const int vcf) const;
  /*!
  \brief Use the sorted elements at two neighbouring grid points to determine
         the sorting permutation of the elements at a third neighbouring grid
         point using a finite differences approximation to the derivative.

  Given two grid points which have been sorted and lie along a line one and two
  steps away, respectively, from a to-be-sorted central grid point, estimate
  the values of each sorted element at the central point and use the estimate to
  determine an assignment permutation for the central grid point.

  In general one can predict the value of a function at x if its value and
  slope are known at a point x₀ to by y₀ and m, respectively, by

      y = m (x-x₀) + y₀

  The slope is the first derivative of the function and can be approximated by
  the difference in the value of the function at two points

      y = [(y₀-y₁)/(x₀-x₁)](x-x₀) + y₀

  If the difference between x₀ and x₁ is the same as the difference between x
  and x₀, as it is in our gridded system, then this simplifies to

      y = 2y₀ - y₁

  where y₀ is the value at the neighbouring point and y₁ is the value at the
  next neighbouring point.

  @param scaleS weight factor for scalar cost
  @param scaleE weight factor for eigenvector cost
  @param scaleV weight factor for vector cost
  @param scaleM weight factor for matrix cost
  @param span the span of the elements of one object, nS*nV*nM*nM
  @param nobj the number of object per grid point
  @param[out] perm the permutation per grid point
  @param cidx the index of the central grid point
  @param nidx the index of the neighbouring grid point, which has been sorted
  @param nnidx the index of the next neighbouring grid point,
               which has also been sorted
  @param ecf which eigenvector cost function to use.
  @param vcf which vector cost function to use.
  */
  template<typename R>
  bool sort_derivative(const R, const R, const R,
                       const size_t, const size_t,
                       ArrayVector<size_t>&,
                       const size_t, const size_t, const size_t, const int vcf) const;
  template<typename R>
  size_t sort_recursion(const size_t centre,
                        const R wS, const R wV, const R wM, const int vcf,
                        const size_t span, const size_t nobj,
                        ArrayVector<size_t>& perm,
                        std::vector<bool>& sorted,
                        std::vector<bool>& locked) const;
  template<typename R=double>
  ArrayVector<size_t> centre_sort_perm(const R scalar_weight=R(1),
                                       const R vector_weight=R(1),
                                       const R matrix_weight=R(1),
                                       const int vcf=0) const;
  //
  std::tuple<std::vector<size_t>,std::vector<size_t>,std::vector<size_t>>
  classify_neighbours(const std::vector<bool>& sorted, const size_t centre_map_idx) const;
  //
  std::tuple<std::vector<size_t>,std::vector<size_t>>
  classify_sorted_neighbours(const std::vector<bool>& sorted, const size_t centre_map_idx) const;
  //
  template<class R>
  bool multi_sort_derivative(const R scales[3], const int func,
    const size_t spobj[2], ArrayVector<size_t>& perm, std::vector<bool>& sorted,
    const size_t cidx, const std::vector<size_t> nidx,
    const std::vector<size_t> nnidx, const size_t no_pairs) const;
  //
  template<class R>
  bool multi_sort_derivative_all(const R scales[3], const int func,
    const size_t spobj[2], ArrayVector<size_t>& perm, std::vector<bool>& sorted,
    const size_t cidx, const std::vector<size_t> nidx,
    const std::vector<size_t> nnidx/*, const size_t no_pairs*/) const;
  //
  template<class R>
  bool multi_sort_difference(const R weights[3], const int func,
    const size_t spobj[2], ArrayVector<size_t>& perm, std::vector<bool>& sorted,
    const size_t cidx, const std::vector<size_t> nidx) const;
  //
  template<class R> size_t multi_sort_recursion(const size_t centre,
    const R weights[3], const int func, const size_t spobj[2],
    ArrayVector<size_t>& perm, std::vector<bool>& sorted,
    std::vector<bool>& locked, std::vector<size_t>& visited) const;
  //
  template<class R> size_t multi_sort(const size_t centre, const R weights[3],
    const int func, const size_t spobj[2], ArrayVector<size_t>& perm,
    std::vector<bool>& sorted, std::vector<bool>& locked,
    std::vector<size_t>& visited) const;
  //
  template<typename R=double>
  ArrayVector<size_t> multi_sort_perm(
    const R scalar_weight=R(1),
    const R vector_weight=R(1),
    const R matrix_weight=R(1),
    const int vcf=0) const ;
  // Calculate the Debye-Waller factor for the provided Q points and ion masses
  template<template<class> class A>
  ArrayVector<double> debye_waller(const A<double>& Q, const std::vector<double>& M, const double t_K) const{
    return data_.debye_waller(Q,M,t_K);
  }
protected:
  void set_size(const size_t *n);
  void calc_span();
  void instantiate_map();
  bool valid_mapping(const size_t l) const;
  bool valid_mapping(const size_t i, const size_t j, const size_t k) const;
  bool is_inbounds(const size_t i, const size_t j, const size_t k) const;
  bool is_inbounds(const size_t* s) const;
};

// TODO: allow for more-compact data arrays with appropriate reducing maps

const double default_zero[3] = {0.,0.,0.};
const double default_step[3] = {1.,1.,1.};

/*! \brief Extends the MapGrid3 class to have positions of each grid point enabling interpolation between mapped points
*/
template<class T> class InterpolateGrid3: public MapGrid3<T>{
  double zero[3]; //!< the 3-vector position of `map[0]`
  double step[3]; //!< the step size along each direction of the grid
public:
  InterpolateGrid3(const size_t *n=default_n, const double *z=default_zero, const double *s=default_step): MapGrid3<T>(n) { this->set_zero(z); this->set_step(s); }
  InterpolateGrid3(const size_t *n, const ArrayVector<T>& av, const double *z=default_zero, const double *s=default_step): MapGrid3<T>(n,av){ this->set_zero(z); this->set_step(s); }
  InterpolateGrid3(const size_t *n, const slong* inmap, const ArrayVector<T>& av, const double *z=default_zero, const double *s=default_step): MapGrid3<T>(n,inmap,av){ this->set_zero(z); this->set_step(s); }


  void set_zero(const double *newzero){ for(int i=0;i<3;i++) this->zero[i] = newzero[i]; }
  void set_step(const double *newstep){ for(int i=0;i<3;i++) this->step[i] = newstep[i]; }

  /*! \brief Return the positions of all grid points
  @param maxN the maximum number of grid points which can be stored in `xyz`
  @param[out] xyz A location to store the grid point coordinates
  @returns the number of grid points returned
  */
  size_t get_grid_xyz(const size_t maxN, double *xyz) const {
    if ( maxN < this->numel() ) return 0;
    size_t i,j,k, n0=this->size(0), n1=this->size(1), n2=this->size(2);
    size_t cnt = 0;
    double s0=this->step[0], z0=this->zero[0], s1=this->step[1], z1=this->zero[1], s2=this->step[2], z2=this->zero[2];
    for (i=0; i<n0; i++)
      for (j=0; j<n1; j++)
        for (k=0; k<n2; k++){
          xyz[ 3*cnt    ] = z0 + s0*i;
          xyz[ 3*cnt +1 ] = z1 + s1*j;
          xyz[ 3*cnt +2 ] = z2 + s2*k;
          cnt++;
        }
    return cnt;
  }
  /*! \brief Return the x coordiates of the grid points
  @param maxN the maximum number of values which can be stored in `x`
  @param[out] x A location to store the x values
  @returns the number of values stored
  */
  size_t get_grid_x(const size_t maxN, double *x) const {
    if ( maxN < this->size(0) ) return 0;
    size_t cnt = 0;
    double s0=this->step[0], z0=this->zero[0];
    for (size_t i=0; i<this->size(0); i++) x[ cnt++ ] = z0 + s0*i;
    return cnt;
  }
  /*! \brief Return the y coordiates of the grid points
  @param maxN the maximum number of values which can be stored in `y`
  @param[out] y A location to store the y values
  @returns the number of values stored
  */
  size_t get_grid_y(const size_t maxN, double *y) const {
    if ( maxN < this->size(1) ) return 0;
    size_t cnt = 0;
    double s1=this->step[1], z1=this->zero[1];
    for (size_t j=0; j<this->size(1); j++) y[ cnt++ ] = z1 + s1*j;
    return cnt;
  }
  /*! \brief Return the z coordiates of the grid points
  @param maxN the maximum number of values which can be stored in `z`
  @param[out] z A location to store the z values
  @returns the number of values stored
  */
  size_t get_grid_z(const size_t maxN, double *z) const {
    if ( maxN < this->size(2) ) return 0;
    size_t cnt = 0;
    double s2=this->step[2], z2=this->zero[2];
    for (size_t k=0; k<this->size(2); k++) z[ cnt++ ] = z2 + s2*k;
    return cnt;
  }
  //! Return an ArrayVector of the coordinates for all grid points
  ArrayVector<double> get_grid_xyz() const {
    ArrayVector<double> xyz(3u,this->numel());
    size_t i,j,k, n0=this->size(0), n1=this->size(1), n2=this->size(2);
    size_t cnt = 0;
    double s0=this->step[0], z0=this->zero[0], s1=this->step[1], z1=this->zero[1], s2=this->step[2], z2=this->zero[2];
    for (i=0; i<n0; i++)
      for (j=0; j<n1; j++)
        for (k=0; k<n2; k++){
          xyz.insert( z0+s0*i, cnt  ,0);
          xyz.insert( z1+s1*j, cnt  ,1);
          xyz.insert( z2+s2*k, cnt++,2);
        }
    return xyz;
  }
  //! Return an ArrayVector of the coordinates for all mapped grid points
  ArrayVector<double> get_mapped_xyz() const {
    ArrayVector<double> xyz(3u,this->valid_mapping_count());
    size_t i,j,k, n0=this->size(0), n1=this->size(1), n2=this->size(2);
    size_t cnt = 0;
    double s0=this->step[0], z0=this->zero[0], s1=this->step[1], z1=this->zero[1], s2=this->step[2], z2=this->zero[2];
    for (i=0; i<n0; ++i)
      for (j=0; j<n1; ++j)
        for (k=0; k<n2; ++k){
          if ( this->valid_mapping(i,j,k) ) {
            xyz.insert( z0+s0*i, cnt  ,0);
            xyz.insert( z1+s1*j, cnt  ,1);
            xyz.insert( z2+s2*k, cnt++,2);
          } else {
            // printf("Invalid mapping at (%u, %u, %u)\n",i,j,k);
          }
        }
    return xyz;
  }
  //! Return an ArrayVector of the grid x coordinates
  ArrayVector<double> get_grid_x() const {
    ArrayVector<double> x(1u,this->size(0));
    double s0=this->step[0], z0=this->zero[0];
    for (size_t i=0; i<x.size(); i++) x.insert( z0+s0*i, i);
    return x;
  }
  //! Return an ArrayVector of the grid y coordinates
  ArrayVector<double> get_grid_y() const {
    ArrayVector<double> y(1u,this->size(1));
    double s1=this->step[1], z1=this->zero[1];
    for (size_t j=0; j<y.size(); j++) y.insert( z1+s1*j, j);
    return y;
  }
  //! Return an ArrayVector of the grid z coordinates
  ArrayVector<double> get_grid_z() const {
    ArrayVector<double> z(1u,this->size(2));
    double s2=this->step[2], z2=this->zero[2];
    for (size_t k=0; k<z.size(); k++) z.insert( z2+s2*k, k);
    return z;
  }
  /*! Find the subscripted indices of the closest grid point to a specified position
  @param x The test position
  @param[out] ijk A storage location for the subscripted indices
  @returns An integer with detailed information about if and how `x` is out of
           bounds of the grid.
  @note `out` encodes details about how each component of `x` is located in
         relation to the boundaries of the grid, using its bits as flags.
         For a given axis n, a number fₙ takes one of four
         values {0,2ⁿ,2ⁿ⁺³,2ⁿ⁺⁶} to indicate that  `x` is between two grid points,
         within machine precision of a grid point, smaller than the lowest grid
         point, or larger than the highest grid point, respectively.
         `out` is then f₀+f₁+f₂.
  */
  unsigned int nearest_index(const double *x, size_t *ijk, unsigned int mask=0u) const {
    unsigned int out=0;
    slong tmp;
    for (int i=0; i<3; i++){
      tmp = (slong)( round( (x[i] - this->zero[i])/this->step[i] ) );
      if (tmp>=0 && static_cast<size_t>(tmp)<this->size(i)){
        if (approx_scalar(this->step[i]*tmp + this->zero[i],x[i]))
          out += 1<<i; // exact match
      } else {
        if (tmp<0) {
          tmp = 0;
          if ((mask>>i)%2 == 1) // underflow allowed
            out += 1<<i; // exact match to avoid out-of-bounds interpolation
          else
            out += 1<<(3+i); // underflow
        } else {
          tmp = this->size(i)-1;
          if ((mask>>i)%2 == 1) // overflow allowed
            out += 1<<i; // exact match to avoid out-of-bounds interpolation
          else
            out += 1<<(6+i); // overflow
        }
      }
      ijk[i] = (size_t)(tmp);
    }
    return out;
  }
  /*! Find the subscripted indices of the grid point to a specified position
      which is smaller than the position in all coordinates.
  @param x The test position
  @param[out] ijk A storage location for the subscripted indices
  @returns An integer with detailed information about if and how `x` is out of
           bounds of the grid.
  @note `out` encodes details about how each component of `x` is located in
         relation to the boundaries of the grid, using its bits as flags.
         For a given axis n, a number fₙ takes one of four
         values {0,2ⁿ,2ⁿ⁺³,2ⁿ⁺⁶} to indicate that  `x` is between two grid points,
         within machine precision of a grid point, smaller than the lowest grid
         point, or larger than the highest grid point, respectively.
         `out` is then f₀+f₁+f₂.
  */
  unsigned int floor_index(const double *x, size_t *ijk, unsigned int mask=0u) const {
    unsigned int out=0;
    slong tmp;
    for (int i=0; i<3; i++){
      tmp = (slong)( floor( (x[i] - this->zero[i])/this->step[i] ) );
      if (tmp>=0 && static_cast<size_t>(tmp)<this->size(i)){
        if (approx_scalar(this->step[i]*tmp + this->zero[i],x[i]))
          out += 1<<i; // exact match
      } else {
        if (tmp<0) {
          tmp = 0;
          if ((mask>>i)%2 == 1) // underflow allowed
            out += 1<<i; // exact match to avoid out-of-bounds interpolation
          else
            out += 1<<(3+i); // underflow
        } else {
          tmp = this->size(i)-1;
          if ((mask>>i)%2 == 1) // overflow allowed
            out += 1<<i; // exact match to avoid out-of-bounds interpolation
          else
            out += 1<<(6+i); // overflow
        }
      }
      ijk[i] = (size_t)(tmp);
    }
    return out;
  }
  //! Perform sanity checks before attempting to interpolate
  template<typename R> unsigned int check_before_interpolating(const ArrayVector<R>& x) const {
    unsigned int mask=0u;
    if (this->data_.size()==0)
      throw std::runtime_error("The grid must be filled before interpolating!");
    if (x.numel()!=3u)
      throw std::runtime_error("InterpolateGrid3 requires x values which are three-vectors.");
    if (this->size(0)<2||this->size(1)<2||this->size(2)<2){
      size_t ijk[3];
      bool allok = true;
      unsigned int flg;
      std::vector<int> must_be_in_bounds;
      for (int i=0; i<3; ++i){
        if (this->size(i)<2)
          mask += 1<<i;
        else
          must_be_in_bounds.push_back(i);
      }
      for (size_t i=0; i<x.size(); ++i){
        flg = this->nearest_index(x.data(i), ijk, mask);
        allok &= 0==(flg>>3);
        if (!allok) break;
      }
      if (!allok){
        std::string msg = "Interpolation on a grid with shape [";
        for (size_t i=0; i<3; ++i) msg += " " + std::to_string(this->size(i));
        msg += " ] can only be performed for points with dot(x, [";
        for (size_t i=0; i<3; ++i) msg += (this->size(i)<2 ? " 0" : " 1");
        msg += " ]) within the grid limits.";
        throw std::runtime_error(msg);
      }
    }
    return mask;
  }
  //! Perform linear interpolation at the specified points expressed in a Reciprocal lattice
  template<typename R> ArrayVector<T> linear_interpolate_at(const LQVec<R>& x) const {return this->linear_interpolate_at(x.get_xyz());}
  //! Perform linear interpolation at the specified points expressed in a Direct lattice
  template<typename R> ArrayVector<T> linear_interpolate_at(const LDVec<R>& x) const {return this->linear_interpolate_at(x.get_xyz());}
  /*! Perform linear interpolation at the specified points expressed in an orthonormal frame
  @param x The coordinates to interpolate at expressed in the same orthonormal frame as the mapping grid
  @returns An ArrayVector of the itnerpolated values
  @note In the event that one or more coordinates of a vector in `x` is an exact
        match for a grid point this routine will perform a lower-dimensional
        interpolation. For three exact matches the point in x *is* a grid point
        and no interpolation is performed, for two exact matches 1D linear
        interpolation is performed, for one exact match bilinear 2D interpolation
        is used, and for no exact matches the method uses trilinear interpolation.
  */
  template<typename R> ArrayVector<T> linear_interpolate_at(const ArrayVector<R>& x) const {
    unsigned int mask = this->check_before_interpolating(x);
    ArrayVector<T> out(this->data_.numel(), x.size());
    std::vector<size_t> corners;
    std::vector<double> weights;
    size_t ijk[3], cnt;
    unsigned int flg;
    int oob=0;
    std::vector<size_t> dirs, corner_count={1u,2u,4u,8u};
    for (size_t i=0; i<x.size(); i++){
      corners.resize(8u);
      weights.resize(8u);
      // find the closest grid subscripted indices to x[i]
      flg = this->nearest_index(x.data(i), ijk, mask);
      // flg = this->floor_index(x.data(i), ijk, mask );
      cnt = 1u; // will be modified if more than one-point interpolation
      // Alternatively, ignore out-of-bounds information by flg &= 7;
      if (flg > 7){
        std::string msg_flg = "grid.h::linear_interpolate_at";
        msg_flg += " Unsure what to do with flg = " + std::to_string(flg);
        msg_flg += " when mask = " + std::to_string(mask);
        throw std::runtime_error(msg_flg);
      }
      if (7==flg)/*+++*/{
        this->sub2map(ijk,corners[0]); // set the first "corner" to this mapped index
        weights[0] = 1.0; // and the weight to one
      } else {
        if (!flg)/*xxx*/{
          dirs.resize(3);
          dirs[0]=0u; dirs[1]=1u; dirs[2]=2u;
        }
        if ( flg && !(flg&(flg-1)) ) dirs.resize(2); // flg&(flg-1) is zero if flg is a power of 2 (or zero)
        if (1==flg)/*+xx*/{ dirs[0] = 1u; dirs[1] = 2u;}
        if (2==flg)/*x+x*/{ dirs[0] = 0u; dirs[1] = 2u;}
        if (4==flg)/*xx+*/{ dirs[0] = 0u; dirs[1] = 1u;}

        if ( flg && (flg&(flg-1)) ) dirs.resize(1);
        if (3==flg)/*++x*/ dirs[0] = 2u;
        if (5==flg)/*+x+*/ dirs[0] = 1u;
        if (6==flg)/*x++*/ dirs[0] = 0u;

        oob = corners_and_weights(this,this->zero,this->step,ijk,x.data(i),corners.data(),weights.data(),3u,dirs);
        cnt = corner_count[dirs.size()];
        if (oob){
          std::string msg = "Point " + std::to_string(i) + " with x = " + x.to_string(i) + " has " + std::to_string(oob) + " corners out of bounds!";
          throw std::runtime_error(msg);
        }
      }
      // now do the actual interpolation:
      // extract an ArrayVector(data_.numel(),cnt) of the corner Arrays
      // multiply all elements at each corner by the weight for that corner
      // sum over the corners, returning an ArrayVector(data_.numel(),1u)
      // and set that ArrayVector as element i of the output ArrayVector
      corners.resize(cnt);
      weights.resize(cnt);
      this->data_.interpolate_at(corners, weights, out, i);
    }
    return out;
  }
  //! Perform linear interpolation in parallel at the specified points expressed in a Reciprocal lattice
  template<typename R> ArrayVector<T> parallel_linear_interpolate_at(const LQVec<R>& x, const int threads) const {return this->parallel_linear_interpolate_at(x.get_xyz(),threads);}
  //! Perform linear interpolation in parallel at the specified points expressed in a Direct lattice
  template<typename R> ArrayVector<T> parallel_linear_interpolate_at(const LDVec<R>& x, const int threads) const {return this->parallel_linear_interpolate_at(x.get_xyz(),threads);}
  /*! Perform linear interpolation in parallel at the specified points expressed in an orthonormal frame
  @param x The coordinates to interpolate at expressed in the same orthonormal frame as the mapping grid
  @param threads The number of OpenMP threads to use, `omp_get_max_threads()` if `threads`≤0
  @returns An ArrayVector of the interpolated values
  @note In the event that one or more coordinates of a vector in `x` is an exact
        match for a grid point this routine will perform a lower-dimensional
        interpolation. For three exact matches the point in x *is* a grid point
        and no interpolation is performed, for two exact matches 1D linear
        interpolation is performed, for one exact match bilinear 2D interpolation
        is used, and for no exact matches the method uses trilinear interpolation.
  */
  template<typename R> ArrayVector<T> parallel_linear_interpolate_at(const ArrayVector<R>& x, const int threads) const {
    unsigned int mask = this->check_before_interpolating(x);
    ArrayVector<T> out(this->data_.numel(), x.size());

    (threads>0) ? omp_set_num_threads(threads) : omp_set_num_threads(omp_get_max_threads());
    std::vector<size_t> corners(8,0);
    std::vector<double> weights(8,0.);
    size_t ijk[3], cnt=0u;
    unsigned int flg=0;
    int oob=0;
    slong xsize = unsigned_to_signed<slong,size_t>(x.size());
    std::vector<size_t> dirs, corner_count={1u,2u,4u,8u};
    size_t n_oob{0};
#pragma omp parallel for default(none) shared(x,out,corner_count,mask) firstprivate(corners,ijk,weights,xsize) private(flg,oob,cnt,dirs) reduction(+:n_oob) schedule(dynamic)
    for (slong si=0; si<xsize; si++){
      size_t i = signed_to_unsigned<size_t,slong>(si);
      corners.resize(8u);
      weights.resize(8u);
      // find the closest grid subscripted indices to x[i]
      flg = this->nearest_index(x.data(i), ijk, mask );
      cnt = 1u; // will be modified if more than one-point interpolation
      oob = 0;
      if (flg <= 7){
        if (7 == flg)/*+++*/{
          this->sub2map(ijk,corners[0]); // set the first "corner" to this mapped index
          weights[0] = 1.0; // and the weight to one
        } else {
          if (0==flg)/*xxx*/{
            dirs.resize(3);
            dirs[0]=0u; dirs[1]=1u; dirs[2]=2u;
          }
          if (1==flg || 2==flg || 4==flg) dirs.resize(2);
          if (1==flg)/*+xx*/{ dirs[0] = 1u; dirs[1] = 2u;}
          if (2==flg)/*x+x*/{ dirs[0] = 0u; dirs[1] = 2u;}
          if (4==flg)/*xx+*/{ dirs[0] = 0u; dirs[1] = 1u;}

          if (3==flg || 5==flg || 6==flg) dirs.resize(1);
          if (3==flg)/*++x*/ dirs[0] = 2u;
          if (5==flg)/*+x+*/ dirs[0] = 1u;
          if (6==flg)/*x++*/ dirs[0] = 0u;

          oob = corners_and_weights(this,this->zero,this->step,ijk,x.data(i),corners.data(),weights.data(),3u,dirs);
          cnt = corner_count[dirs.size()];
        }
        if (!oob) {
          corners.resize(cnt);
          weights.resize(cnt);
          this->data_.interpolate_at(corners, weights, out, i);
        } else {
          ++n_oob;
        }
      } else {
        ++n_oob;
      }
    }
    if (n_oob > 0){
      std::string msg = "parallel_linear_interpolate_at failed with ";
      msg += std::to_string(n_oob) + " out of bounds points.";
      throw std::runtime_error(msg);
    }
    return out;
  }
};


/*! \brief Type information for MapGrid3, MapGrid4, and their subclasses.

In order to compare the array of information stored at each mapped grid point
in order to, e.g., create a sorting permutation, it is necessary to define the
typename of the comparison criterion. Since MapGrid3 and MapGrid4 objects can
hold complex data but the comparison is always performed on a real valued scalar
or array this struct provides a convenient way of providing the comparison type
for templated functions.

| template typename | type | max |
| T | T | std::numeric_limits<T>::max() |
| std::complex<T> | T | std::numeric_limits<T>::max() |

*/
template<class T> struct GridDiffTraits{
  using type = T;
  constexpr static T max = (std::numeric_limits<T>::max)();
};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T> struct GridDiffTraits<std::complex<T>>{
  using type = T;
  constexpr static T max = (std::numeric_limits<T>::max)();
};
#endif

#include "grid.tpp"
#include "grid_sorting.tpp"

#endif
