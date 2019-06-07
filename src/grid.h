/*! \file */
#ifndef _GRID_H_
#define _GRID_H_
typedef long slong; // ssize_t is only defined for gcc?

#include "latvec.h"
#include "neighbours.h"
#include "interpolation.h"
#include <omp.h>
// #include <complex>
// #include <memory>

#include "unsignedtosigned.h"
#include "munkres.h"

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
  ArrayVector<T> data;         //!< The stored ArrayVector indexed by `map`
  ArrayVector<size_t> shape;   //!< A second ArrayVector to indicate a possible higher-dimensional shape of each `data` array
public:
  // constructors
  MapGrid3(const size_t *n=default_n): map(nullptr), data(0,0), shape(1,0)
    { this->set_size(n); };
  MapGrid3(const size_t *n, const ArrayVector<T>& av): map(nullptr)
    { this->set_size(n); this->replace_data(av); };
  MapGrid3(const size_t *n, const slong *inmap, const ArrayVector<T>& av): map(nullptr)
    { this->set_size(n); this->replace_data(av); this->set_map(inmap,n,3u); };
  // copy constructor
  MapGrid3(const MapGrid3<T>& other): map(nullptr) {
    this->resize(other.size(0),other.size(1),other.size(2)); // sets N, calculates span, frees/allocates map memory if necessary
    for (size_t i=0; i<other.numel(); i++) this->map[i] = other.map[i];
    this->data = other.data;
    this->shape= other.shape;
  }
  // destructor
  ~MapGrid3(){
    if ( numel()>0 && map!=nullptr) delete[] map;
  }; // everything else is handled by itself
  // Copy constructor:
  // Assignment operator:
  MapGrid3<T>& operator=(const MapGrid3<T> &other){
    if (this != &other){
      this->resize(other.size(0),other.size(1),other.size(2)); // sets N, calculates span, frees/allocates map memory if necessary
      for (size_t i=0; i<other.numel(); i++) this->map[i] = other.map[i];
      this->data = other.data;
      this->shape= other.shape;
    }
    return *this;
  };
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
  /*! Replace the data stored in the object
  @param newdata the new ArrayVector of data to be stored
  @param newshape the shape information of each array in `newdata`
  */
  int replace_data(const ArrayVector<T>& newdata, const ArrayVector<size_t>& newshape);
  /*! Replace the data stored in the object
  @param newdata the new ArrayVector of data to be stored
  @note This version of the method assumes each array in `newdata` is a vector
  */
  int replace_data(const ArrayVector<T>& newdata);
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
  int sub2map(const size_t *s, size_t *m) const;
  /*! Find the mapping of a point given an array of its three subscripted indices
  @param s the subscripted indices
  @returns the valid mapping index or one more than the maximum mapping if s is invalid.
  */
  size_t sub2map(const size_t *s) const;
  /*! Find the mapping of a point given its linear index
  @param l the linear index
  @param[out] m a pointer to store the mapping index at
  @returns 0 if successful,
           1 if `l` is not a valid linear index,
          -1 if `l` is not a valid mapping index
  */
  int lin2map(const size_t  l, size_t *m) const;
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
  //! Change the size of the mapping grid given three new sizes
  size_t resize(const size_t n0, const size_t n1, const size_t n2);
  //! Change the size of the mapping grid given an array of the new sizes
  size_t resize(const size_t *n);
  //! Return the number of dimensions of each array in the `data` ArrayVector
  size_t data_ndim(void) const;
  //! Return the number of arrays in the `data` ArrayVector
  size_t num_data(void) const;
  //! Return the `shape` ArrayVector containing information about the arrays in the `data` ArrayVector
  ArrayVector<size_t> data_shape(void) const;
  //! Return the three sizes of the mapping grid as an ArrayVector
  ArrayVector<size_t> get_N(void) const;
  /*! Determine the neighbouring grid points of a given grid linear index
  @param centre The linear index to a point in the mapping grid
  @returns Up to 26 neighbouring linear indices, skipping points that are out of the grid
  */
  ArrayVector<size_t> get_neighbours(const size_t centre) const;
  /*! Attempt to determing a sorting permutation for each mapped point in the grid
  by comparing the difference in values stored at neighbouring points.
  The minimization of differences is performed via the Munkres Assignment Algorithm
  @returns An ArrayVector with `numel()==shape.getvalue(0)` and `size()=data.size()`
           containing a permuted indexing scheme which orders the elements at each grid point.
  */
  template <typename R=double>
  ArrayVector<size_t> sort_perm(const size_t n_scalar=0,
                                const size_t n_vector=0,
                                const size_t n_matrix=0,
                                const R scalar_weight=1,
                                const R vector_weight=1,
                                const R matrix_weight=1) const;
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
  /*!
  \brief Assign the elements at one point to the elements at a second point,
         storing the permutation array.

  When sorting the data stored in a grid, it is necessary to make a sorting
  assignment based soley on the difference in values at two points if a partial
  derivative can not be calculated.
  @param nS the number of scalar elements per object
  @param nV the number of vector elements per object
  @param nM the square root of the number of matrix elements per object
  @param scaleS weight factor for scalar cost
  @param scaleV weight factor for vector cost
  @param scaleM weight factor for matrix cost
  @param span the span of the elements of one object, nS*nV*nM*nM
  @param nobj the number of object per grid point
  @param[out] perm the permutation per grid point
  @param cidx the index of the central grid point
  @param nidx the index of the neighbouring grid point, which has been sorted
  @param vcf which vector cost function to use.
  */
  template<typename R>
  bool sort_difference(const size_t, const size_t, const size_t, const R,
                       const R, const R, const size_t, const size_t,
                       ArrayVector<size_t>&,
                       const size_t, const size_t,
                       const int vcf=0) const;
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

  @param nS the number of scalar elements per object
  @param nV the number of vector elements per object
  @param nM the square root of the number of matrix elements per object
  @param scaleS weight factor for scalar cost
  @param scaleV weight factor for vector cost
  @param scaleM weight factor for matrix cost
  @param span the span of the elements of one object, nS*nV*nM*nM
  @param nobj the number of object per grid point
  @param[out] perm the permutation per grid point
  @param cidx the index of the central grid point
  @param nidx the index of the neighbouring grid point, which has been sorted
  @param nnidx the index of the next neighbouring grid point,
               which has also been sorted
  @param vcf which vector cost function to use.
  */
  template<typename R>
  bool sort_derivative(const size_t, const size_t, const size_t, const R,
                       const R,  const R, const size_t, const size_t,
                       ArrayVector<size_t>&,
                       const size_t, const size_t, const size_t,
                       const int vcf=0) const;
  template<typename R>
  ArrayVector<size_t> new_sort_perm(const size_t, const size_t, const size_t,
                                    const R, const R, const R, const int vcf=0
                                   ) const;
  /*! \brief Sum over the data array

  Either add together the elements of the array stored at each mapped point
  or add together all of the arrays.
  @param axis The axis along which to perform the summation -- 1 adds arrays,
              0 (or not 1, really) adds elements of each array.
  @returns An ArrayVector with `numel()==1` for `axis=0` or `size()==1` for
           `axis=1`
  */
  ArrayVector<T> sum_data(const int axis) const{
    return this->data.sum(axis);
  };
  /*! \brief Return a constant reference to the data ArrayVector
  */
  const ArrayVector<T>& get_data() const {
    return this->data;
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
  InterpolateGrid3(const size_t *n=default_n, const double *z=default_zero, const double *s=default_step): MapGrid3<T>(n) { this->set_zero(z); this->set_step(s); };
  InterpolateGrid3(const size_t *n, const ArrayVector<T>& av, const double *z=default_zero, const double *s=default_step): MapGrid3<T>(n,av){ this->set_zero(z); this->set_step(s); };
  InterpolateGrid3(const size_t *n, const slong* inmap, const ArrayVector<T>& av, const double *z=default_zero, const double *s=default_step): MapGrid3<T>(n,inmap,av){ this->set_zero(z); this->set_step(s); };


  void set_zero(const double *newzero){ for(int i=0;i<3;i++) this->zero[i] = newzero[i]; };
  void set_step(const double *newstep){ for(int i=0;i<3;i++) this->step[i] = newstep[i]; };

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
  };
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
  };
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
  };
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
  };
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
  };
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
  };
  //! Return an ArrayVector of the grid x coordinates
  ArrayVector<double> get_grid_x() const {
    ArrayVector<double> x(1u,this->size(0));
    double s0=this->step[0], z0=this->zero[0];
    for (size_t i=0; i<x.size(); i++) x.insert( z0+s0*i, i);
    return x;
  };
  //! Return an ArrayVector of the grid y coordinates
  ArrayVector<double> get_grid_y() const {
    ArrayVector<double> y(1u,this->size(1));
    double s1=this->step[1], z1=this->zero[1];
    for (size_t j=0; j<y.size(); j++) y.insert( z1+s1*j, j);
    return y;
  };
  //! Return an ArrayVector of the grid z coordinates
  ArrayVector<double> get_grid_z() const {
    ArrayVector<double> z(1u,this->size(2));
    double s2=this->step[2], z2=this->zero[2];
    for (size_t k=0; k<z.size(); k++) z.insert( z2+s2*k, k);
    return z;
  };
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
  unsigned int nearest_index(const double *x, size_t *ijk) const {
    unsigned int out=0;
    slong tmp;
    for (int i=0; i<3; i++){
      tmp = (slong)( round( (x[i] - this->zero[i])/this->step[i] ) );
      if (tmp>=0 && tmp<this->size(i)){
        if (approx_scalar(this->step[i]*tmp + this->zero[i],x[i]))
          out += 1<<i; // exact match
      } else {
        if (tmp<0) {
          tmp = 0;
          out += 1<<(3+i); // underflow
        } else {
          tmp = this->size(i)-1;
          out += 1<<(6+i); // overflow
        }
      }
      ijk[i] = (size_t)(tmp);
    }
    return out;
  };
  // /*! Find the linear index of the closest grid point to a specified position
  // @param x The test position
  // @returns The linear index of the closest grid point -- even if `x` is out of bounds
  // */
  // size_t nearest_linear_index(const double *x) const {
  //   size_t ijk[3];
  //   int ret;
  //   if ( (ret = nearest_index(x,ijk)) ){
  //     // somehow indicate that we've hit one or more boundaries?
  //   }
  //   size_t lidx;
  //   ret += this->sub2lin(ijk,&lidx);
  //   return lidx;
  // }
  //! Perform sanity checks before attempting to interpolate
  template<typename R> int check_before_interpolating(const ArrayVector<R>& x){
    if (this->size(0)<2||this->size(1)<2||this->size(2)<2)
      throw std::runtime_error("Interpolation is only possible on grids with at least two elements in each dimension");
    if (this->data.size()==0)
      throw std::runtime_error("The grid must be filled before interpolating!");
    if (x.numel()!=3u)
      throw std::runtime_error("InterpolateGrid3 requires x values which are three-vectors.");
    return 0;
  }
  //! Perform linear interpolation at the specified points expressed in a Reciprocal lattice
  template<typename R> ArrayVector<T> linear_interpolate_at(const LQVec<R>& x){return this->linear_interpolate_at(x.get_xyz());}
  //! Perform linear interpolation at the specified points expressed in a Direct lattice
  template<typename R> ArrayVector<T> linear_interpolate_at(const LDVec<R>& x){return this->linear_interpolate_at(x.get_xyz());}
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
  template<typename R> ArrayVector<T> linear_interpolate_at(const ArrayVector<R>& x){
    this->check_before_interpolating(x);
    ArrayVector<T> out(this->data.numel(), x.size());
    size_t corners[8], ijk[3], cnt;
    unsigned int flg;
    int oob=0;
    double weights[8];
    std::vector<size_t> dirs, corner_count={1u,2u,4u,8u};

    for (size_t i=0; i<x.size(); i++){
      // std::cout << x.to_string(i) << std::endl;
      // find the closest grid subscripted indices to x[i]
      flg = this->nearest_index(x.datapointer(i), ijk );
      // std::cout << "\t closest index [" << std::to_string(ijk[0])
      //           << " " << std::to_string(ijk[1])
      //           << " " << std::to_string(ijk[2]) << "]";
      cnt = 1u; // will be modified if more than one-point interpolation
      // Alternatively, ignore out-of-bounds information by flg &= 7;
      if (flg > 7){
        std::string msg_flg = "grid.h::linear_interpolate_at Unsure what to do with flg = " + std::to_string(flg);
        throw std::runtime_error(msg_flg);
      }
      if (7==flg)/*+++*/{
        // std::cout << "\t exact match!" << std::endl;
        this->sub2map(ijk,corners); // set the first "corner" to this mapped index
        weights[0] = 1.0; // and the weight to one
      } else {
        // std::cout << "\t inexact match (flg="
        //           << std::to_string(flg) << ")";
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

        oob = corners_and_weights(this,this->zero,this->step,ijk,x.datapointer(i),corners,weights,3u,dirs);
        cnt=corner_count[dirs.size()];
        // std::cout << "\t requiring corners [ ";
        // for (size_t gst=0; gst<cnt; gst++) std::cout << std::to_string(corners[gst]) << " ";
        // std::cout << "] with weights [";
        // for (size_t gst=0; gst<cnt; gst++) std::cout << std::to_string(weights[gst]) << " ";
        // std::cout << "]" << std::endl;
        if (oob){
          std::string msg = "Point " + std::to_string(i) + " with x = " + x.to_string(i) + " has " + std::to_string(oob) + " corners out of bounds!";
          throw std::runtime_error(msg);
        }
      }
      // now do the actual interpolation:
      // extract an ArrayVector(this->data.numel(),cnt) of the corner Arrays
      // multiply all elements at each corner by the weight for that corner
      // sum over the corners, returning an ArrayVector(this->data.numel(),1u)
      // and set that ArrayVector as element i of the output ArrayVector
      unsafe_accumulate_to(this->data,cnt,corners,weights,out,i);
    }
    return out;
  };
  //! Perform linear interpolation in parallel at the specified points expressed in a Reciprocal lattice
  template<typename R> ArrayVector<T> parallel_linear_interpolate_at(const LQVec<R>& x, const int threads){return this->parallel_linear_interpolate_at(x.get_xyz(),threads);}
  //! Perform linear interpolation in parallel at the specified points expressed in a Direct lattice
  template<typename R> ArrayVector<T> parallel_linear_interpolate_at(const LDVec<R>& x, const int threads){return this->parallel_linear_interpolate_at(x.get_xyz(),threads);}
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
  template<typename R> ArrayVector<T> parallel_linear_interpolate_at(const ArrayVector<R>& x, const int threads){
    this->check_before_interpolating(x);
    ArrayVector<T> out(this->data.numel(), x.size());

    (threads>0) ? omp_set_num_threads(threads) : omp_set_num_threads(omp_get_max_threads());

    size_t corners[8], ijk[3], cnt=0u;
    unsigned int flg=0;
    int oob=0;
    double weights[8];
    slong xsize = unsigned_to_signed<slong,size_t>(x.size());
    std::vector<size_t> dirs, corner_count={1u,2u,4u,8u};

#pragma omp parallel for shared(x,out,corner_count) firstprivate(corners,ijk,weights,xsize) private(flg,oob,cnt,dirs)
    for (slong si=0; si<xsize; si++){
      size_t i = signed_to_unsigned<size_t,slong>(si);
      // find the closest grid subscripted indices to x[i]
      flg = this->nearest_index(x.datapointer(i), ijk );
      cnt = 1u; // will be modified if more than one-point interpolation
      if (flg > 7){
        std::string msg_flg = "grid.h::parallel_linear_interpolate_at Unsure what to do with flg = " + std::to_string(flg);
        throw std::runtime_error(msg_flg);
      }
      if (7 == flg)/*+++*/{
        this->sub2map(ijk,corners); // set the first "corner" to this mapped index
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

        oob = corners_and_weights(this,this->zero,this->step,ijk,x.datapointer(i),corners,weights,3u,dirs);
        cnt = corner_count[dirs.size()];
        if (oob){
          std::string msg = "Point " + std::to_string(i) + " with x = " + x.to_string(i) + " has " + std::to_string(oob) + " corners out of bounds!";
          throw std::runtime_error(msg);
        }
      }
      // now do the actual interpolation:
      // extract an ArrayVector(this->data.numel(),cnt) of the corner Arrays
      // multiply all elements at each corner by the weight for that corner
      // sum over the corners, returning an ArrayVector(this->data.numel(),1u)
      // and set that ArrayVector as element i of the output ArrayVector
      unsafe_accumulate_to(this->data,cnt,corners,weights,out,i);
    }
    return out;
  };
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
  constexpr static T max = std::numeric_limits<T>::max();
};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T> struct GridDiffTraits<std::complex<T>>{
  using type = T;
  constexpr static T max = std::numeric_limits<T>::max();
};
#endif

#include "grid.hpp"
#include "grid_sorting.hpp"

#endif
