/*! \file */
#ifndef _GRID4_H_
#define _GRID4_H_

typedef long slong; // ssize_t is only defined for gcc?

#include "interpolation.h"
#include "neighbours.h"

// A grid is a 3 (or 4) dimensional object that for a given index, e.g.,
// [i][j][k], contains the (linear) index into a second ArrayVector object.

const size_t default_n4[4] = { 0,0,0,0 };

/*! \brief A class holding a 4D grid which maps into an internal ArrayVector

The MapGrid4 holds a mapping grid where each entry is either a valid index into
a held ArrayVector or is `-1` to indicate that the grid point does not map.
*/
template<class T> class MapGrid4{
protected:
  size_t N[4];                 //!< The number of points along each axis of the grid
  size_t span[4];              //!< The span along each axis, allowing conversion from subscripted to linear indexing
  slong *map;                  //!< The mapping grid
  ArrayVector<T> data;         //!< The stored ArrayVector indexed by `map`
  ArrayVector<size_t> shape;   //!< A second ArrayVector to indicate a possible higher-dimensional shape of each `data` array
public:
  // constructors
  MapGrid4(const size_t *n=default_n4): map(nullptr), data(0,0), shape(1,0)
    { this->set_size(n); };
  MapGrid4(const size_t *n, const ArrayVector<T>& av): map(nullptr)
    { this->set_size(n); this->replace_data(av); };
  MapGrid4(const size_t *n, const slong *inmap, const ArrayVector<T>& av): map(nullptr)
    { this->set_size(n); this->replace_data(av); this->set_map(inmap,n,4u); };
  // copy constructor
  MapGrid4(const MapGrid4<T>& other): map(nullptr) {
    this->resize(other.size(0),other.size(1),other.size(2),other.size(3)); // sets N, calculates span, frees/allocates map memory if necessary
    for (size_t i=0; i<other.numel(); i++) this->map[i] = other.map[i];
    this->data = other.data;
    this->shape= other.shape;
  }
  // destructor
  ~MapGrid4(){
    if ( numel()>0 && map!=nullptr) delete[] map;
  }; // everything else is handled by itself
  // Copy constructor:
  // Assignment operator:
  MapGrid4<T>& operator=(const MapGrid4<T> &other){
    if (this != &other){
      this->resize(other.size(0),other.size(1),other.size(2),other.size(3)); // sets N, calculates span, frees/allocates map memory if necessary
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
  //
  //! Calculate the linear index of a point given its four subscripted indices
  size_t sub2lin(const size_t i, const size_t j, const size_t k, const size_t l) const;
  //! Calculate the linear index of a point given an array of its four subscripted indices
  int sub2lin(const size_t *s, size_t *l) const;
  //! Calculate the subscripted indices of a point given its linear index
  int lin2sub(const size_t  l, size_t *s) const;
  //! Find the mapping of a point given an array of its three subscripted indices
  int sub2map(const size_t* s, size_t& m) const;
  //! Find the mapping of a point given its linear index
  int lin2map(const size_t  l, size_t& m) const;
  //
  //! Return the total number of elements in the mapping grid
  size_t numel(void) const;
  //! Return the number of elements in the mapping grid along a given dimension
  size_t size(const size_t i) const;
  // //! Return the span of elements neccessary to move between neighbouring points
  // size_t span(const size_t i) const;
  //! Return the stride (in bytes) necessary to move between neighbouring points
  size_t stride(const size_t i) const;
  //! Change the size of the mapping grid given four new sizes
  size_t resize(const size_t n0, const size_t n1, const size_t n2, const size_t n3);
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
protected:
  void set_size(const size_t *n);
  void calc_span();
  void instantiate_map();
  bool valid_mapping(const size_t l) const;
  bool valid_mapping(const size_t i, const size_t j, const size_t k, const size_t l) const;
  bool is_inbounds(const size_t i, const size_t j, const size_t k, const size_t l) const;
  bool is_inbounds(const size_t* s) const;
};

// TODO: allow for more-compact data arrays with appropriate reducing maps

const double default_zero4[4] = {0.,0.,0.,0.};
const double default_step4[4] = {1.,1.,1.,1.};

/*! \brief Extends the MapGrid4 class to have positions of each grid point enabling interpolation between mapped points
*/
template<class T> class InterpolateGrid4: public MapGrid4<T>{
  double zero[4]; //!< the 4-vector position of `map[0]`
  double step[4]; //!< the step size along each direction of the grid
public:
  InterpolateGrid4(const size_t *n=default_n, const double *z=default_zero4, const double *s=default_step4): MapGrid4<T>(n) { this->set_zero(z); this->set_step(s); };
  InterpolateGrid4(const size_t *n, const ArrayVector<T>& av, const double *z=default_zero4, const double *s=default_step4): MapGrid4<T>(n,av){ this->set_zero(z); this->set_step(s); };
  InterpolateGrid4(const size_t *n, const slong* inmap, const ArrayVector<T>& av, const double *z=default_zero4, const double *s=default_step4): MapGrid4<T>(n,inmap,av){ this->set_zero(z); this->set_step(s); };


  void set_zero(const double *newzero){ for(int i=0;i<4;i++) this->zero[i] = newzero[i]; };
  void set_step(const double *newstep){ for(int i=0;i<4;i++) this->step[i] = newstep[i]; };

  //! Return the first three components of the positions of all grid points
  ArrayVector<double> get_grid_xyz() const {
    ArrayVector<double> xyz(3u,this->numel());
    size_t i,j,k,l;
    size_t n0=this->size(0), n1=this->size(1), n2=this->size(2), n3=this->size(3);
    size_t cnt = 0;
    double s0=this->step[0], z0=this->zero[0], s1=this->step[1], z1=this->zero[1], s2=this->step[2], z2=this->zero[2];
    for (i=0; i<n0; i++)
      for (j=0; j<n1; j++)
        for (k=0; k<n2; k++)
          for (l=0; l<n3; l++){ // yes we're looping over the fourth grid variable and not outputting it -- this is intentional.
            xyz.insert( z0+s0*i, cnt  ,0);
            xyz.insert( z1+s1*j, cnt  ,1);
            xyz.insert( z2+s2*k, cnt++,2);
          }
    return xyz;
  };
  //! Return all four components of the positions of all grid points
  ArrayVector<double> get_grid_xyzw() const {
    ArrayVector<double> xyzw(4u,this->numel());
    size_t i,j,k,l;
    size_t n0=this->size(0), n1=this->size(1), n2=this->size(2), n3=this->size(3);
    size_t cnt = 0;
    double s0=this->step[0], z0=this->zero[0], s1=this->step[1], z1=this->zero[1], s2=this->step[2], z2=this->zero[2], s3=this->step[3], z3=this->zero[3];
    for (i=0; i<n0; i++)
      for (j=0; j<n1; j++)
        for (k=0; k<n2; k++)
          for (l=0; l<n3; l++){
            xyzw.insert( z0+s0*i, cnt  ,0);
            xyzw.insert( z1+s1*j, cnt  ,1);
            xyzw.insert( z2+s2*k, cnt  ,2);
            xyzw.insert( z3+s3*l, cnt++,3);
          }
    return xyzw;
  };
  //! Return the first three components of the positions of all mapped grid points
  ArrayVector<double> get_mapped_xyz() const {
    ArrayVector<double> xyz(3u,this->valid_mapping_count());
    size_t i,j,k,l;
    size_t n0=this->size(0), n1=this->size(1), n2=this->size(2), n3=this->size(3);
    size_t cnt = 0;
    double s0=this->step[0], z0=this->zero[0], s1=this->step[1], z1=this->zero[1], s2=this->step[2], z2=this->zero[2];
    for (i=0; i<n0; i++)
      for (j=0; j<n1; j++)
        for (k=0; k<n2; k++)
          for (l=0; l<n3; l++){
            if ( this->valid_mapping(i,j,k,l) ) {
              xyz.insert( z0+s0*i, cnt  ,0);
              xyz.insert( z1+s1*j, cnt  ,1);
              xyz.insert( z2+s2*k, cnt++,2);
            }
          }
    return xyz;
  };
  //! Return all four components of the positions of all mapped grid points
  ArrayVector<double> get_mapped_xyzw() const {
    ArrayVector<double> xyzw(4u,this->valid_mapping_count());
    size_t i,j,k,l;
    size_t n0=this->size(0), n1=this->size(1), n2=this->size(2), n3=this->size(3);
    size_t cnt = 0;
    double s0=this->step[0], z0=this->zero[0], s1=this->step[1], z1=this->zero[1], s2=this->step[2], z2=this->zero[2], s3=this->step[3], z3=this->zero[3];
    for (i=0; i<n0; i++)
      for (j=0; j<n1; j++)
        for (k=0; k<n2; k++)
          for (l=0; l<n3; l++){
            if ( this->valid_mapping(i,j,k,l) ) {
              xyzw.insert( z0+s0*i, cnt  ,0);
              xyzw.insert( z1+s1*j, cnt  ,1);
              xyzw.insert( z2+s2*k, cnt  ,2);
              xyzw.insert( z3+s3*l, cnt++,3);
            }
          }
    return xyzw;
  };
  //! Return an ArrayVector of the grid coordinates along its first dimension
  ArrayVector<double> get_grid_x() const {
    ArrayVector<double> x(1u,this->size(0));
    double s0=this->step[0], z0=this->zero[0];
    for (size_t i=0; i<x.size(); i++) x.insert( z0+s0*i, i);
    return x;
  };
  //! Return an ArrayVector of the grid coordinates along its second dimension
  ArrayVector<double> get_grid_y() const {
    ArrayVector<double> y(1u,this->size(1));
    double s1=this->step[1], z1=this->zero[1];
    for (size_t j=0; j<y.size(); j++) y.insert( z1+s1*j, j);
    return y;
  };
  //! Return an ArrayVector of the grid coordinates along its third dimension
  ArrayVector<double> get_grid_z() const {
    ArrayVector<double> z(1u,this->size(2));
    double s2=this->step[2], z2=this->zero[2];
    for (size_t k=0; k<z.size(); k++) z.insert( z2+s2*k, k);
    return z;
  };
  //! Return an ArrayVector of the grid coordinates along its fourth dimension
  ArrayVector<double> get_grid_w() const {
    ArrayVector<double> w(1u,this->size(3));
    double s3=this->step[3], z3=this->zero[3];
    for (size_t l=0; l<w.size(); l++) w.insert( z3+s3*l, l);
    return w;
  };
  /*! Find the subscripted indices of the closest grid point to a specified position
  @param x The test position
  @param[out] ijkl A storage location for the subscripted indices
  @returns An integer with detailed information about if and how `x` is out of
           bounds of the grid.
  @note `out` encodes details about how each component of `x` is located in
         relation to the boundaries of the grid, using its bits as flags.
         For a given axis n, a number fₙ takes one of four
         values {0,2ⁿ,2ⁿ⁺⁴,2ⁿ⁺⁸} to indicate that  `x` is between two grid points,
         within machine precision of a grid point, smaller than the lowest grid
         point, or larger than the highest grid point, respectively.
         `out` is then f₀ + f₁ + f₂ + f₃.
  */
  unsigned int nearest_index(const double *x, size_t *ijkl) const {
    // out contains detailed information about how x is out of bounds.
    unsigned int out=0;
    slong tmp;
    for (int i=0; i<4; i++){
      tmp = (slong)( round( (x[i] - this->zero[i])/this->step[i] ) );
      if (tmp>=0 && static_cast<size_t>(tmp)<this->size(i)){
        if (approx_scalar(this->step[i]*tmp + this->zero[i],x[i]))
          out += 1<<i; // exact match
      } else {
        if (tmp<0) {
          tmp = 0;
          out += 1<<(4+i); // underflow
        } else {
          tmp = this->size(i)-1;
          out += 1<<(8+i); // overflow
        }
      }
      ijkl[i] = (size_t)(tmp);
    }
    return out;
  };
  //! Perform sanity checks before attempting to interpolate
  template<typename R> int check_before_interpolating(const ArrayVector<R>& x){
    if (this->size(0)<2||this->size(1)<2||this->size(2)<2||this->size(3)<2)
      throw std::runtime_error("Interpolation is only possible on grids with at least two elements in each dimension");
    if (this->data.size()==0)
      throw std::runtime_error("The grid must be filled before interpolating!");
    if (x.numel()!=3u)
      throw std::runtime_error("InterpolateGrid4 requires x values which are three-vectors.");
    return 0;
  };
  /*! Perform linear interpolation at the specified points expressed in an orthonormal frame
  @param x The coordinates to interpolate at expressed in the same orthonormal frame as the mapping grid
  @returns An ArrayVector of the itnerpolated values
  @note In the event that one or more coordinates of a vector in `x` is an exact
        match for a grid point this routine will perform a lower-dimensional
        interpolation. For four exact matches the point in x *is* a grid point
        and no interpolation is performed, for three exact matches 1D linear
        interpolation is performed, for two exact match bilinear 2D interpolation
        is used, for one exact match trilinear interpolation is used, and for
        no exact matches the method uses quadralinear interpolation.
  */
  template<typename R> ArrayVector<T> linear_interpolate_at(const ArrayVector<R>& x){
    this->check_before_interpolating(x);
    ArrayVector<T> out(this->data.numel(), x.size());
    size_t corners[16], ijk[4], cnt;
    unsigned int flg;
    int oob;
    double weights[16];
    std::vector<size_t> dirs, corner_count={1u,2u,4u,8u,16u};
    //TODO: switch this to an omp for loop
    for (size_t i=0; i<x.size(); i++){
      // find the closest grid subscripted indices to x[i]
      flg = this->nearest_index(x.datapointer(i), ijk );
      cnt = 1u;
      if (flg > 16){
        std::string msg_flg = "Unsure what to do with flg = " + std::to_string(flg);
        throw std::runtime_error(msg_flg);
      }
      if (15==flg)/*++++*/{
        this->sub2map(ijk,corners[0]);
        weights[0]=1.0;
      } else {
        if (0==flg)/*xxxx*/{
          dirs.resize(4);
          dirs[0]=0u; dirs[1]=1u; dirs[2]=2u; dirs[3]=3u;
        }
        if ( 1==flg|| 2==flg|| 4==flg|| 8==flg) dirs.resize(3);
        if ( 1==flg)/*+xxx*/{ dirs[0]=1u; dirs[1]=2u; dirs[2]=3u; }
        if ( 2==flg)/*x+xx*/{ dirs[0]=0u; dirs[1]=2u; dirs[2]=3u; }
        if ( 4==flg)/*xx+x*/{ dirs[0]=0u; dirs[1]=1u; dirs[2]=3u; }
        if ( 8==flg)/*xxx+*/{ dirs[0]=0u; dirs[1]=1u; dirs[2]=2u; }

        if ( 3==flg|| 5==flg|| 9==flg|| 6==flg|| 10==flg|| 12==flg) dirs.resize(2);
        if ( 3==flg)/*++xx*/{ dirs[0]=2u; dirs[1]=3u; }
        if ( 5==flg)/*+x+x*/{ dirs[0]=1u; dirs[1]=3u; }
        if ( 9==flg)/*+xx+*/{ dirs[0]=1u; dirs[1]=2u; }
        if ( 6==flg)/*x++x*/{ dirs[0]=0u; dirs[1]=3u; }
        if (10==flg)/*x+x+*/{ dirs[0]=0u; dirs[1]=2u; }
        if (12==flg)/*xx++*/{ dirs[0]=0u; dirs[1]=1u; }

        if ( 7==flg|| 11==flg || 13==flg || 14==flg) dirs.resize(1);
        if ( 7==flg)/*+++x*/ dirs[0]=3u;
        if (11==flg)/*++x+*/ dirs[0]=2u;
        if (13==flg)/*+x++*/ dirs[0]=1u;
        if (14==flg)/*x+++*/ dirs[0]=0u;

        // determine the linear indices for the (up to) 16 grid points
        // surrounding x[i] plus their linear-interpolation weights.
        oob = corners_and_weights(this,this->zero,this->step,ijk,x.datapointer(i),corners,weights,4u,dirs);
        cnt = corner_count[dirs.size()];
        if (oob) {
          std::string msg = "Point " + std::to_string(i) + " with x = " + x.to_string(i) + " has " + std::to_string(oob) + " corners out of bounds!";
          throw std::runtime_error(msg);
        }
      }
      unsafe_accumulate_to(this->data,cnt,corners,weights,out,i);
    }
    return out;
  };
  /*! Perform linear interpolation in parallel at the specified points expressed in an orthonormal frame
  @param x The coordinates to interpolate at expressed in the same orthonormal frame as the mapping grid
  @returns An ArrayVector of the itnerpolated values
  @note In the event that one or more coordinates of a vector in `x` is an exact
        match for a grid point this routine will perform a lower-dimensional
        interpolation. For four exact matches the point in x *is* a grid point
        and no interpolation is performed, for three exact matches 1D linear
        interpolation is performed, for two exact match bilinear 2D interpolation
        is used, for one exact match trilinear interpolation is used, and for
        no exact matches the method uses quadralinear interpolation.
  */
  template<typename R> ArrayVector<T> parallel_linear_interpolate_at(const ArrayVector<R>& x,const int threads){
    this->check_before_interpolating(x);
    ArrayVector<T> out(this->data.numel(), x.size());
    size_t corners[16], ijk[4], cnt;
    unsigned int flg;
    int oob;
    double weights[16];
    std::vector<size_t> dirs, corner_count={1u,2u,4u,8u,16u};

    (threads > 0 ) ? omp_set_num_threads(threads) : omp_set_num_threads(omp_get_max_threads());
    slong xsize = unsigned_to_signed<slong,size_t>(x.size());
#pragma omp parallel for shared(x,out,corner_count) firstprivate(corners,ijk,weights,xsize) private(flg,oob,cnt,dirs)
    for (slong si=0; si<xsize; si++){
      size_t i = signed_to_unsigned<size_t,slong>(si);
      // find the closest grid subscripted indices to x[i]
      flg = this->nearest_index(x.datapointer(i), ijk );
      cnt = 1u;
      if (flg > 16){
        std::string msg_flg = "Unsure what to do with flg = " + std::to_string(flg);
        throw std::runtime_error(msg_flg);
      }
      if (15==flg)/*++++*/{
        this->sub2map(ijk,corners[0]);
        weights[0]=1.0;
      } else {
        if (0==flg)/*xxxx*/{
          dirs.resize(4);
          dirs[0]=0u; dirs[1]=1u; dirs[2]=2u; dirs[3]=3u;
        }
        if ( 1==flg|| 2==flg|| 4==flg|| 8==flg) dirs.resize(3);
        if ( 1==flg)/*+xxx*/{ dirs[0]=1u; dirs[1]=2u; dirs[2]=3u; }
        if ( 2==flg)/*x+xx*/{ dirs[0]=0u; dirs[1]=2u; dirs[2]=3u; }
        if ( 4==flg)/*xx+x*/{ dirs[0]=0u; dirs[1]=1u; dirs[2]=3u; }
        if ( 8==flg)/*xxx+*/{ dirs[0]=0u; dirs[1]=1u; dirs[2]=2u; }

        if ( 3==flg|| 5==flg|| 9==flg|| 6==flg|| 10==flg|| 12==flg) dirs.resize(2);
        if ( 3==flg)/*++xx*/{ dirs[0]=2u; dirs[1]=3u; }
        if ( 5==flg)/*+x+x*/{ dirs[0]=1u; dirs[1]=3u; }
        if ( 9==flg)/*+xx+*/{ dirs[0]=1u; dirs[1]=2u; }
        if ( 6==flg)/*x++x*/{ dirs[0]=0u; dirs[1]=3u; }
        if (10==flg)/*x+x+*/{ dirs[0]=0u; dirs[1]=2u; }
        if (12==flg)/*xx++*/{ dirs[0]=0u; dirs[1]=1u; }

        if ( 7==flg|| 11==flg || 13==flg || 14==flg) dirs.resize(1);
        if ( 7==flg)/*+++x*/ dirs[0]=3u;
        if (11==flg)/*++x+*/ dirs[0]=2u;
        if (13==flg)/*+x++*/ dirs[0]=1u;
        if (14==flg)/*x+++*/ dirs[0]=0u;

        oob = corners_and_weights(this,this->zero,this->step,ijk,x.datapointer(i),corners,weights,4u,dirs);
        cnt = corner_count[dirs.size()];
        if (oob) {
          std::string msg = "Point " + std::to_string(i) + " with x = " + x.to_string(i) + " has " + std::to_string(oob) + " corners out of bounds!";
          throw std::runtime_error(msg);
        }
      }
      unsafe_accumulate_to(this->data,cnt,corners,weights,out,i);
    }
    return out;
  };
  /*! Get the size information about the first three components of the grid,
      suitable for use in creating an idential object.
      @returns An ArrayVector with three 1-element arrays containing (N[0:2]-1)/2
  */
  ArrayVector<size_t> get_halfN(void) const {
    ArrayVector<size_t> out(1u,3u,this->N); // this is the Q part of N
    return (out-1)/2; // and we want just half of it
  };
  /*! Get the fourth-dimension specification of the grid,
      suitable for use in creating an identical object.
      @returns an ArrayVector with three 1-element arrays containing {zero[3],step[3],zero[3]+(N[3]-1)*step[3]}
  */
  ArrayVector<double> get_spec(void) const {
    double spec[3];
    spec[0] = this->zero[3];
    spec[1] = this->step[3];
    spec[2] = this->zero[3] + this->step[3] * (double)(this->size(3)-1);
    ArrayVector<double> out(1u,3u,spec);
    return out;
  };
protected:
  /*! Determine the neighbouring grid points of a given grid linear index
  @param centre The linear index to a point in the mapping grid
  @returns Up to 80 neighbouring linear indices, skipping points that are out of the grid
  */
  ArrayVector<size_t> get_neighbours(const size_t centre) const {
    ArrayVector<int> mzp = make_relative_neighbour_indices4(1); // all combinations of [-1,0,+1] for four dimensions, skipping (0,0,0,0)
    ArrayVector<size_t> ijk(4u,1u);
    this->lin2sub(centre, ijk.datapointer(0)); // get the subscripted indices of the centre position
    bool isz[4], ism[4]; // is the centre index 0 (isz) or the maximum (ism)
    for (size_t i=0; i<4u; ++i) isz[i] = 0==ijk.getvalue(0,i);
    for (size_t i=0; i<4u; ++i) ism[i] = this->size(i)-1 <= ijk.getvalue(0,i);
    ArrayVector<bool> is_valid(1u,mzp.size());
    for (size_t i=0; i<mzp.size(); ++i){
      // keep track of if we *can* (or should) add each mzp vector to the centre index
      is_valid.insert(true,i);
      for (size_t j=0; j<mzp.numel(); ++j){
        if (isz[j] && mzp.getvalue(i,j)<0 ) is_valid.insert(false,i);
        if (ism[j] && mzp.getvalue(i,j)>0 ) is_valid.insert(false,i);
      }
    }
    ArrayVector<size_t> tmp(4u,1u);
    for (size_t i=0; i<mzp.size(); ++i){
      if (is_valid.getvalue(i)){
        for (size_t j=0; j<4u; ++j) tmp.insert( ijk.getvalue(0,j) + mzp.getvalue(i,j), 0, j);
        is_valid.insert( this->is_inbounds(tmp.datapointer(0)) ,i); //ensure we only check in-bounds neighbours
      }
    }
    size_t valid_neighbours = 0;
    for (size_t i=0; i<is_valid.size(); ++i) if (is_valid.getvalue(i)) ++valid_neighbours;
    ArrayVector<size_t> neighbours(1u,valid_neighbours);
    int oob = 0;
    size_t valid_neighbour=0;
    for (size_t i=0; i<mzp.size(); ++i){
      if (is_valid.getvalue(i)){
        // we can't use
        //    tmp = mzp[i] + ijk;
        // because the compiler doesn't know what to do with ArrayVector<int> + ArrayVector<size_t>
        for (size_t j=0; j<4u; ++j) tmp.insert( ijk.getvalue(0,j) + mzp.getvalue(i,j), 0, j);
        oob += this->sub2lin(tmp.datapointer(0),neighbours.datapointer(valid_neighbour++));
      }
    }
    if (oob) throw std::runtime_error("Out-of-bounds points found when there should be none.");
    return neighbours;
  };
};

#include "grid4.hpp"

#endif
