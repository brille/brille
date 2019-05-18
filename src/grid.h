#ifndef _GRID_H_
#define _GRID_H_
// #ifdef _WIN32
  typedef long slong; // ssize_t is only defined for gcc?
// #endif
#include "arrayvector.h"
#include "latvec.h"
#include "neighbours.h"
#include "interpolation.h"
#include <omp.h>
#include <complex>
#include <memory>

#include "unsignedtosigned.h"
#include "munkres.h"

// A grid is a 3 (or 4) dimensional object that for a given index, e.g.,
// [i][j][k], contains the (linear) index into a second ArrayVector object.

const size_t default_n[3] = { 0u,0u,0u };

template<class T> class MapGrid3{
protected:
  size_t N[3];
  size_t span[3];
  slong *map; // replace this with a Blitz++ multidimensional array? https://github.com/blitzpp/blitz
  ArrayVector<T> data;
  ArrayVector<size_t> shape; // when "filled" (by a Python call), what is the shape of the numpy.np array that is passed in
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
  void print_N(const bool nl=false) const;
  void print_span(const bool nl=false) const;
  void print_map(void) const;
  int set_map(void);
  int set_map(const slong* inmap, const size_t* n, const size_t d);
  int unsafe_set_map(slong *inmap);
  size_t unsafe_get_map(slong *outmap) const;
  //
  size_t maximum_mapping(const slong *map2check, const size_t num2check) const;
  size_t maximum_mapping(const slong *map2check) const;
  size_t maximum_mapping(void) const;
  //
  size_t valid_mapping_count(void) const;
  //
  int check_map(const ArrayVector<T>& data2check) const;
  int check_map(void) const;
  //
  int replace_data(const ArrayVector<T>& newdata, const ArrayVector<size_t>& newshape);
  int replace_data(const ArrayVector<T>& newdata);
  //
  size_t sub2lin(const size_t i, const size_t j, const size_t k) const;
  int sub2lin(const size_t *s, size_t *l) const;
  int lin2sub(const size_t  l, size_t *s) const;
  int sub2map(const size_t *s, size_t *m) const;
  size_t sub2map(const size_t *s) const;
  int lin2map(const size_t  l, size_t *m) const;
  //
  size_t numel(void) const;
  size_t size(const size_t i) const;
  //
  size_t resize(const size_t n0, const size_t n1, const size_t n2);
  size_t resize(const size_t *n);
  //
  size_t data_ndim(void) const;
  size_t num_data(void) const;
  ArrayVector<size_t> data_shape(void) const;
  ArrayVector<size_t> get_N(void) const;
  //
  ArrayVector<size_t> get_neighbours(const size_t centre) const;
  ArrayVector<size_t> sort_perm(void) const;
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

template<class T> class InterpolateGrid3: public MapGrid3<T>{
  double zero[3];
  double step[3];
public:
  InterpolateGrid3(const size_t *n=default_n, const double *z=default_zero, const double *s=default_step): MapGrid3<T>(n) { this->set_zero(z); this->set_step(s); };
  InterpolateGrid3(const size_t *n, const ArrayVector<T>& av, const double *z=default_zero, const double *s=default_step): MapGrid3<T>(n,av){ this->set_zero(z); this->set_step(s); };
  InterpolateGrid3(const size_t *n, const slong* inmap, const ArrayVector<T>& av, const double *z=default_zero, const double *s=default_step): MapGrid3<T>(n,inmap,av){ this->set_zero(z); this->set_step(s); };


  void set_zero(const double *newzero){ for(int i=0;i<3;i++) this->zero[i] = newzero[i]; };
  void set_step(const double *newstep){ for(int i=0;i<3;i++) this->step[i] = newstep[i]; };

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
  size_t get_grid_x(const size_t maxN, double *x) const {
    if ( maxN < this->size(0) ) return 0;
    size_t cnt = 0;
    double s0=this->step[0], z0=this->zero[0];
    for (size_t i=0; i<this->size(0); i++) x[ cnt++ ] = z0 + s0*i;
    return cnt;
  };
  size_t get_grid_y(const size_t maxN, double *y) const {
    if ( maxN < this->size(1) ) return 0;
    size_t cnt = 0;
    double s1=this->step[1], z1=this->zero[1];
    for (size_t j=0; j<this->size(1); j++) y[ cnt++ ] = z1 + s1*j;
    return cnt;
  };
  size_t get_grid_z(const size_t maxN, double *z) const {
    if ( maxN < this->size(2) ) return 0;
    size_t cnt = 0;
    double s2=this->step[2], z2=this->zero[2];
    for (size_t k=0; k<this->size(2); k++) z[ cnt++ ] = z2 + s2*k;
    return cnt;
  };

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
  ArrayVector<double> get_grid_x() const {
    ArrayVector<double> x(1u,this->size(0));
    double s0=this->step[0], z0=this->zero[0];
    for (size_t i=0; i<x.size(); i++) x.insert( z0+s0*i, i);
    return x;
  };
  ArrayVector<double> get_grid_y() const {
    ArrayVector<double> y(1u,this->size(1));
    double s1=this->step[1], z1=this->zero[1];
    for (size_t j=0; j<y.size(); j++) y.insert( z1+s1*j, j);
    return y;
  };
  ArrayVector<double> get_grid_z() const {
    ArrayVector<double> z(1u,this->size(2));
    double s2=this->step[2], z2=this->zero[2];
    for (size_t k=0; k<z.size(); k++) z.insert( z2+s2*k, k);
    return z;
  };

  int nearest_index(const double *x, size_t *ijk) const {
    /*
    out contains detailed information about how x is out of bounds.
    out = f0 + f1 + f2, where f0 = {1,0,8}, f1 = {2,0,16}, f2 ={4,0,32}
    -- that is, {2^n,0,2^(n+3)} for n=0,1,2 -- the first value
    indicates an index below 0 in that dimension, the last value
    indicates an index above the maximum in that dimension, and a zero
    for an in-bounds index.
    possible values:
      (000) => 000000 ==  0  (-00) => 000001 ==  1  (0-0) => 000010 ==  2
      (--0) => 000011 ==  3  (00-) => 000100 ==  4  (-0-) => 000101 ==  5
      (0--) => 000110 ==  6  (---) => 000111 ==  7  (+00) => 001000 ==  8
      (+-0) => 001010 == 10  (+0-) => 001100 == 12  (+--) => 001110 == 14
      (0+0) => 010000 == 16  (-+0) => 010001 == 17  (0+-) => 010100 == 20
      (-+-) => 010101 == 21  (++0) => 011000 == 24  (++-) => 011100 == 28
      (00+) => 100000 == 32  (-0+) => 100001 == 33  (0-+) => 100010 == 34
      (--+) => 100011 == 35  (+0+) => 101000 == 40  (+-+) => 101010 == 42
      (0++) => 110000 == 48  (-++) => 110001 == 49  (+++) => 111000 == 56
    */
    int out=0;
    slong tmp;
    for (int i=0; i<3; i++){
      tmp = (slong)( round( (x[i] - this->zero[i])/this->step[i] ) );
      if (tmp < 0) { tmp = 0; out += 1<<i; }
      if (tmp >= this->size(i) ) { tmp = this->size(i)-1; out += 1<<(3+i); }
      if (approx_scalar(this->step[i]*tmp + this->zero[i],x[i])) { out -= 1<<i; }; // signal an exact match
      // ijk[i] = signed_to_unsigned<size_t,slong>(tmp);
      ijk[i] = (size_t)(tmp);
    }
    return out;
  };
  size_t nearest_linear_index(const double *x) const {
    size_t ijk[3];
    int ret;
    if ( (ret = nearest_index(x,ijk)) ){
      // somehow indicate that we've hit one or more boundaries?
    }
    size_t lidx;
    ret += this->sub2lin(ijk,&lidx);
    return lidx;
  }
  template<typename R> int check_before_interpolating(const ArrayVector<R>& x){
    if (this->size(0)<2||this->size(1)<2||this->size(2)<2)
      throw std::runtime_error("Interpolation is only possible on grids with at least two elements in each dimension");
    if (this->data.size()==0)
      throw std::runtime_error("The grid must be filled before interpolating!");
    if (x.numel()!=3u)
      throw std::runtime_error("InterpolateGrid3 requires x values which are three-vectors.");
    return 0;
  }
  template<typename R> ArrayVector<T> linear_interpolate_at(const LQVec<R>& x){return this->linear_interpolate_at(x.get_xyz());}
  template<typename R> ArrayVector<T> linear_interpolate_at(const LDVec<R>& x){return this->linear_interpolate_at(x.get_xyz());}

  template<typename R> ArrayVector<T> linear_interpolate_at(const ArrayVector<R>& x){
    this->check_before_interpolating(x);
    ArrayVector<T> out(this->data.numel(), x.size());
    size_t corners[8], ijk[3], cnt;
    int flg, oob=0;
    double weights[8];

    //TODO: switch this to an omp for loop
    for (size_t i=0; i<x.size(); i++){
      // find the closest grid subscripted indices to x[i]
      flg = this->nearest_index(x.datapointer(i), ijk );
      switch (flg){
        case -7: // [ijk] exact match printf("Exact match!\n");
          this->sub2map(ijk,corners); // set the first "corner" to this mapped index
          weights[0] = 1.0; // and the weight to one
          cnt = 1u;
          break;
        case -6: // [0jk] partial match        printf("Partial (0jk) match\n");
          oob = this->get_corners_and_weights_line(0u,corners,weights,ijk,x.datapointer(i));
          cnt = 2u;
          break;
        case -5: // [i0k] partial match        printf("Partial (i0k) match\n");
          oob = this->get_corners_and_weights_line(1u,corners,weights,ijk,x.datapointer(i));
          cnt = 2u;
          break;
        case -4: // [00k] partial match        printf("Partial (00k) match\n");
          oob = this->get_corners_and_weights_area(0u,1u,corners,weights,ijk,x.datapointer(i));
          cnt = 4u;
          break;
        case -3: // [ij0] partial match        printf("Partial (ij0) match\n");
          oob = this->get_corners_and_weights_line(2u,corners,weights,ijk,x.datapointer(i));
          cnt = 2u;
          break;
        case -2: // [0j0] partial match        printf("Partial (0j0) match\n");
          oob = this->get_corners_and_weights_area(0u,2u,corners,weights,ijk,x.datapointer(i));
          cnt = 4u;
          break;
        case -1: // [i00] partial match        printf("Partial (i00) match\n");
          oob = this->get_corners_and_weights_area(1u,2u,corners,weights,ijk,x.datapointer(i));
          cnt = 4u;
          break;
        case  0: // all neighbours exist and no partial exact match
          oob = this->get_corners_and_weights(corners,weights,ijk,x.datapointer(i));
          cnt = 8u;
          break;
        default:
          std::string msg_flg = "Unsure what to do with flg = " + std::to_string(flg);
          throw std::runtime_error(msg_flg);
      }
      if (oob){
        std::string msg = "Point " + std::to_string(i) + " with x = " + x.to_string(i) + " has " + std::to_string(oob) + " corners out of bounds!";
        throw std::runtime_error(msg);
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

  template<typename R> ArrayVector<T> parallel_linear_interpolate_at(const LQVec<R>& x, const int threads){return this->parallel_linear_interpolate_at(x.get_xyz(),threads);}
  template<typename R> ArrayVector<T> parallel_linear_interpolate_at(const LDVec<R>& x, const int threads){return this->parallel_linear_interpolate_at(x.get_xyz(),threads);}
  template<typename R> ArrayVector<T> parallel_linear_interpolate_at(const ArrayVector<R>& x, const int threads){
    this->check_before_interpolating(x);
    ArrayVector<T> out(this->data.numel(), x.size());

    (threads>0) ? omp_set_num_threads(threads) : omp_set_num_threads(omp_get_max_threads());

    size_t corners[8], ijk[3], cnt=0u;
    int flg=0, oob=0;
    double weights[8];
    slong xsize = unsigned_to_signed<slong,size_t>(x.size());

#pragma omp parallel for shared(x,out) firstprivate(corners,ijk,weights,xsize) private(flg,oob,cnt)
    for (slong si=0; si<xsize; si++){
      size_t i = signed_to_unsigned<size_t,slong>(si);
      // find the closest grid subscripted indices to x[i]
      flg = this->nearest_index(x.datapointer(i), ijk );
      switch (flg){
        case -7: // [ijk] exact match printf("Exact match!\n");
          this->sub2map(ijk,corners); // set the first "corner" to this mapped index
          weights[0] = 1.0; // and the weight to one
          cnt = 1u;
          break;
        case -6: // [0jk] partial match        printf("Partial (0jk) match\n");
          oob = this->get_corners_and_weights_line(0u,corners,weights,ijk,x.datapointer(i));
          cnt = 2u;
          break;
        case -5: // [i0k] partial match        printf("Partial (i0k) match\n");
          oob = this->get_corners_and_weights_line(1u,corners,weights,ijk,x.datapointer(i));
          cnt = 2u;
          break;
        case -4: // [00k] partial match        printf("Partial (00k) match\n");
          oob = this->get_corners_and_weights_area(0u,1u,corners,weights,ijk,x.datapointer(i));
          cnt = 4u;
          break;
        case -3: // [ij0] partial match        printf("Partial (ij0) match\n");
          oob = this->get_corners_and_weights_line(2u,corners,weights,ijk,x.datapointer(i));
          cnt = 2u;
          break;
        case -2: // [0j0] partial match        printf("Partial (0j0) match\n");
          oob = this->get_corners_and_weights_area(0u,2u,corners,weights,ijk,x.datapointer(i));
          cnt = 4u;
          break;
        case -1: // [i00] partial match        printf("Partial (i00) match\n");
          oob = this->get_corners_and_weights_area(1u,2u,corners,weights,ijk,x.datapointer(i));
          cnt = 4u;
          break;
        case  0: // all neighbours exist and no partial exact match
          oob = this->get_corners_and_weights(corners,weights,ijk,x.datapointer(i));
          cnt = 8u;
          break;
        default:
          std::string msg_flg = "Unsure what to do with flg = " + std::to_string(flg);
          throw std::runtime_error(msg_flg);
      }
      if (oob){
        std::string msg = "Point " + std::to_string(i) + " with x = " + x.to_string(i) + " has " + std::to_string(oob) + " corners out of bounds!";
        throw std::runtime_error(msg);
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

//   template<typename R> ArrayVector<T> parallel_linear_interpolate_at(const ArrayVector<R>& x,const int threads){
//     this->check_before_interpolating(x);
//     ArrayVector<T> out(this->data.numel(), x.size());
//     const ArrayVector<R>* const xptr = &x;
//           ArrayVector<T>* const outptr = &out;
//     const ArrayVector<T>* const datptr = &(this->data);
//     const InterpolateGrid3<T>* const thsptr = this;
//
//     (threads > 0 ) ? omp_set_num_threads(threads) : omp_set_num_threads(omp_get_max_threads());
//
//     size_t corners[8], ijk[3];
//     double weights[8];
//     // some versions of OpenMP require that the for loop variable be signed.
//     slong xsize = unsigned_to_signed<slong,size_t>(x.size());
//     // Compared to single-processor code, there is no out explicit out-of-bounds
//     // error checking here (nearest_index and get_corners_and_weights) both do
//     // internal checks, but we ignore their results.
//     // TODO: consider putting the error check back in within a pragma omp critical block
//
// #pragma omp parallel for firstprivate(corners,ijk,weights,xsize)
//     for (slong si=0; si<xsize; si++){
//       // size_t i = (size_t)(i); // all functions called with i expect an unsigned integer.
//       size_t i = signed_to_unsigned<size_t,slong>(si);
//       // find the closest grid subscripted indices to x[i]
//       thsptr->nearest_index(xptr->datapointer(i), ijk );
//       // determine the linear indices for the 8 grid points surrounding x[i]
//       // plus their linear-interpolation weights.
//       thsptr->get_corners_and_weights(corners,weights,ijk,xptr->datapointer(i));
//
//       // now do the actual interpolation:
//       // extract an ArrayVector(this->data.numel(),8u) of the corner Arrays
//       // multiply all elements at each corner by the weight for that corner
//       // sum over the corners, returning an ArrayVector(this->data.numel(),1u)
//       // and set that ArrayVector as element i of the output ArrayVector
//       unsafe_accumulate_to(*datptr, 8u, corners, weights, *outptr, i);
//     }
//     return out;
//   };

protected:
  int get_corners_and_weights(size_t *c, double *w, const size_t *ijk, const double *x) const {
    int d[3], oob=0;
    size_t t[3];
    double p[3], m[3], tmp;
    for (int i=0; i<3; i++){
      t[i]=ijk[i];
      tmp = x[i] - (this->zero[i]+ijk[i]*this->step[i]);
      d[i] = tmp < 0 ? -1 : 1;
      p[i] = abs(tmp/(double)(this->step[i]));
      m[i] = 1.0-p[i];
    }
                  oob +=     this->sub2map(t,c   ); w[0] = m[0]*m[1]*m[2]; // (000)
    t[0] += d[0]; oob +=   2*this->sub2map(t,c+1u); w[1] = p[0]*m[1]*m[2]; // (100)
    t[1] += d[1]; oob +=   4*this->sub2map(t,c+2u); w[2] = p[0]*p[1]*m[2]; // (110)
    t[0] -= d[0]; oob +=   8*this->sub2map(t,c+3u); w[3] = m[0]*p[1]*m[2]; // (010)
    t[2] += d[2]; oob +=  16*this->sub2map(t,c+4u); w[4] = m[0]*p[1]*p[2]; // (011)
    t[0] += d[0]; oob +=  32*this->sub2map(t,c+5u); w[5] = p[0]*p[1]*p[2]; // (111)
    t[1] -= d[1]; oob +=  64*this->sub2map(t,c+6u); w[6] = p[0]*m[1]*p[2]; // (101)
    t[0] -= d[0]; oob += 128*this->sub2map(t,c+7u); w[7] = m[0]*m[1]*p[2]; // (001)

    return oob;
  };
  int get_corners_and_weights_area(const size_t a, const size_t b, size_t *c, double *w, const size_t *ijk, const double *x) const {
    int d[2], oob=0;
    double p[2],m[2];
    p[0] = interpolation_direction_and_distance(this->zero[a],this->step[a],ijk[a],x[a]);
    p[1] = interpolation_direction_and_distance(this->zero[b],this->step[b],ijk[b],x[b]);
    for (int i=0; i<2; ++i){
      d[i] = p[i] < 0 ? -1 : 1;
      p[i] = std::abs(p[i]);
      m[i] = 1-p[i];
    }
    size_t t[3];
    for (int i=0; i<3; i++) t[i]=ijk[i];
                  oob +=   this->sub2map(t,c   ); w[0] = m[0]*m[1]; // (00)
    t[a] += d[0]; oob += 2*this->sub2map(t,c+1u); w[1] = p[0]*m[1]; // (10)
    t[b] += d[1]; oob += 4*this->sub2map(t,c+2u); w[2] = p[0]*p[1]; // (11)
    t[a] -= d[0]; oob += 8*this->sub2map(t,c+3u); w[3] = m[0]*p[1]; // (01)
    return oob;
  };
  int get_corners_and_weights_line(const size_t a, size_t *c, double *w, const size_t *ijk, const double *x) const {
    double p = interpolation_direction_and_distance(this->zero[a],this->step[a],ijk[a],x[a]);
    int oob, d = p < 0 ? -1 : 1;
    p = std::abs(p);
    double m = 1.0-p;
    size_t t[3];
    for (int i=0; i<3; i++) t[i]=ijk[i];
    // (0)
    oob  = this->sub2map(t,c);
    w[0] = m;
    // (1)
    t[a] += d;
    oob += 2*this->sub2map(t,c+1u);
    w[1] = p;
    return oob;
  };
};



template<class T> struct GridDiffTraits{
  using type = T;
  constexpr static T max = std::numeric_limits<T>::max();
};
template<class T> struct GridDiffTraits<std::complex<T>>{
  using type = T;
  constexpr static T max = std::numeric_limits<T>::max();
};

#include "grid.hpp"

#endif
