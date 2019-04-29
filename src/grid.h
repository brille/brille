// #ifdef _WIN32
  typedef long slong; // ssize_t is only defined for gcc?
// #endif
#include "arrayvector.h"
#include "neighbours.h"

#ifndef _GRID_H_
#define _GRID_H_
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

class InterpolateGrid3: public MapGrid3<double>{
  double zero[3];
  double step[3];
public:
  InterpolateGrid3(const size_t *n=default_n, const double *z=default_zero, const double *s=default_step): MapGrid3(n) { this->set_zero(z); this->set_step(s); };
  InterpolateGrid3(const size_t *n, const ArrayVector<double>& av, const double *z=default_zero, const double *s=default_step): MapGrid3(n,av){ this->set_zero(z); this->set_step(s); };
  InterpolateGrid3(const size_t *n, const slong* inmap, const ArrayVector<double>& av, const double *z=default_zero, const double *s=default_step): MapGrid3(n,inmap,av){ this->set_zero(z); this->set_step(s); };


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
    int tmp, out=0;
    for (int i=0; i<3; i++){
      tmp = round( (x[i] - this->zero[i])/this->step[i] );
      if (tmp < 0) { tmp = 0; out += 1<<i; }
      if (tmp >= this->size(i) ) { tmp = this->size(i)-1; out += 1<<(3+i); }
      ijk[i] = (size_t)tmp;
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

  template<typename R> ArrayVector<double> linear_interpolate_at(const ArrayVector<R>& x){
    if (this->data.size()==0)
      throw std::runtime_error("The grid must be filled before interpolating!");
    if (x.numel()!=3u)
      throw std::runtime_error("InterpolateGrid3 requires x values which are three-vectors.");
    ArrayVector<double> out(this->data.numel(), x.size());
    size_t corners[8], ijk[3];
    int flg, oob;
    ArrayVector<double> weights(1u, 8u);
    //TODO: switch this to an omp for loop
    for (size_t i=0; i<x.size(); i++){
      // find the closest grid subscripted indices to x[i]
      flg = this->nearest_index(x.datapointer(i), ijk );
      // determine the linear indices for the 8 grid points surrounding x[i]
      // plus their linear-interpolation weights.
      oob = this->get_corners_and_weights(corners,weights.datapointer(0),ijk,x.datapointer(i));
      if (flg || oob) {
        // flg has detailed information of the nearest grid point to ijk deficiencies
        // oob contains the number of corners which are out of bounds.
        printf("linear_interpolate_at: %d corners are out of bounds!\n",oob);
        throw std::runtime_error("out of bounds corner(s)");
      }
      // now do the actual interpolation:
      // extract an ArrayVector(this->data.numel(),8u) of the corner Arrays
      // multiply all elements at each corner by the weight for that corner
      // sum over the corners, returning an ArrayVector(this->data.numel(),1u)
      // and set that ArrayVector as element i of the output ArrayVector
      out.set( i, (this->data.extract(8u, corners) * weights).sum() );
    }
    return out;
  };
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
                  oob += this->sub2map(t,c   ); w[0] = m[0]*m[1]*m[2]; // (000)
    t[0] += d[0]; oob += this->sub2map(t,c+1u); w[1] = p[0]*m[1]*m[2]; // (100)
    t[1] += d[1]; oob += this->sub2map(t,c+2u); w[2] = p[0]*p[1]*m[2]; // (110)
    t[0] -= d[0]; oob += this->sub2map(t,c+3u); w[3] = m[0]*p[1]*m[2]; // (010)
    t[2] += d[2]; oob += this->sub2map(t,c+4u); w[4] = m[0]*p[1]*p[2]; // (011)
    t[0] += d[0]; oob += this->sub2map(t,c+5u); w[5] = p[0]*p[1]*p[2]; // (111)
    t[1] -= d[1]; oob += this->sub2map(t,c+6u); w[6] = p[0]*m[1]*p[2]; // (101)
    t[0] -= d[0]; oob += this->sub2map(t,c+7u); w[7] = m[0]*m[1]*p[2]; // (001)

    return oob;
  };
  ArrayVector<size_t> get_neighbours(const size_t centre) const {
    ArrayVector<int> mzp = make_relative_neighbour_indices(1); // all combinations of [-1,0,+1] for three dimensions, skipping (0,0,0)
    ArrayVector<size_t> ijk(3u,1u);
    this->lin2sub(centre, ijk.datapointer(0)); // get the subscripted indices of the centre position
    bool isz[3];
    for (size_t i=0; i<3u; ++i) isz[i] = 0==ijk.getvalue(0,i);
    ArrayVector<bool> is_valid(1u,mzp.size());
    for (size_t i=0; i<mzp.size(); ++i){
      is_valid.insert(true,i);
      for (size_t j=0; j<mzp.numel(); ++j)
        if (isz[j] && mzp.getvalue(i,j)<0 ) is_valid.insert(false,i);
    }
    ArrayVector<size_t> tmp(3u,1u);
    for (size_t i=0; i<mzp.size(); ++i){
      if (is_valid.getvalue(i)){
        for (size_t j=0; j<3u; ++j) tmp.insert( ijk.getvalue(0,j) + mzp.getvalue(i,j), 0, j);
        is_valid.insert( is_inbounds(tmp.datapointer(0)) ,i); //ensure we only check in-bounds neighbours
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
        for (size_t j=0; j<3u; ++j) tmp.insert( ijk.getvalue(0,j) + mzp.getvalue(i,j), 0, j);
        oob += this->sub2lin(tmp.datapointer(0),neighbours.datapointer(valid_neighbour++));
      }
    }
    if (oob) throw std::runtime_error("Out-of-bounds points found when there should be none.");
    return neighbours;
  };
};

#include "grid.hpp"

#endif
