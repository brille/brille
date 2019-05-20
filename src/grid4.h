#ifndef _GRID4_H_
#define _GRID4_H_
// #ifdef _WIN32
  typedef long slong; // ssize_t is only defined for gcc?
// #endif
#include "arrayvector.h"
#include "interpolation.h"
#include "neighbours.h"

// A grid is a 3 (or 4) dimensional object that for a given index, e.g.,
// [i][j][k], contains the (linear) index into a second ArrayVector object.

const size_t default_n4[4] = { 0,0,0,0 };

template<class T> class MapGrid4{
protected:
  size_t N[4];
  size_t span[4];
  slong *map; // replace this with a Blitz++ multidimensional array? https://github.com/blitzpp/blitz
  ArrayVector<T> data;
  ArrayVector<size_t> shape; // when "filled" (by a Python call), what is the shape of the numpy.np array that is passed in
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
  size_t sub2lin(const size_t i, const size_t j, const size_t k, const size_t l) const;
  int sub2lin(const size_t *s, size_t *l) const;
  int lin2sub(const size_t  l, size_t *s) const;
  int sub2map(const size_t *s, size_t *m) const;
  int lin2map(const size_t  l, size_t *m) const;
  //
  size_t numel(void) const;
  size_t size(const size_t i) const;
  //
  size_t resize(const size_t n0, const size_t n1, const size_t n2, const size_t n3);
  size_t resize(const size_t *n);
  //
  size_t data_ndim(void) const;
  size_t num_data(void) const;
  ArrayVector<size_t> data_shape(void) const;
  ArrayVector<size_t> get_N(void) const;
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

template<class T> class InterpolateGrid4: public MapGrid4<T>{
  double zero[4];
  double step[4];
public:
  InterpolateGrid4(const size_t *n=default_n, const double *z=default_zero4, const double *s=default_step4): MapGrid4<T>(n) { this->set_zero(z); this->set_step(s); };
  InterpolateGrid4(const size_t *n, const ArrayVector<T>& av, const double *z=default_zero4, const double *s=default_step4): MapGrid4<T>(n,av){ this->set_zero(z); this->set_step(s); };
  InterpolateGrid4(const size_t *n, const slong* inmap, const ArrayVector<T>& av, const double *z=default_zero4, const double *s=default_step4): MapGrid4<T>(n,inmap,av){ this->set_zero(z); this->set_step(s); };


  void set_zero(const double *newzero){ for(int i=0;i<4;i++) this->zero[i] = newzero[i]; };
  void set_step(const double *newstep){ for(int i=0;i<4;i++) this->step[i] = newstep[i]; };

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
  ArrayVector<double> get_grid_w() const {
    ArrayVector<double> w(1u,this->size(3));
    double s3=this->step[3], z3=this->zero[3];
    for (size_t l=0; l<w.size(); l++) w.insert( z3+s3*l, l);
    return w;
  };

  int nearest_index(const double *x, size_t *ijkl) const {
    // out contains detailed information about how x is out of bounds.
    int tmp, out=0;
    for (int i=0; i<4; i++){
      tmp = round( (x[i] - this->zero[i])/this->step[i] );
      if (tmp < 0) { tmp = 0; out += 1<<i; }
      if (tmp >= this->size(i) ) { tmp = this->size(i)-1; out += 1<<(4+i); }
      ijkl[i] = (size_t)tmp;
    }
    return out;
  };
  size_t nearest_linear_index(const double *x) const {
    size_t ijkl[3];
    int ret;
    if ( (ret = nearest_index(x,ijkl)) ){
      // somehow indicate that we've hit one or more boundaries?
    }
    size_t lidx;
    ret += this->sub2lin(ijkl,&lidx);
    return lidx;
  };
  template<typename R> int check_before_interpolating(const ArrayVector<R>& x){
    if (this->size(0)<2||this->size(1)<2||this->size(2)<2||this->size(3)<2)
      throw std::runtime_error("Interpolation is only possible on grids with at least two elements in each dimension");
    if (this->data.size()==0)
      throw std::runtime_error("The grid must be filled before interpolating!");
    if (x.numel()!=3u)
      throw std::runtime_error("InterpolateGrid4 requires x values which are three-vectors.");
    return 0;
  };

  template<typename R> ArrayVector<T> linear_interpolate_at(const ArrayVector<R>& x){
    this->check_before_interpolating(x);
    ArrayVector<T> out(this->data.numel(), x.size());
    size_t corners[16], ijk[4], cnt;
    int flg, oob;
    double weights[16];
    std::vector<size_t> dirs, corner_count={1u,2u,4u,8u,16u};
    //TODO: switch this to an omp for loop
    for (size_t i=0; i<x.size(); i++){
      // find the closest grid subscripted indices to x[i]
      flg = this->nearest_index(x.datapointer(i), ijk );
      cnt = 1u;
      if (flg > 0 ){
        std::string msg_flg = "Unsure what to do with flg = " + std::to_string(flg);
        throw std::runtime_error(msg_flg);
      }
      if (-15 == flg)/*++++*/{
        this->sub2map(ijk,corners);
        weights[0]=1.0;
      } else {
        if (0 == flg)/*xxxx*/{
          dirs.resize(4);
          dirs[0]=0u; dirs[1]=1u; dirs[2]=2u; dirs[3]=3u;
        }
        if ( -1==flg|| -2==flg|| -4==flg|| -8==flg) dirs.resize(3);
        if ( -1==flg)/*+xxx*/{ dirs[0]=1u; dirs[1]=2u; dirs[2]=3u; }
        if ( -2==flg)/*x+xx*/{ dirs[0]=0u; dirs[1]=2u; dirs[2]=3u; }
        if ( -4==flg)/*xx+x*/{ dirs[0]=0u; dirs[1]=1u; dirs[2]=3u; }
        if ( -8==flg)/*xxx+*/{ dirs[0]=0u; dirs[1]=1u; dirs[2]=2u; }

        if ( -3==flg|| -5==flg|| -9==flg|| -6==flg||-10==flg||-12==flg) dirs.resize(2);
        if ( -3==flg)/*++xx*/{ dirs[0]=2u; dirs[1]=3u; }
        if ( -5==flg)/*+x+x*/{ dirs[0]=1u; dirs[1]=3u; }
        if ( -9==flg)/*+xx+*/{ dirs[0]=1u; dirs[1]=2u; }
        if ( -6==flg)/*x++x*/{ dirs[0]=0u; dirs[1]=3u; }
        if (-10==flg)/*x+x+*/{ dirs[0]=0u; dirs[1]=2u; }
        if (-12==flg)/*xx++*/{ dirs[0]=0u; dirs[1]=1u; }

        if ( -7==flg||-11==flg||-13==flg||-14==flg) dirs.resize(1);
        if ( -7==flg)/*+++x*/ dirs[0]=3u;
        if (-11==flg)/*++x+*/ dirs[0]=2u;
        if (-13==flg)/*+x++*/ dirs[0]=1u;
        if (-14==flg)/*x+++*/ dirs[0]=0u;

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
  template<typename R> ArrayVector<T> parallel_linear_interpolate_at(const ArrayVector<R>& x,const int threads){
    this->check_before_interpolating(x);
    ArrayVector<T> out(this->data.numel(), x.size());
    size_t corners[16], ijk[4], cnt;
    int flg, oob;
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
      if (flg > 0 ){
        std::string msg_flg = "Unsure what to do with flg = " + std::to_string(flg);
        throw std::runtime_error(msg_flg);
      }
      if (-15 == flg)/*++++*/{
        this->sub2map(ijk,corners);
        weights[0]=1.0;
      } else {
        if (0 == flg)/*xxxx*/{
          dirs.resize(4);
          dirs[0]=0u; dirs[1]=1u; dirs[2]=2u; dirs[3]=3u;
        }
        if ( -1==flg|| -2==flg|| -4==flg|| -8==flg) dirs.resize(3);
        if ( -1==flg)/*+xxx*/{ dirs[0]=1u; dirs[1]=2u; dirs[2]=3u; }
        if ( -2==flg)/*x+xx*/{ dirs[0]=0u; dirs[1]=2u; dirs[2]=3u; }
        if ( -4==flg)/*xx+x*/{ dirs[0]=0u; dirs[1]=1u; dirs[2]=3u; }
        if ( -8==flg)/*xxx+*/{ dirs[0]=0u; dirs[1]=1u; dirs[2]=2u; }

        if ( -3==flg|| -5==flg|| -9==flg|| -6==flg||-10==flg||-12==flg) dirs.resize(2);
        if ( -3==flg)/*++xx*/{ dirs[0]=2u; dirs[1]=3u; }
        if ( -5==flg)/*+x+x*/{ dirs[0]=1u; dirs[1]=3u; }
        if ( -9==flg)/*+xx+*/{ dirs[0]=1u; dirs[1]=2u; }
        if ( -6==flg)/*x++x*/{ dirs[0]=0u; dirs[1]=3u; }
        if (-10==flg)/*x+x+*/{ dirs[0]=0u; dirs[1]=2u; }
        if (-12==flg)/*xx++*/{ dirs[0]=0u; dirs[1]=1u; }

        if ( -7==flg||-11==flg||-13==flg||-14==flg) dirs.resize(1);
        if ( -7==flg)/*+++x*/ dirs[0]=3u;
        if (-11==flg)/*++x+*/ dirs[0]=2u;
        if (-13==flg)/*+x++*/ dirs[0]=1u;
        if (-14==flg)/*x+++*/ dirs[0]=0u;

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

  ArrayVector<size_t> get_halfN(void) const {
    ArrayVector<size_t> out(1u,3u,this->N); // this is the Q part of N
    return (out-1)/2; // and we want just half of it
  };
  ArrayVector<double> get_spec(void) const {
    double spec[3];
    spec[0] = this->zero[3];
    spec[1] = this->step[3];
    spec[2] = this->zero[3] + this->step[3] * (double)(this->size(3)-1);
    ArrayVector<double> out(1u,3u,spec);
    return out;
  };
protected:
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
