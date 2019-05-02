// #ifdef _WIN32
  typedef long slong; // ssize_t is only defined for gcc?
// #endif

#ifndef _BZ_GRID_
#define _BZ_GRID_

#include <limits>

#include "arrayvector.h"
#include "lattice.h"
#include "latvec.h"
#include "bz.h"
#include "grid.h"
#include "grid4.h"

// TODO: Be clever about the map/grid creation (like spglib)?

template<class T> class BrillouinZoneGrid3: public InterpolateGrid3<T>{
protected:
  BrillouinZone brillouinzone;
public:
  BrillouinZoneGrid3(const BrillouinZone bz, const double *d, const int isrlu=1): brillouinzone(bz) { this->determine_map_size(d,isrlu);};
  BrillouinZoneGrid3(const BrillouinZone bz, const size_t *n): brillouinzone(bz) { this->determine_map_step(n); };
  const BrillouinZone get_brillouinzone() const { return this->brillouinzone; };
  size_t get_grid_hkl(const size_t maxN, double *hkl) const {
    ArrayVector<double> xyz = this->get_grid_xyz();
    if ( xyz.size() > maxN ) {printf("!!!!Grid has %lu points but output only has space for %lu\n",xyz.size(),maxN); return 0;}
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!matrix_inverse(fromxyz,toxyz)) printf("transform matrix toxyz has zero determinant?!\n");
    for (size_t i=0; i<xyz.size(); i++) multiply_matrix_vector<double,double,double,3>(hkl+3*i, fromxyz, xyz.datapointer(i));
    return xyz.size();
  };
  ArrayVector<double> get_grid_hkl() const {
    ArrayVector<double> xyz = this->get_grid_xyz();
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!matrix_inverse(fromxyz,toxyz)) printf("transform matrix toxyz has zero determinant?!\n");
    ArrayVector<double> hkl(3,xyz.size());
    for (size_t i=0; i<xyz.size(); i++) multiply_matrix_vector<double,double,double,3>(hkl.datapointer(i), fromxyz, xyz.datapointer(i));
    return hkl;
  };
  ArrayVector<double> get_mapped_hkl() const {
    ArrayVector<double> xyz = this->get_mapped_xyz();
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!matrix_inverse(fromxyz,toxyz)) printf("transform matrix toxyz has zero determinant?!\n");
    ArrayVector<double> hkl(3,xyz.size());
    for (size_t i=0; i<xyz.size(); i++) multiply_matrix_vector<double,double,double,3>(hkl.datapointer(i), fromxyz, xyz.datapointer(i));
    return hkl;
  };
protected:
  void determine_map_size(const double *d_in, const int isrlu){
    double d[3];
    for (int i=0;i<3;i++) d[i]=2.0*d_in[i]; // ×2 here so we can double N later
    if (isrlu){
      LQVec<int> ei(this->brillouinzone.get_lattice(),1u);
      for (int i=0;i<3;i++){
        for (int j=0;j<3;j++) ei.insert(i==j?1:0,0,j);
        // ei.norm() == ArrayVector<double>(1u,1u) --> .getvalue(0) returns its only element
        // d[i] *= ei.norm().getvalue(0);
        d[i] *= norm(ei).getvalue(0);
      }
    }
    // d is "guaranteed" to now be in units of inverse angstrom.
    ArrayVector<double> vxyz = this->brillouinzone.get_vertices().get_xyz();
    double minx[3], maxx[3];
    for (int i=0; i<3; i++) {minx[i]= std::numeric_limits<double>::max(); maxx[i]= -std::numeric_limits<double>::max();} // +/- max
    for (size_t i=0; i<vxyz.size(); i++){
      for (size_t j=0; j<3u; j++){
        if (vxyz.getvalue(i,j) < minx[j]) minx[j] = vxyz.getvalue(i,j);
        if (vxyz.getvalue(i,j) > maxx[j]) maxx[j] = vxyz.getvalue(i,j);
      }
    }
    size_t n[3];
    for (int i=0; i<3; i++) n[i] = 2*std::ceil( (maxx[i]-minx[i])/d[i] ); // ×2 to ensure we have an even number of steps
    // account for using ceil:
    for (int i=0; i<3; i++) d[i] = (maxx[i]-minx[i])/n[i];
    // for interpolation purposes, we want to make sure we go one-step beyond
    // the zone boundary in each direction:
    for (int i=0; i<3; i++) n[i] += 2;
    // which means we need to shift the origin of the grid as well by one step
    double z[3];
    for (int i=0; i<3; i++) z[i] = 0.0 - d[i]*(n[i]/2 -1); // ensure we hit zero, and that N/2 points are on the postive side

    this->resize(n);
    this->set_step(d);
    this->set_zero(z);
    this->set_map(); // fills the 3D map with linear indices -- eventually replace this with a fancier mapping that (a) respects the zone boundary and (b) sticks to the irreducible wedge
    this->truncate_grid_to_brillouin_zone();
  };

  void determine_map_step(const size_t *in){
    size_t n[3];
    bool iszero[3];
    for (int i=0;i<3;i++) {
      iszero[i] = 0u==in[i];
      n[i]= (iszero[i]) ? in[i]+1 : in[i];
    }
    ArrayVector<double> vxyz = this->brillouinzone.get_vertices().get_xyz();
    double maxx[3];
    for (int i=0; i<3; i++) maxx[i]= -std::numeric_limits<double>::max(); // - max
    for (size_t i=0; i<vxyz.size(); i++)
      for (size_t j=0; j<3u; j++)
        if (vxyz.getvalue(i,j) > maxx[j]) maxx[j] = vxyz.getvalue(i,j);
    double d[3];
    for (int i=0; i<3; i++) d[i] = (n[i]<2) ? (iszero[i] ? 0.0 : maxx[i]) : maxx[i]/(n[i]-1);
    double z[3];
    for (int i=0; i<3; i++) z[i] = 0.0 - d[i]*(n[i]-1); // ensure we hit zero, and that N/2 points are on the postive side

    for (int i=0;i<3;i++) if (!iszero[i]) n[i]*=2; // we actually want 2N unless the input requested zero points, in which case we *need* one

    this->resize(n);
    this->set_step(d);
    this->set_zero(z);
    this->set_map();
    this->truncate_grid_to_brillouin_zone();
  };

  void truncate_grid_to_brillouin_zone(){
    LQVec<double> hkl(this->brillouinzone.get_lattice(),this->get_grid_hkl());
    ArrayVector<bool> inbz = this->brillouinzone.isinside(&hkl);
    ArrayVector<bool> keep = inbz; // keep those points that are inside the 1st Brillouin zone
    for (size_t i=0; i<hkl.size(); ++i)
        if (!inbz.getvalue(i) && inbz.extract(this->get_neighbours(i)).areanytrue())
          keep.insert(true,i); // and any out-of-zone points with at least one neighbour in-zone
    slong kept = 0;
    for (size_t i=0; i<hkl.size(); ++i)
      this->map[i] = keep.getvalue(i) ? kept++ : slong(-1);
  };
};

template<class T> class BrillouinZoneGrid4: public InterpolateGrid4<T>{
protected:
  BrillouinZone brillouinzone;
public:
  BrillouinZoneGrid4(const BrillouinZone bz, const double *spec, const double *d, const int isrlu=1): brillouinzone(bz) { this->determine_map_size(spec,d,isrlu);};
  BrillouinZoneGrid4(const BrillouinZone bz, const double *spec, const size_t *n): brillouinzone(bz) { this->determine_map_step(spec,n); };
  const BrillouinZone get_brillouinzone() const { return this->brillouinzone; };
  //
  ArrayVector<double> get_grid_hkl() const {
    ArrayVector<double> xyz = this->get_grid_xyz();
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!matrix_inverse(fromxyz,toxyz)) printf("transform matrix toxyz has zero determinant?!\n");
    ArrayVector<double> hkl(3,xyz.size());
    for (size_t i=0; i<xyz.size(); i++) multiply_matrix_vector<double,double,double,3>(hkl.datapointer(i), fromxyz, xyz.datapointer(i));
    return hkl;
  };
  ArrayVector<double> get_grid_hkle() const {
    ArrayVector<double> hkl = this->get_grid_hkl();
    ArrayVector<double> xyzw = this->get_grid_xyzw();
    for (size_t idx=0; idx<hkl.size(); ++idx)
      for (size_t i3=0; i3<this->size(3); ++i3)
        for (size_t j=0; j<3u; ++j)
          xyzw.insert( hkl.getvalue(idx,j), idx*hkl.size()+i3 , j );
    return xyzw; // which is now hklw (and w==e)
  };
  ArrayVector<double> get_mapped_hkl() const {
    ArrayVector<double> xyz = this->get_mapped_xyz();
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!matrix_inverse(fromxyz,toxyz)) printf("transform matrix toxyz has zero determinant?!\n");
    ArrayVector<double> hkl(3,xyz.size());
    for (size_t i=0; i<xyz.size(); i++) multiply_matrix_vector<double,double,double,3>(hkl.datapointer(i), fromxyz, xyz.datapointer(i));
    return hkl;
  };
  ArrayVector<double> get_mapped_hkle() const {
    ArrayVector<double> hkl = this->get_mapped_hkl();
    ArrayVector<double> xyzw = this->get_mapped_xyzw();
    // this *WILL NOT* give expected results if we ever set parts of the map invalid as a function of the fourth (energy) dimension!
    for (size_t idx=0; idx<hkl.size(); ++idx)
      for (size_t i3=0; i3<this->size(3); ++i3)
        for (size_t j=0; j<3u; ++j)
          xyzw.insert( hkl.getvalue(idx,j), idx*hkl.size()+i3 , j );
    return xyzw; // which is now hklw (and w==e)
  };
protected:
  // spec should be [Emin, Estep, Emax]
  void determine_map_size(const double *spec, const double *d_in, const int isrlu){
    double d[4];
    for (int i=0;i<3;i++) d[i]=2.0*d_in[i]; // ×2 here so we can double N later
    if (isrlu){
      LQVec<int> ei(this->brillouinzone.get_lattice(),1u);
      for (int i=0;i<3;i++){
        for (int j=0;j<3;j++) ei.insert(i==j?1:0,0,j);
        d[i] *= norm(ei).getvalue(0);
      }
    }
    // d is "guaranteed" to now be in units of inverse angstrom (and meV).
    ArrayVector<double> vxyz = this->brillouinzone.get_vertices().get_xyz();
    double minx[3], maxx[3];
    for (int i=0; i<3; i++) {minx[i]= std::numeric_limits<double>::max(); maxx[i]= -std::numeric_limits<double>::max();} // +/- max
    for (size_t i=0; i<vxyz.size(); i++){
      for (size_t j=0; j<3u; j++){
        if (vxyz.getvalue(i,j) < minx[j]) minx[j] = vxyz.getvalue(i,j);
        if (vxyz.getvalue(i,j) > maxx[j]) maxx[j] = vxyz.getvalue(i,j);
      }
    }
    size_t n[4];
    for (int i=0; i<3; i++) n[i] = 2*std::ceil( (maxx[i]-minx[i])/d[i] ); // ×2 to ensure we have an even number of steps
    n[3] = std::ceil( (spec[2]+spec[1]-spec[0])/spec[1] ); // 2+1-0 to ensure 2 is included
    // account for using ceil:
    for (int i=0; i<3; i++) d[i] = (maxx[i]-minx[i])/n[i];
    d[3] = (spec[2]+spec[1]-spec[0])/(n[3]);
    // for interpolation purposes, we want to make sure we go one-step beyond
    // the zone boundary in each Q direction:
    for (int i=0; i<3; i++) n[i] += 2;
    // which means we need to shift the origin of the grid as well by one step
    double z[4];
    for (int i=0; i<3; i++) z[i] = 0.0 - d[i]*(n[i]/2 -1); // ensure we hit zero, and that N/2 points are on the postive side
    z[3] = spec[0]; // and that we dont offset the energy zero

    this->resize(n);
    this->set_step(d);
    this->set_zero(z);
    this->set_map(); // fills the 3D map with linear indices -- eventually replace this with a fancier mapping that (a) respects the zone boundary and (b) sticks to the irreducible wedge
    this->truncate_grid_to_brillouin_zone();
  };

  void determine_map_step(const double *spec, const size_t *in){
    size_t n[4];
    bool iszero[3];
    for (int i=0;i<3;i++) {
      iszero[i] = 0u==in[i];
      n[i]= (iszero[i]) ? in[i]+1 : in[i];
    }
    n[3] = std::ceil( (spec[2]+spec[1]-spec[0])/spec[1] );
    ArrayVector<double> vxyz = this->brillouinzone.get_vertices().get_xyz();
    double maxx[3];
    for (int i=0; i<3; i++) maxx[i]= -std::numeric_limits<double>::max(); // - max
    for (size_t i=0; i<vxyz.size(); i++)
      for (size_t j=0; j<3u; j++)
        if (vxyz.getvalue(i,j) > maxx[j]) maxx[j] = vxyz.getvalue(i,j);
    double d[4];
    for (int i=0; i<3; i++) d[i] = (n[i]<2) ? (iszero[i] ? 0.0 : maxx[i]) : maxx[i]/(n[i]-1);
    d[3] = (spec[2]+spec[1]-spec[0])/n[3];
    double z[4];
    for (int i=0; i<3; i++) z[i] = 0.0 - d[i]*(n[i]-1); // ensure we hit zero, and that N/2 points are on the postive side
    z[3] = spec[0];

    for (int i=0;i<3;i++) if (!iszero[i]) n[i]*=2; // we actually want 2N unless the input requested zero points, in which case we *need* one

    this->resize(n);
    this->set_step(d);
    this->set_zero(z);
    this->set_map();
    this->truncate_grid_to_brillouin_zone();
  };

  void truncate_grid_to_brillouin_zone(){
    LQVec<double> hkl(this->brillouinzone.get_lattice(),this->get_grid_hkl());
    ArrayVector<bool> inbz = this->brillouinzone.isinside(&hkl);
    ArrayVector<bool> keep = inbz; // keep those points that are inside the 1st Brillouin zone
    for (size_t i=0; i<hkl.size(); ++i)
        if (!inbz.getvalue(i) && inbz.extract(this->get_neighbours(i)).areanytrue())
          keep.insert(true,i); // and any out-of-zone points with at least one neighbour in-zone
    slong kept = 0;
    for (size_t i=0; i<hkl.size(); ++i)
      this->map[i] = keep.getvalue(i) ? kept++ : slong(-1);
  };
};


#endif
