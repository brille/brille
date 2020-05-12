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

typedef long slong; // ssize_t is only defined for gcc?

#ifndef _BZ_GRID_
#define _BZ_GRID_

#include "bz.hpp"
#include "grid.hpp"
#include "grid4.hpp"

// #include "triangulation.hpp"

// TODO: Be clever about the map/grid creation (like spglib)?

/*! \brief A class to hold both a BrillouinZone and InterpolateGrid3

  By adding a BrillouinZone object to an InterpolationGrid3 object, it is
  possible to interpolate at arbitrary points within the first Brillouin zone.
*/
template<class T, class R> class BrillouinZoneGrid3: public InterpolateGrid3<T,R>{
protected:
  BrillouinZone brillouinzone;
public:
  /*! Construct using step sizes
    @param bz The BrillouinZone object
    @param d The three step sizes in Q
    @param isrlu A flag to indicate if d is in units of rlu (default) or inverse angstrom
    @note  If d is in relative lattice units an absolute length for each is
           determined by scaling the lengths of the underlying lattice basis
           vectors; e.g., d[0]*= |(100)|, d[1]*=|(010)|, d[2]*=|(001)|
  */
  BrillouinZoneGrid3(const BrillouinZone& bz, const double *d, const int isrlu=1): brillouinzone(bz) { this->determine_map_size(d,isrlu);}
  /*! Construct using number of steps
    @param bz The BrillouinZone object
    @param n The three number of steps
    @note The number of steps along each dimension of the InterpolationGrid3 is
          a function of n. For a given dimension, i, N[i] = 2*n[i]+1 if n[i]>0
          or N[i] = 1 if n[i]==0
  */
  BrillouinZoneGrid3(const BrillouinZone& bz, const size_t *n): brillouinzone(bz) { this->determine_map_step(n); }
  /*! Construct using a maximum tetrahedron volume -- makes a tetrahedron mesh
      instead of a orthogonal grid.
      @param bz The BrillouinZone object
      @param vol The maximum tetrahedron volume
      @param isrlu A flag to indicate if vol is in units of rlu (isrlu=1) or inverse angstrom (isrlu=0)
      @note If vol is in relative lattice units an absolute volume will be
            determined using the unit cell volume of the underlying lattice.
  */
  // BrillouinZoneGrid3(const BrillouinZone bz, const double vol, const int isrlu=1): brillouinzone(bz) {this->determine_map_tri(vol, isrlu);}
  //! Return the BrillouinZone object contained
  const BrillouinZone get_brillouinzone() const { return this->brillouinzone; }
  /*! Return the grid points of the InterpolateGrid3 object in relative lattice units
    @param maxN The maximum number of 3-vectors for which space has been allocated
    @param[out] hkl Where the vectors will be stored
    @returns the number of 3-vectors placed into hkl
  */
  size_t get_grid_hkl(const size_t maxN, double *hkl) const {
    ArrayVector<double> xyz = this->get_grid_xyz();
    if (xyz.size() > maxN) throw std::out_of_range("Output does not have enough storage space.");
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!matrix_inverse(fromxyz,toxyz)) throw std::runtime_error("transform matrix toxyz has zero determinant");
    for (size_t i=0; i<xyz.size(); i++) multiply_matrix_vector<double,double,double,3>(hkl+3*i, fromxyz, xyz.data(i));
    return xyz.size();
  }
  /*! Return the grid points of the InterpolateGrid3 object in relative lattice units
    @returns An ArrayVector containing the 3-vectors
  */
  ArrayVector<double> get_grid_hkl() const {
    ArrayVector<double> xyz = this->get_grid_xyz();
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!matrix_inverse(fromxyz,toxyz)) throw std::runtime_error("transform matrix toxyz has zero determinant");
    ArrayVector<double> hkl(3,xyz.size());
    for (size_t i=0; i<xyz.size(); i++) multiply_matrix_vector<double,double,double,3>(hkl.data(i), fromxyz, xyz.data(i));
    return hkl;
  }
  /*! Return only the mapped grid points of the InterpolateGrid3 object in relative lattice units
    @returns An ArrayVector containing the mapped 3-vectors
  */
  ArrayVector<double> get_mapped_hkl() const {
    ArrayVector<double> xyz = this->get_mapped_xyz();
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!matrix_inverse(fromxyz,toxyz)) throw std::runtime_error("transform matrix toxyz has zero determinant");
    ArrayVector<double> hkl(3,xyz.size());
    for (size_t i=0; i<xyz.size(); i++) multiply_matrix_vector<double,double,double,3>(hkl.data(i), fromxyz, xyz.data(i));
    return hkl;
  }
  template<typename S>
  std::tuple<ArrayVector<T>,ArrayVector<R>>
  ir_interpolate_at(const LQVec<S>& x, const int nthreads, const bool no_move=false) const{
    LQVec<S> ir_q(x.get_lattice(), x.size());
    LQVec<int> tau(x.get_lattice(), x.size());
    std::vector<size_t> rot(x.size(),0u), invrot(x.size(),0u);
    if (no_move){
      ir_q = x;
    } else if (!brillouinzone.ir_moveinto(x, ir_q, tau, rot, invrot, nthreads)){
      std::string msg;
      msg = "Moving all points into the irreducible Brillouin zone failed.";
      throw std::runtime_error(msg);
    }
    // perform the interpolation within the irreducible Brillouin zone
    ArrayVector<T> vals;
    ArrayVector<R> vecs;
    std::tie(vals,vecs) =
      (nthreads > 1) ? this->InterpolateGrid3<T,R>::parallel_linear_interpolate_at(ir_q.get_xyz(), nthreads)
                     : this->InterpolateGrid3<T,R>::linear_interpolate_at(ir_q.get_xyz());
    // we always need the pointgroup operations to 'rotate'
    PointSymmetry psym = brillouinzone.get_pointgroup_symmetry();
    // and might need the Phonon Gamma table
    GammaTable pgt{GammaTable()};
    if (RotatesLike::Gamma == this->data().vectors().rotateslike()){
      pgt.construct(brillouinzone.get_lattice().star(), brillouinzone.add_time_reversal());
    }
    // actually perform the rotation to Q
    this->data().values() .rotate_in_place(vals, ir_q, pgt, psym, rot, invrot, nthreads);
    this->data().vectors().rotate_in_place(vecs, ir_q, pgt, psym, rot, invrot, nthreads);
    // we're done so bundle the output
    return std::make_tuple(vals, vecs);
  }
protected:
  /*! Determines and sets the properties of the underlying grid from three step sizes
    @param d_in A pointer to the three step sizes
    @param isrlu A flag indicating if the step sizes are in relative lattice units or inverse Angstrom
  */
  void determine_map_size(const double *d_in, const int isrlu){
    double d[3];
    for (int i=0;i<3;i++) d[i]=2.0*d_in[i]; // ×2 here so we can double N later
    if (isrlu){
      LQVec<int> ei(this->brillouinzone.get_lattice(),1u);
      for (int i=0;i<3;i++){
        for (int j=0;j<3;j++) ei.insert(i==j?1:0,0,j);
        d[i] *= norm(ei).getvalue(0);
      }
    }
    // d is "guaranteed" to now be in units of inverse angstrom.
    ArrayVector<double> vxyz = this->brillouinzone.get_vertices().get_xyz();
    double minx[3], maxx[3];
    for (int i=0; i<3; i++) {minx[i]= (std::numeric_limits<double>::max)(); maxx[i]= std::numeric_limits<double>::lowest();} // +/- max
    for (size_t i=0; i<vxyz.size(); i++){
      for (size_t j=0; j<3u; j++){
        if (vxyz.getvalue(i,j) < minx[j]) minx[j] = vxyz.getvalue(i,j);
        if (vxyz.getvalue(i,j) > maxx[j]) maxx[j] = vxyz.getvalue(i,j);
      }
    }
    size_t n[3];
    for (int i=0; i<3; i++)
      n[i] = 2*((size_t)std::ceil( (maxx[i]-minx[i])/d[i] ))+1; // ×2 to ensure we have an odd number of steps
    // account for using ceil:
    for (int i=0; i<3; i++) d[i] = (maxx[i]-minx[i])/n[i];
    // for interpolation purposes, we want to make sure we go one-step beyond
    // the zone boundary in each direction:
    for (int i=0; i<3; i++) n[i] += 2;
    double z[3];
    for (int i=0; i<3; i++) z[i] = 0.0 - d[i]*(n[i]-1)/2; // ensure we hit zero, and that N/2 points are on the postive side

    this->resize(n);
    this->set_step(d);
    this->set_zero(z);
    this->set_map(); // fills the 3D map with linear indices -- eventually replace this with a fancier mapping that (a) respects the zone boundary and (b) sticks to the irreducible wedge
    this->truncate_grid_to_brillouin_zone();
  }
  /*! Determines and sets the properties of the underying grid from three step counts
    @param in The input number of steps in three directions
    @note The number of grid points will not be prod(in) but rather prod(N) where
          N[i] = 2*in[i]+1 for in[i]>0 or N[i] = 1 for in[i]==0
  */
  void determine_map_step(const size_t *in){
    size_t n[3];
    bool iszero[3];
    for (int i=0;i<3;i++) {
      iszero[i] = 0u==in[i];
      n[i]= (iszero[i]) ? in[i]+1 : in[i];
    }
    ArrayVector<double> vxyz = this->brillouinzone.get_vertices().get_xyz();
    double maxx[3];
    for (int i=0; i<3; i++) maxx[i]= std::numeric_limits<double>::lowest(); // - max
    for (size_t i=0; i<vxyz.size(); i++)
      for (size_t j=0; j<3u; j++)
        if (vxyz.getvalue(i,j) > maxx[j]) maxx[j] = vxyz.getvalue(i,j);
    double d[3];
    for (int i=0; i<3; i++) d[i] = (n[i]<2) ? (iszero[i] ? 0.0 : maxx[i]) : maxx[i]/n[i];
    double z[3];
    for (int i=0; i<3; i++) z[i] = 0.0 - d[i]*n[i]; // ensure we hit zero, and that N/2 points are on the postive side

    for (int i=0;i<3;i++) if (!iszero[i]) n[i]=2*n[i]+1; // we actually want 2N+1 unless the input requested zero points, in which case we *need* one

    this->resize(n);
    this->set_step(d);
    this->set_zero(z);
    this->set_map();
    this->truncate_grid_to_brillouin_zone();
  }

  /*! \brief Remove unusable grid points

    Utilizing methods of the BrillouinZone object, determine which grid points
    are inside of the first Brillouin zone or within one-step of such a point
    and assign their map indexes sequentially.
    For all points not inside the first Brillouin zone or within one-step of
    such a point, set the map index to an invalid value.
  */
  void truncate_grid_to_brillouin_zone(){
    LQVec<double> hkl(this->brillouinzone.get_lattice(),this->get_grid_hkl());
    ArrayVector<bool> inbz = this->brillouinzone.isinside(hkl);
    ArrayVector<bool> inir = this->brillouinzone.isinside_wedge(hkl);
    ArrayVector<bool> keep(1u, hkl.size());
    // keep those points that are inside the 1st Brillouin zone and the irreducible wedge
    for (size_t i=0; i<hkl.size(); ++i) keep.insert(inbz.getvalue(i)&&inir.getvalue(i), i);
    ArrayVector<bool> orig(keep); // an unmodified copy so that we don't grow beyond 1-neighbour
    for (size_t i=0; i<hkl.size(); ++i)
        if (!orig.getvalue(i) && orig.extract(this->get_neighbours(i)).any_true())
          keep.insert(true,i); // and any out-of-zone points with at least one neighbour in-zone
    slong kept = 0;
    for (size_t i=0; i<hkl.size(); ++i)
      this->map[i] = keep.getvalue(i) ? kept++ : slong(-1);
  }
};

/*! \brief A class to hold both a BrillouinZone and InterpolateGrid4

  By adding a BrillouinZone object to an InterpolationGrid4 object, it is
  possible to interpolate at arbitrary points within the first Brillouin zone.
*/
template<class T,class S> class BrillouinZoneGrid4: public InterpolateGrid4<T,S>{
protected:
  BrillouinZone brillouinzone;
public:
  /*! Construct using an energy-axis specification and step sizes
    @param spec The energy specification, expected of the form [Emin,Estep,Emax]
    @param bz The BrillouinZone object
    @param d The three step sizes in Q
    @param isrlu A flag to indicate if d is in units of rlu (default) or inverse angstrom
    @note  If d is in relative lattice units an absolute length for each is
           determined by scaling the lengths of the underlying lattice basis
           vectors; e.g., d[0]*= |(100)|, d[1]*=|(010)|, d[2]*=|(001)|
  */
  BrillouinZoneGrid4(const BrillouinZone& bz, const double *spec, const double *d, const int isrlu=1): brillouinzone(bz) { this->determine_map_size(spec,d,isrlu);};
  /*! Construct using an energy-axis specification and number of steps
    @param spec The energy specification, expected of the form [Emin,Estep,Emax]
    @param bz The BrillouinZone object
    @param n The three number of steps
    @note The number of steps along each dimension of the InterpolationGrid3 is
          a function of n. For a given dimension, i, N[i] = 2*n[i]+1 if n[i]>0
          or N[i] = 1 if n[i]==0
  */
  BrillouinZoneGrid4(const BrillouinZone& bz, const double *spec, const size_t *n): brillouinzone(bz) { this->determine_map_step(spec,n); };
  // /*! Construct using a maximum tetrahedron volume -- makes a tetrahedron mesh
  //     instead of a orthogonal grid.
  //     @param spec The energy specification, expected of the form [Emin,Estep,Emax]
  //     @param bz The BrillouinZone object
  //     @param vol The maximum tetrahedron volume
  //     @param isrlu A flag to indicate if vol is in units of rlu (isrlu=1) or inverse angstrom (isrlu=0)
  //     @note If vol is in relative lattice units an absolute volume will be
  //           determined using the unit cell volume of the underlying lattice.
  // */
  // BrillouinZoneGrid4(const BrillouinZone bz, const double *spec, const double vol, const int isrlu=1): brillouinzone(bz) {};
  //! Return the BrillouinZone object contained
  const BrillouinZone get_brillouinzone() const { return this->brillouinzone; };
  /*! Return the Q part of the grid points of the InterpolateGrid4 object in relative lattice units
    @returns An ArrayVector containing the 3-vectors
  */
  ArrayVector<double> get_grid_hkl() const {
    ArrayVector<double> xyz = this->get_grid_xyz();
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!matrix_inverse(fromxyz,toxyz)) throw std::runtime_error("transform matrix toxyz has zero determinant");
    ArrayVector<double> hkl(3,xyz.size());
    for (size_t i=0; i<xyz.size(); i++) multiply_matrix_vector<double,double,double,3>(hkl.data(i), fromxyz, xyz.data(i));
    return hkl;
  };
  /*! Return the full grid points of the InterpolateGrid4 object in relative lattice units
    @returns An ArrayVector containing the (Q,E) 4-vectors
  */
  ArrayVector<double> get_grid_hkle() const {
    ArrayVector<double> hkl = this->get_grid_hkl();
    ArrayVector<double> xyzw = this->get_grid_xyzw();
    for (size_t idx=0; idx<hkl.size(); ++idx)
      for (size_t i3=0; i3<this->size(3); ++i3)
        for (size_t j=0; j<3u; ++j)
          xyzw.insert( hkl.getvalue(idx,j), idx*hkl.size()+i3 , j );
    return xyzw; // which is now hklw (and w==e)
  };
  /*! Return the Q part of the mapped grid points of the InterpolateGrid4 object in relative lattice units
    @returns An ArrayVector containing the 3-vectors
  */
  ArrayVector<double> get_mapped_hkl() const {
    ArrayVector<double> xyz = this->get_mapped_xyz();
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!matrix_inverse(fromxyz,toxyz)) throw std::runtime_error("transform matrix toxyz has zero determinant");
    ArrayVector<double> hkl(3,xyz.size());
    for (size_t i=0; i<xyz.size(); i++) multiply_matrix_vector<double,double,double,3>(hkl.data(i), fromxyz, xyz.data(i));
    return hkl;
  };
  /*! Return the full mapped grid points of the InterpolateGrid4 object in relative lattice units
    @returns An ArrayVector containing the (Q,E) 4-vectors
  */
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
  /*! Determines and sets the properties of the underlying grid from an energy specification and three step sizes
    @param spec The energy specification, of the form [Emin,Estep,Emax]
    @param d_in A pointer to the three step sizes
    @param isrlu A flag indicating if the step sizes are in relative lattice units or inverse Angstrom
  */
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
    for (int i=0; i<3; i++) {minx[i]= (std::numeric_limits<double>::max)(); maxx[i]= std::numeric_limits<double>::lowest();} // +/- max
    for (size_t i=0; i<vxyz.size(); i++){
      for (size_t j=0; j<3u; j++){
        if (vxyz.getvalue(i,j) < minx[j]) minx[j] = vxyz.getvalue(i,j);
        if (vxyz.getvalue(i,j) > maxx[j]) maxx[j] = vxyz.getvalue(i,j);
      }
    }
    size_t n[4];
    for (int i=0; i<3; i++)
      n[i] = 2*((size_t)std::ceil( (maxx[i]-minx[i])/d[i] ))+1; // ×2 to ensure we have an even number of steps
    n[3] = (size_t)std::ceil( (spec[2]+spec[1]-spec[0])/spec[1] ); // 2+1-0 to ensure 2 is included
    // account for using ceil:
    for (int i=0; i<3; i++) d[i] = (maxx[i]-minx[i])/n[i];
    d[3] = (spec[2]+spec[1]-spec[0])/(n[3]);
    // for interpolation purposes, we want to make sure we go one-step beyond
    // the zone boundary in each Q direction:
    for (int i=0; i<3; i++) n[i] += 2;
    // which means we need to shift the origin of the grid as well by one step
    double z[4];
    for (int i=0; i<3; i++) z[i] = 0.0 - d[i]*(n[i]-1)/2; // ensure we hit zero, and that N/2 points are on the postive side
    z[3] = spec[0]; // and that we dont offset the energy zero

    this->resize(n);
    this->set_step(d);
    this->set_zero(z);
    this->set_map(); // fills the 3D map with linear indices -- eventually replace this with a fancier mapping that (a) respects the zone boundary and (b) sticks to the irreducible wedge
    this->truncate_grid_to_brillouin_zone();
  };
  /*! Determines and sets the properties of the underying grid from an energy specification and three step counts
    @param spec The energy specification, of the form [Emin,Estep,Emax]
    @param in The input number of steps in three directions
    @note The number of grid points will not be prod(in) but rather prod(N) where
          N[i] = 2*in[i]+1 for in[i]>0 or N[i] = 1 for in[i]==0
  */
  void determine_map_step(const double *spec, const size_t *in){
    size_t n[4];
    bool iszero[3];
    for (int i=0;i<3;i++) {
      iszero[i] = 0u==in[i];
      n[i]= (iszero[i]) ? in[i]+1 : in[i];
    }
    n[3] = (size_t)std::ceil( (spec[2]+spec[1]-spec[0])/spec[1] );
    ArrayVector<double> vxyz = this->brillouinzone.get_vertices().get_xyz();
    double maxx[3];
    for (int i=0; i<3; i++) maxx[i]= std::numeric_limits<double>::lowest(); // - max
    for (size_t i=0; i<vxyz.size(); i++)
      for (size_t j=0; j<3u; j++)
        if (vxyz.getvalue(i,j) > maxx[j]) maxx[j] = vxyz.getvalue(i,j);
    double d[4];
    for (int i=0; i<3; i++) d[i] = (n[i]<2) ? (iszero[i] ? 0.0 : maxx[i]) : maxx[i]/n[i];
    d[3] = (spec[2]+spec[1]-spec[0])/n[3];
    double z[4];
    for (int i=0; i<3; i++) z[i] = 0.0 - d[i]*n[i]; // ensure we hit zero, and that N/2 points are on the postive side
    z[3] = spec[0];

    for (int i=0;i<3;i++) if (!iszero[i]) n[i]=2*n[i]+1; // we actually want 2N+1 unless the input requested zero points, in which case we *need* one

    this->resize(n);
    this->set_step(d);
    this->set_zero(z);
    this->set_map();
    this->truncate_grid_to_brillouin_zone();
  };
  /*! \brief Remove unusable grid points

    Utilizing methods of the BrillouinZone object, determine which grid points
    are inside of the first Brillouin zone or within one-step of such a point
    and assign their map indexes sequentially.
    For all points not inside the first Brillouin zone or within one-step of
    such a point, set the map index to an invalid value.
  */
  void truncate_grid_to_brillouin_zone(){
    LQVec<double> hkl(this->brillouinzone.get_lattice(),this->get_grid_hkl());
    ArrayVector<bool> inbz = this->brillouinzone.isinside(hkl);
    ArrayVector<bool> keep = inbz; // keep those points that are inside the 1st Brillouin zone
    for (size_t i=0; i<hkl.size(); ++i)
        if (!inbz.getvalue(i) && inbz.extract(this->get_neighbours(i)).any_true())
          keep.insert(true,i); // and any out-of-zone points with at least one neighbour in-zone
    slong kept = 0;
    for (size_t i=0; i<hkl.size(); ++i)
      this->map[i] = keep.getvalue(i) ? kept++ : slong(-1);
  };
};


#endif
