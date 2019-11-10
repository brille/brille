typedef long slong;

#ifndef _BZ_TRELLIS_
#define _BZ_TRELLIS_

#include "bz.h"
#include "trellis.h"

template<class T> class BrillouinZoneTrellis3: public PolyhedronTrellis<T>{
  BrillouinZone brillouinzone;
public:
  template<typename... A>
  BrillouinZoneTrellis3(const BrillouinZone& bz, A... args):
    PolyhedronTrellis<T>(bz.get_ir_polyhedron(), args...),
    brillouinzone(bz) {}
  //! get the BrillouinZone object
  BrillouinZone get_brillouinzone(void) const {return this->brillouinzone;}
  //! get the vertices of the trellis in absolute units
  const ArrayVector<double>& get_xyz(void) const {return this->vertices();}
  //! get the vertices of the inner (cubic) nodes in absolute units
  ArrayVector<double> get_inner_xyz(void) const {return this->cube_vertices(); }
  //! get the vertices of the outer (polyhedron) nodes in absolute units
  ArrayVector<double> get_outer_xyz(void) const {return this->poly_vertices(); }
  //! get the vertices of the trellis in relative lattice units
  ArrayVector<double> get_hkl(void) const { return xyz_to_hkl(brillouinzone.get_lattice(),this->vertices());}
  //! get the vertices of the inner (cubic) nodes in relative lattice units
  ArrayVector<double> get_inner_hkl(void) const {return xyz_to_hkl(brillouinzone.get_lattice(),this->cube_vertices()); }
  //! get the vertices of the outer (polyhedron) nodes in relative lattice units
  ArrayVector<double> get_outer_hkl(void) const {return xyz_to_hkl(brillouinzone.get_lattice(),this->poly_vertices()); }
  //! get the indices forming the faces of the tetrahedra
  std::vector<std::array<index_t,4>> get_vertices_per_tetrahedron(void) const {return this->vertices_per_tetrahedron();}


  template<typename R> ArrayVector<T> interpolate_at(const LQVec<R>& x, const int nthreads, const bool no_move=false) const{
    LQVec<R> ir_q(x.get_lattice(), x.size());
    LQVec<int> tau(x.get_lattice(), x.size());
    std::vector<std::array<int,9>> rots(x.size());
    if (no_move){
      // Special mode for testing where no specified points are moved
      // IT IS IMPERITIVE THAT THE PROVIDED POINTS ARE *INSIDE* THE IRREDUCIBLE
      // POLYHEDRON otherwise the interpolation will fail or give garbage back.
      ir_q = x;
      for (size_t i=0; i<x.size(); ++i) rots[i] = {1,0,0, 0,1,0, 0,0,1};
    } else if (!brillouinzone.ir_moveinto(x, ir_q, tau, rots)){
      std::string msg;
      msg = "Moving all points into the irreducible Brillouin zone failed.";
      throw std::runtime_error(msg);
    }
    ArrayVector<T> ir_result = (nthreads < 2)
      ? this->PolyhedronTrellis<T>::interpolate_at(ir_q.get_xyz())
      : this->PolyhedronTrellis<T>::interpolate_at(ir_q.get_xyz(), nthreads);

    if (nthreads < 2)
      this->data().rotate_in_place(ir_result, rots);
    else
      this->data().rotate_in_place(ir_result, rots, nthreads);

    return ir_result;
  }
};

#endif // _BZ_MESH_
