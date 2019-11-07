typedef long slong;

#ifndef _BZ_NEST_
#define _BZ_NEST_

#include "bz.h"
#include "nest.h"

template<class T> class BrillouinZoneNest3: public Nest<T>{
  BrillouinZone brillouinzone;
public:
  template<typename... A>
  BrillouinZoneNest3(const BrillouinZone& bz, A... args):
    Nest<T>(bz.get_ir_polyhedron(), args...),
    brillouinzone(bz) {}
  //! get the BrillouinZone object
  BrillouinZone get_brillouinzone(void) const {return this->brillouinzone;}
  //! get the vertices of the leaf vertices in inverse Angstrom
  ArrayVector<double> get_xyz(void) const {return this->vertices();}
  //! get the vertices of all vertices in absolute units
  const ArrayVector<double>& get_all_xyz(void) const {return this->all_vertices(); }
  //! get the vertices of the leaf vertices in relative lattice units
  ArrayVector<double> get_hkl(void) const { return xyz_to_hkl(brillouinzone.get_lattice(),this->vertices());}
  //! get the vertices of the inner (cubic) nodes in relative lattice units
  ArrayVector<double> get_all_hkl(void) const {return xyz_to_hkl(brillouinzone.get_lattice(),this->all_vertices()); }
  // //! get the indices forming the faces of the tetrahedra
  // std::vector<std::array<size_t,4>> get_vertices_per_tetrahedron(void) const {return this->tetrahedra();}


  template<typename R> ArrayVector<T> interpolate_at(const LQVec<R>& x, const int nthreads, const bool no_move=false) const{
    LQVec<R> ir_q(x.get_lattice(), x.size());
    LQVec<int> tau(x.get_lattice(), x.size());
    std::vector<std::array<int,9>> rots(x.size());
    const std::array<unsigned,4>& el{this->data().elements()};
    std::string msg;
    if (no_move){
      // Special mode for testing where no specified points are moved
      // IT IS IMPERITIVE THAT THE PROVIDED POINTS ARE *INSIDE* THE IRREDUCIBLE
      // POLYHEDRON otherwise the interpolation will fail or give garbage back.
      ir_q = x;
      for (size_t i=0; i<x.size(); ++i) rots[i] = {1,0,0, 0,1,0, 0,0,1};
    } else if (!brillouinzone.ir_moveinto(x, ir_q, tau, rots)){
      msg = "Moving all points into the irreducible Brillouin zone failed.";
      throw std::runtime_error(msg);
    }
    ArrayVector<T> ir_result = nthreads > 1 ? this->Nest<T>::interpolate_at(ir_q.get_xyz(), nthreads) : this->Nest<T>::interpolate_at(ir_q.get_xyz());
    // ArrayVector<T> ir_result = this->Nest<T>::interpolate_at(ir_q.get_xyz(), nthreads);

    // any eigenvector, vector, and matrix (treated as rank-2 tensor) output of
    // the interpolation needs to be rotated.
    if (el[1] || el[2] || el[3]){
      if (el[1] % 3){
        msg = "Eigenvectors should consist of 3 elements (per ion) for each branch: ";
        msg += std::to_string(el[1]) + "%3 != 0";
        throw std::runtime_error(msg);
      }
      if (el[2] %3){
        msg = "Vectors should consist of 3N elements for each branch: ";
        msg += std::to_string(el[2]) + "%3 != 0";
        throw std::runtime_error(msg);
      }
      if (el[3] != 0u && el[3] != 3u){
        msg = "Matrices should be 3x3 for each branch:";
        std::string m = std::to_string(el[3]);
        msg += m + "x" + m + " != 3x3";
        throw std::runtime_error(msg);
      }
      size_t ne = el[1]/3u;
      size_t nv = el[2]/3u;
      size_t nm = el[3]/3u;
      size_t sp = el[0] + ne*3u + nv*3u + nm*9u;
      T tmp_v[3];
      T tmp_m[9];
      std::vector<std::array<int,9>> invR;
      if (nm){ // only allocate and calculate invR if we need it
        invR.resize(rots.size());
        for (size_t i=0; i<rots.size(); ++i)
          matrix_inverse(invR[i].data(), rots[i].data());
      }

      size_t offset;
      for (size_t i=0; i<ir_result.size(); ++i){
        for (size_t b=0; b<this->data().branches(); ++b){
          // we can skip the scalar elements, as they do not rotate.
          offset = b*sp + el[0];
          // eigenvectors and regular vectors rotate the same way
          for (size_t v=0; v<(ne+nv); ++v){
            mul_mat_vec(tmp_v, 3u, rots[i].data(), ir_result.data(i, offset+v*3u));
            for (size_t j=0; j<3u; ++j) ir_result.insert(tmp_v[j], i, offset+v*3u+j);
          }
          offset += (ne+nv)*3u;
          for (size_t m=0; m<nm; ++m){
            // we want R*M*R⁻¹.
            // first calculate M*R⁻¹, storing in tmp_m
            mul_mat_mat(tmp_m, 3u, ir_result.data(i, offset+m*9u), invR[i].data());
            // next calculate R*tmp_m, storing back in the ir_result array
            mul_mat_mat(ir_result.data(i, offset+m*9u), 3u, rots[i].data(), tmp_m);
          }
        }
      }
    }
    return ir_result;
  }
};

#endif
