typedef long slong;

#ifndef _BZ_MESH_
#define _BZ_MESH_

#include "bz.h"
#include "mesh.h"

template<class T> class BrillouinZoneMesh3: public Mesh3<T>{
protected:
  BrillouinZone brillouinzone;
public:
  /*! Construct using a maximum tetrahedron volume -- makes a tetrahedron mesh
      instead of a orthogonal grid.
      @param bz The BrillouinZone object
      @param vol The maximum tetrahedron volume
      @param isrlu A flag to indicate if vol is in units of rlu (isrlu=1) or inverse angstrom (isrlu=0)
      @note If vol is in relative lattice units an absolute volume will be
            determined using the unit cell volume of the underlying lattice.
  */
  // BrillouinZoneMesh3(const BrillouinZone& bz, const double max_size_invA=-1., const double min_angle=20.0, const double max_angle=-1.0, const double max_ratio=-1., const int max_points=-1):
  //   Mesh3<T>(bz.get_ir_vertices().get_xyz(), bz.get_ir_vertices_per_face(), max_size_invA, min_angle, max_angle, max_ratio, max_points),
  //   brillouinzone(bz) {}
  // BrillouinZoneMesh3(const BrillouinZone& bz) Mesh3<T>(bz.get_ir_vertices().get_xyz(), bz.get_ir_vertices_per_face());
  template<typename... A>
  BrillouinZoneMesh3(const BrillouinZone& bz, A... args):
    Mesh3<T>(bz.get_ir_vertices().get_xyz(), bz.get_ir_vertices_per_face(), args...),
    brillouinzone(bz) {}
  // get the BrillouinZone object
  BrillouinZone get_brillouinzone(void) const {return this->brillouinzone;}
  // get the mesh vertices in relative lattice units
  ArrayVector<double> get_mesh_hkl(void) const {
    ArrayVector<double> xyz = this->get_mesh_xyz();
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!matrix_inverse(fromxyz,toxyz)) throw std::runtime_error("transform matrix toxyz has zero determinant");
    ArrayVector<double> hkl(3,xyz.size());
    for (size_t i=0; i<xyz.size(); i++) multiply_matrix_vector<double,double,double,3>(hkl.data(i), fromxyz, xyz.data(i));
    return hkl;
  }
  template<typename R> ArrayVector<T> interpolate_at(const LQVec<R>& x, const int nthreads, const bool no_move=false) const{
    LQVec<R> ir_q(x.get_lattice(), x.size());
    LQVec<int> tau(x.get_lattice(), x.size());
    std::vector<std::array<int,9>> rots(x.size());
    BrillouinZone bz = this->get_brillouinzone();

    std::string msg;
    if (no_move){
      // Special mode for testing where no specified points are moved
      // IT IS IMPERITIVE THAT THE PROVIDED POINTS ARE *INSIDE* THE IRREDUCIBLE
      // POLYHEDRON otherwise the interpolation will fail or give garbage back.
      ir_q = x;
      for (size_t i=0; i<x.size(); ++i) rots[i] = {1,0,0, 0,1,0, 0,0,1};
    } else if (!bz.ir_moveinto(x, ir_q, tau, rots)){
      msg = "Moving all points into the irreducible Brillouin zone failed.";
      throw std::runtime_error(msg);
    }
    ArrayVector<T> ir_result = (nthreads < 2)
      ? this->Mesh3<T>::interpolate_at(ir_q.get_xyz())
      : this->Mesh3<T>::parallel_interpolate_at(ir_q.get_xyz(), nthreads);

    // any eigenvector, vector, and matrix (treated as rank-2 tensor) output of
    // the interpolation needs to be rotated.
    if (this->elements[1] || this->elements[2] || this->elements[3]){
      if (this->elements[1] % 3){
        msg = "Eigenvectors should consist of 3 elements (per ion) for each branch: ";
        msg += std::to_string(this->elements[1]) + "%3 != 0";
        throw std::runtime_error(msg);
      }
      if (this->elements[2] %3){
        msg = "Vectors should consist of 3N elements for each branch: ";
        msg += std::to_string(this->elements[2]) + "%3 != 0";
        throw std::runtime_error(msg);
      }
      if (this->elements[3] != 0u && this->elements[3] != 3u){
        msg = "Matrices should be 3x3 for each branch:";
        std::string m = std::to_string(this->elements[3]);
        msg += m + "x" + m + " != 3x3";
        throw std::runtime_error(msg);
      }
      size_t ne = this->elements[1]/3u;
      size_t nv = this->elements[2]/3u;
      size_t nm = this->elements[3]/3u;
      size_t sp = this->elements[0] + ne*3u + nv*3u + nm*9u;
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
        for (size_t b=0; b<this->branches; ++b){
          // we can skip the scalar elements, as they do not rotate.
          offset = b*sp + this->elements[0];
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

#endif // _BZ_MESH_
