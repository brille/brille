#ifndef BRILLE_LATTICE_DUAL_HPP_
#define BRILLE_LATTICE_DUAL_HPP_

#include <assert.h>
#include <utility>
// #include <vector>
#include "enums.hpp"
#include "primitive.hpp"
#include "basis.hpp"
#include "array2.hpp"
#include "hdf_interface.hpp"
#include "math.hpp"
#include "hall_symbol.hpp"

namespace brille::lattice {

template<class T>
class Lattice{
public:
  using matrix_t = std::array<T,9>;
  using vector_t = std::array<T,3>;
private:
  matrix_t _vectors; //! Real space *column* basis vectors
  matrix_t _reciprocal;
  matrix_t _metric; //! The metric tensor G
  Bravais _bravais;
  Symmetry _space;
  PointSymmetry _point;
  Basis _basis;
public:
  Lattice(matrix_t v, matrix_t r, matrix_t m, Bravais b, Symmetry s, PointSymmetry p, Basis a)
      : _vectors{std::move(v)}, _reciprocal{std::move(r)}, _metric{std::move(m)}, _bravais(b), _space(s), _point(p), _basis(a)
  {}
//
  Lattice(const vector_t & lengths, const vector_t & angles, const Symmetry s, const LengthUnit lu, const AngleUnit au=AngleUnit::not_provided)
  {
    set_vectors(lengths, angles, lu, au);
    set_metric();
    spacegroup_symmetry(s);
    set_point_symmetry();
  }
  Lattice(const matrix_t & vectors, const Symmetry s, const LengthUnit lu)
  {
    set_vectors(vectors, lu);
    set_metric();
    spacegroup_symmetry(s);
    set_point_symmetry();
  }
  Lattice(const vector_t & lengths, const vector_t & angles, const std::string& s, const std::string& c, const LengthUnit lu, const AngleUnit au=AngleUnit::not_provided)
  {
    set_vectors(lengths, angles, lu, au);
    set_metric();
    set_space_symmetry(s, c);
    set_point_symmetry();
  }
  Lattice(const vector_t & lengths, const vector_t & angles, const std::string& s, const LengthUnit lu, const AngleUnit au=AngleUnit::not_provided)
  {
    set_vectors(lengths, angles, lu, au);
    set_metric();
    set_space_symmetry(s);
    set_point_symmetry();
  }
  Lattice(const matrix_t & vectors, const std::string& s, const std::string& c, const LengthUnit lu)
  {
    set_vectors(vectors, lu);
    set_metric();
    set_space_symmetry(s, c);
    set_point_symmetry();
  }
  Lattice(const matrix_t & vectors, const std::string& s, const LengthUnit lu)
  {
    set_vectors(vectors, lu);
    set_metric();
    set_space_symmetry(s);
    set_point_symmetry();
  }

private:
  void set_vectors(const vector_t& v, const vector_t& a, LengthUnit lu, AngleUnit angle_unit){
    AngleUnit au{angle_unit};

    if (AngleUnit::not_provided == au){
      auto minmax = std::minmax_element(a.begin(), a.end());
      if (*minmax.first < T(0))
        throw std::logic_error("Inter-facial cell angles must be positive but " + std::to_string(*minmax.first) + " found");
      au = (*minmax.second < math::pi) ? AngleUnit::radian : AngleUnit::degree;
    }

    vector_t cos, sin;
    switch (au) {
    case AngleUnit::degree:
      for (int i=0; i<3; ++i) std::tie(cos[i], sin[i]) = math::cos_and_sin_d(a[i]);
      break;
    case AngleUnit::radian:
      for (int i=0; i<3; ++i) std::tie(cos[i], sin[i]) = math::cos_and_sin(a[i]);
      break;
    default:
      for (int i=0; i<3; ++i) std::tie(cos[i], sin[i]) = math::cos_and_sin(math::pi * a[i]);
    }
    // unit-length parallelepiped volume:
    T v_sum{0}, v_prod{2};
    for (int i=0; i<3; ++i){
      v_sum += cos[i] * cos[i];
      v_prod *= cos[i];
    }
    T vol = std::sqrt(1 - v_sum + v_prod); // sqrt(1 - sum() + 2 * prod())
    T a_star, b_star, c_star, cb_star, cg_star, sb_star, sg_star, ca, two_pi_over_c;
    switch (lu) {
    case LengthUnit::angstrom:
      {
      // is this a bad idea?
        a_star = math::two_pi * sin[0] / v[0] / vol;
        b_star = math::two_pi * sin[1] / v[1] / vol;
        c_star = math::two_pi * sin[2] / v[2] / vol;
//        ca_star = (cos[1]*cos[2] - cos[0]) / (sin[1] * sin[2]);
        cb_star = (cos[2]*cos[0] - cos[1]) / (sin[2] * sin[0]);
        cg_star = (cos[0]*cos[1] - cos[2]) / (sin[0] * sin[1]);
        sb_star = std::sin(std::acos(cb_star));
        sg_star = std::sin(std::acos(cg_star));
        ca = cos[0];
        two_pi_over_c = math::two_pi / v[2];
      }
      break;
    case LengthUnit::inverse_angstrom:
      {
        a_star = v[0];
        b_star = v[1];
        c_star = v[2];
//        ca_star = cos[0];
        cb_star = cos[1];
        cb_star = cos[2];
        sb_star = sin[1];
        sg_star = sin[2];
        ca = (cos[1]*cos[2] - cos[0]) / (sin[1] * sin[2]);
        two_pi_over_c = (vol * v[2]) / sin[2];
      }
      break;
    default:
      throw std::logic_error("The length unit must be angstrom or inverse angstrom!");
    }
    // the B-matrix as in Acta Cryst. (1967). 22, 457 [with the 2pi convention]
    // http://dx.doi.org/10.1107/S0365110X67000970
    matrix_t B {
        a_star, b_star * cg_star, c_star * cb_star,
        T(0),   b_star * sg_star, -c_star * sb_star * ca,
        T(0),   T(0),             two_pi_over_c
    };
    set_vectors(B, LengthUnit::inverse_angstrom);
  }
  void set_vectors(const matrix_t& m, LengthUnit lu){
    matrix_t inv_m;
    utils::matrix_inverse(inv_m.data(), m.data());
    for (auto & x: inv_m) x *= math::two_pi;
    inv_m = transpose(inv_m);
    if (LengthUnit::angstrom == lu){
      _vectors = m;
      _reciprocal = inv_m;
    } else if (LengthUnit::inverse_angstrom == lu) {
      _vectors = inv_m;
      _reciprocal = m;
    }
  }
  void set_metric() {
    matrix_t vt = transpose(_vectors);
    utils::multiply_matrix_matrix(_metric.data(), _vectors.data(), vt.data());
  }
  void set_space_symmetry(const std::string& s, const std::string& c=""){
    auto no = string_to_hall_number(s, c);
    if (no > 0 && no < 531){
      auto sp = Spacegroup(no);
      _bravais = sp.get_bravais_type();
      _space = sp.get_spacegroup_symmetry();
      _point = sp.get_pointgroup_symmetry();
    } else {
      auto hs = HallSymbol(s);
      Symmetry gens;
      if (hs.validate()){
        gens = hs.get_generators();
        _bravais = hs.getl();
      } else {
        gens.from_ascii(s);
      }
      spacegroup_symmetry(gens);
    }
  }
  void set_point_symmetry(){
    _point = PointSymmetry(get_unique_rotations(_space.getallr(), 0));
  }

public:

  [[nodiscard]] matrix_t real_basis_vectors() const {return _vectors;}
  [[nodiscard]] matrix_t reciprocal_basis_vectors() const {return _reciprocal;}

  [[nodiscard]] Bravais bravais() const {return _bravais;}
  [[nodiscard]] Symmetry spacegroup_symmetry() const {return _space;}
  [[nodiscard]] PointSymmetry pointgroup_symmetry() const {return _point;}
  [[nodiscard]] Basis basis() const {return _basis;}

  [[nodiscard]] matrix_t metric() const {return _metric;}
  [[nodiscard]] matrix_t covariant_metric() const {return _metric;}
  [[nodiscard]] matrix_t contravariant_metric() const {
    matrix_t c;
    brille::utils::matrix_inverse(c.data(), _metric.data());
    return c;
  }

//  [[nodiscard]] T real_a() const {
//    T a{0};
//    for (int i=0; i<3; ++i) a += _vectors[3*i] * _vectors[3*i];
//    return std::sqrt(a);
//  }
//  [[nodiscard]] T real_b() const {
//    T b{0};
//    for (int i=0; i<3; ++i) b += _vectors[1+3*i] * _vectors[1+3*i];
//    return std::sqrt(b);
//  }
//  [[nodiscard]] T real_c() const {
//    T c{0};
//    for (int i=0; i<3; ++i) c+= _vectors[2+3*i] * _vectors[2+3*i];
//    return std::sqrt(c);
//  }

//  [[nodiscard]] T real_alpha() const;
//  [[nodiscard]] T real_beta() const;
//  [[nodiscard]] T real_gamma() const;
//  [[nodiscard]] T real_volume() const;
//  [[nodiscard]] T star_a() const;
//  [[nodiscard]] T star_b() const;
//  [[nodiscard]] T star_c() const;
//  [[nodiscard]] T star_alpha() const;
//  [[nodiscard]] T star_beta() const;
//  [[nodiscard]] T star_gamma() const;
//  [[nodiscard]] T star_volume() const;
//  [[nodiscard]] vector_t lengths(LengthUnit) const;
//  [[nodiscard]] vector_t angles(LengthUnit, AngleUnit) const;

  template<class R>
  [[nodiscard]] bool is_same(const Lattice<R>& o) const {
    // two lattices are *the* same if their basis vectors are the same
    auto bv = o.real_basis_vectors();
    auto rbv = o.reciprocal_basis_vectors();
    return approx::equal(bv, _vectors) && approx::equal(rbv, _reciprocal);
  }
  template<class R>
  [[nodiscard]] bool is_equivalent(const Lattice<R>& o) const {
    return is_same(o);
  }
  template<class R>
  [[nodiscard]] int is_permuted(const Lattice<R>& o) const {
    return is_same(o);
  }
  template<class R>
  bool operator == (const Lattice<R>& l) const {return is_same(l);}
  template<class R>
  bool operator != (const Lattice<R>& l) const {return !is_same(l);}

  Symmetry spacegroup_symmetry(const Symmetry& gens){
    _space = gens.generate();
    _bravais = _space.getcentring();
    _point = PointSymmetry(get_unique_rotations(_space.getallr(), 0));
    return _space;
  }
  [[nodiscard]] bool has_space_inversion() const {return _point.has_space_inversion();}
  [[nodiscard]] bool is_triclinic() const {return _point.higher(1).size() == 0u;}

  Basis basis(const std::vector<vector_t>& pos, const std::vector<ind_t>& typ){
    _basis = Basis(pos, typ);
    return _basis;
  }
  Basis basis(const Basis& b){
    _basis = b;
    return _basis;
  }

  [[nodiscard]] matrix_t to_xyz(LengthUnit lu){
    switch (lu) {
    case LengthUnit::angstrom: return _vectors;
    case LengthUnit::inverse_angstrom: return _reciprocal;
    default:
      throw std::runtime_error("Not implemented");
    }
  }
  [[nodiscard]] matrix_t from_xyz(LengthUnit lu){
    return transpose(to_xyz(lu));
  }
};


template<class T, class... Args>
Lattice<T> Direct(const std::array<T,3> & v, const std::array<T,3> & a, Args ... args){
  return {v, a, args..., LengthUnit::angstrom};
}
template<class T, class... Args>
Lattice<T> Reciprocal(const std::array<T,3> & v, const std::array<T,3> & a, Args ... args){
  return {v, a, args..., LengthUnit::inverse_angstrom};
}
//template<class T, class... Args>
//Lattice<T> Direct(const std::initializer_list<T> & v, const std::initializer_list<T> & a, Args ... args){
//  return {v, a, args..., LengthUnit::angstrom};
//}
//template<class T, class... Args>
//Lattice<T> Reciprocal(const std::initializer_list<T> & v, const std::initializer_list<T> & a, Args ... args){
//  return {v, a, args..., LengthUnit::inverse_angstrom};
//}


template<class T, class... Args>
Lattice<T> Direct(const std::array<T,9> & m, Args ... args){
  return {m, args..., LengthUnit::angstrom};
}
template<class T, class... Args>
Lattice<T> Reciprocal(const std::array<T,9> & m, Args ... args){
  return {m, args..., LengthUnit::inverse_angstrom};
}

}

#endif