#ifndef BRILLE_LATTICE_DUAL_HPP_
#define BRILLE_LATTICE_DUAL_HPP_

#include <assert.h>
#include <utility>
// #include <vector>
#include "enums.hpp"
#include "primitive.hpp"
#include "basis.hpp"
#include "array_.hpp"
#include "hdf_interface.hpp"
#include "math.hpp"
#include "hall_symbol.hpp"
#include "linear_algebra.hpp"

namespace brille::lattice {

template<class T>
class Lattice{
public:
  using matrix_t = std::array<T,9>;
  using vector_t = std::array<T,3>;
private:
  matrix_t _real_vectors; //! Real space *column* basis vectors
  matrix_t _reciprocal_vectors;
  matrix_t _real_metric; //! The metric tensor G
  matrix_t _reciprocal_metric; //! Avoid continually recalculating inverse(G) * 4 pi^2
  Bravais _bravais;
  Symmetry _space;
  PointSymmetry _point;
  Basis _basis;
public:
  explicit Lattice()
    : _real_vectors({{1, 0, 0, 0, 1, 0, 0, 0, 1}}),
      _reciprocal_vectors({{math::two_pi, 0, 0, 0, math::two_pi, 0, 0, 0, math::two_pi}}),
      _real_metric({{1, 0, 0, 0, 1, 0, 0, 0, 1}}),
      _reciprocal_metric({{math::two_pi * math::two_pi, 0, 0, 0, math::two_pi * math::two_pi, 0, 0, 0, math::two_pi * math::two_pi}}),
      _bravais(Bravais::P),
      _space({Motion<int, double>()}), _point(), _basis()
    {}

  Lattice(matrix_t v, matrix_t r, matrix_t m, matrix_t rm, Bravais b, Symmetry s, PointSymmetry p, Basis a)
      : _real_vectors{std::move(v)},
        _reciprocal_vectors{std::move(r)},
        _real_metric{std::move(m)},
        _reciprocal_metric{std::move(rm)},
        _bravais(b), _space(std::move(s)), _point(std::move(p)), _basis(std::move(a))
  {}
//
  Lattice(const vector_t & lengths, const vector_t & angles, const Symmetry s, const LengthUnit lu, const AngleUnit au=AngleUnit::not_provided)
  {
    set_vectors(lengths, angles, lu, au);
    set_metrics();
    spacegroup_symmetry(s);
    set_point_symmetry();
  }
  Lattice(const matrix_t & vectors, const Symmetry s, const LengthUnit lu)
  {
    set_vectors(vectors, lu);
    set_metrics();
    spacegroup_symmetry(s);
    _bravais = _space.getcentring();
    set_point_symmetry();
  }
  Lattice(const vector_t & lengths, const vector_t & angles, const std::string& s, const std::string& c, const LengthUnit lu, const AngleUnit au=AngleUnit::not_provided)
  {
    set_vectors(lengths, angles, lu, au);
    set_metrics();
    set_space_symmetry(s, c);
    _bravais = _space.getcentring();
    set_point_symmetry();
  }
  Lattice(const vector_t & lengths, const vector_t & angles, const std::string& s, const LengthUnit lu, const AngleUnit au=AngleUnit::not_provided)
  {
    set_vectors(lengths, angles, lu, au);
    set_metrics();
    set_space_symmetry(s);
    _bravais = _space.getcentring();
    set_point_symmetry();
  }
  Lattice(const matrix_t & vectors, const MatrixVectors mv, const std::string& s, const std::string& c, const LengthUnit lu)
  {
    set_vectors(MatrixVectors::column ==mv ? vectors : transpose(vectors), lu);
    set_metrics();
    set_space_symmetry(s, c);
    _bravais = _space.getcentring();
    set_point_symmetry();
  }
  Lattice(const matrix_t & vectors, const MatrixVectors mv, const std::string& s, const LengthUnit lu)
  {
    set_vectors(MatrixVectors::column ==mv ? vectors : transpose(vectors), lu);
    set_metrics();
    set_space_symmetry(s);
    _bravais = _space.getcentring();
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

    auto dual = [&](){
      vector_t dual_v, dual_cos, dual_sin;
      // e.g., a* = 2π b×c / a⋅(b×c) → 2π sin(α) / (a * â(b̂×ĉ))
      for (size_t i=0; i<3; ++i) dual_v[i] = math::two_pi * sin[i] / v[i] / vol;
      for (size_t i=0; i<3; ++i) {
        size_t j{(i+1u)%3u}, k{(i+2u)%3u};
        dual_cos[i] = (cos[j]*cos[k] - cos[i]) / (sin[j] * sin[k]);
        dual_sin[i] = std::sqrt(1 - dual_cos[i] * dual_cos[i]); // ok since angles are less than pi
      }
      return std::make_tuple(dual_v, dual_cos, dual_sin);
    };
    vector_t dv, dcos, dsin, rv, rcos, rsin;
    switch (lu) {
    case LengthUnit::angstrom:
      {
        std::tie(rv, rcos, rsin) = dual();
        dv = v; dcos = cos; dsin = sin;
      }
      break;
    case LengthUnit::inverse_angstrom:
      {
        std::tie(dv, dcos, dsin) = dual();
        rv = v; rcos = cos; rsin = sin;
      }
      break;
    default:
      throw std::logic_error("The length unit must be angstrom or inverse angstrom!");
    }
    // the B-matrix as in Acta Cryst. (1967). 22, 457 [with the 2pi convention]
    // http://dx.doi.org/10.1107/S0365110X67000970
    matrix_t B {
        rv[0], rv[1] * rcos[2],  rv[2] * rcos[1],
        T(0),  rv[1] * rsin[2], -rv[2] * rsin[1] * dcos[0],
        T(0),  T(0),             math::two_pi / dv[2]
    };
    set_vectors(B, LengthUnit::inverse_angstrom);
  }
  /*! \brief Set the real and reciprocal basis vector matrices
   *
   * @param m Either the real space or reciprocal space basis vector matrix,
   *          which has the three basis vectors as its *columns*.
   * @param lu An enumerated length unit to indicate if the passed matrix
   *           represents the real (LengthUnit::angstrom) or reciprocal
   *           (LengthUnit::inverse_angstrom) basis vectors.
   */
  void set_vectors(const matrix_t& m, LengthUnit lu){
    auto inv_m = linear_algebra::mat_inverse(m);
    for (auto & x: inv_m) x *= math::two_pi;
    inv_m = transpose(inv_m);
    if (LengthUnit::angstrom == lu){
      _real_vectors = m;
      _reciprocal_vectors = inv_m;
    } else if (LengthUnit::inverse_angstrom == lu) {
      _real_vectors = inv_m;
      _reciprocal_vectors = m;
    }
  }
  void set_metrics() {
    // the two metrics are mutually inverse with a factor of 4 * pi^2
    _real_metric = linear_algebra::mul_mat_mat(_real_vectors, transpose(_real_vectors));
    _reciprocal_metric = linear_algebra::mul_mat_mat(_reciprocal_vectors, transpose(_reciprocal_vectors));
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
  [[nodiscard]] std::string to_string(const LengthUnit lu, const AngleUnit au=AngleUnit::degree) const {
    std::string rep;
    rep += "(" + my_to_string(lengths(lu)) + ")" + (lu == LengthUnit::angstrom ? "Å " : "Å⁻¹ ");
    rep += "(" + my_to_string(angles(lu, au)) + ")" + (AngleUnit::degree == au ? "°" : AngleUnit::radian ==au ? "rad" : "π");
    return rep;
  }
  [[nodiscard]] std::string to_verbose_string(const AngleUnit au=AngleUnit::degree) const {
    std::string rep;
    rep += bravais_string(_bravais) + " Lattice with " + std::to_string(_space.size()) + " space group elements\n";
    rep += "      Direct{" + to_string(LengthUnit::angstrom, au) + "}\n";
    rep += "  Reciprocal{" + to_string(LengthUnit::inverse_angstrom, au) + "}\n";
    return rep;
  }

  [[nodiscard]] matrix_t real_basis_vectors() const {return _real_vectors;}
  [[nodiscard]] matrix_t reciprocal_basis_vectors() const {return _reciprocal_vectors;}

  [[nodiscard]] Bravais bravais() const {return _bravais;}
  [[nodiscard]] Symmetry spacegroup_symmetry() const {return _space;}
  [[nodiscard]] PointSymmetry pointgroup_symmetry() const {return _point;}
  [[nodiscard]] Basis basis() const {return _basis;}

  [[nodiscard]] matrix_t metric(const LengthUnit lu) const {
    switch(lu){
      case LengthUnit::angstrom: return _real_metric;
      case LengthUnit::inverse_angstrom: return _reciprocal_metric;
      default:
        throw std::logic_error("Not implemented length unit");
    }
  }
//  [[nodiscard]] matrix_t covariant_metric() const {return _real_metric;}
//  [[nodiscard]] matrix_t contravariant_metric() const {
//    return linear_algebra::mat_inverse(_real_metric);
//  }

  template<class I>
  [[nodiscard]] std::enable_if_t<std::is_integral_v<I>, vector_t>
  vector(LengthUnit lu, I i) const {
    assert(0u <= i && i <= 3u);
    switch (lu) {
    case LengthUnit::angstrom:
      return {_real_vectors[i], _real_vectors[i + 3], _real_vectors[i + 6]};
    case LengthUnit::inverse_angstrom:
      return {_reciprocal_vectors[i], _reciprocal_vectors[i + 3], _reciprocal_vectors[i + 6]};
    default:
      throw std::logic_error("Not implemented length unit");
    }
  }
  template<class I>
  [[nodiscard]] std::enable_if_t<std::is_integral_v<I>, T>
  length(LengthUnit lu, I i) const {
    return linear_algebra::norm(vector(lu, i));
  }
  template<class I>
  [[nodiscard]] std::enable_if_t<std::is_integral_v<I>, T>
  angle(LengthUnit lu, I i, AngleUnit au = AngleUnit::radian) const {
    auto vj = vector(lu, (i+1) % 3u);
    auto vk = vector(lu, (i+2) % 3u);
    auto cos_angle_i = linear_algebra::dot(vj, vk) / linear_algebra::norm(vj) / linear_algebra::norm(vk);
    switch (au) {
    case AngleUnit::radian:
      return std::acos(cos_angle_i);
    case AngleUnit::degree:
      return std::acos(cos_angle_i) / math::pi * T(180);
    case AngleUnit::pi:
      return std::acos(cos_angle_i) / math::pi;
    default:
      throw std::logic_error("Not implemented angle unit");
    }
  }
  [[nodiscard]] vector_t lengths(LengthUnit lu) const {
    return {length(lu, 0u), length(lu, 1u), length(lu, 2u)};
  }
  [[nodiscard]] vector_t angles(LengthUnit lu, AngleUnit au = AngleUnit::radian) const {
    // this could be more efficient by reusing vi, vj, vk; but maybe its ok.
    return {angle(lu, 0u, au), angle(lu, 1u, au), angle(lu, 2u, au)};
  }
  [[nodiscard]] T volume(LengthUnit lu) const {
    switch (lu){
    case LengthUnit::angstrom:
      return utils::matrix_determinant(_real_vectors.data());
    case LengthUnit::inverse_angstrom:
      return utils::matrix_determinant(_reciprocal_vectors.data());
    default:
      throw std::logic_error("Not implemented length unit for volume");
    }
  }

  template<class R>
  [[nodiscard]] bool is_same(const Lattice<R>& o) const {
    // two lattices are *the* same if their basis vectors are the same
    auto bv = o.real_basis_vectors();
    auto rbv = o.reciprocal_basis_vectors();
    return approx::equal(bv, _real_vectors) && approx::equal(rbv, _reciprocal_vectors);
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

  [[nodiscard]] matrix_t to_xyz(LengthUnit lu) const {
    /*  A vector expressed in the orthonormal coordinate system, e⃗, is related to a vector expressed as coordinates
     *  of the real lattice basis, x⃗, by e⃗ = A x⃗; where A is the column-vector matrix made of the real basis vectors
     *  expressed in e⃗.
     *  Similarly, for a vector expressed as coordinates of the reciprocal lattice basis, q⃗, we have e⃗ = B q⃗.
     *
     *  This method returns the stored basis vectors matrix to take x⃗ or q⃗ to e⃗.
     * */
    switch (lu) {
    case LengthUnit::angstrom: return _real_vectors;
    case LengthUnit::inverse_angstrom: return _reciprocal_vectors;
    default: throw std::runtime_error("Not implemented");
    }
  }
  [[nodiscard]] matrix_t from_xyz(LengthUnit lu){
    /* Similarly to `to_xyz` we make use of the relationships e⃗ = A x⃗ and e⃗ = B q⃗ to enable finding x⃗ or q⃗ from e⃗.
     * Inverting A or B and multiplying from the left we find,  x⃗ = A⁻¹ e⃗   and   q⃗ = B⁻¹ e⃗.
     * However, we can avoid taking the inversion of A or B since we remember that the two sets of basis vectors
     * are related by B = 2π (A⁻¹)ᵀ. [To find B⁻¹, B⁻¹=(2π (A⁻¹)ᵀ)⁻¹=(Aᵀ/2π) --> BB⁻¹= 2π (A⁻¹)ᵀ Aᵀ / 2π= 𝟙]
     *    x⃗ = A⁻¹ e⃗ = (Bᵀ/2π) e⃗   and   q⃗ = B⁻¹ e⃗ = (Aᵀ/2π) e⃗
     * */
    matrix_t mat;
    switch (lu) {
      case LengthUnit::angstrom: mat = _reciprocal_vectors; break;
      case LengthUnit::inverse_angstrom: mat = _real_vectors; break;
      default: throw std::runtime_error("Not implemented");
    }
    for (auto & x: mat) x /= math::two_pi;
    return transpose(mat);
  }

#ifdef USE_HIGHFIVE
  template<class H>
  std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, bool>
  to_hdf(H& obj, const std::string& entry) const {
    auto group = overwrite_group(obj, entry);
    group.createAttribute("real_space_vectors", _real_vectors);
    group.createAttribute("reciprocal_space_vectors", _reciprocal_vectors);
    group.createAttribute("real_metric", _real_metric);
    group.createAttribute("reciprocal_metric", _reciprocal_metric);
    group.createAttribute("bravais", _bravais);
    bool ok{true};
    ok &= _space.to_hdf(group, "spacegroup_symmetry");
    ok &= _point.to_hdf(group, "pointgroup_symmetry");
    ok &= _basis.to_hdf(group, "basis");
    return ok;
  }
  template<class H>
  static std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, Lattice<T>>
  from_hdf(H& obj, const std::string& entry){
    auto group = obj.getGroup(entry);
    matrix_t real, reciprocal, real_metric, reciprocal_metric;
    group.getAttribute("real_space_vectors").read(real);
    group.getAttribute("reciprocal_space_vectors").read(reciprocal);
    group.getAttribute("real_metric").read(real_metric);
    group.getAttribute("reciprocal_metric").read(reciprocal_metric);
    Bravais L;
    group.getAttribute("bravais").read(L);
    auto spg = Symmetry::from_hdf(group, "spacegroup_symmetry");
    auto ptg = PointSymmetry::from_hdf(group, "pointgroup_symmetry");
    auto bas = Basis::from_hdf(group, "basis");
    return {real, reciprocal, real_metric, reciprocal_metric, L, spg, ptg, bas};
  }
#endif

  Lattice<T> primitive() const {
    PrimitiveTransform P(_bravais);
    if (P.does_anything()){
      // The transformation matrix P gives us the primitive basis column-vector
      // matrix Aₚ from the standard basis column-vector matrix Aₛ by
      // Aₚ = Aₛ P.
      // The PrimitiveTransform object contains 6*P (so that it's integer valued)
      auto pv = linear_algebra::mul_mat_mat(_real_vectors, P.get_6P());
      for (auto & x: pv) x /= T(6);
      // calculating the new metric tensor could be done from _real_vectors, but
      // we might as well use the new values in pv
      auto pvm = linear_algebra::mul_mat_mat(pv, transpose(pv));

      /* We have a column-vector matrix of reciprocal lattice basis vectors
       * related to the real space matrix by
       *            Bₛ= 2π(Aₛ⁻¹)ᵀ
       *      (Bₛ/2π)ᵀ= Aₛ⁻¹
       *    (Bₛᵀ/2π)⁻¹= Aₛ
       *     2π(Bₛᵀ)⁻¹= Aₛ, which is of course true for the primitive lattice
       *
       *        2π(Bₚᵀ)⁻¹ = 2π(Bₛᵀ)⁻¹P
       *          (Bₚ⁻¹)ᵀ = (Bₛ⁻¹)ᵀP
       *          Bₚ⁻¹ = Pᵀ Bₛ⁻¹
       *          Bₚ = Bₛ (Pᵀ)⁻¹                                              */
      auto pr = linear_algebra::mul_mat_mat(_reciprocal_vectors, P.get_invPt());
      auto prm = linear_algebra::mul_mat_mat(pr, transpose(pr));

      // strictly we should change the spacegroup and pointgroup symmetry
      // representations, plus the atom basis representation... but that seems
      // overkill for now. FIXME think about this more.
      return {pv, pr, pvm, prm, _bravais, _space, _point, _basis};
    }
    return *this;
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