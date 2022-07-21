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
#include "pointgroup.hpp"

namespace brille::lattice {
template<class T, class S=std::common_type_t<T, double>> std::array<S, 9>
standard_orientation_matrix(const std::array<T, 9>& bv){
  /* for arbitrary basis vectors stored as a row-flattened column vector matrix
   *  [[ax, bx, cx],
   *   [ay, by, cy],
   *   [az, bz, cz]] --> {ax, bx, cx, ay, by, cy, az, bz, cz}
   *
   *   Find the rotation matrix which takes us to the standard upper triangular
   *   form with flattened matrix {a, bx, cx, 0, by, cy, 0, 0, cz}
   * */
	// Avoid (0, -1) -> pi in std::atan2:
	auto arctan = [](T x, T y){return std::abs(x) ? std::atan2(x, y) : 0.;};
  // tan(mu) = -ay / az
  S mu = arctan(-bv[3], bv[6]);
  auto [c_mu, s_mu] = math::cos_and_sin(mu);

  // tan(nu) = -az / (ax * cos(mu))
  S nu = arctan(-bv[6], bv[0] * c_mu);
  auto [c_nu, s_nu] = math::cos_and_sin(nu);

  // tan(xi) = (bx sin(nu) - by cos(nu) sin(mu) + bz cos(nu) cos(mu)) / (by cos(mu) + bz sin(mu))
  S num = bv[1] * s_nu - bv[4] * c_nu * s_mu + bv[7] * c_nu * c_mu;
  S den = bv[4] * c_mu + bv[7] * s_mu;
  S xi = arctan(num, den);
  auto [c_xi, s_xi] = math::cos_and_sin(xi);

  // Rotation matrix, verified 'by hand' multiple times.
  std::array<S,9> rotation_matrix{
      {c_nu,                              s_nu * s_mu,                      -s_nu * c_mu,
       s_xi * s_nu,  c_xi * c_mu - s_xi * c_nu * s_mu,  c_xi * s_mu + s_xi * c_nu * c_mu,
       c_xi * s_nu, -s_xi * c_mu - c_xi * c_nu * s_mu, -s_xi * s_mu + c_xi * c_nu * c_mu}
  };
  return rotation_matrix;
}

template<class T, class S=std::common_type_t<T, double>> S unit_parallelpiped_volume(const std::array<T,3>& cosines){
  // unit-length parallelepiped volume:
  S cos_sum{0}, cos_prod{2};
  for (int i=0; i<3; ++i){
    cos_sum += cosines[i] * cosines[i];
    cos_prod *= cosines[i];
  }
  S unit_volume_squared = S(1) - cos_sum + cos_prod; // 1 - sum() + 2 * prod()
  if (unit_volume_squared <= S(0)) {
    std::string msg = "A lattice with interfacial cosines {";
    for (const auto & c: cosines) msg += " " + std::to_string(c);
    msg += "} has ";
    msg += (unit_volume_squared < 0 ? "imaginary" : "no");
    msg += " volume and can therefore not represent a physical lattice.";
    throw std::invalid_argument(msg);
  }
  S unit_volume = std::sqrt(unit_volume_squared); // sqrt(1 - sum() + 2 * prod())
  return unit_volume;
}

template<class T, class S=std::common_type_t<T, double>> std::tuple<std::array<S,3>, std::array<S,3>, std::array<S,3>>
dual_lattice_parameters(const std::array<T,3>& lengths, const std::array<T,3>& cosines, const std::array<T,3>& sines){
  auto unit_volume = unit_parallelpiped_volume(cosines);

  std::array<S,3> dual_len, dual_cos, dual_sin;
  for (size_t i=0; i<3; ++i) {
    size_t j{(i+1u)%3u}, k{(i+2u)%3u}; // cyclical indexing for angles
    // e.g., a* = 2π b×c / a⋅(b×c) → 2π sin(α) / (a * â(b̂×ĉ))
    dual_len[i] = math::two_pi * sines[i] / lengths[i] / unit_volume;
    dual_cos[i] = (cosines[j]*cosines[k] - cosines[i]) / (sines[j] * sines[k]);
    dual_sin[i] = std::sqrt(S(1) - dual_cos[i] * dual_cos[i]); // ok since angles are less than pi
  }
  return std::make_tuple(dual_len, dual_cos, dual_sin);
}

template<class T, class S=std::common_type_t<T, double>> std::tuple<std::array<S,3>, std::array<S,3>>
inter_facial_angles_to_cosines_sines(const std::array<T,3>& angles, AngleUnit au){
  if (AngleUnit::not_provided == au){
    auto minmax = std::minmax_element(angles.begin(), angles.end());
    if (*minmax.first < T(0))
      throw std::logic_error("Inter-facial cell angles must be positive but " + std::to_string(*minmax.first) + " found");
    au = (*minmax.second < math::pi) ? AngleUnit::radian : AngleUnit::degree;
  }
  std::array<S,3> cos, sin;
  switch (au) {
  case AngleUnit::degree:
    for (int i=0; i<3; ++i) std::tie(cos[i], sin[i]) = math::cos_and_sin_d(angles[i]);
    break;
  case AngleUnit::radian:
    for (int i=0; i<3; ++i) std::tie(cos[i], sin[i]) = math::cos_and_sin(angles[i]);
    break;
  default:
    for (int i=0; i<3; ++i) std::tie(cos[i], sin[i]) = math::cos_and_sin(math::pi * angles[i]);
  }
  return std::make_tuple(cos, sin);
}

template<class T> std::tuple<size_t, size_t> not_right_angles(const std::array<T,3>& cos, const std::array<T,3>& sin){
  using approx_float::scalar;
  size_t index{0}, count{0};
  for (size_t i=0; i<3u; ++i) if (!scalar(T(1), cos[i]) || !scalar(T(0), sin[i])) {
    index = i;
    ++count;
  }
  return std::make_tuple(index, count);
}

template<class T> size_t count_matching(const std::array<T,3>& len){
  using approx_float::scalar;
  auto m01 = scalar(len[0], len[1]);
  auto m02 = scalar(len[0], len[2]);
  if (m01 && m02) return 3u;
  if (m01 || m02 || scalar(len[1], len[2])) return 2u;
  return 0u;
}

template<class T> std::array<T,9> metric_from_column_vectors(const std::array<T,9>& vectors){
  // The column vector matrix A has an associated space metric given by AᵀA
  //
  //    [[ ax bx cx ]
  // A = [ ay by cy ]
  //     [ az bz cz ]]
  //
  //      [[ ax ay az ] [[ ax bx cx ]  [[ a²   a⋅b  a⋅c ]
  // AᵀA = [ bx by bz ]  [ ay by cy ] = [ b⋅a  b²   b⋅c ]
  //       [ cx cy cz ]] [ az bz cz ]   [ c⋅a  c⋅b  c²  ]
  return linear_algebra::mul_mat_mat(transpose(vectors), vectors);
}

template<class T>
class Impl{
public:
  using matrix_t = std::array<T,9>;
  using vector_t = std::array<T,3>;
private:
  matrix_t _real_vectors; //! Real space *column* basis vectors
  matrix_t _reciprocal_vectors;
  matrix_t _real_metric; //! The metric tensor G
  matrix_t _reciprocal_metric; //! Avoid continually recalculating inverse(G) * 4 pi^2
  Bravais _bravais = Bravais::_;
  Symmetry _space;
  PointSymmetry _point;
  Basis _basis;
public:
  Impl() = default;

  /*! \brief Fully specified construction
   * @note Intended (mostly) to be used by the static HDF from_file reader
   *
   * @param v Real space flattened column-basis-vector matrix in Å
   * @param r Reciprocal space flattened column-basis-vector matrix in Å⁻¹
   * @param m Real space flattened metric in Å²
   * @param rm Reciprocal space flattened metric in Å⁻²
   * @param b Centring information (extraneous, since its in the Symmetry too)
   * @param s Real space lattice symmetry operations
   * @param p Real space point symmetry operations
   * @param a Real space lattice atom basis
   */
  Impl(matrix_t v, matrix_t r, matrix_t m, matrix_t rm, Bravais b, Symmetry s, PointSymmetry p, Basis a, bool snap_to_symmetry=true)
      : _real_vectors{std::move(v)},
        _reciprocal_vectors{std::move(r)},
        _real_metric{std::move(m)},
        _reciprocal_metric{std::move(rm)},
        _bravais(b), _space(std::move(s)), _point(std::move(p)), _basis(std::move(a))
  {
    snap_basis_to_symmetry(snap_to_symmetry);
  }
  /*! \brief Lattice parameters and Symmetry constructor
   *
   * @param lu Indicates if the lengths are for the real or reciprocal space basis vectors
   * @param lengths Basis vector lengths in units given by `lu`
   * @param angles Inter-basis-vector angles in units given by `au`
   * @param s The Symmetry operation of the real space lattice
   * @param au Indicates if the angles are expressed in units of degrees, radians, or fractions of pi-radians
   */
  Impl(const LengthUnit lu, const vector_t & lengths, const vector_t & angles, const Symmetry & s, Basis basis=Basis(), bool snap_to_symmetry=true, const AngleUnit au=AngleUnit::not_provided)
  : _basis(std::move(basis))
  {
    spacegroup_symmetry(s);
    _bravais = _space.getcentring();
    set_point_symmetry();
    set_vectors(lengths, angles, lu, au, snap_to_symmetry);
    set_metrics();
    snap_basis_to_symmetry(snap_to_symmetry);
  }
  /*! \brief Lattice basis vectors and Symmetry constructor
   *
   * @param lu Indicates if the vectors are that of the real or reciprocal space
   * @param vectors A flattened row-ordered 3x3 matrix of the basis vectors
   * @param mv Indicates if the pre-flattened matrix was row- or column-vectors
   * @param s The Symmetry operation of the real space lattice
   * @param basis Real space lattice atom basis
   * @param snap_to_symmetry Whether to enforce the symmetry operations on the basis and lattice parameters
   */
  Impl(const LengthUnit lu, const matrix_t & vectors, const MatrixVectors mv, const Symmetry & s, Basis basis=Basis(), bool snap_to_symmetry=true)
  : _basis(std::move(basis))
  {
    spacegroup_symmetry(s);
    _bravais = _space.getcentring();
    set_point_symmetry();
    set_vectors(MatrixVectors::column ==mv ? vectors : transpose(vectors), lu, snap_to_symmetry);
    set_metrics();
    snap_basis_to_symmetry(snap_to_symmetry);
  }
  /*! \brief Lattice parameters and Hermann-Maunguin spacegroup information constructor
   *
   * @param lu Indicates if the lengths are for the real or reciprocal space basis vectors
   * @param lengths Basis vector lengths in units given by `lu`
   * @param angles Inter-basis-vector angles in units given by `au`
   * @param s short or long Hermann-Mauguin group name, as used in the International Tables of Crystallography
   * @param c Hermann-Mauguin centering/axis 'choice', only required if non-default
   * @param basis Real space lattice atom basis
   * @param snap_to_symmetry Whether to enforce the symmetry operations on the basis and lattice parameters
   * @param au Indicates if the angles are expressed in units of degrees, radians, or fractions of pi-radians
   */
  Impl(const LengthUnit lu, const vector_t & lengths, const vector_t & angles, const std::string& s, const std::string& c, Basis basis=Basis(), bool snap_to_symmetry=true, const AngleUnit au=AngleUnit::not_provided)
  : _basis(std::move(basis))
  {
    set_space_symmetry(s, c);
    _bravais = _space.getcentring();
    set_point_symmetry();
    set_vectors(lengths, angles, lu, au, snap_to_symmetry);
    set_metrics();
    snap_basis_to_symmetry(snap_to_symmetry);
  }
  /*! \brief Lattice parameters and string-encoded spacegroup information constructor
   *
  * @param lu Indicates if the lengths are for the real or reciprocal space basis vectors
   *@param lengths Basis vector lengths in units given by `lu`
  * @param angles Inter-basis-vector angles in units given by `au`
  * @param s Hall symbol, CIF xyz operations, or International Table group name
  * @param basis Real space lattice atom basis
  * @param snap_to_symmetry Whether to enforce the symmetry operations on the basis and lattice parameters
  * @param au Indicates if the angles are expressed in units of degrees, radians, or fractions of pi-radians
   */
  Impl(const LengthUnit lu, const vector_t & lengths, const vector_t & angles, const std::string& s, Basis basis=Basis(), bool snap_to_symmetry=true, const AngleUnit au=AngleUnit::not_provided)
  : _basis(std::move(basis))
  {
    set_space_symmetry(s);
    _bravais = _space.getcentring();
    set_point_symmetry();
    set_vectors(lengths, angles, lu, au, snap_to_symmetry);
    set_metrics();
    snap_basis_to_symmetry(snap_to_symmetry);
  }
  /*! \brief Lattice basis vectors and Hermann-Maunguin spacegroup information constructor
   *
   * @param lu Indicates if the lengths are for the real or reciprocal space basis vectors
   * @param vectors A flattened row-ordered 3x3 matrix of the basis vectors
   * @param mv Indicates if the pre-flattened matrix was row- or column-vectors
   * @param s short or long Hermann-Mauguin group name, as used in the International Tables of Crystallography
   * @param c Hermann-Mauguin centering/axis 'choice', only required if non-default
   * @param basis Real space lattice atom basis
   * @param snap_to_symmetry Whether to enforce the symmetry operations on the basis and lattice parameters
   */
  Impl(const LengthUnit lu, const matrix_t & vectors, const MatrixVectors mv, const std::string& s, const std::string& c, Basis basis=Basis(), bool snap_to_symmetry=true)
  : _basis(std::move(basis))
  {
    set_space_symmetry(s, c);
    _bravais = _space.getcentring();
    set_point_symmetry();
    set_vectors(MatrixVectors::column ==mv ? vectors : transpose(vectors), lu, snap_to_symmetry);
    set_metrics();
    snap_basis_to_symmetry(snap_to_symmetry);
  }
  /*! \brief Lattice basis vectors and string-encoded spacegroup information constructor
   *
   * @param lu Indicates if the vectors are of the real or reciprocal basis
   * @param vectors A flattened row-ordered 3x3 matrix of the basis vectors
   * @param mv Indicates if the pre-flattened matrix was row- or column-vectors
   * @param s Hall symbol, CIF xyz operations, or International Table group name
   * @param basis Real space lattice atom basis
   * @param snap_to_symmetry Whether to enforce the symmetry operations on the basis and lattice parameters
   */
  Impl(const LengthUnit lu, const matrix_t & vectors, const MatrixVectors mv, const std::string& s, Basis basis=Basis(), bool snap_to_symmetry=true)
  : _basis(std::move(basis))
  {
    set_space_symmetry(s);
    _bravais = _space.getcentring();
    set_point_symmetry();
    set_vectors(MatrixVectors::column ==mv ? vectors : transpose(vectors), lu, snap_to_symmetry);
    set_metrics();
    snap_basis_to_symmetry(snap_to_symmetry);
  }


private:
  void snap_basis_to_symmetry(const bool snap_to_symmetry){
    if (snap_to_symmetry) {
      auto success = _basis.snap_to(_space.getallm());
      if (!success) {
        std::string msg("Requested snap_to_symmetry failed with result \n");
        msg += _basis.to_string();
        throw std::runtime_error(msg);
      }
    }
  }
  bool snap_to_symmetry(bool is_real, std::array<T, 3>& v, std::array<T, 3>& c, std::array<T, 3>& s){
    /* Averaging is likely to reduce precision, so only do so if necessary */
    auto equal2 = [](const T x, const T y){ return (x == y) ? x : (x + y) / T(2); };
    auto equal3 = [](const T x, const T y, const T z){ return (x == y && y == z) ? x : (x + y + z) / T(3); };
    // check if centring can give us a hint:
    /* If the Lattice has spacegroup operations with centering is it *not*
     * primitive. If the centering indicates that it is Rhombohedral then
     * it is actually the Hexagonal setting of a rhombohedral lattice, so
     * two angles must be 90 degrees and the third 120 degrees.
     * If it is A,B,C Face centred, I Body centred, F face centred, then
     * the symmetry-constrained angles *must* be 90 degrees.
     * */
    // don't use the stored Bravais value in case of user error?
    bool centring_hint{_space.getcentring() != Bravais::P};
    bool second_hint{centring_hint && _space.getcentring() != Bravais::R};
    // determine if basis vectors are connected by symmetry operations
    std::vector<std::array<int,3>> eis {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    // (a-b, b-c, c-a), (a*-b*, b*-c*, c*-a*)
    bool main[3]{false, false, false}, dual[3]{false, false, false};
    for (size_t i=0; i<3u; ++i){
      const auto j{(i+1)%3u};
      // Real lattice vectors rotate like v' = R v
      main[i] = _point.mat_vec_links(eis[i], eis[j]);
      // Reciprocal lattice vectors rotate like q' = tr(R) q = tr(tr(q) R)
      dual[i] = _point.vec_mat_links(eis[i], eis[j]);
    }
    // if we're dealing with Reciprocal lattice parameters, then main and dual
    // need to be swapped!
    if (!is_real) std::swap(main, dual);
    // stash the original lengths and cosines to check for changes
    std::array<T,3> orig_v{v}, orig_cos{c};
    if (main[0] && main[1]) {
      /* All 'main' basis vectors are mapped to each other; so they must all  *
       * have the same length:                                                */
      v[0] = v[1] = v[2] = equal3(v[0], v[1], v[2]);
      if (dual[0] && dual[1]) {
        /* All 'dual' basis vectors are mapped, so the 'main' basis vector    *
         * angles must be the same as well. If the lattice is centred then
         * they all must be exactly 90 degrees, so the cos(angle)=0.0; !      */
        c[0] = c[1] = c[2] = centring_hint ? T(0) : equal3(c[0], c[1], c[2]);
        s[0] = s[1] = s[2] = centring_hint ? T(1) : std::sqrt(1 - c[0] * c[0]);
      } else for (size_t i=0; i<3u; ++i) if (dual[i]) {
        const auto j{(i + 1) % 3u};
        c[i] = c[j] = second_hint ? T(0) : equal2(c[i], c[j]);
        s[i] = s[j] = second_hint ? T(1) : std::sqrt(1 - c[i] * c[i]);
      }
    } else for (size_t i=0; i<3u; ++i) if (main[i]) {
      const auto j{(i + 1) % 3u};
      v[i] = v[j] = equal2(v[i], v[j]);
      if (dual[0] && dual[1]){
        throw std::runtime_error("How can more than two dual-vectors be mapped when only two real-vectors are mapped?");
      }
      for (size_t k=0; k<3u; ++k) if(k!=i && dual[k]){
        throw std::runtime_error("How can a dual-vector be mapped if its real-vector is not mapped?");
      }
      if (dual[i]) {
        if (centring_hint && _space.getcentring() == Bravais::R) {
          const auto k{(i + 2) % 3u};
          std::tie(c[k], s[k]) = math::cos_and_sin_d(is_real ? T(120) : T(60));
        }
        c[i] = c[j] = centring_hint ? T(0) : equal2(c[i], c[j]);
        s[i] = s[j] = centring_hint ? T(1) : std::sqrt(1 - c[i] * c[i]);
      }
    }

    bool changed{false};
    if (main[0] || main[1] || main[2]){
      for (size_t i=0; i<3u; ++i) {
        if (!approx_float::scalar(orig_v[i] - v[i], 0.)) changed = true;
      }
      for(size_t i=0; i<3u; ++i){
        if (!approx_float::scalar(orig_cos[i] - c[i], 0.)) changed = true;
      }
    }
    return changed; // signal if we did anything
  }
  bool snap_parameters_to_symmetry(bool is_real, std::array<T, 3>& v, std::array<T, 3>& c, std::array<T, 3>& s) {
    std::array<T,3> orig_v{v}, orig_cos{c};
    if (snap_to_symmetry(is_real, v, c, s)){
      std::ostringstream msg;
      msg << "Basis vector lengths, angles changed by (";
      for (size_t i=0; i<3u; ++i) {
        msg << " " << v[i] - orig_v[i];
      }
      msg << " )" << (is_real ? u8"Å" : u8"Å⁻¹") << ", (";
      for(size_t i=0; i<3u; ++i){
        msg << " " << (std::acos(c[i]) - std::acos(orig_cos[i])) / math::pi * 180.;
      }
      msg << u8")°";
      info_update(msg.str());
      return true;
    }
    return false;
  }
  void snap_basis_vectors_to_symmetry(){
    // calculate lattice parameters from already-set _reciprocal_vectors
		bool use_real{false};
		auto lu = use_real ? LengthUnit::angstrom : LengthUnit::inverse_angstrom;
    auto v = lengths(lu);
    auto [c, s] = inter_facial_angles_to_cosines_sines(angles(lu), AngleUnit::radian);
    if (snap_to_symmetry(use_real, v, c, s)) {
      // determine the re-orientation matrix necessary to align the real basis
      // vectors with the 'standard' orientation of a || x, b⋅z == 0, b⋅y > 0
      auto R = standard_orientation_matrix(use_real ? _real_vectors : _reciprocal_vectors);
      // build the upper-triangular basis vector column-vector matrix
      matrix_t ut {{
          v[0],  v[1] * c[2], v[2] * c[1],
          T(0) , v[1] * s[2], v[2] * (c[0] - c[2]*c[1])/s[2],
          T(0) , T(0)       , v[2] * unit_parallelpiped_volume(c) / s[2]
      }};
      // rotate the matrix back from the standard orientation to the user orientation
			auto invR_ut = linear_algebra::mul_mat_mat(transpose(R), ut);
			// provide a message to the user if the basis vectors have changed outside of approximate tolerance
			if (!approx_float::matrix(3u, use_real ? _real_vectors.data() : _reciprocal_vectors.data(), invR_ut.data())){
				std::ostringstream msg;
				msg << "Basis vectors changed by [";
				for (size_t i=0; i<3u; ++i){
					msg << "(";
					for (size_t j=0; j<3u; ++j) msg << " " << invR_ut[i + j*3u] -  _reciprocal_vectors[i + j*3u];
					msg << " ), ";
				}
				msg << "] " << (use_real ? u8"Å" : u8"Å⁻¹");
				info_update(msg.str());
			}
			// And actually update the stored parameters to the new reciprocal and real vectors
			_reciprocal_vectors = invR_ut;
			_real_vectors = transpose(linear_algebra::mat_inverse(_reciprocal_vectors));
      for (auto & x: _real_vectors) x *= math::two_pi;
    }
  }
//  void snap_basis_vectors_to_symmetry(){
//    using approx_float::scalar;
//    using linear_algebra::dot;
//    using linear_algebra::norm, linear_algebra::mul_mat_mat, linear_algebra::mul_mat_vec;
//    matrix_t vectors{{0,0,0, 0,0,0, 0,0,0}}, deviations{{0,0,0, 0,0,0, 0,0,0}};
//    size_t counts[3]{0,0,0};
//
//    matrix_t _bv{_real_vectors}, _inv_bv{linear_algebra::mat_inverse(_real_vectors)};
//
//    std::vector<std::array<int,3>> e_i_list {{1, 0, 0}, {0, 1, 0}, {0, 0, 1},
//                                             {-1, 0, 0}, {0, -1, 0}, {0, 0, -1},
//                                             {1, 1, 0}, {1, 0, 1}, {0, 1, 1},
//                                             {-1, 1, 0}, {-1, 0, 1}, {0, -1, 1},
//                                             {1, -1, 0}, {1, 0, -1}, {0, 1, -1},
//                                             {-1, -1, 0}, {-1, 0, -1}, {0, -1, -1},
//                                             {1, 1, 1}, {-1, 1, 1}, {1, -1, 1},
//                                             {1, 1, -1}, {-1, -1, 1}, {-1, 1, -1},
//                                             {1, -1, -1}, {-1, -1, -1}};
//    auto vectors_match = [](const auto & a, const auto & b){
//      if (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) return true;
//      return false;
//    };
//    std::vector<vector_t> vs;
//    for (size_t i=0; i<3u; ++i) vs.push_back(vector(LengthUnit::angstrom, i));
//
//    std::array<std::vector<vector_t>, 3> all_e_j;
//
//    for (const auto & e_i_lu: e_i_list){
//      vector_t e_i_c{0,0,0};
//      for (size_t i=0; i<3u; ++i) if (e_i_lu[i]) {
//        for (size_t j=0; j<3u; ++j) e_i_c[j] += e_i_lu[i] * vs[i][j];
//      }
//      // loop over all output real space basis vectors
//      for (size_t e_j =0; e_j <3u; ++e_j){
//        std::array<int,3> e_j_lu{{0,0,0}}; e_j_lu[e_j] = 1;
//        if (!vectors_match(e_i_lu, e_j_lu)) {
//          info_update("Operations which map ", e_i_lu," to ", e_j_lu);
//          // and for all point group operations that link e_i to e_j
//          for (const auto &r : _point.mat_vec_linking_operations(e_i_lu, e_j_lu)) {
//            // The rotation matrix in cartesian coordinates
//            auto r_c = mul_mat_mat(_bv, mul_mat_mat(r, _inv_bv));
//            // now rotate the cartesian basis vector to point along e_j
//            auto e_j_c = mul_mat_vec(r_c, e_i_c);
//            //          info_update(e_j_c);
//            // add the result to the output basis vector at 'e_j'
//            for (size_t i = 0; i < 3u; ++i) vectors[e_j + 3u * i] += e_j_c[i];
//            // add the difference to the output deviations at 'e_j'
//            for (size_t i = 0; i < 3u; ++i) {
//              auto x = e_j_c[i] - _bv[e_j + 3u * i];
//              deviations[e_j + 3u * i] += x * x;
//            }
//            ++counts[e_j];
//            all_e_j[e_j].push_back(e_j_c);
//          }
//        }
//      }
//    }
//    // Show all found 'matching' basis vectors:
//    for (size_t i=0; i<3; ++i){
//      std::cout << i << ". basis vectors" << std::endl;
//      for (const auto & x: all_e_j[i]){
//        for (const auto & z: x) std::cout << std::fixed << std::setprecision(16) << std::setw(20) << z << " ";
//        std::cout << std::endl;
//      }
//      std::cout << std::endl;
//    }
//
//    // complete finding the average output basis vectors
//    for (size_t z=0; z<9u; ++z) vectors[z] /= static_cast<T>(counts[z % 3u]);
//    for (size_t z=0; z<9u; ++z) deviations[z] /= static_cast<T>(counts[z % 3u]);
//    auto greatest_deviation = std::sqrt(*std::max_element(deviations.begin(), deviations.end()));
//    for (size_t z=0; z<9u; ++z) deviations[z] = std::abs(vectors[z] - _real_vectors[z]);
//    auto greatest_change = *std::max_element(deviations.begin(), deviations.end());
//    info_update("Deviations ", greatest_deviation, " changes ", greatest_change);
//
//    // TODO actually check if the basis vectors changed?
//    _real_vectors = vectors;
//    _reciprocal_vectors = transpose(linear_algebra::mat_inverse(_real_vectors));
//    for (auto & x: _reciprocal_vectors) x *= math::two_pi;
//  }
  void set_vectors(const vector_t& v, const vector_t& a, LengthUnit lu, AngleUnit angle_unit, bool snap_to_symmetry=false){
    auto [cos, sin] = inter_facial_angles_to_cosines_sines(a, angle_unit);
    // arrays are small; so copying shouldn't hurt
    vector_t dv{v}, dcos{cos}, dsin{sin}, rv{v}, rcos{cos}, rsin{sin};
    switch (lu) {
    case LengthUnit::angstrom:
      std::tie(rv, rcos, rsin) = dual_lattice_parameters(v, cos, sin); break;
    case LengthUnit::inverse_angstrom:
      std::tie(dv, dcos, dsin) = dual_lattice_parameters(v, cos, sin); break;
    default:
      throw std::logic_error("The length unit must be angstrom or inverse angstrom!");
    }
    if (snap_to_symmetry && snap_parameters_to_symmetry(false, rv, rcos, rsin)){
      std::tie(dv, dcos, dsin) = dual_lattice_parameters(rv, rcos, rsin);
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
   * @param AorB Either the real space or reciprocal space basis vector matrix,
   *          which has the three basis vectors as its *columns*.
   * @param lu An enumerated length unit to indicate if the passed matrix
   *           represents the real (LengthUnit::angstrom) or reciprocal
   *           (LengthUnit::inverse_angstrom) basis vectors.
   */
  void set_vectors(const matrix_t& AorB, LengthUnit lu, bool snap_to_symmetry=false){
    auto BorA_transposed = linear_algebra::mat_inverse(AorB);
    for (auto & x: BorA_transposed) x *= math::two_pi;
    if (LengthUnit::angstrom == lu){
      _real_vectors = AorB;
      _reciprocal_vectors = transpose(BorA_transposed);
    } else if (LengthUnit::inverse_angstrom == lu) {
      _real_vectors = transpose(BorA_transposed);
      _reciprocal_vectors = AorB;
    }
    if (snap_to_symmetry){
      snap_basis_vectors_to_symmetry();
    }
  }
  void set_metrics() {
    // the two metrics are mutually inverse with a factor of 4 * pi^2
    // and of course A and B are *both* made of column vectors, so their
    // metrics are AᵀA and BᵀB *not AAᵀ and BBᵀ!
    _real_metric = metric_from_column_vectors(_real_vectors);
    _reciprocal_metric = metric_from_column_vectors(_reciprocal_vectors);
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
    if (Bravais::P != _bravais) rep += this->primitive().to_verbose_string();
    return rep;
  }

  [[nodiscard]] matrix_t real_basis_vectors() const {return _real_vectors;}
  [[nodiscard]] matrix_t reciprocal_basis_vectors() const {return _reciprocal_vectors;}

  [[nodiscard]] Bravais bravais() const {return _bravais;}
  [[nodiscard]] Symmetry spacegroup_symmetry() const {return _space;}
  [[nodiscard]] PointSymmetry pointgroup_symmetry() const {return _point;}
  [[nodiscard]] Basis basis() const {return _basis;}

  [[nodiscard]] matrix_t vectors(const LengthUnit lu) const {
    switch (lu){
      case LengthUnit::angstrom: return _real_vectors;
      case LengthUnit::inverse_angstrom: return _reciprocal_vectors;
      default:
        throw std::logic_error("Not implemented length unit");
    }
  }

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
    assert(I(0) <= i && i <= I(3));
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
  [[nodiscard]] bool is_same(const Impl<R>& o) const {
    // two lattices are *the* same if their basis vectors are the same
    auto bv = o.real_basis_vectors();
    auto rbv = o.reciprocal_basis_vectors();
    return approx_float::equal(bv, _real_vectors) && approx_float::equal(rbv, _reciprocal_vectors);
  }
  template<class R>
  [[nodiscard]] bool is_equivalent(const Impl<R>& o) const {
    return is_same(o);
  }
  template<class R>
  [[nodiscard]] int is_permuted(const Impl<R>& o) const {
    return is_same(o);
  }
  template<class R>
  bool operator == (const Impl<R>& l) const {return is_same(l);}
  template<class R>
  bool operator != (const Impl<R>& l) const {return !is_same(l);}

  Symmetry spacegroup_symmetry(const Symmetry& gens){
    _space = gens.generate();
    _bravais = _space.getcentring();
    _point = PointSymmetry(get_unique_rotations(_space.getallr(), 0));
    return _space;
  }
  [[nodiscard]] bool has_space_inversion() const {return _point.has_space_inversion();}
  [[nodiscard]] bool is_triclinic() const {return _point.higher(1).size() == 0u;}

  Basis basis(const std::vector<vector_t>& pos, const std::vector<ind_t>& typ, bool snap_to_symmetry=true){
    _basis = Basis(pos, typ);
    snap_basis_to_symmetry(snap_to_symmetry);
    return _basis;
  }
  Basis basis(const Basis& b, bool snap_to_symmetry=true){
    _basis = b;
    snap_basis_to_symmetry(snap_to_symmetry);
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
  [[nodiscard]] matrix_t from_xyz(LengthUnit lu) const {
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
  static std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, Impl<T>>
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

  Impl<T> primitive() const {
    PrimitiveTransform P(_bravais);
    if (P.does_anything()){
      // The transformation matrix P gives us the primitive basis column-vector
      // matrix Aₚ from the standard basis column-vector matrix Aₛ by
      // Aₚ = Aₛ P.
      // The PrimitiveTransform object contains 6*P (so that it's integer valued)
//      auto P6 = P.get_6P();
//      std::array<T, 9> tP;
//      std::transform(P6.begin(), P6.end(), tP.begin(), [](const auto & x){return static_cast<T>(x)/T(6);});
      auto pv = linear_algebra::mul_mat_mat(_real_vectors, P.get_6P());
      for (auto & x: pv) x /= T(6);
      // calculating the new metric tensor could be done from _real_vectors, but
      // we might as well use the new values in pv
      auto pvm = metric_from_column_vectors(pv);

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
      auto prm = metric_from_column_vectors(pr);

      // strictly we should change the spacegroup and pointgroup symmetry
      // representations, plus the atom basis representation... but that seems
      // overkill for now. FIXME think about this more.
      return {pv, pr, pvm, prm, Bravais::P, _space, _point, _basis};
    }
    return *this;
  }

}; // end of Impl

#define LATTICE_FORWARD_METHOD(RETURNED, NAMED) template<class... A> [[nodiscard]] RETURNED NAMED(A... a) const {return ptr->NAMED(a...);}

template<class T>
class Lattice{
public:
  using matrix_t = typename Impl<T>::matrix_t;
  using vector_t [[maybe_unused]] = typename Impl<T>::vector_t;
private:
  std::shared_ptr<Impl<T>> ptr;
protected:
  std::shared_ptr<Impl<T>> pointer() const {return ptr;}
public:
  Lattice() = default;

  // Define copy constructors to avoid forwarding to Impl
  // Which is only a problem for musl libc, somehow.
  Lattice(Lattice<T>& lat): ptr(lat.pointer()) {}
  Lattice(const Lattice<T>& lat): ptr(lat.pointer()) {}
  Lattice(Lattice<T>&& lat) noexcept : ptr(std::move(lat.ptr)) {}
  // ^^^ copy constructors require we define assignment too:
  Lattice<T>& operator =(const Lattice<T>& other){
    ptr = other.pointer();
    return *this;
  }

  explicit Lattice(std::shared_ptr<Impl<T>> ptr_): ptr(ptr_) {}
  template<class... Args>
  explicit Lattice(Args&&... args) { ptr = std::make_shared<Impl<T>>(std::forward<Args>(args)...);}
  template<class... Args>
  explicit Lattice(Args&... args) {ptr = std::make_shared<Impl<T>>(args...);}

  template<class... Args>
  Lattice(LengthUnit lu, vector_t&& v, vector_t&& a, Args... args): ptr(std::make_shared<Impl<T>>(lu, std::move(v), std::move(a), args...)) {}

  LATTICE_FORWARD_METHOD(std::string, to_string)
  LATTICE_FORWARD_METHOD(std::string, to_verbose_string)
  LATTICE_FORWARD_METHOD(matrix_t, real_basis_vectors)
  LATTICE_FORWARD_METHOD(matrix_t, reciprocal_basis_vectors)
  LATTICE_FORWARD_METHOD(Bravais, bravais)
  LATTICE_FORWARD_METHOD(Symmetry, spacegroup_symmetry)
  LATTICE_FORWARD_METHOD(PointSymmetry, pointgroup_symmetry)
  LATTICE_FORWARD_METHOD(Basis, basis)
  LATTICE_FORWARD_METHOD(matrix_t, vectors)
  LATTICE_FORWARD_METHOD(matrix_t, metric)
  LATTICE_FORWARD_METHOD(vector_t, vector)
  LATTICE_FORWARD_METHOD(T, length)
  LATTICE_FORWARD_METHOD(T, angle)
  LATTICE_FORWARD_METHOD(vector_t, lengths)
  LATTICE_FORWARD_METHOD(vector_t, angles)
  LATTICE_FORWARD_METHOD(T, volume)
  template<class R> bool is_same(const Lattice<R>& l) const {
    if (ptr == l.pointer()) return true;
    return ptr->is_same(*l.pointer());
  }
  template<class R> bool is_equivalent(const Lattice<R>& l) const {
    if (ptr == l.pointer()) return true;
    return ptr->is_equivalent(*l.pointer());
  }
  template<class R> bool is_permuted(const Lattice<R>& l) const {
    if (ptr == l.pointer()) return true;
    return ptr->is_permuted(*l.pointer());
  }
  template<class R> bool operator == (const Lattice<R>& l) const {
    if (ptr == l.pointer()) return true;
    return (*ptr) == (*l.pointer());
  }
  template<class R> bool operator != (const Lattice<R>& l) const {
    if (ptr == l.pointer()) return false;
    return (*ptr) != (*l.pointer());
  }
  LATTICE_FORWARD_METHOD(bool, has_space_inversion)
  LATTICE_FORWARD_METHOD(bool, is_triclinic)
  LATTICE_FORWARD_METHOD(matrix_t, to_xyz)
  LATTICE_FORWARD_METHOD(matrix_t, from_xyz)
#ifdef USE_HIGHFIVE
  LATTICE_FORWARD_METHOD(bool, to_hdf)
  template<class... A> static Lattice<T> from_hdf(A... a){
    auto li = Impl<T>::from_hdf(a...);
    return Lattice<T>(std::make_shared<Impl<T>>(li));
  }
#endif
  Lattice<T> primitive() const {
    if (PrimitiveTransform(ptr->bravais()).does_anything()) {
      return Lattice<T>(std::make_shared<Impl<T>>(ptr->primitive()));
    }
    return Lattice<T>(ptr);
  }
};


template<class T, class... Args>
Lattice<T> Direct(const std::array<T,3> & v, const std::array<T,3> & a, Args ... args){
  return Lattice<T>(LengthUnit::angstrom, v, a, args...);
}
template<class T, class... Args>
Lattice<T> Reciprocal(const std::array<T,3> & v, const std::array<T,3> & a, Args ... args){
  return Lattice<T>(LengthUnit::inverse_angstrom, v, a, args...);
}
//template<class T, class... Args>
//Impl<T> Direct(const std::initializer_list<T> & v, const std::initializer_list<T> & a, Args ... args){
//  return {v, a, args..., LengthUnit::angstrom};
//}
//template<class T, class... Args>
//Impl<T> Reciprocal(const std::initializer_list<T> & v, const std::initializer_list<T> & a, Args ... args){
//  return {v, a, args..., LengthUnit::inverse_angstrom};
//}


template<class T, class... Args>
Lattice<T> Direct(const std::array<T,9> & m, Args ... args){
  return Lattice<T>(LengthUnit::angstrom, m, args...);
}
template<class T, class... Args>
Lattice<T> Reciprocal(const std::array<T,9> & m, Args ... args){
  return Lattice<T>(LengthUnit::inverse_angstrom, m, args...);
}

}

#endif
