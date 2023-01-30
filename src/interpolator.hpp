/* This file is part of brille.

Copyright © 2019,2020 Greg Tucker <greg.tucker@stfc.ac.uk>

brille is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

brille is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with brille. If not, see <https://www.gnu.org/licenses/>.            */
#ifndef BRILLE_INTERPOLATOR_HPP_
#define BRILLE_INTERPOLATOR_HPP_
/*! \file
    \author Greg Tucker
    \brief Class for linear interpolation of arbitrary data
*/
// #include <vector>
// #include <array>
#include <utility>
#include <mutex>
// #include <cassert>
// #include <functional>
#include <omp.h>
#include "phonon.hpp"
#include "permutation.hpp"
#include "permutation_table.hpp"
#include "rotates.hpp"
namespace brille {

/*! \brief A function to calculate a scalar property of two arrays

The `brille::Interpolator` class holds a constant type of information at each
of its points. This information is held in a contiguous array and might
represent some number of sub-arrays each with the same character. In order to
facilitate characterising similarities between sub-arrays at one or different
points it is useful to define a general function signature, the `CostFunction`.

A `CostFunction` takes the number of array elements which should be compared and
two constant pointers to the first element of each (sub)array. Ideally the
function performs some calculation based on the contents of the two arrays and
returns a single double scalar.

Although not a requirement, it is anticipated that the returned value of a
`CostFunction` will be zero if the two arrays are identical and positive if the
arrays differ.
*/
template<class T>
using CostFunction = std::function<double(brille::ind_t, const T*, const T*)>;
//! A template helper to differentiate complex valued containers
template<class T> struct is_complex {enum{value = false};};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T> struct is_complex<std::complex<T>> {enum {value=true};};
#endif

/*! \brief A class to hold and linearly interpolate arbitrary data

The arbitrary data is held in a contiguous two-dimensional Array2 object. Its
first dimension is typically indexed by the point-index of a connected graph
and its second dimension contains a single format of data for all points.
Often one `Interpolator` contains all of the eigenvalues or eigenvectors of a
system for the indexed point, and two are used in concert by a
`DualInterpolator`.

The format of the data for each point is defined by the size of the array and
a user specified sub-array character listing the number of scalar, vector, and
matrix elements. The total number of sub-arrays is then the size of the array
divided by the sum of the character.
If there are \f$i\f$ scalars, \f$j\f$ vector elements, and \f$k\f$ matrix
elements, then each sub-array character is
\f[ c = \left[ s_0, s_1, \ldots, s_i, v_0, v_1, \ldots, v_j, m_0, m_1, \ldots, m_k \right]\f]
and the format of the array for each point is a tiling of the sub-arrays
\f[\left\{c_0, c_1, \ldots, c_b\right\}\f]
for \f$b\f$ total sub-arrays.

The Array2 can often be interpreted as a higher-dimensional Array which has been
flattened beyond the second dimension. The Interpolator object contains a
property to hold this higher-dimensional shape information.
*/
template<class T>
class Interpolator{
public:
  template<class Z> using element_t =std::array<Z,3>;
  using costfun_t = CostFunction<T>;
  using shape_t = std::vector<ind_t>;
  using interpolated_out_t = brille::Array<T>;
  using interpolated_in_t = brille::Array2<T>;
protected:
  bArray<T> data_;      //!< The stored Array2 of points indexed like the holding-Object's vertices
  shape_t shape_;       //!< The shape of the input Array (or Array2)
  element_t<ind_t> _elements; //!< The number of each element type per point and per mode
  RotatesLike rotlike_;   //!< How the elements of `data_` rotate
  LengthUnit lenunit_;   //!< The units of `data_`
  element_t<double> _costmult; //!< The relative (multiplicative) cost for differences in each element type
  element_t<ind_t> _funtype;
  costfun_t _scalarfun; //!< A function to calculate differences between the scalars at two stored points
  costfun_t _vectorfun; //!< A function to calculate differences between the vectors at two stored points
  //costfun_t _matrixfun; //!< A function to calculate the differences between matrices at two stored points
public:
  bool operator!=(const Interpolator<T>& other) const {
    if (data_ != other.data_) return true;
    if (shape_ != other.shape_) return true;
    if (_elements != other._elements) return true;
    if (rotlike_ != other.rotlike_) return true;
    if (lenunit_ != other.lenunit_) return true;
    if (_costmult != other._costmult) return true;
    if (_funtype != other._funtype) return true;
    return false;
  }
  /*! \brief Constructor without data and with optional cost function types

  \param scf_type The scalar cost function *type*
  \param vcf_type The vector cost function *type*
  \see set_cost_info
  */
  explicit Interpolator(int scf_type=0, int vcf_type=0)
  : data_(0,0), _elements({{0,0,0}}), rotlike_{RotatesLike::vector}, lenunit_{LengthUnit::real_lattice}, _costmult({{1,1,1}})
  {
    this->set_cost_info(scf_type, vcf_type);
  }
  /*! \brief Constructor without data specifying cost functions

  \param scf The scalar cost function
  \param vcf The vector cost function
  \see CostFunction
  */
  Interpolator(costfun_t scf, costfun_t vcf)
  : data_(0,0), _elements({{0,0,0}}), rotlike_{RotatesLike::vector}, lenunit_{LengthUnit::real_lattice},
    _costmult({{1,1,1}}), _scalarfun(scf), _vectorfun(vcf)
  {}
  /*! \brief Partial constructor with default cost functions

  \param d The data to be interpolated
  \param sh The \f$N\f$-dimensional shape of the data array
  \param el The sub-array character specifier
  \param rl How any vectors or tensors transform under application of a symmetry operation
  \param lu Units of vectors/tensors
  \see RotatesLike, LengthUnit
  */
  Interpolator(bArray<T>& d, shape_t sh, element_t<ind_t> el, RotatesLike rl, LengthUnit lu)
  : data_(d), shape_(sh), _elements(el), rotlike_{rl}, lenunit_{lu}, _costmult({{1,1,1}})
  {
    this->set_cost_info(0,0);
    this->check_elements();
  }
  /*! \brief Constructor specifying cost function types and weights

  \param d The data to be interpolated
  \param sh The \f$N\f$-dimensional shape of the data array
  \param el The sub-array character specifier
  \param rl How any vectors or tensors transform under application of a symmetry operation
  \param lu Units of vectors/tensors
  \param csf The scalar cost function type
  \param cvf The vector cost function type
  \param wg The relative scalar, vector, and matrix cost scaling values
  \see RotatesLike, LengthUnit, set_cost_info
  */
  Interpolator(bArray<T>& d, shape_t sh, element_t<ind_t> el, RotatesLike rl, LengthUnit lu, int csf, int cvf, element_t<double> wg)
  : data_(d), shape_(sh), _elements(el), rotlike_{rl}, lenunit_{lu}, _costmult(wg)
  {
    this->set_cost_info(csf, cvf);
    this->check_elements();
  }
  /*! \brief Partial constructor with default cost functions

  \param d The \f$N\f$-dimensional data to be interpolated, will be flattend to 2-D
  \param el The sub-array character specifier
  \param rl How any vectors or tensors transform under application of a symmetry operation
  \param lu Units of vectors/tensors
  \see RotatesLike, LengthUnit
  */
  // use the Array2<T>(const Array<T>&) constructor
  Interpolator(brille::Array<T>& d, element_t<ind_t> el, RotatesLike rl, LengthUnit lu)
  : data_(d), shape_(d.shape()), _elements(el), rotlike_{rl}, lenunit_{lu}, _costmult({{1,1,1}})
  {
    this->set_cost_info(0,0);
    this->check_elements();
  }
  /*! \brief Constructor specifying cost function types and weights

  \param d The \f$N\f$-dimensional data to be interpolated, will be flattend to 2-D
  \param el The sub-array character specifier
  \param rl How any vectors or tensors transform under application of a symmetry operation
  \param lu Units of vectors/tensors
  \param csf The scalar cost function type
  \param cvf The vector cost function type
  \param wg The relative scalar, vector, and matrix cost scaling values
  \see RotatesLike, LengthUnit, set_cost_info
  */
  // use the Array2<T>(const Array<T>&) constructor
  Interpolator(brille::Array<T>& d, element_t<ind_t> el, RotatesLike rl, LengthUnit lu, int csf, int cvf, element_t<double> wg)
  : data_(d), shape_(d.shape()), _elements(el), rotlike_{rl}, lenunit_{lu}, _costmult(wg)
  {
    this->set_cost_info(csf, cvf);
    this->check_elements();
  }
  /*! \brief Setup a fake dataset for interpolation

  \param sz The number of points covered by the interpolation object
  \param br The number of sub-arrays per point
  \note This method allows a brille::DualInterpolator to be constructed with
        empty Interpolator objects and then have each replaced with the intended
        interpolation data, one at a time. When the first Interpolator data has
        been replaced the second must contain data for the same number of points
        and the same number of sub-arrays, so this method produces a suitable
        fake dataset for the second Interpolator.
  */
  void setup_fake(const ind_t sz, const ind_t br){
    data_ = bArray<T>(sz, br);
    shape_ = {sz, br};
    _elements = {1u,0u,0u};
  }
  /*! \brief Select scalar and vector cost functions by their type indices

  \param scf The scalar cost function type
  \param vcf The vector cost function type

  A limited set of cost functions are allowed to be set by type, if the function
  you would like to use is not available you can provide a `CostFunction` object
  instead when creating the `Interpolator` object.


  | scalar index | cost function description   | implemented using          |
  | :----------: | :-------------------------- | :------------------------- |
  | any value    | magnitude of the difference | `brille::utils::magnitude` |

  The scalar magnitude function is \f$\sum_k\left|s_{ik} - s_{jk}\right|\f$.

  | vector index | cost function description          | implemented using                |
  | :----------: | :--------------------------------- | :------------------------------- |
  | 1            | magnitude of the difference vector | `brille::utils::vector_distance` |
  | 2            | one minus the vector product       | `brille::utils::vector_product`  |
  | 3            | Euclidean angle between vectors    | `brille::utils::vector_angle`    |
  | 4            | Hermitian angle between vectors    | `brille::utils::hermitian_angle` |
  | other values | sine squared of Hermitian angle    | ^                                |

  The vector difference function is \f$ \left|\mathbf{v}_i - \mathbf{v}_j\right| \f$.
  The vector product function is \f$ 1 - \mathbf{v}^*_i \cdot \mathbf{v}_j \f$.
  \see brille::CostFunction
  */
  void set_cost_info(const int scf, const int vcf){
    _funtype = element_t<ind_t>({0,0,0});
    switch (scf){
      case 0:
      default:
      this->_scalarfun = [](ind_t n, const T* i, const T* j){
        double s{0};
        for (ind_t z=0; z<n; ++z) s += brille::utils::magnitude(i[z]-j[z]);
        return s;
      };
    }
    switch (vcf){
      case 1:
      debug_update("selecting brille::utils::vector_distance");
      this->_vectorfun = [](ind_t n, const T* i, const T* j){
        return brille::utils::vector_distance(n, i, j);
      };
      _funtype[0] = 1;
      break;
      case 2:
      debug_update("selecting 1-brille::utils::vector_product");
      this->_vectorfun = [](ind_t n, const T* i, const T* j){
        return 1-brille::utils::vector_product(n, i, j);
      };
      _funtype[0] = 2;
      break;
      case 3:
      debug_update("selecting brille::utils::vector_angle");
      this->_vectorfun = [](ind_t n, const T* i, const T* j){
        return brille::utils::vector_angle(n, i, j);
      };
      _funtype[0] = 3;
      break;
      case 4:
      debug_update("selecting brille::utils::hermitian_angle");
      this->_vectorfun = [](ind_t n, const T* i, const T* j){
         return brille::utils::hermitian_angle(n,i,j);
      };
      _funtype[0] = 4;
      break;
      default:
      debug_update("selecting sin**2(brille::utils::hermitian_angle)");
      // this->_vectorfun = [](ind_t n, T* i, T* j){return std::abs(std::sin(brille::utils::hermitian_angle(n, i, j)));};
      this->_vectorfun = [](ind_t n, const T* i, const T* j){
        auto sin_theta_H = std::sin(brille::utils::hermitian_angle(n, i, j));
        return sin_theta_H*sin_theta_H;
      };
    }
  }
  /*! \brief Set the cost functions by type and the relative cost weights

  \param scf The scalar cost function type
  \param vcf The vector cost function type
  \param elcost The relative scalar, vector, and matrix cost scaling values

  */
  void set_cost_info(const int scf, const int vcf, const element_t<double>& elcost){
    _costmult = elcost;
    this->set_cost_info(scf, vcf);
  }
  //! Return the number of points the Interpolator covers
  ind_t size(void) const {return data_.size(0);}
  //! Return the number of sub-arrays per point
  ind_t branches(void) const {
    /*
    The data Array is allowed to be anywhere from 1 to 5 dimensional.
    Possible shapes are:
      (P,)            - one branch with one scalar per point
      (P, X)          - one branch with X elements per point
      (P, B, Y)       - B branches with Y elements per branch per point
      (P, B, V, 3)    - B branches with V 3-vectors per branch per point
      (P, B, M, 3, 3) - B branches with M (3,3) matrices per branch per point
    In the two-dimensional case there is a potential ambiguity in that X counts
    both the number of branches and the number of elements.
      X = B*Y
    where Y is the sum of scalar, vector, and matrix elements per branch
      Y = S + 3*V + 9*M
    */
    ind_t nd = static_cast<ind_t>(shape_.size());
    ind_t b = nd>1 ? shape_[1] : 1u;
    if (2 == nd){
      ind_t y = std::accumulate(_elements.begin(), _elements.end(), ind_t(0));
      // zero-elements is the special (initial) case → 1 scalar per branch
      if (y > 0) b /= y;
    }
    return b;
  }
  //! \brief Determine if the sub-arrays are only vectors or matrices
  bool only_vector_or_matrix(void) const {
    // if V or M can not be deduced directly from the shape of data_
    // then this Interpolator represents mixed data
    // purely-scalar data is also classed as 'mixed' for our purposes
    ind_t nd = shape_.size();
    if (5u == nd && 3u == shape_[4] && 3u == shape_[3]) return true;
    if (4u == nd && 3u == shape_[2]) return true;
    if (nd < 4u) return false; // (P,), (P,X), (P,B,Y)
    std::string msg = "Interpolator can not handle a {";
    for (auto x: shape_) msg += " " + std::to_string(x);
    msg += " } data array";
    throw std::runtime_error(msg);
  }
  //! Obtain a reference to the interpolation data Array2
  const bArray<T>& data(void) const {return data_;}
  //! Return the \f$N\f$-dimensional interpolation data Array shape
  shape_t shape(void) const {return shape_;};
  //! Return the \f$N\f$-dimensional interpolation data Array
  brille::Array<T> array(void) const {return brille::Array<T>(data_,shape_);}
  //! Return the sub-array character specifier
  element_t<ind_t> elements(void) const {return _elements;}
  /*! \brief Perform linear interpolation between a set of stored points with
             provided weights

  \note This gateway exists in case it becomes necessary or advantageous to
        handle the 'mixed' character sub-array case separately from the pure
        vector or pure matrix cases.
  \see Interpolator::interpolate_at_mix
  */
  template<class... Args>
  void interpolate_at(Args... args) const {
      this->interpolate_at_mix(args...);
  }
  /*! \brief Apply a symmetry operation to each interpolation result

  \param x    The pre-determined interpolation result
  \param q    The reciprocal lattice point for each interpolation result,
              only used if this Interpolator holds phonon eigenvectors
  \param rt   A helper rotation table, e.g., `brille::GammaTable`,
              only used if this Interpolator holds phonon eigenvectors
  \param ps   All pointgroup symmetry operations for the lattice
  \param r    An index into `ps` for the operation required for each
              interpolation result
  \param invr An index into `ps` for the inverse of the operation required
              for each interpolation result
  \param nth  The number of OpenMP threads that this method should utilize

  This method is a gateway switching on how the stored data transforms under
  application of a symmetry operation.
  \see Interpolator::rip_real, Interpolator::rip_axial, Interpolator::rip_recip
  \see Interpolator::rip_gamma
  */
  template<class R, class RotT>
  bool rotate_in_place(bArray<T>& x,
                       const lattice::LVec<R>& q,
                       const RotT& rt,
                       const PointSymmetry& ps,
                       const std::vector<size_t>& r,
                       const std::vector<size_t>& invr,
                       const int nth=0) const {
    switch (lenunit_) {
      case LengthUnit::real_lattice: {
        switch (rotlike_) {
            case RotatesLike::vector:       return this->rip_real(x,ps,r,invr,nth);
            case RotatesLike::pseudovector: return this->rip_axial(x,ps,r,invr,nth);
            case RotatesLike::Gamma:        return this->rip_gamma(x,q,rt,ps.getall(),r,invr,nth);
            default: throw std::runtime_error("Impossible RotatesLike value!");
        }
      }
      case LengthUnit::reciprocal_lattice: {
        switch (rotlike_) {
            case RotatesLike::vector: return this->rip_recip(x,ps,r,invr,nth);
            default: throw std::runtime_error("LengthUnit, RotatesLike combination not implemented");
        }
      }
      case LengthUnit::angstrom: {
        switch (rotlike_) {
            case RotatesLike::Gamma: {
              // Convert rotation matrix from fractional to Cartesian
              std::vector<std::array<double,9>> rot_cart(ps.size());
              std::array<double,9> tdbl0;
              for (size_t j=0; j<ps.size(); j++){
                // tdbl0 = lattice*rot
                brille::utils::mul_mat_mat(tdbl0.data(), 3u, rt.lattice().to_xyz(LengthUnit::angstrom).data(), ps.data(j));
                // rot_cart = tdbl0*inv(lattice) = lattice*rot*inv(lattice)
                brille::utils::mul_mat_mat(rot_cart[j].data(), 3u, tdbl0.data(), rt.lattice().from_xyz(LengthUnit::angstrom).data());
              }
              return this->rip_gamma(x,q,rt,rot_cart,r,invr,nth);
            }
            default: throw std::runtime_error("LengthUnit, RotatesLike combination not implemented");
        }
      }
      default: throw std::runtime_error("LengthUnit, RotatesLike combination not implemented");
    }
  }
  //! Return the transform type of the stored data under symmetry operation appliation.
  RotatesLike rotateslike() const { return rotlike_; }
  //! Set the transform type of the stored data under symmetry operation appliation.
  RotatesLike rotateslike(const RotatesLike a) {
    rotlike_ = a;
    return rotlike_;
  }
  /*! \brief Replace the data within this object.

  \param nd The new data to be interpolated
  \param sh The \f$N\f$-dimensional shape of the data to be interpolated
  \param ne The new sub-array character specification
  \param rl How the new data transforms under application of a symmetry operation
  \param lu Units of new data
  */
  template<class I>
  void replace_data(
      const bArray<T>& nd,
      const shape_t sh,
      const std::array<I,3>& ne,
      const RotatesLike rl = RotatesLike::vector,
      const LengthUnit lu = LengthUnit::real_lattice)
  {
    data_ = nd;
    shape_ = sh;
    rotlike_ = rl;
    lenunit_ = lu;
    // convert the elements datatype as necessary
    if (ne[1]%3)
      throw std::logic_error("Vectors must have 3N elements per branch");
    if (ne[2]%9)
      throw std::logic_error("Matrices must have 9N elements per branch");
    for (size_t i=0; i<3u; ++i) _elements[i] = static_cast<ind_t>(ne[i]);
    this->check_elements();
  }
  /*! \brief Replace the data within this object.

  \param nd The new \f$N\f$-dimensional data to be interpolated
  \param ne The new sub-array character specification
  \param rl How the new data transforms under application of a symmetry operation
  \param lu Units of new data

  \note The provided data will only be used in-place if it is a contiguous
        row-ordered array. In all other cases a copy will be made before
        flattening into two-dimensions.
  */
  template<class I>
  void replace_data(
      const brille::Array<T>& nd,
      const std::array<I,3>& ne,
      const RotatesLike rl = RotatesLike::vector,
      const LengthUnit lu = LengthUnit::real_lattice)
  {
    data_ = bArray<T>(nd);
    shape_ = nd.shape();
    rotlike_ = rl;
    lenunit_ = lu;
    // convert the elements datatype as necessary
    if (ne[1]%3)
      throw std::logic_error("Vectors must have 3N elements per branch");
    if (ne[2]%9)
      throw std::logic_error("Matrices must have 9N elements per branch");
    for (size_t i=0; i<3u; ++i) _elements[i] = static_cast<ind_t>(ne[i]);
    this->check_elements();
  }
  // Replace the data in this object without specifying the data shape or its elements
  // this variant is necessary since the template specialization above can not have a default value for the elements
  template<template<class> class A>
  void replace_data(const A<T>& nd){
    return this->replace_data(nd, element_t<ind_t>({{0,0,0}}));
  }
  //! Determine the total size of the stored sub-arrays from their character specification
  ind_t branch_span() const { return this->branch_span(_elements);}
  //! Produce an informational string describing the object
  std::string to_string() const {
    std::string str= "{ ";
    for (auto s: shape_) str += std::to_string(s) + " ";
    str += "} data";
    auto b = this->branches();
    if (b){
      str += " with " + std::to_string(b) + " mode";
      if (b>1) str += "s";
    }
    auto n = std::count_if(_elements.begin(), _elements.end(), [](ind_t a){return a>0;});
    if (n){
      str += " of ";
      std::array<std::string,3> types{"scalar", "vector", "matrix"};
      for (size_t i=0; i<3u; ++i) if (_elements[i]) {
        str += std::to_string(_elements[i]) + " " + types[i];
        if (--n>1) str += ", ";
        if (1==n) str += " and ";
      }
      str += " element";
      if (this->branch_span()>1) str += "s";
    }
    return str;
  }

  template<class S>
  void add_cost(const ind_t, const ind_t, std::vector<S>&, bool) const;

  template<typename I>
  bool any_equal_modes(const I idx) const {
    return this->any_equal_modes_(static_cast<ind_t>(idx), this->branches(), this->branch_span());
  }
  /*! \brief Calculate the memory requirements per interpolation result

  Linear interpolation of the stored data requires an output array with the same
  first dimension size as the number of points at which interpolation is to be
  performed. This method can be used by driving programs to estimate how much
  memory such an array will occupy to avoid causing out-of-memory errors at
  runtime.
  */
  size_t bytes_per_point() const {
    size_t n_elements = data_.numel() > 0u ? data_.numel()/data_.size(0) : 0u;
    return n_elements * sizeof(T);
  }
private:
  void check_elements(void){
    // check the input for correctness
    ind_t x = this->branch_span(_elements);
    switch (shape_.size()) {
      case 0u:
        break;
      case 1u: // 1 scalar per branch per point
        if (0u == x) x = _elements[0] = 1u;
        if (x > 1u) throw std::runtime_error("1-D data must represent one scalar per point!") ;
        break;
      case 2u: // (P, B*Y)
        if (0u == x) x = _elements[0] = shape_[1]; // one branch with y scalars per point
        if (shape_[1] % x)
          throw std::runtime_error("2-D data requires an integer number of branches!");
        break;
      case 3u: // (P, B, Y)
        if (0u == x) x = _elements[0] = shape_[2];
        if (shape_[2] != x)
          throw std::runtime_error("3-D data requires that the last dimension matches the specified number of elements!");
        break;
      case 4u: // (P, B, V, 3)
        if (3u != shape_[3])
          throw std::runtime_error("4-D data can only be 3-vectors");
        if (0u == x) x = _elements[1] = shape_[2]*3u;
        if (shape_[2]*3u != x)
          throw std::runtime_error("4-D data requires that the last two dimensions match the specified number of vector elements!");
        break;
      case 5u: // (P, B, M, 3, 3)
        if (3u != shape_[3] || 3u != shape_[4])
          throw std::runtime_error("5-D data can only be matrices");
        if (0u == x) x = _elements[2] = shape_[2]*9u;
        if (shape_[2]*9u != x)
          throw std::runtime_error("5-D data requires the last three dimensions match the specified number of matrix elements!");
        break;
      default: // higher dimensions not (yet) supported
        throw std::runtime_error("Interpolator data is expected to be 1- to 5-D");
    }
  }
  bool any_equal_modes_(const ind_t idx, const ind_t b_, const ind_t s_) {
    // since we're probably only using this when the data is provided and
    // most eigenproblem solvers sort their output by eigenvalue magnitude it is
    // most-likely for mode i and mode i+1 to be equal.
    // ∴ search (i,j), (i+1,j+1), (i+2,j+2), ..., i ∈ (0,N], j ∈ (1,N]
    // for each j = i+1, i+2, i+3, ..., i+N-1
    if (b_ < 2) return false;
    // data_ is always 2D: (N,1), (N,B), or (N,Y)
    for (ind_t offset=1; offset < b_; ++offset)
    for (ind_t i=0, j=offset; j < b_; ++i, ++j)
    if (brille::approx_float::vector(s_, data_.ptr(idx, i*s_), data_.ptr(idx, j*s_))) return true;
    // no matches
    return false;
  }
  template<typename I> ind_t branch_span(const std::array<I,3>& e) const {
    return static_cast<ind_t>(e[0])+static_cast<ind_t>(e[1])+static_cast<ind_t>(e[2]);
  }
  element_t<ind_t> count_scalars_vectors_matrices(void) const {
    element_t<ind_t> no{_elements[0], _elements[1]/3u, _elements[2]/9u};
    return no;
  }
  // the 'mixed' variants of the rotate_in_place implementations
  bool rip_real(bArray<T>&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  bool rip_recip(bArray<T>&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  bool rip_axial(bArray<T>&, const PointSymmetry&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  template<class R, class U>
  bool rip_gamma_complex(bArray<T>&, const lattice::LVec<R>&, const GammaTable&, const std::vector<std::array<U,9>>&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const;
  template<class R, class U, class S=T>
  enable_if_t<is_complex<S>::value, bool>
  rip_gamma(bArray<T>& x, const lattice::LVec<R>& q, const GammaTable& gt, const std::vector<std::array<U,9>>& rot, const std::vector<size_t>& r, const std::vector<size_t>& ir, const int nth) const{
    return rip_gamma_complex(x, q, gt, rot, r, ir, nth);
  }
  template<class R, class U, class S=T>
  enable_if_t<!is_complex<S>::value, bool>
  rip_gamma(bArray<T>&, const lattice::LVec<R>&, const GammaTable&, const std::vector<std::array<U,9>>&, const std::vector<size_t>&, const std::vector<size_t>&, const int) const{
    throw std::runtime_error("RotatesLike == Gamma requires complex valued data!");
  }

  // interpolate_at_*
  void interpolate_at_mix(const std::vector<std::vector<ind_t>>&, const std::vector<ind_t>&, const std::vector<double>&, bArray<T>&, const ind_t, const bool) const;
  void interpolate_at_mix(const std::vector<std::vector<ind_t>>&, const std::vector<std::pair<ind_t,double>>&, bArray<T>&, const ind_t, const bool) const;

#ifdef USE_HIGHFIVE
  public:
  template<class HF> std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, bool>
  to_hdf(HF& object, const std::string& entry) const {
    auto group = overwrite_group(object, entry);
    bool ok{true};
    ok &= data_.to_hdf(group, "data");
    group.createDataSet("shape", shape_);
    group.createDataSet("elements", _elements);
    group.createDataSet("rotlike", rotlike_);
    group.createDataSet("lenunit", lenunit_);
    group.createDataSet("costmult", _costmult);
    group.createDataSet("funtype", _funtype);
    return ok;
  }
  [[nodiscard]] bool to_hdf(const std::string& f, const std::string& d, const unsigned p=HighFive::File::OpenOrCreate) const {
    HighFive::File file(f, p);
    return this->to_hdf(file, d);
  }
  template<class HF> static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, Interpolator<T>>
  from_hdf(HF& object, const std::string& entry){
    auto group = object.getGroup(entry);
    auto d = bArray<T>::from_hdf(group, "data");
    //
    shape_t s;
    element_t<ind_t> e, f;
    RotatesLike r;
    LengthUnit l;
    element_t<double> c;
    //
    group.getDataSet("shape").read(s);
    group.getDataSet("elements").read(e);
    group.getDataSet("rotlike").read(r);
    group.getDataSet("lenunit").read(l);
    group.getDataSet("costmult").read(c);
    group.getDataSet("funtype").read(f);
    //
    return {d, s, e, r, l, static_cast<int>(f[0]), static_cast<int>(f[1]), c};
  }
  static Interpolator<T> from_hdf(const std::string& filename, const std::string& dataset){
    HighFive::File file(filename, HighFive::File::ReadOnly);
    return Interpolator<T>::from_hdf(file, dataset);
  }
#endif // USE_HIGHFIVE
};

#include "interpolator_at.tpp"
#include "interpolator_axial.tpp"
#include "interpolator_cost.tpp"
#include "interpolator_gamma.tpp"
#include "interpolator_real.tpp"
#include "interpolator_recip.tpp"

} // namespace brille
#endif
