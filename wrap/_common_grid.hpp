
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <thread>
#include <tuple>

#include "_array.hpp"
#include "_c_to_python.hpp"
#include "_interpolator.hpp"

#include "bz.hpp"
#include "utilities.hpp"

#ifndef WRAP_BRILLE_COMMON_GRID_HPP_
#define WRAP_BRILLE_COMMON_GRID_HPP_
namespace py = pybind11;
namespace br = brille;

template<template<class, class, class> class Grid, class T, class R, class S>
void def_grid_fill(py::class_<Grid<T,R,S>>& cls){
  using namespace pybind11::literals;
  using namespace brille;
  using Class = Grid<T,R,S>;

  cls.def("fill",
  [](Class& cobj,
    py::array_t<T> pyvals, py::array_t<int, py::array::c_style> pyvalelrl,
    py::array_t<R> pyvecs, py::array_t<int, py::array::c_style> pyvecelrl,
    bool sort
  ){
    using namespace brille;
    profile_update("Start of 'fill' operation");
    size_t count = cobj.vertex_count();
    Interpolator<T> vals = fill_check(pyvals,pyvalelrl,count);
    Interpolator<R> vecs = fill_check(pyvecs,pyvecelrl,count);
    cobj.replace_data(vals, vecs);
    profile_update("  End of 'fill' operation");
    if (sort){
      profile_update("Start of 'sort' operation");
      cobj.sort();
      profile_update("  End of 'sort' operation");
    }
  },
  "values_data"_a, "values_elements"_a,
  "vectors_data"_a, "vectors_elements"_a,
  "sort"_a=false,
R"pbdoc(
Provide data required for interpolation to the grid without cost information.

.. Note
.. ----
.. This method should probably be followed by :py:meth:`set_cost_info` prior to
.. any attempt to interpolate the data in the grid.

Parameters
----------
values_data : :py:class:`numpy.ndarray`
  The eigenvalue data to be stored in the grid. The first dimension must be
  equal in size to the number of grid-vertices. If two dimensional the second
  dimension is interpreted as all information for a single mode flattened
  and concatenated into (scalars, vectors, matrices) -- in that order.
  If more than two dimensional, the second dimension indexes modes and
  higher dimensions will be flattened *as if row ordered* and must flatten into
  a concatenated list of (scalars, vectors, matrices).
  If the provided array can be interpreted as a contiguous row-ordered two
  dimensional array it will be used in place, otherwise a copy will be made.
values_elements: integer vector-like
  A multi-purpose vector containing, in order:

  * the number of scalar-like eigenvalue elements,
  * the number of vector-like eigenvalue *elements* (must be :math:`3\times N`),
  * the number of matrix-like eigenvalue *elements* (must be :math:`9\times N`),
  * an integer :py:class:`RotatesLike` value denoting
    *how* the vector-like and matrix-like parts transform under application
    of a symmetry operation (see note below).
  * an integer :py:class:`LengthUnit` value denoting what units
    the vector-like and matrix-like parts are in (see note below).

vectors_data : :py:class:`numpy.ndarray`
  The eigenvector data to be stored in the grid. Same shape restrictions as
  ``values_data``
vectors_elements:
  Like ``values_elements`` but for the eigenvectors
sort : logical (default ``False``)
  Whether the equivalent-mode permutations should be (re)determined following
  the update to the flags and weights.


Note
----
  Mapping of integers to :py:class:`RotatesLike` values:

  +-------+------------------------------------+
  | value | :py:class:`RotatesLike`            |
  +=======+====================================+
  |   0   |               `vector`             |
  +-------+------------------------------------+
  |   1   |            `pseudovector`          |
  +-------+------------------------------------+
  |   2   |               `Gamma`              |
  +-------+------------------------------------+

  Integer values outside of the mapped range (or missing) are replaced by 0.

  Mapping of integers to :py:class:`LengthUnit` values:

  +-------+------------------------------------+
  | value | :py:class:`LengthUnit`             |
  +=======+====================================+
  |   0   |               `none`               |
  +-------+------------------------------------+
  |   1   |             `angstrom`             |
  +-------+------------------------------------+
  |   2   |         `inverse_angstrom`         |
  +-------+------------------------------------+
  |   3   |           `real_lattice`           |
  +-------+------------------------------------+
  |   4   |        `reciprocal_lattice`        |
  +-------+------------------------------------+

  Integer values outside of the mapped range (or missing) are replaced by 3.

)pbdoc");


  cls.def("fill",
  [](Class& cobj,
    py::array_t<T> pyvals, py::array_t<int> pyvalel, py::array_t<double> pyvalwght,
    py::array_t<R> pyvecs, py::array_t<int> pyvecel, py::array_t<double> pyvecwght,
    bool sort
  ){
    profile_update("Start of 'fill' operation with cost information");
    size_t count = cobj.vertex_count();
    Interpolator<T> vals = fill_check(pyvals,pyvalel,pyvalwght,count);
    Interpolator<R> vecs = fill_check(pyvecs,pyvecel,pyvecwght,count);
    cobj.replace_data(vals, vecs);
    profile_update("  End of 'fill' operation with cost information");
    if (sort){
      profile_update("Start of 'sort' operation");
      cobj.sort();
      profile_update("  End of 'sort' operation");
    }
  },
  "values_data"_a, "values_elements"_a, "values_weights"_a,
  "vectors_data"_a, "vectors_elements"_a, "vectors_weights"_a,
  "sort"_a=false,
R"pbdoc(
Provide all data required for interpolation to the grid at once

Parameters
----------
values_data : :py:class:`numpy.ndarray`
  The eigenvalue data to be stored in the grid. The first dimension must be
  equal in size to the number of grid-vertices. If two dimensional the second
  dimension is interpreted as all information for a single mode flattened
  and concatenated into (scalars, vectors, matrices) -- in that order.
  If more than two dimensional, the second dimension indexes modes and
  higher dimensions will be flattened *as if row ordered* and must flatten into
  a concatenated list of (scalars, vectors, matrices).
  If the provided array can be interpreted as a contiguous row-ordered two
  dimensional array it will be used in place, otherwise a copy will be made.
values_elements: integer vector-like
  A multi-purpose vector containing, in order:

  * the number of scalar-like eigenvalue elements,
  * the number of vector-like eigenvalue *elements* (must be :math:`3\times N`),
  * the number of matrix-like eigenvalue *elements* (must be :math:`9\times N`),
  * an integer :py:class:`RotatesLike` value denoting
    *how* the vector-like and matrix-like parts transform under application
    of a symmetry operation (see note below),
  * an integer :py:class:`LengthUnit` value denoting what units
    the vector-like and matrix-like parts are in (see note below)
  * which scalar cost function should be used (see below),
  * which vector cost function should be used (see below).

  See the note below for the meaning of the last three values.
values_weights : float, vector-like
  The relative cost weights between scalar-, vector-, and matrix- like
  eigenvalue elements stored in the grid
vectors_data : :py:class:`numpy.ndarray`
  The eigenvector data to be stored in the grid. Same shape restrictions as
  **values_data**
vectors_elements:
  Like **values_elements** but for the eigenvectors
vectors_weights : float, vector-like
  The relative cost weights between scalar-, vector-, and matrix- like
  eigenvector elements stored in the grid
sort : logical (default ``False``)
  Whether the equivalent-mode permutations should be (re)determined following
  the update to the flags and weights.


Note
----
  Mapping of integers to :py:class:`RotatesLike` values:

  +-------+------------------------------------+
  | value | :py:class:`RotatesLike`            |
  +=======+====================================+
  |   0   |               `vector`             |
  +-------+------------------------------------+
  |   1   |            `pseudovector`          |
  +-------+------------------------------------+
  |   2   |               `Gamma`              |
  +-------+------------------------------------+

  Mapping of integers to :py:class:`LengthUnit` values:

  +-------+------------------------------------+
  | value | :py:class:`LengthUnit`             |
  +=======+====================================+
  |   0   |               `none`               |
  +-------+------------------------------------+
  |   1   |             `angstrom`             |
  +-------+------------------------------------+
  |   2   |         `inverse_angstrom`         |
  +-------+------------------------------------+
  |   3   |           `real_lattice`           |
  +-------+------------------------------------+
  |   4   |        `reciprocal_lattice`        |
  +-------+------------------------------------+

  Integer values outside of the mapped range (or missing) are replaced by 3.

  Mapping of integers to scalar cost function:

  +-------+------------------------------------+
  | value | function(x,y)                      |
  +=======+====================================+
  |   0   | magnitude(x-y)                     |
  +-------+------------------------------------+

  Mapping of integers to vector cost function:

  +-------+------------------------------------+
  | value | function(vec_x, vec_y)             |
  +=======+====================================+
  |   0   | sin(hermitian_angle(vec_x, vec_y)) |
  +-------+------------------------------------+
  |   1   | vector_distance(vec_x, vec_y)      |
  +-------+------------------------------------+
  |   2   | 1 - vector_product(vec_x, vec_y)   |
  +-------+------------------------------------+
  |   3   | vector_angle(vec_x, vec_y)         |
  +-------+------------------------------------+
  |   4   | hermitian_angle(vec_x, vec_y)      |
  +-------+------------------------------------+

  Integer values outside of the mapped range (or missing) are replaced by 0.

)pbdoc");

  cls.def_property_readonly("bytes_per_point", &Class::bytes_per_point,R"pbdoc(
    Return the memory required per interpolation point *result* in bytes
  )pbdoc");

  cls.def_property_readonly("values",[](Class& cobj){
    // extract the values data wrapped with shape information as an Array
    return a2py(cobj.data().values().array());
  },R"pbdoc(
    Return a shared view of the stored eigenvalues
  )pbdoc");

  cls.def_property_readonly("vectors",[](Class& cobj){
    return a2py(cobj.data().vectors().array());
  },R"pbdoc(
    Return a shared view of the stored eigenvectors
  )pbdoc");
}

template<template<class, class, class> class Grid, class T, class R, class S>
void def_grid_ir_interpolate(py::class_<Grid<T,R,S>>& cls){
  using namespace pybind11::literals;
  using Class = Grid<T,R,S>;
  cls.def("ir_interpolate_at",
  [](
    Class& cobj, py::array_t<double> pyX, const bool& useparallel,
    const int& threads, const bool& no_move
  ){
    using namespace brille;
    using namespace brille::lattice;
    profile_update("Start of 'ir_interpolate_at' operation");
    brille::Array2<double> bX = brille::py2a2(pyX);
    auto qv = LQVec<double>(cobj.get_brillouinzone().get_lattice(), bX);
    profile_update("Q array wrapped for C++ use");
    if (qv.size(qv.ndim()-1) != 3)
      throw std::runtime_error("Interpolation requires one or more 3-vectors");
    // perform the interpolation and rotate and vectors/tensors afterwards
    const int maxth(static_cast<int>(std::thread::hardware_concurrency()));
    int nthreads = (useparallel) ? ((threads < 1) ? maxth : threads) : 1;
    if (no_move) {
      auto [val, vec] = cobj.template ir_interpolate_at<true>(qv, nthreads);
      profile_update("  End of 'ir_interpolate_at' operation");
      return std::make_tuple(brille::a2py(val), brille::a2py(vec));
    } else {
      auto [val, vec] = cobj.ir_interpolate_at(qv, nthreads);
      profile_update("  End of 'ir_interpolate_at' operation");
      return std::make_tuple(brille::a2py(val), brille::a2py(vec));
    }
  },"Q"_a,"useparallel"_a=false,"threads"_a=-1,"do_not_move_points"_a=false,
R"pbdoc(
  Perform linear interpolation of the stored data at irreducible equivalent points

  The irreducible first Brillouin zone is the part of reciprocal space which is
  invariant under application of the integer translations *and* the pointgroup
  operations of a reciprocal space lattice. This method finds points equivalent
  to the input within the irreducible first Brillouin zone and then interpolates
  pre-stored information to provide an estimate at the found positions.

  Parameters
  ----------
  Q : :py:class:`numpy.ndarray`
    A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
    which an interpolated result is required, expressed in units of the
    reciprocal lattice.
  useparallel : bool, optional
    Whether a serial or parallel code should be utilised
  threads : int, optional
    How many OpenMP workers should be utilised; if this value is less than one
    the environment variable ``OMP_NUM_THREADS`` will be used.
  do_not_move_points: bool, optional
    If ``True`` the provided **Q** points must already lie within the first Brillouin
    zone. No check is made to verify this requirement and if any **Q** lie outside
    of the gridded volume out-of-bounds errors may result in bad data or runtime
    errors.

  Returns
  -------
  tuple
      The interpolated eigenvalues and eigenvectors at the equivalent
      irreducible first Brillouin zone points.
      The shape of each output will depend on the shape of the data provided to
      the :py:meth:`~brille._brille.BZTrellisQdc.fill` method. i
      If the filled eigenvalues were of shape
      ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
      ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
      ``[N_Q_points, 3]`` then the output shapes will be
      ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
      for the eigenvalues and eigenvectors, respectively.
)pbdoc");

//  cls.def("ir_interpolate_at_dw",
//  [](
//    Class& cobj, py::array_t<double> pyX, py::array_t<double> pyM,
//    double temp_k, const bool& useparallel, const int& threads, const bool& no_move
//  ){
//    using namespace brille;
//    using namespace brille::lattice;
//    profile_update("Start of 'ir_interpolate_at_dw' operation");
//    brille::Array2<double> bX = brille::py2a2(pyX);
//    auto qv = LQVec<double>(cobj.get_brillouinzone().get_lattice(), bX);
//    profile_update("Q array wrapped for C++ use");
//    if (qv.size(qv.ndim()-1) != 3)
//      throw std::runtime_error("Interpolation requires one or more 3-vectors");
//    // perform the interpolation and rotate and vectors/tensors afterwards
//    const int maxth(static_cast<int>(std::thread::hardware_concurrency()));
//    int nthreads = (useparallel) ? ((threads < 1) ? maxth : threads) : 1;
//    auto [val, vec] = cobj.ir_interpolate_at(qv, nthreads, no_move);
//    profile_update("Interpolated values found");
//    // calculate the Debye-Waller factor
//    auto Wd = brille::a2py(cobj.debye_waller(qv, np2vec(pyM), temp_k));
//    profile_update("  End of 'ir_interpolate_at_dw' operation");
//    return std::make_tuple(brille::a2py(val), brille::a2py(vec), Wd);
//  },"Q"_a,"M/amu"_a,"temperature/K"_a,"useparallel"_a=false,"threads"_a=-1,"do_not_move_points"_a=false,
//R"pbdoc(
//  Perform both the linear interpolation and Debye Waller calculation
//
//  Parameters
//  ----------
//  Q : :py:class:`numpy.ndarray`
//    A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
//    which an interpolated result is required, expressed in units of the
//    reciprocal lattice.
//  M : vector like
//    The masses of atoms in the lattice basis in Atomic Mass Units (amu)
//  temperature : float
//    The temperature at which to perform the Debye-Waller calculation
//  useparallel : bool, optional
//    Whether a serial or parallel code should be utilised
//  threads : int, optional
//    How many OpenMP workers should be utilised; if this value is less than one
//    the environment variable ``OMP_NUM_THREADS`` will be used.
//  do_not_move_points: bool, optional
//    If ``True`` the provided **Q** points must already lie within the first Brillouin
//    zone. No check is made to verify this requirement and if any **Q** lie outside
//    of the gridded volume out-of-bounds errors may result in bad data or runtime
//    errors.
//
//  Returns
//  -------
//  tuple
//      The interpolated eigenvalues and eigenvectors at the equivalent
//      irreducible first Brillouin zone points, and the result of the
//      Debye-Waller calculation at the input Q points.
//      The shape of each output will depend on the shape of the data provided to
//      the :py:meth:`~brille._brille.BZTrellisQdc.fill` method.
//      If the filled eigenvalues were of shape
//      ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
//      ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
//      ``[N_Q_points, 3]`` then the output shapes will be
//      ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
//      for the eigenvalues and eigenvectors, respectively.
//      The Debye Waller result is always a ``[N_Q_points]`` vector.
//)pbdoc");
}

template<template<class, class, class> class Grid, class T, class R, class S>
void def_grid_interpolate(py::class_<Grid<T,R,S>>& cls){
  using namespace pybind11::literals;
  using Class = Grid<T,R,S>;
  cls.def("interpolate_at",
  [](
    Class& cobj, py::array_t<double> pyX, const bool& useparallel,
    const int& threads, const bool& no_move
  ){
    using namespace brille;
    using namespace brille::lattice;
    profile_update("Start of 'interpolate_at' operation");
    brille::Array2<double> bX = brille::py2a2(pyX);
    auto qv = LQVec<double>(cobj.get_brillouinzone().get_lattice(), bX);
    profile_update("Q array wrapped for C++ use");
    if (qv.size(qv.ndim()-1) != 3)
      throw std::runtime_error("Interpolation requires one or more 3-vectors");
    // perform the interpolation and rotate and vectors/tensors afterwards
    const int maxth(static_cast<int>(std::thread::hardware_concurrency()));
    int nthreads = (useparallel) ? ((threads < 1) ? maxth : threads) : 1;
    if (no_move) {
      auto [val, vec] = cobj.template interpolate_at<true>(qv, nthreads);
      profile_update("  End of 'interpolate_at' operation");
      return std::make_tuple(brille::a2py(val), brille::a2py(vec));
    } else {
      auto [val, vec] = cobj.interpolate_at(qv, nthreads);
      profile_update("  End of 'interpolate_at' operation");
      return std::make_tuple(brille::a2py(val), brille::a2py(vec));
    }
      },"Q"_a,"useparallel"_a=false,"threads"_a=-1,"do_not_move_points"_a=false,
R"pbdoc(
  Perform linear interpolation of the stored data at equivalent points

  The first Brillouin zone is the part of reciprocal space which is invariant
  under application of the integer translations of a reciprocal space lattice.
  This method finds points equivalent to the input within the first Brillouin
  zone and then interpolates pre-stored information to provide an estimate at
  the found positions.

  Parameters
  ----------
  Q : :py:class:`numpy.ndarray`
    A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
    which an interpolated result is required, expressed in units of the
    reciprocal lattice.
  useparallel : bool, optional
    Whether a serial or parallel code should be utilised
  threads : int, optional
    How many OpenMP workers should be utilised; if this value is less than one
    the environment variable ``OMP_NUM_THREADS`` will be used.
  do_not_move_points: bool, optional
    If ``True`` the provided **Q** points must already lie within the first Brillouin
    zone. No check is made to verify this requirement and if any **Q** lie outside
    of the gridded volume out-of-bounds errors may result in bad data or runtime
    errors.

  Returns
  -------
  tuple
      The interpolated eigenvalues and eigenvectors at the equivalent
      first Brillouin zone points.
      The shape of each output will depend on the shape of the data provided to
      the :py:meth:`~brille._brille.BZTrellisQdc.fill` method.
      If the filled eigenvalues were of shape
      ``[N_grid_points, N_modes, A, ..., B]``, the eigenvectors were of shape
      ``[N_grid_points, N_modes, C, ..., D]``, and the provided points of shape
      ``[N_Q_points, 3]`` then the output shapes will be
      ``[N_Q_points, N_modes, A, ..., B]`` and ``[N_Q_points, N_modes, C, ..., D]``
      for the eigenvalues and eigenvectors, respectively.
)pbdoc");
}

template<template<class, class, class> class Grid, class T, class R, class S>
void def_grid_sort(py::class_<Grid<T,R,S>>& cls){
  using namespace pybind11::literals;
  using namespace brille;
  using Class = Grid<T,R,S>;

  cls.def("sort",&Class::sort);

  cls.def("set_flags_weights",
  [](Class& cobj,
    py::array_t<int> pyvalfnc, py::array_t<double> pyvalwght,
    py::array_t<int> pyvecfnc, py::array_t<double> pyvecwght, bool sort){
    std::array<double,3> val_wght{{1,1,1}}, vec_wght{{1,1,1}};
    RotatesLike val_rl, vec_rl;
    LengthUnit val_lu, vec_lu;
    int val_sf{0}, val_vf{0}, vec_sf{0}, vec_vf{0};
    std::tie(val_rl,val_lu,val_sf,val_vf,val_wght)=set_check(pyvalfnc,pyvalwght);
    std::tie(vec_rl,vec_lu,vec_sf,vec_vf,vec_wght)=set_check(pyvecfnc,pyvecwght);
    cobj.set_value_cost_info(val_sf, val_vf, val_wght);
    cobj.set_vector_cost_info(vec_sf, vec_vf, vec_wght);
    if (sort) cobj.sort();
  },
  "values_flags"_a, "values_weights"_a, "vectors_flags"_a, "vectors_weights"_a,
  "sort"_a=false,
  R"pbdoc(
  Set :py:class:`~brille._brille.RotatesLike`, `~brille._brille.LengthUnit`
  and cost functions plus relative cost weights for the values and vectors
  stored in the object

  Parameters
  ----------
  values_flags : integer, vector-like
    One or more values indicating the :py:class:`~brille._brille.RotatesLike`
    value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
    value, plus which cost function to use when comparing stored eigenvalues at
    neighbouring grid points for scalar- and vector-like eigenvalues.
  values_weights : float, vector-like
    The relative cost weights between scalar-, vector-, and matrix- like
    eigenvalue elements stored in the grid
  vectors_flags : integer, vector-like
    One or more values indicating the :py:class:`~brille._brille.RotatesLike`
    value for the eigenvalues stored in the object, the `~brille._brille.LengthUnit`
    value, plus which cost function to use when comparing stored eigenvectors at
    neighbouring grid points for scalar- and vector-like eigenvectors.
  vectors_weights : float, vector-like
    The relative cost weights between scalar-, vector-, and matrix- like
    eigenvector elements stored in the grid
  sort : bool, optional
    Whether the equivalent-mode permutations should be (re)determined following
    the update to the flags and weights.


  Note
  ----
    Mapping of integers to :py:class:`~brille._brille.RotatesLike` values:

    +-------+------------------------------------+
    | value | :py:class:`RotatesLike`            |
    +=======+====================================+
    |   0   |               `vector`             |
    +-------+------------------------------------+
    |   1   |            `pseudovector`          |
    +-------+------------------------------------+
    |   2   |               `Gamma`              |
    +-------+------------------------------------+
  
    Mapping of integers to :py:class:`LengthUnit` values:

    +-------+------------------------------------+
    | value | :py:class:`LengthUnit`             |
    +=======+====================================+
    |   0   |               `none`               |
    +-------+------------------------------------+
    |   1   |             `angstrom`             |
    +-------+------------------------------------+
    |   2   |         `inverse_angstrom`         |
    +-------+------------------------------------+
    |   3   |           `real_lattice`           |
    +-------+------------------------------------+
    |   4   |        `reciprocal_lattice`        |
    +-------+------------------------------------+

    Mapping of integers to scalar cost function:
  
    +-------+------------------------------------+
    | value | function(x,y)                      |
    +=======+====================================+
    |   0   | magnitude(x-y)                     |
    +-------+------------------------------------+
  
    Mapping of integers to vector cost function:
  
    +-------+------------------------------------+
    | value | function(vec_x, vec_y)             |
    +=======+====================================+
    |   0   | sin(hermitian_angle(vec_x, vec_y)) |
    +-------+------------------------------------+
    |   1   | vector_distance(vec_x, vec_y)      |
    +-------+------------------------------------+
    |   2   | 1 - vector_product(vec_x, vec_y)   |
    +-------+------------------------------------+
    |   3   | vector_angle(vec_x, vec_y)         |
    +-------+------------------------------------+
    |   4   | hermitian_angle(vec_x, vec_y)      |
    +-------+------------------------------------+

    Integer values outside of the mapped range (or missing) are replaced by 0.
  )pbdoc");
}
//
//template<template<class, class, class> class Grid, class T, class R, class S>
//void def_grid_debye_waller(py::class_<Grid<T,R,S>>& cls){
//  using namespace pybind11::literals;
//  using namespace brille;
//  using namespace brille::lattice;
//  using Class = Grid<T,R,S>;
//
//  cls.def("debye_waller",[](Class& cobj, py::array_t<double> pyX, py::array_t<double> pyM, double temp_k){
//    // handle Q
//    brille::Array2<double> bX = brille::py2a2(pyX);
//    auto qv = LQVec<double>(cobj.get_brillouinzone().get_lattice(), bX);
//    if (qv.size(qv.ndim()-1) != 3)
//      throw std::runtime_error("Interpolation requires one or more 3-vectors");
//    return brille::a2py(cobj.debye_waller(qv, np2vec(pyM), temp_k));
//  }, "Q"_a, "masses"_a, "Temperature_in_K"_a,
//  R"pbdoc(
//  Calculate the Debye-Waller factor by integrating over the Brillouin zone grid
//
//  Parameters
//  ----------
//  pyX : float, array_like
//    A two dimensional array with ``Q.shape[1] == 3`` containing the positions at
//    which to calculate the Debye-Waller factor
//  pyM : float, array_like
//    The atom masses in Atomic Mass Units, matching the atom types provided to
//    the :py:class:`Direct` constructor for the crystal basis
//  temp_K : float
//    The temperature at which to calculate the Debye-Waller factor in Kelvin
//
//
//  Returns
//  -------
//  :py:class:`numpy.ndarray`
//    A one-dimensional array of the Debye-Waller factor evaluated at each input
//
//  )pbdoc");
//}

#ifdef USE_HIGHFIVE
template<template<class, class, class> class Grid, class T, class R, class S>
void def_grid_hdf_interface(py::class_<Grid<T,R,S>>& cls, const std::string& default_entry){
  using namespace pybind11::literals;
  using namespace brille;
  using namespace HighFive;
  using Class = Grid<T,R, S>;
  const std::string default_flag("ac");

  cls.def("to_file",[](Class& cobj, const std::string& filename, const std::string& entry, const std::string& flags){
    unsigned flag{0u};
    if (flags.find('r') != std::string::npos) flag |= File::ReadOnly;
    if (flags.find('x') != std::string::npos) flag |= File::Excl;
    if (flags.find('a') != std::string::npos) flag |= File::ReadWrite;
    if (flags.find('c') != std::string::npos) flag |= File::Create;
    if (flags.find('t') != std::string::npos) flag |= File::Truncate;
    info_update("Provided flags ",flags," is translated to ",flag);
    return cobj.to_hdf(filename, entry, flag);
  }, "filename"_a, "entry"_a=default_entry, "flags"_a=default_flag.c_str(),
  R"pbdoc(
  Save the object to an HDF5 file

  Parameters
  ----------
  filename : str
    The full path specification for the file to write into
  entry: str
    The group path, e.g., "my/cool/grid", where to write inside the file,
    with a default equal to the object Class name
  flags: str
    The HDF5 permissions to use when opening the file. Default 'a' writes to an
    existing file -- if `entry` exists in the file it is overwritten.

  Note
  ----
  Possible `flags` are:

  +---------+-------------------------+----------------+
  | `flags` | meaning                 | HDF equivalent |
  +=========+=========================+================+
  | 'r'     | read                    | H5F_ACC_RDONLY |
  +---------+-------------------------+----------------+
  | 'x'     | write, error if exists  | H5F_ACC_EXCL   |
  +---------+-------------------------+----------------+
  | 'a'     | write, append to file   | H5F_ACC_RDWR   |
  +---------+-------------------------+----------------+
  | 'c'     | write, error if exists  | H5F_ACC_CREAT  |
  +---------+-------------------------+----------------+
  | 't'     | write, replace existing | H5F_ACC_TRUNC  |
  +---------+-------------------------+----------------+


  Returns
  -------
  bool
    Indication of writing success.

  )pbdoc");

    // how do we define this static?
  cls.def_static("from_file",[](const std::string& filename, const std::string& entry){
        return Class::from_hdf(filename, entry);
      }, "filename"_a, "entry"_a=default_entry,
      R"pbdoc(
  Save the object to an HDF5 file

  Parameters
  ----------
  filename : str
    The full path specification for the file to read from
  entry: str
    The group path, e.g., "my/cool/grid", where to read from inside the file,
    with a default equal to the object Class name

  Returns
  -------
  clsObj

  )pbdoc");
}
#endif // USE_HIGHFIVE

#endif
