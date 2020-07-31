#ifndef __COMMON_GRID_HPP
#define __COMMON_GRID_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <thread>
#include <tuple>

#include "_c_to_python.hpp"
#include "_interpolation_data.hpp"

#include "bz.hpp"
#include "utilities.hpp"

namespace py = pybind11;
typedef long slong;

template<template<class, class> class Grid, class T, class R>
void def_grid_fill(py::class_<Grid<T,R>>& cls){
  using namespace pybind11::literals;
  using Class = Grid<T,R>;

  cls.def("fill",[](Class& cobj,
    py::array_t<T> pyvals, py::array_t<int, py::array::c_style> pyvalelrl,
    py::array_t<R> pyvecs, py::array_t<int, py::array::c_style> pyvecelrl,
    bool sort
  ){
    ArrayVector<T> vals;
    ArrayVector<R> vecs;
    std::vector<size_t> val_sh, vec_sh;
    std::array<element_t, 3> val_el{{0,0,0}}, vec_el{{0,0,0}};
    RotatesLike val_rl, vec_rl;
    size_t count = cobj.vertex_count();
    std::tie(vals,val_sh,val_el,val_rl)=fill_check(pyvals,pyvalelrl,count);
    std::tie(vecs,vec_sh,vec_el,vec_rl)=fill_check(pyvecs,pyvecelrl,count);

    cobj.replace_value_data(vals, val_sh, val_el, val_rl);
    cobj.replace_vector_data(vecs, vec_sh, vec_el, vec_rl);
    if (sort) cobj.sort();
  },
  "values_data"_a, "values_elements"_a,
  "vectors_data"_a, "vectors_elements"_a,
  "sort"_a=false,
R"pbdoc(
  Provide data required for interpolation to the grid without cost information

  This method should probably be followed by `set_cost_info` prior to any
  attempt to interpolate the data in the grid.

  Parameters
  ----------
  values_data : :py:class:`numpy.ndarray`
    The eigenvalue data to be stored in the grid. The first dimension must be
    equal in size to the number of grid-vertices. If two dimensional the second
    dimension is interpreted as all information for a single mode flattened
    and concatenated into (scalars, vectors, matrices) -- in that order.
    If more than two dimensional, the second dimension indexes modes and
    higher dimensions will be flattend *as if row ordered* and must flatten into
    a concatenated list of (scalars, vectors, matrices)
  values_elements: integer vector-like
    A multi-purpose vector containing, in order,
      - the number of scalar-like eigenvalue elements
      - the number of vector-like eigenvalue *elements* (must be :math:`3\times N`)
      - the number of matrix-like eigenvalue *elements* (must be :math:`9\times N`)
      - *how* the vector-like and matrix-like parts transform under application
        of a symmetry operation, see the note about `RotatesLike` below
    Extra entries are ignored
  vectors_data : :py:class:`numpy.ndarray`
    The eigenvector data to be stored in the grid. Same shape restrictions as
    `values_data`
  vectors_elements:
    Like `values_elements` but for the eigenvectors
  sort : logical (default False)
    Whether the equivalent-mode permutations should be (re)determined following
    the update to the flags and weights.


  Note
  ----
  Mapping of integers to :py:mod:`~brille._brille.RotatesLike` values:
    {0: Real, 1: Reciprocal, 2: Axial, 3: Gamma}
  Integer values outside of the mapped range (or missing) are replaced by 0.
)pbdoc");


  cls.def("fill",[](Class& cobj,
    py::array_t<T> pyvals, py::array_t<int> pyvalel, py::array_t<double> pyvalwght,
    py::array_t<R> pyvecs, py::array_t<int> pyvecel, py::array_t<double> pyvecwght,
    bool sort
  ){
    ArrayVector<T> vals;
    ArrayVector<R> vecs;
    std::vector<size_t> val_sh, vec_sh;
    std::array<element_t, 3> val_el{{0,0,0}}, vec_el{{0,0,0}};
    std::array<double,3> val_wght{{1,1,1}}, vec_wght{{1,1,1}};
    RotatesLike val_rl, vec_rl;
    int val_sf{0}, val_vf{0}, vec_sf{0}, vec_vf{0};
    size_t count = cobj.vertex_count();
    std::tie(vals,val_sh,val_el,val_rl,val_sf,val_vf,val_wght)=fill_check(pyvals,pyvalel,pyvalwght,count);
    std::tie(vecs,vec_sh,vec_el,vec_rl,vec_sf,vec_vf,vec_wght)=fill_check(pyvecs,pyvecel,pyvecwght,count);

    cobj.replace_value_data(vals, val_sh, val_el, val_rl);
    cobj.replace_vector_data(vecs, vec_sh, vec_el, vec_rl);
    cobj.set_value_cost_info(val_sf, val_vf, val_wght);
    cobj.set_vector_cost_info(vec_sf, vec_vf, vec_wght);
    if (sort) cobj.sort();
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
    higher dimensions will be flattend *as if row ordered* and must flatten into
    a concatenated list of (scalars, vectors, matrices)
  values_elements: integer vector-like
    A multi-purpose vector containing, in order,
      - the number of scalar-like eigenvalue elements
      - the number of vector-like eigenvalue *elements* (must be :math:`3\times N`)
      - the number of matrix-like eigenvalue *elements* (must be :math:`9\times N`)
      - *how* the vector-like and matrix-like parts transform under application
        of a symmetry operation, see the note about `RotatesLike` below
      - which scalar cost function should be used, see note below
      - which vector cost function should be used, see note below
  values_weights : float, vector-like
    The relative cost weights between scalar-, vector-, and matrix- like
    eigenvalue elements stored in the grid
  vectors_data : :py:class:`numpy.ndarray`
    The eigenvector data to be stored in the grid. Same shape restrictions as
    `values_data`
  vectors_elements:
    Like `values_elements` but for the eigenvectors
  vectors_weights : float, vector-like
    The relative cost weights between scalar-, vector-, and matrix- like
    eigenvector elements stored in the grid
  sort : logical (default False)
    Whether the equivalent-mode permutations should be (re)determined following
    the update to the flags and weights.


  Note
  ----
  Mapping of integers to :py:mod:`~brille._brille.RotatesLike` values:
    {0: Real, 1: Reciprocal, 2: Axial, 3: Gamma}
  Mapping of integers to scalar cost function:
    {0: magnitude(x-y), }
  Mapping of integers to vector cost function:
    {0: sin(hermitian_angle(x, y)),  1: vector_distance(x, y),
     2: 1-vector_product(x, y),      3: vector_angle(x, y),
     4: hermitian_angle(x, y), }
  Integer values outside of the mapped range (or missing) are replaced by 0.
)pbdoc");

  cls.def_property_readonly("bytes_per_point", &Class::bytes_per_point,R"pbdoc(
    Return the memory required per interpolation point *result* in bytes
  )pbdoc");

  cls.def_property_readonly("values",[](Class& cobj){
    const auto & v{cobj.data().values()};
    return av2np_shape(v.data(), v.shape(), false);
  },R"pbdoc(
    Return the stored eigenvalues
  )pbdoc");

  cls.def_property_readonly("vectors",[](Class& cobj){
    const auto & v{cobj.data().vectors()};
    return av2np_shape(v.data(), v.shape(), false);
  },R"pbdoc(
    Return the stored eigenvectors
  )pbdoc");
}

template<template<class, class> class Grid, class T, class R>
void def_grid_ir_interpolate(py::class_<Grid<T,R>>& cls){
  using namespace pybind11::literals;
  using Class = Grid<T,R>;
  cls.def("ir_interpolate_at",[](Class& cobj,
                           py::array_t<double> pyX,
                           const bool& useparallel,
                           const int& threads, const bool& no_move){
    py::buffer_info bi = pyX.request();
    if ( bi.shape[bi.ndim-1] !=3 )
      throw std::runtime_error("Interpolation requires one or more 3-vectors");
    // store shape of X before three-vector dimension for shaping output
    std::vector<ssize_t> preshape;
    for (ssize_t i=0; i < bi.ndim-1; ++i) preshape.push_back(bi.shape[i]);
    // copy the Python X array
    BrillouinZone b = cobj.get_brillouinzone();
    Reciprocal lat = b.get_lattice();
    LQVec<double> qv(lat,(double*)bi.ptr, bi.shape, bi.strides); //memcopy
    // perform the interpolation and rotate and vectors/tensors afterwards
    const int maxth(static_cast<int>(std::thread::hardware_concurrency()));
    int nthreads = (useparallel) ? ((threads < 1) ? maxth : threads) : 1;
    ArrayVector<T> valres;
    ArrayVector<R> vecres;
    std::tie(valres, vecres) = cobj.ir_interpolate_at(qv, nthreads, no_move);
    // copy results to Python arrays and return
    py::array_t<T, py::array::c_style> valout = iid2np(valres, cobj.data().values(),  preshape);
    py::array_t<R, py::array::c_style> vecout = iid2np(vecres, cobj.data().vectors(), preshape);
    return std::make_tuple(valout, vecout);
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
    A 2+ dimensional array with `Q.shape[-1] == 3` of the positions at which an
    interpolated result is required, expressed in units of the reciprocal lattice.
  useparallel : logical (default: False)
    Whether a serial or parallel code should be utilised
  threads : int (default: -1)
    How many OpenMP workers should be utilised; if this value is less than one
    the environment variable `OMP_NUM_THREADS` will be used.
  do_not_move_points: logical (default: false)
    If `True` the provided `Q` points must already lie within the first Brillouin
    zone. No check is made to verify this requirement and if any `Q` lie outside
    of the gridded volume out-of-bounds errors may result in bad data or runtime
    errors.

  Returns
  -------
  A `tuple` of interpolated eigenvalues and eigenvectors at the equivalent
  irreducible first Brillouin zone points.
  The shape of each output will depend on the shape of the data provided to the
  `fill` method. If the filled eigenvalues were of shape
  `[N_grid_points, N_modes, A, ..., B]`, the eigenvectors were of shape
  `[N_grid_points, N_modes, C, ..., D]`, and the provided points of shape
  `[x, y, ..., z, 3]` then the output shapes will be `[x, y, ..., z, A, ..., B]`
  and `[x, y, ..., z, C, ..., D]` for the eigenvalues and eigenvectors,
  respectively.
)pbdoc");

  cls.def("ir_interpolate_at_dw",[](Class& cobj,
                           py::array_t<double> pyX,
                           py::array_t<double> pyM, double temp_k,
                           const bool& useparallel,
                           const int& threads, const bool& no_move){
    py::buffer_info bi = pyX.request();
    if ( bi.shape[bi.ndim-1] !=3 )
      throw std::runtime_error("Interpolation requires one or more 3-vectors");
    // store shape of X before three-vector dimension for shaping output
    std::vector<ssize_t> preshape;
    for (ssize_t i=0; i < bi.ndim-1; ++i) preshape.push_back(bi.shape[i]);
    // copy the Python X array
    BrillouinZone b = cobj.get_brillouinzone();
    Reciprocal lat = b.get_lattice();
    LQVec<double> qv(lat,(double*)bi.ptr, bi.shape, bi.strides); //memcopy
    // perform the interpolation and rotate and vectors/tensors afterwards
    const int maxth(static_cast<int>(std::thread::hardware_concurrency()));
    int nthreads = (useparallel) ? ((threads < 1) ? maxth : threads) : 1;
    ArrayVector<T> valres;
    ArrayVector<R> vecres;
    std::tie(valres, vecres) = cobj.ir_interpolate_at(qv, nthreads, no_move);
    // copy results to Python arrays and return
    py::array_t<T, py::array::c_style> valout = iid2np(valres, cobj.data().values(),  preshape);
    py::array_t<R, py::array::c_style> vecout = iid2np(vecres, cobj.data().vectors(), preshape);
    // calculate the Debye-Waller factor
    auto Wd = av2np_squeeze(cobj.debye_waller(qv, np2vec(pyM), temp_k));
    return std::make_tuple(valout, vecout, Wd);
  },"Q"_a,"M/amu"_a,"temperature/K"_a,"useparallel"_a=false,"threads"_a=-1,"do_not_move_points"_a=false,
R"pbdoc(
  Perform both the linear interpolation and Debye Waller calculation

  Parameters
  ----------
  Q : :py:class:`numpy.ndarray`
    A 2+ dimensional array with `Q.shape[-1] == 3` of the positions at which an
    interpolated result is required, expressed in units of the reciprocal lattice.
  M : vector like
    The masses of atoms in the lattice basis in Atomic Mass Units (amu)
  temperature : float
    The temperature at which to perform the Debye Waller calculation
  useparallel : logical (default: False)
    Whether a serial or parallel code should be utilised
  threads : int (default: -1)
    How many OpenMP workers should be utilised; if this value is less than one
    the environment variable `OMP_NUM_THREADS` will be used.
  do_not_move_points: logical (default: false)
    If `True` the provided `Q` points must already lie within the first Brillouin
    zone. No check is made to verify this requirement and if any `Q` lie outside
    of the gridded volume out-of-bounds errors may result in bad data or runtime
    errors.

  Returns
  -------
  A `tuple` of interpolated eigenvalues and eigenvectors at the equivalent
  irreducible first Brillouin zone points, plus the result of the Debye Waller
  calculation at the input Q points.
  The shape of each output will depend on the shape of the data provided to the
  `fill` method. If the filled eigenvalues were of shape
  `[N_grid_points, N_modes, A, ..., B]`, the eigenvectors were of shape
  `[N_grid_points, N_modes, C, ..., D]`, and the provided points of shape
  `[x, y, ..., z, 3]` then the output shapes will be `[x, y, ..., z, A, ..., B]`
  and `[x, y, ..., z, C, ..., D]` for the eigenvalues and eigenvectors,
  respectively.
)pbdoc");
}

template<template<class, class> class Grid, class T, class R>
void def_grid_interpolate(py::class_<Grid<T,R>>& cls){
  using namespace pybind11::literals;
  using Class = Grid<T,R>;
  cls.def("interpolate_at",[](Class& cobj,
                           py::array_t<double> pyX,
                           const bool& useparallel,
                           const int& threads, const bool& no_move){
    py::buffer_info bi = pyX.request();
    if ( bi.shape[bi.ndim-1] !=3 )
      throw std::runtime_error("Interpolation requires one or more 3-vectors");
    // store shape of X before three-vector dimension for shaping output
    std::vector<ssize_t> preshape;
    for (ssize_t i=0; i < bi.ndim-1; ++i) preshape.push_back(bi.shape[i]);
    // copy the Python X array
    BrillouinZone b = cobj.get_brillouinzone();
    Reciprocal lat = b.get_lattice();
    LQVec<double> qv(lat,(double*)bi.ptr, bi.shape, bi.strides); //memcopy
    // perform the interpolation and rotate and vectors/tensors afterwards
    const int maxth(static_cast<int>(std::thread::hardware_concurrency()));
    int nthreads = (useparallel) ? ((threads < 1) ? maxth : threads) : 1;
    ArrayVector<T> valres;
    ArrayVector<R> vecres;
    std::tie(valres, vecres) = cobj.interpolate_at(qv, nthreads, no_move);
    // copy the results to Python arrays and return
    auto valout = iid2np(valres, cobj.data().values(),  preshape);
    auto vecout = iid2np(vecres, cobj.data().vectors(), preshape);
    return std::make_tuple(valout, vecout);
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
    A 2+ dimensional array with `Q.shape[-1] == 3` of the positions at which an
    interpolated result is required, expressed in units of the reciprocal lattice.
  useparallel : logical (default: False)
    Whether a serial or parallel code should be utilised
  threads : int (default: -1)
    How many OpenMP workers should be utilised; if this value is less than one
    the environment variable `OMP_NUM_THREADS` will be used.
  do_not_move_points: logical (default: false)
    If `True` the provided `Q` points must already lie within the first Brillouin
    zone. No check is made to verify this requirement and if any `Q` lie outside
    of the gridded volume out-of-bounds errors may result in bad data or runtime
    errors.

  Returns
  -------
  A `tuple` of interpolated eigenvalues and eigenvectors at the equivalent first
  Brillouin zone points. The shape of each output will depend on the shape of
  the data provided to the `fill` method. If the filled eigenvalues were of
  shape `[N_grid_points, N_modes, A, ..., B]`, the eigenvectors were of shape
  `[N_grid_points, N_modes, C, ..., D]`, and the provided points of shape
  `[x, y, ..., z, 3]` then the output shapes will be `[x, y, ..., z, A, ..., B]`
  and `[x, y, ..., z, C, ..., D]` for the eigenvalues and eigenvectors,
  respectively.
)pbdoc");
}

template<template<class, class> class Grid, class T, class R>
void def_grid_sort(py::class_<Grid<T,R>>& cls){
  using namespace pybind11::literals;
  using Class = Grid<T,R>;

  cls.def("sort",&Class::sort);

  cls.def("set_flags_weights",
  [](Class& cobj,
    py::array_t<int> pyvalfnc, py::array_t<double> pyvalwght,
    py::array_t<int> pyvecfnc, py::array_t<double> pyvecwght, bool sort){
    std::array<double,3> val_wght{{1,1,1}}, vec_wght{{1,1,1}};
    RotatesLike val_rl, vec_rl;
    int val_sf{0}, val_vf{0}, vec_sf{0}, vec_vf{0};
    std::tie(val_rl,val_sf,val_vf,val_wght)=set_check(pyvalfnc,pyvalwght);
    std::tie(vec_rl,vec_sf,vec_vf,vec_wght)=set_check(pyvecfnc,pyvecwght);
    cobj.set_value_cost_info(val_sf, val_vf, val_wght);
    cobj.set_vector_cost_info(vec_sf, vec_vf, vec_wght);
    if (sort) cobj.sort();
  },
  "values_flags"_a, "values_weights"_a, "vectors_flags"_a, "vectors_weights"_a,
  "sort"_a=false,
  R"pbdoc(
  Set RotatesLike and cost functions plus relative cost weights for the
  values and vectors stored in the object

  Parameters
  ----------
  values_flags : integer, vector-like
    One or more values indicating the :py:mod:`~brille._brille.RotatesLike`
    value for the eigenvalues stored in the object, plus which cost function
    to use when comparing stored eigenvalues at neighbouring grid points for
    scalar- and vector-like eigenvalues.
  values_weights : float, vector-like
    The relative cost weights between scalar-, vector-, and matrix- like
    eigenvalue elements stored in the grid
  vectors_flags : integer, vector-like
    One or more values indicating the :py:mod:`~brille._brille.RotatesLike`
    value for the eigenvectors stored in the object, plus which cost function
    to use when comparing stored eigenvectors at neighbouring grid points for
    scalar- and vector-like eigenvectors.
  vectors_weights : float, vector-like
    The relative cost weights between scalar-, vector-, and matrix- like
    eigenvector elements stored in the grid
  sort : logical (default False)
    Whether the equivalent-mode permutations should be (re)determined following
    the update to the flags and weights.


  Note
  ----
  Mapping of integers to :py:mod:`~brille._brille.RotatesLike` values:
    {0: Real, 1: Reciprocal, 2: Axial, 3: Gamma}
  Mapping of integers to scalar cost function:
    {0: magnitude(x-y), }
  Mapping of integers to vector cost function:
    {0: sin(hermitian_angle(x, y)),  1: vector_distance(x, y),
     2: 1-vector_product(x, y),      3: vector_angle(x, y),
     4: hermitian_angle(x, y), }
  Integer values outside of the mapped range (or missing) are replaced by 0.
  )pbdoc");
}

template<template<class, class> class Grid, class T, class R>
void def_grid_debye_waller(py::class_<Grid<T,R>>& cls){
  using namespace pybind11::literals;
  using Class = Grid<T,R>;

  cls.def("debye_waller",[](Class& cobj, py::array_t<double> pyQ, py::array_t<double> pyM, double temp_k){
    // handle Q
    py::buffer_info bi = pyQ.request();
    if ( bi.shape[bi.ndim-1] !=3 )
      throw std::runtime_error("debye_waller requires one or more 3-vectors");
    BrillouinZone b = cobj.get_brillouinzone();
    Reciprocal lat = b.get_lattice();
    LQVec<double> cQ(lat, (double*)bi.ptr, bi.shape, bi.strides); //memcopy
    return av2np_squeeze(cobj.debye_waller(cQ, np2vec(pyM), temp_k));
  }, "Q"_a, "masses"_a, "Temperature_in_K"_a);
}


#endif
