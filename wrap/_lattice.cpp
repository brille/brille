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
#include <pybind11/pybind11.h>
#include "_lattice.hpp"

namespace py = pybind11;

void wrap_lattice(py::module &m){
  using namespace pybind11::literals;
  using namespace brille;
  // disable pybind11 default doc strings:
  py::options options;
  options.disable_function_signatures();
  // Declare the interface to the superclass Lattice
  py::class_<Lattice> cls(m,"Lattice",R"pbdoc(
  A space-spanning lattice in three dimensions

  A space-spanning lattice in :math:`N` dimensions has :math:`N` basis vectors
  which can be described fully by their :math:`N` lengths and the
  :math:`\sum_1^{N-1} 1` angles between each set of basis vectors, or
  :math:`\sum_1^N 1 = \frac{1}{2}N(N+1)` scalars in total.
  Storing only their lengths and angles is therefore always more efficient than
  storing all :math:`N^2` components of the basis vectors in *an* orthonormal
  frame.
  This class stores the 3 lengths and 3 angles required to describe a
  3 dimensional space-spanning lattice, plus the volume of the lattice unit cell
  and optional symmetry attributes.

  Attributes
  ----------
  a,b,c : float
        The basis vector lengths
  alpha,beta,gamma : float
        The angles between the basis vectors, internally always in radian
  volume : float
        The volume of the lattice unit cell in units of length cubed.
  bravais : :py:class:`~brille._brille.Bravais`
        The centring type of the lattice
  spacegroup : :py:class:`~brille._brille.Symmetry`
        The Spacegroup symmetry operations of the lattice
  pointgroup : :py:class:`~brille._brille.PointSymmetry`
        The Pointgroup symmetry operations of the lattice
  basis : :py:class:`~brille._brille.Basis`
        The positions of all atoms within the lattice unit cell
  )pbdoc");

  cls.def(py::init<double,double,double,double,double,double,int>(),
          "a"_a,"b"_a,"c"_a,"alpha"_a=brille::halfpi,"beta"_a=brille::halfpi,"gamma"_a=brille::halfpi,
          "HallNumber"_a=1, R"pbdoc(
  Initialize a Lattice from its scalar basis vector lengths and angles

  Parameters
  ----------
  a,b,c : float
        The basis vector lengths
  alpha,beta,gamma : float
        The basis vector angles in radians *or* degrees with values
        :math:`\alpha, \beta, \gamma \in (0,2\pi)` treated as radian
        and all others degrees.
  HallNumber : int, optional
        An integer indication of the Hall symmetry.

  Examples
  --------
  >>> Lattice(a, b, c, alpha, beta, gamma, HallNumber)

  )pbdoc");
  //
  cls.def(py::init([](py::array_t<double> lens, py::array_t<double> angs, int hall) {
    py::buffer_info linfo = lens.request(), ainfo = angs.request();
    if ( linfo.ndim!=1 || ainfo.ndim!=1)
      throw std::runtime_error("Number of dimensions must be one");
    if ( linfo.shape[0] < 3 || ainfo.shape[0] < 3 )
      throw std::runtime_error("(At least) three lengths and angles required.");
    return Lattice((double*) linfo.ptr, linfo.strides,
                   (double*) ainfo.ptr, ainfo.strides, hall);
  }),"lengths"_a,"angles"_a,"HallNumber"_a=1,R"pbdoc(
  Initialize a Lattice from its basis vector lengths and angles

  Parameters
  ----------
  lengths : array_like
        The basis vector lengths
  angles : array_like
        The basis vector angles in radians *or* degrees
  HallNumber : int, optional

  Examples
  --------
  >>> Lattice((a, b, c), [alpha, beta, gamma], HallNumber)
  >>> Lattice([a, b, c], (alpha, beta, gamma), HallNumber)
  )pbdoc");
  //
  cls.def(py::init([](py::array_t<double> vecs, int hall) {
    py::buffer_info info = vecs.request();
    if ( info.ndim!=2 )
      throw std::runtime_error("Number of dimensions must be two");
    if ( info.shape[0] != 3 || info.shape[1] != 3 )
      throw std::runtime_error("Three three-vectors required.");
    return Lattice((double*) info.ptr, info.strides, hall);
  }),"vectors"_a,"HallNumber"_a=1,R"pbdoc(
  Initialize a Lattice from its basis vectors

  Parameters
  ----------
  vectors : matrix_like
        The basis vectors as the rows of a 2-D matrix
  HallNumber : int, optional

  Examples
  --------
  >>> avec = np.array([3.95, 0, 0])
  >>> bvec = np.array([0, 3.95, 0])
  >>> cvec = np.array([0, 0, 12.9])
  >>> Lattice(np.array([avec, bvec, cvec]))

  >>> Lattice(np.array([[3.95, 0, 0], [0, 3.95, 0], [0, 0, 12.9]]))
  )pbdoc");

  // accessors
  cls.def_property_readonly("a",     &Lattice::get_a, R"pbdoc()pbdoc");
  cls.def_property_readonly("b",     &Lattice::get_b);
  cls.def_property_readonly("c",     &Lattice::get_c);
  cls.def_property_readonly("alpha", &Lattice::get_alpha);
  cls.def_property_readonly("beta",  &Lattice::get_beta);
  cls.def_property_readonly("gamma", &Lattice::get_gamma);
  cls.def_property_readonly("volume",&Lattice::get_volume);
  cls.def_property_readonly("bravais",&Lattice::get_bravais_type);
  //cls.def_property("spacegroup",&Lattice::get_spacegroup_symmetry,&Lattice::set_spacegroup_symmetry);
  cls.def_property("spacegroup",
    [](const Lattice& l){
    return l.get_spacegroup_symmetry();
    },
    [](Lattice& l, const Symmetry& s){
    return l.set_spacegroup_symmetry(s);
  });
  //cls.def_property_readonly("pointgroup",&Lattice::get_pointgroup_symmetry);
  cls.def_property_readonly("pointgroup",
    [](const Lattice& l){
    return l.get_pointgroup_symmetry();
    });
  // The next line would require that the Basis type is exposed to Python
  //cls.def_property_readonly("basis",&Lattice::get_basis);

  cls.def("get_covariant_metric_tensor",[](Lattice &l){
    auto result = py::array_t<double, py::array::c_style >({3,3});
    py::buffer_info bi = result.request();
    double *cmt = (double *) bi.ptr;
    l.get_covariant_metric_tensor( cmt );
    return result;
  },R"pbdoc(
  Calculate the covariant metric tensor of the lattice

  Returns
  -------
  matrix_like
        The metric of the lattice

        .. math::
            g_{ij} =
            \begin{pmatrix}
            a^2 & ab\cos\gamma & ac\cos\beta \\
            ab\cos\gamma & b^2 & bc\cos\alpha \\
            ac\cos\beta & bc\cos\alpha & c^2
            \end{pmatrix}


  )pbdoc");

  cls.def("get_contravariant_metric_tensor",[](Lattice &l){
    auto result = py::array_t<double, py::array::c_style >({3,3});
    py::buffer_info bi = result.request();
    l.get_contravariant_metric_tensor((double *) bi.ptr );
    return result;
  },R"pbdoc(
  Calculate the contravariant metric tensor of the lattice

  Returns
  -------
  matrix_like
        The inverse of the metric of the lattice

        .. math::
            g^{ij} =
            \begin{pmatrix}
            a^2 & ab\cos\gamma & ac\cos\beta \\
            ab\cos\gamma & b^2 & bc\cos\alpha \\
            ac\cos\beta & bc\cos\alpha & c^2
            \end{pmatrix}^{-1}

  )pbdoc");

  cls.def_property_readonly("star",[](Lattice &){throw std::runtime_error("Bare Lattices do not have a reciprocal!");});
  cls.def("issame",&Lattice::issame);
  cls.def("__eq__",&Lattice::issame);
  cls.def("isapprox",&Lattice::isapprox);
  cls.def("isstar",[](Lattice &, Lattice ){throw std::runtime_error("Bare Lattices do not have a reciprocal!");});
  cls.def("__repr__",&Lattice::string_repr);

  // Inform pybind11 that the specializations Direct and Reciprocal exist
  py::class_<Direct,Lattice> direct(m,"Direct");
  py::class_<Reciprocal,Lattice> reciprocal(m,"Reciprocal");
  // and define their common extensions
  declare_lattice_methods<Direct>(direct,"Å");
  declare_lattice_methods<Reciprocal>(reciprocal,"Å⁻¹");

  reciprocal.def_property_readonly("B",[](Reciprocal &r){
    auto result = py::array_t<double, py::array::c_style>({3,3});
    py::buffer_info bi = result.request();
    r.get_B_matrix((double *)bi.ptr);
    return result;
  },R"pbdoc(
    Calculate the B matrix

    Equation 3 of `Busing and Levy, Acta Cryst. (1967). 22, 457 <https://doi.org/10.1107/S0365110X67000970>`_
    *Angle Calculations for 3- and 4- Circle X-ray and Neutron Diffractometers*

    Returns
    -------
    matrix_like

      .. math::
          B =
          \begin{pmatrix}
          b_1 & b_2\cos\beta_3 & b_3\cos\beta_2 \\
          0 & b_2\sin\beta_3 & -b_3\sin\beta_2\cos\alpha_1 \\
          0 & 0 & \frac{2\pi}{a_3}
          \end{pmatrix}

      For a lattice with basis vector lengths
      :math:`(a_1, a_2, a_3)` and basis vector angles
      :math:`(\alpha_1, \alpha_2, \alpha_3)` and its reciprocal lattice with
      basis vector lengths :math:`(b_1, b_2, b_3)` and basis vector angles
      :math:`(\beta_1, \beta_2, \beta_3)` and the relationship

      .. math::
          \mathbf{b}_i = 2\pi \frac{\mathbf{a}_j\times \mathbf{a}_k}{\mathbf{a}_i\cdot\left(\mathbf{a}_j\times \mathbf{a}_k\right)}

      for :math:`(i,j,k) \in \left\{(1,2,3),(2,3,1),(3,1,2)\right\}`


    )pbdoc");
}
