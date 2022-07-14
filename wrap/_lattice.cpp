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

#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <thread>
#include <utility>

#include "_array.hpp"
#include "_c_to_python.hpp"
#include "lattice_dual.hpp"
#include "utilities.hpp"

namespace py = pybind11;

static brille::LengthUnit lat_type(bool direct){
  return direct ? brille::LengthUnit::angstrom : brille::LengthUnit::inverse_angstrom;
}
static brille::MatrixVectors mat_type(bool row){
  return row ? brille::MatrixVectors::row : brille::MatrixVectors::column;
}

void wrap_lattice(py::module &m){
  using namespace pybind11::literals;
  using namespace brille;
  using namespace brille::lattice;
  // disable pybind11 default doc strings:
  py::options options;
  options.disable_function_signatures();
  // Declare the interface to the superclass Lattice
  py::class_<Lattice<double>> cls(m,"Lattice",R"pbdoc(
  A space-spanning lattice in three dimensions

  A space-spanning lattice in :math:`N` dimensions has :math:`N` basis vectors
  which can be described fully by their :math:`N` lengths and the
  :math:`\sum_1^{N-1} 1` angles between each set of basis vectors, or
  :math:`\sum_1^N 1 = \frac{1}{2}N(N+1)` scalars in total.
  This class stores the basis vectors of the lattice described in an orthonormal space,
  plus the metric of the space, and the equivalent information for the dual of the lattice.

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

  cls.def(py::init(
    [](const py::array_t<double>& lengths,
       const py::array_t<double>& angles,
       const Symmetry& sym,
       const Basis& bas,
       const bool snap,
       const bool dir){
    auto len = np2sa<double,3>(lengths);
    auto ang = np2sa<double,3>(angles);
    return Lattice<double>(lat_type(dir), len, ang, sym, bas, snap);
  }),"basis_vector_lengths"_a, "basis_vector_angles"_a, "symmetry"_a, "basis"_a=Basis(), py::kw_only(), "snap_to_symmetry"_a=true, "real_space"_a=true);
  cls.def(py::init(
    [](const py::array_t<double>& lengths,
       const py::array_t<double>& angles,
       const std::string& sym,
       const Basis& bas,
       const bool snap,
       const bool dir){
    auto len = np2sa<double,3>(lengths);
    auto ang = np2sa<double,3>(angles);
    return Lattice<double>(lat_type(dir),len, ang, sym, bas, snap);
  }),"basis_vector_lengths"_a, "basis_vector_angles"_a, "symmetry_information"_a="P 1", "basis"_a=Basis(), py::kw_only(), "snap_to_symmetry"_a=true, "real_space"_a=true);
  cls.def(py::init(
    [](const py::array_t<double>& lengths,
       const py::array_t<double>& angles,
       const std::string& n,
       const std::string& c,
       const Basis& bas,
       const bool snap,
       const bool dir) {
    auto len = np2sa<double,3>(lengths);
    auto ang = np2sa<double,3>(angles);
    return Lattice<double>(lat_type(dir), len, ang, n, c, bas, snap);
  }), "basis_vector_lengths"_a, "basis_vector_angles"_a, "HM_name"_a, "HM_choice"_a, "basis"_a=Basis(), py::kw_only(), "snap_to_symmetry"_a=true, "real_space"_a=true);
  cls.def(py::init(
    [](const py::array_t<double>& vectors,
       const Symmetry& sym,
       const Basis& bas,
       const bool snap,
       const bool dir,
       const bool row) {
    auto mat = np2sa<double,9>(vectors);
    return Lattice<double>(lat_type(dir), mat, mat_type(row), sym, bas, snap);
  }),"basis_vectors"_a, "symmetry"_a, "basis"_a=Basis(), py::kw_only(), "snap_to_symmetry"_a=true, "real_space"_a=true, "row_vectors"_a=true);
  cls.def(py::init(
    [](const py::array_t<double>& vectors,
       const std::string& sym,
       const Basis& bas,
       const bool snap,
       const bool dir,
       const bool row) {
    auto mat = np2sa<double,9>(vectors);
    return Lattice<double>(lat_type(dir), mat, mat_type(row), sym, bas, snap);
  }),"basis_vectors"_a, "symmetry_information"_a="P 1", "basis"_a=Basis(), py::kw_only(), "snap_to_symmetry"_a=true, "real_space"_a=true, "row_vectors"_a=true);
  cls.def(py::init(
    [](const py::array_t<double>& vectors,
       const std::string& n,
       const std::string& c,
       const Basis& bas,
       const bool snap,
       const bool dir,
       const bool row) {
    auto mat = np2sa<double,9>(vectors);
    return Lattice<double>(lat_type(dir), mat, mat_type(row), n, c, bas, snap);
  }),"basis_vectors"_a, "HM_name"_a, "HM_choice"_a, "basis"_a=Basis(), py::kw_only(), "snap_to_symmetry"_a=true, "real_space"_a=true, "row_vectors"_a=true);

  // accessors
  cls.def_property_readonly("real_vectors",[](const Lattice<double>& lat){
    return sa2np<int>({3,3}, lat.real_basis_vectors());
  });
  cls.def_property_readonly("reciprocal_vectors",[](const Lattice<double>& lat){
    return sa2np<int>({3,3}, lat.reciprocal_basis_vectors());
  });
  cls.def_property_readonly("a",     [](const Lattice<double>& lat){return lat.length(LengthUnit::angstrom, 0);});
  cls.def_property_readonly("b",     [](const Lattice<double>& lat){return lat.length(LengthUnit::angstrom, 1);});
  cls.def_property_readonly("c",     [](const Lattice<double>& lat){return lat.length(LengthUnit::angstrom, 2);});
  cls.def_property_readonly("alpha", [](const Lattice<double>& lat){return lat.angle(LengthUnit::angstrom, 0);});
  cls.def_property_readonly("beta",  [](const Lattice<double>& lat){return lat.angle(LengthUnit::angstrom, 1);});
  cls.def_property_readonly("gamma", [](const Lattice<double>& lat){return lat.angle(LengthUnit::angstrom, 2);});
  cls.def_property_readonly("volume",[](const Lattice<double>& lat){return lat.volume(LengthUnit::angstrom);});

  cls.def_property_readonly("a_star",     [](const Lattice<double>& lat){return lat.length(LengthUnit::inverse_angstrom, 0);});
  cls.def_property_readonly("b_star",     [](const Lattice<double>& lat){return lat.length(LengthUnit::inverse_angstrom, 1);});
  cls.def_property_readonly("c_star",     [](const Lattice<double>& lat){return lat.length(LengthUnit::inverse_angstrom, 2);});
  cls.def_property_readonly("alpha_star", [](const Lattice<double>& lat){return lat.angle(LengthUnit::inverse_angstrom, 0);});
  cls.def_property_readonly("beta_star",  [](const Lattice<double>& lat){return lat.angle(LengthUnit::inverse_angstrom, 1);});
  cls.def_property_readonly("gamma_star", [](const Lattice<double>& lat){return lat.angle(LengthUnit::inverse_angstrom, 2);});
  cls.def_property_readonly("volume_star",[](const Lattice<double>& lat){return lat.volume(LengthUnit::inverse_angstrom);});

  cls.def_property_readonly("bravais",[](const Lattice<double>& l){return l.bravais();});
  cls.def_property_readonly("basis",[](const Lattice<double>& l){return l.basis();});
  cls.def_property("spacegroup",
    [](const Lattice<double>& l){
    return l.spacegroup_symmetry();
    },
    [](Lattice<double>& l, const Symmetry& s){
    return l.spacegroup_symmetry(s);
  });
  cls.def_property_readonly("pointgroup",
    [](const Lattice<double>& l){
    return l.pointgroup_symmetry();
    });
  // The next line would require that the Basis type is exposed to Python
  //cls.def_property_readonly("basis",&Lattice::get_basis);

  cls.def("get_covariant_metric_tensor",[](Lattice<double> &l){
    return sa2np<int>({3,3}, l.metric(LengthUnit::angstrom)); // copies data :/
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

  cls.def("get_contravariant_metric_tensor",[](Lattice<double> &l){
    return sa2np<int>({3,3}, linear_algebra::mat_inverse(l.metric(LengthUnit::angstrom)));
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

  cls.def("vectors",[](Lattice<double>&l, LengthUnit lu){
    return sa2np<int>({3,3}, l.vectors(lu));
  });
  cls.def("vector",[](Lattice<double>&l, LengthUnit lu, const int index){
    return sa2np<int>({3}, l.vector(lu, index));
  });
//  cls.def("angles",[](Lattice<double>&l, LengthUnit lu, AngleUnit au){
//    return sa2np<int>({3}, l.angles(lu, au));
//  });
  cls.def("metric",[](Lattice<double>&l, LengthUnit lu){
    return sa2np<int>({3,3}, l.metric(lu));
  });

  cls.def("__eq__",[](const Lattice<double>& a, const Lattice<double>&b){return a == b;});
  cls.def("__repr__",[](const Lattice<double>& l) {return l.to_verbose_string(AngleUnit::degree);});
  cls.def("str",[](const Lattice<double>& l, LengthUnit lu){return l.to_string(lu);});
}
