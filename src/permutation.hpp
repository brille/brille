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

/* This file holds implementations of Linear Assignment Problem solvers which
   can be used to find the sorting permutations at grid points.              */
#ifndef BRILLE_PERMUTATION_H
#define BRILLE_PERMUTATION_H
#include <complex>
#include <cstdint>
#include "array.hpp"
#include "array_latvec.hpp" // defines bArray
#include "sorting_status.hpp"
#include "munkres.hpp"
#include "lapjv.hpp"
#include "smp.hpp"
namespace brille {

/*! \brief Type information for cost-matrix elements used by Linear Assignment
     Problem solvers.

Linear Assignment Problem algorithms can only handle real-valued costs.
Since it is desirable to make assignments of complex-valued data we
need means by which to identify the underlying real data type.

| template typename | type | max |
| T | T | std::numeric_limits<T>::max() |
| std::complex<T> | T | std::numeric_limits<T>::max() |
*/
template<class T> struct CostTraits{
  using type = T;
  constexpr static T max = (std::numeric_limits<T>::max)();
};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T> struct CostTraits<std::complex<T>>{
  using type = T;
  constexpr static T max = (std::numeric_limits<T>::max)();
};
#endif



/*! \brief Use Munkres' Assignment algorithm to determine a permutation

For two arrays of data, located in memory at `centre` and `neighbour`, and each
representing `Nobj` sets of  `Nel[0]` scalars, `Nel[1]` vector elements, and `Nel[2]`
matrix elements, determine the sorting permutation that maps the elements at
`centre` onto the same global mapping as the sorting permutation already stored
at `permutations[neighbour_idx]`. The resultant global permutation is stored
at `permutations[centre_idx]`.

Each array of data to be compared must be formatted as:

      [ 0{(0…Nel[0]-1)(0…Nel[1]-1)(0…Nel[2]-1)}
        1{(0…Nel[0]-1)(0…Nel[1]-1)(0…Nel[2]-1)}
        ⋮
        Nobj-1{(0…Nel[0]-1)(0…Nel[1]-1)(0…Nel[2]-1)} ]

The function constructs an `Nobj`×`Nobj` cost matrix, where each element is
given by

      Cᵢⱼ = Wscl×∑ₖ₋₀ᴺˢᶜˡ√(centre[i,k]-neighbour[j,k])
          + Wvec×𝔣ᵥ(centre_vec[i],neighbour_vec[j])
          + Wmat×𝔣ₘ(centre_mat[i],neighbour_mat[j])

where `Wscl`, `Wvec`, and `Wmat` are weight factors for adjusting the relative
cost of the scalar, vector, and matrix differences, respectively;
𝔣ᵥ is the vector cost function; and 𝔣ₘ is the matrix cost function.

The vector cost function is selected by `vec_cost_func` to be one of:
the angle between the vectors, the length of the difference between the vectors,
or one minus the inner product of the vectors.

The matrix cost function is the Frobenius norm of the difference between the
two matrices.

@param centre The values to be sorted by the determined permutation
@param neighbour The values against which the `centre` values are compared
@param Nel The number of scalars, vector elements, and matrix elements per object
@param Wscl The cost weight for scalar elements
@param Wvec The cost weight for vectors
@param Wmat The cost weight for matrices
@param span The number of elements per object, `Nel[0]+Nel[1]+Nel[2]`
@param Nobj The number of objects at each of `centre` and `neighbour`
@param[out] permutations Contains the global permutation sorting for the
                         neighbour and is where the centre permutation is stored
@param centre_idx The index into `permutations` where the output is stored
@param neighbour_idx The index into `permutations` to find the neighbour permutation
@param vec_cost_func Used to select the eigenvector cost function:
                      0 --> `vector_angle`
                      1 --> `vector_distance`
                      2 --> `1-brille::utils::vector_product`
                      3 --> `abs(sin(brille::utils::hermitian_angle))`
@returns `true` if the permutation was assigned successfully, otherwise `false`
*/
template<class T, class R, class I,
          typename=typename std::enable_if<std::is_same<typename CostTraits<T>::type, R>::value>::type
        >
bool munkres_permutation(const T* centre, const T* neighbour, const std::array<I,3>& Nel,
                         const R Wscl, const R Wvec, const R Wmat,
                         const size_t span, const size_t Nobj,
                         bArray<size_t>& permutations,
                         const size_t centre_idx, const size_t neighbour_idx,
                         const int vec_cost_func = 0
                       ){
/* An earlier version of this function took `centre` and `neighbour` arrays
   which were effectively [span, Nobj] 2D arrays. This has the unfortunate
   property that individual objects were not contiguous in memory but instead
   had a stride equal to Nobj.
   Now this function requires arrays which are [Nobj, span], placing each
   object in contiguous memory and eliminating the need to copy sub-object
   vector/matrices into contiguous memory before calling subroutines.
*/
// initialize variables
R s_cost{0}, v_cost{0}, m_cost{0};
brille::assignment::Munkres<R> munkres(Nobj);
size_t *assignment = new size_t[Nobj]();
const T *c_i, *n_j;
// calculate costs and fill the Munkres cost matrix
for (size_t i=0; i<Nobj; ++i){
  for (size_t j=0; j<Nobj; ++j){
    c_i = centre+i*span;
    n_j = neighbour+j*span;
    s_cost = R(0);
    if (Nel[0]){
      for (size_t k=0; k<Nel[0]; ++k)
        s_cost += magnitude(c_i[k] - n_j[k]);
      c_i += Nel[0];
      n_j += Nel[0];
    }
    if (Nel[1]){
      /* As long as the eigenvectors at each grid point are the eigenvectors of
         a Hermitian matrix (for which H=H†, † ≡ complex conjugate transpose)
         then the *entire* Nel[1] dimensional eigenvectors at a given grid
         point *are* orthogonal, and equivalent-mode eigenvectors at neighbouring
         grid points should be least-orthogonal.

         If ⃗a and ⃗b are eigenvectors of D with eigenvalues α and β,
         respectively, then (using Einstein summation notation)
         <D ⃗a, ⃗b > = (Dᵢⱼaⱼ)† bᵢ = Dᵢⱼ* aⱼ* bᵢ = aⱼ* Dᵢⱼ* bᵢ = < ⃗a, D† ⃗b >
         and, since D ⃗a = α ⃗a and D† ⃗b = D ⃗b = β ⃗b,
         < α ⃗a, ⃗b > = < ⃗a, β ⃗b > → (α-β)< ⃗a, ⃗b > = 0.
         ∴ < ⃗a, ⃗b > = 0 for non-degenerate solutions to D ϵ = ω² ϵ.
         For degenerate solutions, eigenvalue solvers still tend to return an
         arbitrary linear combination c₀ ⃗a + s₁ ⃗b and s₀ ⃗a + c₁ ⃗b which are
         still orthogonal.

         It stands to reason that at nearby grid points the orthogonality will
         only be approximate, and we can instead try to identify the any
         least-orthogonal modes as being equivalent.
      */
      switch(vec_cost_func){
        case 0: v_cost = std::abs(std::sin(brille::utils::hermitian_angle(Nel[1], c_i, n_j))); break;
        case 1: v_cost = brille::utils::vector_distance(Nel[1], c_i, n_j); break;
        case 2: v_cost = 1-brille::utils::vector_product(Nel[1], c_i, n_j); break;
        case 3: v_cost = brille::utils::vector_angle(Nel[1], c_i, n_j); break;
      }
      c_i += Nel[1];
      n_j += Nel[1];
    }
    if (Nel[2]){
      I nel2 = std::sqrt(Nel[2]);
      if (nel2*nel2 != Nel[2])
        throw std::runtime_error("Non-square matrix in munkres_permutation");
      m_cost = brille::utils::frobenius_distance(nel2, c_i, n_j);
    }
    // for each i we want to determine the cheapest j
    munkres.get_cost()[i*Nobj+j] = Wscl*s_cost + Wvec*v_cost + Wmat*m_cost;
  }
}
// use the Munkres' algorithm to determine the optimal assignment
munkres.run_assignment();
if (!munkres.get_assignment(assignment)){
  delete[] assignment;
  return false;
}
/* use the fact that the neighbour objects have already had their global
   permutation saved into `permutations` to determine the global permuation
   for the centre objects too; storing the result into `permutations` as well.
*/
brille::ind_t nind = static_cast<brille::ind_t>(neighbour_idx);
brille::ind_t cind = static_cast<brille::ind_t>(centre_idx);
for (brille::ind_t i=0; i<Nobj; ++i)
  for (brille::ind_t j=0; j<Nobj; ++j)
    if (permutations.val(nind, i) == assignment[j]) permutations.val(cind, i) = j;
delete[] assignment;
return true;
}


/*! \brief Use Junker-Volgenant algorithm to determine a permutation

For two arrays of data, located in memory at `centre` and `neighbour`, and each
representing `Nobj` sets of  `Nel[0]` scalars, `Nel[1]` eigenvector elements,
`Nel[1]` vector elements, and `Nel[2]` matrix elements,
determine the sorting permutation that maps the elements at `centre` onto the
same global mapping as the sorting permutation already stored at
`permutations[neighbour_idx]`. The resultant global permutation is stored at
`permutations[centre_idx]`.

Each array of data to be compared must be formatted as:

      [ 0{(0…Nel[0]-1)(0…Nel[1]-1)(0…Nel[2]-1)}
        1{(0…Nel[0]-1)(0…Nel[1]-1)(0…Nel[2]-1)}
        ⋮
        Nobj-1{(0…Nel[0]-1)(0…Nel[1]-1)(0…Nel[2]-1)} ]

The function constructs an `Nobj`×`Nobj` cost matrix, where each element is
given by

      Cᵢⱼ = Wscl×∑ₖ₋₀ᴺˢᶜˡ√(centre[i,k]-neighbour[j,k])
          + Wvec×𝔣ᵥ(centre_vec[i],neighbour_vec[j])
          + Wmat×𝔣ₘ(centre_mat[i],neighbour_mat[j])

where `Wscl`, `Wvec`, and `Wmat` are weight factors for adjusting the
relative cost of the scalar, vector, and matrix differences,
respectively; 𝔣ᵥ is the vector cost
function; and 𝔣ₘ is the matrix cost function.

The vector cost function is selected by `vec_cost_func` to be one of:
the absolute value of the sine of the hermitian angle between eigenvectors,
the distance between eigenvectors, or one minus the inner product between
eigenvectors.

The matrix cost function is the Frobenius norm of the difference between the
two matrices.

@param centre The values to be sorted by the determined permutation
@param neighbour The values against which the `centre` values are compared
@param Nel[0] The number of scalars per object
@param Nel[1] The number of vector elements per object
@param Nel[2] The square root of the number of matrix elements per object
@param Wscl The cost weight for scalar elements
@param Wvec The cost weight for vectors
@param Wmat The cost weight for matrices
@param span The number of elements per object, `Nel[0]+Nel[1]+Nel[2]`
@param Nobj The number of objects at each of `centre` and `neighbour`
@param[out] permutations Contains the global permutation sorting for the
                         neighbour and is where the centre permutation is stored
@param centre_idx The index into `permutations` where the output is stored
@param neighbour_idx The index into `permutations` to find the neighbour permutation
@param vec_cost_func Used to select the eigenvector cost function:
                      0 --> `abs(sin(brille::utils::hermitian_angle))`
                      1 --> `vector_distance`
                      2 --> `1-brille::utils::vector_product`
                      3 --> `vector_angle`
@returns `true` if the permutation was assigned successfully, otherwise `false`
*/
template<class T, class R, class I,
          typename=typename std::enable_if<std::is_same<typename CostTraits<T>::type, R>::value>::type
        >
bool jv_permutation(const T* centre, const T* neighbour, const std::array<I,3>& Nel,
                    const R Wscl, const R Wvec, const R Wmat,
                    const size_t span, const size_t Nobj,
                    bArray<size_t>& permutations,
                    const size_t centre_idx, const size_t neighbour_idx,
                    const int vec_cost_func = 0
                   ){
/* An earlier version of this function took `centre` and `neighbour` arrays
   which were effectively [span, Nobj] 2D arrays. This has the unfortunate
   property that individual objects were not contiguous in memory but instead
   had a stride equal to Nobj.
   Now this function requires arrays which are [Nobj, span], placing each
   object in contiguous memory and eliminating the need to copy sub-object
   vector/matrices into contiguous memory before calling subroutines.
*/
// initialize variables
R s_cost{0}, v_cost{0}, m_cost{0};
R *cost=nullptr, *usol=nullptr, *vsol=nullptr;
cost = new R[Nobj*Nobj];
usol = new R[Nobj];
vsol = new R[Nobj];
int *rowsol=nullptr, *colsol=nullptr;
rowsol = new int[Nobj];
colsol = new int[Nobj];

const T *c_i, *n_j;
// calculate costs and fill the Munkres cost matrix
for (size_t i=0; i<Nobj; ++i){
  for (size_t j=0; j<Nobj; ++j){
    c_i = centre+i*span;
    n_j = neighbour+j*span;
    s_cost = R(0);
    if (Nel[0]){
      for (size_t k=0; k<Nel[0]; ++k)
        s_cost += magnitude(c_i[k] - n_j[k]);
      c_i += Nel[0];
      n_j += Nel[0];
    }
    if (Nel[1]){
      switch(vec_cost_func){
        case 0: v_cost = std::abs(std::sin(brille::utils::hermitian_angle(Nel[1], c_i, n_j))); break;
        case 1: v_cost = brille::utils::vector_distance(Nel[1], c_i, n_j); break;
        case 2: v_cost = 1-brille::utils::vector_product(Nel[1], c_i, n_j); break;
        case 3: v_cost =     brille::utils::vector_angle(Nel[1], c_i, n_j); break;
      }
      c_i += Nel[1];
      n_j += Nel[1];
    }
    if (Nel[2]){
      I nel2 = std::sqrt(Nel[2]);
      if (nel2*nel2 != Nel[2])
        throw std::runtime_error("Non-square matrix in jv_permutation");
      m_cost = brille::utils::frobenius_distance(nel2, c_i, n_j);
    }
    // for each i we want to determine the cheapest j
    // cost[i*Nobj+j] = std::log(Wscl*s_cost + Wvec*v_cost + Wmat*m_cost);
    cost[i*Nobj+j] = Wscl*s_cost + Wvec*v_cost + Wmat*m_cost;
  }
}

// use the Jonker-Volgenant algorithm to determine the optimal assignment
/*
There might be a hidden problem here.
As discussed in the README at https://github.com/hrldcpr/pyLAPJV

Supposedly, if two costs are equally smallest (in a row) to machine precision
then the Jonker-Volgenant algorithm enters an infinite loop.
The version in lapjv.h has a check to avoid this but it might still be a
problem.
*/
brille::assignment::lapjv((int)Nobj, cost, false, rowsol, colsol, usol, vsol);
/* use the fact that the neighbour objects have already had their global
   permutation saved into `permutations` to determine the global permuation
   for the centre objects too; storing the result into `permutations` as well.
*/
brille::ind_t nind = static_cast<brille::ind_t>(neighbour_idx);
brille::ind_t cind = static_cast<brille::ind_t>(centre_idx);
for (brille::ind_t i=0; i<Nobj; ++i)
  for (brille::ind_t j=0; j<Nobj; ++j)
    if (permutations.val(nind, i) == rowsol[j]) permutations.val(cind, i) = static_cast<size_t>(j);

delete[] cost;
delete[] usol;
delete[] vsol;
delete[] rowsol;
delete[] colsol;
return true;
}

template<class S, class T, class R, class I,
          typename=typename std::enable_if<std::is_same<typename CostTraits<T>::type, R>::value>::type
        >
bool jv_permutation(const S* centre_vals, const T* centre_vecs,
                    const S* neighbour_vals, const T* neighbour_vecs,
                    const std::array<I,3>& vals_Nel, const std::array<I,3>& vecs_Nel,
                    const R Wscl, const R Wvec, const R Wmat,
                    const size_t vals_span, const size_t vecs_span, const size_t Nobj,
                    bArray<size_t>& permutations,
                    const size_t centre_idx, const size_t neighbour_idx,
                    const int vec_cost_func = 0
                   ){
/* An earlier version of this function took `centre` and `neighbour` arrays
   which were effectively [span, Nobj] 2D arrays. This has the unfortunate
   property that individual objects were not contiguous in memory but instead
   had a stride equal to Nobj.
   Now this function requires arrays which are [Nobj, span], placing each
   object in contiguous memory and eliminating the need to copy sub-object
   vector/matrices into contiguous memory before calling subroutines.
*/
// initialize variables
R s_cost{0}, v_cost{0}, m_cost{0};
R *cost=nullptr, *usol=nullptr, *vsol=nullptr;
cost = new R[Nobj*Nobj];
usol = new R[Nobj];
vsol = new R[Nobj];
int *rowsol=nullptr, *colsol=nullptr;
rowsol = new int[Nobj];
colsol = new int[Nobj];

const S *vals_c_i, *vals_n_j;
const T *vecs_c_i, *vecs_n_j;
// calculate costs and fill the Munkres cost matrix
for (size_t i=0; i<Nobj; ++i){
  for (size_t j=0; j<Nobj; ++j){
    vals_c_i = centre_vals+i*vals_span;
    vals_n_j = neighbour_vals+j*vals_span;
    vecs_c_i = centre_vecs+i*vecs_span;
    vecs_n_j = neighbour_vecs+j*vecs_span;
    s_cost = R(0);
    if (vals_Nel[0]){
      for (size_t k=0; k<vals_Nel[0]; ++k)
        s_cost += magnitude(vals_c_i[k] - vals_n_j[k]);
      vals_c_i += vals_Nel[0];
      vals_n_j += vals_Nel[0];
    }
    if (vecs_Nel[0]){
      for (size_t k=0; k<vecs_Nel[0]; ++k)
        s_cost += magnitude(vecs_c_i[k] - vecs_n_j[k]);
      vecs_c_i += vecs_Nel[0];
      vecs_n_j += vecs_Nel[0];
    }
    v_cost = R(0);
    if (vals_Nel[1]){
      switch(vec_cost_func){
        case 0: v_cost += std::abs(std::sin(brille::utils::hermitian_angle(vals_Nel[1], vals_c_i, vals_n_j))); break;
        case 1: v_cost +=  brille::utils::vector_distance(vals_Nel[1], vals_c_i, vals_n_j); break;
        case 2: v_cost += 1-brille::utils::vector_product(vals_Nel[1], vals_c_i, vals_n_j); break;
        case 3: v_cost +=     brille::utils::vector_angle(vals_Nel[1], vals_c_i, vals_n_j); break;
      }
      vals_c_i += vals_Nel[1];
      vals_n_j += vals_Nel[1];
    }
    if (vecs_Nel[1]){
      switch(vec_cost_func){
        case 0: v_cost += std::abs(std::sin(brille::utils::hermitian_angle(vecs_Nel[1], vecs_c_i, vecs_n_j))); break;
        case 1: v_cost +=  brille::utils::vector_distance(vecs_Nel[1], vecs_c_i, vecs_n_j); break;
        case 2: v_cost += 1-brille::utils::vector_product(vecs_Nel[1], vecs_c_i, vecs_n_j); break;
        case 3: v_cost +=     brille::utils::vector_angle(vecs_Nel[1], vecs_c_i, vecs_n_j); break;
      }
      vecs_c_i += vecs_Nel[1];
      vecs_n_j += vecs_Nel[1];
    }
    m_cost = R(0);
    if (vals_Nel[2]){
      I nel2 = static_cast<I>(std::sqrt(vals_Nel[2]));
      if (nel2*nel2 != vals_Nel[2])
        throw std::runtime_error("Non-square matrix in jv_permutation");
      m_cost = brille::utils::frobenius_distance(nel2, vals_c_i, vals_n_j);
    }
    if (vecs_Nel[2]){
      I nel2 = static_cast<I>(std::sqrt(vecs_Nel[2]));
      if (nel2*nel2 != vecs_Nel[2])
        throw std::runtime_error("Non-square matrix in jv_permutation");
      m_cost = brille::utils::frobenius_distance(nel2, vecs_c_i, vecs_n_j);
    }
    // for each i we want to determine the cheapest j
    // cost[i*Nobj+j] = std::log(Wscl*s_cost + Wvec*v_cost + Wmat*m_cost);
    cost[i*Nobj+j] = Wscl*s_cost + Wvec*v_cost + Wmat*m_cost;
  }
}
// use the Jonker-Volgenant algorithm to determine the optimal assignment
/*
There might be a hidden problem here.
As discussed in the README at https://github.com/hrldcpr/pyLAPJV

Supposedly, if two costs are equally smallest (in a row) to machine precision
then the Jonker-Volgenant algorithm enters an infinite loop.
The version in lapjv.h has a check to avoid this but it might still be a
problem.
*/
brille::assignment::lapjv((int)Nobj, cost, false, rowsol, colsol, usol, vsol);
/* use the fact that the neighbour objects have already had their global
   permutation saved into `permutations` to determine the global permuation
   for the centre objects too; storing the result into `permutations` as well.
*/
brille::ind_t nind = static_cast<brille::ind_t>(neighbour_idx);
brille::ind_t cind = static_cast<brille::ind_t>(centre_idx);
for (brille::ind_t i=0; i<Nobj; ++i)
  for (brille::ind_t j=0; j<Nobj; ++j)
    if (permutations.val(nind, i) == rowsol[j]) permutations.val(cind, i) = static_cast<size_t>(j);

delete[] cost;
delete[] usol;
delete[] vsol;
delete[] rowsol;
delete[] colsol;
return true;
}

/*! \brief Use Junker-Volgenant algorithm to determine a permutation from a cost matrix

@param cost A square matrix held in a std::vector of the costs for each assignment
@returns The permutation as a std::vector<size_t>
*/
template<class I, class T>
std::vector<I> jv_permutation(const std::vector<T>& cost){
  I Nobj = static_cast<I>(std::sqrt(cost.size()));
  assert( cost.size() == Nobj*Nobj);
  std::vector<T> usol(Nobj,T(0)), vsol(Nobj,T(0));
  std::vector<I> rows(Nobj,I(0)), cols(Nobj,I(0));

  brille::assignment::lapjv(Nobj, cost.data(), false, rows.data(), cols.data(), usol.data(), vsol.data());
  return rows;
}
template<class T>
std::vector<int> jv_permutation(const std::vector<T>& cost){
  return jv_permutation<int,T>(cost);
}

template<class T, class I>
bool jv_permutation_fill(const std::vector<T>& cost, std::vector<I>& row){
  I Nobj = static_cast<I>(std::sqrt(cost.size()));
  assert(cost.size() == static_cast<size_t>(Nobj)*static_cast<size_t>(Nobj));
  std::vector<T> u(Nobj,T(0)), v(Nobj,T(0));
  std::vector<I> col(Nobj,0);
  row.resize(Nobj);
  brille::assignment::lapjv(Nobj, cost.data(), false, row.data(), col.data(), u.data(), v.data());
  return true;
}
template<class T, class I>
bool jv_permutation_fill(const std::vector<T>& cost, std::vector<I>& row, std::vector<I>& col){
  I Nobj = static_cast<I>(std::sqrt(cost.size()));
  assert(cost.size() == static_cast<size_t>(Nobj)*static_cast<size_t>(Nobj));
  std::vector<T> u(Nobj,T(0)), v(Nobj,T(0));
  row.resize(Nobj);
  col.resize(Nobj);
  brille::assignment::lapjv(Nobj, cost.data(), false, row.data(), col.data(), u.data(), v.data());
  return true;
}

/*! \brief Use a Stable Matching algorithm to determine a permutation

For two arrays of data, located in memory at `centre` and `neighbour`, and each
representing `Nobj` sets of  `Nel[0]` scalars, `Nel[1]` vector elements, and `Nel[2]` matrix elements,
determine the sorting permutation that maps the elements at `centre` onto the
same global mapping as the sorting permutation already stored at
`permutations[neighbour_idx]`. The resultant global permutation is stored at
`permutations[centre_idx]`.

Each array of data to be compared must be formatted as:

      [ 0{(0…Nel[0]-1)(0…Nel[1]-1)(0…Nel[2]-1)}
        1{(0…Nel[0]-1)(0…Nel[1]-1)(0…Nel[2]-1)}
        ⋮
        Nobj-1{(0…Nel[0]-1)(0…Nel[1]-1)(0…Nel[2]-1)} ]

The function constructs an `Nobj`×`Nobj` cost matrix, where each element is
given by

      Cᵢⱼ = Wscl×∑ₖ₋₀ᴺˢᶜˡ√(centre[i,k]-neighbour[j,k])
          + Wvec×𝔣ᵥ(centre_vec[i],neighbour_vec[j])
          + Wmat×𝔣ₘ(centre_mat[i],neighbour_mat[j])

where `Wscl`, `Wvec`, and `Wmat` are weight factors for adjusting the
relative cost of the scalar, vector, and matrix differences,
respectively; 𝔣ᵥ is the vector cost function; and 𝔣ₘ is the matrix cost function.

The vector cost function is selected by `vec_cost_func` to be one of:
the absolute value of the sine of the hermitian angle between eigenvectors,
the distance between eigenvectors, or one minus the inner product between
eigenvectors.

The matrix cost function is the Frobenius norm of the difference between the
two matrices.

@param centre The values to be sorted by the determined permutation
@param neighbour The values against which the `centre` values are compared
@param Nel[0] The number of scalars per object
@param Nel[1] The number of vector elements per object
@param Nel[2] The square root of the number of matrix elements per object
@param Wscl The cost weight for scalar elements
@param Wvec The cost weight for vectors
@param Wmat The cost weight for matrices
@param span The number of elements per object, `Nel[0]+Nel[1]+Nel[2]`
@param Nobj The number of objects at each of `centre` and `neighbour`
@param[out] permutations Contains the global permutation sorting for the
                         neighbour and is where the centre permutation is stored
@param centre_idx The index into `permutations` where the output is stored
@param neighbour_idx The index into `permutations` to find the neighbour permutation
@param vec_cost_func Used to select the eigenvector cost function:
                      0 --> `abs(sin(brille::utils::hermitian_angle))`
                      1 --> `vector_distance`
                      2 --> `1-brille::utils::vector_product`
                      3 --> `vector_angle`
@returns `true` if the permutation was assigned successfully, otherwise `false`
*/
template<class T, class R, class I,
          typename=typename std::enable_if<std::is_same<typename CostTraits<T>::type, R>::value>::type
        >
bool sm_permutation(const T* centre, const T* neighbour, const std::array<I,3>& Nel,
                    const R Wscl, const R Wvec, const R Wmat,
                    const size_t span, const size_t Nobj,
                    bArray<size_t>& permutations,
                    const size_t centre_idx, const size_t neighbour_idx,
                    const int vec_cost_func = 0
                   ){
/* An earlier version of this function took `centre` and `neighbour` arrays
   which were effectively [span, Nobj] 2D arrays. This has the unfortunate
   property that individual objects were not contiguous in memory but instead
   had a stride equal to Nobj.
   Now this function requires arrays which are [Nobj, span], placing each
   object in contiguous memory and eliminating the need to copy sub-object
   vector/matrices into contiguous memory before calling subroutines.
*/
// initialize variables
R s_cost{0}, v_cost{0}, m_cost{0};
R* cost = new R[Nobj*Nobj];
size_t* rowsol = new size_t[Nobj];
size_t* colsol = new size_t[Nobj];

const T *c_i, *n_j;
// calculate costs and fill the Munkres cost matrix
for (size_t i=0; i<Nobj; ++i){
  for (size_t j=0; j<Nobj; ++j){
    c_i = centre+i*span;
    n_j = neighbour+j*span;
    s_cost = R(0);
    if (Nel[0]){
      for (size_t k=0; k<Nel[0]; ++k)
        s_cost += magnitude(c_i[k] - n_j[k]);
      c_i += Nel[0];
      n_j += Nel[0];
    }
    if (Nel[1]){
      switch(vec_cost_func){
        case 0: v_cost = std::abs(std::sin(brille::utils::hermitian_angle(Nel[1], c_i, n_j))); break;
        case 1: v_cost =  brille::utils::vector_distance(Nel[1], c_i, n_j); break;
        case 2: v_cost = 1-brille::utils::vector_product(Nel[1], c_i, n_j); break;
        case 3: v_cost =     brille::utils::vector_angle(Nel[2], c_i, n_j); break;
        default: std::cout << "Unknown vector cost function. None used." << std::endl;
      }
      c_i += Nel[1];
      n_j += Nel[1];
    }
    if (Nel[2]){
      I nel2 = std::sqrt(Nel[2]);
      if (nel2*nel2 != Nel[2])
        throw std::runtime_error("Non-square matrix in sm_permutation");
      m_cost = brille::utils::frobenius_distance(nel2, c_i, n_j);
    }
    // for each i we want to determine the cheapest j
    cost[i*Nobj+j] = Wscl*s_cost + Wvec*v_cost + Wmat*m_cost;
  }
}

brille::assignment::smp(Nobj, cost, rowsol, colsol, false);
/* use the fact that the neighbour objects have already had their global
   permutation saved into `permutations` to determine the global permuation
   for the centre objects too; storing the result into `permutations` as well.
*/
brille::ind_t nind = static_cast<brille::ind_t>(neighbour_idx);
brille::ind_t cind = static_cast<brille::ind_t>(centre_idx);
for (brille::ind_t i=0; i<Nobj; ++i)
  for (brille::ind_t j=0; j<Nobj; ++j)
    if (permutations.val(nind, i) == rowsol[j]) permutations.val(cind, i) = static_cast<size_t>(j);


delete[] cost;
delete[] rowsol;
delete[] colsol;
return true;
}

// The following apply_permutation has been adapted from
//  https://devblogs.microsoft.com/oldnewthing/20170104-00/?p=95115
/*! \brief Apply a permutation to a random access iterable object

For a permutation held in one random access iterable object of length `N`
comprised of the integers `(0,N]` in some order, and a second iterable
containing `N` objects, find the permuted iterable such that

  permuted_object_iterable[i] = object_iterable[index_iterable[i]]

for all `i` by swapping elements in place and using the permutation iterable
as a scratch workspace.

Upon completion the permutation iterable will be ordered 0:N-1.

@param objects An iterator pointing to the first element of the object iterable
@param end An iterator pointing to the end of the object iterable
@param indices An iterator pointing to the first element of the permutation iterable
*/
template<typename ObjItr, typename PermItr>
void
apply_permutation(ObjItr objects, ObjItr end, PermItr indices){
  using Obj = typename std::iterator_traits<ObjItr>::value_type;
  using Dif = typename std::iterator_traits<PermItr>::value_type;
  Dif numel = end - objects;
  for (Dif i=0; i<numel; ++i) if (i != indices[i]) {
    // move the object to a temporary location (fallsback to copy)
    Obj obji{std::move(objects[i])};
    // keep track of where we are in the swap loop
    Dif current = i;
    while (i != indices[current]){
      Dif next = indices[current];
      objects[current] = std::move(objects[next]);
      indices[current] = current;
      current = next;
    }
    objects[current] = std::move(obji);
    indices[current] = current;
  }
}

/*! \brief Apply an inverse permutation to a random access iterable object

For an inverse permutation held in one random access iterable object of length
`N` comprised of the integers `(0,N]` in some order, and a second iterable
containing `N` objects, find the inverse permuted iterable such that

  permuted_object_iterable[index_iterable[i]] = object_iterable[i]

for all `i` by swapping elements in place and using the inverse permutation
iterable as a scratch workspace.

@param objects An iterator pointing to the first element of the object iterable
@param end An iterator pointing to the end of the object iterable
@param indices An iterator pointing to the first element of the inverse permutation iterable
*/
template<typename ObjItr, typename PermItr>
void
apply_inverse_permutation(ObjItr objects, ObjItr end, PermItr indices){
  using Obj = typename std::iterator_traits<ObjItr>::value_type;
  using Dif = typename std::iterator_traits<PermItr>::value_type;
  Dif numel = end - objects;
  for (Dif i=0; i<numel; ++i) while (i != indices[i]) {
    // pop-out the targeted object and its index
    Dif pop_idx = indices[indices[i]];
    Obj pop_obj{std::move(objects[indices[i]])};
    // move the current object and index to the target
    objects[indices[i]] = std::move(objects[i]);
    indices[indices[i]] = indices[i];
    // and put the popped object and index back in here
    objects[i] = std::move(pop_idx);
    indices[i] = pop_idx;
  }
}

} // namespace brille
#endif
