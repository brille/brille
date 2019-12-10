/* Copyright 2019 Greg Tucker
//
// This file is part of brille.
//
// brille is free software: you can redistribute it and/or modify it under the
// terms of the GNU Affero General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// brille is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with brille. If not, see <https://www.gnu.org/licenses/>.            */

/* This file holds implementations of Linear Assignment Problem solvers which
   can be used to find the sorting permutations at grid points.              */
#ifndef _PERMUTATION_H
#define _PERMUTATION_H

#include "munkres.hpp"
#include "lapjv.hpp"
#include "smp.hpp"

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
@param Nel The number of scalars, eigenvector elements, vector elements, and matrix elements per object
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
                      2 --> `1-vector_product`
                      3 --> `abs(sin(hermitian_angle))`
@returns `true` if the permutation was assigned successfully, otherwise `false`
*/
template<class T, class R, class I,
          typename=typename std::enable_if<std::is_same<typename CostTraits<T>::type, R>::value>::type
        >
bool munkres_permutation(const T* centre, const T* neighbour, const std::array<I,3>& Nel,
                         const R Wscl, const R Wvec, const R Wmat,
                         const size_t span, const size_t Nobj,
                         ArrayVector<size_t>& permutations,
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
Munkres<R> munkres(Nobj);
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
        case 0: v_cost = std::abs(std::sin(hermitian_angle(Nel[1], c_i, n_j))); break;
        case 1: v_cost = vector_distance(Nel[1], c_i, n_j); break;
        case 2: v_cost = 1-vector_product(Nel[1], c_i, n_j); break;
        case 3: v_cost = vector_angle(Nel[1], c_i, n_j); break;
      }
      c_i += Nel[1];
      n_j += Nel[1];
    }
    if (Nel[2]){
      m_cost = frobenius_distance(Nel[2], c_i, n_j);
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
for (size_t i=0; i<Nobj; ++i)
  for (size_t j=0; j<Nobj; ++j)
    if (permutations.getvalue(neighbour_idx,i)==assignment[j])
      permutations.insert(j,centre_idx,i);
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
                      0 --> `abs(sin(hermitian_angle))`
                      1 --> `vector_distance`
                      2 --> `1-vector_product`
                      3 --> `vector_angle`
@returns `true` if the permutation was assigned successfully, otherwise `false`
*/
template<class T, class R, class I,
          typename=typename std::enable_if<std::is_same<typename CostTraits<T>::type, R>::value>::type
        >
bool jv_permutation(const T* centre, const T* neighbour, const std::array<I,3>& Nel,
                    const R Wscl, const R Wvec, const R Wmat,
                    const size_t span, const size_t Nobj,
                    ArrayVector<size_t>& permutations,
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
        case 0: v_cost = std::abs(std::sin(hermitian_angle(Nel[1], c_i, n_j))); break;
        case 1: v_cost = vector_distance(Nel[1], c_i, n_j); break;
        case 2: v_cost = 1-vector_product(Nel[1], c_i, n_j); break;
        case 3: v_cost =     vector_angle(Nel[1], c_i, n_j); break;
      }
      c_i += Nel[1];
      n_j += Nel[1];
    }
    if (Nel[2]){
      m_cost = frobenius_distance(Nel[2], c_i, n_j);
    }
    // for each i we want to determine the cheapest j
    cost[i*Nobj+j] = std::log(Wscl*s_cost + Wvec*v_cost + Wmat*m_cost);
  }
}
// use the Jonker-Volgenant algorithm to determine the optimal assignment
/*
There might be a hidden problem here.

Supposedly, if two costs are equally smallest (in a row) to machine precision
then the Jonker-Volgenant algorithm enters an infinite loop.
The version in lapjv.h has a check to avoid this but it might still be a
problem.
*/
lapjv((int)Nobj, cost, false, rowsol, colsol, usol, vsol);
/* use the fact that the neighbour objects have already had their global
   permutation saved into `permutations` to determine the global permuation
   for the centre objects too; storing the result into `permutations` as well.
*/
for (size_t i=0; i<Nobj; ++i)
  for (size_t j=0; j<Nobj; ++j)
    if (permutations.getvalue(neighbour_idx,i)==static_cast<size_t>(rowsol[j])) // or should this be colsol?
      permutations.insert(j,centre_idx,i);

delete[] cost;
delete[] usol;
delete[] vsol;
delete[] rowsol;
delete[] colsol;
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
                      0 --> `abs(sin(hermitian_angle))`
                      1 --> `vector_distance`
                      2 --> `1-vector_product`
                      3 --> `vector_angle`
@returns `true` if the permutation was assigned successfully, otherwise `false`
*/
template<class T, class R, class I,
          typename=typename std::enable_if<std::is_same<typename CostTraits<T>::type, R>::value>::type
        >
bool sm_permutation(const T* centre, const T* neighbour, const std::array<I,3>& Nel,
                    const R Wscl, const R Wvec, const R Wmat,
                    const size_t span, const size_t Nobj,
                    ArrayVector<size_t>& permutations,
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
        case 0: v_cost = std::abs(std::sin(hermitian_angle(Nel[1], c_i, n_j))); break;
        case 1: v_cost =  vector_distance(Nel[1], c_i, n_j); break;
        case 2: v_cost = 1-vector_product(Nel[1], c_i, n_j); break;
        case 3: v_cost =     vector_angle(Nel[2], c_i, n_j); break;
        default: std::cout << "Unknown vector cost function. None used." << std::endl;
      }
      c_i += Nel[1];
      n_j += Nel[1];
    }
    if (Nel[2]){
      m_cost = frobenius_distance(Nel[2], c_i, n_j);
    }
    // for each i we want to determine the cheapest j
    cost[i*Nobj+j] = Wscl*s_cost + Wvec*v_cost + Wmat*m_cost;
  }
}

smp(Nobj, cost, rowsol, colsol, false);
/* use the fact that the neighbour objects have already had their global
   permutation saved into `permutations` to determine the global permuation
   for the centre objects too; storing the result into `permutations` as well.
*/
for (size_t i=0; i<Nobj; ++i)
  for (size_t j=0; j<Nobj; ++j)
    if (permutations.getvalue(neighbour_idx,i)==rowsol[j]) // or should this be colsol?
      permutations.insert(j,centre_idx,i);

delete[] cost;
delete[] rowsol;
delete[] colsol;
return true;
}

#endif
