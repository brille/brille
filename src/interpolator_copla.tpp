/* This file is part of brille.

Copyright © 2021 Greg Tucker <gregory.tucker@ess.eu>

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

/* A note about the use of the vectors in a GammaTable object:

For speed purposes the GammaTable vectors do not contain lattice information
when accessed by .vector, as below. There is another accessor .ldvector which
adds the lattice, but this significantly slows down the dot product due to
equivalent/inverse lattice checks reconstructing the Direct/Reciprocal objects
multiple times *for every q point*(!).

The per-q-point check can be avoided as long as we check *before* using the
GammTable vectors through their .vector accessor. Then we know that the dot
product of q = (h,k,l) and d = (a,b,c) is q⋅d = 2π (h*a + k*b + l*c).
*/

// Implement in-place 'rotation' of eigenvectors of Copla's quadratic boson Hamiltonian
// There is (hopefully) only a label permutation without a phase-factor
// For now, we can make use of a GammaTable for the permutations and live with
// the extraneous displacement vectors being tabulated as well
template<class T>
template<class R>
bool Interpolator<T>::rip_copla_complex(
  bArray<T>& x, const LQVec<R>& q, const GammaTable& pgt,
  const PointSymmetry&, const std::vector<size_t>& ridx, const std::vector<size_t>& invRidx,
  const int nthreads
) const {
  profile_update("Start Interpolator::rip_copla_complex method");
  if (! pgt.lattice().isstar(q.get_lattice()))
    throw std::runtime_error("The q points and GammaTable must be in mutually reciprocal lattices");
  verbose_update("Interpolator::rip_copla_complex called with ",nthreads," threads");
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  // Copla's matrix is 2Nx2N for a system with N atoms, so its eigenvectors
  // must have 2N elements
  auto one_n = pgt.atom_count();
  if (_elements[1] != 2*one_n){
    std::cout << "Copla matrix eigenvector transformation requires 2N elements" << std::endl;
    return false;
  }
  // Copla's method does not support tensors, maybe it can be extended.
  if (_elements[2] > 0){
    std::cout << "Copla's method does not support tensors" << std::endl;
    return false;
  }
  const ind_t b_{this->branches()}, s_{this->branch_span()};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
  // Start the parallel region so we caninitialize per-thread temporary memory
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
  #pragma omp parallel default(none) shared(x,q,pgt,ridx,invRidx,one_n,xsize)
#else
  #pragma omp parallel default(none) shared(b_,s_,x,q,pgt,ridx,invRidx,one_n,xsize)
#endif
  {
  std::vector<T> tA(2*one_n);
  #pragma omp for schedule(static)
  for (long long si=0; si<xsize; ++si){
    ind_t i = brille::utils::s2u<ind_t, long long>(si);
    T * xi = x.ptr(i);
    for (ind_t b=0; b<b_; ++b){
      // scalar elements do not need to be rotated, so skip them
      ind_t o = b*s_ + _elements[0];
      // the *full 2N* eigenvector for mode j transforms as
      //  ϵʲ(Rq)  =  Γ(q;R) eʲ(q)
      // where Γ(q;R) is a 2N×2N matrix with elements given by
      // the submatrices Γₖᵥᵅᵝ(q;R) = Rᵅᵝ δ(v,F₀⁻¹(k,R))
      // (without a phase factor, like exp{i q⋅[R⁻¹xₖ - xᵥ]} )
      // The GammaTable contains predetermined F₀(k,R) and R⁻¹xₖ - xᵥ
      // to make this calculation more straightforward.
      //
      // Permute the atom indexes for the eigenvector, which is doubled
      for (ind_t k=0; k<one_n; ++k){
        auto v0idx = pgt.F0(k, invRidx[i]);
        tA[v0idx] = xi[o+k];
        tA[v0idx+one_n] = xi[o+k+one_n];
      }
      for (ind_t k=0; k<2*one_n; ++k) xi[o+k] = tA[k];
      //
      // if (_elements[2]>0){
      //   // do something, maybe follow interpolator_gamma
      // }
    }
  }
  }
  profile_update("  End Interpolator::rip_copla_complex method");
  return true;
}
