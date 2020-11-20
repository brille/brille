/* This file is part of brille.

Copyright © 2020 Greg Tucker <greg.tucker@stfc.ac.uk>

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

template<class T, class R,
class S = typename std::common_type<T,R>::type,
class I = typename brille::ind_t
>
static std::complex<S>
e_iqd(const LQVec<T>& q, const I i, const bArray<R>& d, const size_t j)
{
  // const double pi = 3.14159265358979323846;
  S dotqd{0};
  I jj = static_cast<I>(j);
  for (I k=0; k<3u; ++k)
    dotqd += q.val(i,k)*d.val(jj,k);
  std::complex<S> i_2pi_x(0, 2*brille::pi*dotqd);
  return std::exp(i_2pi_x);
}

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

template<class T>
template<class R>
bool Interpolator<T>::rip_gamma_complex(
  bArray<T>& x, const LQVec<R>& q, const GammaTable& pgt,
  const PointSymmetry& ptsym, const std::vector<size_t>& ridx, const std::vector<size_t>& invRidx,
  const int nthreads
) const {
  profile_update("Start Interpolator::rip_gamma_complex method");
  // construct a lambda to calculate the phase for qᵢ, and [R⁻¹xₖ - xᵥ]
  auto e_iqd_gt = [q,pgt](ind_t i, ind_t k, size_t r){
    return e_iqd(q, i, pgt.vectors(), pgt.vector_index(k,r));
  };
  if (! pgt.lattice().isstar(q.get_lattice()))
    throw std::runtime_error("The q points and GammaTable must be in mutually reciprocal lattices");
  verbose_update("Interpolator::rip_gamma_complex called with ",nthreads," threads");
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  auto no = this->count_scalars_vectors_matrices();
  if (!std::any_of(no.begin()+1, no.end(), [](ind_t n){return n>0;}))
    return false;
  // for Γ transformations of tensors, there *must* be N×N in total:
  ind_t Nmat = static_cast<ind_t>(std::sqrt(no[2]))/3;
  if (no[2] != 9*Nmat*Nmat){
    std::cout << "Atomic displacement Gamma transformation requires NxN 3x3 tensors!" << std::endl;
    return false;
  }
  const ind_t b_{this->branches()}, s_{this->branch_span()};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
  std::vector<T> tA;
  T t0[9], t1[9];
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
#pragma omp parallel for default(none) shared(x,q,pgt,ptsym,ridx,invRidx,e_iqd_gt) private(t0,t1,tA) firstprivate(no,Nmat,xsize) schedule(static)
#else
#pragma omp parallel for default(none) shared(b_,s_,x,q,pgt,ptsym,ridx,invRidx,e_iqd_gt) private(t0,t1,tA) firstprivate(no,Nmat,xsize) schedule(static)
#endif
  for (long long si=0; si<xsize; ++si){
    ind_t i = brille::utils::s2u<ind_t, long long>(si);
    T * xi = x.ptr(i);
    size_t Rii = ridx[i], iRii = invRidx[i];
    for (ind_t b=0; b<b_; ++b){
      // scalar elements do not need to be rotated, so skip them
      ind_t o = b*s_ + no[0];
      // the *full 3N* eigenvector for mode j transforms as
      //  ϵʲ(Rq)  =  Γ(q;R) eʲ(q)
      // where Γ(q;R) is a 3N×3N matrix with elements given by
      // the submatrices Γₖᵥᵅᵝ(q;R) = Rᵅᵝ δ(v,F₀⁻¹(k,R)) exp{i q⋅[R⁻¹xₖ - xᵥ]}
      // -- I *think* that for a given k and R ∑ᵥ δ(v,F₀⁻¹(k,R)) ≡ 1
      //    that is, all equivalent-atom mappings are singular. THIS SHOULD BE VERIFIED!
      // The GammaTable contains predetermined F₀(k,R) and R⁻¹xₖ - xᵥ
      // to make this calculation more straightforward.
      //
      // Rotate, permute, and apply the phase factor simultaneously into a temporary array
      if (no[1]>0){
        ind_t o0 = o;
        if (tA.size() < no[1]*3u) tA.resize(no[1]*3u);
        for (ind_t k=0; k<no[1]; ++k){
          // is q already expressed in the right lattice? hopefully!
          // if the vector rotates by Γ is *must* be complex, so T is a complex type
          // use its constructor to make i q⋅[R⁻¹xₖ - xᵥ] and calculate the (k,R) phase
          brille::utils::mul_mat_vec(t0, 3u, ptsym.get(iRii).data(), xi+o);
          auto v0idx = 3u*pgt.F0(k, iRii); // ×3u to account for stride
          T phase = e_iqd_gt(i, k, iRii);
          for (int j=0; j<3; ++j) tA[v0idx+j] = phase*t0[j];
          o += 3u; // increase the offset for each vector
        }
        for (ind_t j=0; j<no[1]*3u; ++j) xi[o0+j] = tA[j];
      }
      if (no[2]>0){
        if (tA.size() < no[2]*9u) tA.resize(no[2]*9u);
        for (ind_t n=0; n<Nmat; ++n){
          T Rph = e_iqd_gt(i, n, Rii);
          ind_t v = static_cast<ind_t>(pgt.F0(n, Rii));
          for (ind_t m=0; m<Nmat; ++m){
            T iRph = e_iqd_gt(i, m, iRii);
            ind_t k = static_cast<ind_t>(pgt.F0(m, iRii));
            // Calculate R⁻¹*M*R in two steps
            // first calculate M*R, storing in t0
            brille::utils::mul_mat_mat(t0, 3u, xi+o+9u*(n*Nmat+m), ptsym.get(Rii).data());
            // next calculate R⁻¹*t0, storing in the temporary all matrix array
            brille::utils::mul_mat_mat(t1, 3u, ptsym.get(iRii).data(), t0);
            // include the R R⁻¹ phase factor
            for (int j=0; j<9; ++j) tA[(v*Nmat+k)*9u+j] = Rph*iRph*t1[j];
          }
        }
        for (ind_t j=0; j<no[2]*9u; ++j) xi[o+j] = tA[j];
      }
    }
  }
  profile_update("  End Interpolator::rip_gamma_complex method");
  return true;
}
