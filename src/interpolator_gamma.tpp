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
  data_t<T>& x, const LQVec<R>& q, const GammaTable& pgt,
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
  const ind_t b_{this->branches()};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
#pragma omp parallel for default(none) shared(x,q,pgt,ptsym,ridx,invRidx,e_iqd_gt) firstprivate(no,Nmat,xsize) schedule(static)
#else
#pragma omp parallel for default(none) shared(b_,x,q,pgt,ptsym,ridx,invRidx,e_iqd_gt) firstprivate(no,Nmat,xsize) schedule(static)
#endif
  for (long long si=0; si<xsize; ++si){
    ind_t i = brille::utils::s2u<ind_t, long long>(si);
    auto xi = x.slice(i); // (1,), (B,), (Y,) or (B,Y)
    size_t Rii = ridx[i], iRii = invRidx[i];
    for (ind_t b=0; b<b_; ++b){
      auto xib = (b_>1 ? xi.slice(b) : xi).contiguous_copy(); // (1,) or (Y,)
      // scalar elements do not need to be rotated, so skip them
      ind_t o = no[0];
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
        T tmp_v[3];
        std::vector<T> tmpvecs(no[1]*3u);
        for (ind_t k=0; k<no[1]; ++k){
          // is q already expressed in the right lattice? hopefully!
          // if the vector rotates by Γ is *must* be complex, so T is a complex type
          // use its constructor to make i q⋅[R⁻¹xₖ - xᵥ] and calculate the (k,R) phase
          brille::utils::mul_mat_vec(tmp_v, 3u, ptsym.get(iRii).data(), xib.ptr()+o);
          auto v0idx = 3u*pgt.F0(k, iRii);
          T phase = e_iqd_gt(i, k, iRii);
          for (int j=0; j<3; ++j) tmpvecs[v0idx+j] = phase*tmp_v[j];
          o += 3u;
        }
        for (ind_t j=0; j<no[1]*3u; ++j) xib.val(j) = tmpvecs[j];
      }
      if (no[2]>0){
        T tmp_m[9], tmp_m2[9];
        std::vector<T> tmpmats(no[2]*9u);
        for (ind_t n=0; n<Nmat; ++n){
          T Rph = e_iqd_gt(i, n, Rii);
          ind_t v = static_cast<ind_t>(pgt.F0(n, Rii));
          for (ind_t m=0; m<Nmat; ++m){
            T iRph = e_iqd_gt(i, m, iRii);
            ind_t k = static_cast<ind_t>(pgt.F0(m, iRii));
            // Calculate R⁻¹*M*R in two steps
            // first calculate M*R, storing in tmp_m
            brille::utils::mul_mat_mat(tmp_m, 3u, xib.ptr()+o+(n*Nmat+m)*9u, ptsym.get(Rii).data());
            // next calculate R⁻¹*tmp_m, storing in the temporary all matrix array
            brille::utils::mul_mat_mat(tmp_m2, 3u, ptsym.get(iRii).data(), tmp_m);
            // include the R R⁻¹ phase factor
            for (int j=0; j<9; ++j) tmpmats[(v*Nmat+k)*9u+j] = Rph*iRph*tmp_m2[j];
          }
        }
        for (ind_t j=0; j<no[2]*9u; ++j) xib.val(j) = tmpmats[j];
      }
      // we may have made a copy when we produced the contiguous xib
      if (!xi.shares_with(xib)){
        for (ind_t z=0; z<xib.size(0); ++z) xi.val(b,z) = xib.val(z);
      }
    }
  }
  profile_update("  End Interpolator::rip_gamma_complex method");
  return true;
}

template<class T>
template<class R>
bool Interpolator<T>::rip_gamma_vec_complex(
  data_t<T>& x, const LQVec<R>& q, const GammaTable& pgt,
  const PointSymmetry& ptsym, const std::vector<size_t>&, const std::vector<size_t>& invRidx,
  const int nthreads
) const {
  profile_update("Start Interpolator::rip_gamma_vec_complex method");
  // construct a lambda to calculate the phase for qᵢ, and [R⁻¹xₖ - xᵥ]
  auto e_iqd_gt = [q,pgt](ind_t i, ind_t k, size_t r){
    return e_iqd(q, i, pgt.vectors(), pgt.vector_index(k,r));
  };
  if (! pgt.lattice().isstar(q.get_lattice()))
    throw std::runtime_error("The q points and GammaTable must be in mutually reciprocal lattices");
  verbose_update("Interpolator::rip_gamma_vec_complex called with ",nthreads," threads");
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  // for Γ transformations of tensors, there *must* be N×N in total:
  const ind_t b_{data_.size(1)}, v_{data_.size(2)};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
  T tmp_v[3]; // allocated private variable possible as it's statically sized
  std::vector<T> tmpvecs(v_*3u); // private possible but will not be pre-allocated in workers unless *first*private
  #if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
  #pragma omp parallel for default(none) shared(x,q,pgt,ptsym,invRidx,e_iqd_gt) firstprivate(xsize,tmp_v,tmpvecs) schedule(static)
  #else
  #pragma omp parallel for default(none) shared(b_,v_,x,q,pgt,ptsym,invRidx,e_iqd_gt) firstprivate(xsize,tmp_v,tmpvecs) schedule(static)
  #endif
  for (long long si=0; si<xsize; ++si){
    ind_t i = brille::utils::s2u<ind_t, long long>(si);
    auto xi = x.slice(i); // (B,V,3)
    // size_t Rii = ridx[i], iRii = invRidx[i]; // ridx[i] is only used with the matrix form:
    auto iRii = invRidx[i];
    for (ind_t b=0; b<b_; ++b){
      auto xib = xi.slice(b); // (V,3)
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
      for (ind_t k=0; k<v_; ++k){
        auto xibk = xib.slice(k); // required to be contiguous (3,)
        // is q already expressed in the right lattice? hopefully!
        // if the vector rotates by Γ is *must* be complex, so T is a complex type
        // use its constructor to make i q⋅[R⁻¹xₖ - xᵥ] and calculate the (k,R) phase
        brille::utils::mul_mat_vec(tmp_v, 3u, ptsym.get(iRii).data(), xibk.ptr());
        auto v0idx = 3u*pgt.F0(k, iRii);
        T phase = e_iqd_gt(i, k, iRii);
        for (int j=0; j<3; ++j) tmpvecs[v0idx+j] = phase*tmp_v[j];
      }
      for (ind_t j=0; j<v_; ++j)
        for(ind_t k=0; k<3u; ++k)
          xib.val(j,k) = tmpvecs[j*3u + k];
    }
  }
  profile_update("  End Interpolator::rip_gamma_vec_complex method");
  return true;
}


template<class T>
template<class R>
bool Interpolator<T>::rip_gamma_mat_complex(
  data_t<T>& x, const LQVec<R>& q, const GammaTable& pgt,
  const PointSymmetry& ptsym, const std::vector<size_t>& ridx, const std::vector<size_t>& invRidx,
  const int nthreads
) const {
  profile_update("Start Interpolator::rip_gamma_mat_complex method");
  // construct a lambda to calculate the phase for qᵢ, and [R⁻¹xₖ - xᵥ]
  auto e_iqd_gt = [q,pgt](ind_t i, ind_t k, size_t r){
    return e_iqd(q, i, pgt.vectors(), pgt.vector_index(k,r));
  };
  if (! pgt.lattice().isstar(q.get_lattice()))
    throw std::runtime_error("The q points and GammaTable must be in mutually reciprocal lattices");
  verbose_update("Interpolator::rip_gamma_mat_complex called with ",nthreads," threads");
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  const ind_t b_{data_.size(1)}, m_{data_.size(2)};
  // for Γ transformations of tensors, there *must* be N×N in total:
  ind_t Nmat = static_cast<ind_t>(std::sqrt(m_));
  if (m_ != Nmat*Nmat){
    std::cout << "Atomic displacement Gamma transformation requires NxN 3x3 tensors!" << std::endl;
    return false;
  }
  T tmp_m[9], tmp_m2[9];
  std::vector<T> tmpmats(m_*9u); // these won't be allocated in the different threads :(
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
#pragma omp parallel for default(none) shared(x,q,pgt,ptsym,ridx,invRidx,e_iqd_gt) firstprivate(Nmat,xsize,tmp_m,tmp_m2,tmpmats) schedule(static)
#else
#pragma omp parallel for default(none) shared(b_,m_,x,q,pgt,ptsym,ridx,invRidx,e_iqd_gt) firstprivate(Nmat,xsize,tmp_m,tmp_m2,tmpmats) schedule(static)
#endif
  for (long long si=0; si<xsize; ++si){
    ind_t i = brille::utils::s2u<ind_t, long long>(si);
    auto xi = x.slice(i); // (B,M,3,3)
    size_t Rii = ridx[i], iRii = invRidx[i];
    for (ind_t b=0; b<b_; ++b){
      auto xib = xi.slice(b); // (M,3,3) with each (3,3) matrix required to be row-ordered contiguous
      for (ind_t n=0; n<Nmat; ++n){
        T Rph = e_iqd_gt(i, n, Rii);
        ind_t v = static_cast<ind_t>(pgt.F0(n, Rii));
        for (ind_t m=0; m<Nmat; ++m){
          T iRph = e_iqd_gt(i, m, iRii);
          ind_t k = static_cast<ind_t>(pgt.F0(m, iRii));
          // Calculate R⁻¹*M*R in two steps
          // first calculate M*R, storing in tmp_m
          brille::utils::mul_mat_mat(tmp_m, 3u, xib.ptr(n*Nmat+m), ptsym.get(Rii).data());
          // next calculate R⁻¹*tmp_m, storing in the temporary all matrix array
          brille::utils::mul_mat_mat(tmp_m2, 3u, ptsym.get(iRii).data(), tmp_m);
          // include the R R⁻¹ phase factor
          for (int j=0; j<9; ++j) tmpmats[(v*Nmat+k)*9u+j] = Rph*iRph*tmp_m2[j];
        }
      }
      for (ind_t m=0; m<m_; ++m){
        // loop over the m_ matrices since we only require that the (3,3) part
        // be contigous and row-ordered.
        T* pibm = xib.ptr(m);
        for (ind_t z=0; z<9u; ++z) pibm[z] = tmpmats[m*9u+z];
      }
    }
  }
  profile_update("  End Interpolator::rip_gamma_mat_complex method");
  return true;
}
