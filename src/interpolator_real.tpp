template<class T>
bool Interpolator<T>::rip_real(
  data_t<T>& x, const PointSymmetry& ptsym, const std::vector<size_t>& r, const std::vector<size_t>& invR, const int nthreads
) const {
  profile_update("Start Interpolator::rip_real method");
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  auto no = this->count_scalars_vectors_matrices();
  if (!std::any_of(no.begin()+1, no.end(), [](ind_t n){return n>0;}))
    return false;
  T tmp_v[3], tmp_m[9];
  const ind_t b_{this->branches()};
  std::array<int,9> ident = {1,0,0, 0,1,0, 0,0,1};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
// otherwise gcc complains that const b_ is *already* shared
#pragma omp parallel for default(none) shared(x,ptsym,r,invR) private(tmp_v,tmp_m) firstprivate(ident,no,xsize) schedule(static)
#else
#pragma omp parallel for default(none) shared(b_,x,ptsym,r,invR) private(tmp_v,tmp_m) firstprivate(ident,no,xsize) schedule(static)
#endif
  for (long long si=0; si<xsize; ++si){
    ind_t i = brille::utils::s2u<ind_t, long long>(si);
    auto xi = x.slice(i); // (1,), (B,), (Y,), (B, Y)
    if (!brille::approx::matrix(3, ident.data(), ptsym.get(r[i]).data())){
      for (ind_t b=0; b<b_; ++b){
        auto xib = (b_>1 ? xi.slice(b) : xi).contiguous_copy();  // (1,) or (Y,)
        // scalar elements do not need to be rotated, so skip them
        ind_t o = no[0];
        // rotate real vectors: since Q = Rᵀq + τ → Rv
        for (ind_t v=0; v<no[1]; ++v){
          T *xptr = xib.ptr(o);
          brille::utils::mul_mat_vec(tmp_v, 3u, ptsym.get(r[i]).data(), xptr);
          for (int j=0; j<3; ++j) xptr[j] = tmp_v[j];
          o += 3u; // shift 3 for each vector
        }
        for (ind_t m=0; m<no[2]; ++m){
          T *xptr = xib.ptr(o);
          // Calculate R*M*R⁻¹ in two steps
          // first calculate M*R⁻¹, storing in tmp_m
          brille::utils::mul_mat_mat(tmp_m, 3u, xptr, ptsym.get(invR[i]).data());
          // next calculate R*tmp_m, storing back in the x array
          brille::utils::mul_mat_mat(xptr, 3u, ptsym.get(r[i]).data(), tmp_m);
          o += 9u; // shift 9 for each matrix
        }
        // we may have made a copy when we produced the contiguous xib
        if (!xi.shares_with(xib)){
          for (ind_t z=0; z<xib.size(0); ++z) xi.val(b,z) = xib.val(z);
        }
      }
    }
  }
  profile_update("  End Interpolator::rip_real method");
  return true;
}

template<class T>
bool
Interpolator<T>::rip_real_vec(data_t<T>& x, const PointSymmetry& ptsym, const std::vector<size_t>& r, const std::vector<size_t>& invR, const int nthreads) const {
  profile_update("Start Interpolator::rip_real_vec method");
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  if (data_.ndim() != 4) return false;
  T tmp_v[3];
  const ind_t b_{data_.size(1)}, v_{data_.size(2)};
  std::array<int,9> ident = {1,0,0, 0,1,0, 0,0,1};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
// otherwise gcc complains that const b_ is *already* shared
#pragma omp parallel for default(none) shared(x,ptsym,r,invR) private(tmp_v) firstprivate(ident,xsize) schedule(static)
#else
#pragma omp parallel for default(none) shared(b_,v_,x,ptsym,r,invR) private(tmp_v) firstprivate(ident,xsize) schedule(static)
#endif
  for (long long si=0; si<xsize; ++si){
    ind_t i = brille::utils::s2u<ind_t, long long>(si);
    auto xi = x.slice(i); // (B, V, 3)
    if (!brille::approx::matrix(3, ident.data(), ptsym.get(r[i]).data())){
      for (ind_t b=0; b<b_; ++b){
        auto xib = xi.slice(b);  // (V,3)
        // we chose R such that Q = Rᵀq + τ, but are rotating real vectors
        for (ind_t v=0; v<v_; ++v){
          auto xibv = xib.slice(b); // (3,) *required* to be contiguous!
          brille::utils::mul_mat_vec(tmp_v, 3u, ptsym.get(r[i]).data(), xibv.ptr());
          for (ind_t j=0; j<3u; ++j) xibv.val(j) = tmp_v[j];
        }
      }
    }
  }
  profile_update("  End Interpolator::rip_real_vec method");
  return true;
}

template<class T>
bool
Interpolator<T>::rip_real_mat(data_t<T>& x, const PointSymmetry& ptsym, const std::vector<size_t>& r, const std::vector<size_t>& invR, const int nthreads) const {
  profile_update("Start Interpolator::rip_real_mat method");
  omp_set_num_threads( (nthreads>0) ? nthreads : omp_get_max_threads() );
  if (data_.ndim() != 5u) return false;
  T tmp_m[9];
  const ind_t b_{data_.size(1)}, m_{data_.size({2})};
  std::array<int,9> ident = {1,0,0, 0,1,0, 0,0,1};
  // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
  long long xsize = brille::utils::u2s<long long, ind_t>(x.size(0));
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
// otherwise gcc complains that const b_ is *already* shared
#pragma omp parallel for default(none) shared(x,ptsym,r,invR) private(tmp_m) firstprivate(ident,xsize) schedule(static)
#else
#pragma omp parallel for default(none) shared(b_,m_,x,ptsym,r,invR) private(tmp_m) firstprivate(ident,xsize) schedule(static)
#endif
  for (long long si=0; si<xsize; ++si){
    ind_t i = brille::utils::s2u<ind_t, long long>(si);
    auto xi = x.slice(i); // (B, M, 3, 3)
    if (!brille::approx::matrix(3, ident.data(), ptsym.get(r[i]).data())){
      for (ind_t b=0; b<b_; ++b){
        auto xib = xi.slice(b);  // (M, 3, 3)
        for (ind_t m=0; m<m_; ++m){
          auto xibm = xib.slice(m); // (3,3) *required* to be row-ordered contiguous
          // Calculate R*M*R⁻¹ in two steps
          // first calculate M*R⁻¹, storing in tmp_m
          brille::utils::mul_mat_mat(tmp_m, 3u, xibm.ptr(), ptsym.get(invR[i]).data());
          // next calculate R*tmp_m, storing back in the x array
          brille::utils::mul_mat_mat(xibm.ptr(), 3u, ptsym.get(r[i]).data(), tmp_m);
        }
      }
    }
  }
  profile_update("  End Interpolator::rip_real_mat method");
  return true;
}
