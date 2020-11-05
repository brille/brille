// template<class T>
// void Interpolator<T>::interpolate_at_mix(
//   const std::vector<std::vector<ind_t>>& permutations,
//   const std::vector<ind_t>& indices,
//   const std::vector<double>& weights,
//   data_t<T>& out,
//   const ind_t to,
//   const bool arbitrary_phase_allowed
// ) const {
//   if (indices.size()==0 || weights.size()==0)
//     throw std::logic_error("Interpolation requires input data!");
//   ind_t span = this->branch_span();
//   ind_t b_ = this->branches();
//   verbose_update("Combining\n",data_.extract(indices).to_string(),"with weights ", weights);
//   auto ox = out.slice(to);
//   if (arbitrary_phase_allowed){
//     auto d0 = data_.slice(indices[0]); // a reference point for the phase calculation
//     for (size_t x=0; x<indices.size(); ++x){
//       auto dx = data_.slice(x); // (1,), (B,), (Y,), or (B,Y)
//       for (ind_t b=0; b < b_; ++b){
//         auto d0b = b_>1 ? d0.slice(b) : d0; // (1,) or (Y,)
//         auto oxb = b_>1 ? ox.slice(b) : ox; // (1,) or (Y,)
//         auto dxp = b_>1 ? dx.slice(permutations[x][b]): dx; // (1,) or (Y,)
//         T *pxb = oxb.ptr();
//         const T *pxp = dxp.ptr();
//         T eith = brille::utils::antiphase(span, d0b.ptr(), pxp);
//         for (ind_t s=0; s<span; ++s) pxb[s] += eith*weights[x]*pxp[s];
//       }
//     }
//   } else {
//     for (size_t x=0; x<indices.size(); ++x){
//       auto dx = data_.slice(x);
//       for (ind_t b=0; b < b_; ++b){
//         auto dxp = b_>1 ? dx.slice(permutations[x][b]) : dx;
//         auto oxb = b_>1 ? ox.slice(b) : ox;
//         T *pxb = oxb.ptr();
//         const T *pxp = dxp.ptr();
//         for (ind_t s=0; s<span; ++s) pxb[s] += weights[x]*dxp[s];
//       }
//     }
//   }
// }


template<class T>
void Interpolator<T>::interpolate_at_mix(
  const std::vector<std::vector<ind_t>>& permutations,
  const std::vector<std::pair<ind_t,double>>& indices_weights,
  data_t<T>& out,
  const ind_t to,
  const bool arbitrary_phase_allowed
) const {
  if (indices_weights.size()==0)
    throw std::logic_error("Interpolation requires input data!");
  ind_t span = this->branch_span();
  std::vector<int> dummy;
  ind_t b_ = this->branches();
  auto ox = out.slice(to);
  if (arbitrary_phase_allowed){
    auto d0 = data_.slice(indices_weights[0].first); // reference point for phase calculation
    std::transform(permutations.begin(), permutations.end(), indices_weights.begin(), std::back_inserter(dummy),
    [&](const std::vector<ind_t>& perm, const std::pair<ind_t,double>& iw){
      auto dx = data_.slice(iw.first);
      for (ind_t b=0; b<b_; ++b)
      {
        auto oxb = (b_>1 ? ox.slice(b) : ox);
        auto d0b = (b_>1 ? d0.slice(b) : d0);
        auto dxp = (b_>1 ? dx.slice(perm[b]) : dx);
        T eith = brille::utils::antiphase(span, d0b.ptr(), dxp.ptr());
        for (ind_t s=0; s<span; ++s) oxb.val(s) += iw.second*eith*dxp.val(s);
      }
      return 1;
    });
  } else {
    std::transform(permutations.begin(), permutations.end(), indices_weights.begin(), std::back_inserter(dummy),
    [&](const std::vector<ind_t>& perm, const std::pair<ind_t,double>& iw){
      auto dx = data_.slice(iw.first);
      for (ind_t b=0; b<b_; ++b)
      {
        auto oxb = (b_>1 ? ox.slice(b) : ox);
        auto dxp = (b_>1 ? dx.slice(perm[b]) : dx);
        T *pxb = oxb.ptr();
        const T *pxp = dxp.ptr();
        for (ind_t s=0; s<span; ++s) pxb[s] += iw.second*pxp[s];
      }
      return 1;
    });
  }
}


// template<class T>
// void Interpolator<T>::interpolate_at_vec(
//   const std::vector<std::vector<ind_t>>& permutations,
//   const std::vector<ind_t>& indices,
//   const std::vector<double>& weights,
//   data_t<T>& out,
//   const ind_t to,
//   const bool arbitrary_phase_allowed
// ) const {
//   if (indices.size()==0 || weights.size()==0)
//     throw std::logic_error("Interpolation requires input data!");
//   // data_ *has* shape (n, B, V, 3) with at least the last dimension contiguous
//   const ind_t b_{data_.size(1)}, v_{data_.size(2)};
//   verbose_update("Combining\n",data_.extract(indices).to_string(),"with weights ", weights);
//   auto ox = out.slice(to);
//   if (arbitrary_phase_allowed){
//     auto d0 = data_.slice(indices[0]); // (B,V,3) a reference point for the phase calculation
//     for (size_t x=0; x<indices.size(); ++x){
//       auto dx = data_.slice(x); // (B,V,3)
//       for (ind_t b=0; b<b_; ++b){
//         auto oxb = ox.slice(b); // (v,3)
//         auto dxp = dx.slice(permutations[x][b]); // (v,3)
//         T eith = brille::utils::unsafe_antiphase(d0.slice(b), dxp);
//         for (ind_t v=0; v<v_; ++v){
//           T *oxbv = oxb.ptr(v), *dxpv = dxp.ptr(v);
//           for (ind_t i=0; i<3u; ++i) oxbv[i] += eith*weights[x]*dxpv[i];
//         }
//       }
//     }
//   } else {
//     for (size_t x=0; x<indices.size(); ++x){
//       auto dx = data_.slice(x); // (B,V,3)
//       for (ind_t b=0; b<b_; ++b){
//         auto oxb = ox.slice(b); // (v,3)
//         auto dxp = dx.slice(permutations[x][b]); // (v,3)
//         for (ind_t v=0; v<v_; ++v){
//           T *oxbv = oxb.ptr(v), *dxpv = dxp.ptr(v);
//           for (ind_t i=0; i<3u; ++i) oxbv[i] += weights[x]*dxpv[i];
//         }
//       }
//     }
//   }
// }


template<class T>
void Interpolator<T>::interpolate_at_vec(
  const std::vector<std::vector<ind_t>>& permutations,
  const std::vector<std::pair<ind_t,double>>& indices_weights,
  data_t<T>& out,
  const ind_t to,
  const bool arbitrary_phase_allowed
) const {
  if (indices_weights.size()==0)
    throw std::logic_error("Interpolation requires input data!");
  std::vector<int> dummy;
  const ind_t b_{data_.size(1)}, v_{data_.size(2)};
  auto ox = out.slice(to);
  if (arbitrary_phase_allowed){
    auto d0 = data_.slice(indices_weights[0].first); // reference point for phase calculation
    std::transform(permutations.begin(), permutations.end(), indices_weights.begin(), std::back_inserter(dummy),
    [&](const std::vector<ind_t>& perm, const std::pair<ind_t,double>& iw){
      auto dx = data_.slice(iw.first);
      for (ind_t b=0; b<b_; ++b)
      {
        auto oxb = ox.slice(b);
        auto dxp = dx.slice(perm[b]);
        T eith = brille::utils::unsafe_antiphase(d0.slice(b), dxp);
        for (ind_t v=0; v<v_; ++v){
          T *oxbv = oxb.ptr(v), *dxpv = dxp.ptr(v);
          for (ind_t i=0; i<3u; ++i) oxbv[i] += eith*iw.second*dxpv[i];
        }
      }
      return 1;
    });
  } else {
    std::transform(permutations.begin(), permutations.end(), indices_weights.begin(), std::back_inserter(dummy),
    [&](const std::vector<ind_t>& perm, const std::pair<ind_t,double>& iw){
      auto dx = data_.slice(iw.first);
      for (ind_t b=0; b<b_; ++b) {
        auto oxb = ox.slice(b);
        auto dxp = dx.slice(perm[b]);
        for (ind_t v=0; v<v_; ++v){
          T *oxbv = oxb.ptr(v), *dxpv = dxp.ptr(v);
          for (ind_t i=0; i<3u; ++i) oxbv[i] +=iw.second*dxpv[i];
        }
      }
      return 1;
    });
  }
}

// template<class T>
// void Interpolator<T>::interpolate_at_mat(
//   const std::vector<std::vector<ind_t>>& permutations,
//   const std::vector<ind_t>& indices,
//   const std::vector<double>& weights,
//   data_t<T>& out,
//   const ind_t to,
//   const bool arbitrary_phase_allowed
// ) const {
//   if (indices.size()==0 || weights.size()==0)
//     throw std::logic_error("Interpolation requires input data!");
//   // data_ *has* shape (n, B, M, 3, 3) with at least the last dimension contiguous
//   const ind_t b_{data_.size(1)}, m_{data_.size(2)};
//   verbose_update("Combining\n",data_.extract(indices).to_string(),"with weights ", weights);
//   auto ox = out.slice(to);
//   if (arbitrary_phase_allowed){
//     auto d0 = data_.slice(indices[0]); // (B,M,3.3) a reference point for the phase calculation
//     for (size_t x=0; x<indices.size(); ++x){
//       auto dx = data_.slice(x); // (B,M,3,3)
//       for (ind_t b=0; b<b_; ++b){
//         auto oxb = ox.slice(b); // (M,3,3)
//         auto dxp = dx.slice(permutations[x][b]); // (M,3,3)
//         T eith = brille::utils::unsafe_antiphase(d0.slice(b), dxp);
//         for (ind_t m=0; m<m_; ++m){
//           T *oxbm = oxb.ptr(m), *dxpm = dxp.ptr(m);
//           for (ind_t i=0; i<9u; ++i) oxbm[i] += eith*weights[x]*dxpm[i];
//         }
//       }
//     }
//   } else {
//     for (size_t x=0; x<indices.size(); ++x){
//       auto dx = data_.slice(x); // (B,M,3,3)
//       for (ind_t b=0; b<b_; ++b){
//         auto oxb = ox.slice(b); // (M,3,3)
//         auto dxp = dx.slice(permutations[x][b]); // (M,3,3)
//         for (ind_t m=0; m<m_; ++m){
//           T *oxbm = oxb.ptr(m), *dxpm = dxp.ptr(m);
//           for (ind_t i=0; i<9u; ++i) oxbm[i] += weights[x]*dxpm[i];
//         }
//       }
//     }
//   }
// }


template<class T>
void Interpolator<T>::interpolate_at_mat(
  const std::vector<std::vector<ind_t>>& permutations,
  const std::vector<std::pair<ind_t,double>>& indices_weights,
  data_t<T>& out,
  const ind_t to,
  const bool arbitrary_phase_allowed
) const {
  if (indices_weights.size()==0)
    throw std::logic_error("Interpolation requires input data!");
  std::vector<int> dummy;
  const ind_t b_{data_.size(1)}, m_{data_.size(2)};
  auto ox = out.slice(to);
  if (arbitrary_phase_allowed){
    auto d0 = data_.slice(indices_weights[0].first);
    std::transform(permutations.begin(), permutations.end(), indices_weights.begin(), std::back_inserter(dummy),
    [&](const std::vector<ind_t>& perm, const std::pair<ind_t,double>& iw){
      auto dx = data_.slice(iw.first);
      for (ind_t b=0; b<b_; ++b) {
        auto oxb = ox.slice(b); // (M,3,3)
        auto dxp = dx.slice(perm[b]); // (M,3,3)
        T eith = brille::utils::unsafe_antiphase(d0.slice(b), dxp);
        for (ind_t m=0; m<m_; ++m){
          T *oxbm = oxb.ptr(m), *dxpm = dxp.ptr(m);
          for (ind_t i=0; i<9u; ++i) oxbm[i] += eith*iw.second*dxpm[i];
        }
      }
      return 1;
    });
  } else {
    std::transform(permutations.begin(), permutations.end(), indices_weights.begin(), std::back_inserter(dummy),
    [&](const std::vector<ind_t>& perm, const std::pair<ind_t,double>& iw){
      auto dx = data_.slice(iw.first);
      for (ind_t b=0; b<b_; ++b) {
        auto oxb = ox.slice(b);
        auto dxp = dx.slice(perm[b]);
        for (ind_t m=0; m<m_; ++m){
          T *oxbm = oxb.ptr(m), *dxpm = dxp.ptr(m);
          for (ind_t i=0; i<9u; ++i) oxbm[i] += iw.second*dxpm[i];
        }
      }
      return 1;
    });
  }
}
