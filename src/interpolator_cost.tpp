template<class T>
template<typename S>
void
Interpolator<T>::add_cost_mix(const ind_t i0, const ind_t i1, std::vector<S>& cost, const bool arbitrary_phase_allowed) const {
  auto x0 = data_.slice(i0); // → (1,), (B,), (Y,), or (B,Y)
  auto x1 = data_.slice(i1);
  // can this be extended to account for the arbitrary phase too?
  S s_cost{0}, v_cost{0}, m_cost{0};
  ind_t span = this->branch_span();
  ind_t b_ = this->branches();
  auto e_ = _elements;
  if (arbitrary_phase_allowed){ // if the _vectorfun uses the Hermitian angle, e^iθ *never* matters.
    auto phased = std::unique_ptr<T[]>(new T[span]);
    for (ind_t i=0; i<b_; ++i){
      auto b0i = (b_>1? x0.slice(i) : x0).contiguous_copy(); // (1,) or (Y,)
      const T* p0i = b0i.ptr(); // don't combine with prior line in case a copy was required
      for (ind_t j=0; j<b_; ++j){
        auto b1j = (b_>1 ? x1.slice(j) : x1).contiguous_copy();
        brille::utils::inplace_antiphase(span, p0i, b1j.ptr(), phased.get());
        if (e_[0]) s_cost = this->_scalarfun(e_[0], p0i, phased.get());
        if (e_[1]) v_cost = this->_vectorfun(e_[1], p0i+e_[0], phased.get()+e_[0]);
        if (e_[2]){
          m_cost = 0;
          ind_t o{e_[0]+e_[1]};
          for (ind_t m=0; m<e_[2]/3; ++m)
            m_cost += brille::utils::frobenius_distance(3u, p0i+o+9u*m, phased.get()+o+9u*m);
        }
        cost[i*b_+j] += _costmult[0]*s_cost + _costmult[1]*v_cost + _costmult[2]*m_cost;
      }
    }
  } else {
    for (ind_t i=0; i<b_; ++i){
      auto b0i = (b_>1 ? x0.slice(i) : x0).contiguous_copy(); // (1,) or (Y,)
      const T* p0i = b0i.ptr(); // don't combine with prior line in case a copy was required
      for (ind_t j=0; j<b_; ++j){
        auto b1j = (b_>1 ? x1.slice(j) : x1).contiguous_copy();
        const T* p1j = b1j.ptr();
        if (e_[0]) s_cost = this->_scalarfun(e_[0], p0i, p1j);
        if (e_[1]) v_cost = this->_vectorfun(e_[1], p0i+e_[0], p1j+e_[0]);
        if (e_[2]){
          m_cost = 0;
          ind_t o{e_[0]+e_[1]};
          for (ind_t m=0; m<e_[2]/3; ++m)
            m_cost += brille::utils::frobenius_distance(3u, p0i+o+9u*m, p1j+o+9u*m);
        }
        cost[i*b_+j] += _costmult[0]*s_cost + _costmult[1]*v_cost + _costmult[2]*m_cost;
      }
    }
  }
}

template<class T>
template<typename S>
void
Interpolator<T>::add_cost_vec(const ind_t i0, const ind_t i1, std::vector<S>& cost, const bool arbitrary_phase_allowed) const {
  auto x0 = data_.slice(i0); // → (B, V, 3) Array
  auto x1 = data_.slice(i1); // → (B, V, 3) Array
  ind_t b_ = x0.size(1);
  ind_t v_ = x0.size(2)*3u;
  if (arbitrary_phase_allowed){ // if the _vectorfun uses the Hermitian angle, e^iθ *never* matters.
    auto phased = std::unique_ptr<T[]>(new T[v_]);
    for (ind_t i=0; i<b_; ++i){
      auto b0i = x0.slice(i).contiguous_copy(); // a copy is only made if the slice *is not already* contiguous
      for (ind_t j=0; j<b_; ++j){
        auto b1j = x1.slice(j).contiguous_copy();
        brille::utils::inplace_antiphase(v_, b0i.ptr(), b1j.ptr(), phased.get());
        cost[i*b_+j] += _costmult[1]*this->_vectorfun(v_, b0i.ptr(), phased.get());
      }
    }
  } else {
    for (ind_t i=0; i<b_; ++i){
      auto b0i = x0.slice(i).contiguous_copy();
      for (ind_t j=0; j<b_; ++j){
        auto b1j = x1.slice(j).contiguous_copy();
        cost[i*b_+j] += _costmult[1]*this->_vectorfun(v_, b0i.ptr(), b1j.ptr());
      }
    }
  }
}

template<class T>
template<typename S>
void
Interpolator<T>::add_cost_mat(const ind_t i0, const ind_t i1, std::vector<S>& cost, const bool arbitrary_phase_allowed) const {
  auto x0 = data_.slice(i0); // → (B, M, 3, 3) Array
  auto x1 = data_.slice(i1); // → (B, M, 3, 3) Array
  ind_t b_ = x0.size(1);
  ind_t m_ = x0.size(2);
  if (arbitrary_phase_allowed){
    ind_t allm = m_ * 9u;
    T* phased = new T[allm]();
    for (ind_t i=0; i<b_; ++i){
      auto b0i = x0.slice(i).contiguous_copy(); // row-ordered contiguous (M, 3, 3)
      for (ind_t j=0; j<b_; ++j){
        auto b1j = x1.slice(j).contiguous_copy();
        brille::utils::inplace_antiphase(allm, b0i.ptr(), b1j.ptr(), phased);
        // the Frobenius distance for a N by M matrix is the same as the vector distance
        // of the N*M array treated as a single vector.
        // This means we can skip-over
        S fd{0};
        for (ind_t m = 0; m < m_; ++m)
          fd += brille::utils::frobenius_distance(3u, b0i.ptr()+9u*m, phased+9u*m);
        cost[i*b_+j] += _costmult[2]*fd;
      }
    }
    delete[] phased;
  } else {
    for (ind_t i=0; i<b_; ++i){
      auto b0i = x0.slice(i).contiguous_copy(); // row-ordered contiguous (M, 3, 3)
      for (ind_t j=0; j<b_; ++j){
        auto b1j = x1.slice(j).contiguous_copy();
        S fd{0};
        for (ind_t m = 0; m < m_; ++m)
          fd += brille::utils::frobenius_distance(3u, b0i.ptr()+9u*m, b1j.ptr()+9u*m);
        cost[i*b_+j] += _costmult[2]*fd;
      }
    }
  }
}
