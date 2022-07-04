#include "bz.hpp"

using namespace brille;
using namespace brille::lattice;

void BrillouinZone::_moveinto_prim(const LVec<double>& Q, LVec<double>& q, LVec<int>& tau, const LVec<double>& pa, const LVec<double>& pb, const LVec<double>& pc) const {
  // Q, q, tau *must* be in the primitive lattice already

  // the face centre points and normals in the primitive lattice
  auto normals = this->get_primitive_normals();
  normals = normals/norm(normals); // ensure they're normalised
  auto taus = (2.0*this->get_primitive_points()).round(); // the points *must* be the face center vectors!
  auto tau_lens = norm(taus);
  size_t max_count = taus.size(0);
  // ensure that q and tau can hold each q_i and tau_i
  q.resize(Q.size(0));
  tau.resize(Q.size(0));
  auto snQ = utils::u2s<long long, ind_t>(Q.size(0));
#pragma omp parallel for default(none) shared(Q, tau, q, pa, pb, pc, normals, taus, tau_lens, snQ, max_count) schedule(dynamic)
  for (long long si=0; si<snQ; si++){
    auto i = utils::s2u<ind_t, long long>(si);
    auto tau_i = Q.view(i).round();
    auto q_i = Q.view(i) - tau_i;
    auto last_shift = tau_i;
    size_t count{0};
    while (count++ < max_count && !point_inside_all_planes(pa, pb, pc, q_i)){
      auto qi_dot_normals = dot(q_i , normals);
      auto N_hkl = (qi_dot_normals/ tau_lens).round().to_std();
      if (std::any_of(N_hkl.begin(), N_hkl.end(), [](int a){return a > 0;})){
        int max_nm{0};
        ind_t max_at{0};
        for (ind_t j=0; j<N_hkl.size(); ++j) {
          // protect against oscillating by ±τ
          if (N_hkl[j] > 0 && N_hkl[j] >= max_nm &&
              (0 == max_nm ||
               (norm(taus.view(j) + last_shift).all(brille::cmp::gt, 0.) &&
                qi_dot_normals[j] > qi_dot_normals[max_at]))
          ) {
            max_nm = N_hkl[max_at = j];
          }
        }
        q_i -= taus.view(max_at) * static_cast<double>(max_nm); // ensure we subtract LVec<double>
        tau_i += taus.view(max_at) * max_nm; // but add LVec<int>
        last_shift = taus.view(max_at) * max_nm;
      }
    }
    q.set(i, q_i);
    tau.set(i, tau_i);
  }
}

bool BrillouinZone::moveinto(const LVec<double>& Q, LVec<double>& q, LVec<int>& tau, const int threads) const {
  profile_update("BrillouinZone::moveinto called with ",threads," threads");
  omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
  bool already_same = _inner.is_same(Q.lattice());
  LVec<double> Qprim(Q.type(), _inner);
  LVec<double> qprim(q.type(), _inner);
  LVec<int> tauprim(tau.type(), _inner);
  PrimitiveTransform PT(_outer.bravais());
  bool transform_needed = PT.does_anything() && _outer.is_same(Q.lattice());
  if (!(already_same || transform_needed)){
    std::string msg = "Q points provided to BrillouinZone::moveinto must be ";
    msg += "in the standard or primitive lattice used to define ";
    msg += "the BrillouinZone object";
    throw std::runtime_error(msg);
  }
  if (transform_needed)  Qprim = transform_to_primitive(_outer,Q);
  const auto & Qsl = transform_needed ? Qprim : Q;
  auto & qsl = transform_needed ? qprim : q;
  auto & tausl = transform_needed? tauprim : tau;

  auto [a, b, c] = _first.planes();
  auto pa = transform_to_primitive(_outer, a);
  auto pb = transform_to_primitive(_outer, b);
  auto pc = transform_to_primitive(_outer, c);

  this->_moveinto_prim(Qsl, qsl, tausl, pa, pb, pc);
  if (transform_needed){ // then we need to transform back q and tau
    q   = transform_from_primitive(_outer,qsl);
    tau = transform_from_primitive(_outer,tausl);
  }
  auto allinside = this->isinside(q);
  if (std::count(allinside.begin(), allinside.end(), false) > 0){
    std::transform(allinside.begin(), allinside.end(), allinside.begin(), [](const auto & x){return !x;});
    info_update(Q.extract(allinside).size(0), " of ", Q.size(0), " still outside?");
    info_update("outside Q:\nnp.array(\n", Q.extract(allinside).to_string(), ")");
    info_update("outside tau:\nnp.array(\n", tau.extract(allinside).to_string(), ")");
    info_update("outside q\nnp.array(\n", q.extract(allinside).to_string(), ")");
    info_update("outside q(xyz)\nnp.array(\n", q.extract(allinside).xyz().to_string(), ")");
    throw std::runtime_error("Not all points inside Brillouin zone");
    // return false;
  }
  return true; // otherwise, an error has been thrown
}

bool BrillouinZone::ir_moveinto(const LVec<double>& Q, LVec<double>& q, LVec<int>& tau, std::vector<size_t>& Ridx, std::vector<size_t>& invRidx, const int threads) const {
  profile_update("BrillouinZone::ir_moveinto called with ",threads," threads");
  omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
  /* The Point group symmetry information has all rotation matrices defined
     * in the conventional unit cell -- which is our `_outer`.
     * Consequently, we must work in the outer lattice here.  */
  if (!_outer.is_same(Q.lattice()))
    throw std::runtime_error("Q points provided to ir_moveinto must be in the standard lattice used to define the BrillouinZone object");
  // ensure q, tau, and Rm can hold one for each Q.
  ind_t nQ = Q.size(0);
  auto Qshape = Q.shape();
  q.resize(Qshape);
  tau.resize(Qshape);
  Ridx.resize(nQ);
  invRidx.resize(nQ);
  // find q₁ₛₜ in the first Brillouin zone and τ ∈ [reciprocal lattice vectors]
  // such that Q = q₁ₛₜ + τ
  this->moveinto(Q, q, tau, threads);
  auto lat = Q.lattice();
  // OpenMP 2 (VS) doesn't like unsigned loop counters
  size_t n_outside{0};
  auto snQ = utils::u2s<long long, ind_t>(nQ);
#pragma omp parallel default(none) shared(Ridx, invRidx, q, lat, snQ) reduction(+:n_outside)
  {
    // get the PointSymmetry object, containing all operations
    PointSymmetry psym = this->get_pointgroup_symmetry();
    auto eidx = psym.find_identity_index();
    std::array<double,3> q_j{0,0,0}; // temporary result storage
    std::vector<std::array<int, 9>> r_transpose;
    for (const auto& r: psym.getall()) r_transpose.push_back(transpose(r));
#pragma omp for schedule(dynamic)
    for (long long si = 0; si < snQ; ++si) {
      auto i = utils::s2u<ind_t, long long>(si);
      bool inside{_inside_wedge_outer(q.view(i))};
      if (inside){
        // any q already in the irreducible zone need no rotation → identity
        invRidx[i] = Ridx[i] = eidx;
      } else {
        // find the jᵗʰ operation which moves qᵢ into the irreducible zone
        for (ind_t j = 0; j < psym.size(); ++j) if (inside) break; else {
            // The point symmetry matrices relate *real space* vectors!
            // We must use their transposes' to rotate reciprocal space vectors.
            utils::multiply_matrix_vector(q_j.data(), r_transpose[j].data(), q.ptr(i));
            auto lq_j = from_std_like(q, q_j);
            if (_inside_wedge_outer(lq_j)) {
              /* store the result */
              // and (Rⱼᵀ)⁻¹ ∈ G, such that Qᵢ = (Rⱼᵀ)⁻¹⋅qᵢᵣ + τᵢ.
              q.set(i, lq_j);   // keep Rⱼᵀ⋅qᵢ as qᵢᵣ
              invRidx[i] = j; // Rⱼ *is* the inverse of what we want for output
              Ridx[i] = psym.get_inverse_index(j); // find the index of Rⱼ⁻¹
              inside = true;
            }
          }
      }
      if (!inside) ++n_outside;
    }
  }
  if (n_outside) for (ind_t i=0; i<nQ; ++i) if (!_inside_wedge_outer(q.view(i))){
        std::string msg = "Q = " + Q.to_string(i);
        msg += " is outside of the irreducible BrillouinZone ";
        msg += " : tau = " + tau.to_string(i) + " , q = " + q.to_string(i);
        throw std::runtime_error(msg);
        return false;
      }
  return true; // otherwise we hit the runtime error above
}


bool BrillouinZone::ir_moveinto_wedge(const LVec<double>& Q, LVec<double>& q, std::vector<size_t>& R, const int threads) const {
  omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
  /* The Pointgroup symmetry information comes from, effectively, spglib which
  has all rotation matrices defined in the conventional unit cell -- which is
  our `_outer`. Consequently we must work in the _outer here.  */
  if (!_outer.is_same(Q.lattice()))
    throw std::runtime_error("Q points provided to ir_moveinto must be in the standard lattice used to define the BrillouinZone object");
  // ensure q and R can hold one for each Q.
  ind_t nQ = Q.size(0);
  auto Qshape = Q.shape();
  q.resize(Qshape);
  R.resize(nQ);
  auto lat = Q.lattice();
  // OpenMP 2 (VS) doesn't like unsigned loop counters
  size_t n_outside{0};
  auto snQ = brille::utils::u2s<long long, ind_t>(nQ);
#pragma omp parallel default(none) shared(R, q, Q, lat, snQ) reduction(+:n_outside)
  {
    // get the PointSymmetry object, containing all operations
    auto psym = this->_outer.pointgroup_symmetry();
    if (time_reversal) psym = psym.add_space_inversion();
    auto eidx = psym.find_identity_index();
    std::array<double, 3> q_j{0,0,0}; // temporary result storage
    std::vector<std::array<int, 9>> r_transpose;
    for (const auto& r: psym.getall()) r_transpose.push_back(transpose(r));
#pragma omp for schedule(dynamic)
    for (long long si = 0; si < snQ; ++si) {
      auto i = brille::utils::s2u<ind_t, long long>(si);
      // any q already in the irreducible zone need no rotation → identity
      bool inside{_inside_wedge_outer(Q.view(i))};
      if (inside){
        q.set(i, Q.view(i));
        R[i] = eidx;
      } else {
        // for others find the jᵗʰ operation which moves qᵢ into the irreducible zone
        for (ind_t j = 0; j < psym.size(); ++j) if (inside) break; else {
            // The point symmetry matrices relate *real space* vectors! We must use their transposes' to rotate reciprocal space vectors.
            brille::utils::multiply_matrix_vector(q_j.data(), r_transpose[j].data(), Q.ptr(i));
            auto lq_j = from_std_like(Q, q_j);
            if (_inside_wedge_outer(lq_j)) { /* store the result */
              q.set(i, lq_j); // keep Rⱼᵀ⋅Qᵢ as qᵢᵣ
              R[i] = psym.get_inverse_index(j); // and (Rⱼᵀ)⁻¹ ∈ G, such that Q = (Rⱼᵀ)⁻¹⋅qᵢᵣ
              inside = true;
            }
          }
      }
      if (!inside) ++n_outside;
    }
  }
  if (n_outside > 0) for (ind_t i=0; i<nQ; ++i) if (!_inside_wedge_outer(q.view(i))){
        std::string msg = "Q = " + Q.to_string(i);
        msg += " is outside of the irreducible reciprocal space wedge ";
        msg += " , irQ = " + q.to_string(i);
        throw std::runtime_error(msg);
        return false;
      }
  return true; // otherwise we hit the runtime error above
}
