#include "bz.hpp"

using namespace brille;

void BrillouinZone::_moveinto_prim(const LQVec<double>& Q, LQVec<double>& q, LQVec<int>& tau, const LQVec<double>& pa, const LQVec<double>& pb, const LQVec<double>& pc) const {
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
// #pragma omp parallel for default(none) shared(Q, tau, q, pa, pb, pc, points, normals, taus, tau_lens, snQ, max_count) schedule(dynamic)
  for (long long si=0; si<snQ; si++){
    auto i = utils::s2u<ind_t, long long>(si);
    auto tau_i = Q.view(i).round();
    auto q_i = Q.view(i) - tau_i;
    auto last_shift = tau_i;
    size_t count{0};
    while (count++ < max_count && !point_inside_all_planes(pa, pb, pc, q_i)){
      auto qi_dot_normals = dot(q_i , normals);
      auto N_hkl = (qi_dot_normals/ tau_lens).round().to_std();
      auto q_i_d_n = qi_dot_normals.to_std();
      if (std::any_of(N_hkl.begin(), N_hkl.end(), [](int a){return a > 0;})){
        int max_nm{0};
        ind_t max_at{0};
        for (ind_t j=0; j<N_hkl.size(); ++j) {
          if (qi_dot_normals[j] != q_i_d_n[j]) {
            info_update("Why is ", qi_dot_normals[j], " different from ",
                        q_i_d_n[j]);
          }
          // protect against oscillating by ±τ
          if (N_hkl[j] > 0 && N_hkl[j] >= max_nm &&
              (0 == max_nm ||
               (norm(taus.view(j) + last_shift).all(brille::cmp::gt, 0.) &&
                qi_dot_normals[j] > qi_dot_normals[max_at]))
              //              (0 == max_nm ||
              //              (norm(taus.view(j)+last_shift).all(brille::cmp::gt,
              //              0.) && q_i_d_n[j] > q_i_d_n[max_at]))
          ) {
            max_nm = N_hkl[max_at = j];
          }
        }
        q_i -= taus.view(max_at) * static_cast<double>(max_nm); // ensure we subtract LQVec<double>
        tau_i += taus.view(max_at) * max_nm; // but add LQVec<int>
        last_shift = taus.view(max_at) * max_nm;
      }
    }
    q.set(i, q_i);
    tau.set(i, tau_i);
  }
}

bool BrillouinZone::moveinto(const LQVec<double>& Q, LQVec<double>& q, LQVec<int>& tau, const int threads) const {
  profile_update("BrillouinZone::moveinto called with ",threads," threads");
  omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
  bool already_same = lattice.issame(Q.get_lattice());
  LQVec<double> Qprim(lattice);
  LQVec<double> qprim(lattice);
  LQVec<int> tauprim(lattice);
  PrimitiveTransform PT(outerlattice.get_bravais_type());
  bool transform_needed = PT.does_anything() && outerlattice.issame(Q.get_lattice());
  if (!(already_same || transform_needed)){
    std::string msg = "Q points provided to BrillouinZone::moveinto must be ";
    msg += "in the standard or primitive lattice used to define ";
    msg += "the BrillouinZone object";
    throw std::runtime_error(msg);
  }
  if (transform_needed)  Qprim = transform_to_primitive(outerlattice,Q);
  const LQVec<double> & Qsl = transform_needed ? Qprim : Q;
  LQVec<double> & qsl = transform_needed ? qprim : q;
  LQVec<int> & tausl = transform_needed? tauprim : tau;

  LQVec<double> a, b, c, pa, pb, pc;
  std::tie(a,b,c) = first_poly.planes();
  pa = transform_to_primitive(outerlattice, a);
  pb = transform_to_primitive(outerlattice, b);
  pc = transform_to_primitive(outerlattice, c);

//
//  // the face centre points and normals in the primitive lattice:
//  auto points = this->get_primitive_points();
//  auto normals = this->get_primitive_normals();
//  normals = normals/norm(normals); // ensure they're normalised
//  auto taus = (2.0*points).round();
//  auto taulen = norm(taus);
//  size_t max_count = taus.size(0);
//  // ensure that qsl and tausl can hold each qi and taui
//  qsl.resize(Qsl.size(0));
//  tausl.resize(Qsl.size(0));
//  auto snQ = utils::u2s<long long, ind_t>(Qsl.size(0));
//#pragma omp parallel for default(none) shared(Qsl, tausl, qsl, points, normals, taus, taulen, snQ, max_count) schedule(dynamic)
//  for (long long si=0; si<snQ; si++){
//    auto i = utils::s2u<ind_t, long long>(si);
//    auto taui = Qsl.view(i).round();
//    auto qi = Qsl.view(i) - taui;
//    auto last_shift = taui;
//    size_t count{0};
//    // FIXME insert tolerance here
//    while (count++ < max_count && dot(normals, qi-points).any(brille::cmp::gt,0.)){
//      auto qi_dot_normals = dot(qi , normals);
//      auto Nhkl = (qi_dot_normals/taulen).round().to_std();
//      auto qidn = qi_dot_normals.to_std();
//      if (std::any_of(Nhkl.begin(), Nhkl.end(), [](int a){return a > 0;})){
//        int maxnm{0};
//        ind_t maxat{0};
//        for (ind_t j=0; j<Nhkl.size(); ++j)
//          // protect against oscillating by ±τ
//          if (Nhkl[j]>0 && Nhkl[j]>=maxnm && (0==maxnm || (norm(taus.view(j)+last_shift).all(brille::cmp::gt, 0.) && qidn[j]>qidn[maxat]))){
//            maxnm = Nhkl[maxat=j];
//          }
//        qi -= taus.view(maxat) * static_cast<double>(maxnm); // ensure we subtract LQVec<double>
//        taui += taus.view(maxat) * maxnm; // but add LQVec<int>
//        last_shift = taus.view(maxat) * maxnm;
//      }
//    }
//    qsl.set(i, qi);
//    tausl.set(i, taui);
//  }
  this->_moveinto_prim(Qsl, qsl, tausl, pa, pb, pc);
  if (transform_needed){ // then we need to transform back q and tau
    q   = transform_from_primitive(outerlattice,qsl);
    tau = transform_from_primitive(outerlattice,tausl);
  }
  auto allinside = this->isinside(q);
  if (std::count(allinside.begin(), allinside.end(), false) > 0){

//    info_update("inside Q:\nnp.array(", Q.extract(allinside).to_string(), ")");
//    info_update("inside tau:\nnp.array(", tau.extract(allinside).to_string(), ")");
//    info_update("inside q\nnp.array(", q.extract(allinside).to_string(), ")");
//    info_update("inside Qsl:\nnp.array(", Qsl.extract(allinside).to_string(), ")");
//    info_update("inside tausl:\nnp.array(", tausl.extract(allinside).to_string(), ")");
//    info_update("inside qsl\nnp.array(", qsl.extract(allinside).to_string(), ")");
//    for (ind_t i=0; i<Q.size(0); ++i) if (!allinside[i]){
//      info_update("Q  =",Q.to_string(i)  ," tau  =",tau.to_string(i)  ," q  =",q.to_string(i));
//      info_update("Qsl=",Qsl.to_string(i)," tausl=",tausl.to_string(i)," qsl=",qsl.to_string(i),"\n");
//    }
    std::transform(allinside.begin(), allinside.end(), allinside.begin(), [](const auto & x){return !x;});
    info_update(Q.extract(allinside).size(0), " of ", Q.size(0), " still outside?");
    info_update("outside Q:\nnp.array(\n", Q.extract(allinside).to_string(), ")");
    info_update("outside tau:\nnp.array(\n", tau.extract(allinside).to_string(), ")");
    info_update("outside q\nnp.array(\n", q.extract(allinside).to_string(), ")");
    info_update("outside Qsl:\nnp.array(\n", Qsl.extract(allinside).to_string(), ")");
    info_update("outside tausl:\nnp.array(\n", tausl.extract(allinside).to_string(), ")");
    info_update("outside qsl\nnp.array(\n", qsl.extract(allinside).to_string(), ")");
    throw std::runtime_error("Not all points inside Brillouin zone");
    // return false;
  }
  return true; // otherwise, an error has been thrown
}

bool BrillouinZone::ir_moveinto(const LQVec<double>& Q, LQVec<double>& q, LQVec<int>& tau, std::vector<size_t>& Ridx, std::vector<size_t>& invRidx, const int threads) const {
  profile_update("BrillouinZone::ir_moveinto called with ",threads," threads");
  omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
  /* The Point group symmetry information has all rotation matrices defined
     * in the conventional unit cell -- which is our `outerlattice`.
     * Consequently, we must work in the outer lattice here.  */
  if (!outerlattice.issame(Q.get_lattice()))
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
  auto lat = Q.get_lattice();
  // OpenMP 2 (VS) doesn't like unsigned loop counters
  size_t n_outside{0};
  auto snQ = utils::u2s<long long, ind_t>(nQ);
#pragma omp parallel default(none) shared(Ridx, invRidx, q, lat, snQ) reduction(+:n_outside)
  {
    // get the PointSymmetry object, containing all operations
    PointSymmetry psym = this->get_pointgroup_symmetry();
    auto eidx = psym.find_identity_index();
    LQVec<double> qj(lat, 1u); // a place to hold the multiplication result
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
            utils::multiply_matrix_vector(qj.ptr(0), r_transpose[j].data(), q.ptr(i));
            if (_inside_wedge_outer(qj)) {
              /* store the result */
              // and (Rⱼᵀ)⁻¹ ∈ G, such that Qᵢ = (Rⱼᵀ)⁻¹⋅qᵢᵣ + τᵢ.
              q.set(i, qj);   // keep Rⱼᵀ⋅qᵢ as qᵢᵣ
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


bool BrillouinZone::ir_moveinto_wedge(const LQVec<double>& Q, LQVec<double>& q, std::vector<size_t>& R, const int threads) const {
  omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
  /* The Pointgroup symmetry information comes from, effectively, spglib which
  has all rotation matrices defined in the conventional unit cell -- which is
  our `outerlattice`. Consequently we must work in the outerlattice here.  */
  if (!outerlattice.issame(Q.get_lattice()))
    throw std::runtime_error("Q points provided to ir_moveinto must be in the standard lattice used to define the BrillouinZone object");
  // ensure q and R can hold one for each Q.
  ind_t nQ = Q.size(0);
  auto Qshape = Q.shape();
  q.resize(Qshape);
  R.resize(nQ);
  auto lat = Q.get_lattice();
  // OpenMP 2 (VS) doesn't like unsigned loop counters
  size_t n_outside{0};
  auto snQ = brille::utils::u2s<long long, ind_t>(nQ);
#pragma omp parallel default(none) shared(R, q, Q, lat, snQ) reduction(+:n_outside)
  {
    // get the PointSymmetry object, containing all operations
    PointSymmetry psym = this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
    auto eidx = psym.find_identity_index();
    LQVec<double> qj(lat, 1u);
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
            brille::utils::multiply_matrix_vector(qj.ptr(0), r_transpose[j].data(), Q.ptr(i));
            if (_inside_wedge_outer(qj)) { /* store the result */
              q.set(i, qj); // keep Rⱼᵀ⋅Qᵢ as qᵢᵣ
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
