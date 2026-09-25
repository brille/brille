#include "bz.hpp"
#include "plane_set.hpp"
#include <atomic>

using namespace brille;
using namespace brille::lattice;

/*!\brief Run the primitive lattice move into routine
 *
 * Initially wrwitten to operate on a subset of the provided Q array
 * from `start` to `stop` (exclusive) such that it could be called in a parallel
 * loop. However, something makes serial execution faster than even 12-core
 * parallel execution -- likely something with memory access in the underpinning
 * `Array2` class?
 * */

namespace {
/* The first Brillouin zone as plain numbers, for moving many points into it.

The previous implementation built lattice-aware temporaries (views,
differences, cross and dot products) for every plane and point: profiling
showed ~46 % of its time in shared_ptr reference counting and ~11 % in
malloc/free, which also stopped it from scaling with threads. This does the
same arithmetic, in the same order and with the same helpers, on plain arrays,
so the results are unchanged.
*/
struct FirstZone {
  using v3 = std::array<double, 3>;
  PlaneSet planes;                // the zone's faces, for point_inside_all_planes
  std::vector<v3> normals;        // unit face normals (hkl)
  std::vector<std::array<int, 3>> taus; // the face-centre reciprocal lattice vectors
  std::vector<double> tau_lens;
  std::array<double, 9> metric{}; // of the points' LengthUnit, for dot(q, normals)

  FirstZone(const LVec<double>& pa, const LVec<double>& pb, const LVec<double>& pc,
            const LVec<double>& n, const LVec<int>& t, const Array2<double>& tl, double f, int at)
  : planes(pa, pb, pc, f, at), metric(pa.lattice().metric(pa.type())) {
    for (ind_t j=0; j<n.size(0); ++j){
      normals.push_back({n.val(j,0), n.val(j,1), n.val(j,2)});
      taus.push_back({t.val(j,0), t.val(j,1), t.val(j,2)});
      tau_lens.push_back(tl.val(j,0));
    }
  }
  [[nodiscard]] bool inside(const v3& q) const { return planes.inside(q); }
  // move one point into the first Brillouin zone
  void move(const v3& Q, double* q_out, int* tau_out) const {
    std::array<int, 3> tau, last_shift;
    v3 q;
    for (int k=0; k<3; ++k){
      tau[k] = static_cast<int>(std::round(Q[k]));
      q[k] = Q[k] - tau[k];
    }
    last_shift = tau;
    const auto max_count = static_cast<ind_t>(taus.size());
    ind_t count{0};
    std::vector<double> qdn(normals.size());
    std::vector<int> N(normals.size());
    while (count++ < max_count && !inside(q)){
      v3 mq;
      brille::utils::mul_mat_vec(mq.data(), 3u, metric.data(), q.data()); // dot(q, normals)
      for (size_t j=0; j<normals.size(); ++j){
        double d{0};
        for (int k=0; k<3; ++k) d += mq[k] * normals[j][k];
        qdn[j] = d;
        N[j] = static_cast<int>(std::round(qdn[j] / tau_lens[j]));
      }
      if (std::any_of(N.begin(), N.end(), [](int n){return n > 0;})) {
        int max_nm{0};
        size_t max_at{0};
        for (size_t j=0; j<N.size(); ++j){
          const bool shift_nonzero = taus[j][0] + last_shift[0] != 0 || taus[j][1] + last_shift[1] != 0
                                  || taus[j][2] + last_shift[2] != 0;
          // protect against oscillating by ±τ
          if (N[j] > 0 && N[j] >= max_nm && (0 == max_nm || (shift_nonzero && qdn[j] > qdn[max_at]))) {
            max_nm = N[max_at = j];
          }
        }
        for (int k=0; k<3; ++k){
          q[k] -= static_cast<double>(taus[max_at][k]) * static_cast<double>(max_nm);
          tau[k] += taus[max_at][k] * max_nm;
          last_shift[k] = taus[max_at][k] * max_nm;
        }
      }
    }
    for (int k=0; k<3; ++k){ q_out[k] = q[k]; tau_out[k] = tau[k]; }
  }
};

/* The irreducible wedge test, _inside_wedge_outer, as plain numbers.

dot(normals, p).all(c, 0, ftol, atol) on lattice vectors, repeated for every
point and every trial operation, was ~99 % of ir_moveinto after the first-zone
step was made fast. This keeps its arithmetic and comparisons, including that
the le_ge case compares without tolerance (Array2::all's two-argument form).
*/
struct IrWedge {
  using v3 = std::array<double, 3>;
  std::vector<v3> weighted;       // metric · normal, as same_lattice_dot computes it
  std::vector<v3> plain;          // the normals, for the star-lattice case
  bool star{false};
  brille::cmp expr;
  double ftol;
  int atol;

  IrWedge(const LVec<double>& normals, const LVec<double>& like, brille::cmp c, double f, int at)
  : expr(c), ftol(f), atol(at) {
    if (normals.size(0) == 0) return;
    star = normals.star_lattice(like);
    const auto metric = normals.lattice().metric(normals.type());
    for (ind_t i=0; i<normals.size(0); ++i){
      v3 n{normals.val(i,0), normals.val(i,1), normals.val(i,2)}, mn;
      brille::utils::mul_mat_vec(mn.data(), 3u, metric.data(), n.data());
      weighted.push_back(mn);
      plain.push_back(n);
    }
  }
  [[nodiscard]] bool inside(const v3& p) const {
    if (weighted.empty()) return true;
    std::array<double, 64> small{};
    std::vector<double> large;
    double* d = small.data();
    if (weighted.size() > small.size()) { large.resize(weighted.size()); d = large.data(); }
    for (size_t i=0; i<weighted.size(); ++i){
      double out{0};
      if (star) {
        for (int k=0; k<3; ++k) out += plain[i][k] * p[k];
        out *= brille::math::two_pi;
      } else {
        for (int k=0; k<3; ++k) out += weighted[i][k] * p[k];
      }
      d[i] = out;
    }
    auto all = [&](brille::cmp e, double t, int n){
      brille::Comparer<double,double> op(e, t, t, n);
      for (size_t i=0; i<weighted.size(); ++i) if (!op(d[i], 0.)) return false;
      return true;
    };
    if (brille::cmp::le_ge == expr) return all(brille::cmp::le, 0., 1) || all(brille::cmp::ge, 0., 1);
    return all(expr, ftol, atol);
  }
};

std::pair<LVec<double>, LVec<int>>
fast_moveinto_prim(const LVec<double>& Q, const LVec<double>& normals, const LVec<int>& taus,
                   const Array2<double>& tau_lens, double ftol, int atol,
                   const LVec<double>& pa, const LVec<double>& pb, const LVec<double>& pc, int threads) {
  const FirstZone zone(pa, pb, pc, normals, taus, tau_lens, ftol, atol);
  LVec<double> q(Q.type(), Q.lattice(), Q.size(0));
  LVec<int> tau(Q.type(), Q.lattice(), Q.size(0));
  const auto pool = ThreadPool::getInstance();
  if (threads > 0) pool->resize(threads); else pool->resize();
  const auto workers = pool->size();
  const ind_t n = Q.size(0);
  for (size_t worker=0; worker<workers; ++worker){
    auto [first, last] = thread_slice(n, workers, worker);
    pool->enqueue([&, first=first, last=last](){
      for (size_t i=first; i<last; ++i){
        const auto ii = static_cast<ind_t>(i);
        zone.move({Q.val(ii,0), Q.val(ii,1), Q.val(ii,2)}, q.ptr(ii), tau.ptr(ii));
      }
    });
  }
  pool->wait();
  return std::make_pair(q, tau);
}
} // namespace


bool BrillouinZone::moveinto(const LVec<double>& Q, LVec<double>& q, LVec<int>& tau, const int threads) const {
  profile_update("BrillouinZone::moveinto called with ",threads," threads");
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
  if (transform_needed)  Qprim = parallel_transform_to_primitive(_outer, Q, threads);
  const auto & Qsl = transform_needed ? Qprim : Q;
//  auto & qsl = transform_needed ? qprim : q;
//  auto & tausl = transform_needed? tauprim : tau;

  auto [a, b, c] = _first.planes();
  auto pa = parallel_transform_to_primitive(_outer, a, threads);
  auto pb = parallel_transform_to_primitive(_outer, b, threads);
  auto pc = parallel_transform_to_primitive(_outer, c, threads);

//  // more than one thread == slow??
//  this->_moveinto_prim(Qsl, qsl, tausl, pa, pb, pc, threads);
//  if (transform_needed){ // then we need to transform back q and tau
//    q   = parallel_transform_from_primitive(_outer, qsl, threads);
//    tau = parallel_transform_from_primitive(_outer, tausl, threads);
//  }

  // single threaded is faster?!
  // the face centre points and normals in the primitive lattice
  auto normals = this->get_primitive_normals();
  normals = normals/norm(normals); // ensure they're normalised
  auto taus = (2.0*this->get_primitive_points()).round(); // the points *must* be the face center vectors!
  auto tau_lens = norm(taus);
  auto q_tau = fast_moveinto_prim(Qsl, normals, taus, tau_lens, float_tolerance, approx_tolerance, pa, pb, pc, threads);
  if (transform_needed){ // then we need to transform back q and tau
    q   = parallel_transform_from_primitive(_outer, q_tau.first, threads);
    tau = parallel_transform_from_primitive(_outer, q_tau.second, threads);
  } else {
    q = q_tau.first;
    tau = q_tau.second;
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
  profile_update("BrillouinZone::moveinto finished with ",threads," threads");
  return true; // otherwise, an error has been thrown
}

bool BrillouinZone::ir_moveinto(const LVec<double>& Q, LVec<double>& q, LVec<int>& tau, std::vector<size_t>& Ridx, std::vector<size_t>& invRidx, const int threads) const {
  profile_update("BrillouinZone::ir_moveinto called with ",threads," threads");
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
//  auto snQ = utils::u2s<long long, ind_t>(nQ);
  /*FIXME The following implementation is *SLOWER* than single threaded!*/
//#pragma omp parallel default(none) shared(Ridx, invRidx, q, lat, snQ) reduction(+:n_outside)
//  {
//    // get the PointSymmetry object, containing all operations
//    PointSymmetry psym = this->get_pointgroup_symmetry();
//    auto eidx = psym.find_identity_index();
//    std::array<double,3> q_j{0,0,0}; // temporary result storage
//    std::vector<std::array<int, 9>> r_transpose;
//    for (const auto& r: psym.getall()) r_transpose.push_back(transpose(r));
//#pragma omp for
//    for (long long si = 0; si < snQ; ++si) {
//      auto i = utils::s2u<ind_t, long long>(si);
//      bool inside{_inside_wedge_outer(q.view(i))};
//      if (inside){
//        // any q already in the irreducible zone need no rotation → identity
//        invRidx[i] = Ridx[i] = eidx;
//      } else {
//        // find the jᵗʰ operation which moves qᵢ into the irreducible zone
//        for (ind_t j = 0; j < psym.size(); ++j) if (inside) break; else {
//            // The point symmetry matrices relate *real space* vectors!
//            // We must use their transposes' to rotate reciprocal space vectors.
//            utils::multiply_matrix_vector(q_j.data(), r_transpose[j].data(), q.ptr(i));
//            auto lq_j = from_std_like(q, q_j);
//            if (_inside_wedge_outer(lq_j)) {
//              /* store the result */
//              // and (Rⱼᵀ)⁻¹ ∈ G, such that Qᵢ = (Rⱼᵀ)⁻¹⋅qᵢᵣ + τᵢ.
//              q.set(i, lq_j);   // keep Rⱼᵀ⋅qᵢ as qᵢᵣ
//              invRidx[i] = j; // Rⱼ *is* the inverse of what we want for output
//              Ridx[i] = psym.get_inverse_index(j); // find the index of Rⱼ⁻¹
//              inside = true;
//            }
//          }
//      }
//      if (!inside) ++n_outside;
//    }
//  }

//  /*FIXME As is the following implementation*/
//  PointSymmetry psym = this->get_pointgroup_symmetry();
//  auto eidx = psym.find_identity_index();
//  std::vector<std::array<int, 9>> r_transpose;
//  for (const auto& r: psym.getall()) r_transpose.push_back(transpose(r));
//#pragma omp parallel for default(none) shared(Ridx, invRidx, q, lat, snQ, eidx, r_transpose, psym) reduction(+:n_outside)
//  for (long long si = 0; si < snQ; ++si) {
//    auto i = utils::s2u<ind_t, long long>(si);
//    bool inside{_inside_wedge_outer(q.view(i))};
//    if (inside){
//      // any q already in the irreducible zone need no rotation → identity
//      invRidx[i] = Ridx[i] = eidx;
//    } else {
//      std::array<double,3> q_j{0,0,0}; // temporary result storage
//      // find the jᵗʰ operation which moves qᵢ into the irreducible zone
//      for (ind_t j = 0; j < psym.size(); ++j) if (inside) break; else {
//          // The point symmetry matrices relate *real space* vectors!
//          // We must use their transposes' to rotate reciprocal space vectors.
//          utils::multiply_matrix_vector(q_j.data(), r_transpose[j].data(), q.ptr(i));
//          auto lq_j = from_std_like(q, q_j);
//          if (_inside_wedge_outer(lq_j)) {
//            /* store the result */
//            // and (Rⱼᵀ)⁻¹ ∈ G, such that Qᵢ = (Rⱼᵀ)⁻¹⋅qᵢᵣ + τᵢ.
//            q.set(i, lq_j);   // keep Rⱼᵀ⋅qᵢ as qᵢᵣ
//            invRidx[i] = j; // Rⱼ *is* the inverse of what we want for output
//            Ridx[i] = psym.get_inverse_index(j); // find the index of Rⱼ⁻¹
//            inside = true;
//          }
//        }
//    }
//    if (!inside) ++n_outside;
//  }

  PointSymmetry psym = this->get_pointgroup_symmetry();
  const auto eidx = psym.find_identity_index();
  std::vector<std::array<int, 9>> r_transpose;
  std::vector<size_t> inverse_index;
  for (size_t j=0; j<psym.size(); ++j){
    r_transpose.push_back(transpose(psym.get(j)));
    inverse_index.push_back(psym.get_inverse_index(j));
  }
  const IrWedge wedge(get_ir_wedge_normals(), q, no_ir_mirroring ? brille::cmp::ge : brille::cmp::le_ge,
                      float_tolerance, approx_tolerance);
  std::atomic<size_t> outside{0};
  const auto pool = ThreadPool::getInstance();
  if (threads > 0) pool->resize(threads); else pool->resize();
  const auto workers = pool->size();
  for (size_t worker=0; worker<workers; ++worker){
    auto [first, last] = thread_slice(nQ, workers, worker);
    pool->enqueue([&, first=first, last=last](){
      for (size_t si=first; si<last; ++si){
        const auto i = static_cast<ind_t>(si);
        const std::array<double,3> qi{q.val(i,0), q.val(i,1), q.val(i,2)};
        bool inside = wedge.inside(qi);
        if (inside){
          invRidx[i] = Ridx[i] = eidx;
        } else {
          std::array<double,3> q_j{0,0,0};
          for (size_t j = 0; j < r_transpose.size() && !inside; ++j) {
            utils::multiply_matrix_vector(q_j.data(), r_transpose[j].data(), qi.data());
            if (wedge.inside(q_j)) {
              for (int k=0; k<3; ++k) q.ptr(i)[k] = q_j[k]; // keep Rⱼᵀ⋅qᵢ as qᵢᵣ
              invRidx[i] = j; // Rⱼ *is* the inverse of what we want for output
              Ridx[i] = inverse_index[j]; // find the index of Rⱼ⁻¹
              inside = true;
            }
          }
        }
        if (!inside) ++outside;
      }
    });
  }
  pool->wait();
  n_outside = outside;

  profile_update("BrillouinZone::ir_moveinto finished with ",threads," threads");
  if (n_outside) for (ind_t i=0; i<nQ; ++i) if (!_inside_wedge_outer(q.view(i))){
        std::string msg = "Q = " + Q.to_string(i);
        msg += " is outside of the irreducible BrillouinZone ";
        msg += " : tau = " + tau.to_string(i) + " , q = " + q.to_string(i);
        throw std::runtime_error(msg);
        return false;
      }
  return true; // otherwise we hit the runtime error above
}


// bool BrillouinZone::ir_moveinto_wedge(const LVec<double>& Q, LVec<double>& q, std::vector<size_t>& R, const int threads) const {
//   omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
//   /* The Pointgroup symmetry information comes from, effectively, spglib which
//   has all rotation matrices defined in the conventional unit cell -- which is
//   our `_outer`. Consequently we must work in the _outer here.  */
//   if (!_outer.is_same(Q.lattice()))
//     throw std::runtime_error("Q points provided to ir_moveinto must be in the standard lattice used to define the BrillouinZone object");
//   // ensure q and R can hold one for each Q.
//   ind_t nQ = Q.size(0);
//   auto Qshape = Q.shape();
//   q.resize(Qshape);
//   R.resize(nQ);
//   auto lat = Q.lattice();
//   // OpenMP 2 (VS) doesn't like unsigned loop counters
//   size_t n_outside{0};
//   auto snQ = brille::utils::u2s<long long, ind_t>(nQ);
// #pragma omp parallel default(none) shared(R, q, Q, lat, snQ) reduction(+:n_outside)
//   {
//     // get the PointSymmetry object, containing all operations
//     auto psym = this->_outer.pointgroup_symmetry();
//     if (time_reversal) psym = psym.add_space_inversion();
//     auto eidx = psym.find_identity_index();
//     std::array<double, 3> q_j{0,0,0}; // temporary result storage
//     std::vector<std::array<int, 9>> r_transpose;
//     for (const auto& r: psym.getall()) r_transpose.push_back(transpose(r));
// #pragma omp for schedule(dynamic)
//     for (long long si = 0; si < snQ; ++si) {
//       auto i = brille::utils::s2u<ind_t, long long>(si);
//       // any q already in the irreducible zone need no rotation → identity
//       bool inside{_inside_wedge_outer(Q.view(i))};
//       if (inside){
//         q.set(i, Q.view(i));
//         R[i] = eidx;
//       } else {
//         // for others find the jᵗʰ operation which moves qᵢ into the irreducible zone
//         for (ind_t j = 0; j < psym.size(); ++j) if (inside) break; else {
//             // The point symmetry matrices relate *real space* vectors! We must use their transposes' to rotate reciprocal space vectors.
//             brille::utils::multiply_matrix_vector(q_j.data(), r_transpose[j].data(), Q.ptr(i));
//             auto lq_j = from_std_like(Q, q_j);
//             if (_inside_wedge_outer(lq_j)) { /* store the result */
//               q.set(i, lq_j); // keep Rⱼᵀ⋅Qᵢ as qᵢᵣ
//               R[i] = psym.get_inverse_index(j); // and (Rⱼᵀ)⁻¹ ∈ G, such that Q = (Rⱼᵀ)⁻¹⋅qᵢᵣ
//               inside = true;
//             }
//           }
//       }
//       if (!inside) ++n_outside;
//     }
//   }
//   if (n_outside > 0) for (ind_t i=0; i<nQ; ++i) if (!_inside_wedge_outer(q.view(i))){
//         std::string msg = "Q = " + Q.to_string(i);
//         msg += " is outside of the irreducible reciprocal space wedge ";
//         msg += " , irQ = " + q.to_string(i);
//         throw std::runtime_error(msg);
//         return false;
//       }
//   return n_outside == 0;
// }


bool BrillouinZone::ir_moveinto_wedge(const LVec<double>& Q, LVec<double>& q, std::vector<size_t>& R, const int threads) const {
  /* The Pointgroup symmetry information comes from, effectively, spglib which
  has all rotation matrices defined in the conventional unit cell -- which is
  our `_outer`. Consequently we must work in the _outer here.  */
  if (!_outer.is_same(Q.lattice()))
    throw std::runtime_error("Q points provided to ir_moveinto must be in the standard lattice used to define the BrillouinZone object");
  // ensure q and R can hold one for each Q.
  const ind_t nQ = Q.size(0);
  const auto Qshape = Q.shape();
  q.resize(Qshape);
  R.resize(nQ);
  auto lat = Q.lattice();

  const auto pool = ThreadPool::getInstance();
  if (threads > 0) pool->resize(threads); else pool->resize();
  const auto workers = pool->size();
  std::vector outside_counts(workers, 0u);
  auto make_task = [&](const size_t thread) {
    // setup for this thread:
    auto psym = this->_outer.pointgroup_symmetry();
    if (time_reversal) psym = psym.add_space_inversion();
    const auto eidx = psym.find_identity_index();
    std::array<double, 3> q_j{0,0,0}; // temporary result storage
    std::vector<std::array<int, 9>> r_transpose;
    for (const auto& r: psym.getall()) r_transpose.push_back(transpose(r));
    const auto [first, last] = thread_slice(nQ, workers, thread);
    // the actual task that the thread should execute // capture the thread number to avoid all threads sharing it
    // psym, eidx and r_transpose are local to make_task, so the task must own copies
    auto task = [&,frst=first,lst=last,iam=thread,psym=std::move(psym),eidx,r_transpose=std::move(r_transpose)]() {
      for (size_t i=frst; i<lst; ++i) {
        bool inside{_inside_wedge_outer(Q.view(i))};
        if (inside){
          q.set(i, Q.view(i));
          R[i] = eidx;
        } else {
          // for others find the jᵗʰ operation which moves qᵢ into the irreducible zone
          for (ind_t j = 0; j < psym.size(); ++j) {
            if (inside) break;
            // The point symmetry matrices relate *real space* vectors! We must use their transposes' to rotate reciprocal space vectors.
            utils::multiply_matrix_vector(q_j.data(), r_transpose[j].data(), Q.ptr(i));
            if (const auto lq_j = from_std_like(Q, q_j); _inside_wedge_outer(lq_j)) {
              /* store the result */
              q.set(i, lq_j); // keep Rⱼᵀ⋅Qᵢ as qᵢᵣ
              R[i] = psym.get_inverse_index(j); // and (Rⱼᵀ)⁻¹ ∈ G, such that Q = (Rⱼᵀ)⁻¹⋅qᵢᵣ
              inside = true;
            }
          }
        }
        if (!inside) ++outside_counts[iam];
      }
    };
    return task;
  };
  for (size_t thread=0; thread < workers; ++thread) {
    pool->enqueue(make_task(thread));
  }
  pool->wait();
  const auto n_outside = std::accumulate(outside_counts.begin(), outside_counts.end(), 0u);

  if (n_outside > 0) for (ind_t i=0; i<nQ; ++i) if (!_inside_wedge_outer(q.view(i))){
        std::string msg = "Q = " + Q.to_string(i);
        msg += " is outside of the irreducible reciprocal space wedge ";
        msg += " , irQ = " + q.to_string(i);
        throw std::runtime_error(msg);
        return false;
      }
  return n_outside == 0;
}