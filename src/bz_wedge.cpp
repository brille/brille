/* This file is part of brille.

Copyright © 2019,2020 Greg Tucker <greg.tucker@stfc.ac.uk>

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
#include "bz.hpp"
using namespace brille;
using namespace brille::lattice;

bArray<bool> keep_if(const LVec<double>& normals, const LVec<double>& points, const double ftol, const int tol){
  // determine whether we should keep points based on the provided normals.
  // We want to only keep those points in the positive half-space for any one
  // normal, but if all normals contribute to reduce the remaining points to
  // lie in a single plane we instead want to keep all points.
  std::vector<size_t> nop(normals.size(0), 0); // number of not-on-plane points
  for (ind_t i=0; i<normals.size(0); ++i)
    nop[i] = dot(normals.view(i), points).is(brille::cmp::gt, 0., ftol, tol).count();
  // If there are no planes with 0 off-plane points, divide the space
  bArray<bool> keep(points.size(0), 1u, true);
  if (std::find(nop.begin(),nop.end(),0u)==nop.end())
    for (ind_t i=0; i<points.size(0); ++i)
      keep[i] = dot(normals, points.view(i)).all(brille::cmp::ge, 0., ftol, tol);
  return keep;
}

bool have_volume(const bArray<double> & p, const double ftol, const int tol){
  double origin[3]{0,0,0};
  double total{0};
  for (ind_t i=0; i<p.size(0) - 2; ++i) total += std::abs(orient3d(origin, p.ptr(i), p.ptr(i+1), p.ptr(i+2)));
  return !approx_float::scalar(total, 0., ftol, ftol, tol);
}

bool BrillouinZone::wedge_brute_force(const bool special_2_folds, const bool special_mirrors, const bool sort_by_length, const bool sort_one_sym){
  debug_update("Starting from first Brillouin zone\n", _first.python_string());
  std::string pmsg = "BrillouinZone::wedge_brute_force(";
  if (special_2_folds) pmsg += "2-folds";
  if (special_2_folds && special_mirrors) pmsg += ",";
  if (special_mirrors) pmsg += "mirrors";
  pmsg += ")";
  profile_update("Start ",pmsg);
  debug_exec(std::string msg;)
  // Grab the pointgroup symmetry operations
  auto full_ps = _outer.pointgroup_symmetry();
  if (time_reversal) full_ps = full_ps.add_space_inversion();
  // Now restrict the symmetry operations to those with order > 1.
  auto ps = full_ps.higher(1);
  // if the full pointgroup contains only 𝟙 (and -𝟙) then ps is empty and  this algorithm can not be used.
  // (triclinic systems are handled outside of this method)
  assert(ps.size() > 0);

  // Get and combine the characteristic points of the first Brillouin zone:
  // The face centres, face corners, and mid-face-edge points.
  auto special = cat(0, this->get_points(), this->get_vertices(), this->get_half_edges());
//  info_update("First Bz\n", _first.python_string(), "with special points\n,np.array(", get_xyz(special).to_string(), ")");

  std::vector<size_t> perm(ps.size());
  std::iota(perm.begin(), perm.end(), 0u); // 0u, 1u, 2u,...
  if (sort_by_length){
    // ps is sorted in increasing rotation-order order, but we want to sort
    // by increasing stationary-axis length (so that, e.g, [100] is dealt with
    // before [111]).
    std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b){
        LVec<int> lq(LengthUnit::angstrom, _outer, 2u); // FIXME stationary axes are real-space vectors. Fixed?
        lq.set(0u, ps.axis(a));
        lq.set(1u, ps.axis(b));
        return lq.dot(0u,0u) < lq.dot(1u,1u);
    });
  } else {
    // ps is sorted in increasing rotation-order order, but we want to sort
    // by decreasing isometry so that we deal with rotoinversions last
    std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b){
      return ps.isometry(a) > ps.isometry(b);
    });
  }
  ps.permute(perm); // ps now sorted

  debug_update("Sorted symmetry operations:");
  for (size_t i=0; i<ps.size(); ++i){
    debug_update(ps.get(i), " ", ps.axis(i), " ", ps.isometry(i));
  }

  bArray<bool> keep;

  // Keep track of *which* normals go into determining the irreducible wedge:
  // count how many normals we might need
  // 2-folds divide space with one normal
  // 3-, 4-, 6-folds need two normals
  size_t n_expected{0};
  for (size_t i=0; i<ps.size(); ++i)
    n_expected += (ps.order(i)==2) ? 1u : 2u;
  LVec<double> cutting_normals(LengthUnit::inverse_angstrom, _outer, n_expected);
  ind_t n_cut{0};

  std::vector<bool> sym_unused(ps.size(), true);
  // deal with any two-fold axes along êᵢ first:
  // The stationary vector of each rotation is a real space vector!
  auto vec = LDVec<int>(_outer, 1u);// must be int since ps.axis returns array<int,3>
  auto nrm = LQVec<double>(_outer, 1u);

  std::vector<std::array<double,3>> eiv{{{1,0,0}},{{0,1,0}},{{0,0,1}},{{1,1,0}},{{1,-1,0}},{{1,0,1}},{{0,1,1}},{{1,0,-1}},{{0,1,-1}},{{1,1,1}}};
  auto eiv_array = bArray<double>::from_std(eiv);
  auto eis = LQVec<double>(_outer, eiv_array);
  auto reis = LDVec<double>(_outer, eiv_array);
  size_t is_nth_ei;

  debug_update_if(special_2_folds,"Deal with 2-fold rotations with axes along highest-symmetry directions first");
  if (special_2_folds) for (size_t i=0; i<ps.size(); ++i) if (ps.isometry(i)==2){
    vec.set(0, ps.axis(i));
    // First check if this stationary axis is along a reciprocal space vector
    is_nth_ei = norm(cross(eis, vec.star())).is(brille::cmp::eq, 0., float_tolerance, approx_tolerance).first();
    if (is_nth_ei < 9 /* This is less than great practice */){
      debug_update("2-fold axis ",i," is ei* No. ",is_nth_ei);
      debug_update(ps.get(i)," along ",ps.axis(i), " is ", vec.star(),"  == ", eis.to_string(is_nth_ei));
      ind_t e1, e2;
      switch (is_nth_ei){
        case 0: /* (100)⋆ */ e1=2; e2=0; /* n = (001)×(100)⋆ */ break;
        case 1: /* (010)⋆ */ e1=0; e2=1; /* n = (100)×(010)⋆ */ break;
        case 2: /* (001)⋆ */ e1=1; e2=2; /* n = (010)×(001)⋆ */ break;
        case 3: /* (110)⋆ */ e1=3; e2=2; /* n = (110)×(001)⋆ */ break;
        case 4: /* (1̄10)⋆ */ e1=2; e2=4; /* n = (001)×(1̄10)⋆ */ break;
        case 5: /* (101)⋆ */ e1=1; e2=5; /* n = (010)×(101)⋆ */ break;
        case 6: /* (011)⋆ */ e1=6; e2=0; /* n = (011)×(100)⋆ */ break;
        case 7: /* (10̄1)⋆ */ e1=7; e2=1; /* n = (10̄1)×(010)⋆ */ break;
        case 8: /* (01̄1)⋆ */ e1=0; e2=8; /* n = (100)×(01̄1)* */ break;
        default: e1=0; e2=0;
      }
      // The plane normal is the cross product of the first real space vector
      // (expressed in units of the reciprocal lattice) and the second
      // reciprocal space vector.
      nrm.set(0, cross(reis.view(e1).star(), eis.view(e2)));
      nrm /= norm(nrm);
      auto n_count = norm(cross(eis, nrm)).is(brille::cmp::eq, 0., float_tolerance, approx_tolerance).count();
      if (n_count == 1){
        // keep any special points beyond the bounding plane
        keep = dot(nrm, special).is(brille::cmp::ge, 0., float_tolerance, approx_tolerance);
        verbose_update("1 Keeping special points with",nrm.to_string(0)," dot p >= 0:\n",cat(1, special ,1.0 * keep).to_string());
        special = special.extract(keep);
        debug_update("Retained special points\nnp.array(", get_xyz(special).to_string(), ")");
        sym_unused[i] = false;
        cutting_normals.set(n_cut++, nrm);
      } else {
        debug_update("can not be used here since there are ", n_count, " matches to its normal, ",nrm.to_string(0)," in eis");
      }
    }
    // Stationary axis along real space basis vector
    is_nth_ei = norm(cross(reis, vec)).is(brille::cmp::eq, 0., float_tolerance, approx_tolerance).first();
    if (sym_unused[i] && is_nth_ei < 2){
      debug_update("2-fold axis ",i," is ei No. ",is_nth_ei);
      switch (is_nth_ei){
        case 0: nrm.set(0, eiv[1]); break; /* (100) → n = (010)* */
        case 1: nrm.set(0, eiv[0]); break; /* (010) → n = (100)* */
        default: throw std::runtime_error("Unreachable path reached!");
      }
      nrm /= norm(nrm);
      // keep any special points beyond the bounding plane
      keep = dot(nrm, special).is(brille::cmp::ge, 0., float_tolerance, approx_tolerance);
      verbose_update("Keeping special (LVec) points p, with (LVec)",nrm.to_string(0)," dot p >= 0:\n",cat(1, special ,1.0 * keep).to_string());
      special = special.extract(keep);
      debug_update("Retained special points\nnp.array(", get_xyz(special).to_string(), ")");
      sym_unused[i] = false;
      cutting_normals.set(n_cut++, nrm);
    }
  }
  debug_update_if(special_mirrors,"Deal with mirror planes");
  if (special_mirrors) for (size_t i=0; i<ps.size(); ++i) if (ps.isometry(i)==-2){
    vec.set(0, ps.axis(i)); // the mirror plane normal is in the direct lattice
    nrm.set(0, vec.star()); // and we want the normal in the reciprocal lattice
    nrm /= norm(nrm);
    keep = dot(nrm, special).is(brille::cmp::ge, 0., float_tolerance, approx_tolerance);
    // we need at least three points (plus Γ) to define a polyhedron
    // If we are not keeping three points, check if applying the mirror plane
    // pointing the other way works for us:
    if (keep.count() < 3 || !have_volume(special.extract(keep), float_tolerance, approx_tolerance)){
      nrm = -1*nrm; // - change nrm since we save it for later
      keep = dot(nrm, special).is(brille::cmp::ge, 0., float_tolerance, approx_tolerance);
    }
    if (keep.count() > 2 && have_volume(special.extract(keep), float_tolerance, approx_tolerance)){
      verbose_update("3 Keeping special points with ",nrm.to_string(0)," dot p >= 0:\n",cat(1, special ,1.0 * keep).to_string());
      special = special.extract(keep);
      debug_update("Retained special points\nnp.array(", get_xyz(special).to_string(), ")");
      sym_unused[i] = false;
      cutting_normals.set(n_cut++, nrm);
    }
  }
  debug_update("Now figure out how all special points are related for each symmetry operation");
  // Find which points are mapped onto equivalent points by each symmetry operation
  std::vector<std::vector<ind_t>> one_sym;
  std::vector<ind_t> one_type, type_order;
  std::vector<bool> unfound(special.size(0), true), type_unfound;
  for (ind_t i=0; i<ps.size(); ++i) if (sym_unused[i]){
    debug_update("Unused symmetry ",i," with order ", ps.order(i));
    debug_update(ps.get(i));
    one_sym.clear();
    for (auto b: unfound) b = true;
    for (ind_t j=0; j<special.size(0); ++j) if (unfound[j]){
      one_type.clear();
      one_type.push_back(j);
      unfound[j] = false;
      for (ind_t k=j+1; k<special.size(0); ++k) {
        if (unfound[k] && special.match(k, j, transpose(ps.get(i)), -ps.order(i), float_tolerance,
                                        approx_tolerance)) { // -order checks all possible rotations
          one_type.push_back(k);
          unfound[k] = false;
        }
      }
      debug_update("Points equivalent to ",j," for symmetry ",i,":",one_type);
      // sort the equivalent points by their relative order for this operation
      // such that Rⁱj ≡ type_order[i]
      type_order.clear();
      type_order.insert(type_order.begin(), ps.order(i), special.size(0)); // set default to (non-indexable) size of the special array
      type_unfound.clear(); type_unfound.insert(type_unfound.begin(), one_type.size(), true);
      for (int o=0; o<ps.order(i); ++o)
      for (ind_t k=0; k<one_type.size(); ++k) if (type_unfound[k]){
        if (special.match(one_type[k], j, transpose(ps.get(i)), o, float_tolerance, approx_tolerance)){
        // if (special.match(one_type[k], j, transpose(ps.get_proper(i)), o)){
          type_order[o] = one_type[k];
          type_unfound[k] = false;
        }
      }
      // and store the sorted equivalent indices
      debug_update("Which are sorted by their rotation order:",type_order);
      one_sym.push_back(type_order);
    }
    // sort one_sym by the number of valid (indexable) equivalent points
    if (sort_one_sym){
      std::sort(one_sym.begin(), one_sym.end(),
        [&](std::vector<ind_t>& a, std::vector<ind_t>& b){
          size_t as = a.size() - std::count(a.begin(), a.end(), special.size(0));
          size_t bs = b.size() - std::count(b.begin(), b.end(), special.size(0));
          return as > bs;
        }
      );
    }
    debug_update("one_sym planes:\nnp.array(",get_xyz(special).to_string(),"),",one_sym);
    size_t keep_count;
//    LVec<double> pt0(_outer, 1u), pt1(_outer, 1u);
//    for (size_t s=0; s<one_sym.size(); ++s) if (sym_unused[i]/*always true?*/){
//      // we need at least two equivalent points, ideally there will be the same number as the order
//      if (one_sym[s].size() > static_cast<size_t>(std::count(one_sym[s].begin(), one_sym[s].end(), special.size(0))+1)){
//        type_order = one_sym[s];
//        debug_update("Highest-multiplicity distinct point type:",type_order);
//      // auto loc = std::find_if(one_sym.begin(), one_sym.end(), [&](std::vector<size_t>x){return x.size()==static_cast<size_t>(ps.order(i));});
//      // if (loc != one_sym.end()){
//      //   size_t idx = loc - one_sym.begin();
//        // type_order = one_sym[idx];
//        // grab the rotation stationary axis
//        vec.set(0, ps.axis(i));
//        for (size_t j=0, k=1; j<type_order.size(); ++j, k=(j+1)%type_order.size())
//        if (sym_unused[i] && type_order[j]<special.size(0) && type_order[k]<special.size(0)){
//          if (ps.order(i)>2){
//            // hold the two special points in their own LVec<double> (!not <int>!)
//            pt0 = special.view(type_order[j]);
//            pt1 = special.view(type_order[k]);
//            // we have two plane normals to worry about:
//            debug_update("Stationary vector",vec.to_string(0)," and special points",pt0.to_string(0)," and",pt1.to_string(0));
//            // find both cross products, remembering that we want normals
//            // pointing *into* the wedge.
//            nrm.resize(2);
//            if ( dot(pt1, cross(vec.star(), pt0)).all(brille::cmp::lt,0., approx_tolerance) ){
//              // the rotation is left handed, so swap the special points
//              nrm.set(0, cross(vec.star(), pt1));
//              nrm.set(1, cross(pt0, vec.star()));
//            } else {
//              nrm.set(0, cross(vec.star(), pt0));
//              nrm.set(1, cross(pt1, vec.star()));
//            }
//            debug_update("give normals:", nrm.to_string(0), " and", nrm.to_string(1));
//            // now check that all special points are inside of the wedge defined by the normals
//          } else {
//            // order == 2, so only one normal to worry about:
//            nrm = cross(vec.star(), special.view(type_order[j]));
//            // make sure we don't remove all points out of the plane containing
//            // the rotation axis and the two special points
//            if (dot(nrm, special).is(brille::cmp::gt,0., approx_tolerance).count() == 0)
//              nrm *= -1; // switch the cross product
//          }
//          // check to make sure that using these normals do not remove all but a plane of points
//          keep = keep_if(nrm, special, approx_tolerance);
//          keep_count = keep.count();
//          // We need at least three points (plus Γ) to define a polyhedron.
//          // Also skip the extraction if we are keeping all points
//          if (keep_count > 2 && keep_count < keep.size(0)){
//            debug_update("Keeping special points (keep, h, k, l):\n",cat(1, 1.0 * keep, special).to_string());
//            special = special.extract(keep);
//            sym_unused[i]=false;
//            for (ind_t nc=0; nc<nrm.size(0); ++nc)
//                cutting_normals.set(n_cut++, nrm.view(nc));
//          }
//        }
//      }
//    }
    for (auto & s : one_sym) {
      if (sym_unused[i]/*always true?*/) {
        // we need at least two equivalent points, ideally there will be the same number as the order
        if (s.size() > static_cast<size_t>(std::count(s.begin(), s.end(), special.size(0)) + 1)) {
          type_order = s;
          debug_update("Highest-multiplicity distinct point type:", s);
          // grab the rotation stationary axis
          vec.set(0, ps.axis(i));
          auto qvec = vec.star();
          auto not_parallel = [&](const LVec<double> &x) {
            return dot(x, qvec).abs().is(brille::cmp::neq, norm(x) * norm(qvec), float_tolerance, float_tolerance,
                                         approx_tolerance).all();
          };
          for (size_t j = 0, k = 1; j < s.size(); ++j, k = (j + 1) % s.size()) {
            if (sym_unused[i] && s[j] < special.size(0) && s[k] < special.size(0)) {
              bool non_parallel = true;
              if (ps.order(i) > 2) {
                // hold the two special points in their own LVec<double> (!not <int>!)
                const auto pt0{special.view(s[j])};
                const auto pt1{special.view(s[k])};
                non_parallel = not_parallel(pt0) && not_parallel(pt1);
                // we have two plane normals to worry about:
                debug_update("Stationary reciprocal vector", qvec.to_string(0), " and special points", pt0.to_string(0),
                             " and", pt1.to_string(0));
                // find both cross products, remembering that we want normals
                // pointing *into* the wedge.
                nrm.resize(2);
                if (dot(pt1, cross(qvec, pt0)).all(brille::cmp::lt, 0., float_tolerance, approx_tolerance)) {
                  // the rotation is left handed, so swap the special points
                  nrm.set(0, cross(qvec, pt1));
                  nrm.set(1, cross(pt0, qvec));
                } else {
                  nrm.set(0, cross(qvec, pt0));
                  nrm.set(1, cross(pt1, qvec));
                }
                debug_update("give normals:", nrm.to_string(0), " and", nrm.to_string(1));
                // now check that all special points are inside of the wedge defined by the normals
              } else {
                non_parallel = not_parallel(special.view(s[j]));
                // order == 2, so only one normal to worry about:
                nrm = cross(qvec, special.view(s[j]));
                // make sure we don't remove all points out of the plane containing
                // the rotation axis and the two special points
                if (dot(nrm, special).is(brille::cmp::gt, 0., float_tolerance, approx_tolerance).count() == 0)
                  nrm *= -1; // switch the cross product
              }
              if (non_parallel) {
                nrm /= norm(nrm);
                // check to make sure that using these normals do not remove all but a plane of points
                keep = keep_if(nrm, special, float_tolerance, approx_tolerance);
                keep_count = keep.count();
                // We need at least three points (plus Γ) to define a polyhedron.
                // Also skip the extraction if we are keeping all points
                if (keep_count > 2 && keep_count < keep.size(0)) {
                  verbose_update("Keeping special points (keep, h, k, l):\n", cat(1, 1.0 * keep, special).to_string());
                  special = special.extract(keep);
                  debug_update("Retained special points\nnp.array(", get_xyz(special).to_string(), ")");
                  sym_unused[i] = false;
                  for (ind_t nc = 0; nc < nrm.size(0); ++nc)
                    cutting_normals.set(n_cut++, nrm.view(nc));
                }
              }
            }
          }
        }
      }
    }
  }

  // debug_update("Remaining special points\n", special.to_string());
  if (n_cut > 0) { /*protect against view(0,0)*/
    auto cn = cutting_normals.view(0, n_cut); // the cutting direction is opposite the normal
    auto [ca, cb, cc] = plane_points_from_normal(-1.0 * cn, 0.0 * cn);
    debug_update("Now cut the first Brillouin zone by normal(s)\n", get_xyz(-1.0 * cn).to_string());
    debug_update("With on-plane points,\n",get_xyz(ca).to_string(),get_xyz(cb).to_string(),get_xyz(cc).to_string());
    _irreducible = _first.cut(ca, cb, cc, float_tolerance, approx_tolerance);
  } else {
    info_update("No cutting normals in wedge search");
  }
  // copy functionality of set_ir_vertices, which set the normals as well
  if (this->check_ir_polyhedron()){
    verbose_update("Irreducible Polyhedron check succeeded. Set irreducible wedge normals");
    this->set_ir_wedge_normals(this->get_ir_polyhedron_wedge_normals());
  }

//  // append the Γ point onto the list of special points:
//  special.resize(special.size()+1u);
//  for (size_t i=0; i<3u; ++i) special.insert(0, special.size()-1, i);
//  // and then form the irreducible polyhedron by finding the convex hull
//  // of these points: (which sets the wedge normals too)
//  this->set_ir_vertices(special);

  /*If the full pointsymmetry has space inversion *and* there is only a single
    stationary rotation/rotoinversion axis the thus-far-found polyhedron has
    twice the volume of the true irreducible polyhedron.
    We need to split the polyhedron in half by a plane perpendicular to *one of*
    the full-Brillouin zone normal directions.
    This is easier for spacegroups with orthogonal basis vectors -- we can
    usually pick the stationary axis -- but is trickier if the basis vectors
    are not orthogonal.
  */
  if (full_ps.has_space_inversion()){
    // ps only contains operations *with* a stationary axis (described in real space)
    auto all_axes = LDVec<int>(_outer, bArray<int>::from_std(ps.axes()));
    auto ir_volume = _irreducible.volume();
    auto goal_volume = _first.volume()/static_cast<double>(full_ps.size());
    auto aaiu = all_axes.is_unique();
    auto cnt = all_axes.unique().size(0);
    if ((ir_volume > goal_volume && cnt == 1) || brille::approx_float::scalar(ir_volume, 2.0*goal_volume) ){
      debug_update("Deal with -1 since there is only one stationary axis (or doubled volume for some other reason)");
      auto bz_n = this->get_normals();
      // auto wg_n = this->get_ir_wedge_normals(); // !! This is empty if the thus-far-found volume is wrong
      auto wg_n = this->get_ir_polyhedron_wedge_normals(); // This pulls directly from the ir_polyhedron
      auto gamma = bz_n.view(0) * 0.;
      for (ind_t i=0; i<bz_n.size(0); ++i){
        auto [a, b, c] = plane_points_from_normal(bz_n.view(i), gamma);
        auto div = _irreducible.one_cut(a, b, c, float_tolerance, approx_tolerance);
        if (approx_float::scalar(div.volume(), goal_volume, float_tolerance, float_tolerance, approx_tolerance)){
          // set div to be the ir_polyhedron
          _irreducible = div;
          // add the new normal to wedge normals list
          // extract instead of view to avoid copying the whole bz_n Array
          // when the inversion happens.
          wg_n.append(/*dimension to expand*/ 0u, -bz_n.extract(i));
          this->set_ir_wedge_normals(wg_n);
          break;
        }
      }
      if (approx_float::scalar(_irreducible.volume(), ir_volume, float_tolerance, float_tolerance, approx_tolerance)){
        debug_update("Polyhedron volume still double expected.");
        bool proceed=true;
        // check for other dividing planes :/
        //LVec<double> cij(this->_outer, 1u);
        for (ind_t i=0; proceed && i<bz_n.size(0)-1; ++i)
        for (ind_t j=i+1; proceed && j<bz_n.size(0); ++j){
          auto [a, b, c] = plane_points_from_normal(bz_n.cross(i, j), gamma);
          auto div = _irreducible.one_cut(a, b, c, float_tolerance, approx_tolerance);
          if (approx_float::scalar(div.volume(), goal_volume, float_tolerance, float_tolerance, approx_tolerance)){
            _irreducible = div;
            wg_n.append(0u, -bz_n.cross(i,j));
            this->set_ir_wedge_normals(wg_n);
            proceed = false;
          }
        }
      }
    }
  }
  bool success = this->check_ir_polyhedron();
  if (!success) pmsg += " failed";
  profile_update("  End ",pmsg);
  return success;
}


bool BrillouinZone::wedge_triclinic(){
  /* Assuming that this is spacegroup P -1 we have inversion symmetry only.
     We always want to find a set of planes that divides symmetry equivalent
     regions of space. With only -1 symmetry we have to make a choice as to
     *which* space inversion we care about.
     Any single 'wedge' normal will divide the first Brillouin zone -- why not
     choose one perpendicular to a first Brillouin zone face? (or, failing that,
     one of [100], [010], [001], or [111]).
  */
  using namespace brille;
  auto normals = this->get_normals(); // the first Brillouin zone normals
  auto origin = 0 * normals.view(0);
  for (ind_t i=0; i<normals.size(0); ++i){
    this->set_ir_wedge_normals(normals.view(i));
    auto [a, b, c] = plane_points_from_normal(-1*normals.view(i), origin);
    _irreducible = _first.one_cut(a, b, c);
    if (this->check_ir_polyhedron()) return true;
  }
  std::vector<std::array<double,3>> sv{{{1,0,0}},{{0,1,0}},{{0,0,1}},{{1,1,1}}};
  auto nrm = LQVec<double>(_outer, bArray<double>::from_std(sv));
  for (ind_t i=0; i<nrm.size(0); ++i){
    this->set_ir_wedge_normals(nrm.view(i));
    auto [a, b, c] = plane_points_from_normal(-1*nrm.view(i), origin);
    _irreducible = _first.one_cut(a, b, c);
    if (this->check_ir_polyhedron()) return true;
  }
  return false;
}
