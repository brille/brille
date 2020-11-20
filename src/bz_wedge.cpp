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

void BrillouinZone::wedge_search(const bool pbv, const bool pok){
  debug_exec(std::string update_msg;)
  // Get the full pointgroup symmetry information
  PointSymmetry fullps = this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
  // And use it to find only the highest-order rotation operation along each
  // unique stationary axis.
  // PointSymmetry rotps = fullps.nfolds(1); // 1 to request only orders>1
  PointSymmetry rotps = fullps.higher(1); // 1 to request only orders>1
  // Get the vectors pointing to each full Brillouin zone facet cetre
  auto xyz = cat(0,this->get_points(), this->get_vertices(), this->get_half_edges());

  // debug_update("xyz=\n", xyz.to_string());

  // if rotps is empty, there are no 2+-fold rotations, act like we have ̄1:
  if (rotps.size()==0){
    debug_update("No 2+-fold rotation operations");
    this->ir_wedge_is_ok(xyz.view(0)); // for ̄1, assigns this->ir_polyhedron
    return;
  }

  int max_order=rotps.order(rotps.size()-1); // since the rotations are sorted

  // rotps is sorted in increasing rotation-order order, but we want to sort
  // by increasing stationary-axis length (so that, e.g, [100] is dealt with
  // before [111]). We could ammend the PointSymmetry.sort method but since the
  // rotations are defined in a lattice that would only work for cubic systems;
  // e.g., in a hexagonal lattice [100], [010], and [1-10] are all the same
  // length. So we need to create a sorting permutation here for rotps to use.
  std::vector<size_t> perm(rotps.size());
  std::iota(perm.begin(), perm.end(), 0u); // 0u, 1u, 2u,...
  std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b){
    LQVec<int> lq(this->outerlattice, 2u);
    lq.set(0u, rotps.axis(a));
    lq.set(1u, rotps.axis(b));
    return lq.dot(0u,0u) < lq.dot(1u,1u);
  });
  // with the permutation found, permute rotps:
  rotps.permute(perm);

  debug_update("Rotations:\n", rotps.getall());

  // make lattice vectors from the stationary axes
  LQVec<int> z(this->outerlattice, bArray<int>::from_std(rotps.axes()));
  // Find ̂zᵢ⋅ẑⱼ
  std::vector<std::vector<int>> dotij(z.size(0));
  double dottmp;
  for (size_t i=0; i<z.size(0); ++i) for (size_t j=0; j<z.size(0); ++j){
    dottmp = std::abs(z.dot(i,j)/z.norm(i)/z.norm(j));
    dotij[i].push_back( brille::approx::scalar(dottmp,0.) ? -1: brille::approx::scalar(dottmp, 1.) ? 1 : 0);
  }
  debug_update("dot(i,j) flag\n", dotij);

  // Find a suitable in-rotation-plane vector for each stationary axis
  // LQVec<double> x(this->outerlattice, bArray<double>(rotps.perpendicular_axes()));
  LQVec<double> x(this->outerlattice, bArray<int>::from_std(rotps.perpendicular_axes()));
  // or pick one of the BZ facet-points if the basis vector preffered flag
  if (pbv) for (size_t i=0; i<x.size(0); ++i) for (size_t j=0; j<xyz.size(0); ++j)
  if (dot(xyz.view(j), z.view(i)).all(brille::cmp::eq, 0.)){
    x.set(i, xyz.view(j));
    break;
  }

  /* Two rotations with ẑᵢ⋅ẑⱼ≠0 should have (Rⁿx̂ᵢ)⋅(Rᵐx̂ⱼ)=1 for some n,m.
     To start, assume that n and m are 0 and find the x̂ᵢ such that x̂ᵢ⋅(ẑᵢ×ẑⱼ)=1
  */
  std::vector<bool> handled(z.size(0), false);
  bool flag;
  for (size_t i=0; i<z.size(0)-1; ++i) for (size_t j=i+1; j<z.size(0); ++j)
  // if zᵢ and zⱼ are neither parallel or perpendicular
  if (0 == dotij[i][j]){
    if (!handled[i]){
      if (!pbv /*basis vectors not preferred*/) x.set(i, z.cross(i, j));
      if (handled[j] && !brille::approx::scalar(x.dot(i,j), 0.) && x.dot(i,j)<0) x.set(i, -x.extract(i));
      handled[i] = true;
    }
    if (!handled[j]){
      if (!pbv /*basis vectors not preferred*/){
        flag = norm(cross(x.view(i), z.view(j))).all(brille::cmp::gt,1e-10);
        // if both or neither parallel is ok (pok) and zⱼ∥xᵢ, xⱼ=zᵢ×zⱼ; otherwise xⱼ=xᵢ
        // x.set(j, (pok^u_parallel_v) ? z.cross(i,j) : x.view(i));
        x.set(j, (pok||flag) ? x.view(i) : z.cross(i,j));
      }
      if (!brille::approx::scalar(x.dot(i,j), 0.) && x.dot(i,j)<0) x.set(j, -x.extract(j));
      handled[j] = true;
    }
  }

  auto y = cross(z, x); // complete a right-handed coordinate system
  x= x/norm(x);
  y= y/norm(y);
  // debug_update("    z                x                        y           ");
  // debug_update("--------- -----------------------  -----------------------";)
  // debug_update(z.append(1,x).append(1,y).to_string()); // z, x&y dont have the same type so can't be appended :(

  /* Each symmetry operation in rotps is guaranteed to be a proper rotation,
     but there is no guarantee that it represents a *right handed* rotation.
     For the normal vectors to point the right way, we need to know if it is. */
  std::vector<bool> is_right_handed(z.size(0), true);
  auto Rv = y.extract(0); // to ensure we have the same lattice, make a copy
  for (size_t i=0; i<z.size(0); ++i) if (rotps.order(i)>2){
    brille::utils::multiply_matrix_vector(Rv.ptr(0), rotps.data(i), x.ptr(i));
    if (dot(y.view(i), Rv).all(brille::cmp::lt, 0.))
      is_right_handed[i] = false;
  }
  debug_update("Right handed:", is_right_handed);

  /* We now have for every symmetry operation a consistent vector in the
     rotation plane. We can now use z, x, and rotps to find and add the wedge
     normals. We should find one normal for each 2-fold axis and two normals
     for each 3-, 4-, and 6-fold axis. Plus we need space for one extra normal
     if there is only one rotation axis and the pointgroup has inversion.     */
  size_t found=0, expected = (rotps.size()==1) ? 1 : 0;
  for (size_t i=0; i<rotps.size(); ++i) expected += (rotps.order(i)>2) ? 2 : 1;
  LQVec<double> normals(this->outerlattice, expected);
  if (rotps.size() == 1){
    // We are assuming the system has inversion symmetry. We handle the case
    // where it doesn't elsewhere.
    this->wedge_normal_check(z.view(0), normals, found);
    if (found < 1)
      throw std::runtime_error("About to view normals 0:0 which is not allowed");
    this->ir_wedge_is_ok(normals.view(0,found)); // updates this->wedge_normals
  }
  bool accepted;
  int order;
  LQVec<double> vi(this->outerlattice, max_order), zi(this->outerlattice, 1u);
  LQVec<double> vxz, zxv;
  for (size_t i=0; i<rotps.size(); ++i){
    zi = is_right_handed[i] ? z.view(i) : -z.extract(i);
    order = rotps.order(i);
    debug_update("\nOrder ",order,", z=",zi.to_string());
    accepted = false;
    vi.set(0, x.view(i)); // do something better here?
    for (int j=1; j<order; ++j) brille::utils::multiply_matrix_vector(vi.ptr(j), rotps.data(i), vi.ptr(j-1));
    zxv = cross(zi, vi.view(0,order)); // order is guaranteed > 0
    zxv /= norm(zxv);
    debug_update("          R^n v                 z x (R^n v)      ");
    debug_update("------------------------ ------------------------");
    debug_update( vi.append(1,zxv).to_string() );
    if (2==order){
      // one-normal version of wedge_normal_check allows for either ±n
      // → no need to check z×Rv = z×(-x) = -(z×x)
      accepted = this->wedge_normal_check(zxv.view(0), normals, found);
      /* if we couldn't add a 2-fold normal, we have more work to do. */
      if (!accepted){
        // for now, hope that 90° away is good enough
        this->wedge_normal_check(vi.view(0), normals, found);
      }
    } else {
      /* Consecutive acceptable normals *must* point into the irreducible wedge
         otherwise they destroy the polyhedron.
         Check between each pair of in-plane vectors (Rⁿx, Rⁿ⁺¹x),
         for n=[0,order) with the last check (Rᵒ⁻¹x,R⁰x)≡(Rᵒ⁻¹x,x)            */
      // if (this->isinside_wedge(zi, /*constructing=*/false).getvalue(0)){
        for (int j=0; j<order; ++j){
          if (accepted) break;
          accepted = this->wedge_normal_check(zxv.view(j), -zxv.extract((j+1)%order), normals, found);
        }
      // } // isinside_wedge
    } // order>2
  }
  if (found < 1)
    throw std::runtime_error("about to view normals 0:0, which is not allowed");
  this->ir_wedge_is_ok(normals.view(0,found));
  debug_update("wedge_search finished");
}

bArray<bool> keep_if(const LQVec<double>& normals, const LQVec<double>& points){
  // determine whether we should keep points based on the provided normals.
  // We want to only keep those points in the positive half-space for any one
  // normal, but if all normals contribute to reduce the remaining points to
  // lie in a single plane we instead want to keep all points.
  std::vector<size_t> nop(normals.size(0), 0); // number of not-on-plane points
  for (size_t i=0; i<normals.size(0); ++i)
    nop[i] = dot(normals.view(i), points).is(brille::cmp::gt,0.).count();
  // If there are no planes with 0 off-plane points, divide the space
  bArray<bool> keep(points.size(0), 1u, true);
  if (std::find(nop.begin(),nop.end(),0u)==nop.end())
    for (size_t i=0; i<points.size(0); ++i)
      keep[i] = dot(normals, points.view(i)).all(brille::cmp::ge,0.);
  return keep;
}

bool BrillouinZone::wedge_brute_force(const bool special_2_folds, const bool special_mirrors, const bool sort_by_length, const bool sort_one_sym){
  std::string pmsg = "BrillouinZone::wedge_brute_force(";
  if (special_2_folds) pmsg += "2-folds";
  if (special_2_folds && special_mirrors) pmsg += ",";
  if (special_mirrors) pmsg += "mirrors";
  pmsg += ")";
  profile_update("Start ",pmsg);
  debug_exec(std::string msg;)
  // Grab the pointgroup symmetry operations
  PointSymmetry fullps = this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
  // Now restrict the symmetry operations to those with order > 1.
  PointSymmetry ps = fullps.higher(1);
  // if the full pointgroup contains only 𝟙 (and -𝟙) then ps is empty and there
  // this algorithm can not be used:
  if (0 == ps.size()){
    if (fullps.size() > 1) // 𝟙 & -𝟙 → P -1 triclinic spacegroup
      this->wedge_triclinic();
    return this->check_ir_polyhedron();
  }

  // Get and combine the characteristic points of the first Brillouin zone:
  // The face centres, face corners, and mid-face-edge points.
  auto special = cat(0, this->get_points(), this->get_vertices(), this->get_half_edges());
  std::vector<size_t> perm(ps.size());
  std::iota(perm.begin(), perm.end(), 0u); // 0u, 1u, 2u,...
    // ps is sorted in increasing rotation-order order, but we want to sort
  if (sort_by_length){
      // by increasing stationary-axis length (so that, e.g, [100] is dealt with
      // before [111]).
      std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b){
          LQVec<int> lq(this->outerlattice, 2u);
          lq.set(0u, ps.axis(a));
          lq.set(1u, ps.axis(b));
          return lq.dot(0u,0u) < lq.dot(1u,1u);
      });
  } else {
      // by decreasing isometry so that we deal with rotoinversions last
      std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b){
        return ps.isometry(a) > ps.isometry(b);
      });
  }
  ps.permute(perm); // ps now sorted

  bArray<bool> keep;

  // Keep track of *which* normals go into determining the irreducible wedge:
  // count how many normals we might need
  // 2-folds divide space with one normal
  // 3-, 4-, 6-folds need two normals
  size_t n_expected{0};
  for (size_t i=0; i<ps.size(); ++i)
    n_expected += (ps.order(i)==2) ? 1u : 2u;
  LQVec<double> cutting_normals(this->outerlattice, n_expected);
  size_t n_cut{0};

  std::vector<bool> sym_unused(ps.size(), true);
  // deal with any two-fold axes along êᵢ first:
  // The stationary vector of each rotation is a real space vector!
  LDVec<int> vec(this->outerlattice.star(), 1u);// must be int since ps.axis returns array<int,3>
  LQVec<double> nrm(this->outerlattice, 1u);
  std::vector<std::array<double,3>> eiv{{{1,0,0}},{{0,1,0}},{{0,0,1}},{{1,1,0}},{{1,-1,0}},{{1,0,1}},{{0,1,1}},{{1,0,-1}},{{0,1,-1}},{{1,1,1}}};
  auto eiv_array = bArray<double>::from_std(eiv);
  LQVec<double> eis(this->outerlattice, eiv_array);
  LDVec<double> reis(this->outerlattice.star(), eiv_array);
  size_t is_nth_ei;

  debug_update_if(special_2_folds,"Deal with 2-fold rotations with axes along highest-symmetry directions first");
  if (special_2_folds) for (size_t i=0; i<ps.size(); ++i) if (ps.isometry(i)==2){
    vec.set(0, ps.axis(i));
    // First check if this stationary axis is along a reciprocal space vector
    is_nth_ei = norm(cross(eis, vec.star())).is(brille::cmp::eq, 0.).first();
    if (is_nth_ei < 9 /* This is less than great practice */){
      debug_update("2-fold axis ",i," is ei* No. ",is_nth_ei);
      size_t e1, e2;
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
      if (norm(cross(eis, nrm)).is(brille::cmp::eq, 0.).count() == 1){
        // keep any special points beyond the bounding plane
        keep = dot(nrm, special).is(brille::cmp::ge, 0.);
        debug_update("Keeping special points with ",nrm.to_string(0)," dot p >= 0:\n",special.to_string(),keep.to_string());
        special = special.extract(keep);
        sym_unused[i] = false;
        cutting_normals.set(n_cut++, nrm);
      }
    }
    // Stationary axis along real space basis vector
    is_nth_ei = norm(cross(reis, vec)).is(brille::cmp::eq, 0.).first();
    if (sym_unused[i] && is_nth_ei < 2){
      debug_update("2-fold axis ",i," is ei No. ",is_nth_ei);
      switch (is_nth_ei){
        case 0: nrm.set(0, eiv[1]); break; /* (100) → n = (010)* */
        case 1: nrm.set(0, eiv[0]); break; /* (010) → n = (100)* */
      }
      // keep any special points beyond the bounding plane
      keep = dot(nrm, special).is(brille::cmp::ge, 0.);
      debug_update("Keeping special points with ",nrm.to_string(0)," dot p >= 0:\n", special.to_string(), keep.to_string());
      special = special.extract(keep);
      sym_unused[i] = false;
      cutting_normals.set(n_cut++, nrm);
    }
  }
  debug_update_if(special_mirrors,"Deal with mirror planes");
  if (special_mirrors) for (size_t i=0; i<ps.size(); ++i) if (ps.isometry(i)==-2){
    vec.set(0, ps.axis(i)); // the mirror plane normal is in the direct lattice
    nrm.set(0, vec.star()); // and we want the normal in the reciprocal lattice
    keep = dot(nrm, special).is(brille::cmp::ge, 0.);
    // we need at least three points (plus Γ) to define a polyhedron
    // If we are not keeping three points, check if applying the mirror plane
    // pointing the other way works for us:
    if (keep.count() < 3){
      nrm = -1*nrm; // - change nrm since we save it for later
      keep = dot(nrm, special).is(brille::cmp::ge, 0.);
    }
    if (keep.count() > 2){
      debug_update("Keeping special points with\n",nrm.to_string()," dot p >=0:\n", special.to_string(), keep.to_string());
      special = special.extract(keep);
      sym_unused[i] = false;
      cutting_normals.set(n_cut++, nrm);
    }
  }
  debug_update("Now figure out how all special points are related for each symmetry operation");
  // Find which points are mapped onto equivalent points by each symmetry operation
  std::vector<std::vector<size_t>> one_sym;
  std::vector<size_t> one_type, type_order;
  std::vector<bool> unfound(special.size(0), true), type_unfound;
  for (size_t i=0; i<ps.size(); ++i) if (sym_unused[i]){
    debug_update("Unused symmetry ",i," with order ", ps.order(i));
    debug_update(ps.get(i));
    one_sym.clear();
    for (auto b: unfound) b = true;
    for (size_t j=0; j<special.size(0); ++j) if (unfound[j]){
      one_type.clear();
      one_type.push_back(j);
      unfound[j] = false;
      for (size_t k=j+1; k<special.size(0); ++k)
      if (unfound[k] && special.match(k, j, transpose(ps.get(i)), -ps.order(i))){ // -order checks all possible rotations
      // if (unfound[k] && special.match(k, j, transpose(ps.get_proper(i)), -ps.order(i))){ // -order checks all possible rotations
        one_type.push_back(k);
        unfound[k] = false;
      }
      debug_update("Point equivalent to ",j," for symmetry ",i,":",one_type);
      // sort the equivalent points by their relative order for this operation
      // such that Rⁱj ≡ type_order[i]
      type_order.clear();
      type_order.insert(type_order.begin(), ps.order(i), special.size(0)); // set default to (non-indexable) size of the special array
      type_unfound.clear(); type_unfound.insert(type_unfound.begin(), one_type.size(), true);
      for (int o=0; o<ps.order(i); ++o)
      for (size_t k=0; k<one_type.size(); ++k) if (type_unfound[k]){
        if (special.match(one_type[k], j, transpose(ps.get(i)), o)){
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
        [&](std::vector<size_t>& a, std::vector<size_t>& b){
          size_t as = a.size() - std::count(a.begin(), a.end(), special.size(0));
          size_t bs = b.size() - std::count(b.begin(), b.end(), special.size(0));
          return as > bs;
        }
      );
    }
    size_t keep_count;
    LQVec<double> pt0(this->outerlattice, 1u), pt1(this->outerlattice, 1u);
    for (size_t s=0; s<one_sym.size(); ++s) if (sym_unused[i]/*always true?*/){
      // we need at least two equivalent points, ideally there will be the same number as the order
      if (one_sym[s].size() > static_cast<size_t>(std::count(one_sym[s].begin(), one_sym[s].end(), special.size(0))+1)){
        type_order = one_sym[s];
        debug_update("Highest-multiplicity distinct point type:",type_order);
      // auto loc = std::find_if(one_sym.begin(), one_sym.end(), [&](std::vector<size_t>x){return x.size()==static_cast<size_t>(ps.order(i));});
      // if (loc != one_sym.end()){
      //   size_t idx = loc - one_sym.begin();
        // type_order = one_sym[idx];
        // grab the rotation stationary axis
        vec.set(0, ps.axis(i));
        for (size_t j=0, k=1; j<type_order.size(); ++j, k=(j+1)%type_order.size())
        if (sym_unused[i] && type_order[j]<special.size(0) && type_order[k]<special.size(0)){
          if (ps.order(i)>2){
            // hold the two special points in their own LQVec<double> (!not <int>!)
            pt0 = special.view(type_order[j]);
            pt1 = special.view(type_order[k]);
            // we have two plane normals to worry about:
            debug_update("Stationary vector",vec.to_string(0)," and special points",pt0.to_string(0)," and",pt1.to_string(0));
            // find both cross products, remembering that we want normals
            // pointing *into* the wedge.
            nrm.resize(2);
            if ( dot(pt1, cross(vec.star(), pt0)).all(brille::cmp::lt,0.) ){
              // the rotation is left handed, so swap the special points
              nrm.set(0, cross(vec.star(), pt1));
              nrm.set(1, cross(pt0, vec.star()));
            } else {
              nrm.set(0, cross(vec.star(), pt0));
              nrm.set(1, cross(pt1, vec.star()));
            }
            debug_update("give normals:", nrm.to_string(0), " and", nrm.to_string(1));
            // now check that all special points are inside of the wedge defined by the normals
          } else {
            // order == 2, so only one normal to worry about:
            nrm = cross(vec.star(), special.view(type_order[j]));
            // make sure we don't remove all points out of the plane containing
            // the rotation axis and the two special points
            if (dot(nrm, special).is(brille::cmp::gt,0.).count() == 0)
              nrm *= -1; // switch the cross product
          }
          // check to make sure that using these normals do not remove all but a plane of points
          keep = keep_if(nrm, special);
          keep_count = keep.count();
          // We need at least three points (plus Γ) to define a polyhedron.
          // Also skip the extraction if we are keeping all points
          if (keep_count > 2 && keep_count < keep.size(0)){
            debug_update("Keeping ",keep.to_string()," of special points:\n",special.to_string());
            special = special.extract(keep);
            sym_unused[i]=false;
            for (size_t nc=0; nc<nrm.size(0); ++nc)
                cutting_normals.set(n_cut++, nrm.view(nc));
          }
        }
      }
    }
  }

  // debug_update("Remaining special points\n", special.to_string());
  bArray<double> cn(0,3);
  if (n_cut > 0) /*protect against view(0,0)*/
    cn = cutting_normals.view(0,n_cut).get_xyz(); // the cutting direction is opposite the normal
  bArray<double> cp(cn.shape(), 0.);
  this->ir_polyhedron = Polyhedron::bisect(this->polyhedron, -1*cn, cp);
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
  if (fullps.has_space_inversion()){
    // ps only contains operations *with* a stationary axis (described in real space)
    LDVec<int> all_axes(this->outerlattice.star(), bArray<int>::from_std(ps.axes()));
    double ir_volume = this->ir_polyhedron.get_volume();
    double goal_volume = this->polyhedron.get_volume()/static_cast<double>(fullps.size());
    auto aaiu = all_axes.is_unique();
    if (std::count(aaiu.begin(), aaiu.end(), true) == 1u || brille::approx::scalar(ir_volume, 2.0*goal_volume) ){
      debug_update("Deal with -1 since there is only one stationary axis (or doubled volume for some other reason)");
      Polyhedron div;
      auto bz_n = this->get_normals();
      // auto wg_n = this->get_ir_wedge_normals(); // !! This is empty if the thus-far-found volume is wrong
      auto wg_n = this->get_ir_polyhedron_wedge_normals(); // This pulls directly from the ir_polyhedron
      bArray<double> gamma({1u,3u}, 0.);
      for (size_t i=0; i<bz_n.size(0); ++i){
        div = this->ir_polyhedron.divide(bz_n.view(i).get_xyz(), gamma);
        if (brille::approx::scalar(div.get_volume(), goal_volume)){
          // set div to be the ir_polyhedron
          this->ir_polyhedron = div;
          // add the new normal to wedge normals list
          // extract instead of view to avoid copying the whole bz_n Array
          // when the inversion happens.
          wg_n.append(/*dimension to expand*/ 0u, -bz_n.extract(i));
          this->set_ir_wedge_normals(wg_n);
          break;
        }
      }
      if (brille::approx::scalar(this->ir_polyhedron.get_volume(), ir_volume)){
        debug_update("Polyhedron volume still double expected.");
        bool proceed=true;
        // check for other dividing planes :/
        //LQVec<double> cij(this->outerlattice, 1u);
        for (size_t i=0; proceed && i<bz_n.size(0)-1; ++i)
        for (size_t j=i+1; proceed && j<bz_n.size(0); ++j){
          div = this->ir_polyhedron.divide(bz_n.cross(i,j).get_xyz(), gamma);
          if (brille::approx::scalar(div.get_volume(), goal_volume)){
            this->ir_polyhedron = div;
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


void BrillouinZone::wedge_triclinic(void){
  /* Assuming that this is spacegroup P -1 we have inversion symmetry only.
     We always want to find a set of planes that divides symmetry equivalent
     regions of space. With only -1 symmetry we have to make a choice as to
     *which* space inversion we care about.
     We could restrict ourselves to {x≥0, y, z}, {x, y≥0, z}, {x, y, z≥0} or
     their opposites, or any other subspace for which xyz ≥ 0.
     Equivalently then, we can resstrict ourselves to the subspace where
     x̂⋅(111) ≥ 0.
  */
  using namespace brille;
  bArray<double> vec(1u, 3u, 1.);
  LQVec<double> nrm(this->outerlattice, vec);
  this->set_ir_wedge_normals(nrm);
  this->irreducible_vertex_search();
}
