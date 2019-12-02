/* Copyright 2019 Greg Tucker
//
// This file is part of brille.
//
// brille is free software: you can redistribute it and/or modify it under the
// terms of the GNU Affero General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// brille is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with brille. If not, see <https://www.gnu.org/licenses/>.            */

#include "bz.hpp"

void BrillouinZone::wedge_search(const bool pbv, const bool pok){
  debug_exec(std::string update_msg;)
  // Get the full pointgroup symmetry information
  PointSymmetry fullps = this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
  // And use it to find only the highest-order rotation operation along each
  // unique stationary axis.
  // PointSymmetry rotps = fullps.nfolds(1); // 1 to request only orders>1
  PointSymmetry rotps = fullps.higher(1); // 1 to request only orders>1
  // Get the vectors pointing to each full Brillouin zone facet cetre
  LQVec<double> xyz = cat(cat(this->get_points(), this->get_vertices()), this->get_half_edges());

  // debug_update("xyz=\n", xyz.to_string());

  // if rotps is empty, there are no 2+-fold rotations, act like we have ̄1:
  if (rotps.size()==0){
    debug_update("No 2+-fold rotation operations");
    this->ir_wedge_is_ok(xyz.first(1u)); // for ̄1, assigns this->ir_polyhedron
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
  LQVec<int> z(this->outerlattice, rotps.axes());
  // Find ̂zᵢ⋅ẑⱼ
  std::vector<std::vector<int>> dotij(z.size());
  double dottmp;
  for (size_t i=0; i<z.size(); ++i) for (size_t j=0; j<z.size(); ++j){
    dottmp = std::abs(z.dot(i,j)/z.norm(i)/z.norm(j));
    dotij[i].push_back( approx_scalar(dottmp,0.) ? -1: approx_scalar(dottmp, 1.) ? 1 : 0);
  }
  debug_update("dot(i,j) flag\n", dotij);

  // Find a suitable in-rotation-plane vector for each stationary axis
  LQVec<double> x(this->outerlattice, rotps.perpendicular_axes());
  // or pick one of the BZ facet-points if the basis vector preffered flag
  if (pbv) for (size_t i=0; i<x.size(); ++i) for (size_t j=0; j<xyz.size(); ++j)
  if (dot(xyz.get(j), z.get(i)).all_approx(0.)){
    x.set(i, xyz.get(j));
    break;
  }

  /* Two rotations with ẑᵢ⋅ẑⱼ≠0 should have (Rⁿx̂ᵢ)⋅(Rᵐx̂ⱼ)=1 for some n,m.
     To start, assume that n and m are 0 and find the x̂ᵢ such that x̂ᵢ⋅(ẑᵢ×ẑⱼ)=1
  */
  std::vector<bool> handled(z.size(), false);
  bool flag;
  for (size_t i=0; i<z.size()-1; ++i) for (size_t j=i+1; j<z.size(); ++j)
  // if zᵢ and zⱼ are neither parallel or perpendicular
  if (0 == dotij[i][j]){
    if (!handled[i]){
      if (!pbv /*basis vectors not preferred*/) x.set(i, z.cross(i, j));
      if (handled[j] && !approx_scalar(x.dot(i,j), 0.) && x.dot(i,j)<0) x.set(i, -x.get(i));
      handled[i] = true;
    }
    if (!handled[j]){
      if (!pbv /*basis vectors not preferred*/){
        flag = norm(cross(x.get(i), z.get(j))).all_approx(Comp::gt,1e-10);
        // if both or neither parallel is ok (pok) and zⱼ∥xᵢ, xⱼ=zᵢ×zⱼ; otherwise xⱼ=xᵢ
        // x.set(j, (pok^u_parallel_v) ? z.cross(i,j) : x.get(i));
        x.set(j, (pok||flag) ? x.get(i) : z.cross(i,j));
      }
      if (!approx_scalar(x.dot(i,j), 0.) && x.dot(i,j)<0) x.set(j, -x.get(j));
      handled[j] = true;
    }
  }

  LQVec<double> y = cross(z, x); // complete a right-handed coordinate system
  x= x/norm(x);
  y= y/norm(y);
  debug_update("    z                x                        y           ");
  debug_exec(update_msg = "--------- -----------------------  -----------------------\n";)
  for (size_t i=0; i<z.size(); ++i){
    debug_exec(update_msg += z.to_string(i) +" "+  x.to_string(i) +" "+ y.to_string(i) +"\n";)
  }
  debug_update(update_msg);

  /* Each symmetry operation in rotps is guaranteed to be a proper rotation,
     but there is no guarantee that it represents a *right handed* rotation.
     For the normal vectors to point the right way, we need to know if it is. */
  std::vector<bool> is_right_handed(z.size(), true);
  LQVec<double> Rv(this->outerlattice, 1u);
  for (size_t i=0; i<z.size(); ++i) if (rotps.order(i)>2){
    multiply_matrix_vector(Rv.data(0), rotps.data(i), x.data(i));
    if (dot(y.get(i), Rv.get(0)).getvalue(0) < 0)
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
    this->wedge_normal_check(z.get(0), normals, found);
    this->ir_wedge_is_ok(normals.first(found)); // updates this->wedge_normals
  }
  bool accepted;
  int order;
  LQVec<double> vi(this->outerlattice, max_order), zi(this->outerlattice, 1u);
  LQVec<double> vxz, zxv;
  for (size_t i=0; i<rotps.size(); ++i){
    zi = is_right_handed[i] ? z.get(i) : -z.get(i);
    order = rotps.order(i);
    debug_update("\nOrder ", order, ", z=", zi.to_string(0));
    accepted = false;
    vi.set(0, x.get(i)); // do something better here?
    for (int j=1; j<order; ++j) multiply_matrix_vector(vi.data(j), rotps.data(i), vi.data(j-1));
    zxv = cross(zi, vi.first(order));
    zxv = zxv/norm(zxv);
    debug_exec(update_msg ="          R^n v                 z x (R^n v)      \n";)
    debug_exec(update_msg+="------------------------ ------------------------\n";)
    for (int j=0; j<order; ++j){
      debug_exec(update_msg += vi.to_string(j) +" "+ zxv.to_string(j)+"\n");
    }
    debug_update(update_msg);
    if (2==order){
      // one-normal version of wedge_normal_check allows for either ±n
      // → no need to check z×Rv = z×(-x) = -(z×x)
      accepted = this->wedge_normal_check(zxv.get(0), normals, found);
      /* if we couldn't add a 2-fold normal, we have more work to do. */
      if (!accepted){
        // for now, hope that 90° away is good enough
        this->wedge_normal_check(vi.get(0), normals, found);
      }
    } else {
      /* Consecutive acceptable normals *must* point into the irreducible wedge
         otherwise they destroy the polyhedron.
         Check between each pair of in-plane vectors (Rⁿx, Rⁿ⁺¹x),
         for n=[0,order) with the last check (Rᵒ⁻¹x,R⁰x)≡(Rᵒ⁻¹x,x)            */
      // if (this->isinside_wedge(zi, /*constructing=*/false).getvalue(0)){
        for (int j=0; j<order; ++j){
          if (accepted) break;
          accepted = this->wedge_normal_check(zxv.get(j), -zxv.get((j+1)%order), normals, found);
        }
      // } // isinside_wedge
    } // order>2
  }
  this->ir_wedge_is_ok(normals.first(found));
  debug_update("wedge_search finished");
}


ArrayVector<bool> keep_if(const LQVec<double>& normals, const LQVec<double>& points){
  // determine whether we should keep points based on the provided normals.
  // We want to only keep those points in the positive half-space for any one
  // normal, but if all normals contribute to reduce the remaining points to
  // lie in a single plane we instead want to keep all points.
  ArrayVector<bool> keep(1u, points.size(), true);
  std::vector<size_t> nop(normals.size(), 0); // number of not-on-plane points
  for (size_t i=0; i<normals.size(); ++i)
    nop[i] = dot(normals.extract(i), points).is_approx(Comp::gt,0.).count_true();
  // If there are no planes with 0 off-plane points, divide the space
  if (std::find(nop.begin(),nop.end(),0u)==nop.end())
    for (size_t i=0; i<points.size(); ++i)
      keep.insert(dot(normals, points.extract(i)).all_approx(Comp::ge,0.), i);
  return keep;
}

void BrillouinZone::wedge_brute_force(const bool special_2_folds, const bool special_mirrors, const bool sort_by_length){
  debug_exec(std::string msg;)
  Spacegroup sg = this->outerlattice.get_spacegroup_object();
  // For P1 the irreducible Brillouin zone *is* the first Brillouin zone
  // which has already been set
  if (/*P 1*/sg.get_international_table_number() < 2) return;
  // The other triclinic spacegroup, P -1, is too complicated for this algorithm
  if (/*P 1 or P-1*/sg.get_international_table_number() < 3){
    this->wedge_triclinic();
    return;
  }
  // Get and combine the characteristic points of the first Brillouin zone:
  // The face centres, face corners, and mid-face-edge points.
  LQVec<double> special = cat(this->get_points(), this->get_vertices(), this->get_half_edges());

  // Grab the pointgroup symmetry operations
  PointSymmetry fullps = this->outerlattice.get_pointgroup_symmetry(this->time_reversal);
  // Now restrict the symmetry operations to those with order > 1.
  PointSymmetry ps = fullps.higher(1);
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

  ArrayVector<bool> keep(1u, 0u);

  // Keep track of *which* normals go into determining the irreducible wedge:
  LQVec<double> cutting_normals(this->outerlattice, ps.size());
  size_t n_cut{0};

  std::vector<bool> sym_unused(ps.size(), true);
  // deal with any two-fold axes along êᵢ first:
  // The stationary vector of each rotation is a real space vector!
  LDVec<int> vec(this->outerlattice.star(), 1u);// must be int since ps.axis returns array<int,3>
  LQVec<double> nrm(this->outerlattice, 1u);
  std::vector<std::array<double,3>> eiv{{{1,0,0}},{{0,1,0}},{{0,0,1}},{{1,1,0}},{{1,-1,0}},{{1,0,1}},{{0,1,1}},{{1,0,-1}},{{0,1,-1}},{{1,1,1}}};
  LQVec<double> eis(this->outerlattice, eiv.size());
  for (size_t i=0; i<eis.size(); ++i) eis.set(i, eiv[i]);
  LDVec<double> reis(this->outerlattice.star(), eiv.size());
  for (size_t i=0; i<reis.size(); ++i) reis.set(i, eiv[i]);
  size_t is_nth_ei;

  debug_update_if(special_2_folds,"Deal with 2-fold rotations with axes along highest-symmetry directions first");
  if (special_2_folds) for (size_t i=0; i<ps.size(); ++i) if (ps.isometry(i)==2){
    vec.set(0, ps.axis(i));
    // First check if this stationary axis is along a reciprocal space vector
    is_nth_ei = norm(cross(eis, vec.star())).is_approx(Comp::eq, 0.).first_true();
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
      nrm.set(0, cross(reis.extract(e1).star(), eis.extract(e2)));
      if (norm(cross(eis, nrm)).is_approx(Comp::eq, 0.).count_true() == 1){
        // keep any special points beyond the bounding plane
        keep = dot(nrm, special).is_approx(Comp::ge, 0.);
        debug_update("Keeping special points with\n",nrm.to_string(0)," dot p >= 0:\n", special.to_string(keep));
        special = special.extract(keep);
        sym_unused[i] = false;
        cutting_normals.set(n_cut++, nrm);
      }
    }
    // Stationary axis along real space basis vector
    is_nth_ei = norm(cross(reis, vec)).is_approx(Comp::eq, 0.).first_true();
    if (sym_unused[i] && is_nth_ei < 2){
      debug_update("2-fold axis ",i," is ei No. ",is_nth_ei);
      switch (is_nth_ei){
        case 0: nrm.set(0, eiv[1]); break; /* (100) → n = (010)* */
        case 1: nrm.set(0, eiv[0]); break; /* (010) → n = (100)* */
      }
      // keep any special points beyond the bounding plane
      keep = dot(nrm, special).is_approx(Comp::ge, 0.);
      debug_update("Keeping special points with\n",nrm.to_string(0)," dot p >= 0:\n", special.to_string(keep));
      special = special.extract(keep);
      sym_unused[i] = false;
      cutting_normals.set(n_cut++, nrm);
    }
  }
  debug_update_if(special_mirrors,"Deal with mirror planes");
  if (special_mirrors) for (size_t i=0; i<ps.size(); ++i) if (ps.isometry(i)==-2){
    vec.set(0, ps.axis(i)); // the mirror plane normal is in the direct lattice
    nrm.set(0, vec.star()); // and we want the normal in the reciprocal lattice
    keep = dot(nrm, special).is_approx(Comp::ge, 0.);
    // we need at least three points (plus Γ) to define a polyhedron
    // If we are not keeping three points, check if applying the mirror plane
    // pointing the other way works for us:
    if (keep.count_true() < 3){
      nrm = -1*nrm; // - change nrm since we save it for later
      keep = dot(nrm, special).is_approx(Comp::ge, 0.);
    }
    if (keep.count_true() > 2){
      debug_update("Keeping special points with\n",nrm.to_string(0)," dot p >=0:\n", special.to_string(keep));
      special = special.extract(keep);
      sym_unused[i] = false;
      cutting_normals.set(n_cut++, nrm);
    }
  }
  debug_update("Now figure out how all special points are related for each symmetry operation");
  // Find which points are mapped onto equivalent points by each symmetry operation
  std::vector<std::vector<size_t>> one_sym;
  std::vector<size_t> one_type, type_order;
  std::vector<bool> unfound(special.size(), true), type_unfound;
  for (size_t i=0; i<ps.size(); ++i) if (sym_unused[i]){
    debug_update("Unused symmetry ",i," with order ", ps.order(i));
    debug_update(ps.get(i));
    one_sym.clear();
    for (auto b: unfound) b = true;
    for (size_t j=0; j<special.size(); ++j) if (unfound[j]){
      one_type.clear();
      one_type.push_back(j);
      unfound[j] = false;
      for (size_t k=j+1; k<special.size(); ++k)
      if (unfound[k] && special.rotate_approx(k, j, transpose(ps.get(i)), -ps.order(i))){ // -order checks all possible rotations
      // if (unfound[k] && special.rotate_approx(k, j, transpose(ps.get_proper(i)), -ps.order(i))){ // -order checks all possible rotations
        one_type.push_back(k);
        unfound[k] = false;
      }
      debug_update("Point equivalent to ",j," for symmetry ",i,":",one_type);
      // sort the equivalent points by their relative order for this operation
      // such that Rⁱj ≡ type_order[i]
      type_order.clear();
      type_order.insert(type_order.begin(), ps.order(i), special.size()); // set default to (non-indexable) size of the special array
      type_unfound.clear(); type_unfound.insert(type_unfound.begin(), one_type.size(), true);
      for (int o=0; o<ps.order(i); ++o)
      for (size_t k=0; k<one_type.size(); ++k) if (type_unfound[k]){
        if (special.rotate_approx(one_type[k], j, transpose(ps.get(i)), o)){
        // if (special.rotate_approx(one_type[k], j, transpose(ps.get_proper(i)), o)){
          type_order[o] = one_type[k];
          type_unfound[k] = false;
        }
      }
      // and store the sorted equivalent indices
      debug_update("Which are sorted by their rotation order:",type_order);
      one_sym.push_back(type_order);
    }
    // sort one_sym by the number of valid (indexable) equivalent points
    std::sort(one_sym.begin(), one_sym.end(),
      [&](std::vector<size_t>& a, std::vector<size_t>& b){
        size_t as = a.size() - std::count(a.begin(), a.end(), special.size());
        size_t bs = b.size() - std::count(b.begin(), b.end(), special.size());
        return as > bs;
      }
    );
    size_t keep_count;
    LQVec<double> pt0(this->outerlattice, 1u), pt1(this->outerlattice, 1u);
    for (size_t s=0; s<one_sym.size(); ++s) if (sym_unused[i]/*always true?*/){
      // we need at least two equivalent points, ideally there will be the same number as the order
      if (one_sym[s].size() > static_cast<size_t>(std::count(one_sym[s].begin(), one_sym[s].end(), special.size())+1)){
        type_order = one_sym[s];
        debug_update("Highest-multiplicity distinct point type:",type_order);
      // auto loc = std::find_if(one_sym.begin(), one_sym.end(), [&](std::vector<size_t>x){return x.size()==static_cast<size_t>(ps.order(i));});
      // if (loc != one_sym.end()){
      //   size_t idx = loc - one_sym.begin();
        // type_order = one_sym[idx];
        // grab the rotation stationary axis
        vec.set(0, ps.axis(i));
        for (size_t j=0, k=1; j<type_order.size(); ++j, k=(j+1)%type_order.size())
        if (sym_unused[i] && type_order[j]<special.size() && type_order[k]<special.size()){
          if (ps.order(i)>2){
            // hold the two special points in their own LQVec<double> (!not <int>!)
            pt0 = special.extract(type_order[j]);
            pt1 = special.extract(type_order[k]);
            // we have two plane normals to worry about:
            debug_update("Stationary vector",vec.to_string(0),"and special points",pt0.to_string(0),"and",pt1.to_string(0));
            // find both cross products, remembering that we want normals
            // pointing *into* the wedge.
            nrm.resize(2);
            if ( dot(pt1, cross(vec.star(), pt0)).all_approx(Comp::lt,0.) ){
              // the rotation is left handed, so swap the special points
              nrm.set(0, cross(vec.star(), pt1));
              nrm.set(1, cross(pt0, vec.star()));
            } else {
              nrm.set(0, cross(vec.star(), pt0));
              nrm.set(1, cross(pt1, vec.star()));
            }
            debug_update("give normals:", nrm.to_string(";"));
            // now check that all special points are inside of the wedge defined by the normals
          } else {
            // order == 2, so only one normal to worry about:
            nrm = cross(vec.star(), special.extract(type_order[j]));
            // make sure we don't remove all points out of the plane containing
            // the rotation axis and the two special points
            if (dot(nrm, special).is_approx(Comp::gt,0.).count_true() == 0)
              nrm *= -1; // switch the cross product
          }
          // check to make sure that using these normals do not remove all but a plane of points
          keep = keep_if(nrm, special);
          keep_count = keep.count_true();
          // We need at least three points (plus Γ) to define a polyhedron.
          // Also skip the extraction if we are keeping all points
          if (keep_count > 2 && keep_count < keep.size()){
            debug_update("Keeping special points:\n",special.to_string(keep));
            special = special.extract(keep);
            sym_unused[i]=false;
            for (size_t nc=0; nc<nrm.size(); ++nc)
                cutting_normals.set(n_cut++, nrm.extract(nc));
          }
        }
      }
    }
  }

  // debug_update("Remaining special points\n", special.to_string());
  ArrayVector<double> cn = cutting_normals.first(n_cut).get_xyz(); // the cutting direction is opposite the normal
  ArrayVector<double> cp(3u, cn.size(), 0.);
  this->ir_polyhedron = Polyhedron::bisect(this->polyhedron, -1*cn, cp);
  // copy functionality of set_ir_vertices, which set the normals as well
  if (this->check_ir_polyhedron())
    this->set_ir_wedge_normals(this->get_ir_polyhedron_wedge_normals());

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
    LDVec<int> all_axes(this->outerlattice.star(), ps.axes());
    double ir_volume = this->ir_polyhedron.get_volume();
    double goal_volume = this->polyhedron.get_volume()/static_cast<double>(fullps.size());
    if (all_axes.is_unique().count_true() == 1u || approx_scalar(ir_volume, 2.0*goal_volume) ){
      debug_update("Deal with -1 since there is only one stationary axis (or doubled volume for some other reason)");
      Polyhedron div;
      LQVec<double> bz_n = this->get_normals();
      // LQVec<double> wg_n = this->get_ir_wedge_normals(); // !! This is empty if the thus-far-found volume is wrong
      LQVec<double> wg_n = this->get_ir_polyhedron_wedge_normals(); // This pulls directly from the ir_polyhedron
      ArrayVector<double> gamma(3u, 1u);
      for (size_t i=0; i<3u; ++i) gamma.insert(0, 0, i);
      for (size_t i=0; i<bz_n.size(); ++i){
        div = this->ir_polyhedron.divide(bz_n.extract(i).get_xyz(), gamma);
        if (approx_scalar(div.get_volume(), goal_volume)){
          // set div to be the ir_polyhedron
          this->ir_polyhedron = div;
          // add the new normal to wedge normals list
          wg_n.resize(wg_n.size()+1);
          wg_n.set(wg_n.size()-1, -bz_n.extract(i));
          this->set_ir_wedge_normals(wg_n);
          break;
        }
      }
      if (approx_scalar(this->ir_polyhedron.get_volume(), ir_volume)){
        debug_update("Polyhedron volume still double expected.");
        bool proceed=true;
        // check for other dividing planes :/
        LQVec<double> cij(this->outerlattice, 1u);
        for (size_t i=0; proceed && i<bz_n.size()-1; ++i)
        for (size_t j=i+1; proceed && j<bz_n.size(); ++j){
          div = this->ir_polyhedron.divide(bz_n.cross(i,j).get_xyz(), gamma);
          if (approx_scalar(div.get_volume(), goal_volume)){
            this->ir_polyhedron = div;
            wg_n.resize(wg_n.size()+1);
            wg_n.set(wg_n.size()-1, -bz_n.cross(i,j));
            this->set_ir_wedge_normals(wg_n);
            proceed = false;
          }
        }
      }
    }
  }

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
  LQVec<double> nrm(this->outerlattice, 1u);
  nrm.set(0, std::array<double,3>({{1,1,1}}));
  this->set_ir_wedge_normals(nrm);
  this->irreducible_vertex_search();
}
