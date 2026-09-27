/* This file is part of brille.

Copyright © 2026 Greg Tucker <gregory.tucker@ess.eu>

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
#ifndef BRILLE_LATTICE_MESH_HPP_
#define BRILLE_LATTICE_MESH_HPP_
/*! \file
    \brief The structured mesh of the irreducible Brillouin zone, as Mesh3 uses it

`LatticeMesh::from_zone` builds a `latticetri::LatticeTri` for a `BrillouinZone`'s
own irreducible zone, and keeps its vertices (Cartesian, Å⁻¹) and tetrahedra with
a bucket grid that finds a point's tetrahedron in constant time.
*/
#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>
#if defined(__GLIBC__)
#include <malloc.h>
#endif
#include "array_.hpp"
#include "bz.hpp"
#include "hdf_interface.hpp"
#include "lattice_tri.hpp"
#include "permutation_table.hpp"

namespace brille {

class LatticeMesh {
  using vec3 = std::array<double, 3>;
  bArray<double> positions_;                 // (vertices, 3), Cartesian
  bArray<ind_t> tetrahedra_;                 // (tetrahedra, 4)
  bool limited_{false};                      // max_points made the grid coarser than asked
  long long divisions_{0};                   // the grid lattice is Λ* / divisions_
  // location: each tetrahedron's origin and the inverse of its edge matrix, and
  // the tetrahedra whose bounding box meets each cubic bucket
  std::vector<std::array<double, 12>> inverse_;
  vec3 lo_{};
  double side_{1};
  std::array<long long, 3> dims_{1, 1, 1};
  std::vector<size_t> offsets_;
  std::vector<ind_t> members_;

public:
  //! How far, in barycentric weight, a point may be outside a tetrahedron and still be located in it
  static constexpr double surface_tolerance{1e-10};

  LatticeMesh() = default;
  LatticeMesh(bArray<double> positions, bArray<ind_t> tetrahedra, const bool limited = false, const long long divisions = 0)
      : positions_(std::move(positions)), tetrahedra_(std::move(tetrahedra)), limited_(limited), divisions_(divisions) {
    index();
  }

  //! What LatticeTri needs to mesh a zone's irreducible part, in the primitive reciprocal basis
  struct Inputs {
    std::array<double, 9> metric{};              //!< the reciprocal metric (Å⁻²)
    std::vector<latticetri::mat3i> ops;          //!< the point group, with time reversal's inversion
    std::vector<latticetri::int3> cone;          //!< the zone's wedge: c·x >= 0 inside
    std::array<double, 9> basis{};               //!< columns: the primitive reciprocal vectors (Å⁻¹)
  };

private:
  // What refinement needs: the construction's inputs, every refinement so far (the
  // marked tetrahedra, by their vertices, and the edge limit), and the LatticeTri
  // they give, rebuilt from them on first use so a mesh that is never refined
  // doesn't hold it.
  std::optional<Inputs> inputs_;
  std::vector<std::vector<std::array<ind_t, 4>>> history_marked_;
  std::vector<double> history_min_edge_;
  mutable std::shared_ptr<latticetri::LatticeTri> tri_;

public:
  //! A refinement worked out on a copy of the triangulation, ready to `commit`
  struct Refinement {
    std::shared_ptr<latticetri::LatticeTri> tri;
    std::vector<std::array<ind_t, 4>> marked;
    double min_edge{0};
    bArray<double> points;   //!< the new vertices, Cartesian (Å⁻¹)
  };
  [[nodiscard]] bool refinable() const { return inputs_.has_value(); }
  //! Whether the triangulation that refinement works on is held (rather than rebuilt when next needed)
  [[nodiscard]] bool holds_triangulation() const { return static_cast<bool>(tri_); }
  /*! \brief Free the triangulation that refinement works on

  The mesh is unchanged and can still be refined: the triangulation is rebuilt from
  what built the mesh and the refinements since, when next needed, which costs a
  build and a replay of those refinements.
  */
  void release_triangulation() { tri_.reset(); }
  /*! \brief Work out refining the tetrahedra `tets` (by index), without changing the mesh

  Each marked tetrahedron is bisected once, with closure; boundary edges are split
  with their equivalents, so paired zone faces keep matching. A tetrahedron whose
  longest edge is at most twice `min_edge` (Å⁻¹) is not split. The new vertices are
  those that `commit` would add, in the order it would add them.
  */
  [[nodiscard]] Refinement plan(const std::vector<ind_t> & tets, const double min_edge) const {
    Refinement out;
    for (const auto t: tets) {
      if (t >= number_of_tetrahedra()) throw std::out_of_range("tetrahedron index " + std::to_string(t) + " is not in the mesh");
      out.marked.push_back({tetrahedra_.val(t, 0), tetrahedra_.val(t, 1), tetrahedra_.val(t, 2), tetrahedra_.val(t, 3)});
    }
    out.min_edge = min_edge;
    out.tri = std::make_shared<latticetri::LatticeTri>(triangulation());
    out.tri->refine(as_tets(out.marked), min_edge);
    const auto & xp = out.tri->vertices();
    const auto first = static_cast<size_t>(number_of_vertices());
    out.points = bArray<double>(static_cast<ind_t>(xp.size() - first), 3u);
    for (size_t i = first; i < xp.size(); ++i) {
      const auto x = apply(inputs_->basis, xp[i]);
      for (int k = 0; k < 3; ++k) out.points.val(static_cast<ind_t>(i - first), k) = x[k];
    }
    return out;
  }
  //! Apply a refinement from `plan`: the new vertices are appended, so existing vertices keep their indices
  void commit(Refinement && r) {
    positions_ = cat(0, positions_, r.points);
    tri_ = std::move(r.tri);
    const auto & tets = tri_->tetrahedra();
    tetrahedra_ = bArray<ind_t>(static_cast<ind_t>(tets.size()), 4u);
    for (size_t i = 0; i < tets.size(); ++i) for (int k = 0; k < 4; ++k) tetrahedra_.val(static_cast<ind_t>(i), k) = static_cast<ind_t>(tets[i][k]);
    history_marked_.push_back(std::move(r.marked));
    history_min_edge_.push_back(r.min_edge);
    index();
  }

private:
  static std::vector<std::array<size_t, 4>> as_tets(const std::vector<std::array<ind_t, 4>> & marked) {
    std::vector<std::array<size_t, 4>> out;
    for (const auto & t: marked) out.push_back({t[0], t[1], t[2], t[3]});
    return out;
  }
  //! The triangulation, rebuilt from the inputs and history if not held
  const latticetri::LatticeTri & triangulation() const {
    if (!tri_) {
      if (!inputs_)
        throw std::runtime_error("This mesh can't be refined: it was read from a file without what built it (e.g. a TetGen mesh)");
      auto tri = std::make_shared<latticetri::LatticeTri>(inputs_->metric, inputs_->ops, divisions_, inputs_->cone);
      for (size_t i = 0; i < history_marked_.size(); ++i) tri->refine(as_tets(history_marked_[i]), history_min_edge_[i]);
      if (tri->vertices().size() != static_cast<size_t>(number_of_vertices()) || tri->tetrahedra().size() != static_cast<size_t>(number_of_tetrahedra()))
        throw std::runtime_error("Rebuilding the mesh for refinement gave a different mesh");
      tri_ = std::move(tri);
    }
    return *tri_;
  }

public:
  //! The primitive metric, the point group as integer matrices, and the zone's own wedge
  static Inputs inputs(const BrillouinZone & bz) {
    Inputs in;
    const auto outer = bz.get_lattice();
    const auto Bp = outer.primitive().reciprocal_basis_vectors();   // columns: primitive reciprocal vectors
    in.basis = Bp;
    const auto A = outer.real_basis_vectors();                      // columns: conventional direct vectors
    const auto Bi = inverse(Bp);
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) in.metric[3 * i + j] += Bp[3 * k + i] * Bp[3 * k + j];
    // the point group (with time reversal's inversion) on primitive reciprocal coordinates
    const auto Ai = inverse(A);
    const auto ps = bz.get_pointgroup_symmetry();
    for (size_t j = 0; j < ps.size(); ++j) {
      const auto W = ps.get(j);
      std::array<double, 9> w{};
      for (int i = 0; i < 9; ++i) w[i] = static_cast<double>(W[i]);
      // Cartesian rotation A W A⁻¹, then Bp⁻¹ (that) Bp
      const auto R = product(Bi, product(product(A, w), product(Ai, Bp)));
      latticetri::mat3i r{};
      for (int i = 0; i < 9; ++i) {
        r[i] = std::llround(R[i]);
        if (std::abs(R[i] - static_cast<double>(r[i])) > 1e-6)
          throw std::runtime_error("a point group operation is not integer on the primitive reciprocal lattice");
      }
      if (std::find(in.ops.begin(), in.ops.end(), r) == in.ops.end()) in.ops.push_back(r);
    }
    // the zone's own wedge, as integer c with c·x >= 0 inside
    const auto ir_xyz = bz.get_ir_polyhedron().vertices().xyz();
    vec3 centre{0, 0, 0};
    for (ind_t i = 0; i < ir_xyz.size(0); ++i) for (int k = 0; k < 3; ++k) centre[k] += ir_xyz.val(i, k) / static_cast<double>(ir_xyz.size(0));
    const auto centre_p = apply(Bi, centre);
    // A wedge plane passes through Γ, so its normal is that of an irreducible face
    // through Γ, computed from the face's vertices. brille's wedge normals are only
    // used to pick those faces: for zones with faces far smaller than themselves it
    // can take a sliver face on the zone boundary for part of the wedge, and a
    // sliver's normal is not a lattice direction.
    const auto gamma_faces = faces_through_gamma(bz);
    const auto normals = bz.get_ir_wedge_normals();
    if (normals.size(0)) {
      const auto xyz = normals.xyz();
      for (ind_t i = 0; i < xyz.size(0); ++i) {
        const vec3 given{xyz.val(i, 0), xyz.val(i, 1), xyz.val(i, 2)};
        const auto face = parallel_to(given, gamma_faces);
        if (!face) continue;
        vec3 c{};
        for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) c[k] += Bp[3 * l + k] * (*face)[l];   // Bpᵀ N
        auto n = integer_direction(c);
        double side{0};
        for (int k = 0; k < 3; ++k) side += static_cast<double>(n[k]) * centre_p[k];
        if (side < 0) for (auto & x: n) x = -x;
        if (std::find(in.cone.begin(), in.cone.end(), n) == in.cone.end()) in.cone.push_back(n);
      }
      if (in.cone.empty()) throw std::runtime_error("none of the zone's wedge normals is that of an irreducible face through Γ");
    }
    return in;
  }

  /*! \brief The structured mesh of a zone's irreducible part

  \param bz the zone; the mesh fills exactly its irreducible polyhedron
  \param max_volume the largest grid tetrahedron volume (Å⁻³), which sets the grid
         spacing; not positive for the coarsest grid, Λ* itself
  \param max_points if positive, the grid is coarsened until the estimated vertex
         count is at most this, and `refinement_limited` says so
  */
  static LatticeMesh from_zone(const BrillouinZone & bz, const double max_volume = -1, const int max_points = -1) {
    auto mesh = build(bz, max_volume, max_points);
#if defined(__GLIBC__)
    // Building frees many small allocations, which glibc keeps mapped (about 240 MiB
    // for 5×10⁴ vertices, four times what the mesh holds): return them
    malloc_trim(0);
#endif
    return mesh;
  }

private:
  static LatticeMesh build(const BrillouinZone & bz, const double max_volume, const int max_points) {
    const auto in = inputs(bz);
    const auto & Bp = in.basis;
    const auto ir = bz.get_ir_polyhedron();
    // the grid spacing: a Kuhn cell of Λ*/n holds six tetrahedra
    const double cell = std::abs(determinant(Bp));
    long long n{1};
    if (max_volume > 0) n = std::max(1LL, static_cast<long long>(std::ceil(std::cbrt(cell / (6 * max_volume)))));
    bool limited{false};
    if (max_points > 0) {
      const double fraction = ir.volume() / cell;
      auto estimate = [&](const long long m) {
        const double k = static_cast<double>(m);
        return fraction * k * k * k + 6 * std::cbrt(fraction) * k * k + 8;   // interior plus a boundary layer
      };
      while (n > 1 && estimate(n) > max_points) { --n; limited = true; }
    }
    latticetri::LatticeTri tri(in.metric, in.ops, n, in.cone);
    const auto & xp = tri.vertices();
    const auto & tets = tri.tetrahedra();
    bArray<double> positions(xp.size(), 3u);
    for (size_t i = 0; i < xp.size(); ++i) {
      const auto x = apply(Bp, xp[i]);
      for (int k = 0; k < 3; ++k) positions.val(i, k) = x[k];
    }
    bArray<ind_t> tetrahedra(tets.size(), 4u);
    for (size_t i = 0; i < tets.size(); ++i) for (int k = 0; k < 4; ++k) tetrahedra.val(i, k) = static_cast<ind_t>(tets[i][k]);
    LatticeMesh mesh(positions, tetrahedra, limited, n);
    mesh.inputs_ = in;
    // the mesh must fill the zone that ir_moveinto moves points into
    const double volume = mesh.volume(), expected = ir.volume();
    if (std::abs(volume - expected) > 1e-8 * expected)
      throw std::runtime_error("the structured mesh volume " + std::to_string(volume) + " differs from the irreducible zone's "
                               + std::to_string(expected));
    return mesh;
  }

public:
  [[nodiscard]] bool refinement_limited() const { return limited_; }
  [[nodiscard]] long long divisions() const { return divisions_; }
  [[nodiscard]] ind_t number_of_vertices() const { return positions_.size(0); }
  [[nodiscard]] ind_t number_of_tetrahedra() const { return tetrahedra_.size(0); }
  [[nodiscard]] const bArray<double> & get_vertex_positions() const { return positions_; }
  [[nodiscard]] const bArray<ind_t> & get_vertices_per_tetrahedron() const { return tetrahedra_; }
  [[nodiscard]] double volume() const {
    double v{0};
    for (const auto & m: inverse_) v += 1 / (6 * std::abs(determinant({m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11]})));
    return v;
  }
  [[nodiscard]] std::set<size_t> collect_keys() const {
    std::set<size_t> keys;
    const auto n = number_of_vertices();
    for (ind_t i = 0; i < number_of_tetrahedra(); ++i) {
      const auto v = tetrahedra_.view(i).to_std();
      const auto k = permutation_table_keys_from_indicies(v.begin(), v.end(), n);
      keys.insert(k.begin(), k.end());
    }
    return keys;
  }
  [[nodiscard]] std::string to_string() const {
    return "LatticeMesh(" + std::to_string(number_of_vertices()) + " vertices, " + std::to_string(number_of_tetrahedra())
           + " tetrahedra, grid Λ*/" + std::to_string(divisions_) + ")";
  }

  /*! \brief The vertices of the tetrahedron holding x (a single 3-vector), with their weights

  Vertices with zero weight are left out. A point at most `surface_tolerance`
  outside the mesh, in barycentric weight, is located in the nearest tetrahedron;
  a point further outside gives an empty result.
  */
  [[nodiscard]] std::vector<std::pair<ind_t, double>> locate(const bArray<double> & x) const {
    if (x.ndim() != 2u || x.size(0) != 1u || x.size(1) != 3u)
      throw std::runtime_error("locate requires a single 3-element vector.");
    const vec3 p{x.val(0, 0), x.val(0, 1), x.val(0, 2)};
    std::vector<std::pair<ind_t, double>> vw;
    std::array<long long, 3> b{};
    for (int k = 0; k < 3; ++k) b[k] = std::clamp(static_cast<long long>(std::floor((p[k] - lo_[k]) / side_)), 0LL, dims_[k] - 1);
    std::array<double, 4> w{};
    for (size_t j = offsets_[bucket(b)]; j < offsets_[bucket(b) + 1]; ++j) {
      weights(members_[j], p, w);
      if (*std::min_element(w.begin(), w.end()) >= 0) return result(members_[j], w);
    }
    // on a face shared with a tetrahedron in no bucket searched yet, or outside by round-off:
    // take the best tetrahedron in this and the neighbouring buckets
    double best{-surface_tolerance};
    ind_t found{number_of_tetrahedra()};
    std::array<double, 4> found_w{};
    for (long long i = std::max(0LL, b[0] - 1); i <= std::min(dims_[0] - 1, b[0] + 1); ++i)
      for (long long j = std::max(0LL, b[1] - 1); j <= std::min(dims_[1] - 1, b[1] + 1); ++j)
        for (long long k = std::max(0LL, b[2] - 1); k <= std::min(dims_[2] - 1, b[2] + 1); ++k) {
          const auto c = bucket({i, j, k});
          for (size_t m = offsets_[c]; m < offsets_[c + 1]; ++m) {
            weights(members_[m], p, w);
            const double least = *std::min_element(w.begin(), w.end());
            if (least >= best) { best = least; found = members_[m]; found_w = w; }
          }
        }
    if (found == number_of_tetrahedra()) return vw;
    double total{0};
    for (auto & v: found_w) { v = std::max(v, 0.0); total += v; }
    for (auto & v: found_w) v /= total;
    return result(found, found_w);
  }

  template<class HF>
  std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, bool>
  to_hdf(HF & obj, const std::string & entry) const {
    auto group = overwrite_group(obj, entry);
    group.createAttribute("kind", std::string("LatticeMesh"));
    group.createAttribute("refinement_limited", limited_ ? 1 : 0);
    group.createAttribute("divisions", divisions_);
    bool ok{true};
    ok &= positions_.to_hdf(group, "positions");
    ok &= tetrahedra_.to_hdf(group, "vertices_per_tetrahedron");
    if (inputs_) {
      // what refinement needs: the construction's inputs and the refinements so far
      auto r = group.createGroup("refinement");
      r.createDataSet("metric", std::vector<double>(inputs_->metric.begin(), inputs_->metric.end()));
      r.createDataSet("basis", std::vector<double>(inputs_->basis.begin(), inputs_->basis.end()));
      std::vector<long long> ops, cone;
      for (const auto & g: inputs_->ops) ops.insert(ops.end(), g.begin(), g.end());
      for (const auto & c: inputs_->cone) cone.insert(cone.end(), c.begin(), c.end());
      r.createAttribute("operations", inputs_->ops.size());
      r.createAttribute("cone_planes", inputs_->cone.size());
      if (!ops.empty()) r.createDataSet("operations", ops);
      if (!cone.empty()) r.createDataSet("cone", cone);
      std::vector<unsigned long long> marked, counts;
      for (const auto & step: history_marked_) {
        counts.push_back(step.size());
        for (const auto & t: step) marked.insert(marked.end(), t.begin(), t.end());
      }
      r.createAttribute("steps", history_marked_.size());
      if (!history_marked_.empty()) {
        r.createDataSet("marked_counts", counts);
        if (!marked.empty()) r.createDataSet("marked_vertices", marked);
        r.createDataSet("min_edge", history_min_edge_);
      }
    }
    return ok;
  }
  /*! Read a mesh written by `to_hdf`, or the finest layer of a TetGen mesh written by
  an earlier brille (a `TetTri` of layers) */
  template<class HF>
  static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, LatticeMesh>
  from_hdf(HF & obj, const std::string & entry) {
    auto group = obj.getGroup(entry);
    if (group.hasAttribute("layers")) {
      size_t layers;
      group.getAttribute("layers").read(layers);
      auto finest = group.getGroup("layers").getGroup(std::to_string(layers - 1));
      return {bArray<double>::from_hdf(finest, "positions"), bArray<ind_t>::from_hdf(finest, "vertices_per_tetrahedron")};
    }
    int limited{0};
    long long divisions{0};
    group.getAttribute("refinement_limited").read(limited);
    group.getAttribute("divisions").read(divisions);
    LatticeMesh mesh(bArray<double>::from_hdf(group, "positions"), bArray<ind_t>::from_hdf(group, "vertices_per_tetrahedron"), limited != 0, divisions);
    if (group.exist("refinement")) {
      auto r = group.getGroup("refinement");
      Inputs in;
      std::vector<double> metric, basis;
      r.getDataSet("metric").read(metric);
      r.getDataSet("basis").read(basis);
      std::copy(metric.begin(), metric.end(), in.metric.begin());
      std::copy(basis.begin(), basis.end(), in.basis.begin());
      size_t n_ops{0}, n_cone{0}, steps{0};
      r.getAttribute("operations").read(n_ops);
      r.getAttribute("cone_planes").read(n_cone);
      r.getAttribute("steps").read(steps);
      std::vector<long long> ops, cone;
      if (n_ops) r.getDataSet("operations").read(ops);
      if (n_cone) r.getDataSet("cone").read(cone);
      for (size_t i = 0; i < n_ops; ++i) { latticetri::mat3i g{}; std::copy(ops.begin() + 9 * i, ops.begin() + 9 * i + 9, g.begin()); in.ops.push_back(g); }
      for (size_t i = 0; i < n_cone; ++i) in.cone.push_back({cone[3 * i], cone[3 * i + 1], cone[3 * i + 2]});
      mesh.inputs_ = in;
      if (steps) {
        std::vector<unsigned long long> counts, marked;
        r.getDataSet("marked_counts").read(counts);
        if (r.exist("marked_vertices")) r.getDataSet("marked_vertices").read(marked);
        r.getDataSet("min_edge").read(mesh.history_min_edge_);
        size_t at{0};
        for (const auto c: counts) {
          std::vector<std::array<ind_t, 4>> step;
          for (size_t k = 0; k < c; ++k, at += 4)
            step.push_back({static_cast<ind_t>(marked[at]), static_cast<ind_t>(marked[at + 1]), static_cast<ind_t>(marked[at + 2]), static_cast<ind_t>(marked[at + 3])});
          mesh.history_marked_.push_back(std::move(step));
        }
      }
    }
    return mesh;
  }

private:
  static std::array<double, 9> product(const std::array<double, 9> & a, const std::array<double, 9> & b) {
    std::array<double, 9> c{};
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) c[3 * i + j] += a[3 * i + k] * b[3 * k + j];
    return c;
  }
  static vec3 apply(const std::array<double, 9> & m, const vec3 & x) {
    vec3 y{0, 0, 0};
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) y[i] += m[3 * i + k] * x[k];
    return y;
  }
  static double determinant(const std::array<double, 9> & m) {
    return m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
  }
  static std::array<double, 9> inverse(const std::array<double, 9> & m) {
    const double d = determinant(m);
    if (d == 0) throw std::runtime_error("singular matrix");
    return {(m[4] * m[8] - m[5] * m[7]) / d, (m[2] * m[7] - m[1] * m[8]) / d, (m[1] * m[5] - m[2] * m[4]) / d,
            (m[5] * m[6] - m[3] * m[8]) / d, (m[0] * m[8] - m[2] * m[6]) / d, (m[2] * m[3] - m[0] * m[5]) / d,
            (m[3] * m[7] - m[4] * m[6]) / d, (m[1] * m[6] - m[0] * m[7]) / d, (m[0] * m[4] - m[1] * m[3]) / d};
  }
  //! The smallest integer vector along c, which must be (a multiple of) a rational direction
  static latticetri::int3 integer_direction(const vec3 & c) {
    double least{0};
    for (const auto x: c) if (std::abs(x) > 1e-12 * std::max({std::abs(c[0]), std::abs(c[1]), std::abs(c[2])}) && (least == 0 || std::abs(x) < least)) least = std::abs(x);
    for (long long k = 1; k <= 100000; ++k) {
      latticetri::int3 n{};
      bool whole{true};
      for (int i = 0; i < 3; ++i) {
        const double y = static_cast<double>(k) * c[i] / least;
        n[i] = std::llround(y);
        whole &= std::abs(y - static_cast<double>(n[i])) <= 1e-7 * std::max(1.0, std::abs(y));
      }
      if (whole) {
        const long long g = std::gcd(std::gcd(std::llabs(n[0]), std::llabs(n[1])), std::llabs(n[2]));
        for (auto & x: n) x /= g;
        return n;
      }
    }
    throw std::runtime_error("an irreducible wedge normal, (" + std::to_string(c[0]) + ", " + std::to_string(c[1]) + ", "
                             + std::to_string(c[2]) + ") on the primitive reciprocal lattice, is not a lattice direction");
  }
  //! Unit normals, Cartesian, of the irreducible polyhedron's faces that pass through Γ
  static std::vector<vec3> faces_through_gamma(const BrillouinZone & bz) {
    const auto ir = bz.get_ir_polyhedron();
    const auto xyz = ir.vertices().xyz();
    double size{0};
    for (ind_t i = 0; i < xyz.size(0); ++i) for (int k = 0; k < 3; ++k) size = std::max(size, std::abs(xyz.val(i, k)));
    std::vector<vec3> out;
    for (const auto & f: ir.faces().faces()) {
      // Newell's normal and the vertex centre
      vec3 n{0, 0, 0}, p{0, 0, 0};
      for (size_t j = 0; j < f.size(); ++j) {
        const auto a = f[j], b = f[(j + 1) % f.size()];
        for (int k = 0; k < 3; ++k) {
          const int k1 = (k + 1) % 3, k2 = (k + 2) % 3;
          n[k] += (xyz.val(a, k1) - xyz.val(b, k1)) * (xyz.val(a, k2) + xyz.val(b, k2));
          p[k] += xyz.val(a, k) / static_cast<double>(f.size());
        }
      }
      const double length = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
      if (length == 0) continue;
      for (auto & x: n) x /= length;
      if (std::abs(n[0] * p[0] + n[1] * p[1] + n[2] * p[2]) <= 1e-9 * size) out.push_back(n);
    }
    return out;
  }
  //! The face normal parallel (either way) to n, if there is one
  static std::optional<vec3> parallel_to(const vec3 & n, const std::vector<vec3> & faces) {
    const double length = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (const auto & f: faces)
      if (std::abs(f[0] * n[0] + f[1] * n[1] + f[2] * n[2]) >= (1 - 1e-6) * length) return f;
    return std::nullopt;
  }
  [[nodiscard]] size_t bucket(const std::array<long long, 3> & b) const {
    return static_cast<size_t>((b[0] * dims_[1] + b[1]) * dims_[2] + b[2]);
  }
  void weights(const ind_t t, const vec3 & p, std::array<double, 4> & w) const {
    const auto & m = inverse_[t];
    const vec3 d{p[0] - m[0], p[1] - m[1], p[2] - m[2]};
    w[1] = m[3] * d[0] + m[4] * d[1] + m[5] * d[2];
    w[2] = m[6] * d[0] + m[7] * d[1] + m[8] * d[2];
    w[3] = m[9] * d[0] + m[10] * d[1] + m[11] * d[2];
    w[0] = 1 - w[1] - w[2] - w[3];
  }
  [[nodiscard]] std::vector<std::pair<ind_t, double>> result(const ind_t t, const std::array<double, 4> & w) const {
    std::vector<std::pair<ind_t, double>> vw;
    for (int i = 0; i < 4; ++i) if (!brille::approx_float::scalar(w[i], 0.)) vw.emplace_back(tetrahedra_.val(t, i), w[i]);
    return vw;
  }
  void index() {
    const ind_t nt = number_of_tetrahedra();
    inverse_.resize(nt);
    std::vector<std::array<vec3, 2>> boxes(nt);
    double total{0};
    for (ind_t t = 0; t < nt; ++t) {
      std::array<vec3, 4> p;
      for (int k = 0; k < 4; ++k) for (int i = 0; i < 3; ++i) p[k][i] = positions_.val(tetrahedra_.val(t, k), i);
      // columns: the edges from p0; its inverse maps x - p0 to weights 1..3
      const std::array<double, 9> e{p[1][0] - p[0][0], p[2][0] - p[0][0], p[3][0] - p[0][0],
                                    p[1][1] - p[0][1], p[2][1] - p[0][1], p[3][1] - p[0][1],
                                    p[1][2] - p[0][2], p[2][2] - p[0][2], p[3][2] - p[0][2]};
      const auto ei = inverse(e);
      inverse_[t] = {p[0][0], p[0][1], p[0][2], ei[0], ei[1], ei[2], ei[3], ei[4], ei[5], ei[6], ei[7], ei[8]};
      total += std::abs(determinant(e)) / 6;
      boxes[t] = {p[0], p[0]};
      for (int k = 1; k < 4; ++k) for (int i = 0; i < 3; ++i) {
        boxes[t][0][i] = std::min(boxes[t][0][i], p[k][i]);
        boxes[t][1][i] = std::max(boxes[t][1][i], p[k][i]);
      }
    }
    offsets_.assign(2, 0);
    members_.clear();
    if (nt == 0) return;
    lo_ = boxes[0][0];
    vec3 hi = boxes[0][1];
    for (const auto & b: boxes) for (int i = 0; i < 3; ++i) { lo_[i] = std::min(lo_[i], b[0][i]); hi[i] = std::max(hi[i], b[1][i]); }
    // buckets about four mean tetrahedra in volume, and not too many of them
    side_ = std::cbrt(4 * total / static_cast<double>(nt));
    for (;;) {
      for (int i = 0; i < 3; ++i) dims_[i] = std::max(1LL, static_cast<long long>(std::ceil((hi[i] - lo_[i]) / side_)));
      if (dims_[0] * dims_[1] * dims_[2] <= 4 * static_cast<long long>(nt) + 64) break;
      side_ *= 1.25;
    }
    auto range = [&](const std::array<vec3, 2> & box, const int i) {
      const auto a = std::clamp(static_cast<long long>(std::floor((box[0][i] - lo_[i]) / side_)), 0LL, dims_[i] - 1);
      const auto b = std::clamp(static_cast<long long>(std::floor((box[1][i] - lo_[i]) / side_)), 0LL, dims_[i] - 1);
      return std::make_pair(a, b);
    };
    std::vector<size_t> count(static_cast<size_t>(dims_[0] * dims_[1] * dims_[2]) + 1, 0);
    for (const auto & box: boxes) {
      const auto [a0, a1] = range(box, 0);
      const auto [b0, b1] = range(box, 1);
      const auto [c0, c1] = range(box, 2);
      for (auto a = a0; a <= a1; ++a) for (auto b = b0; b <= b1; ++b) for (auto c = c0; c <= c1; ++c) ++count[bucket({a, b, c}) + 1];
    }
    std::partial_sum(count.begin(), count.end(), count.begin());
    offsets_ = count;
    members_.resize(offsets_.back());
    std::vector<size_t> next(offsets_.begin(), offsets_.end() - 1);
    for (ind_t t = 0; t < nt; ++t) {
      const auto [a0, a1] = range(boxes[t], 0);
      const auto [b0, b1] = range(boxes[t], 1);
      const auto [c0, c1] = range(boxes[t], 2);
      for (auto a = a0; a <= a1; ++a) for (auto b = b0; b <= b1; ++b) for (auto c = c0; c <= c1; ++c) members_[next[bucket({a, b, c})]++] = t;
    }
  }
};

}
#endif
