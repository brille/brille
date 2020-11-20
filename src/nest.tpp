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

template<class T, class S>
void Nest<T,S>::construct(const Polyhedron& poly, const size_t max_branchings, const double max_volume){
  SimpleTet root_tet(poly);
  double exponent;
  exponent = std::log(root_tet.maximum_volume()/max_volume)/std::log(static_cast<double>(max_branchings));
  // we want more than one tetrahedra at the root node
  if (root_tet.number_of_tetrahedra()<2 && max_branchings > 0){
    root_tet = SimpleTet(poly, max_volume*std::pow(static_cast<double>(max_branchings)-1, exponent));
    exponent = std::log(root_tet.maximum_volume()/max_volume)/std::log(static_cast<double>(max_branchings));
  }
  verbose_update("Largest root tetrahedron ",root_tet.maximum_volume()," and maximum leaf tetrahedron ",max_volume," give exponent ",exponent);
  // copy-over the vertices of the root node
  vertices_ = root_tet.get_vertices();
  // we need to make a guess about how many vertices we'll need:
  size_t number_density = static_cast<size_t>(poly.get_volume()/max_volume); // shockingly, this is about right for total number of vertices!
  size_t nVerts = vertices_.size(0);
  vertices_.resize(number_density + nVerts);
  // copy over the per-tetrahedron vertex indices to the root's branches
  const auto& tvi{root_tet.get_vertices_per_tetrahedron()};
  for (brille::ind_t i=0; i<tvi.size(0); ++i)
  {
    std::array<brille::ind_t,4> single;
    for (brille::ind_t j=0; j<4u; ++j)
      single[j] = tvi.val(i,j); // no need to adjust indices at this stage
    // create a branch for this tetrahedron
    NestNode branch(single, root_tet.circumsphere_info(i), root_tet.volume(i));
    // and subdivide it if necessary
    if (max_branchings > 0 && branch.volume() > max_volume)
      this->subdivide(branch, 1u, max_branchings, max_volume, exponent, nVerts);
    // storing the resulting branch/leaf at this root
    root_.branches().push_back(branch);
  }
  // ensure that we only keep actual vertices:
  if (vertices_.size(0) > nVerts) vertices_.resize(nVerts);
}

template<class T,class S>
void Nest<T,S>::subdivide(
  NestNode& node, const size_t nBr, const size_t maxBr,
  const double max_volume, const double exp, size_t& nVerts
){
  if (!node.is_leaf()) return; // return node; // we can only branch un-branched nodes
  Polyhedron poly(vertices_.extract(node.boundary().vertices()));
  double mult = (maxBr > nBr) ? std::pow(static_cast<double>(maxBr-nBr),exp) : 1.0;
  SimpleTet node_tet(poly, max_volume*mult);
  // add any new vertices to the object's array, keep a mapping for all:
  std::vector<size_t> map;
  const auto& ntv{node_tet.get_vertices()};
  for (size_t i=0; i<ntv.size(0); ++i){
    const auto vi{ntv.view(i)};
    // If we're clever we can possibly simplify this
    bool notfound{true};
    if (nVerts > 0){
      auto idx = norm(vertices_.view(0,nVerts)-vi).find(brille::cmp::eq,0.);
      if (idx.size()>1)
        throw std::runtime_error("Multiple matching vertices?!");
      if (idx.size()==1){
        map.push_back(idx[0]);
        notfound = false;
      }
    }
    if (notfound) {
      // new vertex. make sure we have room to store it:
      if (vertices_.size(0) < nVerts+1) vertices_.resize(nVerts*2);
      vertices_.set(nVerts, vi);
      map.push_back(nVerts++);
    }
  }
  // copy over the per-tetrahedron vertex indices, applying the mapping as we go
  const auto& tvi{node_tet.get_vertices_per_tetrahedron()};
  for (brille::ind_t i=0; i<tvi.size(0); ++i)
  {
    std::array<brille::ind_t,4> single;
    for (brille::ind_t j=0; j<4u; ++j)
      single[j] = map[tvi.val(i,j)];
    // create a branch for this tetrahedron
    NestNode branch(single, node_tet.circumsphere_info(i), node_tet.volume(i));
    // and subdivide it if necessary
    if (nBr < maxBr && branch.volume() > max_volume)
      this->subdivide(branch, nBr+1u, maxBr, max_volume, exp, nVerts);
    // storing the resulting branch/leaf at this node
    node.branches().push_back(branch);
  }
}
