/*
// Copyright © 2019,2020 Greg Tucker <greg.tucker@stfc.ac.uk>
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
// along with brille. If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef BRILLE_TRIANGULATION_HPP_
#define BRILLE_TRIANGULATION_HPP_
#include <vector>
#include <array>
#include <cassert>
#include <algorithm>
#include "array_latvec.hpp" // defines bArray
#include "tetgen.h"
#include "debug.hpp"
#include "balltrellis.hpp"
#include "approx.hpp"
namespace brille {

template<class T, size_t N> static size_t find_first(const std::array<T,N>& x, const T val){
  auto at = std::find(x.begin(), x.end(), val);
  if (at == x.end()) throw std::logic_error("Value not found?!");
  return std::distance(x.begin(), at);
}

// // In case we need to store the input polyhedron and mesh-refining parameters:
// class TetgenInput{
//   bArray<double> _vertices;
//   std::vector<std::vector<int>> _vertices_per_facet;
//   double _max_cell_size;
//   double _min_dihedral_angle;
//   double _max_dihedral_angle;
//   double _radius_edge_ratio;
//   int _max_mesh_points;
// public:
//   TetgenInput(
//     const bArray<double>& v, const std::vector<std::vector<int>>& f,
//     const double mcs=-1.0, const double nda=-1.0, const double xda=-1.0,
//     const double rer=-1.0, const int mmp=-1):
//       _vertices(v), _vertices_per_facet(f), _max_cell_size(mcs),
//       _min_dihedral_angle(nda), _max_dihedral_angle(xda),
//       _radius_edge_ratio(rer), _max_mesh_points(mmp) {}
//   const bArray<double>& vertices() const { return _vertices; }
//   const std::vector<std::vector<int>>& vertices_per_facet() const {return _vertices_per_facet; }
//   double max_cell_size() const { return _max_cell_size; }
//   double max_dihedral_angle() const { return _max_dihedral_angle; }
//   double min_dihedral_angle() const { return _min_dihedral_angle; }
//   double radius_edge_ratio() const { return _radius_edge_ratio; }
//   int max_mesh_points() const { return _max_mesh_points; }
// };

// template<class T, template<class> class L, typename=typename std::enable_if<std::is_base_of<bArray<T>,L<T>>::value>::type>
class TetTri{
  // L<T> vertex_positions; // (nVertices, 3), bArray<T>, LQVec<T>, or LDVec<T>
  size_t nVertices;
  size_t nTetrahedra;
  bArray<double> vertex_positions; // (nVertices, 3)
  bArray<size_t> vertices_per_tetrahedron; // (nTetrahedra, 4)
  std::vector<std::vector<size_t>> tetrahedra_per_vertex; // (nVertices,)(1+,)
  std::vector<std::vector<size_t>> neighbours_per_tetrahedron; // (nTetrahedra,)(1+,)
  //tetgenio tgsource; // we need to store the output of tetgen so that we can refine the mesh
  // BallNode tetrahedraTree;
  // std::vector<BallLeaf> leaves;
  Trellis tetrahedraTrellis;
  std::vector<TrellisLeaf> leaves;
public:
  size_t number_of_vertices(void) const {return nVertices;}
  size_t number_of_tetrahedra(void) const {return nTetrahedra;}
  const bArray<double>& get_vertex_positions(void) const {return vertex_positions;}
  const bArray<size_t>& get_vertices_per_tetrahedron(void) const {return vertices_per_tetrahedron;}

  TetTri(void): nVertices(0), nTetrahedra(0), vertex_positions({0u,3u}), vertices_per_tetrahedron({0u,4u}){}
  TetTri(const tetgenio& tgio, const double fraction): vertex_positions({0u,3u}), vertices_per_tetrahedron({0u,4u}){ //, tgsource(tgio){
    nVertices = static_cast<size_t>(tgio.numberofpoints);
    nTetrahedra = static_cast<size_t>(tgio.numberoftetrahedra);
    // copy-over all vertex positions:
    vertex_positions.resize(nVertices);
    for (size_t i=0; i<nVertices; ++i)
    for (size_t j=0; j<3u; ++j)
    vertex_positions.val(i,j) = tgio.pointlist[3*i+j];
    // copy-over all tetrahedron vertex indices
    vertices_per_tetrahedron.resize(nTetrahedra);
    for (size_t i=0; i<nTetrahedra; ++i)
    for (size_t j=0; j<4u; ++j)
    vertices_per_tetrahedron.val(i,j) = static_cast<size_t>(tgio.tetrahedronlist[i*tgio.numberofcorners+j]);
    // for (size_t i=0; i<nTetrahedra; ++i)
    // Construct the tetrahedra per vertex vector of vectors
    tetrahedra_per_vertex.resize(nVertices);
    for (size_t i=0; i<nVertices; ++i)
    for (size_t j=0; j<nTetrahedra; ++j)
    for (size_t k=0; k<4u; ++k)
    if (vertices_per_tetrahedron.val(j,k)==i)
      tetrahedra_per_vertex[i].push_back(j);
    // Construct the neighbours per tetrahedron vector of vectors
    neighbours_per_tetrahedron.resize(nTetrahedra);
    for (size_t i=0; i<nTetrahedra; ++i)
    for (size_t j=0; j<4u; ++j)
    if (tgio.neighborlist[i*4u+j] >= 0)
      neighbours_per_tetrahedron[i].push_back(static_cast<size_t>(tgio.neighborlist[i*4u+j]));
    // ensure that all tetrahedra have positive (orient3d) volume
    this->correct_tetrahedra_vertex_ordering();
    // construct the tree for faster locating:
    this->make_balltree(fraction);
    // Create a string full of object information:
  }
  std::string to_string(void) const {
    std::string str;
    str  = std::to_string(nVertices) + " vertices";
    str += " in " + std::to_string(nTetrahedra) + " tetrahedra";
    str += " with a Trellis[" + tetrahedraTrellis.to_string() + "]";
    return str;
  }
  // emulate the locate functionality of CGAL -- for which we need an enumeration
  enum Locate_Type{ VERTEX, EDGE, FACET, CELL, OUTSIDE_CONVEX_HULL };
  /*! \brief Locate which tetrahedron contains a specified point

  By comparing distances between all tetrehedron vertices and a test point find
  the vertex closest to the point. Then check all tetrahedra containing the
  closest vertex to see whether the test point is within or on the surface of
  the tetrahedron. The second check is performed by constructing four new
  tetrahedra by replacing each of the source tetrahedron's vertices by the test
  point in turn and then comparing their volumes. If all four tetrahedra have
  positive volumes the test point is inside the source tetrahedron. If any one
  has negative volume the test point is outside of the source tetrahedron. If
  one has zero volume the test point lies on one of the faces of the source
  tetrahedron. If two have zero volume the test point lies on one of the source
  tetrahedron edges. And if three have zero volume the test point is one of the
  source tetrahedron's vertices.
  This function indicates which type of relationship the test point has to the
  returned tetrahedron (index) by use of an enumeration. Depending on which
  relationship is present further information is conveyed by up to two
  unsigned integers.

  @param x A single three-element Array that is the test point
  @param [out] type The relationship between the test point and tetrahedron
  @param [out] v0 The first integer conveying relational information
  @param [out] v1 The second integer conveying relational information
  @returns The index of the found tetrahedron

  | type                | returned index                              | v0                        | v1                          |
  |---------------------|---------------------------------------------|---------------------------|-----------------------------|
  | VERTEX              | arbitrary tetrahedron containing the vertex | vertex index              | undefined                   |
  | EDGE                | arbitrary tetrahedron containing the edge   | vertex at one end of edge | vertex at other end of edge |
  | FACET               | one of the two tetrahedra with this face    | vertex opposite face      | undefined                   |
  | CELL                | the tetrahedron containing the point        | undefined                 | undefined                   |
  | OUTSIDE_CONVEX_HULL | zero                                        | undefind                  | undefined                   |

  @note When this method becomes a bottleneck (when the number of vertices
        in the mesh becomes too large) a large performance benefit can be had
        at the expense of preprocessing and memory use for metadata.
        As written now the method requires nVertices distance calculations to
        locate the vertex closest to the supplied point. One way to reduce this
        is to be able to quickly narrow-down the number of vertices which must
        be checked by, e.g., dividing the vertices into 3D bins. Bin dimesnions
        can be chosen to keep a similar number of vertices in each bin.
        The metadata we need to maintain are the bin boundaries, and either a
        linked list of the vertices in each bin or (after sorting the vertices)
        a list of the bin boundaries in the vertices (first vertex in each bin),
        e.g., [0, first vertex in second bin, first vertex in third bin, ...].
        Then, to locate a point's closest vertex, we first find which bin it
        falls into and then only check its distance to vertices in that bin and
        its neighbours, 3³ = 27 bins in total. If we choose an 8×8×8 binning we
        have 512 total bins to keep track of and would only need to check 5.3%
        of the total vertices; for 16×16×16 we would have 4096 bins and would
        only need to check 0.66% of all vertices.
  */
  size_t old_locate(const bArray<double>& x, Locate_Type& type, size_t& v0, size_t& v1) const{
    std::vector<size_t> v;
    std::vector<double> w;
    size_t idx = this->locate(x, v, w);
    switch (v.size()){
      case 0:
        type = OUTSIDE_CONVEX_HULL;
        return 0;
      case 1:
        type = VERTEX;
        v0 = v[0];
        break;
      case 2:
        type = EDGE;
        v0 = v[0];
        v1 = v[1];
        break;
      case 3:
        type = FACET;
        for (size_t i=0; i<4u; ++i){
          v0 = vertices_per_tetrahedron.val(idx, i);
          if (std::find(v.begin(), v.end(), v0) == v.end()) break;
        }
        break;
      case 4:
        type = CELL;
        break;
      default:
        throw std::logic_error("Intentionally unreachable switch statement..");
    }
    return idx;
  }
  // Make a new locator which slots into the interpolation routine more easily
  // size_t locate(const bArray<double>& x, std::vector<size_t>& v, std::vector<double>& w) const {
  //   if (x.numel() != 3u || x.size() != 1u)
  //     throw std::runtime_error("locate requires a single 3-element vector.");
  //   std::array<double,4> ws;
  //   v.clear();
  //   w.clear(); // make sure w is back to zero-elements
  //
  //   size_t tet_idx;
  //   for (tet_idx=0; tet_idx < nTetrahedra; ++tet_idx) if (this->might_contain(tet_idx, x)){
  //     this->weights(tet_idx, x, ws);
  //     // if all weights are greater or equal to ~zero, we can use this tetrahedron
  //     if (std::all_of(ws.begin(), ws.end(), [](double z){return z>0. || brille::approx::scalar(z,0.);})){
  //       for (size_t i=0; i<4u; ++i) if (!brille::approx::scalar(ws[i], 0.)){
  //         v.push_back(vertices_per_tetrahedron.getvalue(tet_idx, i));
  //         w.push_back(ws[i]);
  //       }
  //       break;
  //     }
  //   }
  //   return tet_idx;
  // }
  size_t locate(const bArray<double>& x, std::vector<size_t>& v, std::vector<double>& w) const {
    if (x.ndim()!=2u || x.size(0)!=1u || x.size(1)!=3u)
      throw std::runtime_error("locate requires a single 3-element vector.");
    std::array<double,4> ws;
    v.clear();
    w.clear(); // make sure w is back to zero-elements

    // verbose_update("Find which tetrahedra might contain the point ",x.to_string(0));
    // std::vector<size_t> tets_to_check = tetrahedraTree.all_containing_leaf_indexes(x);
    // verbose_update("The tree would have us search ",tets_to_check);
    // for (size_t tet_idx: tets_to_check) if (this->unsafe_might_contain(tet_idx, x)){

    for (size_t node: tetrahedraTrellis.nodes_to_search(x)) // returns a distance-sorted list of nodes which might have a containing leaf
    for (auto leaf: tetrahedraTrellis.node_leaves(node))
    if (this->unsafe_might_contain(leaf.index(), x)){
      this->weights(leaf.index(), x, ws);
      // if all weights are greater or equal to ~zero, we can use this tetrahedron
      if (std::all_of(ws.begin(), ws.end(), [](double z){return z>0. || brille::approx::scalar(z,0.);})){
        for (size_t i=0; i<4u; ++i) if (!brille::approx::scalar(ws[i], 0.)){
          v.push_back(vertices_per_tetrahedron[leaf.index(), i]);
          w.push_back(ws[i]);
        }
        return leaf.index();
      }
    }
    throw std::runtime_error("The containing tetrahedra was not found?!");
    return nTetrahedra;
  }
  // and a special version which doesn't return the weights
  size_t locate(const bArray<double>& x, std::vector<size_t>& v) const{
    std::vector<double> w;
    return this->locate(x, v, w);
  }

  // /* If the following function is more useful than the preceeding, it could be
  //    advantageous to replicate the above code in this function's for loop.
  //    Either way, this should probably be parallelised with OpenMP.
  // */
  // std::vector<std::vector<size_t>> locate_all_for_interpolation(const bArray<double>& x) const {
  //   if (x.numel()!=3u){
  //     std::string msg = "locate_all requires 3-element vector(s)";
  //     throw std::runtime_error(msg);
  //   }
  //   std::vector<std::vector<size_t>> all_idx(x.size());
  //   for (size_t i=0; i<x.size(); ++i)
  //     all_idx[i] = this->locate_for_interpolation(x.extract(i));
  //   return all_idx;
  // }

  /* Given a vertex in the mesh, return a vector of all of the other vertices to
     which it is connected.
  */
  std::vector<size_t> neighbours(const bArray<double>& x) const {
    std::vector<size_t> v;
    this->locate(x, v);
    if (v.size() != 1u){
      std::string msg = "The provided point is not a mesh vertex.";
      throw std::runtime_error(msg);
    }
    return this->neighbours(v[0]);
  }
  std::vector<size_t> neighbours(const size_t vert) const {
    if (vert >= this->nVertices){
      std::string msg = "The provided vertex index is out of bounds";
      throw std::out_of_range(msg);
    }
    std::vector<size_t> n;
    size_t v;
    for (size_t t: this->tetrahedra_per_vertex[vert])
    // for (size_t v: this->vertices_per_tetrahedron[t]) // would work if vertices_per_tetrahedron was a vector<array<size_t,4>>
    for (size_t j=0; j<4u; ++j){
      v = this->vertices_per_tetrahedron.val(t, j);
      if ( v!= vert && std::find(n.begin(), n.end(), v) == n.end() ) n.push_back(v);
    }
    return n;
  }
  double volume(const size_t tet) const {
    double v;
    const ind_t* i = vertices_per_tetrahedron.ptr(tet);
    v = orient3d(
      vertex_positions.ptr(i[0]), vertex_positions.ptr(i[1]),
      vertex_positions.ptr(i[2]), vertex_positions.ptr(i[3]))/6.0;return v;
  }
  bool might_contain(const size_t tet, const bArray<double>& x) const {
    if (x.ndim()!=2u || x.size(0)!=1u || x.size(1)!=3u)
      throw std::runtime_error("x must be a single 3-vector");
    if (tet >= nTetrahedra) return false;
    return this->unsafe_might_contain(tet, x);
  }
protected:
  bool unsafe_might_contain(const size_t tet, const bArray<double>& x) const {
    return leaves[tet].fuzzy_contains(x);
  }
  void weights(const size_t tet, const bArray<double>& x, std::array<double,4>& w) const {
    double vol6 = 6.0*this->volume(tet);
    const ind_t* i = vertices_per_tetrahedron.ptr(tet);
    const auto& p{vertex_positions};
    w[0] = orient3d(x.ptr(0),    p.ptr(i[1]), p.ptr(i[2]), p.ptr(i[3]) )/vol6;
    w[1] = orient3d(p.ptr(i[0]), x.ptr(0),    p.ptr(i[2]), p.ptr(i[3]) )/vol6;
    w[2] = orient3d(p.ptr(i[0]), p.ptr(i[1]), x.ptr(0),    p.ptr(i[3]) )/vol6;
    w[3] = orient3d(p.ptr(i[0]), p.ptr(i[1]), p.ptr(i[2]), x.ptr(0)    )/vol6;
  }
  void correct_tetrahedra_vertex_ordering(void){
    for (size_t i=0; i<nTetrahedra; ++i)
    if (std::signbit(this->volume(i))) // the volume of tetrahedra i is negative
    vertices_per_tetrahedron.swap(i, 0,1); // swap two vertices to switch sign
  }
  void make_balltree(const double fraction){
    // Construct a vector of BallLeaf objects for each tetrahedra:
    // std::vector<BallLeaf> leaves;
    leaves.clear();
    std::array<double,3> centre;
    double radius;
    verbose_update("Pull together the circumsphere information for all tetrahedra");
    tetgenmesh tgm; // to get access to circumsphere
    const auto& p{vertex_positions};
    for (size_t i=0; i<nTetrahedra; ++i){
      const ind_t* v = vertices_pertetrahedron.ptr(i);
      // use tetgen's circumsphere to find the centre and radius for each tetrahedra
      tgm.circumsphere(p.ptr(v[0]), p.ptr(v[1]), p.ptr(v[2]), p.ptr(v[3]), centre.data(), &radius);
      // leaves.push_back(BallLeaf(centre, radius, i));
      leaves.push_back(TrellisLeaf(centre, radius, i));
    }
    // // construct the full tree structure at once using an algorithm similar to
    // // Omohundro's Kd algorithm in 'Five Balltree Construction Algorithms'
    // if (nTetrahedra > 0){
    //   // we will construct a binary tree. If we have N leaves to place at the
    //   // end of the branches, they can occupy as many as N branches.
    //   // A tree with N final branches has log₂(N) levels, so we will never need
    //   // to exceed this.
    //   // If we reduce the number of branchings we increase the number of leaves
    //   // per terminal branch. There must be some trade-off between tree-granularity
    //   // and tree-size for how quickly a containing tetrahedra can be located.
    //   size_t maximum_branchings = static_cast<size_t>(std::log2(nTetrahedra))-2;
    //   maximum_branchings = 3;
    //   verbose_update("Construct a tree for tetrahedra location with up to ",maximum_branchings," branchings.");
    //   tetrahedraTree = construct_balltree(leaves, maximum_branchings);
    //   verbose_update("Tree for locating tetrahedra now exists");
    //   // info_update("Tetrahedra tree:\n",tetrahedraTree.to_string());
    // }
    tetrahedraTrellis = construct_trellis(leaves, fraction);
  }
};

template <class T>
TetTri
triangulate(
  const bArray<T>& verts,
  const std::vector<std::vector<int>>& vpf,
  const double max_cell_size=-1.0,
  const double min_dihedral=-1.0,
  const double max_dihedral=-1.0,
  const double radius_edge_ratio=-1.0,
  const int max_mesh_points=-1,
  const double fraction=1.0
) {
  assert(verts.ndim()==2 && verts.size(1)==3);// otherwise we can't make a 3-D triangulation
  // create the tetgenbehavior object which contains all options/switches for tetrahedralize
  verbose_update("Creating `tetgenbehavior` object");
  tetgenbehavior tgb;
  tgb.plc = 1; // we will always tetrahedralize a piecewise linear complex
  tgb.quality = 1; // we will (almost) always improve the tetrahedral mesh
  tgb.neighout = 1; // we *need* the neighbour information to be stored into tgo.
  if (max_cell_size > 0){
    // volume constraint with a specified size
    tgb.fixedvolume = 1;
    tgb.maxvolume = max_cell_size;
  } else{
    //volume constraint without a specified size?
    tgb.varvolume = 1;
  }
  if (max_mesh_points>0) tgb.steinerleft = max_mesh_points;
  if (radius_edge_ratio>0) tgb.minratio = radius_edge_ratio;
  if (min_dihedral>0) tgb.mindihedral = min_dihedral;
  if (max_dihedral>0) tgb.optmaxdihedral = max_dihedral;
  #ifndef VERBOSE_MESHING
  tgb.quiet = 1;
  #endif
  #ifdef VERBOSE_MESHING
  tgb.verbose = 10000;
  #endif

  // make the input and output tetgenio objects and fill the input with our polyhedron
  verbose_update("Creating input and output `tetgenio` objects");
  tetgenio tgi, tgo;
  // we have to handle initializing points/facets, but tetgenio has a destructor
  // which handles deleting all non-NULL fields.
  verbose_update("Initialize and fill the input object's pointlist parameter");
  tgi.numberofpoints = static_cast<int>(verts.size(0));
  tgi.pointlist = new double[3*tgi.numberofpoints];
  tgi.pointmarkerlist = new int[tgi.numberofpoints];
  //tgi.point2tetlist = new int[tgi.numberofpoints];
  int idx=0;
  for (size_t i=0; i<verts.size(0); ++i){
    tgi.pointmarkerlist[i] = static_cast<int>(i);
    for (size_t j=0; j<verts.size(1); ++j)
      tgi.pointlist[idx++] = verts.val(i,j);
  }
  verbose_update("Initialize and fill the input object's facetlist parameter");
  tgi.numberoffacets = static_cast<int>(vpf.size());
  tgi.facetlist = new tetgenio::facet[tgi.numberoffacets];
  tgi.facetmarkerlist = new int[tgi.numberoffacets];
  for (size_t i=0; i<vpf.size(); ++i){
    tgi.facetmarkerlist[i] = static_cast<int>(i);
    tgi.facetlist[i].numberofpolygons = 1;
    tgi.facetlist[i].polygonlist = new tetgenio::polygon[1];
    tgi.facetlist[i].polygonlist[0].numberofvertices = static_cast<int>(vpf[i].size());
    tgi.facetlist[i].polygonlist[0].vertexlist = new int[tgi.facetlist[i].polygonlist[0].numberofvertices];
    for (size_t j=0; j<vpf[i].size(); ++j)
      tgi.facetlist[i].polygonlist[0].vertexlist[j] = vpf[i][j];
  }
  // The input is now filled with the piecewise linear complex information.
  // so we can call tetrahedralize:
  verbose_update("Calling tetgen::tetrahedralize");
  try {
      tetrahedralize(&tgb, &tgi, &tgo);
  } catch (const std::logic_error& e) {
    std::string msg = "tetgen threw a logic error with message\n" + std::string(e.what());
    throw std::runtime_error(msg);
  } catch (const std::runtime_error& e) {
    std::string msg = "tetgen threw a runtime_error with message\n" + std::string(e.what());
    throw std::runtime_error(msg);
  } catch (...) {
    std::string msg = "tetgen threw an undetermined error";
    throw std::runtime_error(msg);
  }
  verbose_update("Constructing TetTri object");
  return TetTri(tgo, fraction);
}

} // end namespace brille
#endif // _TRIANGULATION_H_
