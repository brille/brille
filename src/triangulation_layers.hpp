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

#ifndef BRILLE_TRIANGULATION_HPP_
#define BRILLE_TRIANGULATION_HPP_
#include <set>
#include <vector>
#include <array>
#include <omp.h>
#include <cassert>
#include <algorithm>
#include "array_latvec.hpp" // defines bArray
#include "tetgen.h"
#include "debug.hpp"
#include "polyhedron.hpp"
#include "approx.hpp"
namespace brille {

template<class T, size_t N> static size_t find_first(const std::array<T,N>& x, const T val){
  auto at = std::find(x.begin(), x.end(), val);
  if (at == x.end()) throw std::logic_error("Value not found?!");
  return std::distance(x.begin(), at);
}

class TetTriLayer{
  using ind_t = brille::ind_t;
  ind_t nVertices;
  ind_t nTetrahedra;
  bArray<double> vertex_positions; // (nVertices, 3)
  bArray<ind_t> vertices_per_tetrahedron; // (nTetrahedra, 4)
  std::vector<std::vector<ind_t>> tetrahedra_per_vertex; // (nVertices,)(1+,)
  std::vector<std::vector<ind_t>> neighbours_per_tetrahedron; // (nTetrahedra,)(1+,)
  bArray<double> circum_centres; // (nTetrahedra, 3);
  std::vector<double> circum_radii; // (nTetrahedra,)
public:
  ind_t number_of_vertices(void) const {return nVertices;}
  ind_t number_of_tetrahedra(void) const {return nTetrahedra;}
  const bArray<double>& get_vertex_positions(void) const {return vertex_positions;}
  const bArray<ind_t>& get_vertices_per_tetrahedron(void) const {return vertices_per_tetrahedron;}
  const bArray<double>& get_circum_centres(void) const {return circum_centres;}
  const std::vector<double>& get_circum_radii(void) const {return circum_radii;}
  Polyhedron get_tetrahedron(const ind_t idx) const {
    if (nTetrahedra <= idx)
      throw std::out_of_range("The requested tetrahedron does not exist.");
    auto tet = vertices_per_tetrahedron.view(idx);
    std::vector<std::vector<int>> vpf{{0,1,2},{0,2,3},{1,0,3},{1,3,2}};
    return Polyhedron(vertex_positions.extract(tet), vpf);
  }

  explicit TetTriLayer()
  : nVertices(0), nTetrahedra(0), vertex_positions(0u,3u),
  vertices_per_tetrahedron(0u,4u), circum_centres(0u,3u)
  {}
  TetTriLayer(const tetgenio& tgio)
  : vertex_positions(0u,3u), vertices_per_tetrahedron(0u,4u), circum_centres(0u,3u)
  {
    nVertices = static_cast<ind_t>(tgio.numberofpoints);
    nTetrahedra = static_cast<ind_t>(tgio.numberoftetrahedra);
    // copy-over all vertex positions:
    vertex_positions.resize(nVertices);
    for (ind_t i=0; i<nVertices; ++i) for (ind_t j=0; j<3u; ++j)
      vertex_positions.val(i,j)=tgio.pointlist[3*i+j];
    // copy-over all tetrahedron vertex indices
    vertices_per_tetrahedron.resize(nTetrahedra);
    for (ind_t i=0; i<nTetrahedra; ++i) for (ind_t j=0; j<4u; ++j)
      vertices_per_tetrahedron.val(i,j) = static_cast<ind_t>(tgio.tetrahedronlist[i*tgio.numberofcorners+j]);
    // Construct the tetrahedra per vertex vector of vectors
    tetrahedra_per_vertex.resize(nVertices);
    for (ind_t i=0; i<nVertices; ++i) for (ind_t j=0; j<nTetrahedra; ++j) for (ind_t k=0; k<4u; ++k)
      if (vertices_per_tetrahedron.val(j,k)==i) tetrahedra_per_vertex[i].push_back(j);

    // Construct the neighbours per tetrahedron vector of vectors
    neighbours_per_tetrahedron.resize(nTetrahedra);
    for (ind_t i=0; i<nTetrahedra; ++i)
    for (ind_t j=0; j<4u; ++j)
    if (tgio.neighborlist[i*4u+j] >= 0)
      neighbours_per_tetrahedron[i].push_back(static_cast<ind_t>(tgio.neighborlist[i*4u+j]));
    // ensure that all tetrahedra have positive (orient3d) volume
    this->correct_tetrahedra_vertex_ordering();
    // Calculate the circumsphere information:
    this->determine_circumspheres();
  }
  // Create a string full of object information:
  std::string to_string(void) const {
    std::string str;
    str  = std::to_string(nVertices) + " vertices";
    str += " in " + std::to_string(nTetrahedra) + " tetrahedra";
    return str;
  }
  // ind_t locate(const bArray<double>& x, std::vector<ind_t>& v, std::vector<double>& w) const {
  //   // no specified tetrahedra to check against, so check them all
  //   std::vector<ind_t> tosearch(nTetrahedra);
  //   std::iota(tosearch.begin(), tosearch.end(), 0u);
  //   return locate(tosearch, x, v, w);
  // }
  // ind_t locate(const std::vector<ind_t>& tosearch, const bArray<double>& x, std::vector<ind_t>& v, std::vector<double>& w) const {
  //   if (x.ndim()!=2u || x.size(0)!=1u || x.size(1)!=3u)
  //     throw std::runtime_error("locate requires a single 3-element vector.");
  //   if (std::any_of(tosearch.begin(), tosearch.end(), [this](ind_t a){return a>=this->nTetrahedra;}))
  //     throw std::domain_error("Out-of-bounds tetrahedra index to search");
  //   ind_t found = this->unsafe_locate(tosearch, x, v, w);
  //   if (found >= nTetrahedra)
  //     throw std::runtime_error("The point was not located!");
  //   return found;
  // }
  // ind_t unsafe_locate(const std::vector<ind_t>& tosearch, const bArray<double>& x, std::vector<ind_t>& v, std::vector<double>& w) const {
  //   std::array<double,4> ws;
  //   v.clear();
  //   w.clear(); // make sure w is back to zero-elements
  //   for (ind_t idx: tosearch){
  //     if (this->unsafe_might_contain(idx, x) && this->unsafe_contains(idx, x, ws)){
  //       // unsafe_contains sets the weights in ws
  //       for (ind_t i=0; i<4u; ++i) if (!brille::approx::scalar(ws[i], 0.)){
  //         v.push_back(vertices_per_tetrahedron.val(idx,i));
  //         w.push_back(ws[i]);
  //       }
  //       return idx;
  //     }
  //   }
  //   return nTetrahedra;
  // }
  ind_t locate(const bArray<double>& x, std::vector<std::pair<ind_t,double>>& vw) const {
    if (x.ndim()!=2u || x.size(0)!=1u || x.size(1)!=3u)
      throw std::runtime_error("locate requires a single 3-element vector.");
    return unsafe_locate(x, vw);
  }
  ind_t unsafe_locate(const bArray<double>& x, std::vector<std::pair<ind_t,double>>& vw) const {
    // no specified tetrahedra to check against, so check them all
    std::vector<ind_t> tosearch(nTetrahedra);
    std::iota(tosearch.begin(), tosearch.end(), 0u);
    return unsafe_locate(tosearch, x, vw);
  }
  ind_t locate(const std::vector<ind_t>& tosearch, const bArray<double>& x, std::vector<std::pair<ind_t,double>>& vw) const {
    if (x.ndim()!=2u || x.size(0)!=1u || x.size(1)!=3u)
      throw std::runtime_error("locate requires a single 3-element vector.");
    if (std::any_of(tosearch.begin(), tosearch.end(), [this](ind_t a){return a>=this->nTetrahedra;}))
      throw std::domain_error("Out-of-bounds tetrahedra index to search");
    ind_t found = this->unsafe_locate(tosearch, x, vw);
    if (found >= nTetrahedra)
      throw std::runtime_error("The point was not located!");
    return found;
  }
  ind_t unsafe_locate(const std::vector<ind_t>& tosearch, const bArray<double>& x, std::vector<std::pair<ind_t,double>>& vw) const {
    std::array<double,4> ws;
    vw.clear();// make sure w is back to zero-elements
    for (ind_t idx: tosearch){
      if (this->unsafe_might_contain(idx, x) && this->unsafe_contains(idx, x, ws)){
        // unsafe_contains sets the weights in ws
        for (ind_t i=0; i<4u; ++i) if (!brille::approx::scalar(ws[i], 0.))
          vw.push_back(std::make_pair(vertices_per_tetrahedron.val(idx,i), ws[i]));
        return idx;
      }
    }
    return nTetrahedra;
  }
  std::vector<ind_t> neighbours(const ind_t vert) const {
    if (vert >= this->nVertices){
      std::string msg = "The provided vertex index is out of bounds";
      throw std::out_of_range(msg);
    }
    std::vector<ind_t> n;
    ind_t v;
    for (ind_t t: this->tetrahedra_per_vertex[vert]) for (ind_t j=0; j<4u; ++j)
    {
      v = this->vertices_per_tetrahedron.val(t,j);
      if ( v!= vert && std::find(n.begin(), n.end(), v) == n.end() ) n.push_back(v);
    }
    return n;
  }
  double volume(const ind_t tet) const {
    const ind_t* i = vertices_per_tetrahedron.ptr(tet,0);
    double v;
    v = orient3d(
      vertex_positions.ptr(i[0],0),
      vertex_positions.ptr(i[1],0),
      vertex_positions.ptr(i[2],0),
      vertex_positions.ptr(i[3],0) )/6.0;
    return v;
  }
  std::array<double,3> volume_statistics() const {
    // total volume, minimum volume, maximum volume
    std::array<double,3> vs{0.,(std::numeric_limits<double>::max)(),std::numeric_limits<double>::lowest()};
    for (ind_t i=0; i<nTetrahedra; ++i){
      double v = this->volume(i);
      vs[0] += v; // add the volume to the total
      if (v<vs[1]) vs[1] = v; // new minimum
      if (v>vs[2]) vs[2] = v; // new maximum
    }
    return vs;
  }
  bool might_contain(const ind_t tet, const bArray<double>& x) const {
    if (x.ndim()!=2u || x.size(0)!=1u || x.size(1)!=3u)
      throw std::runtime_error("x must be a single 3-vector");
    if (tet >= nTetrahedra) return false;
    return this->unsafe_might_contain(tet, x);
  }
  bool contains(const ind_t tet, const bArray<double>& x) const{
    if (x.ndim()!=2u || x.size(0)!=1u || x.size(1)!=3u)
      throw std::runtime_error("x must be a single 3-vector");
    if (tet >= nTetrahedra) return false;
    return this->unsafe_contains(tet, x);
  }
  std::set<size_t> collect_keys() const {
    long long ntets = brille::utils::u2s<long long, ind_t>(this->number_of_tetrahedra());
    ind_t nvert = this->number_of_vertices();
    std::set<size_t> keys;
    #pragma omp parallel for default(none) shared(ntets, keys, nvert)
    for (long long si=0; si<ntets; ++si){
      ind_t i = brille::utils::s2u<ind_t, long long>(si);
      auto v = vertices_per_tetrahedron.view(i).to_std();
      std::set<size_t> t = permutation_table_keys_from_indicies(v.begin(), v.end(), nvert);
      #pragma omp critical
      {
        keys.insert(t.begin(), t.end());
      }
    }
    return keys;
  }
protected:
  bool unsafe_might_contain(const ind_t tet, const bArray<double>& x) const {
    return norm(x-circum_centres.view(tet)).all(brille::cmp::le, circum_radii[tet]);
  }
  bool unsafe_contains(const ind_t tet, const bArray<double>& x) const {
    std::array<double,4> w{0.,0.,0.,0.};
    return this->unsafe_contains(tet,x,w);
  }
  bool unsafe_contains(const ind_t tet, const bArray<double>& x, std::array<double,4>& w) const {
    this->weights(tet, x, w);
    return std::all_of(w.begin(), w.end(), [](double z){ return (z>0.||brille::approx::scalar(z,0.)); });
  }
  void weights(const ind_t tet, const bArray<double>& x, std::array<double,4>& w) const {
    const ind_t* i = vertices_per_tetrahedron.ptr(tet);
    double vol6 = 6.0*this->volume(tet);
    w[0] = orient3d(
      x.ptr(0),
      vertex_positions.ptr(i[1]),
      vertex_positions.ptr(i[2]),
      vertex_positions.ptr(i[3]) )/vol6;
    w[1] = orient3d(
      vertex_positions.ptr(i[0]),
      x.ptr(0),
      vertex_positions.ptr(i[2]),
      vertex_positions.ptr(i[3]) )/vol6;
    w[2] = orient3d(
      vertex_positions.ptr(i[0]),
      vertex_positions.ptr(i[1]),
      x.ptr(0),
      vertex_positions.ptr(i[3]) )/vol6;
    w[3] = orient3d(
      vertex_positions.ptr(i[0]),
      vertex_positions.ptr(i[1]),
      vertex_positions.ptr(i[2]),
      x.ptr(0)                   )/vol6;
  }
  void correct_tetrahedra_vertex_ordering(void){
    for (ind_t i=0; i<nTetrahedra; ++i)
    if (std::signbit(this->volume(i))) // the volume of tetrahedra i is negative
    vertices_per_tetrahedron.swap(i, 0,1); // swap two vertices to switch sign
  }
  void determine_circumspheres(void){
    // ensure that the properties can hold all data
    circum_centres.resize(nTetrahedra);
    circum_radii.resize(nTetrahedra);
    verbose_update("Pull together the circumsphere information for all tetrahedra");
    tetgenmesh tgm; // to get access to circumsphere
    for (ind_t i=0; i<nTetrahedra; ++i){
      const ind_t* v = vertices_per_tetrahedron.ptr(i);
      // use tetgen's circumsphere to find the centre and radius for each tetrahedra
      tgm.circumsphere(
        vertex_positions.ptr(v[0]),
        vertex_positions.ptr(v[1]),
        vertex_positions.ptr(v[2]),
        vertex_positions.ptr(v[3]),
        circum_centres.ptr(i), circum_radii.data()+i);
    }
  }
};

// If we ever get around to making the tetrahedra mesh refinable, then we will
// likely only ever need to refine the lowest layer or a TetTri object, followed
// by re-running find-connections from the next-highest layer.

// If we ever get around to *saving* the mesh to file, we will only need to save
// the *lowest* layer since we can produce a (possibly suboptimal) strata of
// meshes by taking convex-hull of the saved mesh and working down from there.
class TetTri{
  using ind_t = brille::ind_t;
  typedef std::vector<ind_t> TetSet;
  // typedef std::map<size_t, TetSet> TetMap;
  typedef std::vector<TetSet> TetMap;
  //
  std::vector<TetTriLayer> layers;
  std::vector<TetMap> connections;
public:
  explicit TetTri() {};
  TetTri(const std::vector<TetTriLayer>& l)
  : layers(l)
  {
    this->find_connections();
  }
  void find_connections(const size_t highest=0){
    if (highest < layers.size()-1)
    for (size_t i=highest; i<layers.size()-1; ++i)
    connections.push_back(this->connect(i,i+1));
  }
  // make general locate and neighbour methods which drill-down through the layers
  // ind_t locate(const bArray<double>& x, std::vector<ind_t>& v, std::vector<double>& w) const {
  //   if (x.ndim()!=2u || x.size(0)!=1u || x.size(1)!=3u)
  //     throw std::runtime_error("locate requires a single 3-element vector.");
  //   if (layers.size() < 1)
  //     throw std::runtime_error("Can not locate without triangulation");
  //   // find the point within the highest-layer tetrahedra:
  //   ind_t idx = layers[0].locate(x,v,w);
  //   // use the layer-connection map to restrict the search in the next layer's tetrahedra
  //   for (size_t i=1; i<layers.size(); ++i){
  //     const TetSet& tosearch = connections[i-1][idx];
  //     idx = layers[i].unsafe_locate(tosearch, x, v, w);
  //   }
  //   return idx;
  // }
  // ind_t locate(const bArray<double>& x, std::vector<ind_t>& v) const{
  //   std::vector<double> w;
  //   return this->locate(x, v, w);
  // }
  std::vector<std::pair<ind_t,double>>
  locate(const bArray<double>& x) const {
    if (x.ndim()!=2u || x.size(0)!=1u || x.size(1)!=3u)
      throw std::runtime_error("locate requires a single 3-element vector.");
    if (layers.size() < 1)
      throw std::runtime_error("Can not locate without triangulation");
    // setup temporary vertex(index) and weight vectors
    std::vector<std::pair<ind_t,double>> vw;
    // find the point within the highest-layer tetrahedra:
    ind_t idx = layers[0].unsafe_locate(x,vw);
    // use the layer-connection map to restrict the search in the next layer's tetrahedra
    for (size_t i=1; i<layers.size(); ++i){
      const TetSet& tosearch = connections[i-1][idx];
      idx = layers[i].unsafe_locate(tosearch, x, vw);
    }
    return vw;
  }
  // return the neighbouring vertices to a provided mesh-vertex in the lowest layer.
  std::vector<ind_t> neighbours(const bArray<double>& x) const {
    std::vector<std::pair<ind_t,double>> vw = this->locate(x);
    if (vw.size() != 1u){
      std::string msg = "The provided point is not a mesh vertex.";
      throw std::runtime_error(msg);
    }
    return layers.back().neighbours(vw[0].first);
  }
  std::vector<ind_t> neighbours(const ind_t v) const {
    return layers.back().neighbours(v);
  }
  std::string to_string() const {
    std::string str="";
    if (layers.size()>1){
      str += "Layers|";
      for (auto layer: layers) str += layer.to_string() + "|";
    } else {
      str += layers[0].to_string();
    }
    return str;
  }
  // provide convenience functions which pass-through to the lowest layer
  ind_t number_of_tetrahedra() const {return layers.back().number_of_tetrahedra(); }
  ind_t number_of_vertices() const {return layers.back().number_of_vertices(); }
  const bArray<double>& get_vertex_positions() const {return layers.back().get_vertex_positions(); }
  const bArray<ind_t>& get_vertices_per_tetrahedron() const { return layers.back().get_vertices_per_tetrahedron(); }
  std::set<size_t> collect_keys() const {return layers.back().collect_keys();}
private:
  TetMap connect(const size_t high, const size_t low) const{
    omp_set_num_threads(omp_get_max_threads());
    Stopwatch<> stopwatch;
    stopwatch.tic();
    TetMap map(layers[high].number_of_tetrahedra());
    long mapsize = brille::utils::u2s<long, ind_t>(map.size());
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
// this version is necessary with g++ <= 8.3.0
#pragma omp parallel for default(none) shared(map, mapsize) schedule(dynamic)
#else
// this version is necessary with g++ == 9.2.0
#pragma omp parallel for default(none) shared(map, mapsize, high, low) schedule(dynamic)
#endif
    for (long ui=0; ui<mapsize; ++ui){
      ind_t i = brille::utils::s2u<ind_t, long>(ui);
      // initialize the map
      map[i] = TetSet();
      auto cchi = layers[high].get_circum_centres().view(i);
      // get a Polyhedron object for the ith higher-tetrahedra in case we need it
      Polyhedron tethi = layers[high].get_tetrahedron(i);
      std::vector<double> sumrad;
      for (double r: layers[low].get_circum_radii()) sumrad.push_back(layers[high].get_circum_radii()[i]+r);
      // if two circumsphere centers are closer than the sum of their radii
      // they are close enough to possibly overlap:
      auto close_enough = norm(layers[low].get_circum_centres() - cchi).is(brille::cmp::le, sumrad);
      for (ind_t j=0; j < close_enough.size(); ++j) if (close_enough[j])
      {
        bool add = false;
        // check if any vertex of the jth lower-tetrahedra is inside of the ith higher-tetrahedra
        for (ind_t k=0; k<4u; ++k)
        {
          if (!add && layers[high].contains(i, layers[low].get_vertex_positions().view(layers[low].get_vertices_per_tetrahedron().val(j,k))))
            add = true;
        }
        // even if no vertex is inside of the ith higher-tetrahedra, the two tetrahedra
        // can overlap -- and checking for this overlap is complicated.
        // make the Polyhedron class do the heavy lifting.
        // if (add || tethi.intersects(ll.get_tetrahedron(j))) map[i].push_back(j);
        if (add || layers[low].get_tetrahedron(j).intersects(tethi)) map[i].push_back(j);
      }
    }
    stopwatch.toc();
    info_update("Connect ",layers[high].number_of_tetrahedra()," to ",layers[low].number_of_tetrahedra()," completed in ",stopwatch.elapsed()," ms");
    // we now have a TetMap which contains, for every tetrahedral index of the
    // higher level, all tetrahedral indices of the lower level which touch the
    // higher tetrahedron or share some part of its volume.
    return map;
  }
};

template <typename T>
TetTriLayer
triangulate_one_layer(const bArray<T>& verts,
                     const std::vector<std::vector<int>>& vpf,
                     const double max_cell_size=-1.0,
                     const int max_mesh_points=-1)
{
  assert(verts.ndim()==2 && verts.size(1)==3); // otherwise we can't make a 3-D triangulation
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
  tgi.pointlist = nullptr;
  tgi.pointlist = new double[3*tgi.numberofpoints];
  tgi.pointmarkerlist = nullptr;
  tgi.pointmarkerlist = new int[tgi.numberofpoints];
  //tgi.point2tetlist = new int[tgi.numberofpoints];
  using ind_t = brille::ind_t;
  int idx=0;
  for (ind_t i=0; i<verts.size(0); ++i){
    tgi.pointmarkerlist[i] = static_cast<int>(i);
    for (ind_t j=0; j<verts.size(1); ++j) tgi.pointlist[idx++] = verts.val(i,j);
  }
  verbose_update("Initialize and fill the input object's facetlist parameter");
  tgi.numberoffacets = static_cast<int>(vpf.size());
  tgi.facetlist = new tetgenio::facet[tgi.numberoffacets];
  tgi.facetmarkerlist = nullptr;
  tgi.facetmarkerlist = new int[tgi.numberoffacets];
  for (ind_t i=0; i<vpf.size(); ++i){
    tgi.facetmarkerlist[i] = static_cast<int>(i);
    tgi.facetlist[i].numberofpolygons = 1;
    tgi.facetlist[i].polygonlist = new tetgenio::polygon[1];
    tgi.facetlist[i].polygonlist[0].numberofvertices = static_cast<int>(vpf[i].size());
    tgi.facetlist[i].polygonlist[0].vertexlist = nullptr;
    tgi.facetlist[i].polygonlist[0].vertexlist = new int[tgi.facetlist[i].polygonlist[0].numberofvertices];
    for (ind_t j=0; j<vpf[i].size(); ++j)
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
  verbose_update("Constructing TetTriLayer object");
  return TetTriLayer(tgo);
}

template <typename T>
TetTri
triangulate(const bArray<T>& verts,
            const std::vector<std::vector<int>>& vpf,
            const double max_cell_size=-1.0,
            const int layer_count=5,
            const int max_mesh_points=-1)
{
//
  profile_update("Create layered triangulation");
  assert(verts.ndim()==2 && verts.size(1)==3);
  std::vector<TetTriLayer> layers;
  if (layer_count < 2){
    profile_update("Less than two layers requested");
    layers.push_back(triangulate_one_layer(verts, vpf, max_cell_size, max_mesh_points));
  } else {
    profile_update(layer_count," layers requested");
    // the highest-layer is the most-basic tetrahedral triangulation for the input
    layers.push_back(triangulate_one_layer(verts, vpf, -1, max_mesh_points));
    if (max_cell_size > 0.0){
      std::array<double,3> layer_vol_stats = layers[0].volume_statistics();
      profile_update("Highest layer has volume total, min, max: ", layer_vol_stats);
      double highest_maxvol = layer_vol_stats[2];
      TetTriLayer lowest = triangulate_one_layer(verts, vpf, max_cell_size, max_mesh_points);
      layer_vol_stats = lowest.volume_statistics();
      profile_update(" Lowest layer has volume total, min, max: ", layer_vol_stats);
      double lowest_maxvol = layer_vol_stats[2];

      double exponent = std::log(highest_maxvol/lowest_maxvol)/std::log(static_cast<double>(layer_count));
      profile_update("So an exponent of ",exponent," will be used.");
      for (int i=layer_count-1; i>1; --i){ // start from second highest layer, finish at second lowest
        double layer_maxvol = lowest_maxvol * std::pow(static_cast<double>(i), exponent);
        layers.push_back(triangulate_one_layer(verts, vpf, layer_maxvol, max_mesh_points));
      }
      // keep the lowest layer as well
      layers.push_back(lowest);
    }

  }
  return TetTri(layers);
}

} // end namespace brille
#endif // _TRIANGULATION_H_
