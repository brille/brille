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

#ifndef _TRIANGULATION_H_
#define _TRIANGULATION_H_
#include <set>
#include <vector>
#include <array>
#include <omp.h>
#include <cassert>
#include <algorithm>
#include "tetgen.h"
#include "debug.hpp"
#include "polyhedron.hpp"

template<class T, size_t N> static size_t find_first(const std::array<T,N>& x, const T val){
  auto at = std::find(x.begin(), x.end(), val);
  if (at == x.end()) throw std::logic_error("Value not found?!");
  return std::distance(x.begin(), at);
}

class TetTriLayer{
  size_t nVertices;
  size_t nTetrahedra;
  ArrayVector<double> vertex_positions; // (nVertices, 3)
  ArrayVector<size_t> vertices_per_tetrahedron; // (nTetrahedra, 4)
  std::vector<std::vector<size_t>> tetrahedra_per_vertex; // (nVertices,)(1+,)
  std::vector<std::vector<size_t>> neighbours_per_tetrahedron; // (nTetrahedra,)(1+,)
  ArrayVector<double> circum_centres; // (nTetrahedra, 3);
  std::vector<double> circum_radii; // (nTetrahedra,)
public:
  size_t number_of_vertices(void) const {return nVertices;}
  size_t number_of_tetrahedra(void) const {return nTetrahedra;}
  const ArrayVector<double>& get_vertex_positions(void) const {return vertex_positions;}
  const ArrayVector<size_t>& get_vertices_per_tetrahedron(void) const {return vertices_per_tetrahedron;}
  const ArrayVector<double>& get_circum_centres(void) const {return circum_centres;}
  const std::vector<double>& get_circum_radii(void) const {return circum_radii;}
  Polyhedron get_tetrahedron(const size_t idx) const {
    if (nTetrahedra <= idx)
      throw std::out_of_range("The requested tetrahedron does not exist.");
    ArrayVector<size_t> tet = vertices_per_tetrahedron.extract(idx);
    std::vector<std::vector<int>> vpf{{0,1,2},{0,2,3},{1,0,3},{1,3,2}};
    return Polyhedron(vertex_positions.extract(tet), vpf);
  }

  TetTriLayer(void): nVertices(0), nTetrahedra(0), vertex_positions({3u,0u}), vertices_per_tetrahedron({4u,0u}), circum_centres({3u,0u}){}
  TetTriLayer(const tetgenio& tgio): vertex_positions({3u,0u}), vertices_per_tetrahedron({4u,0u}), circum_centres({3u,0u}){
    nVertices = static_cast<size_t>(tgio.numberofpoints);
    nTetrahedra = static_cast<size_t>(tgio.numberoftetrahedra);
    // copy-over all vertex positions:
    vertex_positions.resize(nVertices);
    for (size_t i=0; i<nVertices; ++i)
    for (size_t j=0; j<3u; ++j)
    vertex_positions.insert(tgio.pointlist[3*i+j], i,j);
    // copy-over all tetrahedron vertex indices
    vertices_per_tetrahedron.resize(nTetrahedra);
    for (size_t i=0; i<nTetrahedra; ++i)
    for (size_t j=0; j<4u; ++j)
    vertices_per_tetrahedron.insert(static_cast<size_t>(tgio.tetrahedronlist[i*tgio.numberofcorners+j]),i,j);
    for (size_t i=0; i<nTetrahedra; ++i)
    // Construct the tetrahedra per vertex vector of vectors
    tetrahedra_per_vertex.resize(nVertices);
    for (size_t i=0; i<nVertices; ++i)
    for (size_t j=0; j<nTetrahedra; ++j)
    for (size_t k=0; k<4u; ++k)
    if (vertices_per_tetrahedron.getvalue(j,k)==i)
      tetrahedra_per_vertex[i].push_back(j);
    // Construct the neighbours per tetrahedron vector of vectors
    neighbours_per_tetrahedron.resize(nTetrahedra);
    for (size_t i=0; i<nTetrahedra; ++i)
    for (size_t j=0; j<4u; ++j)
    if (tgio.neighborlist[i*4u+j] >= 0)
      neighbours_per_tetrahedron[i].push_back(static_cast<size_t>(tgio.neighborlist[i*4u+j]));
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
  size_t locate(const ArrayVector<double>& x, std::vector<size_t>& v, std::vector<double>& w) const {
    // no specified tetrahedra to check against, so check them all
    std::vector<size_t> tosearch(nTetrahedra);
    std::iota(tosearch.begin(), tosearch.end(), 0u);
    return locate(tosearch, x, v, w);
  }
  size_t locate(const std::vector<size_t>& tosearch, const ArrayVector<double>& x, std::vector<size_t>& v, std::vector<double>& w) const {
    if (x.numel() != 3u || x.size() != 1u)
      throw std::runtime_error("locate requires a single 3-element vector.");
    if (std::any_of(tosearch.begin(), tosearch.end(), [this](size_t a){return a>=this->nTetrahedra;}))
      throw std::domain_error("Out-of-bounds tetrahedra index to search");
    size_t found = this->unsafe_locate(tosearch, x, v, w);
    if (found >= nTetrahedra)
      throw std::runtime_error("The point was not located!");
    return found;
  }
  size_t unsafe_locate(const std::vector<size_t>& tosearch, const ArrayVector<double>& x, std::vector<size_t>& v, std::vector<double>& w) const {
    std::array<double,4> ws;
    v.clear();
    w.clear(); // make sure w is back to zero-elements

    for (size_t idx: tosearch)
    if (this->unsafe_might_contain(idx, x) && this->unsafe_contains(idx, x, ws)){
      // unsafe_contains sets the weights in ws
      for (size_t i=0; i<4u; ++i) if (!approx_scalar(ws[i], 0.)){
        v.push_back(vertices_per_tetrahedron.getvalue(idx, i));
        w.push_back(ws[i]);
      }
      return idx;
    }
    return nTetrahedra;
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
      v = this->vertices_per_tetrahedron.getvalue(t, j);
      if ( v!= vert && std::find(n.begin(), n.end(), v) == n.end() ) n.push_back(v);
    }
    return n;
  }
  double volume(const size_t tet) const {
    double v;
    v = orient3d(
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 0u)),
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 1u)),
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 2u)),
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 3u)) )/6.0;
    return v;
  }
  std::array<double,3> volume_statistics() const {
    // total volume, minimum volume, maximum volume
    std::array<double,3> vs{0.,(std::numeric_limits<double>::max)(),std::numeric_limits<double>::lowest()};
    for (size_t i=0; i<nTetrahedra; ++i){
      double v = this->volume(i);
      vs[0] += v; // add the volume to the total
      if (v<vs[1]) vs[1] = v; // new minimum
      if (v>vs[2]) vs[2] = v; // new maximum
    }
    return vs;
  }
  bool might_contain(const size_t tet, const ArrayVector<double>& x) const {
    if (x.size() < 1u || x.numel()<3u) throw std::runtime_error("x must be a single 3-vector");
    if (tet >= nTetrahedra) return false;
    return this->unsafe_might_contain(tet, x);
  }
  bool contains(const size_t tet, const ArrayVector<double>& x) const{
    if (x.size() < 1u || x.numel()<3u) throw std::runtime_error("x must be a single 3-vector");
    if (tet >= nTetrahedra) return false;
    return this->unsafe_contains(tet, x);
  }
  std::set<size_t> collect_keys() const {
    long long ntets = unsigned_to_signed<long long, size_t>(this->number_of_tetrahedra());
    size_t nvert = this->number_of_vertices();
    std::set<size_t> keys;
    #pragma omp parallel for default(none) shared(ntets, keys, nvert)
    for (long long si=0; si<ntets; ++si){
      size_t i = signed_to_unsigned<size_t, long long>(si);
      auto v = vertices_per_tetrahedron.extract(i).to_std();
      std::set<size_t> t = permutation_table_keys_from_indicies(v.begin(), v.end(), nvert);
      #pragma omp critical
      {
        keys.insert(t.begin(), t.end());
      }
    }
    return keys;
  }
protected:
  bool unsafe_might_contain(const size_t tet, const ArrayVector<double>& x) const {
    return norm(x - circum_centres.extract(tet)).all_approx(Comp::le, circum_radii[tet]);
  }
  bool unsafe_contains(const size_t tet, const ArrayVector<double>& x) const {
    std::array<double,4> w{0.,0.,0.,0.};
    return this->unsafe_contains(tet,x,w);
  }
  bool unsafe_contains(const size_t tet, const ArrayVector<double>& x, std::array<double,4>& w) const {
    this->weights(tet, x, w);
    return std::all_of(w.begin(), w.end(), [](double z){ return (z>0.||approx_scalar(z,0.)); });
  }
  void weights(const size_t tet, const ArrayVector<double>& x, std::array<double,4>& w) const {
    double vol6 = 6.0*this->volume(tet);
    w[0] = orient3d(
      x.data(),
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 1u)),
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 2u)),
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 3u)) )/vol6;
    w[1] = orient3d(
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 0u)),
      x.data(),
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 2u)),
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 3u)) )/vol6;
    w[2] = orient3d(
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 0u)),
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 1u)),
      x.data(),
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 3u)) )/vol6;
    w[3] = orient3d(
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 0u)),
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 1u)),
      vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 2u)),
      x.data()                                                          )/vol6;
  }
  void correct_tetrahedra_vertex_ordering(void){
    for (size_t i=0; i<nTetrahedra; ++i)
    if (std::signbit(this->volume(i))) // the volume of tetrahedra i is negative
    vertices_per_tetrahedron.swap(i, 0,1); // swap two vertices to switch sign
  }
  void determine_circumspheres(void){
    // ensure that the properties can hold all data
    circum_centres.resize(nTetrahedra);
    circum_radii.resize(nTetrahedra);
    verbose_update("Pull together the circumsphere information for all tetrahedra");
    tetgenmesh tgm; // to get access to circumsphere
    for (size_t i=0; i<nTetrahedra; ++i){
      // use tetgen's circumsphere to find the centre and radius for each tetrahedra
      tgm.circumsphere(
        vertex_positions.data(vertices_per_tetrahedron.getvalue(i, 0u)),
        vertex_positions.data(vertices_per_tetrahedron.getvalue(i, 1u)),
        vertex_positions.data(vertices_per_tetrahedron.getvalue(i, 2u)),
        vertex_positions.data(vertices_per_tetrahedron.getvalue(i, 3u)),
        circum_centres.data(i), circum_radii.data()+i);
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
  typedef std::vector<size_t> TetSet;
  // typedef std::map<size_t, TetSet> TetMap;
  typedef std::vector<TetSet> TetMap;
  //
  std::vector<TetTriLayer> layers;
  std::vector<TetMap> connections;
public:
  TetTri() {};
  TetTri(const std::vector<TetTriLayer>& l): layers(l) {
    this->find_connections();
  }
  void find_connections(const size_t highest=0){
    if (highest < layers.size()-1)
    for (size_t i=highest; i<layers.size()-1; ++i)
    connections.push_back(this->connect(i,i+1));
  }
  // make general locate and neighbour methods which drill-down through the layers
  size_t locate(const ArrayVector<double>& x, std::vector<size_t>& v, std::vector<double>& w) const {
    if (layers.size() < 1)
      throw std::runtime_error("Can not locate without triangulation");
    // find the point within the highest-layer tetrahedra:
    size_t idx = layers[0].locate(x,v,w);
    // use the layer-connection map to restrict the search in the next layer's tetrahedra
    for (size_t i=1; i<layers.size(); ++i){
      const TetSet& tosearch = connections[i-1][idx];
      idx = layers[i].locate(tosearch, x, v, w); // switch this to unsafe_locate
    }
    return idx;
  }
  size_t locate(const ArrayVector<double>& x, std::vector<size_t>& v) const{
    std::vector<double> w;
    return this->locate(x, v, w);
  }
  // return the neighbouring vertices to a provided mesh-vertex in the lowest layer.
  std::vector<size_t> neighbours(const ArrayVector<double>& x) const {
    std::vector<size_t> v;
    this->locate(x, v);
    if (v.size() != 1u){
      std::string msg = "The provided point is not a mesh vertex.";
      throw std::runtime_error(msg);
    }
    return layers.back().neighbours(v[0]);
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
  size_t number_of_tetrahedra() const {return layers.back().number_of_tetrahedra(); }
  size_t number_of_vertices() const {return layers.back().number_of_vertices(); }
  const ArrayVector<double>& get_vertex_positions() const {return layers.back().get_vertex_positions(); }
  const ArrayVector<size_t>& get_vertices_per_tetrahedron() const { return layers.back().get_vertices_per_tetrahedron(); }
  std::set<size_t> collect_keys() const {return layers.back().collect_keys();}
private:
  // TetMap connect(const size_t high, const size_t low) const{
  //   Stopwatch<> stopwatch;
  //   stopwatch.tic();
  //   const TetTriLayer& hl{layers[high]}, ll{layers[low]};
  //   const ArrayVector<double>& cch = hl.get_circum_centres();
  //   const ArrayVector<double>& ccl = ll.get_circum_centres();
  //   const std::vector<double>& crh = hl.get_circum_radii();
  //   const std::vector<double>& crl = ll.get_circum_radii();
  //   const ArrayVector<double>& lvrt = ll.get_vertex_positions();
  //   const ArrayVector<size_t>& ltet = ll.get_vertices_per_tetrahedron();
  //   TetMap map(hl.number_of_tetrahedra());
  //   for (size_t i=0; i<map.size(); ++i){
  //     // initialize the map
  //     map[i] = TetSet();
  //     ArrayVector<double> cchi = cch.extract(i);
  //     // get a Polyhedron object for the ith higher-tetrahedra in case we need it
  //     Polyhedron tethi = hl.get_tetrahedron(i);
  //     std::vector<double> sumrad;
  //     for (double r: crl) sumrad.push_back(crh[i]+r);
  //     // if two circumsphere centers are closer than the sum of their radii
  //     // they are close enough to possibly overlap:
  //     for (size_t j: find(norm(ccl - cchi).is_approx(Comp::le, sumrad))){
  //       bool add = false;
  //       // check if any vertex of the jth lower-tetrahedra is inside of the ith higher-tetrahedra
  //       for (size_t k=0; k<4u; ++k) if (!add && hl.contains(i, lvrt.extract(ltet.getvalue(j, k)))) add = true;
  //       // even if no vertex is inside of the ith higher-tetrahedra, the two tetrahedra
  //       // can overlap -- and checking for this overlap is complicated.
  //       // make the Polyhedron class do the heavy lifting.
  //       // if (add || tethi.intersects(ll.get_tetrahedron(j))) map[i].push_back(j);
  //       if (add || ll.get_tetrahedron(j).intersects(tethi)) map[i].push_back(j);
  //       // if (add || ll.get_tetrahedron(j).fuzzy_intersects(tethi)) map[i].push_back(j);
  //     }
  //   }
  //   stopwatch.toc();
  //   info_update("Connect ",hl.number_of_tetrahedra()," to ",ll.number_of_tetrahedra()," completed in ",stopwatch.elapsed()," ms");
  //   // we now have a TetMap which contains, for every tetrahedral index of the
  //   // higher level, all tetrahedral indices of the lower level which touch the
  //   // higher tetrahedron or share some part of its volume.
  //   return map;
  // }
TetMap connect(const size_t high, const size_t low) const{
  omp_set_num_threads(omp_get_max_threads());
  Stopwatch<> stopwatch;
  stopwatch.tic();
  TetMap map(layers[high].number_of_tetrahedra());
  long mapsize = unsigned_to_signed<long, size_t>(map.size());
#if defined(__GNUC__) && !defined(__llvm__) && __GNUC__ < 9
// this version is necessary with g++ <= 8.3.0
#pragma omp parallel for default(none) shared(map, mapsize) schedule(dynamic)
#else
// this version is necessary with g++ == 9.2.0
#pragma omp parallel for default(none) shared(map, mapsize, high, low) schedule(dynamic)
#endif
  for (long i=0; i<mapsize; ++i){
    // initialize the map
    map[i] = TetSet();
    ArrayVector<double> cchi = layers[high].get_circum_centres().extract(i);
    // get a Polyhedron object for the ith higher-tetrahedra in case we need it
    Polyhedron tethi = layers[high].get_tetrahedron(i);
    std::vector<double> sumrad;
    for (double r: layers[low].get_circum_radii()) sumrad.push_back(layers[high].get_circum_radii()[i]+r);
    // if two circumsphere centers are closer than the sum of their radii
    // they are close enough to possibly overlap:
    for (size_t j: find(norm(layers[low].get_circum_centres() - cchi).is_approx(Comp::le, sumrad))){
      bool add = false;
      // check if any vertex of the jth lower-tetrahedra is inside of the ith higher-tetrahedra
      for (size_t k=0; k<4u; ++k) if (!add && layers[high].contains(i, layers[low].get_vertex_positions().extract(layers[low].get_vertices_per_tetrahedron().getvalue(j, k)))) add = true;
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
TetTri triangulate(const ArrayVector<T>& verts,
                   const std::vector<std::vector<int>>& vpf,
                   const double max_cell_size=-1.0,
                   const int layer_count=5,
                   const int max_mesh_points=-1){
//
  info_update("Create layered triangulation");
  assert(verts.numel() == 3);
  std::vector<TetTriLayer> layers;
  if (layer_count < 2){
    info_update("Less than two layers requested");
    layers.push_back(triangulate_one_layer(verts, vpf, max_cell_size, max_mesh_points));
  } else {
    info_update(layer_count," layers requested");
    // the highest-layer is the most-basic tetrahedral triangulation for the input
    layers.push_back(triangulate_one_layer(verts, vpf, -1, max_mesh_points));
    if (max_cell_size > 0.0){
      std::array<double,3> layer_vol_stats = layers[0].volume_statistics();
      info_update("Highest layer has volume total, min, max: ", layer_vol_stats);
      double highest_maxvol = layer_vol_stats[2];
      TetTriLayer lowest = triangulate_one_layer(verts, vpf, max_cell_size, max_mesh_points);
      layer_vol_stats = lowest.volume_statistics();
      info_update(" Lowest layer has volume total, min, max: ", layer_vol_stats);
      double lowest_maxvol = layer_vol_stats[2];

      double exponent = std::log(highest_maxvol/lowest_maxvol)/std::log(static_cast<double>(layer_count));
      info_update("So an exponent of ",exponent," will be used.");
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
template <typename T>
TetTriLayer triangulate_one_layer(const ArrayVector<T>& verts,
                   const std::vector<std::vector<int>>& vpf,
                   const double max_cell_size=-1.0,
                   const int max_mesh_points=-1
) {
  assert(verts.numel() == 3); // otherwise we can't make a 3-D triangulation
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
  tgi.numberofpoints = static_cast<int>(verts.size());
  tgi.pointlist = nullptr;
  tgi.pointlist = new double[3*tgi.numberofpoints];
  tgi.pointmarkerlist = nullptr;
  tgi.pointmarkerlist = new int[tgi.numberofpoints];
  //tgi.point2tetlist = new int[tgi.numberofpoints];
  int idx=0;
  for (size_t i=0; i<verts.size(); ++i){
    tgi.pointmarkerlist[i] = static_cast<int>(i);
    for (size_t j=0; j<verts.numel(); ++j)
      tgi.pointlist[idx++] = verts.getvalue(i,j);
  }
  verbose_update("Initialize and fill the input object's facetlist parameter");
  tgi.numberoffacets = static_cast<int>(vpf.size());
  tgi.facetlist = new tetgenio::facet[tgi.numberoffacets];
  tgi.facetmarkerlist = nullptr;
  tgi.facetmarkerlist = new int[tgi.numberoffacets];
  for (size_t i=0; i<vpf.size(); ++i){
    tgi.facetmarkerlist[i] = static_cast<int>(i);
    tgi.facetlist[i].numberofpolygons = 1;
    tgi.facetlist[i].polygonlist = new tetgenio::polygon[1];
    tgi.facetlist[i].polygonlist[0].numberofvertices = static_cast<int>(vpf[i].size());
    tgi.facetlist[i].polygonlist[0].vertexlist = nullptr;
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
  verbose_update("Constructing TetTriLayer object");
  return TetTriLayer(tgo);
}

#endif // _TRIANGULATION_H_
