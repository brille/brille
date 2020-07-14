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

#ifndef _TRIANGULATION_SIMPLE_H_
#define _TRIANGULATION_SIMPLE_H_
#include <vector>
#include <array>
#include <cassert>
#include <algorithm>
#include "tetgen.h"
#include "debug.hpp"

class SimpleTet{
  ArrayVector<double> vertex_positions; // (nVertices, 3)
  ArrayVector<size_t> vertices_per_tetrahedron; // (nTetrahedra, 4)
public:
  SimpleTet(void): vertex_positions({3u,0u}), vertices_per_tetrahedron({4u,0u}){}
  SimpleTet(const Polyhedron& poly, const double max_volume=-1, const bool addGamma=false): vertex_positions({3u,0u}), vertices_per_tetrahedron({4u,0u}){
    const ArrayVector<double>& verts{poly.get_vertices()};
    const std::vector<std::vector<int>>& vpf{poly.get_vertices_per_face()};
    // create the tetgenbehavior object which contains all options/switches for tetrahedralize
    verbose_update("Creating `tetgenbehavior` object");
    tetgenbehavior tgb;
    tgb.plc = 1; // we will always tetrahedralize a piecewise linear complex
    if (max_volume > 0){
      tgb.quality = 1; // we will (almost) always improve the tetrahedral mesh
    // tgb.neighout = 1; // we *need* the neighbour information to be stored into tgo.
    // // tgb.mindihedral = 20.; // degrees, avoid very accute edges
      tgb.fixedvolume = 1;
      tgb.maxvolume = max_volume;
    // } else{
    //   tgb.varvolume = 1;
    }
    #ifdef VERBOSE_MESHING
    tgb.verbose = 10000;
    #else
    tgb.quiet = 1;
    #endif
    // make the input and output tetgenio objects and fill the input with our polyhedron
    verbose_update("Creating input and output `tetgenio` objects");
    tetgenio tgi, tgo;
    // we have to handle initializing points/facets, but tetgenio has a destructor
    // which handles deleting all non-NULL fields.
    verbose_update("Initialize and fill the input object's pointlist parameter");
    tgi.numberofpoints = static_cast<int>(verts.size());
    tgi.pointlist = new double[3*tgi.numberofpoints];
    tgi.pointmarkerlist = new int[tgi.numberofpoints];
    size_t numel = verts.numel();
    for (size_t i=0; i<verts.size(); ++i){
      tgi.pointmarkerlist[i] = static_cast<int>(i);
      for (size_t j=0; j<numel; ++j)
        tgi.pointlist[i*numel+j] = verts.getvalue(i,j);
    }
    verbose_update_if(addGamma,"Check whether the Gamma point is present");
    bool gammaPresent{false};
    if (addGamma) {
      for (size_t i=0; i<verts.size(); ++i){
        bool isGamma{true};
	for (size_t j=0; j<numel; ++j)
          isGamma &= approx_scalar(tgi.pointlist[i*numel+j], 0.);
        if (isGamma) for (size_t j=0; j<numel; ++j) tgi.pointlist[i*numel+j] = 0.;
	gammaPresent |= isGamma;
      }
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
      if (addGamma && !gammaPresent){
        tgb.insertaddpoints = 1;
        tetgenio addin;
        addin.numberofpoints = 1;
        addin.pointlist = new double[3];
        addin.pointmarkerlist = new int[1];
        for (int j=0; j<3; ++j) addin.pointlist[j] = 0.;
        addin.pointmarkerlist[0] = verts.size();
        tetrahedralize(&tgb, &tgi, &tgo, &addin);
      } else {
        tetrahedralize(&tgb, &tgi, &tgo);
      }
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
    verbose_update("Copy generated tetgen vertices to SimpleTet object");
    vertex_positions.resize(tgo.numberofpoints);
    for (size_t i=0; i<vertex_positions.size(); ++i) for (size_t j=0; j<3u; ++j)
    vertex_positions.insert(tgo.pointlist[3*i+j], i,j);
    verbose_update("Copy generated tetgen indices to SimpleTet object");
    vertices_per_tetrahedron.resize(tgo.numberoftetrahedra);
    for (size_t i=0; i<vertices_per_tetrahedron.size(); ++i) for (size_t j=0; j<4u; ++j)
    vertices_per_tetrahedron.insert(static_cast<size_t>(tgo.tetrahedronlist[i*tgo.numberofcorners+j]),i,j);
    // ensure that all tetrahedra have positive (orient3d) volume
    this->correct_tetrahedra_vertex_ordering();
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
  double maximum_volume(void) const {
    double vol{0}, maxvol{0};
    for (size_t i=0; i<this->number_of_tetrahedra(); ++i){
      vol = this->volume(i);
      if (vol > maxvol) maxvol = vol;
    }
    return maxvol;
  }
  size_t number_of_vertices(void) const {return vertex_positions.size();}
  size_t number_of_tetrahedra(void) const {return vertices_per_tetrahedron.size();}
  const ArrayVector<double>& get_vertices(void) const {return vertex_positions;}
  const ArrayVector<double>& get_vertex_positions(void) const {return vertex_positions;}
  const ArrayVector<size_t>& get_vertices_per_tetrahedron(void) const {return vertices_per_tetrahedron;}
  std::vector<std::array<size_t,4>> std_vertices_per_tetrahedron(void) const {
    std::vector<std::array<size_t,4>> stdvpt;
    for (size_t i=0; i<this->number_of_tetrahedra(); ++i){
      ArrayVector<size_t> itet = vertices_per_tetrahedron.extract(i);
      stdvpt.push_back({itet.getvalue(0,0), itet.getvalue(0,1), itet.getvalue(0,2), itet.getvalue(0,3)});
    }
    return stdvpt;
  }
  std::array<double,4> circumsphere_info(const size_t tet) const {
    if (tet < vertices_per_tetrahedron.size()){
      std::array<double,4> centre_radius;
      tetgenmesh tgm; // to get access to circumsphere
      tgm.circumsphere(
        vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 0u)),
        vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 1u)),
        vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 2u)),
        vertex_positions.data(vertices_per_tetrahedron.getvalue(tet, 3u)),
        centre_radius.data(),centre_radius.data()+3);
        return centre_radius;
    }
    std::string msg = "The provided tetrahedra index ";
    msg += std::to_string(tet) + " is out of range for ";
    msg += std::to_string(vertices_per_tetrahedron.size());
    msg += " tetrahedra.";
    throw std::out_of_range(msg);
  }
protected:
  void correct_tetrahedra_vertex_ordering(void){
    for (size_t i=0; i<this->number_of_tetrahedra(); ++i)
    if (std::signbit(this->volume(i))) // the volume of tetrahedra i is negative
    vertices_per_tetrahedron.swap(i, 0,1); // swap two vertices to switch sign
  }
};
#endif // _TRIANGULATION_H_
