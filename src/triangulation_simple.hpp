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

#ifndef BRILLE_TRIANGULATION_SIMPLE_HPP_
#define BRILLE_TRIANGULATION_SIMPLE_HPP_
/*! \file
    \author Greg Tucker
    \brief A class to interact with TetGen in the simplest case
*/
#include "array_.hpp" // defines bArray
#include "tetgen.h"
namespace brille {

/*! \brief A class to handle interaction with `TetGen`

`TetGen` is much more flexible than `brille` requires. This class exists to
drive triangulation of a single convex polyhedron with limited control over the
triangulation parameters.
Following the `TetGen` `tetrahedralize` call, the resulting vertices and
tetrahedra indexing is collected before freeing the TetGet structure.
*/
template<class T, template<class> class A>
class SimpleTet{
  using vert_t = A<T>;
  using poly_t = polyhedron::Poly<T,A>;
  vert_t vertex_positions;        /*!< The triangulated vertex positions */ // (nVertices, 3)
  bArray<ind_t> vertices_per_tetrahedron; /*!< The vertices in each of the triangulated tetrahedra */ // (nTetrahedra, 4)
public:
  //! Explicit empty constructor
  explicit SimpleTet()
  : vertex_positions(0u,3u), vertices_per_tetrahedron(0u,4u)
  {}
  /*! \brief Triangulate a convex polyhedron

  \param poly       the Polyhedron to triangulate
  \param max_volume The maximum tetrahedron volume allowed within the
                    triangulation, in the same units as the polyhedron volume
  \param addGamma   If true the point (0,0,0) will be ensured to exist in the
                    triangulated vertrices
  */
  explicit SimpleTet(const poly_t& poly, const double max_volume=-1, const bool addGamma=false)
  : vertex_positions(0u,3u), vertices_per_tetrahedron(0u,4u)
  {
    const auto& verts{poly.vertices()};
    const auto& vpf{poly.faces().faces()};
    debug_update("Create SimpleTet from vertices, and facets\n",verts.to_string(),",\n", vpf);
    // create the tetgenbehavior object which contains all options/switches for tetrahedralize
    debug_update("Creating `tetgenbehavior` object");
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
    debug_update("Creating input and output `tetgenio` objects");
    tetgenio tgi, tgo;
    // we have to handle initializing points/facets, but tetgenio has a destructor
    // which handles deleting all non-NULL fields.
    debug_update("Initialize and fill the input object's pointlist parameter");
    tgi.numberofpoints = static_cast<int>(verts.size(0));
    tgi.pointlist = new double[3*tgi.numberofpoints];
    tgi.pointmarkerlist = new int[tgi.numberofpoints];
    ind_t numel = verts.size(1);
    for (ind_t i=0; i<verts.size(0); ++i){
      tgi.pointmarkerlist[i] = static_cast<int>(i);
      for (ind_t j=0; j<numel; ++j)
        tgi.pointlist[i*numel+j] = verts.val(i,j);
    }
    debug_update_if(addGamma,"Check whether the Gamma point is present");
    bool gammaPresent{false};
    if (addGamma) {
      for (ind_t i=0; i<verts.size(0); ++i){
        bool isGamma{true};
        for (ind_t j=0; j<numel; ++j)
          isGamma &= brille::approx_float::scalar(tgi.pointlist[i*numel+j], 0.);
        if (isGamma) for (ind_t j=0; j<numel; ++j) tgi.pointlist[i*numel+j] = 0.;
	      gammaPresent |= isGamma;
      }
    }
    debug_update("Initialize and fill the input object's facetlist parameter");
    tgi.numberoffacets = static_cast<int>(vpf.size());
    tgi.facetlist = new tetgenio::facet[tgi.numberoffacets];
    tgi.facetmarkerlist = new int[tgi.numberoffacets];
    for (size_t i=0; i<vpf.size(); ++i){
      tgi.facetmarkerlist[i] = static_cast<int>(i);
      tgi.facetlist[i].numberofpolygons = 1;
      tgi.facetlist[i].polygonlist = new tetgenio::polygon[1];
      tgi.facetlist[i].numberofholes = 0;
      tgi.facetlist[i].holelist = nullptr;
      tgi.facetlist[i].polygonlist[0].numberofvertices = static_cast<int>(vpf[i].size());
      tgi.facetlist[i].polygonlist[0].vertexlist = new int[tgi.facetlist[i].polygonlist[0].numberofvertices];
      for (size_t j=0; j<vpf[i].size(); ++j)
        tgi.facetlist[i].polygonlist[0].vertexlist[j] = static_cast<int>(vpf[i][j]);
    }
    // The input is now filled with the piecewise linear complex information.
    // so we can call tetrahedralize:
    debug_update("Calling tetgen::tetrahedralize");
    try {
      if (addGamma && !gammaPresent){
        tgb.insertaddpoints = 1;
        tetgenio addin;
        addin.numberofpoints = 1;
        addin.pointlist = new double[3];
        addin.pointmarkerlist = new int[1];
        for (int j=0; j<3; ++j) addin.pointlist[j] = 0.;
        addin.pointmarkerlist[0] = verts.size(0);
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
    debug_update("Copy generated tetgen vertices to SimpleTet object");
    vertex_positions.resize(tgo.numberofpoints);
    for (ind_t i=0; i<vertex_positions.size(0); ++i)
      for (ind_t j=0; j<3u; ++j)
        vertex_positions.val(i,j) = tgo.pointlist[3*i+j];
    debug_update("Copy generated tetgen indices to SimpleTet object");
    vertices_per_tetrahedron.resize(tgo.numberoftetrahedra);
    for (ind_t i=0; i<vertices_per_tetrahedron.size(0); ++i)
      for (ind_t j=0; j<4u; ++j)
        vertices_per_tetrahedron.val(i,j) = static_cast<ind_t>(tgo.tetrahedronlist[i*tgo.numberofcorners+j]);
    // ensure that all tetrahedra have positive (orient3d) volume
    this->correct_tetrahedra_vertex_ordering();
  }
  //! Return the volume of the indexed tetrahedron
  [[nodiscard]] double volume(const ind_t tet) const {
    double v;
    const ind_t* i = vertices_per_tetrahedron.ptr(tet);
    v = orient3d(
      vertex_positions.ptr(i[0]),
      vertex_positions.ptr(i[1]),
      vertex_positions.ptr(i[2]),
      vertex_positions.ptr(i[3]) )/6.0;
    return v;
  }
  //! Return the volume of the largest tetrahedron
  [[nodiscard]] double maximum_volume() const {
    double vol{0}, maxvol{0};
    for (ind_t i=0; i<this->number_of_tetrahedra(); ++i){
      vol = this->volume(i);
      if (vol > maxvol) maxvol = vol;
    }
    return maxvol;
  }
  //! Return the number of triangulated vertices
  [[nodiscard]] ind_t number_of_vertices() const {return vertex_positions.size(0);}
  //! Return the number of triangulated tetrahedra
  [[nodiscard]] ind_t number_of_tetrahedra() const {return vertices_per_tetrahedron.size(0);}
  //! Return a constant reference to the triangulated vertices
  [[nodiscard]] const vert_t& get_vertices() const {return vertex_positions;}
  //! Return a constant reference to the triangulated vertices
  [[nodiscard]] const vert_t& get_vertex_positions() const {return vertex_positions;}
  //! Return a constant reference to the triangulated tetrahedra vertex indices
  [[nodiscard]] const bArray<ind_t>& get_vertices_per_tetrahedron() const {return vertices_per_tetrahedron;}
  //! Convert the tetrahedra vertex indices to a nested standard container
  [[nodiscard]] std::vector<std::array<ind_t,4>> std_vertices_per_tetrahedron() const {
    std::vector<std::array<ind_t,4>> stdvpt;
    for (ind_t i=0; i<this->number_of_tetrahedra(); ++i){
      const ind_t* v = vertices_per_tetrahedron.ptr(i);
      stdvpt.push_back({{v[0], v[1], v[2], v[3]}});
    }
    return stdvpt;
  }
  /*! \brief Use TetGen to determine circumsphere information

  For any tetrahedron there exists a sphere which contains on its surface the
  four tetrahedron vertices. The circumscribing sphere is characterised by its
  centre and its radius, which can be used to quickly decide if a fifth point
  in space has any chance of being *inside* the tetrahedron.

  \param tet the index of the tetrahedron within all triangulated tetrahedra
  \returns the three elements of the circumsphere centre and its radius
  */
  [[nodiscard]] std::array<double,4> circumsphere_info(const ind_t tet) const {
    if (tet < vertices_per_tetrahedron.size(0)){
      const ind_t* i = vertices_per_tetrahedron.ptr(tet);
      std::array<double,4> centre_radius{0,0,0,0};
      tetgenmesh tgm; // to get access to circumsphere
      tgm.circumsphere(
        vertex_positions.ptr(i[0]), vertex_positions.ptr(i[1]),
        vertex_positions.ptr(i[2]), vertex_positions.ptr(i[3]),
        centre_radius.data(),centre_radius.data()+3);
        return centre_radius;
    }
    std::string msg = "The provided tetrahedra index ";
    msg += std::to_string(tet) + " is out of range for ";
    msg += std::to_string(vertices_per_tetrahedron.size(0));
    msg += " tetrahedra.";
    throw std::out_of_range(msg);
  }
protected:
  void correct_tetrahedra_vertex_ordering(){
    for (ind_t i=0; i<this->number_of_tetrahedra(); ++i)
    if (std::signbit(this->volume(i))) // the volume of tetrahedra i is negative
    vertices_per_tetrahedron.swap(i, 0,1); // swap two vertices to switch sign
  }
};

} // end namespace brille
#endif // BRILLE_TRIANGULATION_SIMPLE_HPP_
