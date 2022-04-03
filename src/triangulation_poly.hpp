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

#ifndef BRILLE_TRIANGULATION_POLY_HPP_
#define BRILLE_TRIANGULATION_POLY_HPP_
/*! \file
    \author Greg Tucker
    \brief A class to interact with TetGen in the simplest case
*/

#include "polyhedron_flex.hpp"
#include "array_tetgen.hpp"

namespace brille::polyhedron {

template<class T, template<class> class A>
std::enable_if_t<isBareArray<T,A>, std::tuple<int, A<T>, Array2<ind_t>>>
triangulate(const T max_volume, const bool addGamma, const A<T>& points, const Faces& faces){
  const auto& f{faces.faces()};
  debug_update("Triangulate vertices and facets\n",points.to_string(),",\n", f);
  debug_update("Creating `tetgenbehavior` object with options/switches for tetrahedralize");
  tetgenbehavior tgb;
  tgb.plc = 1; // we will always tetrahedralize a piecewise linear complex
  if (max_volume > 0){
    tgb.quality = 1; // we will (almost) always improve the tetrahedral mesh
    tgb.fixedvolume = 1;
    tgb.maxvolume = max_volume;
  }
  // if max_volume < 0, force the number of Steiner points to be zero?
  // this would make triangulating un-triangulatable regions break, probably.
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
  tgi.numberofpoints = static_cast<int>(points.size(0));
  tgi.pointlist = new double[3*tgi.numberofpoints];
  tgi.pointmarkerlist = new int[tgi.numberofpoints];
  for (ind_t i=0; i<points.size(0); ++i){
    tgi.pointmarkerlist[i] = static_cast<int>(i);
    for (ind_t j=0; j<3u; ++j)
      tgi.pointlist[i*3u+j] = points.val(i,j);
  }
  debug_update_if(addGamma,"Check whether the Gamma point is present");
  bool gammaPresent{false};
  if (addGamma){
    gammaPresent = points.row_is(brille::cmp::eq, 0 * points.view(0)).any();
  }
  debug_update("Initialize and fill the input object's facetlist parameter");
  tgi.numberoffacets = static_cast<int>(f.size());
  tgi.facetlist = new tetgenio::facet[tgi.numberoffacets];
  tgi.facetmarkerlist = new int[tgi.numberoffacets];
  for (size_t i=0; i<f.size(); ++i){
    tgi.facetmarkerlist[i] = static_cast<int>(i);
    tgi.facetlist[i].numberofpolygons = 1;
    tgi.facetlist[i].polygonlist = new tetgenio::polygon[1];
    tgi.facetlist[i].numberofholes = 0;
    tgi.facetlist[i].holelist = nullptr;
    tgi.facetlist[i].polygonlist[0].numberofvertices = static_cast<int>(f[i].size());
    tgi.facetlist[i].polygonlist[0].vertexlist = new int[tgi.facetlist[i].polygonlist[0].numberofvertices];
    for (size_t j=0; j<f[i].size(); ++j)
      tgi.facetlist[i].polygonlist[0].vertexlist[j] = static_cast<int>(f[i][j]);
  }
  // The input is now filled with the piecewise linear complex information.
  // so we can call tetrahedralize:
  debug_update("Calling tetgen::tetrahedralize");
  try {
    if (addGamma && !gammaPresent){
      tgb.insertaddpoints = 1;
      tetgenio add_in;
      add_in.numberofpoints = 1;
      add_in.pointlist = new double[3];
      add_in.pointmarkerlist = new int[1];
      for (int j=0; j<3; ++j) add_in.pointlist[j] = 0.;
      add_in.pointmarkerlist[0] = points.size(0);
      tetrahedralize(&tgb, &tgi, &tgo, &add_in);
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
  Array2<T> position(tgo.numberofpoints, 3u);
  for (ind_t i=0; i<position.size(0); ++i)
    for (ind_t j=0; j<3u; ++j)
      position.val(i,j) = tgo.pointlist[3*i+j];
  debug_update("Copy generated tetgen indices to SimpleTet object");
  Array2<ind_t> tetrahedra(tgo.numberoftetrahedra, 4u);
  for (ind_t i=0; i<tetrahedra.size(0); ++i) {
    for (ind_t j = 0; j < 4u; ++j) {
      // tetgen returns vertex order different than I expect, so we can swap any two indices
      // or permute their order by one
      tetrahedra.val(i, (j + 1) % 4) = static_cast<ind_t>(tgo.tetrahedronlist[i * tgo.numberofcorners + j]);
    }
  }
  auto extra = tgo.numberofpoints - tgi.numberofpoints;
  if (addGamma && !gammaPresent) extra -= 1; // we know about Gamma already
  return std::make_tuple(extra, position, tetrahedra);
}
template<class T, template<class> class A>
std::enable_if_t<isLatVec<T,A>, std::tuple<int, A<T>, Array2<ind_t>>>
triangulate(const T max_volume, const bool addGamma, const A<T>& points, const Faces& faces){
  auto [n, v, t] = triangulate(max_volume, addGamma, points.xyz(), faces);
  return std::make_tuple(n, from_xyz_like(points, v), t);
}

/*! \brief A class to handle interaction with `TetGen`

`TetGen` is much more flexible than `brille` requires. This class exists to
drive triangulation of a single convex polyhedron with limited control over the
triangulation parameters.
Following the `TetGen` `tetrahedralize` call, the resulting vertices and
tetrahedra indexing is collected before freeing the TetGet structure.
*/
template<class T, template<class> class A>
class LQPolyTet{
  int extra_;
  A<T> points_;   /*!< The triangulated vertex positions */ // (nVertices, 3)
  bArray<ind_t> tetrahedra_; /*!< The vertices in each of the triangulated tetrahedra */ // (nTetrahedra, 4)
public:
  explicit LQPolyTet() = default;
  /*! \brief Triangulate a convex polyhedron

  \param poly       the Polyhedron to triangulate
  \param max_volume The maximum tetrahedron volume allowed within the
                    triangulation, in the same units as the polyhedron volume
  \param addGamma   If true the point (0,0,0) will be ensured to exist in the
                    triangulated vertrices
  */
  explicit LQPolyTet(const Poly<T,A>& poly, const bool addGamma=false, const T max_volume=-1)
  {
    std::tie(extra_, points_, tetrahedra_) = triangulate(max_volume, addGamma, poly.vertices(), poly.faces());
    // ensure that all tetrahedra have positive (orient3d) volume
    this->correct_tetrahedra_vertex_ordering();
  }
  [[nodiscard]] bool any() const {return extra_ > 0;}
  [[nodiscard]] int count() const {return extra_;}
  //! Return the volume of the indexed tetrahedron
  T volume(const ind_t tet) const {
    return triple_product(points_, tetrahedra_.view(tet)).sum() / T(6);
  }
  //! Return the volume of the largest tetrahedron
  T maximum_volume() const {
    return triple_product(points_, tetrahedra_).max();
  }
  //! Return the number of triangulated vertices
  [[nodiscard]] ind_t number_of_vertices() const {return points_.size(0);}
  //! Return the number of triangulated tetrahedra
  [[nodiscard]] ind_t number_of_tetrahedra() const {return tetrahedra_.size(0);}
  //! Return a constant reference to the triangulated vertices
  const A<T>& get_vertices() const {return points_;}
  //! Return a constant reference to the triangulated tetrahedra vertex indices
  [[nodiscard]] const bArray<ind_t>& get_vertices_per_tetrahedron() const {return tetrahedra_;}
  //! Convert the tetrahedra vertex indices to a nested standard container
  [[nodiscard]] std::vector<std::array<ind_t,4>> std_vertices_per_tetrahedron() const {
    std::vector<std::array<ind_t,4>> stdvpt;
    for (ind_t i=0; i<this->number_of_tetrahedra(); ++i){
      const ind_t* v = tetrahedra_.ptr(i);
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
  std::array<T,4> circumsphere_info(const ind_t tet) const {
    if (tet < tetrahedra_.size(0)){
      return four_point_sphere(points_, tetrahedra_.view(tet));
    }
    std::string msg = "The provided tetrahedra index ";
    msg += std::to_string(tet) + " is out of range for ";
    msg += std::to_string(tetrahedra_.size(0));
    msg += " tetrahedra.";
    throw std::out_of_range(msg);
  }
protected:
  void correct_tetrahedra_vertex_ordering(){
    for (ind_t i=0; i<this->number_of_tetrahedra(); ++i) {
      if (std::signbit(this->volume(i))) {
        tetrahedra_.swap(i, 0, 1);
      }
    }
  }
};

} // end namespace brille
#endif // BRILLE_TRIANGULATION_POLY_HPP_
