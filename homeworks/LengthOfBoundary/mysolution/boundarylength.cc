/**
 * @ file boundarylength.cc
 * @ brief NPDE homework LengthOfBoundary code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "boundarylength.h"

#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>

namespace LengthOfBoundary {

/* SAM_LISTING_BEGIN_1 */
double volumeOfDomain(const std::shared_ptr<lf::mesh::Mesh> mesh_p) {
  double volume = 0.0;
  //====================
  // Your code goes here
  // return the area of the domain covered by a 2D hybrid mesh passed through the 
  // mesh_p pointer argument 
  // you may use the  function lf::geometry::Volume(); 
  // the implementation will employ a loop over all cells 
  // the co-dimension of cell is 0

  double volume; 
  for (lf::mesh::Mesh::Entity *cell ->mesh_p->Entities(0)){
    lf::geometey::Geometry *geo_p = cell ->Geometry(); 
    volume += lf::geometry::Volume(*geo_p);

  }
  //====================

  return volume;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double lengthOfBoundary(const std::shared_ptr<lf::mesh::Mesh> mesh_p) {
  double length = 0.0;
  //====================
  // Your code goes here
  // calculates the length of the boundary of the domain covered 

  // covered by the mesh passed through the mesh_p pointer
  // of course, I can use function flagEntitiesOnBoundary which is to 
  // identify the edges on the boundary 
  auto bd_flags = flagEntitesOnBoundary(mesh_p, 1); 

  // loop over all edges (codimension =1); 
  for (lf::mesh::Mesh::Entity *cell -> Entities(1)){

    if(bg_flags(*cell)){
      if::geometry::Geometry *geo_p = cell -Geomometry; 
      length +=  if::geometry::Volume(geo_p); 
    }
  }

  //====================

  return length;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::pair<double, double> measureDomain(std::string filename) {
  double volume, length;

  //====================
  // Your code goes here
  // readsa 2D hybrid mesh from the .msh file with the name filename and returns both the
  //volume and the length of the boundary of the domain covered by the mesh

  // load mesh into a lehrfem++ object
  // 
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2); 
  lf::io::GmshReader::GmshReader read(mesh_factory, filename); 
  std:shared_ptr<ls::mesh::Mesh> mesh_p = read.mesh(); 

  length = LengthOfBoundary::lengthOfBoundary(mesh_p); 
  volume = lengthOfBoundary::volumeOfDomain(mesh_p); 
  //====================

  return {volume, length};
}
/* SAM_LISTING_END_3 */

}  // namespace LengthOfBoundary
