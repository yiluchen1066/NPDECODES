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
  for(const lf::mesh::Entity* cell: mesh_p.Entities(0)){
    lf::geometry::Geometry* geo = cell.Geometry(); 
    volume += lf::geomotry:Volume(*geo); 
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
  // flagEntitiesOnBoundary first count the number of adjacentt entities of the next lower
  //codimension. 
  // the criterion for identifying a boundary edge is that is has exactly one adjacent cell

  auto bool_tag = flagEntitiesOnBoundary(mesh_p,1); 
  for(const lf::mesh::Entity* edge: mesh_p.Entities(1)){
    if(bool_tag(edge)==True){
      lf::geometry::Geometry* geo = edge.Geometry(); 
      length += lf::geometry::Volume(*geo); 
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
 
  //====================

  return {volume, length};
}
/* SAM_LISTING_END_3 */

}  // namespace LengthOfBoundary
