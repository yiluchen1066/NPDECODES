/**
 * @file
 * @brief NPDE homework "Handling degrees of freedom (DOFs) in LehrFEM++"
 * @author Julien Gacon
 * @date March 1st, 2019
 * @copyright Developed at ETH Zurich
 */

#include "lfppdofhandling.h"

#include <Eigen/Dense>
#include <array>
#include <memory>

#include "lf/assemble/assemble.h"
#include "lf/base/base.h"
#include "lf/geometry/geometry.h"
#include "lf/mesh/mesh.h"
#include "lf/mesh/utils/utils.h"

namespace LFPPDofHandling {

/* SAM_LISTING_BEGIN_1 */
std::array<std::size_t, 3> countEntityDofs(
    const lf::assemble::DofHandler &dofhandler) {
  std::array<std::size_t, 3> entityDofs;
  //====================
  // Your code goes here
  // returns the number of global shape functions (managed by the local-> global index ing mapping encoded
  // in the dofhander) associated with mesh entities of co-domensioon 0, 1, and 2 respectively
  // the co-dimension also serves as index for the returned array
  
  // loop over all mesh entitesm get the inices of the global shape functions
  // associated with them via DofHandker::InteriorGlobalDofIndices()
  // and count them 
  // this number is alo available through DofHandler::NumInteriorDofs()

  // iterate over entities in the mesh and get interior number of dofs for each
  std::shared_ptr<lf::mesh::Mesh> mesh = dofhandler.Mesh(); 
  for (int co_dim =0 ; co_dim <=2; co_dim++){
    entityDofs[co_dim] =0; 
    for (auto el : mesh->Entities(co_dim)){
      entity_Dofs[co_dim] += dofhandler.interiorGlobalDofIndices(el); 
    }
  }
  

  //====================
  return entityDofs;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::size_t countBoundaryDofs(const lf::assemble::DofHandler &dofhandler) {
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  // given an entity, bd\_flags(entity) == true, if the entity is on the
  // boundary
  lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags(
      lf::mesh::utils::flagEntitiesOnBoundary(mesh));
  std::size_t no_dofs_on_bd = 0;
  //====================
  // Your code goes here
  // tells the number of global shape functions associated with mesh entities located on the 
  // boundary 
  // use the function flagEntitiesonBoundary to obtain the array of flags whether an entity is located on the boundary 
  // edges and nodes can be on the boundary; 
  for (auto edge: mesh->Entity(1)){
    if (bd_flags(edge) == TRUE){
      no_dofs_on_bd += dofhandler.InteriorGlobalDofIndices(edge); 
    }
  }

  for (auto node: mesh->Entity(2)){
    if (bd_flags(node) == TRUE){
      no_dofs_on_bd += dofhandler.InteriorGlobalDofIndices(node); 
    }
  }

  //====================
  return no_dofs_on_bd;
}
/* SAM_LISTING_END_2 */

// clang-format off
/* SAM_LISTING_BEGIN_3 */
double integrateLinearFEFunction(
    const lf::assemble::DofHandler& dofhandler,
    const Eigen::VectorXd& mu) {
  double I = 0;
  //====================
  // Your code goes here
  // computes the ...
  // whose nodal bases coefficients are passed in mu 
  // the dofhandler arguments provides the local-> global mapping 
  // check whether the lf::assemble::DofHandler object dofhandler really fits the 
  // lagrangian finite element space 
  // run over all cells
  // get the indices of the tent functions covering them via 
  //DofHandler::GlobalDofIndices, 
  // sum the corresponding entries of the coefficient vectors
  // in linear Lagarangian FE, the integral over the basis functions 
  // over a triangle K is 1/3*vol(k)

  std::shared_ptr<lf::mesh::Mesh> mesh = dofhandler.Mesh(); 


  for(auto cell: mesh->Entities(0)){
    // lf::assemble::DofHandler::GlobalDofIndices 
    // access to indices of global dof's belonging to an entity 
    auto glob_dof_indices = dofhandler.GlobalDofIndices(*cell); 

    double I_dof =0; 
    for(auto dof_indice = glob_dof_indices.start(); dof_indice =glob_dof_indices.end(); dof_indice++){
      I_dof += mu(dof_indice); 
    }
    lf::geometry::Geometry *cell->Geometry(); 
    I_dof *= 1/3*lf::geometry::Volume(cell); 
    I+=I_dof; 
  }
  
  //====================
  return I;
}
/* SAM_LISTING_END_3 */
// clang-format on

/* SAM_LISTING_BEGIN_4 */
double integrateQuadraticFEFunction(const lf::assemble::DofHandler &dofhandler,
                                    const Eigen::VectorXd &mu) {
  double I = 0;
  //====================
  // Your code goes here
  //====================
  return I;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd convertDOFsLinearQuadratic(
    const lf::assemble::DofHandler &dofh_Linear_FE,
    const lf::assemble::DofHandler &dofh_Quadratic_FE,
    const Eigen::VectorXd &mu) {
  if (dofh_Linear_FE.Mesh() != dofh_Quadratic_FE.Mesh()) {
    throw "Underlying meshes must be the same for both DOF handlers!";
  }
  std::shared_ptr<const lf::mesh::Mesh> mesh =
      dofh_Linear_FE.Mesh();                          // get the mesh
  Eigen::VectorXd zeta(dofh_Quadratic_FE.NumDofs());  // initialise empty zeta
  // safety guard: always set zero if you're not sure to set every entry later
  // on for us this shouldn't be a problem, but just to be sure
  zeta.setZero();

  for (const auto *cell : mesh->Entities(0)) {
    // check if the spaces are actually linear and quadratic
    //====================
    // Your code goes here
    //====================
    // get the global dof indices of the linear and quadratic FE spaces, note
    // that the vectors obey the LehrFEM++ numbering, which we will make use of
    // lin\_dofs will have size 3 for the 3 dofs on the nodes and
    // quad\_dofs will have size 6, the first 3 entries being the nodes and
    // the last 3 the edges
    //====================
    // Your code goes here
    // assign the coefficients of mu to the correct entries of zeta, use
    // the previous subproblem 2-9.a
    //====================
  }
  return zeta;
}
/* SAM_LISTING_END_5 */

}  // namespace LFPPDofHandling
