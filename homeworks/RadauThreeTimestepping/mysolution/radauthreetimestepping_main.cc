/**
 * @file radauthreetimestepping_main.cc
 * @brief NPDE homework RadauThreeTimestepping
 * @author Erick Schulz
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <iostream>
#include <memory>

#include "radauthreetimestepping.h"
#include "radauthreetimesteppingode.h"

using namespace RadauThreeTimestepping;

int main(int /*argc*/, char ** /*argv*/) {
  /* Solving the ODE problem */
  // This function prints to the terminal the convergence rates and average rate
  // of a convergence study performed for the ODE (d/dt)y = -y.
  testConvergenceTwoStageRadauLinScalODE();

  /* Solving the parabolic heat equation */
  // Create a Lehrfem++ square tensor product mesh
  lf::mesh::utils::TPTriagMeshBuilder builder(
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
  // Set mesh parameters following the Builder pattern
  // Domain is the unit square

  builder.setBottomLeftCorner(Eigen::Vector2d{-1.0, -1.0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNumXCells(50)
      .setNumYCells(50);
  auto mesh_p = builder.Build();

  /* SAM_LISTING_BEGIN_1 */
  //====================
  // Your code goes here
  // augment the mian function that it calls the function solveHeatEvolution() from subproblem k
  // with m =50 and outputs the final temperature 

  // finite element space: 
  auto fe_space = std::make_shared<lf::uscalfe::FeSpacelangrangeO1<double>>(mesh_p); 

  // obtain local-> global index mapping for the current finite elemet space
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()}; 
  const lf::uscalfe::size_type N_dof(dofh.NumDofs()); 



  int m = 50; 
  double final_time = 1.0; 

  Eigen::VectorXd discrete_sol = solveHeatEvolution(dofh, m, final_time); 

  
  //====================

  return 0;
}
