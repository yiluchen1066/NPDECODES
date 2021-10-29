/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 06.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "solve.h"

#include <Eigen/Core>
#include <iostream>

#include "mylinearfeelementmatrix.h"
#include "mylinearloadvector.h"

namespace ElementMatrixComputation {

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solvePoissonBVP() {
  // Convert the globally defined function f to a LehrFEM++ mesh function object
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};

  // The basis expansion coefficient vector for the finite-element solution
  Eigen::VectorXd solution = Eigen::VectorXd::Zero(1);

  //====================
  // Your code goes here
  // that is to call solve() in order to compute the basis expansion coefficients of a solution 
  // of the Neumann boundary value problem
  // make use of the LehrFEM++'s local PROVIDER types for element matrices for the bilinear form 
  //for the weak-laplacian
  // and element vectors for the linear form 
  // initialize the lf::uscalfe::LinearFELocalLoadVector appropriately so that the right-hand side source function 
  // provided by the funtion f(Eigen::Vector2d x); 

  // define the element matrix and element vector builders and solve the system 
  lf::uscalfe::LinearFELaplaceElementMatrix elmat_builder; 
  lf::uscalfe::LinearFELocalLoadVector<double, decitype(mf_f)> elvec_builder(mf_f); 
  Eigen::VectorXd solve(elmat_builder, elvec_provider)； 
  //====================

  return solution;
}
/* SAM_LISTING_END_2 */

Eigen::VectorXd solveNeumannEq() {
  // Define the solution vector
  Eigen::VectorXd solution;

  //====================
  // Your code goes here
  //====================

  return solution;
}

}  // namespace ElementMatrixComputation
