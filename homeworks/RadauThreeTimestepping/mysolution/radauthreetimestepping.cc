/**
 * @file radauthreetimestepping.cc
 * @brief NPDE homework RadauThreeTimestepping
 * @author Erick Schulz, edited by Oliver Rietmann
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include "radauthreetimestepping.h"

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <cmath>
#include <iostream>
#include <unsupported/Eigen/KroneckerProduct>

namespace RadauThreeTimestepping {

/**
 * @brief Implementation of the right hand side (time dependent) source vector
 * for the parabolic heat equation
 * @param dofh A reference to the DOFHandler
 * @param time The time at which to evaluate the source vector
 * @returns The source vector at time `time`
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd rhsVectorheatSource(const lf::assemble::DofHandler &dofh,
                                    double time) {
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  // Right-hand side vector has to be set to zero initially
  Eigen::VectorXd phi(N_dofs);
  //====================
  // Your code goes here
  // conpute the right-hand-side vector arising from the Galerkin finite element discrtizaition of the 
  // linear functional on the right hand side of the spatial variational formulation 
  // the argument dofh supplies information about the mesh and the local-to-global index mapping
  auto f [time] (Eigen::Vector2d x){
    const double pi = 3.1415926; 
    Eigen::Vector2d PI (std::cos(pi*time),std::sin(pi*time)); 
    return (((x-0.5*PI).norm()<0.5)? 1;0); 
  }

  // pointer to current mesh; 
  auto mesh_p= dofh.mesh(); 
  phi.setZero(); 


  //====================
  return phi;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Heat evolution solver: the solver obtains the
 * discrete evolution operator from the Radau3MOLTimestepper class and
 * repeatedly iterates its applicaiton starting from the initial condition
 * @param dofh The DOFHandler object
 * @param m is total number of steps until final time final_time (double)
 * @param final_time The duration for which to solve the PDE
 * @returns The solution at the final timestep
 */
/* SAM_LISTING_BEGIN_6 */
Eigen::VectorXd solveHeatEvolution(const lf::assemble::DofHandler &dofh,
                                   unsigned int m, double final_time) {
  Eigen::VectorXd discrete_heat_sol(dofh.NumDofs());
  //====================
  // Your code goes here
  // solves the IBVP 9.1.1 ober the time interval. using m equidistant timesteps and the finite element 
  // spaces as encoded in the argument dofh; 
  double tau = final_time/m; 
  const ls::uscalfe::size_type N_dof(dofh.NumDofs());//N

  Radau3MOLTimestepper radau_solver(dofh); 
  // starting with the zero initial condition; 
  Eigen::VectorXd descrete_heat_cur = radau_solver.discreteEvolutionOperator(0.0, tau, const Eigen::VectorXd::Zero(N_dofs)); 

  Eigen::VectorXd discrete_heat_next; 
  for (int i =1; i<m; i++){
    discrete_heat_next = radar_solver.discreteEvolutionOperator(i*tau, tau, discrete_heat_cur); 
    discrete_heat_cur = discrete_heat_next; 
  }
  discrete_heat_sol = discrete_heat_cur; 
  //====================
  return discrete_heat_sol;
}
/* SAM_LISTING_END_6 */

/* Implementing member function Eval of class LinFEMassMatrixProvider*/
Eigen::Matrix<double, 3, 3> LinFEMassMatrixProvider::Eval(
    const lf::mesh::Entity &tria) {
  Eigen::Matrix<double, 3, 3> elMat;
  //====================
  // Your code goes here
  //====================
  return elMat;  // return the local mass element matrix
}

/* Implementing constructor of class Radau3MOLTimestepper */
/* SAM_LISTING_BEGIN_4 */
Radau3MOLTimestepper::Radau3MOLTimestepper(const lf::assemble::DofHandler &dofh)
    : dofh_(dofh) {
  //====================
  // Your code goes here
  // Add any additional members you need in the header file
  //====================
}
/* SAM_LISTING_END_4 */

/* Implementation of Radau3MOLTimestepper member functions */
// The function discreteEvolutionOperator() returns the discretized evolution
// operator as obtained from the Runge-Kutta Radau IIA 2-stages method using the
// Butcher table as stored in the Radau3MOLTimestepper class
/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd Radau3MOLTimestepper::discreteEvolutionOperator(
    double time, double tau, const Eigen::VectorXd &mu) const {
  Eigen::VectorXd discrete_evolution_operator(dofh_.NumDofs());
  //====================
  // Your code goes here
  // 
  //firstly we need to compute the increments ki i=1,2 by solving 9.1.15 using 
  //Eigen's direct sparse elimation solver 
  // to that end, you have to build the system matrix as an object of type 
  // Eigen::SparseMatrix<double> by using the Kronecker-product formula 
  // the matrices A and M should have been initialized in the constructor 
  // althernatively, you could already have initialized the two 2N *2N matrices 
  

  //dimension of finite element space: 
  const lf::uscalfe::size_type N_dof(dofh_.NumDofs()); 
  
  // building the linear system for the implicitely defined increments 
  
  // first assembly the right hand side using block initialization; 
  Eigen::VectorXd linearSystem_rhs(2*N_dof);
  // A_, c_ have been initialized in the constructor 
  Eigen::VectorXd rhs_abstraction = A_*mu; 
  linearSystem_rhs << rhsVectorheatSource(dofh_, time+c_(0)*tau)-rhs_abstraction, rhsVectorheatSource(dofh_, time+tau)-rhs_abstraction; 

  // then implicit Runge-Kutta methods leads to systems of equations that must be solved in order to 
  // obtained the increments 

  Eigen::SparseMatrix<double> linSys_mat; 
  linSys_mat = M_Kp_ +tau*A_kp_; 
  Eigen::SparseLU<Eigen::SpaseMatrix<double>> solver; 
  solver.compute(linSys_mat); 
  Eigen::VectorXd k = solver.solve(linearSystem_rhs); 

  discrete_evolution_operator = mu + tau*(k.topRows(N_dof)*b_[0] + b_[1] *k.bottomRows(N_dof)); 



  // 
  //====================
  return discrete_evolution_operator;
}
/* SAM_LISTING_END_5 */

}  // namespace RadauThreeTimestepping
