/**
 * @file systemode_main.cc
 * @brief NPDE homework SystemODE
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "../../../lecturecodes/helperfiles/polyfit.h"
#include "systemode.h"

/* SAM_LISTING_BEGIN_0 */
int main() {
  // PARAMETERS
  double T = 1;
  int n = 5;

  // INITIAL VALUE
  Eigen::VectorXd y0(2 * n);
  for (int i = 0; i < n; ++i) {
    y0(i) = (i + 1.) / n;
    y0(i + n) = -1;
  }

  // SETUP
  double conv_rate = 0;
  std::cout << std::setw(8) << "M" << std::setw(20) << "Error" << std::endl;

  //====================
  // Your code goes here
  // apply the classical Runge-Kutta method of order 4 to solve a particular initial value problem 
  // for the ODE derivide in problem a. 


  // build the tridiagonal C matrix 
  Eigen::SparseMatrix<double> C(n,n); 
  C.reserve(Eigen::VectorXi::Contant(n,3)); 
  C.insert(0,0)=2; 
  for (int i=1; i<n; i++){
    C.insert(i,i)=2; 
    C.insert(i-1,i)=-1;
    C.insert(i,i-1)=-1; 
  }
  C.compress();


  // compute the right-hand side f 
  auto f [n,C] (Eigen::VectorXi y) {

    Eigen::VectorXd fy(2*n); 
    fy.head(n) = y.tail(n); 
    //compute r
    Eigen::VectorXd r(n); 
    r.insert(0) = y(0) *(y(0)+y(1)); 
    r.insert(n-1) = y(n-1)*(y(n-1)+y(n-1)); 
    for (int i=1; i<n-1; i++){
      r.insert(i) = y(i)*(y(i-1)+y(i+1)); 
    }
    Eigen::SparseLU <Eigen::SparseMatrix <double>>csolver; 
    csolver.compute(C); 
    fy.tail(n) = csolver.solve(C); 
    return fy; 
  }

  // compute the exact solution 
  //use N = 2^12 steps to calculate an approximate exact solution 
  int stepNum = std::pow(2,12); 
  double h = T/stepNum; 
  Eigen::VectorXd y_init; 
  Eigen::VectorXd y_exact; 
  y_init=y0; 

  for (int i= 1; i<stepNum; i++){
    y_exact = SystemODE::rk4step(f, h, y_init); 
    y_init = y_exact; 
  }

  // calculate solution using N = 2, 2^2, 2^3, ....., 2^10 timesteps 
  kmax = 10; 
  Eigen::VectorXd Error(kmax); 
  for (int k=1;k< kmax+1; k++){
    int stepNum2 = std::pow(2,k); 
    doubel h_2 = T/stepNum2; 
    Eigen::VectorXd Y_init; 
    Eigen::VectorXd Y_exact; 
    Y_init=y0; 

    for (int i=1; i<stepNum2; i++){
      Y_exact = SystemODE::rk4step(f,h_2, Y_init); 
      Y_init = Y_exact; 
    }

    Error(k) = (Y_exact-y_exact).norm();
  }

  // 
  //====================

  std::cout << "Convergence rate: " << std::round(std::abs(conv_rate))
            << std::endl;

  return 0;
}
/* SAM_LISTING_END_0 */
