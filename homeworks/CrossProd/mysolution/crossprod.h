#ifndef CROSSPROD_H_
#define CROSSPROD_H_

/**
 * @file crossprod.h
 * @brief NPDE homework CrossProd code
 * @author Unknown, Oliver Rietmann
 * @date 31.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/LU>
#include <utility>
#include <vector>

#include <iomanip>
#include <iostream>

#include "implicitrkintegrator.h"

namespace CrossProd {

/* SAM_LISTING_BEGIN_1 */
template <class Function, class Jacobian>
std::vector<Eigen::VectorXd> solve_imp_mid(Function &&f, Jacobian &&Jf,
                                           double T, const Eigen::VectorXd &y0,
                                           unsigned int M) {
  std::vector<Eigen::VectorXd> res(M + 1);
  // Construct the implicit mid-point method with the class
  // implicit_RKIntegrator and execute the .solve() method.
  // Return the vector containing all steps including initial and final value.
  //====================
  // Your code goes here

  int s =1; 
  Eigen::Matrix A(s,s); 
  Eigen::VectorXd b(s); 
  CrossProd::implicitRKIntegrator RK(A,b); 
  res = RK.solve(std::forward<Function>(f), std::forward<Jacobian> (Jf), T, y0, M);

  //====================
  return res;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <class Function, class Jacobian>
std::vector<Eigen::VectorXd> solve_lin_mid(Function &&f, Jacobian &&Jf,
                                           double T, const Eigen::VectorXd &y0,
                                           unsigned int M) {
  std::vector<Eigen::VectorXd> res;
  // Implement the linear implicit mid-point method for
  // an autonomous ODE y' = f(y), y(0) = y0. Return the vector containing
  // all steps including initial and final value.
  //====================
  // Your code goes here

  //initial step size 
  double h= T/M; 
  // store initial data 
  Eigen::VectorXd ytemp1 = y0; 
  Eigen::VectorXd ytemp2 = y1; 
  res.push_back(y0); 

  //pointers to efficient swapping of state vectors 
  Eigen::VectorXd *ynew = ytemp1; 
  Eigen::VectorXd *yold = ytemp2; 
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3,3); 

  for(int i=1; i<M+1; i++){
    *ynew = *yold + h*(I-0.5*h*Jf(yold).lu().solve(f(*yold))); 
    res.push_back(*ynew); 
    std::swap(yold, ynew); 
  }
  //
  //====================

  return res;
}
/* SAM_LISTING_END_2 */

void tab_crossprod();

}  // namespace CrossProd

#endif  // #define CROSSPROD
