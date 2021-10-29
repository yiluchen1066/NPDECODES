#ifndef SEMIMPRK_H_
#define SEMIMPRK_H_

/**
 * @file semimprk.h
 * @brief NPDE homework SemImpRK code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>
#include <vector>

namespace SemImpRK {

// Solve the autonomous IVP y' = f(y), y(0) = y0 using the Rosenbrock method
/* SAM_LISTING_BEGIN_0 */
template <class Func, class Jac>
std::vector<Eigen::VectorXd> SolveRosenbrock(Func &&f, Jac &&df,
                                             const Eigen::VectorXd &y0,
                                             unsigned int M, double T) {
  // Will contain all states computed by the ROW-SSM
  std::vector<Eigen::VectorXd> res(M + 1);
  //====================
  // Your code goes here
  // that aapplies the Rosenbrock method 7.4.2 for solving an initial value problem for an 
  res[0]=y0; 
  double h = T/M; 
  int size = y0.size(); // N 
  double a = 1./(std::sqrt(2)+2.); 

  Eigen::VectorXd k1(size); 
  Eigen::VectorXd k2(size); 

  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(size, size);
  Eigen::MatrixXd J(size,size); 
  Eigen::MatrixXd W(size,size); 

  for (int iter =1; iter < M+1; iter++){
    Eigen::VectorXd &yold = res[i-1]; 
    J=df(yold); 
    W=I-a*h*J; 
    k1 = W.lu().solve(f(yold));
    k2 = W.lu().solve(f(yold+0.5*h*k1)-a*h*J*k1); 
    res[i] = yold + h*k2; 
  }

  //====================
  return res;
}
/* SAM_LISTING_END_0 */

double CvgRosenbrock();

}  // namespace SemImpRK

#endif  // #ifndef SEMIMPRK_H_
