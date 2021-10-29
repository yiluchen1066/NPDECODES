/**
 * @file systemode.h
 * @brief NPDE homework SystemODE
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>

namespace SystemODE {

// Single step of RK4 for the ODE y' = f(y)
/* SAM_LISTING_BEGIN_1 */
template <typename Function>
Eigen::VectorXd rk4step(Function &&f, double h, Eigen::VectorXd &y0) {
  Eigen::VectorXd eval(y0.size());
  //====================
  // Your code goes here
  // performs a single stepsize h of the classical Runge-Kutta method of order 4 
  // for the first-order ODE obtained in problem a (associated right hand side function passed via the functor f
  // ) here y0 contains the vector before the step and the fucntion returns the vector y1 after the time step
  Eigen::VectorXd k1 = f(y0);
  Eigen::VectorXd k2 = f(y0+0.5*h*k1); 
  Eigen::VectorXd k3 = f(y0+0.5*h*k2); 
  Eigen::VectorXd k4 = f(y0+0.5*h*k3); 
  eval = y0 + h*(1/6*k1+2/6*k2+2/6*k3+1/6*k4); 
  //====================
  return eval;
}
/* SAM_LISTING_END_1 */

}  // namespace SystemODE
