/**
 * @file burgersequation.cc
 * @brief NPDE homework BurgersEquation code
 * @author Oliver Rietmann
 * @date 15.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "burgersequation.h"

#include <Eigen/Core>
#include <cmath>

namespace BurgersEquation {
/* SAM_LISTING_BEGIN_1 */
constexpr double PI = 3.14159265358979323846;

double Square(double x) { return x * x; }

double f(double x) { return 2.0 / 3.0 * std::sqrt(x * x * x); }

Eigen::VectorXd solveBurgersGodunov(double T, unsigned int N) {
  double h = 5.0 / N;           // meshwidth
  double tau = h;               // timestep = meshwidth by CFL condition
  int m = std::round(T / tau);  // no. of timesteps

  // initialize vector with initial nodal values
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N + 1, -1.0, 4.0);
  Eigen::VectorXd mu =
      x.unaryExpr([](double x) {
         return 0.0 <= x && x <= 1.0 ? Square(std::sin(PI * x)) : 0.0;
       }).eval();

  for (int i=0; i<m; m++){
    for (int int j= N; j>0; j--){
      mu(j) = mu(j) -tau/h*(f(mu(j))-f(mu(j-1))); 
    }
    mu(0)=0.0; 
  }
  //====================
  // Your code goes here

  // implement the function that solves 11.1.2 based on the smi discretization developed in sub-problem 11-1g
  // using N equispaced nodes xj in ipace interval -1, 4. 
  // fpr temporal discretization use explicit Euler timestepping with timestep tau=h, 
  //where h is the mesh width in space. 
  // the function is supposed to return nodal values (averages over dual cells) of the
  // the solution at time T. 
  // for computingthe flux you may assume that the solution is zero outside the domian, uss the initial condition 
  // and find the mu(0) by simply sampling w0 at the nodes of mesh, outside -1,4 the solution is set to zero

  //====================
  return mu;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Converts a large vector on  a grid to a smaller vector correponding to
 * a sub-grid.
 *
 * @param mu vector of function values on a spacial grid of size N_large
 * @param N divides the size N_large of mu
 * @return a vector mu_sub of size N, that represents mu on a sub-grid of size N
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd reduce(const Eigen::VectorXd &mu, unsigned int N) {
  Eigen::VectorXd mu_sub(N + 1);
  int fraction = mu.size() / N;
  for (int j = 0; j < N + 1; ++j) {
    mu_sub(j) = mu(j * fraction);
  }
  return mu_sub;
}

Eigen::Matrix<double, 3, 4> numexpBurgersGodunov() {
  const unsigned int N_large = 3200;
  Eigen::Vector2d T{0.3, 3.0};
  Eigen::Vector4i N{5 * 10, 5 * 20, 5 * 40, 5 * 80};
  Eigen::Vector4d h;
  for (int i = 0; i < 4; ++i) h(i) = 5.0 / N(i);

  Eigen::Matrix<double, 3, 4> result;
  result.row(0) = h.transpose();

  //====================
  // Your code goes here
  //tabulate the discrete error norm as a
  for (int k=0; k<1; k++){
    Eigen::VectorXd mu_ref = solveBurgerGodunov(T(k), 3200); 
    Eigen::Vector4d error; 
    for (int i=0; i<4; i++){
      Eigen::VectorXd mu_ref_sub = reduce(mu_ref, N(i)); 
      Eigen::VectorXd mu = solveBurgerGodunov(T(k), h(i)); 
      error(i) = (mu-mu_ref_sub).lpNorm<1>(); 
ß
    }
    result.row(k+1) = error.transpose();
  }


  //====================

  return result;
}
/* SAM_LISTING_END_2 */

}  // namespace BurgersEquation
