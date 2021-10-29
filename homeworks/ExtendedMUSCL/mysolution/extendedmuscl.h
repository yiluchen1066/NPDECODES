#ifndef EXTENDEDMUSCL_H_
#define EXTENDEDMUSCL_H_

/**
 * @file extendedmuscl.h
 * @brief NPDE homework ExtendedMUSCL code
 * @author Oliver Rietmann
 * @date 04.08.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "slopelimfluxdiff.h"

namespace ExtendedMUSCL {

/**
 * @brief Computes the Godunov numerical Flux for
 * flux function f(u) = u(log(u) - 1) and u > 0.
 *
 * @param v left state, v > 0
 * @param w right state, w > 0
 * @return F_GD(v, w).
 */
double logGodunovFlux(double v, double w);

/**
 * @brief Computes the scaled slopes for linear reconstruction,
 * so that it serves as argument for slopelimfluxdiff(per).
 *
 * @param mu_left $\mu_{j-1}$
 * @param mu_center $\mu_{j}$
 * @param mu_right $\mu_{j+1}$
 * @return MC slope limiter as in the problem description
 */
double limiterMC(double mu_left, double mu_center, double mu_right);

/**
 * @brief Solves the ODE $\dot{y} = f(y)$ using the SSP
 * method given in the problem description.
 *
 * @param f right-hand side of ODE modeling std::function<State(State)>
 * @param y inital data, State supports addition and multiplication
 * by a double (from the left). Examples: Eigen::VectorXd or double
 * @param tau small timestep tau > 0
 * @return solution y after time tau.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR, typename State>
State sspEvolop(FUNCTOR &&f, State y, double tau) {
  State y_tau;
  //====================
  // Your code goes here
  // implement the discrete operator \chi according to 11.7.14. 
  // the type State must support basic vector operations like addition and scalar multiplication
  // with a real number
  // the f must supply the right hand-side function f of the autonomous ODE and is expected to have the signature...
  State k; 
  k= y+tau*f(y); 
  k= 3/4*y+1/4*k+1/4*tau*f(k); 
  y_tau = 1/3*y+2/3*k+2/3*tau*f(k); 
  return y_tau;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Solves $\dot{u} + \partial_xf(u) = 0$ for
 * $f(u)=u(\log(u)-1)$ using the finite volume method
 * from the problem description.
 *
 * @param u0 inital data modeling std::function<double(double)>
 * @param T final time, T > 0
 * @param n number of finite volume cells
 * @return approximate solution u(x, T).
 */
/* SAM_LISTING_BEGIN_4 */
template <typename U0_FUNCTOR>
Eigen::VectorXd solveClaw(U0_FUNCTOR &&u0, double T, unsigned int n) {
  // Set up spacial mesh and inital data.
  double a = 0.0;
  double b = 1.0;
  double h = (b - a) / n; // h spatialstep
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n, a + 0.5 * h, b - 0.5 * h);

  // Approximate dual cell averages at t=0
  Eigen::VectorXd mu = x.unaryExpr(u0);

  double alpha = mu.minCoeff();  // lower bound for initial data
  double beta = mu.maxCoeff();   // upper bound for initial data
  assert(alpha > 0.0 && beta > 0.0);
  assert(alpha > 0.0 && beta > 0.0); 

  //====================
  // Your code goes here
  //set timestep tau according to the CFL confition
  double tau = h/(std::max(std::abs(std::log(alpha)),std::abs(std::log(beta)))); 

  auto semi_discrete_rhs = [h_inv] (const Eigen::VectorXd &mu) -> Eigen::VectorXd{
    return -h_inv * slopelimfluxdiffper(mu, &logGodunovFlux, &limiterMC); 
  }
  // timestepping: solve tge semi-discrete ODE
  int N = T/(tau+0.5); // tau timestep
  for(int i =0; i<N; i++){
    mu = sspEvolop(semidiscrete_rhs, mu, tau); 
  }
ß

  //====================

  return mu;
}
/* SAM_LISTING_END_4 */

/* Stores the cell averages of a FV solution obtained by the MUSCL scheme in a
 * 1-periodic setting to file. The filename is passed as an argument. The other
 * arguments are the same as for soveClaw(), which is called by this function.
 */
template <typename U0_FUNCTOR>
void storeMUSCLSolution(const std::string &filename, U0_FUNCTOR &&u0, double T,
                        unsigned int n) {
  const Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols,
                                  ", ", "\n");
  Eigen::VectorXd mu = solveClaw(std::forward<U0_FUNCTOR>(u0), T, n);
  std::ofstream file;
  file.open(filename.c_str());
  file << mu.transpose().format(CSVFormat) << std::endl;
  file.close();
}

/* The function takess two sequences (usually of reals), whose values are read
 * as nodal values on the cell midpoints of an equidistant grid on [0,1].
 * 1-piodic continuation outside [0,1] is assumed. The values on the source grid
 * are linearly interpolated to the cell midpoints of the fine grid.
 */
/* SAM_LISTING_BEGIN_5 */
template <typename VECSOURCE, typename VECDEST>
void interpolate(const VECSOURCE &s, VECDEST &d) {
  // Determine number of cells
  const std::size_t n = s.size();
  const std::size_t N = d.size();
  // Mesh widths/cell sizes
  const double H = 1.0 / n;
  const double h = 1.0 / N;
//====================
// Your code goes here
//loop over nodes of destination mesh 
for(int j=0; j<N; j++){
  // need xd xsl xsr xd
  // need svl svr
  double xd = j*h+0.5; 
  std::size_t k= (std::size_t)std::roundxd/H); 
  double xsl = k*H-0.5; 
  double xsl = xsr + H; 

  double svl=s[(k<1)?n-1:k-1]; 
  double svr=s[(k>=n-1)?0:k]; 

  d[j]= ((xd-xsr)*svl+(xsl-xd)*svr)/H; 

}


//====================
}
/* SAM_LISTING_END_5 */

/* Function conducting a convergence study for the MUSCL scheme. The solution a
 * prescribed final time T is computed for different spatial and temporal
 * resolutions adn then compared with a highly resolved solution obtain on a
 * mesh with 8192 cells. The discrete L1 and L\infty norms of the errors are
 * written to a table u0 passes the initial value, and T the final time.
 */
/* SAM_LISTING_BEGIN_6 */
template <typename U0_FUNCTOR>
void studyCvgMUSCLSolution(U0_FUNCTOR &&u0, double T) {
  // For temporarily storing the number of cells and the associated error norms
  std::vector<std::tuple<std::size_t, double, double>> result{};
  constexpr int n_ref = 8192;
  // Compute reference solution
  std::cout << "Computing reference solution at T = " << T << " with " << n_ref
            << " cells" << std::endl;
  Eigen::VectorXd u_ref{solveClaw(std::forward<U0_FUNCTOR>(u0), T, n_ref)};
  //====================
  // Your code goes here
  // compute solutions on different meshes 
  int min =4; 
  int max =10; 
  for (int i=4; i<=10; i++){
    n = std::pow(2,n); 
    Eigen::VectorXd u ={solveClaw(std::forward<U0_FUNCTOR>(u0), T, n)}; 
    Eigen::VectorXd u_prop = Eigen::VectorXd::Zero(u_ref.size()); 
    interpolate(u,u_prop); 
    double linf_error = (u_prop-u_ref).lpNorm(); 
    double l1_error = (u_prop-u_ref).lpNorm<1>()/n_ref; 

  }
  //====================
  std::cout << "n \t linf error \t l1 error" << std::endl;
  for (auto &data : result) {
    std::cout << std::get<0>(data) << " \t " << std::get<1>(data) << " \t "
              << std::get<2>(data) << std::endl;
  }
}
/* SAM_LISTING_END_6 */

}  // namespace ExtendedMUSCL

#endif  // EXTENDEDMUSCL_H_
