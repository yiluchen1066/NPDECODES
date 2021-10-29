/**
 * @file semimprk.cc
 * @brief NPDE homework SemImpRK code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "semimprk.h"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../../lecturecodes/helperfiles/polyfit.h"

namespace SemImpRK {

/* SAM_LISTING_BEGIN_0 */
double CvgRosenbrock() {
  double cvgRate = 0.0;
  // Use polyfit to estimate the rate of convergence
  // for SolveRosenbrock.
  //====================
  // Your code goes here
  // explore the order of method 7.4.2 empirically by applying it to the IVP for the limit cycle 
  // use uniform timesteps for size h=2^(-k), k= 4,...,10, and compute a reference solution yref with 
  // timestep size h =2^(-12). Monitor the maximal error on the temporal mesh 
  // max||yj-yref(tj)||2 and use it to estimate a rate of algebraic convergence by means of linear regression. 

  Eigen::Matrix2d R; 
  R << 0.0, -1.0, 1.0, 0.0; 
  int lambda =1; 
  auto f [R,lambda](Eigen::VectorXd y) -> Eigen::VectorXd{

    return R*y+lambda*(1-y.squaredNorm())*y; 
  }

  auto df [R, lambda](Eigen::VectorXd y) {
    Eigen::Matrix2d J; 
    J << lambda*(1-y.squaredNorm())-2*lambda*y(0)*y(0), -1-2*lambda*y(0)*y(1), 1-2*lambda*y(0)*y(1), lambda*(1-y.squaredNorm())-2*lambda*y(0)**2; 
    return J; 
  }

  int M_ref = std::pow(2,12); 
  std::vector<Eigen::VectorXd> res_ref = SolveRosenbrock(f, df, y0, M_ref, T); 

  Eigen::ArrayXd K= Eigen::ArrayXd::linspace(7,4,10); 

  Eigen::ArrayXd Error(K.size()); 
  Eigen::ArrayXd M(K.size()); 

  for (int i=0; i<K.size(); i++){
    M[i] = std::pow(2,K[i]); 
    std::vector<Eigen::VectorXd> res = SolveRosenbrock(f,df, y0, M[i], T); 
    
    double maxerr =0; 

    for(int j=0; j<res.size(); j++){
      maxerr = std::max(maxerr, (res[j]-res_ref[j*M_ref/M[i]].norm()); 
    }
    Error[i] = maxerr; 
  }


  SolveRosenbrock(Func &&f, Jac &&df, const Eigen::VectorXd &y0, unsigned int M, double T)

  //====================
  return cvgRate;
}
/* SAM_LISTING_END_0 */

}  // namespace SemImpRK
