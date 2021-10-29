/**
 * @file extendedmuscl.cc
 * @brief NPDE homework ExtendedMUSCL code
 * @author Oliver Rietmann
 * @date 04.08.2019
 * @copyright Developed at ETH Zurich
 */

#include "extendedmuscl.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>

namespace ExtendedMUSCL {

/* SAM_LISTING_BEGIN_1 */
double logGodunovFlux(double v, double w) {
  double godunov_numerical_flux;
  auto f = [](double u) { return u * (std::log(u) - 1.0); };
  auto df = [](double u) { return std::log(u); };
  //====================
  // Your code goes here

  if (v >= w)
  {
    godunov_numerical_flux=std::max(f(v), f(w));
  }else
  {
    if (v<w && w<=1)
    {
      godunov_numerical_flux =f(w); 
    } else if (v<=1 && w>1){
      godunov_numerical_flux =f(1.0); 
    } else if (v>=1 && w>v){
      godunov_numerical_flux = f(v); 
    }
  }
  //====================
  return godunov_numerical_flux;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_4 */
double limiterMC(double mu_left, double mu_center, double mu_right) {
  double scaled_slope;

  //====================
  // Your code goes here
  double slope_c = (mu_right-mu_left)/2; 
  double slope_l = (mu_center-mu_left)*2; 
  double slope_r = (mu_right-mu_center)*2; 

  if (slope_c >0 && slope_l >0 && slope_r >0){
    scaled_slope = std::min(slope_c, std::min(slope_l,slope_r)); 
  } else if (slope_c <0 && slope_c <0 && slope_r <0){
    scaled_slope = std::max(slope_c, std::max(slope_l,slope_r)); 
  } else {
    scaled_slope = 0; 
  }


  //====================

  return scaled_slope;
}
/* SAM_LISTING_END_4 */

}  // namespace ExtendedMUSCL
