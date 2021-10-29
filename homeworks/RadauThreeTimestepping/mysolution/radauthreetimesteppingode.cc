/**
 * @file radauthreetimesteppingode.cc
 * @brief NPDE homework RadauThreeTimestepping
 * @author Erick Schulz
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include "radauthreetimesteppingode.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace RadauThreeTimestepping {

/* SAM_LISTING_BEGIN_1 */
std::vector<double> twoStageRadauTimesteppingLinScalODE(unsigned int m) {
  std::vector<double> sol_vec;
  //====================
  // Your code goes here
  //that uses m equidistant steps of the 2-stage Radau method to solve the linear scalar initial
  // value problem 9.1.10
  // on the time interval [0,5] and returns the generated sequence 
  sol_vec.push_back(1.0); 
  double tau = 5/m; // m steps; tau stepsizes; 

  double evolution_operator = (1-tau*(1+tau/6)/((1+5/12*tau)*(1+tau/4)+tau**2/16)); 

  for (int i = 1; i<m+1; i++){
    sol_vec.push_back(evolution_operator*sol_vec.at(i));  
  }

  //====================
  return sol_vec;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void testConvergenceTwoStageRadauLinScalODE() {
  constexpr int nIter = 10;       // total number of iterations
  double max_norm_errors[nIter];  // errors vector for all approx. sols
  double rates[nIter - 1];        // The rates of convergence
  double avg_rate = 0.0;  // The average rate of convergence over all iterations

  //====================
  // Your code goes here
  // produces an output that permits tou qualitatively and quantitively accesss the convergence of the 
  // 2-stage radau implicit runge-kutta single step methods for the initial value problem 
  // one way to proceed is to produce a sequence of N solutions with m =10*2**k eduidistant timesteps
  // what is the observed order of this single-steo methods 
  //====================
  /* SAM_LISTING_END_2 */

  // Printing results
  std::cout << "\n" << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "         Convergence of two-stage Radau Method           "
            << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "--------------------- RESULTS ---------------------------"
            << std::endl;
  std::cout << "Iteration"
            << "\t| Nsteps"
            << "\t| error"
            << "\t\t| rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < nIter; k++) {
    std::cout << k << "\t"
              << "\t|" << 10 * std::pow(2, k) << "\t\t|" << max_norm_errors[k];
    if (k > 0) {
      std::cout << "\t|" << rates[k - 1];
    }
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "Average rate of convergence: " << avg_rate << "\n" << std::endl;
}

}  // namespace RadauThreeTimestepping
