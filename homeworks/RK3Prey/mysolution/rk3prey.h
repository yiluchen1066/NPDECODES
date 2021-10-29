/**
 * @file rk3prey_main.cc
 * @brief NPDE homework RK3Prey code
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>
#include <vector>

namespace RK3Prey {

// Butcher Tableau based Runge-Kutta explicit solver for autonomous ODEs
/* SAM_LISTING_BEGIN_0 */
class RKIntegrator {
 public:
  RKIntegrator(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) : A_(A), b_(b), s_(b.size()){
    //====================
    // Your code goes here
    //====================
    assert(A.cols() == A.rows() && "Matrix must be square."); 
    assert(A.cols() == b.size() && "Incompatible matrix/vector size.")
  }

  // Explicit Runge-Kutta numerical integrator
  template <class Function>
  std::vector<Eigen::VectorXd> solve(Function &&f, double T,
                                     const Eigen::VectorXd &y0, int M) const;

 private:
  //====================
  // Your code goes here
  const Eigen::MatrixXd A_; 
  const Eigen::VectorXd b_; 
  int s_; // size of butcher tableau
  //====================
};
/* SAM_LISTING_END_0 */

/* Solves an autonomous ODE y' = f(y), y(0) = y0, using a
 * RK scheme from the Butcher tableau provided by the
 * constructor. Performs N equidistant steps up to time T */
/* SAM_LISTING_BEGIN_1 */
template <typename Function>
std::vector<Eigen::VectorXd> RKIntegrator::solve(Function &&f, double T,
                                                 const Eigen::VectorXd &y0,
                                                 int M) const { 
  int dim = y0.size();  // dimension
  double h = T / M;     // step size
  std::vector<Eigen::VectorXd> sol;
  sol.reserve(M + 1);
  sol.push_back(y0); 

  //====================
  // Your code goes here
  // initialize the vector 

  Eigen::VectorXd incr(dim);
  std::vector<Eigen::VectorXd> k; 
  k.reserve(s_);  
  

  Eigen::VectorXd step(dim); 

  for (int iter =0; iter<M; iter++){
    step.setZero(); 
    k.clear(); 
    k.push_back(f(sol.at(iter))); 
    step = step + b_(0) + k(0); 
    for (int i=1; i<s_+1; i++){
      for (int j=1; j<i; j++){
        incr = incr + A_(i,j)*k(j); 
      }
      k.push_back(f(sol.at(iter)+h*incr)); 
      step = step+h*b_(i)*k.back(); 
    }
    sol.push_back(f(sol.at(iter)+h*step)); 
  }


  //Rugge-Kutta looping tools
  
  //stepping 
  // M number of equidistant timesteps 
  // f the right hand sirde field 
  // the final tiem T 
   
  

  //====================
  return sol;
}
/* SAM_LISTING_END_1 */

}  // namespace RK3Prey
