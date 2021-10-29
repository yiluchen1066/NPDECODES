/**
 * @file matode.cc
 * @brief NPDE homework MatODE code
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>

namespace MatODE {

/* SAM_LISTING_BEGIN_3 */
Eigen::MatrixXd eeulstep(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Y0,
                         double h) {
  //====================
  // Your code goes here
  return Y0+h*A*Y0; 
  //====================
  
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::MatrixXd ieulstep(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Y0,
                         double h) {
  //====================
  // Your code goes here
  int n=A.rows();
  return (Eigen::MatrixXd::Identity(n,n)-h*A).partialPivLu().solve(Y0);
  //====================
  
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
Eigen::MatrixXd impstep(const Eigen::MatrixXd& A, const Eigen::MatrixXd& Y0,
                        double h) {
  //====================
  // Your code goes here
  //====================
  int n=A.nrows();
  return (Eigen::MatrixXd::Identity(n,n)-0.5*h*A).partialPivLu().solve(Y0+0.5*h*A*Y0);
  
}
/* SAM_LISTING_END_5 */

}  // namespace MatODE
