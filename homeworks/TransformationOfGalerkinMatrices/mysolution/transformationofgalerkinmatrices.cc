/**
 * @file
 * @brief NPDE homework TransformationOfGalerkinMatrices code
 * @author Erick Schulz
 * @date 01/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "transformationofgalerkinmatrices.h"

#include <cassert>

namespace TransformationOfGalerkinMatrices {

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::Triplet<double>> transformCOOmatrix(
    const std::vector<Eigen::Triplet<double>> &A) {
  std::vector<Eigen::Triplet<double>> A_t{};  // return value

  // First step: find the size of the matrix by searching the maximal
  // indices. Depends on the assumption that no zero rows/columns occur.
  int rows_max_idx = 0, cols_max_idx = 0;
  for (const Eigen::Triplet<double> &triplet : A) {
    rows_max_idx =
        (triplet.row() > rows_max_idx) ? triplet.row() : rows_max_idx;
    cols_max_idx =
        (triplet.col() > cols_max_idx) ? triplet.col() : cols_max_idx;
  }
  int n_rows = rows_max_idx + 1;
  int n_cols = cols_max_idx + 1;

  // Make sure we deal with a square matrix
  assert(n_rows == n_cols);
  // The matrix size must have even parity
  assert(n_cols % 2 == 0);

  int N = n_cols;      // Size of (square) matrix
  int M = n_cols / 2;  // Half the size
  //====================
  // Your code goes here
  // implement a C++ function 
  // that takes the Galerkin matrix A in triplet formatee and returns A~ also in triplet formate. 
  // it can be taken granted that the argument matrix has no rows or columns with all zero entries, 
  // which makes it possible to determine the matrix size from the triplet information. 
  // the emplace_back() method adds an element to a vector 
  for (const Eigen::Triplet<double>&it :A){
    i = it.row()+1; 
    j = it.col()+1; 
    if (i%2==0 && j%2==0){
      A_t.emplace_back(i/2-1, j/2-1, it.value());
      A_t.emplace_back(i/2-1, j/2+M-1, -it.value()); 
      A_t.emplace_back(i/2+M-1, j/2-1, -it.value()); 
      A_t.emplace_back(i/2+M-1, j/2+M-1, it.value()); 
    } else if (i%2==0 && j%2!=0){
      A_t.emplace_back(i/2-1, (j+1)/2-1, it.value()); 
      A_t.emplace_back(i/2-1, (j+1)/2+M-1, it.value()); 
      A_t.emplace_back(i/2+M-1, (j+1)/2+M-1, -it.value()); 
      A_t.emplace_back(i/2+M-1, (j+1)/2-1, -it.value());
    } else if (i%2!=0 && j%2==0){
      A_t.emplace_back((i+1)/2-1, j/2-1, it.value()); 
      A_t.emplace_back((i+1)/2+M-1, j/2+M-1, -it.value()); 
      A_t.emplace_back((i+1)/2+M-1, j/2-1,it.value()); 
      A_t.emplace_back((i+1)/2-1, j/2+M-1, -it.value()); 
      
    } else if (i%2!=0 && j%2!=0){
      A_t.emplace_back((i+1)/2-1, (j+1)/2-1, it.value()); 
      A_t.emplace_back((i+1)/2+M-1, (j+1)/2+M-1, it.value()); 
      A_t.emplace_back((i+1)/2-1, (j+1)/2+M-1, it.value()); 
      A_t.emplace_back((i+1)/2+M-1, (j+1)/2-1, it.value());
    }
  }



  //====================
  return A_t;
}
/* SAM_LISTING_END_1 */

}  // namespace TransformationOfGalerkinMatrices
