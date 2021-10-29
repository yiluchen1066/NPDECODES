/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "mylinearfeelementmatrix.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>

namespace ElementMatrixComputation {

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 4, 4> MyLinearFEElementMatrix::Eval(
    const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());
  // Matrix for returning element matrix
  Eigen::Matrix<double, 4, 4> elem_mat;

  //====================
  // Your code goes here
  // eval member function returns the lememt matrix for the finite elemet space and the bilinear form 
  // in the weak formulation of 2.8.4
  Eigen::Matrix<double, 4,4> elment_mat_lap; 
  Eigen::Matrix<double,4,4> element_mat_mass; 

  // define the class
  lf::uscalfe::LinearFELaplaceElementMatrix laplace_elem_builder; 
  element_mat_lap = laplace_elem_builder.Eval(*cell); 

  // computations differ depending on the type of the cell
  case lf::base:RefEl::kTria():{
    double area = 0.5*((vertices(0,1)-vertices(0,0))*(vertices(1,2)-vertices(1,0)
    )-(vertices(1,1)-vertices(1,0))*(vertices(0,2)-vertices(0,0))); 
    element_mat_mass >> 2.,1.,1.,0.,
     1.,2.,1.,0., 
     1.,1.,2.,0.,
     0.,0.,0.,0.; 
    element_mat_mass = area/12*elemet_mat_mass; 
    break; 

  } case lf::base::RefEl::Kquad():{
    double area = (vertices(0,1)-vertices(0,0))*(vertices(1,3)-vertices(1,0)); 
    element_mat_mass >> 4.,2.,1.,2.,
    2.,4.,2.,1.,
    1.,2.,4.,2.,
    2.,1.,2.,4.; 
    element_mat_mass = area/36*element_mat_mass;
    break; 
  }

  elem_mat = element_mat_mass + element_mat_lap; 

  // 
  //====================

  return elem_mat;
}
/* SAM_LISTING_END_1 */
}  // namespace ElementMatrixComputation
