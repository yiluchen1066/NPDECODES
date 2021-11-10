#include "incidencematrices.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <array>
#include <memory>

namespace IncidenceMatrices {

/** @brief Create the mesh consisting of a triangle and quadrilateral
 *         from the exercise sheet.
 * @return Shared pointer to the hybrid2d mesh.
 */
std::shared_ptr<lf::mesh::Mesh> createDemoMesh() {
  // builder for a hybrid mesh in a world of dimension 2
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);

  // Add points
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{0, 0});    // (0)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{1, 0});    // (1)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{1, 1});    // (2)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{0, 1});    // (3)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{0.5, 1});  // (4)

  // Add the triangle
  // First set the coordinates of its nodes:
  Eigen::MatrixXd nodesOfTria(2, 3);
  nodesOfTria << 1, 1, 0.5, 0, 1, 1;
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),  // we want a triangle
      std::array<lf::mesh::Mesh::size_type, 3>{
          {1, 2, 4}},  // indices of the nodes
      std::make_unique<lf::geometry::TriaO1>(nodesOfTria));  // node coords

  // Add the quadrilateral
  Eigen::MatrixXd nodesOfQuad(2, 4);
  nodesOfQuad << 0, 1, 0.5, 0, 0, 0, 1, 1;
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kQuad(),
      std::array<lf::mesh::Mesh::size_type, 4>{{0, 1, 4, 3}},
      std::make_unique<lf::geometry::QuadO1>(nodesOfQuad));

  std::shared_ptr<lf::mesh::Mesh> demoMesh_p = mesh_factory_ptr->Build();

  return demoMesh_p;
}

/** @brief Compute the edge-vertex incidence matrix G for a given mesh
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *        such as lf::mesh::hybrid2d::Mesh)
 * @return The edge-vertex incidence matrix as Eigen::SparseMatrix<int>
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<int> computeEdgeVertexIncidenceMatrix(
    const lf::mesh::Mesh &mesh) {
  // Store edge-vertex incidence matrix here
  Eigen::SparseMatrix<int, Eigen::RowMajor> G;
  using size_type = lf::mesh::Mesh::size_type; 

  //====================
  // Your code goes here
  // index() method of lf:mesh:Mesh provides the numbering of entities which underlies 
  // indexing of the entries of the incidence matrixes. 
  // lf::mesh::Mesh::Index() provides a consecutive numbering of all mesh entities of a specific co-dimension 
  // lf::mesh::Entity::SubEntities() returns an array of subentities and fix their ordering. 
  const std::size_t nnz_row=2; 
  const size_type num_edge = mesh.Numentities(1); 
  const size_type num_node = mesh.Numentities(2); 
  
  Eigen::SparseMatrix<int, Eigen::RowMajor> G(num_edge, num_node); 
  G.reserve(Eigen::VectorXi::Constant(num_edge,nnz_row)); 

  // to compute G efficiently, we iterate over all edges and 
  // check the index of the nodes at its end. 
  // this is the efficient way to do the assembly, introduced as distribute scheme
  // in class
  for(lf::mesh::Entity *edge: mesh.Entities(1)){
    nonstd::span<const lf::mesh::Entity *const> node{edge.SubEntities(1)}; // the relative codimension is 1
    size_type edge_index = mesh.Index(*edge); 

    size_type node_start_index = mesh.Index(*node[0]); 
    size_type node_end_index = mesh.Index(*node[1]); //
    // subentities returned from SubEntities() can be accessed through [] operator using their local index 
    G.coeffRef(edge_index, node_start_index) +=1; 
    G.coeffRef(edge_index, node_end_index) -=1; 
  }
  
  //====================

  return G;
}
/* SAM_LISTING_END_1 */

/** @brief Compute the cell-edge incidence matrix D for a given mesh
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *        such as lf::mesh::hybrid2d::Mesh)
 * @return The cell-edge incidence matrix as Eigen::SparseMatrix<int>
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::SparseMatrix<int> computeCellEdgeIncidenceMatrix(
    const lf::mesh::Mesh &mesh) {
  // Store cell-edge incidence matrix here
  Eigen::SparseMatrix<int, Eigen::RowMajor> D;

  //====================
  // Your code goes here
  using size_type = lf::mesh::Mesh::size_type; 
  const std::size_t nnz = 4; 
  const size_type num_cell = mesh.NumEntities(0); 
  const size_type num_edge = mesh.NumEntities(1); 
  Eigen::SparseMatrix<int, Eigen::RowMajor> D(num_cell,num_edge); 
  D.reserve(Eigen::VectorXi::Constant(num_cell,nnz)); 
  for(lf::mesh::Entity *cell: mesh.Entities(0)){
    // get edges and their orientations of a cell 
    nonstd::span<const Entity* const> edges = cell->SubEntities(1); 
    nonstd::span<const Entity* const> orientations = cell->RelativeOrientations(); 
    size_type cell_index = mesh.Index(*cell); 
    auto edge_start = edges.begin(); 
    auto orientation_start = orientation.begin(); 
    // get the index of each edge and add their orientations to D
    for(; edge_start !=edges.end() && orientation_start != orientation.end(); edge_start++, orientation++){
      size_type edge_index=mesh.Index(*edge[edge_start]);
      D.coeffRef(cell_index,edge_index) += lf::mesh::to_sign(orientation_start);  
    }
  }


  //====================

  return D;
}
/* SAM_LISTING_END_2 */

/** @brief For a given mesh test if the product of cell-edge and edge-vertex
 *        incidence matrix is zero: D*G == 0?
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *             such as lf::mesh::hybrid2d::Mesh)
 * @return true, if the product is zero and false otherwise
 */
/* SAM_LISTING_BEGIN_3 */
bool testZeroIncidenceMatrixProduct(const lf::mesh::Mesh &mesh) {
  bool isZero = false;

  //====================
  // Your code goes here
  // returns true whenever the two incidence matrices of the 2D hybrid mesh satisfy 
  // the relationship asserted in 2.6.6 
  Eigen::SparseMatrix<int> G = computeEdgeVertexIncidenceMatrix(*mesh); 
  Eigen::SparseMatrix<int> D = computeCellEdgeIncidenceMatrix(*mesh); 
  auto O = G*D; 
  
  isZero = O.norm()==0; 
  }
  //====================
  return isZero;
}
/* SAM_LISTING_END_3 */

}  // namespace IncidenceMatrices
