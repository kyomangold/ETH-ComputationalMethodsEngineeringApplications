#pragma once
#include <Eigen/Core>
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <igl/setdiff.h>
#include <igl/colon.h>
#include "dofs.hpp"


//! Finds the boundary and interior vertices, and evaluates 
//! the function g on the boundary vertices.
//!
//! @param[out] u should be of size <number of vertices>. At the end:
//!             u(index) = g(vertices(index))  for each boundary index.
//! @param[out] interiorDofs the list of dofs that are not on the boundary.
//! @param[in]  quadraticDofs 
//! @param[in]  g the boundary value function g.
void setDirichletBoundary(Eigen::VectorXd& u,
                          Eigen::VectorXi& interiorDofs,
                          const QDofs& quadraticDofs,
                          const std::function<double(double, double)>& g)
{
  //get dofs
  Eigen::MatrixXi qdof;
  int N = quadraticDofs.get_dofs(qdof);
  
  // Find boundary edges
  Eigen::MatrixXi boundaryEdges;
  Eigen::MatrixXi triangles = qdof.block(0,0, qdof.rows(), 3);
  igl::boundary_facets(triangles, boundaryEdges);
  
  // Find dof indices corresponding to boundary edges
  Eigen::VectorXi boundaryEdgeIndices(boundaryEdges.rows());
  auto edgemap = quadraticDofs.vertex2edge;
  for(int k=0; k<boundaryEdges.rows(); k++){
    auto edge = boundaryEdges.row(k);
    auto v0 = edge(0);
    auto v1 = edge(1);
    // get edges connected to these vertices
    auto edgesv0 = (*edgemap.find(v0)).second;
    auto edgesv1 = (*edgemap.find(v1)).second;
    // traverse them and find edge (intersection)
    for(auto e0 : edgesv0){
      for(auto e1: edgesv1){
	if(e0==e1){//found edge!
	  // mark as boundary
	  boundaryEdgeIndices(k) = e0;
	  break; break;
	}
      }
    }
  }  
 
  // Find boundary vertices
  Eigen::VectorXi boundaryVertexIndices, IA, IC;
  igl::unique(boundaryEdges, boundaryVertexIndices, IA, IC);

  // set boundary data
  auto nodes = quadraticDofs.nodes;
  for ( int i = 0; i < boundaryVertexIndices.size(); ++i) {
    const int index = boundaryVertexIndices(i);
    const auto& x = nodes.row(index);
    u (index) = g(x(0), x(1));
  }
    
  // Combine this information to get boundary dofs
  Eigen::VectorXi boundaryDofs(boundaryEdges.rows()+boundaryVertexIndices.rows() );
  boundaryDofs<< boundaryVertexIndices, boundaryEdgeIndices;  
  // Get interior dofs
  Eigen::VectorXi allIndices;
  igl::colon<int>(0, N - 1, allIndices);
  igl::setdiff(allIndices, boundaryDofs, interiorDofs, IA);  
  
}
