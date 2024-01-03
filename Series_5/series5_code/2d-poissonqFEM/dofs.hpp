#ifndef QDOFS_HPP
#define QDOFS_HPP

#pragma once
#include <set>
#include <Eigen/Core>
#include <igl/edges.h>
#include <igl/vertex_triangle_adjacency.h>

class QDofs{
public:
  // list of edges connected to given vertex
  std::map<int, std::set<int>> vertex2edge;
  Eigen::MatrixXd nodes;

private:
  const Eigen::MatrixXd& vertices_;
  const Eigen::MatrixXi& triangles_;
  Eigen::MatrixXi edges_;
  Eigen::MatrixXi qdof_;


public:  
  //! QDofs constructor
  //!
  //! @param[in] vertices a list of triangle vertices
  //! @param[in] triangles a list of triangles
  QDofs(const Eigen::MatrixXd& vertices,
	const Eigen::MatrixXi& triangles)
    : vertices_(vertices)
    , triangles_(triangles){

    set_quadraticDofs();
  }

  int get_dofs(Eigen::MatrixXi& qdof) const{
    qdof = qdof_;
    return edges_.rows()+vertices_.rows();
  }

  Eigen::MatrixXd get_vertices() const{
    return vertices_;
  }

private:  
  //! Create matrix containing the global indices of the degrees of freedom
  //! associated to quadratic FEM. First numbers vertices and then edges.
  //! First three columns are vertices, the last three columns are the edges.
  void set_quadraticDofs(){

    // get list of mesh edges
    igl::edges(triangles_, edges_);

    // resize dof_ so it has right length and 6 columns
    qdof_.resize(triangles_.rows(), 6);
    // set frist 3 columns to match the vertices defining the triangles
    qdof_.block(0,0, triangles_.rows(), 3) = triangles_;

    int Ne = edges_.rows();
    int Nv = vertices_.rows();

    // resize matrix with nodes coordinates
    nodes.resize(Ne+Nv,3);
    // set first Nv rows to be the vertices
    nodes.block(0,0, vertices_.rows(), 3) = vertices_;

    // get list of elements to which each vertex is connected
    std::vector<std::vector<int> > VF;
    std::vector<std::vector<int> > VFi;
    igl::vertex_triangle_adjacency(Nv,triangles_, VF, VFi);
    
    // traverse the edges
    for(int i=0; i<Ne; i++){
      // get vertices attached to current edge
      auto edgevert = edges_.row(i);

      // add midpoint to matrix of nodes coordinates
      nodes.row(i+Nv) = (vertices_.row(edgevert(0)) + vertices_.row(edgevert(1)) )/2.;
      
      // get list of elements connect to those vertices
      auto tria_vert0 = VF[edgevert(0)];
      auto tria_vert1 = VF[edgevert(1)];
      // get list of local index of those vertices in said elements
      auto tria_vertid0 = VFi[edgevert(0)];
      auto tria_vertid1 = VFi[edgevert(1)];
      // traverse elements connected to first vertex
      for(int k0 = 0; k0< tria_vert0.size(); k0++){
	auto tria0 = tria_vert0[k0];
	// traverse elements connected to second element
	for(int k1 = 0; k1< tria_vert1.size(); k1++){
	  auto tria1 = tria_vert1[k1];
	  // if the elements match, it means we have found a triangle
	  // that contains the current edge
	  if( tria0 == tria1){
	    int id = tria_vertid0[k0] + tria_vertid1[k1] ;
	    // fill the matrix dof_ accordingly
	    if (id == 1){qdof_(tria0, 3) = Nv+i;} // edge 0, dof 3
	    else if (id == 3){qdof_(tria0, 4) = Nv+i;} // edge 1, dof 4
	    else if (id == 2){qdof_(tria0, 5) = Nv+i;} // edge 2, dof 5
	    else{
	      std::cout << " error!" << std::endl;
	    }
	    
	  }// end found element
	  
	} //looping on second vertex
      } //looping on first vertex 

            // add edge to the vertex2edge vector
      for(int k=0; k<2; k++){
	// get set of each vertex
	auto auxset = vertex2edge.find(edgevert(k));
	// if vertex index is already there 
	if(auxset!= vertex2edge.end()){
	  // add it to its set
	  (*auxset).second.insert((Nv+i));
	}
	else{ // if not, create set
	  std::set<int> vset;
	  vset.insert((Nv+i));
	  vertex2edge.insert(std::pair<int, std::set<int>>(edgevert(k), vset));
	}
      }
      
      
    }// looping on edges  

  }
  
};
#endif




