#pragma once
#include <Eigen/Core>
#include <string>
#include <igl/readSTL.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include "stiffness_matrix_assembly.hpp"
#include "load_vector_assembly.hpp"
#include "dirichlet_boundary.hpp"
#include <chrono>

typedef Eigen::VectorXd Vector;

//----------------solveBegin----------------
//! Solve the FEM system.
//!
//! @param[out] u will at the end contain the FEM solution.
//! @param[in] vertices list of triangle vertices for the mesh
//! @param[in] triangles list of triangles (described by indices)
//! @param[in] f the RHS f (as in the exercise)
//! @param[in] g the boundary value (as in the exercise)
//! @param[in] sigma the function sigma as in the exercise
//! @param[in] r the parameter r from the lecture notes
//! return number of degrees of freedom (without the boundary dofs)
int solveFiniteElement(Vector& u,
		       const Eigen::MatrixXd& vertices,
		       const Eigen::MatrixXi& triangles,
		       const std::function<double(double, double)>& f,
		       const std::function<double(double, double)>& sigma,		       
		       const std::function<double(double, double)>& g,
		       double r)
{
  SparseMatrix A(vertices.rows(), vertices.rows());
  auto startAssembly = std::chrono::high_resolution_clock::now();
    //// CMEA_START_TEMPLATE
    assembleStiffnessMatrix(A, vertices, triangles, sigma, r);
    //// CMEA_END_TEMPLATE
    auto endAssembly = std::chrono::high_resolution_clock::now();
    Vector F;
    //// CMEA_START_TEMPLATE
    assembleLoadVector(F, vertices, triangles, f);
    //// CMEA_END_TEMPLATE
    
    u.resize(vertices.rows());
    u.setZero();
    Eigen::VectorXi interiorVertexIndices;

    // set Dirichlet Boundary conditions
    //// CMEA_START_TEMPLATE
    
    setDirichletBoundary(u, interiorVertexIndices, vertices, triangles, g);
    
    F -= A * u;
    //// CMEA_END_TEMPLATE
    

    std::cout << "Assembly used " << (std::chrono::duration_cast<std::chrono::microseconds>(endAssembly-startAssembly).count() / (1000*1000.)) << " seconds\n";
    SparseMatrix AInterior;

    igl::slice(A, interiorVertexIndices, interiorVertexIndices, AInterior);
    Eigen::SimplicialLDLT<SparseMatrix> solver;

    Vector FInterior;

    igl::slice(F, interiorVertexIndices, FInterior);
    auto startSolve = std::chrono::high_resolution_clock::now();
    //initialize solver for AInterior
    //// CMEA_START_TEMPLATE
    solver.compute(AInterior);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Could not decompose the matrix");
    }
    //// CMEA_END_TEMPLATE

    //solve interior system
    //// CMEA_START_TEMPLATE
    Vector uInterior = solver.solve(FInterior);
    auto endSolve = std::chrono::high_resolution_clock::now();
    igl::slice_into(uInterior, interiorVertexIndices, u);
    //// CMEA_END_TEMPLATE
    std::cout << "Solving used " << (std::chrono::duration_cast<std::chrono::microseconds>(endSolve-startSolve).count() / (1000*1000.0)) << " seconds\n";
    return interiorVertexIndices.size();

}
//----------------solveEnd----------------
