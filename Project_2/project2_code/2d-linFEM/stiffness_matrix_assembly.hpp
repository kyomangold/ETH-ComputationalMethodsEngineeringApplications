#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

#include "stiffness_matrix.hpp"

//! Sparse Matrix type. Makes using this type easier.
typedef Eigen::SparseMatrix<double> SparseMatrix;

//! Used for filling the sparse matrix.
typedef Eigen::Triplet<double> Triplet;

//----------------AssembleMatrixBegin----------------
//! Assemble the stiffness matrix
//! for the linear system
//!
//! @param[out] A will at the end contain the Galerkin matrix
//! @param[in] vertices a list of triangle vertices
//! @param[in] triangles a list of triangles
//! @param[in] sigma the function sigma as in the exercise
//! @param[in] r the parameter r as in the exercise
template<class Matrix>
void assembleStiffnessMatrix(Matrix& A, const Eigen::MatrixXd& vertices,
			     const Eigen::MatrixXi& triangles,
			     const std::function<double(double, double)>& sigma = constantFunction,
			     double r=0)
{
    
    const int numberOfElements = triangles.rows();
    A.resize(vertices.rows(), vertices.rows());
    
    std::vector<Triplet> triplets;

    triplets.reserve(numberOfElements * 3 * 3);
    //// CMEA_START_TEMPLATE
    for (int i = 0; i < numberOfElements; ++i) {
        auto& indexSet = triangles.row(i);

        const auto& a = vertices.row(indexSet(0));
        const auto& b = vertices.row(indexSet(1));
        const auto& c = vertices.row(indexSet(2));

        Eigen::Matrix3d stiffnessMatrix;
        computeStiffnessMatrix(stiffnessMatrix, a, b, c, sigma, r);

        for (int n = 0; n < 3; ++n) {
            for (int m = 0; m < 3; ++m) {
                auto triplet = Triplet(indexSet(n), indexSet(m), stiffnessMatrix(n, m));
                triplets.push_back(triplet);
            }
        }
    }
    //// CMEA_END_TEMPLATE
    A.setFromTriplets(triplets.begin(), triplets.end());
}
//----------------AssembleMatrixEnd----------------
