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
//! @param[in] dofs a list of the dofs' indices in each triangle
template<class Matrix>
void assembleStiffnessMatrix(Matrix& A, const Eigen::MatrixXd& vertices,
			     const Eigen::MatrixXi& dofs, const int& N)
{
    
    const int numberOfElements = dofs.rows();
    A.resize(N,N);
    
    std::vector<Triplet> triplets;

    triplets.reserve(numberOfElements * 6 * 10);
    //// CMEA_START_TEMPLATE
    for (int i = 0; i < numberOfElements; ++i) {
        auto& indexSet = dofs.row(i);

        const auto& a = vertices.row(indexSet(0));
        const auto& b = vertices.row(indexSet(1));
        const auto& c = vertices.row(indexSet(2));

        Eigen::MatrixXd stiffnessMatrix(6,6);
        computeStiffnessMatrix(stiffnessMatrix, a, b, c);

        for (int n = 0; n < 6; ++n) {
            for (int m = 0; m < 6; ++m) {
                auto triplet = Triplet(indexSet(n), indexSet(m), stiffnessMatrix(n, m));
                triplets.push_back(triplet);
            }
        }
    }
    //// CMEA_END_TEMPLATE
    A.setFromTriplets(triplets.begin(), triplets.end());
}
//----------------AssembleMatrixEnd----------------
