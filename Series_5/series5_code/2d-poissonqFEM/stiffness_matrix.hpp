#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include "grad_shape.hpp"
#include "coordinate_transform.hpp"
#include "integrate.hpp"
#include "shape.hpp"

//----------------compMatrixBegin----------------
//! Evaluate the stiffness matrix on the triangle spanned by
//! the points (a, b, c).
//!
//! Here, the stiffness matrix A is a 6x6 matrix
//! 
//! $$A_{ij} = \int_{K} ( \nabla \lambda_i^K(x, y) \cdot  \nabla \lambda_j^K(x, y)\; dV$$
//! 
//! where $K$ is the triangle spanned by (a, b, c).
//!
//! @param[out] stiffnessMatrix should be a 6x6 matrix
//!                        At the end, will contain the integrals above.
//!
//! @param[in] a the first corner of the triangle
//! @param[in] b the second corner of the triangle
//! @param[in] c the third corner of the triangle
template<class MatrixType, class Point>
void computeStiffnessMatrix(MatrixType& stiffnessMatrix,
                            const Point& a,
                            const Point& b,
                            const Point& c)
{
    Eigen::Matrix2d coordinateTransform = makeCoordinateTransform(b - a, c - a);
    double volumeFactor = std::abs(coordinateTransform.determinant());
    Eigen::Matrix2d elementMap = coordinateTransform.inverse().transpose();
    //// CMEA_START_TEMPLATE
    for (int i = 0; i < 6; ++i) {
        for (int j = i; j < 6; ++j) {

            stiffnessMatrix(i, j) = integrate([&](double x, double y) {
                Eigen::Vector2d gradphiI = elementMap * gradientShapefun(i, x, y);
                Eigen::Vector2d gradphiJ = elementMap * gradientShapefun(j, x, y);
               
                return volumeFactor * gradphiI.dot(gradphiJ);

            });
        }
    }
    
    // Make symmetric (we did not need to compute these value above)
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < i; ++j) {
            stiffnessMatrix(i, j) = stiffnessMatrix(j, i);
        }
    }
    //// CMEA_END_TEMPLATE

    
}
//----------------compMatrixEnd----------------
