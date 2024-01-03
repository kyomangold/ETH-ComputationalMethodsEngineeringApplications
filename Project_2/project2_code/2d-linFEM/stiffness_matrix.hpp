#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include "grad_shape.hpp"
#include "coordinate_transform.hpp"
#include "integrate.hpp"
#include "shape.hpp"
#include "constant_function.hpp"

//----------------compMatrixBegin----------------
//! Evaluate the stiffness matrix on the triangle spanned by
//! the points (a, b, c).
//!
//! Here, the stiffness matrix A is a 3x3 matrix
//! 
//! $$A_{ij} = \int_{K} ( sigma(x,y) \nabla \lambda_i^K(x, y) \cdot  \nabla \lambda_j^K(x, y)\; + r \lambda_i^K(x,y)\lambda_j^K(x,y) )dV$$
//! 
//! where $K$ is the triangle spanned by (a, b, c).
//!
//! @param[out] stiffnessMatrix should be a 3x3 matrix
//!                        At the end, will contain the integrals above.
//!
//! @param[in] a the first corner of the triangle
//! @param[in] b the second corner of the triangle
//! @param[in] c the third corner of the triangle
//! @param[in] sigma the function sigma as in the exercise
//! @param[in] r the parameter r as in the exercise
template<class MatrixType, class Point>
void computeStiffnessMatrix(MatrixType& stiffnessMatrix,
                            const Point& a,
                            const Point& b,
                            const Point& c,
			    const std::function<double(double, double)>& sigma = constantFunction,
                            const double r=0)
{
    Eigen::Matrix2d coordinateTransform = makeCoordinateTransform(b - a, c - a);
    double volumeFactor = std::abs(coordinateTransform.determinant());
    
    Eigen::Matrix2d elementMap = coordinateTransform.inverse().transpose();
    //// CMEA_START_TEMPLATE
    for (int i = 0; i < 3; ++i) {
        for (int j = i; j < 3; ++j) {

	  stiffnessMatrix(i, j) = integrate([&](double x, double y) {
	      Eigen::Vector2d gradLambdaI = elementMap * gradientLambda(i, x, y);
	      Eigen::Vector2d gradLambdaJ = elementMap * gradientLambda(j, x, y);

	      auto lambdaI = lambda(i, x, y);
	      auto lambdaJ = lambda(j, x, y);

	      Eigen::Vector2d z = coordinateTransform * Eigen::Vector2d(x, y)
	                          + Eigen::Vector2d(a(0), a(1));
	                     
	      return volumeFactor * (sigma(z(0),z(1))*gradLambdaI.dot(gradLambdaJ)
				     +r*lambdaI*lambdaJ);
	      
            });
        }
    }
    
    // Make symmetric (we did not need to compute these value above)
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < i; ++j) {
            stiffnessMatrix(i, j) = stiffnessMatrix(j, i);
        }
    }
    //// CMEA_END_TEMPLATE
    
}
//----------------compMatrixEnd----------------
