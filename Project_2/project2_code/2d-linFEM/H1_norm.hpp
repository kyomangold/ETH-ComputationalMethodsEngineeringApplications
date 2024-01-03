#pragma once
#include <functional>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "coordinate_transform.hpp"
#include "integrate.hpp"
#include "shape.hpp"

//! Computes the H^1 differences between 
//! u1 (considered as coefficients for the shape functions)
//! and u2grad.
double computeH1Difference(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& triangles,
    const Eigen::VectorXd& u1,
			   const std::function<Eigen::Vector2d(double, double)>& u2grad) 
{
    const int numberOfElements = triangles.rows();

    double error = 0;
    for (int i = 0; i < numberOfElements; ++i) {
        auto& indexSet = triangles.row(i);

        const int i0 = indexSet(0);
        const int i1 = indexSet(1);
        const int i2 = indexSet(2);

        const auto& a = vertices.row(i0);
        const auto& b = vertices.row(i1);
        const auto& c = vertices.row(i2);

        auto coordinateTransform = makeCoordinateTransform(b - a, c - a);
        auto volumeFactor = std::abs(coordinateTransform.determinant());
	Eigen::Matrix2d elementMap = coordinateTransform.inverse().transpose();

	
        error += integrate([&](double x, double y) {
            Eigen::Vector2d z = coordinateTransform * Eigen::Vector2d(x, y) + Eigen::Vector2d(a(0), a(1));

	    Eigen::Vector2d approximate_grad = u1(i0) * elementMap * gradientLambda(0, x, y) + u1(i1) * elementMap * gradientLambda(1, x, y) + u1(i2) * elementMap * gradientLambda(2, x, y);
	    Eigen::Vector2d difference_grad = u2grad(z(0), z(1)) - approximate_grad;
            return difference_grad.dot(difference_grad) * volumeFactor;
        });
    }

    return std::sqrt(error);
}
