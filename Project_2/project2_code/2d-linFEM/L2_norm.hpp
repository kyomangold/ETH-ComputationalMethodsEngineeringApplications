#pragma once
#include <functional>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "coordinate_transform.hpp"
#include "integrate.hpp"
#include "shape.hpp"

//! Computes the L^2 differences between 
//! u1 (considered as coefficients for the shape functions)
//! and u2.
double computeL2Difference(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& triangles,
    const Eigen::VectorXd& u1,
    const std::function<double(double, double)>& u2) 
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

        error += integrate([&](double x, double y) {
            Eigen::Vector2d z = coordinateTransform * Eigen::Vector2d(x, y) + Eigen::Vector2d(a(0), a(1));

            double approximateValue = u1(i0) * lambda(0, x, y) + u1(i1) * lambda(1, x, y) + u1(i2) * lambda(2, x, y);

            return std::pow(std::abs(u2(z(0), z(1)) - approximateValue), 2) * volumeFactor;
        });
    }

    return std::sqrt(error);
}
