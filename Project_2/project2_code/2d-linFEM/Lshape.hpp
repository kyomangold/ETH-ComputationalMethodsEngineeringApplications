#pragma once
#include <iostream>
#include <Eigen/Core>
#include <igl/readMESH.h>

#include "writer.hpp"
#include "fem_solve.hpp"

double f_lshape(double x, double y) {
    // return 0;
	return sin(M_PI * y) * sin(M_PI * y);
}

double g_lshape(double x, double y) {
    /* double r = std::sqrt(x * x + y * y);

    double theta = std::atan2(y, x);

    // Adjust for the region where theta < 0
    if (theta < 0) {
        theta += 2 * M_PI;
    }

    return std::pow(r, 2.0 / 3.0) * std::sin(2 * theta / 3); */
    return x*x*x + y;
}

Eigen::Vector2d g_grad_lshape(double x, double y) {
    /* double r = std::sqrt(x * x + y * y);

    double theta = std::atan2(y, x);

    // Adjust for the region where theta < 0
    if (theta < 0) {
        theta += 2 * M_PI;
    }
    Eigen::Vector2d grad;
    grad << -2.0/3.0*std::pow(r, -1.0 / 3.0) * std::sin(theta / 3), 2.0/3.0*std::pow(r, -1.0 / 3.0) * std::cos(theta / 3);
    return grad; */
    Eigen::Vector2d grad;
    grad << 3*x*x, 1;
    return grad;
}

void solveL(double r, const std::function<double(double, double)>& sigma) {
    std::cout << "Solving L-shape" << std::endl;
    Vector u;

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi triangles;
    Eigen::MatrixXi tetrahedra;

    igl::readMESH(CMEA_DATA_PATH "Lshape_5.mesh", vertices, tetrahedra, triangles);

    solveFiniteElement(u, vertices, triangles, f_lshape, sigma, g_lshape, r);

    writeToFile("Lshape_values.txt", u);
    writeMatrixToFile("Lshape_vertices.txt", vertices);
    writeMatrixToFile("Lshape_triangles.txt", triangles);

}
