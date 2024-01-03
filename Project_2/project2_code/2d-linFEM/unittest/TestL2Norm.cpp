#include <gtest/gtest.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <gtest/gtest.h>
#include "integrate.hpp"
#include "L2_norm.hpp"

TEST(TestL2Norm, ConstantTest) {
    // In this test, we check that the L2-norm function
    // is able to get a simple constant right

    const double constant = 4.0;
    auto f = [=](double, double) {
        return constant;
    };

    // We only add one triangle:

    Eigen::MatrixXd vertices(3, 3);
    vertices << 0, 0, 0,
        1, 0, 0,
        0, 1, 0;

    Eigen::MatrixXi triangles(1, 3);
    triangles << 0, 1, 2;

    Eigen::VectorXd u(3);
    u << constant, constant, constant;

    const double error = computeL2Difference(vertices, triangles, u, f);

    ASSERT_EQ(0.0, error) << "L2-norm does not give 0 = ||C-C||_2";
}


TEST(TestL2Norm, ShapeFunctionTest) {
    // In this test, we check that the L2-norm function
    // is able to destinguish the shape functions


    for (int i = 0; i < 3; ++i) {
        // We only add one triangle:

        Eigen::MatrixXd vertices(3, 3);
        vertices << 0, 0, 0,
            1, 0, 0,
            0, 1, 0;

        Eigen::MatrixXi triangles(1, 3);
        triangles << 0, 1, 2;

        Eigen::VectorXd u(3);
        u << (i == 0) , (i == 1), (i == 2);

        auto f = [i](double x, double y) {
            return lambda(i, x, y);
        };

        const double error = computeL2Difference(vertices, triangles, u, f);

        ASSERT_EQ(0.0, error) << "L2-norm does not give 0 = ||lambda_" << i <<"-lambda_" << i << "||_2";
    }
}



