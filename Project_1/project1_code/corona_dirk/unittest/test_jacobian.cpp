#include <gtest/gtest.h>
#include <Eigen/Core>
#include "../coronaoutbreak.hpp"

#define NEAR_TOLERANCE 1e-7
TEST(TestJacobian, FirstColCheck) {
    CoronaOutbreak outbreak;
    int dim = 4;
    Eigen::MatrixXd JFu(dim,dim);
    Eigen::VectorXd u(dim);
    u << 0, 1, 0, 0;

    outbreak.computeJF(JFu,1,u);

    ASSERT_NEAR(JFu(0,0), -0.0002885, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,0), 0.00029849, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,0), 0.0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(3,0), 0.0, NEAR_TOLERANCE);

    ASSERT_NEAR(JFu(0,1), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,1), -0.500035, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,1), 0.5, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(3,1), 1.5e-5, NEAR_TOLERANCE);
    
    ASSERT_NEAR(JFu(0,2), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,2), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,2), -0.50251999, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(3,2), 0.5, NEAR_TOLERANCE);

    ASSERT_NEAR(JFu(0,3), 0.002, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,3), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,3), 0, NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(3,3), -0.00202, NEAR_TOLERANCE);
}

TEST(TestJacobian, SecondAndThirdColCheck) {
    CoronaOutbreak outbreak;
    int dim = 4;
    Eigen::MatrixXd JFu(dim,dim), Je(dim,dim);
    Eigen::VectorXd u(dim);
    u << 1, 0, 0, 0;

    outbreak.computeJF(JFu,1,u);

    Je.row(0) << 9.9999999999999991e-06, -0.00029849999999999999, -0.0014925, 0.002;
    Je.row(1) <<  0,-0.49973650000000003,0.0014925, 0;
    Je.row(2) << 0, 0.5, -0.50251999999999997, 0;
    Je.row(3) << 0, 1.5e-05, 0.5, -0.00202;

    ASSERT_NEAR(JFu(0,0),Je(0,0), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(0,1),Je(0,1), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(0,2),Je(0,2), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(0,3),Je(0,3), NEAR_TOLERANCE);

    ASSERT_NEAR(JFu(1,0),Je(1,0), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,1),Je(1,1), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,2),Je(1,2), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,3),Je(1,3), NEAR_TOLERANCE);

    ASSERT_NEAR(JFu(2,0),Je(2,0), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,1),Je(2,1), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,2),Je(2,2), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,3),Je(2,3), NEAR_TOLERANCE);

    ASSERT_NEAR(JFu(3,0),Je(3,0), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(3,1),Je(3,1), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(3,2),Je(3,2), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(3,3),Je(3,3), NEAR_TOLERANCE);
}


TEST(TestJacobian, LastColCheck) {
    CoronaOutbreak outbreak(0,0,0.13,0,0,0.1,0.98,0,0);
    int dim = 4;
    Eigen::MatrixXd JFu(dim,dim), Je(dim,dim);
    Eigen::VectorXd u(dim);
    u << 0, 0, 0, 0;

    outbreak.computeJF(JFu,1,u);

    Je.row(0) << -0.88, 0, 0, 0.13;
    Je.row(1) <<  0, -0.98, 0, 0;
    Je.row(2) << 0, 0, -0.98, 0;
    Je.row(3) << 0, 0, 0, -1.11;

    ASSERT_NEAR(JFu(0,0),Je(0,0), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(0,1),Je(0,1), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(0,2),Je(0,2), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(0,3),Je(0,3), NEAR_TOLERANCE);

    ASSERT_NEAR(JFu(1,0),Je(1,0), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,1),Je(1,1), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,2),Je(1,2), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(1,3),Je(1,3), NEAR_TOLERANCE);

    ASSERT_NEAR(JFu(2,0),Je(2,0), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,1),Je(2,1), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,2),Je(2,2), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(2,3),Je(2,3), NEAR_TOLERANCE);

    ASSERT_NEAR(JFu(3,0),Je(3,0), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(3,1),Je(3,1), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(3,2),Je(3,2), NEAR_TOLERANCE);
    ASSERT_NEAR(JFu(3,3),Je(3,3), NEAR_TOLERANCE);
}
