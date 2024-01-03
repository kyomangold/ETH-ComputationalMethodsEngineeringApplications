#include <gtest/gtest.h>
#include <Eigen/Core>
#include "../dirksolver.hpp"
#include <iostream>

#define NEAR_TOLERANCE 1e-7

TEST(TestGFunctions, G1Check) {
    const double dt = 2*sqrt(3)/(1+sqrt(3));
    CoronaOutbreak outbreak(2,4,1,1,1,1,0,0,0,1);
    Eigen::VectorXd y(4), u(4), G(4);

    y << 1, 2, 3, 4;
    u << 4, 8, 15, 19;
    G.setConstant(0.0);

    DIRKSolver solver(outbreak);
    solver.computeG1(G, y, 0, u, dt);

    ASSERT_NEAR(G(0), 1, NEAR_TOLERANCE);
    ASSERT_NEAR(G(1), 11, NEAR_TOLERANCE);
    ASSERT_NEAR(G(2), 11, NEAR_TOLERANCE);
    ASSERT_NEAR(G(3), 14, NEAR_TOLERANCE);
}


TEST(TestGFunctions, G2Check) {
    const double dt = 2*sqrt(3)/(1+sqrt(3));
    CoronaOutbreak outbreak(2,4,1,1,1,1,0,0,0,1);
    Eigen::VectorXd y1(4), u(4), y(4), G(4);

    y1 << 1, 2, 3, 4;
    u << 4, 8, 15, 19;
    y << 5, 0, 0, 4;
    G.setConstant(0.0);

    DIRKSolver solver(outbreak);
    solver.computeG2(G, y, 0, u, dt, y1);
    Eigen::Vector4d Ge(12. - 2.*dt, -2. + 5. * dt, 17 - dt, 13 - dt);

    ASSERT_NEAR(G(0),Ge(0), NEAR_TOLERANCE);
    ASSERT_NEAR(G(1),Ge(1), NEAR_TOLERANCE);
    ASSERT_NEAR(G(2),Ge(2), NEAR_TOLERANCE);
    ASSERT_NEAR(G(3),Ge(3), NEAR_TOLERANCE);
}


TEST(TestNewtonMethod, NewtonY1Check) {
    CoronaOutbreak outbreak;
    Eigen::VectorXd u(4), uOut(4);

    u << 4, 8, 15, 19;

    DIRKSolver solver(outbreak);
	solver.newtonSolveY1(uOut, u, 0.01, 1, 1e-10, 100);

    ASSERT_NEAR(uOut(0), 3.999521057255166, NEAR_TOLERANCE);
    ASSERT_NEAR(uOut(1), 7.9693515344298502, NEAR_TOLERANCE);
    ASSERT_NEAR(uOut(2), 14.972088014857869, NEAR_TOLERANCE);
    ASSERT_NEAR(uOut(3), 19.058737881152023, NEAR_TOLERANCE);
}


TEST(TestNewtonMethod, NewtonY2Check) {
    CoronaOutbreak outbreak;
    Eigen::VectorXd u(4), uOut(4), y(4);
    u << 4, 8, 15, 19;
    y << 16, 23, 42, 10;

    DIRKSolver solver(outbreak);
	solver.newtonSolveY2(uOut, u, y, 0.01, 1, 1e-10, 100);

    ASSERT_NEAR(uOut(0), 4.0058230106186841, NEAR_TOLERANCE);
    ASSERT_NEAR(uOut(1), 8.0290955994451281, NEAR_TOLERANCE);
    ASSERT_NEAR(uOut(2), 15.027563089642353, NEAR_TOLERANCE);
    ASSERT_NEAR(uOut(3), 18.937829649864085, NEAR_TOLERANCE);
}

TEST(TestDirkSolver, DirkSolverCheckSizes) {
    std::vector<std::vector<double> > u(5);
    std::vector<double> time(201, 0);
    
    for(int i = 0; i<u.size(); i++){
        u[i].resize(201,0);
        u[i][0] = 0.;
    }

    CoronaOutbreak outbreak;
    DIRKSolver dirk(outbreak);
    dirk.solve(u, time, 1, 200);

    ASSERT_EQ(u[0].size(), 201);
    ASSERT_EQ(u[1].size(), 201);
    ASSERT_EQ(u[2].size(), 201);
    ASSERT_EQ(u[3].size(), 201);
    ASSERT_EQ(u[4].size(), 201);
    ASSERT_EQ(time.size(), 201);
}

TEST(TestDirkSolver, DirkSolverWroteLast) {
    std::vector<std::vector<double> > u(5);
    std::vector<double> time(201, 0);

    for(int i = 0; i<u.size(); i++){
        u[i].resize(201,0);
        u[i][0] = 1.;
    }

    CoronaOutbreak outbreak;
    DIRKSolver dirk(outbreak);
    dirk.solve(u, time, 1, 200);

    EXPECT_NE(u[0][200], 0);
    EXPECT_NE(u[1][200], 0);
    EXPECT_NE(u[2][200], 0);
    EXPECT_NE(u[3][200], 0);
    EXPECT_NE(u[4][200], 0);
    ASSERT_NE(time[200], 0);
}


TEST(TestDirkSolver, DirkSolverCheck) {
    const double T = 42;
    const int N = 1138;
    std::vector<std::vector<double> > u(5);
    std::vector<double> time(N + 1, 0);

    for(int i = 0; i<u.size(); i++){
        u[i].resize(N + 1,0);
    }

    u[0][0] = 1123;
    u[1][0] = 6536;
    u[2][0] = 5321;
    u[3][0] = 3321;
    u[4][0] = 0;

    CoronaOutbreak outbreak;
    DIRKSolver dirk(outbreak);
    dirk.solve(u, time, T, N);

    ASSERT_NEAR(u[0][N-2], 298.01271226267249, NEAR_TOLERANCE);
    ASSERT_NEAR(u[1][N-2], 93.517705913615629, NEAR_TOLERANCE);
    ASSERT_NEAR(u[2][N-2], 87.660290107132013, NEAR_TOLERANCE);
    ASSERT_NEAR(u[3][N-2], 15740.037176309437, NEAR_TOLERANCE);
    ASSERT_NEAR(u[4][N-2], 82.040496575054078, NEAR_TOLERANCE);
}

