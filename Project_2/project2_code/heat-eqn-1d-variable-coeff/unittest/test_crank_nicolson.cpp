#include <gtest/gtest.h>

#include "crank_nicolson.hpp"
#include "boundaries.hpp"
#include "coefficient_functions.hpp"

//! Change this to see how well your code matches.
#define CRANK_NICOLSON_ABSOLUTE_TOLERANCE 1e-6

TEST(TestCrankNicolson, SizeCheck6x5) {
    const int N = 6;

    const auto T = 0.25;

    const auto dt = T / 4.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    auto solutionPair = crankNicolson(u0,
            dt,
            T,
            N,
            boundaries::zero,
            boundaries::zero,
            coefficient_functions::a2);

    const auto time = solutionPair.second;
    const auto u = solutionPair.first;

    ASSERT_EQ(5, u.cols());
    ASSERT_EQ(1, time.cols());

    ASSERT_EQ(5, time.rows());
    ASSERT_EQ(N + 2, u.rows());



}



TEST(TestCrankNicolson, TimeVectorCheck) {
    const int N = 8;

    const auto T = 0.25;

    const auto dt = T / 16.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    auto solutionPair = crankNicolson(u0,
            dt,
            T,
            N,
            boundaries::zero,
            boundaries::zero,
            coefficient_functions::a2);

    const auto time = solutionPair.second;
    const auto u = solutionPair.first;

    ASSERT_EQ(17, u.cols());
    ASSERT_EQ(1, time.cols());

    ASSERT_EQ(17, time.rows());
    ASSERT_EQ(N + 2, u.rows());



    for (int i = 0; i < time.rows(); ++i) {
        ASSERT_DOUBLE_EQ(dt * i, time(i))
                << "Wrong time at row " << i;
    }


}


TEST(TestCrankNicolson, ZeroInitialDataConstant) {
    const int N = 8;

    const auto T = 0.25;

    const auto dt = T / 8.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    u0.setZero();
    auto solutionPair = crankNicolson(u0,
            dt,
            T,
            N,
            boundaries::zero,
            boundaries::zero,
            coefficient_functions::a2);

    const auto time = solutionPair.second;
    const auto u = solutionPair.first;

    ASSERT_EQ(9, u.cols());
    ASSERT_EQ(1, time.cols());

    ASSERT_EQ(9, time.rows());
    ASSERT_EQ(N + 2, u.rows());



    for (int i = 0; i < u.rows(); ++i) {
        for (int j = 0; j < u.cols(); ++j) {
            ASSERT_DOUBLE_EQ(0.0, u(i, j))
                    << "Solution was given zero initial data, but is not zero at\n"
                        << "(" << i << ", " << j << ")";
        }
    }
}

TEST(TestCrankNicolson, FortyTwoInitialDataConstant) {
    const int N = 8;

    const auto T = 0.25;

    const auto dt = T / 8.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    u0.setOnes();
    u0 = 42.0 * u0;

    // make sure our boundary is 42
    auto g42 = [](double ) {
        return 42;
    };

    auto solutionPair = crankNicolson(u0,
            dt,
            T,
            N,
            g42,
            g42,
            coefficient_functions::a2);

    const auto time = solutionPair.second;
    const auto u = solutionPair.first;

    ASSERT_EQ(9, u.cols());
    ASSERT_EQ(1, time.cols());

    ASSERT_EQ(9, time.rows());
    ASSERT_EQ(N + 2, u.rows());



    for (int i = 0; i < u.rows(); ++i) {
        for (int j = 0; j < u.cols(); ++j) {
            ASSERT_DOUBLE_EQ(42.0, u(i, j))
                    << "Solution was given constant = 42 initial data, but is not zero at\n"
                        << "(" << i << ", " << j << ")";
        }
    }
}



TEST(TestCrankNicolson, FortyTwoInitialDataConstantZeroBoundary) {
    const int N = 8;

    const auto T = 0.25;

    const auto dt = T / 8.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    u0.setOnes();
    u0 = 42.0 * u0;
    u0[0] = 0;
    u0[N + 1] = 0;



    auto solutionPair = crankNicolson(u0,
            dt,
            T,
            N,
            boundaries::zero,
            boundaries::zero,
            coefficient_functions::a2);

    const auto time = solutionPair.second;
    const auto u = solutionPair.first;

    ASSERT_EQ(9, u.cols());
    ASSERT_EQ(1, time.cols());

    ASSERT_EQ(9, time.rows());
    ASSERT_EQ(N + 2, u.rows());




    for (int j = 0; j < u.cols(); ++j) {
        ASSERT_DOUBLE_EQ(0.0, u(0, j))
                << "Solution does not obey boundary at left side at column " << j;

        ASSERT_DOUBLE_EQ(0.0, u(N + 1, j))
                << "Solution does not obey boundary at right side at column " << j;
    }



}

TEST(TestCrankNicolson, FortyTwoInitialDataTimeVaryingBoundary) {
    const int N = 8;

    const auto T = 0.25;

    const auto dt = T / 8.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    u0.setOnes();
    u0 = 42.0 * u0;
    u0[0] = 0;
    u0[N + 1] = 0;


    auto gL = [](double t) {
        return std::sin(15.3 * M_PI * t);
    };

    auto gR = [](double t) {
        return std::cos(15.3 * M_PI * t);
    };

    auto solutionPair = crankNicolson(u0,
            dt,
            T,
            N,
            gL,
            gR,
            coefficient_functions::a2);

    const auto time = solutionPair.second;
    const auto u = solutionPair.first;

    ASSERT_EQ(9, u.cols());
    ASSERT_EQ(1, time.cols());

    ASSERT_EQ(9, time.rows());
    ASSERT_EQ(N + 2, u.rows());




    for (int j = 1; j < u.cols(); ++j) {
        ASSERT_DOUBLE_EQ(gL(time[j]), u(0, j))
                << "Solution does not obey boundary at left side at column " << j;

        ASSERT_DOUBLE_EQ(gR(time[j]), u(N + 1, j))
                << "Solution does not obey boundary at right side at column " << j;
    }



}





TEST(TestCrankNicolson, SquaredInitialDataZeroBoundary) {
    const int N = 3;

    const auto T = 0.25;

    const auto dt = T / 2.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N + 2, 0, 1);

    for (int i = 0; i < N + 2; i++) {
        u0[i] = (x[i] - 0.5) * (x[i] - 0.5) - 0.25;
    }



    auto solutionPair = crankNicolson(u0,
            dt,
            T,
            N,
            boundaries::zero,
            boundaries::zero,
            coefficient_functions::a2);

    const auto time = solutionPair.second;
    const auto u = solutionPair.first;

    ASSERT_EQ(3, u.cols());
    ASSERT_EQ(1, time.cols());

    ASSERT_EQ(3, time.rows());
    ASSERT_EQ(N + 2, u.rows());

    const std::string helpMessage =
        "Precomputed check failed\n. You can not always expect to pass this test, your solution can be correct even if this fails\n (but check that the values are close, if they are not, something is wrong).\nIf this test fails with the values being almost equal, try using the Eigen::SparseLU solver\n\nExpected solution:\n\t[ 0.     -0.1875 -0.25   -0.1875  0.    ]\n[ 0.         -0.04464286 -0.07142857 -0.04464286  0.        ]\n[ 0.         -0.01403061 -0.01530612 -0.01403061  0.        ]\n\nWhole solution:\n\t";
    // time index = 0
    ASSERT_NEAR(0.0000000000000000, u(0, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.1875000000000000, u(1, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.2500000000000000, u(2, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.1875000000000000, u(3, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0000000000000000, u(4, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;

    // time index = 1
    ASSERT_NEAR(0.0000000000000000, u(0, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.0446428571428571, u(1, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.0714285714285714, u(2, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.0446428571428571, u(3, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0000000000000000, u(4, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;

    // time index = 2
    ASSERT_NEAR(0.0000000000000000, u(0, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.0140306122448980, u(1, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.0153061224489796, u(2, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.0140306122448980, u(3, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0000000000000000, u(4, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;

}




TEST(TestCrankNicolson, SquaredInitialDataTimeVaryingBoundary) {
    const int N = 3;

    const auto T = 0.25;

    const auto dt = T / 2.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N + 2, 0, 1);

    for (int i = 0; i < N + 2; i++) {
        u0[i] = (x[i] - 0.5) * (x[i] - 0.5) - 0.25;
    }

    auto gL = [](double t) {
        return t;
    };

    auto gR = [](double t) {
        return t * t;
    };



    auto solutionPair = crankNicolson(u0,
            dt,
            T,
            N,
            gL,
            gR,
            coefficient_functions::a2);

    const auto time = solutionPair.second;
    const auto u = solutionPair.first;

    ASSERT_EQ(3, u.cols());
    ASSERT_EQ(1, time.cols());

    ASSERT_EQ(3, time.rows());
    ASSERT_EQ(N + 2, u.rows());


    const std::string helpMessage =
        "Precomputed check failed\n. You can not always expect to pass this test, your solution can be correct even if this fails\n (but check that the values are close, if they are not, something is wrong).\nIf this test fails with the values being almost equal, try using the Eigen::SparseLU solver\n\nExpected solution:\n\t[ 0.     -0.1875 -0.25   -0.1875  0.    ]\n[ 0.125       0.00372024 -0.05133929 -0.0327381   0.015625  ]\n[0.25       0.12790533 0.06377551 0.04109977 0.0625    ]\n\nWhole solution:\n\t";
    // time index = 0
    ASSERT_NEAR(0.0000000000000000, u(0, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.1875000000000000, u(1, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.2500000000000000, u(2, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.1875000000000000, u(3, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0000000000000000, u(4, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;

    // time index = 1
    ASSERT_NEAR(0.1250000000000000, u(0, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0037202380952381, u(1, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.0513392857142857, u(2, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.0327380952380952, u(3, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0156250000000000, u(4, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;

    // time index = 2
    ASSERT_NEAR(0.2500000000000000, u(0, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.1279053287981859, u(1, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0637755102040816, u(2, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0410997732426304, u(3, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0625000000000000, u(4, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;



}




TEST(TestCrankNicolson,
    SquaredInitialDataTimeVaryingBoundaryVaryingCoefficients) {
    const int N = 3;

    const auto T = 0.25;

    const auto dt = T / 2.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N + 2, 0, 1);

    for (int i = 0; i < N + 2; i++) {
        u0[i] = (x[i] - 0.5) * (x[i] - 0.5) - 0.25;
    }

    auto gL = [](double t) {
        return t;
    };

    auto gR = [](double t) {
        return t * t;
    };


    auto a = [] (double x) {
        return 1 + x;
    };

    auto solutionPair = crankNicolson(u0,
            dt,
            T,
            N,
            gL,
            gR,
            a);

    const auto time = solutionPair.second;
    const auto u = solutionPair.first;

    ASSERT_EQ(3, u.cols());
    ASSERT_EQ(1, time.cols());

    ASSERT_EQ(3, time.rows());
    ASSERT_EQ(N + 2, u.rows());

    const std::string helpMessage =
        "Precomputed check failed\n. You can not always expect to pass this test, your solution can be correct even if this fails\n (but check that the values are close, if they are not, something is wrong).\nIf this test fails with the values being almost equal, try using the Eigen::SparseLU solver\n\nExpected solution:\n\t[ 0.     -0.1875 -0.25   -0.1875  0.    ]\n[0.125      0.03702867 0.00368027 0.01445205 0.015625  ]\n[0.25       0.15465417 0.0987858  0.06220094 0.0625    ]\n\nWhole solution:\n\t";
    // time index = 0
    ASSERT_NEAR(0.0000000000000000, u(0, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.1875000000000000, u(1, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.2500000000000000, u(2, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(-0.1875000000000000, u(3, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0000000000000000, u(4, 0), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;

    // time index = 1
    ASSERT_NEAR(0.1250000000000000, u(0, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0370286673553719, u(1, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0036802685950413, u(2, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0144520488980716, u(3, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0156250000000000, u(4, 1), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;

    // time index = 2
    ASSERT_NEAR(0.2500000000000000, u(0, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.1546541675204335, u(1, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0987858012886187, u(2, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0622009444558280, u(3, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;
    ASSERT_NEAR(0.0625000000000000, u(4, 2), CRANK_NICOLSON_ABSOLUTE_TOLERANCE)
            << helpMessage << u;


}





