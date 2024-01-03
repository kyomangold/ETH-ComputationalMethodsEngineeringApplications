#include <gtest/gtest.h>

#include "forward_euler.hpp"
#include "boundaries.hpp"
#include "coefficient_functions.hpp"

TEST(TestForwardEuler, SizeCheck6x5) {
    const int N = 6;

    const auto T = 0.25;

    const auto dt = T / 4.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    auto solutionPair = forwardEuler(u0,
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



TEST(TestForwardEuler, TimeVectorCheck) {
    const int N = 8;

    const auto T = 0.25;

    const auto dt = T / 16.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    auto solutionPair = forwardEuler(u0,
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


TEST(TestForwardEuler, ZeroInitialDataConstant) {
    const int N = 8;

    const auto T = 0.25;

    const auto dt = T / 8.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    u0.setZero();
    auto solutionPair = forwardEuler(u0,
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

TEST(TestForwardEuler, FortyTwoInitialDataConstant) {
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

    auto solutionPair = forwardEuler(u0,
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


TEST(TestForwardEuler, FortyTwoInitialDataConstantZeroBoundary) {
    const int N = 8;

    const auto T = 0.25;

    const auto dt = T / 8.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    u0.setOnes();
    u0 = 42.0 * u0;
    u0[0] = 0;
    u0[N + 1] = 0;



    auto solutionPair = forwardEuler(u0,
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

TEST(TestForwardEuler, FortyTwoInitialDataTimeVaryingBoundary) {
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

    auto solutionPair = forwardEuler(u0,
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


TEST(TestForwardEuler, SquaredInitialDataZeroBoundary) {
    const int N = 3;

    const auto T = 0.25;

    const auto dt = T / 2.0;


    Eigen::VectorXd u0;
    u0.resize(N + 2);
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N + 2, 0, 1);

    for (int i = 0; i < N + 2; i++) {
        u0[i] = (x[i] - 0.5) * (x[i] - 0.5) - 0.25;
    }



    auto solutionPair = forwardEuler(u0,
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



    // time index = 0
    ASSERT_DOUBLE_EQ(0.0, u(0, 0));
    ASSERT_DOUBLE_EQ(-0.1875, u(1, 0));
    ASSERT_DOUBLE_EQ(-0.25, u(2, 0));
    ASSERT_DOUBLE_EQ(-0.1875, u(3, 0));
    ASSERT_DOUBLE_EQ(0.0, u(4, 0));

    // time index = 1
    ASSERT_DOUBLE_EQ(0.0, u(0, 1));
    ASSERT_DOUBLE_EQ(0.0625, u(1, 1));
    ASSERT_DOUBLE_EQ(0.0, u(2, 1));
    ASSERT_DOUBLE_EQ(0.0625, u(3, 1));
    ASSERT_DOUBLE_EQ(0.0, u(4, 1));

    // time index = 2
    ASSERT_DOUBLE_EQ(0.0, u(0, 2));
    ASSERT_DOUBLE_EQ(-0.1875, u(1, 2));
    ASSERT_DOUBLE_EQ(0.25, u(2, 2));
    ASSERT_DOUBLE_EQ(-0.1875, u(3, 2));
    ASSERT_DOUBLE_EQ(0.0, u(4, 2));
}



TEST(TestForwardEuler, SquaredInitialDataTimeVaryingBoundary) {
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



    auto solutionPair = forwardEuler(u0,
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


    // time index = 0
    ASSERT_DOUBLE_EQ(0.0, u(0, 0));
    ASSERT_DOUBLE_EQ(-0.1875, u(1, 0));
    ASSERT_DOUBLE_EQ(-0.25, u(2, 0));
    ASSERT_DOUBLE_EQ(-0.1875, u(3, 0));
    ASSERT_DOUBLE_EQ(0.0, u(4, 0));

    // time index = 1
    ASSERT_DOUBLE_EQ(0.125, u(0, 1));
    ASSERT_DOUBLE_EQ(0.0625, u(1, 1));
    ASSERT_DOUBLE_EQ(0.0, u(2, 1));
    ASSERT_DOUBLE_EQ(0.0625, u(3, 1));
    ASSERT_DOUBLE_EQ(0.015625, u(4, 1));

    // time index = 2
    ASSERT_DOUBLE_EQ(0.25, u(0, 2));
    ASSERT_DOUBLE_EQ(0.0625, u(1, 2));
    ASSERT_DOUBLE_EQ(0.25, u(2, 2));
    ASSERT_DOUBLE_EQ(-0.15625, u(3, 2));
    ASSERT_DOUBLE_EQ(0.0625, u(4, 2));
}


TEST(TestForwardEuler,
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

    auto solutionPair = forwardEuler(u0,
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

    // time index = 0
    ASSERT_DOUBLE_EQ(0.0, u(0, 0));
    ASSERT_DOUBLE_EQ(-0.1875, u(1, 0));
    ASSERT_DOUBLE_EQ(-0.25, u(2, 0));
    ASSERT_DOUBLE_EQ(-0.1875, u(3, 0));
    ASSERT_DOUBLE_EQ(0.0, u(4, 0));

    // time index = 1
    ASSERT_DOUBLE_EQ(0.125, u(0, 1));
    ASSERT_DOUBLE_EQ(0.125, u(1, 1));
    ASSERT_DOUBLE_EQ(0.125, u(2, 1));
    ASSERT_DOUBLE_EQ(0.25, u(3, 1));
    ASSERT_DOUBLE_EQ(0.015625, u(4, 1));

    // time index = 2
    ASSERT_DOUBLE_EQ(0.25, u(0, 2));
    ASSERT_DOUBLE_EQ(0.125, u(1, 2));
    ASSERT_DOUBLE_EQ(0.5, u(2, 2));
    ASSERT_DOUBLE_EQ(-1.0078125, u(3, 2));
    ASSERT_DOUBLE_EQ(0.0625, u(4, 2));


}
