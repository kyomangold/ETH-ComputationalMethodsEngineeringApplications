#include "create_poisson_matrix.hpp"
#include "boundaries.hpp"
#include "coefficient_functions.hpp"
#include <gtest/gtest.h>

TEST(TestCreatePoissonMatrix, SizeCheck5x5) {
    const auto N = 5;
    const auto a = [](double) {
        return 1.0;
    };

    const auto A = createPoissonMatrix(N, a);

    ASSERT_EQ(N, A.cols())
            << "Poisson Matrix does not have correct number of column. N = " << N
                << ", A.cols() = " << A.cols();
    ASSERT_EQ(N, A.rows())
            << "Poisson Matrix does not have correct number of rows. N = " << N
                << ", A.rows() = " << A.rows();


}

TEST(TestCreatePoissonMatrix, SizeCheck42x42) {
    const auto N = 42;
    const auto a = [](double) {
        return 1.0;
    };

    const auto A = createPoissonMatrix(N, a);

    ASSERT_EQ(N, A.cols())
            << "Poisson Matrix does not have correct number of column. N = " << N
                << ", A.cols() = " << A.cols();
    ASSERT_EQ(N, A.rows())
            << "Poisson Matrix does not have correct number of rows. N = " << N
                << ", A.rows() = " << A.rows();


}

TEST(TestCreatePoissonMatrix, ElementCheckConstantCoefficient) {
    const auto N = 4;
    const auto a = [](double) {
        return 1.0;
    };

    // we will multiply A with h^2 to make the coefficient nicer
    const auto h = 1. / (N + 1);
    const SparseMatrix A = h * h * createPoissonMatrix(N, a);


    // first check the size
    ASSERT_EQ(N, A.cols())
            << "Poisson Matrix does not have correct number of column. N = " << N
                << ", A.cols() = " << A.cols();
    ASSERT_EQ(N, A.rows())
            << "Poisson Matrix does not have correct number of rows. N = " << N
                << ", A.rows() = " << A.rows();

    // Now A should be
    //  2 -1  0  0
    // -1  2 -1  0
    //  0 -1  2 -1
    //  0  0 -1  2
    // then we manually check each element

    // first row
    ASSERT_DOUBLE_EQ(2., A.coeff(0, 0));
    ASSERT_DOUBLE_EQ(-1., A.coeff(0, 1));
    ASSERT_DOUBLE_EQ(0., A.coeff(0, 2));
    ASSERT_DOUBLE_EQ(0., A.coeff(0, 3));

    // second row
    ASSERT_DOUBLE_EQ(-1., A.coeff(1, 0));
    ASSERT_DOUBLE_EQ(2., A.coeff(1, 1));
    ASSERT_DOUBLE_EQ(-1., A.coeff(1, 2));
    ASSERT_DOUBLE_EQ(0., A.coeff(1, 3));

    // third row
    ASSERT_DOUBLE_EQ(0., A.coeff(2, 0));
    ASSERT_DOUBLE_EQ(-1., A.coeff(2, 1));
    ASSERT_DOUBLE_EQ(2., A.coeff(2, 2));
    ASSERT_DOUBLE_EQ(-1., A.coeff(2, 3));

    // fourth row
    ASSERT_DOUBLE_EQ(0., A.coeff(3, 0));
    ASSERT_DOUBLE_EQ(0., A.coeff(3, 1));
    ASSERT_DOUBLE_EQ(-1., A.coeff(3, 2));
    ASSERT_DOUBLE_EQ(2., A.coeff(3, 3));

}

TEST(TestCreatePoissonMatrix, ElementCheckCubedCoefficient) {
    const auto N = 4;
    const auto a = [](double x) {
        return 1.0 + x * x;
    };

    // we will multiply A with h^2 to make the coefficient nicer
    const auto h = 1. / (N + 1);
    SparseMatrix A = h * h * createPoissonMatrix(N, a);


    // first check the size
    ASSERT_EQ(N, A.cols())
            << "Poisson Matrix does not have correct number of column. N = " << N
                << ", A.cols() = " << A.cols();
    ASSERT_EQ(N, A.rows())
            << "Poisson Matrix does not have correct number of rows. N = " << N
                << ", A.rows() = " << A.rows();

    // first row
    ASSERT_DOUBLE_EQ(a(h * 1) * 2., A.coeff(0, 0));
    ASSERT_DOUBLE_EQ(a(h * 1) * -1., A.coeff(0, 1));
    ASSERT_DOUBLE_EQ(a(h * 1) * 0., A.coeff(0, 2));
    ASSERT_DOUBLE_EQ(a(h * 1) * 0., A.coeff(0, 3));

    // second row
    ASSERT_DOUBLE_EQ(a(h * 2) * - 1., A.coeff(1, 0));
    ASSERT_DOUBLE_EQ(a(h * 2) * 2., A.coeff(1, 1));
    ASSERT_DOUBLE_EQ(a(h * 2) * - 1., A.coeff(1, 2));
    ASSERT_DOUBLE_EQ(a(h * 2) * 0., A.coeff(1, 3));

    // third row
    ASSERT_DOUBLE_EQ(a(h * 3) * 0., A.coeff(2, 0));
    ASSERT_DOUBLE_EQ(a(h * 3) * -1., A.coeff(2, 1));
    ASSERT_DOUBLE_EQ(a(h * 3) * 2., A.coeff(2, 2));
    ASSERT_DOUBLE_EQ(a(h * 3) * -1., A.coeff(2, 3));

    // fourth row
    ASSERT_DOUBLE_EQ(a(h * 4) * 0., A.coeff(3, 0));
    ASSERT_DOUBLE_EQ(a(h * 4) * 0., A.coeff(3, 1));
    ASSERT_DOUBLE_EQ(a(h * 4) * -1., A.coeff(3, 2));
    ASSERT_DOUBLE_EQ(a(h * 4) * 2., A.coeff(3, 3));

}


