#include <Eigen/Sparse>
#include <iostream>
#include "writer.hpp"
#include <cmath>
#include <Eigen/SparseCholesky>
#include <stdexcept>

//! Sparse Matrix type. Makes using this type easier.
typedef Eigen::SparseMatrix<double> SparseMatrix;

//! Used for filling the sparse matrix.
typedef Eigen::Triplet<double> Triplet;

//! Vector type
typedef Eigen::VectorXd Vector;


typedef double(*FunctionPointer)(double);

//----------------poissonBegin----------------
//! Create the 1D Poisson matrix
//! @param[out] A will contain the Poisson matrix
//! @param[in] N the number of points
void createPoissonMatrix(SparseMatrix& A, int N) {
    A.resize(N, N);
    std::vector<Triplet> triplets;
    triplets.reserve(N + 2 * N - 2);

    for (int i = 1; i < N-1; ++i) {
    // This is the diagonal
        triplets.push_back(Triplet(i, i, 2));

        //// CMEA_START_TEMPLATE
        triplets.push_back(Triplet(i, i - 1, -1));
        triplets.push_back(Triplet(i, i + 1, -1));
        //// CMEA_END_TEMPLATE
    }
    // Deal with first and last row here
    //// CMEA_START_TEMPLATE
    // first row
    triplets.push_back(Triplet(0,0,2));
    triplets.push_back(Triplet(0,1,-1));
    // last row
    triplets.push_back(Triplet(N-1,N-1,2));
    triplets.push_back(Triplet(N-1,N-2,-1));
    //// CMEA_END_TEMPLATE


    A.setFromTriplets(triplets.begin(), triplets.end());
}
//----------------poissonEnd----------------


//----------------RHSBegin----------------
//! Create the right hand side for the poisson problem.
//! @note This scales the right hand side (ie. $dx^2 * f(x)$)
//!
//! @param[out] rhs will contain the right hand side
//! @param[in] f function pointer to f
//! @param[in] N the number of points to use
//! @param[in] dx the cell length
void createRHS(Vector& rhs, FunctionPointer f, int N, double dx) {
    rhs.resize(N);
    // Set RHS
    //// CMEA_START_TEMPLATE
    for (int i = 0; i < N; ++i) {
        const double x = (i + 1) * dx;
        rhs[i] = dx * dx * f(x);
    }
    //// CMEA_END_TEMPLATE
}
//----------------RHSEnd----------------

//----------------PoissonSolveBegin----------------
//! Solves the Poisson equation
//!
//!   $-u''(x) = f(x) $
//!
//! on [0,1] with boundary values $u(0)=u(1) = 0$.
//!
//! @param[out] u should contain the solution u at the end
//! @param[in] f should be a function pointer to f
//! @param[in] N as in the exercise
void poissonSolve(Vector& u, FunctionPointer f, int N) {
    double dx = 1.0 / (N + 1);

    SparseMatrix A;
    // create the matrix
    //// CMEA_START_TEMPLATE
    createPoissonMatrix(A, N);
    //// CMEA_END_TEMPLATE

    Vector rhs;
    // create RHS
    //// CMEA_START_TEMPLATE
    createRHS(rhs, f, N, dx);
    //// CMEA_END_TEMPLATE

    Eigen::SparseLU<SparseMatrix> solver;

    solver.compute(A);

    if ( solver.info() !=  Eigen::Success) {
        throw std::runtime_error("Could not decompose the matrix");
    }

    // Find u: ....
    //// CMEA_START_TEMPLATE
	u.resize(N);
	u = solver.solve(rhs);
    //// CMEA_END_TEMPLATE
}
//----------------PoissonSolveEnd----------------


//! Test if the Poisson matrix is correctly set up for one case
//! This does NOT guarantee that the code is correct, it is only a small
//! indication
void testPoissonMatrix() {
    SparseMatrix A;
    const int N = 13;
    createPoissonMatrix(A, N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; j++) {

            if (i==j) {
                if ( A.coeff(i,j) != 2) {
                    throw std::runtime_error("Poisson matrix: Wrong Poisson matrix");
                }
            } else if (i == j - 1 || i == j + 1) {
                if (A.coeff(i,j) != -1) {
                    throw std::runtime_error("Poisson matrix: Wrong upper or lower diagonal");
                }
            } else {
                if (A.coeff(i,j) != 0) {
                    throw std::runtime_error("Poisson matrix: Matrix is not band diagonal");
                }
            }
        }
    }
}


double F(double x) {
    return sin(2 * M_PI * x);
}

double exact(double x) {
    return 1.0/ (4 * M_PI * M_PI) * F(x);
}

//----------------convergenceBegin----------------
//! Computes error for a range of cell lengths and stores them to errors
//! @param[out] errors the errors computed
//! @param[out] resolutions the resolutions used
void poissonConvergence( std::vector<double>& errors,
                         std::vector<double>& resolutions) {

    const int startK = 4;
    const int endK = 10;
    errors.resize(endK - startK + 1);
    resolutions.resize(errors.size());
    for (int k = startK; k <= endK; ++k) {
        const int N = (1 << k) - 1;
	    // compute the solution and the error
	    //// CMEA_START_TEMPLATE
        Vector u;
        poissonSolve(u, F, N);

        double maxError = 0;
        const double dx = 1.0 / (N + 1);
        for (int i = 0; i < N; ++i) {
            const double x = (i+1) * dx;
            const double error = std::abs(u[i] - exact(x));
            maxError = std::max(error, maxError);

        }

        errors[k - startK] = maxError;
        resolutions[k - startK] = N;
	    //// CMEA_END_TEMPLATE
    }

}
//----------------convergenceEnd----------------

int main(int argc, char** argv) {

    try {
        testPoissonMatrix();
    } catch (const std::exception& e) {
	    std::cout << "Warning: " << e.what() << " (or not implemented)" << std::endl;
        return -1;
    }

    Vector u;

    poissonSolve(u, F, 200);
    writeToFile("u_fd.txt", u);

    std::vector<double> errors;
    std::vector<double> resolutions;
    poissonConvergence(errors, resolutions);
    writeToFile("errors_fd.txt", errors);
    writeToFile("resolutions_fd.txt", resolutions);
}
