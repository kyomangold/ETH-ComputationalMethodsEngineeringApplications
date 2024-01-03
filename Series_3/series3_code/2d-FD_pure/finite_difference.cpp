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


//! Our function pointer, typedef'd to make it easier to use
typedef double(*FunctionPointer)(double, double);

//----------------poissonBegin----------------
//! Create the Poisson matrix for 2D finite difference.
//! @param[out] A will be the Poisson matrix (as in the exercise)
//! @param[in] N number of elements in the x-direction
void createPoissonMatrix2D(SparseMatrix& A, int N) {
     // Fill the matrix A using setFromTriplets - method (see other exercise
    // for how to use it).
    A.resize(N*N, N*N);
    //// CMEA_START_TEMPLATE
    std::vector<Triplet> triplets;
    triplets.reserve(5*N*N-4*N);
    for (int i = 0; i < N*N; ++i) {
        triplets.push_back(Triplet(i, i, 4));
        if (i % N != 0) {
            triplets.push_back(Triplet(i, i-1, -1));
        }
        if (i % N != N - 1) {
            triplets.push_back(Triplet(i, i+1, -1));
        }
        if (i >= N) {
            triplets.push_back(Triplet(i, i-N, -1));
        }
        if (i < N*N - N ) {
            triplets.push_back(Triplet(i, i+N, -1));
        }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    //// CMEA_END_TEMPLATE
}
//----------------poissonEnd----------------


//----------------RHSBegin----------------
//! Create the Right hand side for the 2D finite difference
//! @param[out] rhs will at the end contain the right hand side
//! @param[in] f the right hand side function f
//! @param[in] N the number of points in the x direction
//! @param[in] dx the cell width
void createRHS(Vector& rhs, FunctionPointer f, int N, double dx) {
    rhs.resize(N * N);
    // fill up RHS
    // remember that the index (i,j) corresponds to i*N+j
    //// CMEA_START_TEMPLATE
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            const double x = (i + 1) * dx;
            const double y = (j + 1) * dx;
            rhs[j * N + i] = dx * dx * f(x, y);
        }
    }
    //// CMEA_END_TEMPLATE
}
//----------------RHSEnd----------------


//----------------solveBegin----------------
//! Solve the Poisson equation in 2D
//! @param[out] u will contain the solution u
//! @param[in] f the function pointer to f
//! @param[in] N the number of points to use (in x direction)
void poissonSolve(Vector& u, FunctionPointer f, int N) {
    // Solve Poisson 2D here
    //// CMEA_START_TEMPLATE
    double dx = 1.0 / (N + 1);

    SparseMatrix A;
    createPoissonMatrix2D(A, N);

    Vector rhs;
    createRHS(rhs, f, N, dx);

    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(A);

    if ( solver.info() !=  Eigen::Success) {
        throw std::runtime_error("Could not decompose the matrix");
    }
    u.resize((N + 2) * (N + 2));
    u.setZero();

    Vector innerU = solver.solve(rhs);

    // Copy vector to inner u.
    for (int i = 1; i < N + 1; ++i) {
        for (int j = 1; j < N + 1; ++j) {
            u[i * (N + 2) + j] = innerU[(i - 1) * N + j - 1];
        }
    }

    //// CMEA_END_TEMPLATE
}
//----------------solveEnd----------------

//! Test if the Poisson matrix is correctly set up for one case
//! This does NOT guarantee that the code is correct, it is only a small
//! indication
void testPoissonMatrix() {
    SparseMatrix A;
    const int N = 13;
    createPoissonMatrix2D(A, N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; j++) {

            if (i==j) {
                if ( A.coeff(i,j) != 4) {
                    throw std::runtime_error("Poisson matrix: Wrong Poisson matrix");
                }
            } else if (i == j - 1 || i == j + 1 || i == j+N || i == j-N) {
                if (A.coeff(i,j) != -1) {
                    throw std::runtime_error("Poisson matrix: Wrong upper or lower diagonals");
                }
            } else {
                if (A.coeff(i,j) != 0) {
                    throw std::runtime_error("Poisson matrix: Matrix is not band diagonal");
                }
            }
        }
    }
}

double F(double x, double y) {
    return 8*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y);
}



int main(int, char**) {
    Vector u;

    try {
        testPoissonMatrix();
    } catch (const std::exception& e) {
        std::cout << "Warning: " << e.what() << " (or not implemented)" << std::endl;
        return -1;
    }

    poissonSolve(u, F, 100);
    writeToFile("u_fd.txt", u);
}




