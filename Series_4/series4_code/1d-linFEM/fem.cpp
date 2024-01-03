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
typedef double(*FunctionPointer)(double);

//----------------StiffnessMatrixBegin----------------
//! Create the Stiffness matrix for 1D with boundary terms
//! @param[out] A will be the Galerkin matrix (as in the exercise)
//! @param[in] N as in the exercise
//! @param[in] dx the cell length

void createStiffnessMatrix(SparseMatrix& A, int N, double dx) {
    std::vector<Triplet> triplets;
    A.resize(N, N);

    // Reserve the space we will allocate
    triplets.reserve(N + 2 * N - 2);

    //// CMEA_START_TEMPLATE

    //first row
    triplets.push_back(Triplet(0,0,2./dx));
    triplets.push_back(Triplet(0,1,-1./dx));

    for (int i = 1; i < N-1; ++i) {
        triplets.push_back(Triplet(i, i, 2./dx)); // diagonal
        triplets.push_back(Triplet(i, i-1, -1./dx)); // subdiagonal
        triplets.push_back(Triplet(i, i+1, -1./dx)); // superdiagonal
    }

    // last row
    triplets.push_back(Triplet(N-1,N-1, 2./dx));
    triplets.push_back(Triplet(N-1,N-2, -1./dx));

    //// CMEA_END_TEMPLATE
    A.setFromTriplets(triplets.begin(), triplets.end());
}
//----------------StiffnessMatrixEnd----------------

//----------------RHSBegin----------------
//! Creates the right hand side for the finite element method in the exercise
//! @param[out] rhs the resulting right hand side (or phi in the exercise)
//! @param[in] f the function pointer to f
//! @param[in] dx the cell length
//! @param[in] x the x points

void createRHS(Vector& rhs, FunctionPointer f, int N, double dx, const Vector& x) {
    rhs.resize(N);

    //// CMEA_START_TEMPLATE
    for (int i = 0; i < N; ++i) {
        rhs[i] = dx * f(x[i]);
    }
    //// CMEA_END_TEMPLATE
}
//----------------RHSEnd----------------


//----------------FEMsolveBegin----------------
//! Solve the equation
//!
//!    $-u''(x) = f(x)$
//!
//! using FEM on the interval $[a,b]$.
//!
//! @param[out] u at the end, will have all the values of u
//! @param[out] x will be the x points
//! @param[in] f should be a pointer to the function f
//! @param[in] N number of inner cells to use
//! @param[in] a endpoint on left side
//! @param[in] b endpoint on right side
//! @param[in] ua boundary value at a
//! @param[in] ub boundary value at b

void femSolve(Vector& u, Vector& x, FunctionPointer f, int N, double a, double b,
              double ua = 0.0, double ub = 0.0) {

    double dx = (b - a) / (N + 1);

    // Fill x vector
    x.setLinSpaced(N+2, a, b);

    SparseMatrix A;
    createStiffnessMatrix(A, N, dx);

    Vector phi;
    createRHS(phi, f, N, dx, x.segment(1,N));

    //// CMEA_START_TEMPLATE
    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(A);

    if ( solver.info() !=  Eigen::Success) {
        throw std::runtime_error("Could not decompose the matrix");
    }
    //// CMEA_END_TEMPLATE

    u.resize((N + 2));
    u.setZero();

    // Do boundary conditions
    //// CMEA_START_TEMPLATE
    // offset function approach:
    // g linear, such that w = u - g verifies -w'' = f, w(a)=w(b)=0
    Vector g(N+2);
    Vector e = Vector::Ones(N+2);
    g = (x-a*e)*(ub-ua)/(b-a) + ua*e;
    //// CMEA_END_TEMPLATE

    //solve inner system
    //// CMEA_START_TEMPLATE
    u.segment(1, N) = solver.solve(phi);
    u += g;
    //// CMEA_END_TEMPLATE

}
//----------------FEMsolveEnd----------------

void computeWithoutBoundaryValues() {
    Vector u;
    Vector x;
    femSolve(u, x, std::sin, 100, -M_PI, M_PI);

    writeToFile("u_wo_bc_fem.txt", u);
    writeToFile("x_wo_bc_fem.txt", x);
}

void computeWithBoundaryValues() {
    Vector u;
    Vector x;
    femSolve(u, x, std::cos, 100, -M_PI, M_PI, -1.0, 1.0/2.0);

    writeToFile("u_w_bc_fem.txt", u);
    writeToFile("x_w_bc_fem.txt", x);
}

int main(int, char**) {
    computeWithoutBoundaryValues();
    computeWithBoundaryValues();
}




