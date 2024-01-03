#include <Eigen/Sparse>
#include <iostream>
#include "writer.hpp"
#include <cmath>
#include <Eigen/SparseCholesky>
#include <stdexcept>
#include <functional>

//! Sparse Matrix type. Makes using this type easier.
typedef Eigen::SparseMatrix<double> SparseMatrix;

//! Used for filling the sparse matrix.
typedef Eigen::Triplet<double> Triplet;

//! Vector type
typedef Eigen::VectorXd Vector;

//! Create the 1D Poisson matrix
//! @param[out] A will contain the Poisson matrix
//! @param[in] N the number of interior points
void createPoissonMatrix(SparseMatrix& A, int N) {
  //// CMEA_START_TEMPLATE
    A.resize(N, N);
    double h=1./(N+1);
    std::vector<Triplet> triplets;
    triplets.reserve(N + 2 * N - 2);
    for (int i = 0; i < N; ++i) {
        triplets.push_back(Triplet(i, i, 2./(h*h)));
        if (i > 0) {
            triplets.push_back(Triplet(i, i - 1, -1./(h*h)));
        }
        if (i < N - 1){
            triplets.push_back(Triplet(i, i + 1, -1./(h*h)));
        }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    //// CMEA_END_TEMPLATE
}

/// Uses the explicit Euler method to compute u from time 0 to time T 
///
/// @param[out] u at all time steps up to time T, each column corresponding to a time step (including the initial condition as first column)
/// @param[out] time the time levels
/// @param[in] u0 the initial data, as column vector
/// @param[in] dt the time step size
/// @param[in] T the final time at which to compute the solution (which we assume to be a multiple of dt)
/// @param[in] N the number of interior grid points
/// @param[in] gL function of time with the Dirichlet condition at left boundary
/// @param[in] gR function of time with the Dirichlet condition at right boundary
///

void explicitEuler(Eigen::MatrixXd & u, Vector & time,
                   const Vector u0, double dt, double T, int N,
                   const std::function<double(double)>& gL,
                   const std::function<double(double)>& gR) {
    const unsigned int nsteps = round(T/dt);
    const double h=1./(N+1);
    u.resize(N,nsteps+1);
    time.resize(nsteps+1);

    //// CMEA_START_TEMPLATE
    /* Initialize A */
    SparseMatrix A;
    createPoissonMatrix(A,N);
    /* Initialize u */
    u.col(0)<<u0;
    time[0]=0.;
    Vector G;
    G.resize(N);
    G.setZero();
    for(unsigned k=0; k<nsteps; k++)
    {
        G[0] = dt*gL(time[k])/(h*h);
        G[N-1] = dt*gR(time[k])/(h*h);

        u.col(k+1)=u.col(k)-dt*A*u.col(k)+G;
        time[k+1]=(k+1)*dt;
    }
    //// CMEA_END_TEMPLATE
}

/// Uses the Crank-Nicolson method to compute u from time 0 to time T 
///
/// @param[out] u at all time steps up to time T, each column corresponding to a time step (including the initial condition as first column)
/// @param[out] time the time levels
/// @param[in] u0 the initial data, as column vector
/// @param[in] dt the time step size
/// @param[in] T the final time at which to compute the solution (which we assume to be a multiple of dt)
/// @param[in] N the number of interior grid points
/// @param[in] gL function of time with the Dirichlet condition at left boundary
/// @param[in] gR function of time with the Dirichlet condition at right boundary
///

void CrankNicolson(Eigen::MatrixXd & u, Vector & time,
                   const Vector u0, double dt, double T, int N,
                   const std::function<double(double)>& gL,
                   const std::function<double(double)>& gR) {
  //// CMEA_START_TEMPLATE
    const unsigned int nsteps = round(T/dt);
    const double h=1./(N+1);
    u.resize(N,nsteps+1);
    time.resize(nsteps+1);
    /* Initialize A */
    SparseMatrix A;
    createPoissonMatrix(A,N);
    SparseMatrix B(N,N);
    B.setIdentity();
    B+=dt/2.*A;
    /* Initialize u */
    u.col(0)<<u0;
    time[0]=0.;
    /* Initialize solver and compute Cholesky decomposition of B (Note: since dt is constant, the matrix B is the same for all timesteps)*/
    Eigen::SimplicialLDLT<SparseMatrix> solver;
    solver.compute(B);
    Vector G1,G2;
    G1.resize(N);
    G1.setZero();
    G2.resize(N);
    G2.setZero();
    for(unsigned k=0; k<nsteps; k++)
    {
        time[k+1]=(k+1)*dt;
        G1[0] = dt*gL(time[k])/(h*h);
        G1[N-1] = dt*gR(time[k])/(h*h);

        G2[0] = dt*gL(time[k+1])/(h*h);
        G2[N-1] = dt*gR(time[k+1])/(h*h);

        u.col(k+1)=solver.solve(u.col(k)-dt/2*A*u.col(k)
				+0.5*(G1+G2));

    }

    //// CMEA_END_TEMPLATE
}

double U0(double x) {
    return 1+std::min(2*x,2-2*x);
}

int main(int, char**) {
    double T = 0.3;
    double dt = 0.0002; // Change this for explicit / implicit time stepping comparison
    int N=40;
    Vector u0(N);
    double h=1./(N+1);
    /* Initialize u0 */
    for(int i=0;i<u0.size();i++)
        u0[i] = U0(h*(i+1));
    auto gR = [](double t) { return std::exp(-10*t); };
    auto gL = [](double t) { return std::exp(-10*t); };
    Eigen::MatrixXd u_euler;
    Vector time;
    explicitEuler(u_euler,time,u0,dt,T,N, gL, gR );
    writeToFile("time.txt",time);
    writeMatrixToFile("u_explEuler.txt", u_euler);
    Eigen::MatrixXd u_cn;
    CrankNicolson(u_cn,time,u0,dt,T,N, gL, gR);
    writeMatrixToFile("u_cn.txt",u_cn);
}




