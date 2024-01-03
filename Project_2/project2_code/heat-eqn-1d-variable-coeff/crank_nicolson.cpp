#include "crank_nicolson.hpp"


//! Uses the Crank-Nicolson method to compute u from time 0 to time T
//!
//! @param[in] u0 the initial data, as column vector (size N+2)
//! @param[in] dt the time step size
//! @param[in] T the final time at which to compute the solution (which we assume to be a multiple of dt)
//! @param[in] N the number of interior grid points
//! @param[in] gL function of time with the Dirichlet condition at left boundary
//! @param[in] gR function of time with the Dirichlet condition at right boundary
//! @param[in] a the coefficient function a
//!
//! @note the vector returned should include the boundary values!
//!
std::pair<Eigen::MatrixXd, Eigen::VectorXd> crankNicolson(
    const Eigen::VectorXd& u0,
    double dt, double T, int N,
    const std::function<double(double)>& gL,
    const std::function<double(double)>& gR,
    const std::function<double(double)>& a) {

    Eigen::VectorXd time;
    Eigen::MatrixXd u;


    //// CMEA_START_TEMPLATE
    const int nsteps = int(round(T / dt));
    const double h = 1. / (N + 1);
    u.resize(N + 2, nsteps + 1);
    time.resize(nsteps + 1);
    /* Initialize A */

    auto A = createPoissonMatrix(N, a);

    SparseMatrix B(N, N);
    B.setIdentity();
    B += dt / 2.*A;


    /* Initialize u */
    u.col(0) << u0;
    // initialize time
    time[0] = 0.;


    /* Initialize solver and compute LU decomposition of B (Note: since dt is constant, the matrix B is the same for all timesteps)*/
    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(B);

    Eigen::VectorXd G1 = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd G2 = Eigen::VectorXd::Zero(N);

    for (int k = 0; k < nsteps; k++) {
        time[k + 1] = (k + 1) * dt;
        G1[0] = a(h) * dt * gL(time[k]) / (h * h);
        G1[N - 1] = a(1 - h) * dt * gR(time[k]) / (h * h);

        G2[0] = a(h) * dt * gL(time[k + 1]) / (h * h);
        G2[N - 1] = a(1 - h) *  dt * gR(time[k + 1]) / (h * h);

        const Eigen::VectorXd rhs =
            u.col(k).segment(1, N) - dt / 2 * A * u.col(k).segment(1, N)
            + 0.5 * (G1 + G2);

        u.col(k + 1).segment(1, N) = solver.solve(rhs);

        u.col(k + 1)[0] = gL(time[k + 1]);

        u.col(k + 1)[N + 1] = gR(time[k + 1]);
    }

    //// CMEA_END_TEMPLATE

    return std::make_pair(u, time);
}
