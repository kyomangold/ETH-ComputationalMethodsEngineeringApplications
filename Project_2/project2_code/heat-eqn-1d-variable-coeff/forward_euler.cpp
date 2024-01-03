#include "forward_euler.hpp"
#include "create_poisson_matrix.hpp"

//! Uses the explicit forward Euler method to compute u from time 0 to time T
//!
//! @param[in] u0 the initial data, as column vector (size N+2)
//! @param[in] dt the time step size
//! @param[in] T the final time at which to compute the solution (which we assume to be a multiple of dt)
//! @param[in] N the number of interior grid points
//! @param[in] gL function of time with the Dirichlet condition at left boundary
//! @param[in] gR function of time with the Dirichlet condition at right boundary
//! @param[in] a the coefficient function a
//!
//! @return u at all time steps up to time T, each column corresponding to a time step (including the initial condition as first column)
//!
//! @note the vector returned should include the boundary values!
std::pair<Eigen::MatrixXd, Eigen::VectorXd> forwardEuler(
    const Eigen::VectorXd& u0,
    double dt,
    double T,
    int N,
    const std::function<double(double)>& gL,
    const std::function<double(double)>& gR,
    const std::function<double(double)>& a) {


    const int nsteps = int(round(T / dt));

    const double h = 1. / (N + 1);

    Eigen::MatrixXd u;
    u.resize(N + 2, nsteps + 1);

    Eigen::VectorXd time;
    time.resize(nsteps + 1);

    //// CMEA_START_TEMPLATE
    // Initialize A
    SparseMatrix A = createPoissonMatrix(N, a);
    // Initialize u
    u.col(0) << u0;
    time[0] = 0.;

    Eigen::VectorXd G = Eigen::VectorXd::Zero(N);


    for (int k = 0; k < nsteps; k++) {
        G[0] = a(h) * dt * gL(time[k]) / (h * h);
        G[N - 1] = a(1 - h) * dt * gR(time[k]) / (h * h);

        u.col(k + 1).segment(1, N) =
            u.col(k).segment(1, N) - dt * A * u.col(k).segment(1, N) + G;

        time[k + 1] = (k + 1) * dt;

        u.col(k + 1)[0] = gL(time[k + 1]);

        u.col(k + 1)[N + 1] = gR(time[k + 1]);
    }


    //// CMEA_END_TEMPLATE
    return std::make_pair(u, time);
}

