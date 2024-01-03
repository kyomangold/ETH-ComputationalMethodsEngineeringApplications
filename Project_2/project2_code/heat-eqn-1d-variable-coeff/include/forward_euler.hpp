#include <Eigen/Dense>
#include <functional>

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
std::pair<Eigen::MatrixXd, Eigen::VectorXd>
forwardEuler(const Eigen::VectorXd& u0,
    double dt,
    double T,
    int N,
    const std::function<double(double)>& gL,
    const std::function<double(double)>& gR,
    const std::function<double(double)>& a);
