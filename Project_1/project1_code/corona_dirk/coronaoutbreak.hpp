#pragma once
#include <functional>
#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <iostream>

class CoronaOutbreak {
public:

    ///
    /// Evaluates right hand side of the system of ODEs U' = F(t,U)
    ///    S' = Pi*S - alpha*S*I - beta*S*E - delta*S + d*R
    ///    E' = alpha*S*I + beta*S*E - e*E - gamma *E - delta*E
    ///    I' = e*E -f*I - (sigma +delta) *I
    ///    R' = f*I + gamma *E - d*R - delta *R
    /// @param[in/out] F evaluation of RHS
    /// @param[in] U VectorXd with [S,E,I,R]
    /// @param[in] t current time
    ///
    //----------------FStart----------------
    void computeF(Eigen::VectorXd& F, double t, Eigen::VectorXd U) {
        double S = U[0];
        double E = U[1];
        double I = U[2];
        double R = U[3];

        F[0] = Pi * S - alpha(t) * S * I - beta(t) * S * E - delta * S + d * R;
        F[1] = alpha(t) * S * I + beta(t) * S * E  - e * E - gamma * E - delta * E;
        F[2] = e * E - f * I - sigma * I - delta * I;
        F[3] = f * I + gamma * E - d * R - delta * R; 
    }
    //----------------FEnd---------------

    ///
    /// Computes Jacobian matrix of F
    /// @param[in/out] J Jacobian matrix
    /// @param[in] t time
    /// @param[in] U VectorXd with [S,E,I,R]
    ///
    //----------------JFStart----------------
    void computeJF(Eigen::MatrixXd& J, double t, Eigen::VectorXd U) {
        //// CMEA_START_TEMPLATE
        double S = U[0];
        double E = U[1];
        double I = U[2];
        double R = U[3];

        J.setConstant(0.0);

        J(0, 0) = Pi - delta - alpha(t) * I - beta(t) * E;
        J(0, 1) = -beta(t) * S;
        J(0, 2) = -alpha(t) * S;
        J(0, 3) = d;

        J(1, 0) = alpha(t) * I + beta(t) * E;
        J(1, 1) = - (e + gamma + delta) + beta(t) * S;
        J(1, 2) = alpha(t) * S;
        
        J(2, 1) = e;
        J(2, 2) = -f - sigma - delta;
        
        J(3, 1) = gamma;
        J(3, 2) = f;
        J(3, 3) = -(d + delta);
        //// CMEA_END_TEMPLATE
    }
    //----------------JFEnd----------------

    ///
    /// Computes Dead class D
    /// @param[in/out] u vector of vectors of variables. Dimensions: 5 x # time steps
    /// @param[in] dt time step, n integer corresponding to t_n+1.
    ///
    //----------------DStart----------------
    void computeD(std::vector< std::vector<double> >& u, const double dt, int n) {
        //// CMEA_START_TEMPLATE
        double Sn = u[0][n-1];
        double En = u[1][n-1];
        double In = u[2][n-1];
        double Rn = u[3][n-1];

        double Snp1 = u[0][n];
        double Enp1 = u[1][n];
        double Inp1 = u[2][n];
        double Rnp1 = u[3][n];
        
        double Hn = delta * (Sn + En + Rn + In) + sigma * In;
        double Hnp1 = delta * (Snp1 + Enp1 + Rnp1 + Inp1) + sigma * Inp1;

        u[4][n] = u[4][n-1] + 0.5 * (Hn + Hnp1) * dt;
        //// CMEA_END_TEMPLATE
    }
    //----------------DEnd----------------

    double alpha(double t) {
        return a * (1. - r * t / (t + 1.));
    }

    double beta(double t) {
        return b * (1. - r * t / (t + 1.));
    }

    CoronaOutbreak() = default;

    CoronaOutbreak(double a_, double b_, double d_, double e_, double f_, 
        double Pi_, double delta_, double gamma_, double sigma_, double r_ = 0) :
        a(a_), b(b_), d(d_), e(e_), f(f_), Pi(Pi_), delta(delta_), gamma(gamma_), sigma(sigma_), r(r_) {
        // Empty
    }


    std::vector<double> computeExactNoCorona(double t, double S0) {
        std::vector<double> result(5);
        result[0] = S0 * std::exp( (Pi - delta) * t);
        result[1] = 0.;
        result[2] = 0.;
        result[3] = 0.;
        result[4] = (delta * S0 / (Pi - delta)) * ( std::exp( (Pi - delta) * t) - 1. );
        std::cout.precision(17);
        std::cout << "Exact solution: " << result[0] << " " << result[1]
            << " " << result[2] << " " << result[3] << " " << result[4] << std::endl;
        return result;
    }


private:

    const double r = 0.01;

    //! Birth rate
    const double Pi = 3e-5;

    //! Death rate
    const double delta = 2e-5;
    const double gamma = 1.5e-5;
    const double sigma = 2.5e-3;

    // The magic Garland constants:
    const double a = 1.5e-3;
    const double b = 0.3e-3;
    const double d = 2e-3;
    const double e = 0.5;
    const double f = 0.5;
};
