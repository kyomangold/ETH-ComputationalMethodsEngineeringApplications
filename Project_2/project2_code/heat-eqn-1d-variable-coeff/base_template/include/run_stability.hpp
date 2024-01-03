#pragma once
#include <Eigen/Dense>
#include <functional>

//! Type of the solver function. The solver function takes:
//! u0, dt, T, N, gL, gR, a
using SolveFunction = std::function<std::pair<Eigen::MatrixXd, Eigen::VectorXd>
    (const Eigen::VectorXd&,
        double,
        double,
        int,
        const std::function<double(double)>&,
        const std::function<double(double)>&,
        const std::function<double(double)>&)>;

//! Runs different timesteps for the given solver
void runStability(const SolveFunction& solver,
    const std::string& solverName);
