#pragma once
#include <functional>
#include <string>
#include <Eigen/Dense>


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

//! Runs different boundaries for the given solver
void runBoundaries(const SolveFunction& solver,
    const std::string& solverName);
