#pragma once
#include <Eigen/Dense>

//! @param[in] N + 1 is the number of cells
//!
//! @return a vector of length N+2 initialized to initial condition u0Function
Eigen::VectorXd computeU0(int N, std::function<double(double)> u0Function);


namespace initial {
//! @return min(2 * x, 2 - 2 * x)
double zeroZero(double x);

//! @return 1 + min(2 * x, 2 - 2 * x)
double oneOne(double x);

//! @return x + min(2 * x, 2 - 2 * x)
double zeroOne(double x);

}
