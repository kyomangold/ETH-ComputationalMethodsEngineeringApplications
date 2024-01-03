#pragma once
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

//! Sparse Matrix type. Makes using this type easier.
using SparseMatrix =  Eigen::SparseMatrix<double>;


//! Create the 1D Poisson matrix
//! @param[in] N the number of interior points
//! @param[in] a the coefficient function a
//!
//! @returns the Poisson matrix.
SparseMatrix createPoissonMatrix(int N, const std::function<double(double)>& a);
