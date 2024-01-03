#include "create_poisson_matrix.hpp"

//! Used for filling the sparse matrix.
using Triplet = Eigen::Triplet<double>;

//! Create the 1D Poisson matrix
//! @param[in] N the number of interior points
//! @param[in] a the coefficient function a
//!
//! @returns the Poisson matrix.
SparseMatrix createPoissonMatrix(int N,
    const std::function<double(double)>& a) {

    SparseMatrix A;
    //// CMEA_START_TEMPLATE
    A.resize(N, N);
    double h = 1. / (N + 1);
    std::vector<Triplet> triplets;
    auto x = Eigen::VectorXd::LinSpaced(N + 2, 0, 1);
    triplets.reserve(size_t(N + 2 * N - 2));

    for (int i = 0; i < N; ++i) {
        triplets.push_back(Triplet(i, i, 2.*a(x[i + 1]) / (h * h)));

        if (i > 0) {
            triplets.push_back(Triplet(i, i - 1, -a(x[i + 1]) / (h * h)));
        }

        if (i < N - 1) {
            triplets.push_back(Triplet(i, i + 1, -a(x[i + 1]) / (h * h)));
        }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    //// CMEA_END_TEMPLATE
    return A;
}

