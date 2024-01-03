#include "initial_data.hpp"

Eigen::VectorXd computeU0(int N, std::function<double(double)> u0Function) {
    Eigen::VectorXd u0Vector(N + 2);

    double h = 1. / (N + 1);

    //auto u0Function = [] (double x) {
    //    return 1 + std::min(2 * x, 2 - 2 * x);
    //};

    for (int i = 0; i < N + 2; i++) {
        u0Vector[i] = u0Function(h * i);
    }

    return u0Vector;
}

namespace initial {
    double zeroZero(double x) {
        return std::min(2 * x, 2 - 2 * x);
    }

    double zeroOne(double x) {
        return x + std::min(2 * x, 2 - 2 * x);
    }

    double oneOne(double x) {
        return 1 + std::min(2 * x, 2 - 2 * x);
    }
}
