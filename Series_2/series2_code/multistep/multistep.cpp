#include <writer.hpp>
#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <iostream>
#include <cmath>

///
/// Computes N iterations of Adam-Bashforth 2
/// @param[out] y will contain the values of y (not y'!) at
///             the end. should have length N+1
/// @param[out] time will contain the times at the end, should have length N+1
/// @param[in] T the end time
/// @param[in] N the number of iterations to perform
/// @param[in] y0 the initial y value y(0)
/// @param[in] dydt0 the initial value y'(0)
///
//----------------AB2Begin----------------
void adamBashforth2(std::vector<double>& y, std::vector<double>& time,
    double T, int N, double y0, double dydt0) {
    const double dt = T / N;
    // First compute u_1 using FE
    // Compute Adam-Bashforth for u_2, u_3, ...

    //// CMEA_START_TEMPLATE
    Eigen::Matrix2d A;
    A << 0, 1,
        -6, -5;

    Eigen::Vector2d uPrevPrev(y0, dydt0);
    Eigen::Vector2d uPrev = uPrevPrev + dt * A * uPrevPrev;

    y[0] = y0;
    y[1] = uPrev[0];
    time[0] = 0;
    for (int i = 2; i < N + 1; ++i) {
        Eigen::Vector2d uNext = uPrev + (dt/2) * A * (-uPrevPrev + 3 * uPrev);
        uPrevPrev = uPrev;
        uPrev = uNext;
        y[i] = uNext[0];
        time[i] = time[i - 1] + dt;
    }
    //// CMEA_END_TEMPLATE
}
//----------------AB2End----------------

double exactSolution(double t) {
    //// CMEA_START_TEMPLATE
    return 5 * exp(-2 * t) - 4 * exp(-3 * t);
    //// CMEA_RETURN_TEMPLATE
    //// CMEA_END_TEMPLATE
}

int main(int, char**) {
    double T = 1;

    std::vector<double> errors;
    std::vector<double> resolutions;
    double error = 0;
    for (int k = 5; k <= 12; ++k) {
        int N = 1 << k;
        std::vector<double> y(N + 1);
        std::vector<double> time(N + 1, 0);

        double y0 = 1;
        double dydt0 = 2;

        adamBashforth2(y, time, T, N, y0, dydt0);

	// Compute error
        //// CMEA_START_TEMPLATE
        double error = std::abs(exactSolution(T) - y.back());
        //// CMEA_END_TEMPLATE
        errors.push_back(error);
        resolutions.push_back(N);

        if (k == 12) {
            writeToFile("y.txt", y);
            writeToFile("time.txt", time);
        }
    }
    writeToFile("errors_ab2.txt", errors);
    writeToFile("resolutions.txt", resolutions);
}




