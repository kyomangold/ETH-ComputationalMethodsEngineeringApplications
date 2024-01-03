#include "run_stability.hpp"
#include <vector>
#include "writer.hpp"
#include "coefficient_functions.hpp"
#include "boundaries.hpp"
#include "initial_data.hpp"

#include <iostream>

void runStability(const SolveFunction& solver,
    const std::string& solverName) {

    std::vector<std::function<double(double)>> coefficients = {
        coefficient_functions::a1,
        coefficient_functions::a2,
        coefficient_functions::a3
    };

    const auto N = 127;
    const auto dxSquared = 1.0 / ((N + 1) * (N + 1));

    std::vector<double> timesteps = {
        128.0 * dxSquared / 2.0,
        8.0 * dxSquared / 2.0,
        1.0 * dxSquared / 2.0
    };



    using namespace initial;
    const auto u0 = computeU0(N, zeroZero);

    const auto T = 0.25;

    for (size_t i = 0; i < coefficients.size(); ++i) {
        for (size_t j = 0; j < timesteps.size(); ++j) {


            const auto& a = coefficients[i];
            auto dt = timesteps[j];
            std::cout << "Running with a_" << i << ", dt =  " << dt << std::endl;

            const auto gL = boundaries::zero;
            const auto gR = boundaries::zero;

            auto solutionPair = solver(u0, dt, T, N, gL, gR, a);

            auto u = solutionPair.first;
            auto time = solutionPair.second;

            writeToFile(solverName + "_stability_t_a" + std::to_string(i)
                + "_dt" + std::to_string(j)
                + ".txt", time);


            writeMatrixToFile(solverName + "_stability_u_a" + std::to_string(i)
                + "_dt" + std::to_string(j)
                + ".txt", u);

        }
    }
}

