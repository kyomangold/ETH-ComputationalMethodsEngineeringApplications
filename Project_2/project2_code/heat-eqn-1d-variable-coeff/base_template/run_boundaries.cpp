#include "run_boundaries.hpp"

#include "coefficient_functions.hpp"
#include "initial_data.hpp"
#include "writer.hpp"
#include "boundaries.hpp"

void runBoundaries(const SolveFunction& solver,
    const std::string& solverName) {

    using namespace boundaries;
    std::vector<BoundaryCondition> boundaryConditions = {
        BoundaryCondition(zero, zero),
        BoundaryCondition(zero, one),
        BoundaryCondition(exponential, exponential)
    };
    
    using namespace initial;
    std::vector<std::function<double(double)>> initialConditions = {
        zeroZero,
        zeroOne,
        oneOne
    };

    const auto N = 63;

//    const auto u0 = computeU0(N);

    const auto T = 0.25;

    const auto dt = 0.5 / ((N + 1) * (N + 1));

    const auto a = coefficient_functions::a1;


    for (size_t k = 0; k < boundaryConditions.size(); ++k) {

        const auto bc = boundaryConditions[k];
        const auto gL = bc.first;
        const auto gR = bc.second;

        const auto u0 = computeU0(N, initialConditions[k]);

        auto solutionPair = solver(u0, dt, T, N, gL, gR, a);

        auto u = solutionPair.first;
        auto time = solutionPair.second;

        writeToFile(solverName + "_boundaries_t_g" + std::to_string(k)
            + ".txt", time);


        writeMatrixToFile(solverName + "_boundaries_u_g" + std::to_string(k)
            + ".txt", u);

    }
}

