#include "writer.hpp"
#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <iostream>
#include <chrono>
#include "coronaoutbreak.hpp"
#include "dirksolver.hpp"


//----------------mainBegin----------------
int main(int argc, char** argv) {

    double T = 101;
    CoronaOutbreak outbreak(0,0,0,0,0,0.03,0.02,0,0);
    std::vector<double> u0(5);
    u0[0] = 500;
    u0[1] = 0;
    u0[2] = 0;
    u0[3] = 0;
    u0[4] = 0;
    
    // Compute the exact solution for the parameters above
    std::vector<double> exact = outbreak.computeExactNoCorona(T, u0[0]);

    // Initialize solver object for the parameters above
    DIRKSolver dirkSolver(outbreak);

    int minExp = 0;
    int maxExp = 8;
    int countExponents = maxExp - minExp +1;
    std::vector<double> numbers(countExponents);
    std::vector<double> walltimes(countExponents);
    std::vector<double> errors(countExponents);

    //// CMEA_START_TEMPLATE
    std::vector<std::vector<double> > u(5);
    int baseN = 200;

    for (int i = 0; i < countExponents; i++) {
        int N = (1 << (i + minExp) ) * baseN;
        std::cout << "Running for N = " << N << "..." << std::endl;

        std::vector<double> time(N + 1, 0);
        u[0].resize(N + 1, 0);
        u[1].resize(N + 1, 0);
        u[2].resize(N + 1, 0);
        u[3].resize(N + 1, 0);
        u[4].resize(N + 1, 0);

        u[0][0] = u0[0];
        u[1][0] = u0[1];
        u[2][0] = u0[2];
        u[3][0] = u0[3];
        u[4][0] = u0[4];

        auto begin = std::chrono::high_resolution_clock::now();
        dirkSolver.solve(u, time, T, N);
        auto end = std::chrono::high_resolution_clock::now();
       
       	int N_metric = u.size();
       	std::cout << "Approx sol: ";
        for (int k = 0; k< N_metric; k++){
            errors[i] += std::abs(u[k].back() - exact[k]);
        	std::cout << u[k].back() << " ";
        }
        std::cout << std::endl;

        numbers[i] = N;
        walltimes[i] = std::chrono::duration_cast<std::chrono::nanoseconds>
            (end - begin).count();
    }
    //// CMEA_END_TEMPLATE

    writeToFile("numbers.txt", numbers);
    writeToFile("errors.txt", errors);
    writeToFile("walltimes.txt", walltimes);

}
//----------------mainEnd----------------
