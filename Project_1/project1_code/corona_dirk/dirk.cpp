#include "writer.hpp"
#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <iostream>
#include "dirksolver.hpp"


int main(int argc, char** argv) {

    double T = 450;
    int N = 100;

    if (argc > 1)  {
        // Read N from command line
        // We use atof because we want to allow things like 1e7
        N = int(atof(argv[1]));
    }

    // We combine our system unknowns U = [S, E, I, R] and D in the
    // problem vector of unknowns u = [S, E, I, R, D].
    std::vector<std::vector<double> > u(5);
    std::vector<double> time(N + 1, 0);
    
    for (int i=0; i<u.size(); ++i){
        u[i].resize(N + 1, 0);
        u[i][0] = 0.;
    }

    u[0][0] = 500;
    u[2][0] = 1;

    CoronaOutbreak outbreak;
    DIRKSolver dirk(outbreak);
    dirk.solve(u, time, T, N);

    writeToFile("S.txt", u[0]);
    writeToFile("E.txt", u[1]);
    writeToFile("I.txt", u[2]);
    writeToFile("R.txt", u[3]);
    writeToFile("D.txt", u[4]);
    writeToFile("time.txt", time);
}




