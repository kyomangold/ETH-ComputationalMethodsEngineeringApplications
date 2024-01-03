#include "crank_nicolson.hpp"
#include "run_boundaries.hpp"
#include <Eigen/Dense>

int main(int, char**) {
    runBoundaries(crankNicolson, "crank_nicolson");
}
