#include "forward_euler.hpp"
#include "run_boundaries.hpp"


int main(int, char**) {
    runBoundaries(forwardEuler, "forward_euler");
}
