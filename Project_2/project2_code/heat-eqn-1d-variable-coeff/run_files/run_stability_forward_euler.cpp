#include "run_stability.hpp"
#include "forward_euler.hpp"
#include "coefficient_functions.hpp"
#include "initial_data.hpp"
#include "writer.hpp"


int main(int, char**) {
    runStability(forwardEuler, "forward_euler");
}
