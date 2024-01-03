#include "run_stability.hpp"
#include "crank_nicolson.hpp"
int main(int, char**) {
    runStability(crankNicolson, "crank_nicolson");
}
