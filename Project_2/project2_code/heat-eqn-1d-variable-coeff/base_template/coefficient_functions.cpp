#include "coefficient_functions.hpp"

namespace coefficient_functions {

double a1(double) {
    return 0.1;
}

double a2(double) {
    return 1;
}

double a3(double x) {
    return 0.5 + 0.25 * std::sin(4 * M_PI * x);
}
}
