#include "boundaries.hpp"

namespace boundaries {
double zero(double) {
    return 0;
}

double one(double) {
    return 1;
}

double exponential(double x) {
    return std::exp(-10 * x);
}

}
