#include "midpoint.hpp"

double midpoint_rule(double a, double b, int n, FunctionPointer f) {
    //// CMEA_START_TEMPLATE
    double sum = 0;
    const double h = (b - a) / n;
    for(int i = 0; i < n; ++i) {
        sum += f(a + (i+0.5) * h);
    }
    return sum * h;
    //// CMEA_RETURN_TEMPLATE
    //// CMEA_END_TEMPLATE
}
