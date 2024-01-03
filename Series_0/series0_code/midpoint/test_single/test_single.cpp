// To use our previously written midpoint rule function
#include "midpoint.hpp"

// For printing to the terminal
#include <iostream>

// On some platforms we need to add this in order
// to get M_PI defined
#define _USE_MATH_DEFINES

// for our usual math functions and constants
#include <math.h>

// We use these two to set the precision of our output. 
#include <limits>
#include <iomanip>

double f(double x) {
    return sin(M_PI * x);
}

int main(int, char**) {
    // TODO: Compute the proper approximation here:
    double In = 0; 
    //// CMEA_START_TEMPLATE
    In = midpoint_rule(0.2, 1.3, 100, f);
    //// CMEA_END_TEMPLATE

    // Set high precision for output, easier to see what is going on:
    std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
    
    // We print out the value of the midpoint rule here:
    std::cout << "In = " << In << std::endl;

    return 0;
}
