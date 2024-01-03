#include "midpoint.hpp" // To use our library
#include "writer.hpp" // This is the output function to write to file

// We store our results in a vector
#include <vector>

// On some platforms we need to add this in order
// to get M_PI defined
#define _USE_MATH_DEFINES

// for our usual math functions and constants
#include <math.h>

double f(double x) {
    return sin(M_PI * x);
}


int main(int, char**) {
    //// CMEA_START_TEMPLATE
    const double a = 0.2;
    const double b = 1.3;
    double exact = (cos(M_PI*a) - cos(M_PI * b)) / M_PI;
    
    // We store the errors in this vector
    std::vector<double> errors;

    // We store the number of subintervals we have used here
    std::vector<int> numberOfSubintervals;

    for(int k = 4; k <= 11; k++) {
        int n = 1 << k;

        // Compute the midpoint rule here.
        double In = midpoint_rule(a, b, n, f);
        // Compute the correct error:
        double error = fabs(In - exact);
        errors.push_back(error);
        numberOfSubintervals.push_back(n);
    }
    
    // Write result to disk
    writeToFile("series0_1_d_errors.txt", errors);
    writeToFile("series0_1_d_numbers.txt", numberOfSubintervals);

    //// CMEA_END_TEMPLATE

    return 0;
}
