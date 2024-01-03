#pragma once
#include <cmath>
namespace coefficient_functions {

//! @return 0.1 for all x
double a1(double);

//! @return 1 for all x
double a2(double);

//! @return 0.5 + 0.25 sin(4*pi*x)
double a3(double x);
}
