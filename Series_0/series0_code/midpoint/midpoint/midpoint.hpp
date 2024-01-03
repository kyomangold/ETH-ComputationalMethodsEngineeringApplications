#pragma once
///
/// This is the type of a function taking as parameter a double, and
/// return a double
///
typedef double(*FunctionPointer)(double);

///
/// Computes the midpoint rule to approximate the integral
///
/// \param a the left endpoint
/// \param b the right endpoint
/// \param n the number of subintervals to use
/// \param f the function to compute the integral over
///
double midpoint_rule(double a, double b, int n, FunctionPointer f);
