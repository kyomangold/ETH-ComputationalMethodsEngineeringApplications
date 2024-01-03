#include <vector>
#include <iostream>
#include <cmath>
#include "writer.hpp"

double f(double t, double u) {
    //// CMEA_START_TEMPLATE
    return std::exp(-2*t) - 2*u;
    //// CMEA_RETURN_TEMPLATE
    //// CMEA_END_TEMPLATE
}


/// Uses the SSP RK3 method to compute u from time 0 to time T 
/// for the ODE $u'=e^{-2t}-2u$
///
/// @param[out] u at the end of the call (solution at all time steps)
/// @param[out] time contains the time levels
/// @param[in] u0 the initial data
/// @param[in] dt the step size
/// @param[in] T the final time up to which to compute the solution.
///

//----------------SSPRK3Begin----------------
void SSPRK3(std::vector<double> & u, std::vector<double> & time,
          const double & u0, double dt, double T) {
    const unsigned int nsteps = std::round(T/dt);
    u.resize(nsteps+1);
    time.resize(nsteps+1);

    // Write your SSPRK3 code here
    //// CMEA_START_TEMPLATE
    time[0] = 0.;
    u[0] = u0;

    for(unsigned int i = 0; i < nsteps; i++)
    {
        time[i+1] = (i+1)*dt;
        double k1 = f(time[i], u[i]);
        double k2 = f(time[i+1], u[i] + dt*k1);
	double k3 = f(time[i] + 0.5*dt, u[i] + 0.25*dt*(k1 + k2));
        u[i+1] = u[i] + dt*( (1./6)*(k1 + k2) + (2./3)*k3);
    }
    //// CMEA_END_TEMPLATE
}
//----------------SSPRK3End----------------

int main(int argc, char** argv) {

    double T = 10.0;

    double dt = 0.5;

    // To make some plotting easier, we take the dt parameter in as an optional
    // parameter.
    if (argc == 2) {
        dt = atof(argv[1]);
    }

    const double u0 = 0.;
    std::vector<double> time;
    std::vector<double> u;
    SSPRK3(u,time,u0,dt,T);

    writeToFile("solution.txt", u);
    writeToFile("time.txt",time);

    return 0;
}
