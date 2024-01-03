//// CMEA_START_TEMPLATE
#include <vector>
#include <iostream>
#include <cmath>
#include "writer.hpp"

/// Uses the SSPRK3's method to compute u from time 0 to time T 
/// for the ODE $u'=e^{-2t}-2u$
///
/// @param[out] u at the end of the call (solution at all time steps)
/// @param[out] time contains the time levels
/// @param[in] u0 the initial data
/// @param[in] dt the step size
/// @param[in] T the final time up to which to compute the solution.
///

double f(double t, double u) {
    return std::exp(-2*t) - 2*u;
}

void SSPRK3(std::vector<double> & u, std::vector<double> & time,
          const double & u0, double dt, double T) {
    const unsigned int nsteps = std::round(T/dt);
    u.resize(nsteps+1);
    time.resize(nsteps+1);
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
}

//----------------mainBegin----------------
int main() {

    double T = 10.0;
    std::vector<double> dt(8);
    std::vector<double> error(8);
    const double u0 = 0.;
    for(int i=0; i<8; i++) {
        dt[i]=std::pow(0.5,i);
        std::vector<double> time;
        std::vector<double> u;
        SSPRK3(u,time,u0,dt[i],T);
        double uex = T*std::exp(-2.*T);
        error[i]=std::abs(u.back()-uex);
    }
    writeToFile("dt.txt", dt);
    writeToFile("error.txt", error);
    return 0;
}
//----------------mainEnd----------------
//// CMEA_BLANKFILE_TEMPLATE
//// CMEA_END_TEMPLATE
