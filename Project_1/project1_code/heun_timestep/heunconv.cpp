//// CMEA_START_TEMPLATE
#include <Eigen/Core>
#include <vector>
#include <iostream>
#include "writer.hpp"

/// Uses the Heun's method to compute u from time 0 to time T
/// for the ODE $u'=e^{-2t}-2u$
///
/// @param[out] u at the end of the call (solution at all time steps)
/// @param[out] time contains the time levels
/// @param[in] u0 the initial data
/// @param[in] dt the step size
/// @param[in] T the final time up to which to compute the solution.
///

void Heun(std::vector<double> & u, std::vector<double> & time,
          const double & u0, double dt, double T) {
    const unsigned int nsteps = round(T/dt);
    u.resize(nsteps+1);
    time.resize(nsteps+1);

    time[0] = 0.;
    u[0] = u0;

    for(unsigned int i = 0; i < nsteps; i++) {
        time[i+1] = (i+1)*dt;
        double k1 = std::exp(-2.*time[i])-2.*u[i];
        double k2 = std::exp(-2.*time[i+1])-2.*(u[i]+dt*k1);
        u[i+1] = u[i] + dt*0.5*(k1+k2);
    }
}

int main() {

    double T = 10.0;
    std::vector<double> dt(8);
    std::vector<double> error(8);
    const double u0 = 0.;

    for(int i=0; i<8; i++) {
        dt[i]=std::pow(0.5,i+1);
        std::vector<double> time;
        std::vector<double> u;
        Heun(u,time,u0,dt[i],T);
        double uex = T*std::exp(-2.*T);
        std::cout << time.back() << std::endl;
        error[i]=std::abs(u.back()-uex);
    }
    writeToFile("dt.txt", dt);
    writeToFile("error.txt", error);
    return 0;
}
//// CMEA_BLANKFILE_TEMPLATE
//// CMEA_END_TEMPLATE
