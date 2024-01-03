#include <Eigen/Core>
#include <iostream>
#include "writer.hpp"
#include <cmath>
#include <stdexcept>
#include <functional>
//! Vector type
typedef Eigen::VectorXd Vector;

//----------------UpwindFDBegin----------------
/// Uses forward Euler and upwind finite differences to compute u from time 0 to time T 
///
/// @param[out] u solution at all time steps up to time T, each column corresponding to a time step (including the initial condition as first column)
/// @param[out] time the time levels
/// @param[in] u0 the initial conditions, as column vector
/// @param[in] dt the time step size
/// @param[in] T the final time at which to compute the solution (which we assume to be a multiple of dt)
/// @param[in] N the number of grid points, INCLUDING BOUNDARY POINTS
/// @param[in] a the advection velocity
/// @param[in] xL left limit of the domain
/// @param[in] xR right limit of the domain
///

void UpwindFD(Eigen::MatrixXd & u,  Vector & time, const Vector u0, double dt, double T, int N,
              const std::function<double(double)>& a, double xL, double xR)
{
    const unsigned int nsteps = round(T/dt);
    const double dx = (xR - xL)/(N-1);
    u.resize(N + 2,nsteps+1);
    time.resize(nsteps+1);

    /* Initialize u */
    //// CMEA_START_TEMPLATE
    u.col(0).segment(1,N) = u0;
	u(0,0) = u0[0];
	u(N+1,0) = u0[N-1];
    //// CMEA_END_TEMPLATE
	
    /* Main loop */
    //// CMEA_START_TEMPLATE
    for(unsigned k=0; k<nsteps; k++)
    {
        for(unsigned j=1; j<N+1; j++) {

            double x = xL + (j-1)*dx;
            u(j,k+1) = u(j,k) - 0.5*dt/dx*( a(x)*(u(j+1,k)-u(j-1,k)) - fabs(a(x))*(u(j+1,k)-2*u(j,k)+u(j-1,k)));

        }
		/* Outflow boundary conditions */
		u(0, k+1) = u(1, k+1);
		u(N+1, k+1) = u(N, k+1);
        time[k+1]=(k+1)*dt;
    }
    //// CMEA_END_TEMPLATE
}
//----------------UpwindFDEnd----------------

//----------------CenteredFDBegin----------------
/// Uses forward Euler and centered finite differences to compute u from time 0 to time T 
///
/// @param[out] u solution at all time steps up to time T, each column corresponding to a time step (including the initial condition as first column)
/// @param[out] time the time levels
/// @param[in] u0 the initial conditions, as column vector
/// @param[in] dt the time step size
/// @param[in] T the final time at which to compute the solution (which we assume to be a multiple of dt)
/// @param[in] N the number of grid points, INCLUDING BOUNDARY POINTS
/// @param[in] a the advection velocity
/// @param[in] xL left limit of the domain
/// @param[in] xR right limit of the domain
///

void CenteredFD(Eigen::MatrixXd & u,  Vector & time, const Vector u0, double dt, double T, int N, const std::function<double(double)>& a, double xL, double xR)
{
    const unsigned int nsteps = round(T/dt);
    const double dx = (xR-xL)/(N-1);
    u.resize(N + 2,nsteps+1);
    time.resize(nsteps+1);

    /* Initialize u */
    //// CMEA_START_TEMPLATE
    u.col(0).segment(1,N) = u0;
	u(0,0) = u0[0];
	u(N+1,0) = u0[N-1];
    //// CMEA_END_TEMPLATE
	
    time[0]=0;
    
    /* Main loop */
    //// CMEA_START_TEMPLATE
    for(unsigned k=0; k<nsteps; k++)
    {
        for(unsigned j=1; j<N+1; j++) {
            double x = xL + (j-1)*dx;
            u(j,k+1) = u(j,k) - dt/(2.*dx)*a(x)*(u(j+1,k)-u(j-1,k));
        }
        
		/* Outflow boundary conditions */
		u(0, k+1) = u(1, k+1);
		u(N+1, k+1) = u(N, k+1);
        time[k+1]=(k+1)*dt;
    }
    //// CMEA_END_TEMPLATE
}
//----------------CenteredFDEnd----------------


//----------------UpwindSrcFDBegin----------------
/// Uses forward Euler and centered finite differences to compute u from time 0 to time T 
///
/// @param[out] u solution at all time steps up to time T, each column corresponding to a time step (including the initial condition as first column)
/// @param[out] time the time levels
/// @param[in] u0 the initial conditions, as column vector
/// @param[in] dt the time step size
/// @param[in] T the final time at which to compute the solution (which we assume to be a multiple of dt)
/// @param[in] N the number of grid points, INCLUDING BOUNDARY POINTS
/// @param[in] a the advection velocity
/// @param[in] xL left limit of the domain
/// @param[in] xR right limit of the domain
/// @param[in] f the source term
void UpwindFDwithSource(Eigen::MatrixXd & u,  Vector & time, const Vector u0, double dt, double T, int N,
              const std::function<double(double)>& a, double xL, double xR, const std::function<double(double,double)>& f)
{
    const unsigned int nsteps = round(T/dt);
    const double dx = (xR-xL)/(N-1);
    u.resize(N + 2,nsteps+1);
    time.resize(nsteps+1);

    /* Initialize u */
    //// CMEA_START_TEMPLATE
    u.col(0).segment(1,N) = u0;
	u(0,0) = u0[0];
	u(N+1,0) = u0[N-1];
    //// CMEA_END_TEMPLATE
	
    time[0]=0;
    
    /* Main loop */
    //// CMEA_START_TEMPLATE
    for(unsigned k=0; k<nsteps; k++)
    {
        for(unsigned j=1; j<N+1; j++) {
            double x = xL + (j-1)*dx;
            u(j,k+1) = dt*f(x,time[k]) + u(j,k) - 0.5*dt/dx*( a(x)*(u(j+1,k)-u(j-1,k)) - fabs(a(x))*(u(j+1,k)-2*u(j,k)+u(j-1,k)));
        }
	
		/* Outflow boundary conditions */
		u(0, k+1) = u(1, k+1);
		u(N+1, k+1) = u(N, k+1);
        time[k+1]=(k+1)*dt;
    } 
    //// CMEA_END_TEMPLATE
}
//----------------UpwindSrcFDEnd----------------


/* Initial condition: rectangle */
double U0(double x) { 
    if(x<0.25 || x> 0.75)
        return 0;
    else
        return 2.;
}

int main(int, char**) {
    double T = 2.;
    double dt = 0.002; // Change this for timestep comparison
    int N=101;
    double xL = 0;
    double xR = 5;

    auto a = [](double x) {
    //// CMEA_START_TEMPLATE
	return 1 + 8*std::exp(-x/2);
    //// CMEA_RETURN_TEMPLATE
    //// CMEA_END_TEMPLATE
    };

    auto f = [&](double x, double t) {
    //// CMEA_START_TEMPLATE
		double peak_t = 0.2;
		double peak_x = 1;
		double ampl = 200;
		return ampl*std::exp(-1000*(x-peak_x)*(x-peak_x))*std::pow(std::sin(t*2*M_PI), 16);
    //// CMEA_RETURN_TEMPLATE
    //// CMEA_END_TEMPLATE
    };

    Vector u0(N);
    double h=(xR-xL)/(N-1);
    /* Initialize u0 */
    for(int i=0;i<u0.size();i++)
        u0[i] = U0(xL + h*(double)i);

    Eigen::MatrixXd u_upwind;
    Vector time_upwind;
    UpwindFD(u_upwind,time_upwind,u0,dt,T,N,a,xL,xR);
    writeToFile("time_upwind.txt",time_upwind);
    writeMatrixToFile("u_upwind.txt", u_upwind);

    Eigen::MatrixXd u_upwind_src;
    Vector time_upwind_src;
    UpwindFDwithSource(u_upwind_src,time_upwind_src,u0,dt,T,N,a,xL,xR,f);
    writeToFile("time_upwind_src.txt",time_upwind_src);
    writeMatrixToFile("u_upwind_src.txt", u_upwind_src);

    Eigen::MatrixXd u_centered;
    Vector time_centered;
    CenteredFD(u_centered,time_centered,u0,dt,T,N,a,xL,xR);
    writeToFile("time_centered.txt",time_centered);
    writeMatrixToFile("u_centered.txt", u_centered);

}
