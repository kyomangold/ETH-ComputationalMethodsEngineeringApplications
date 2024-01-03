#include <Eigen/Core>
#include <Eigen/LU> 
#include <vector>
#include <iostream>
#include "writer.hpp"
#include <stdexcept>


///
/// Uses the Newton method to solve 
///
///    vNext = v + (h/2) [F(v) - F(vNext)]
///
/// for vNext given the initial guess initialGuess.
///
/// The results will be stored to vNext.
///
/// @param[out] vNext at the end of the call, this will have vNext
/// @param[in] initialGuess the initial guess to use for newton
/// @param[in] v the previously obtained v (ie. vn)
/// @param[in] k1 the reaction factor k1 (see exercise)
/// @param[in] k2 the reaction factor k2 (see exercise)
/// @param[in] dt the step size
/// @param[in] tolerance the limit on when to stop the newton iterations
/// @param[in] maxIterations maximum number of iterations to perform.
///
//  HINT:
//  Your code should have the structure:
/* for( max number of iterations ) {
 	compute residual vector
	if( method can stop ) (*)
		return;
	Compute DG (**)
	solve DG deltaX = residual (***)
	vNext = ...
   }
   throw runtime _error (tolerance not reached)

	* depending on residual.norm()
	** you can set the matrix with with 
	   DG << a11, a12, a13, a14,
	         a21, a22, a23, a24,
	         a31, a32, a33, a34,
	         a41, a42, a43, a44;
	   Where aij is the value of the i-th row and j-th column.
	*** you can use: Vector4d deltaX = DG.lu().solve(residual);
*/


//----------------NewtonBegin----------------

/// Computes F(u) where F is as in the exercise
///
void computeF(Eigen::Vector4d& F, const Eigen::Vector4d& u, double k1, double k2) {
    //// CMEA_START_TEMPLATE
    F << -k1 * u[0]*u[1] + k2*u[2]*u[3],
	 -k1 * u[0]*u[1] + k2*u[2]*u[3],
	  k1 * u[0]*u[1] - k2*u[2]*u[3],
	  k1 * u[0]*u[1] - k2*u[2]*u[3];
    //// CMEA_END_TEMPLATE
}

/// Computes DF and stores it to DF
///
void computeDF(Eigen::Matrix4d& DF, const Eigen::Vector4d& v, double k1, double k2) {
    //// CMEA_START_TEMPLATE
    DF << -k1*v[1], -k1*v[0],  k2*v[3],  k2*v[2],
	  -k1*v[1], -k1*v[0],  k2*v[3],  k2*v[2],
	   k1*v[1],  k1*v[0], -k2*v[3], -k2*v[2],
	   k1*v[1],  k1*v[0], -k2*v[3], -k2*v[2];
    //// CMEA_END_TEMPLATE
}


void newtonSolve(Eigen::Vector4d& vNext,
	const Eigen::Vector4d& initialGuess,
	const Eigen::Vector4d& v,
	double k1, double k2, double dt, double tolerance,
	int maxIterations) {
    
	// F(v), notice that we only have to compute Fv once
	Eigen::Vector4d Fv;
	computeF(Fv, v, k1, k2);

	vNext = initialGuess;

	//// CMEA_START_TEMPLATE
	for (int iteration = 0; iteration < maxIterations; ++iteration) {
		// F(vNext)
		Eigen::Vector4d FvNext;
		computeF(FvNext, vNext, k1, k2);

		// Compute the residual, this gives us the error. (residual = G(u))
		Eigen::Vector4d residual = v + (dt / 2.0) * (Fv + FvNext) - vNext;
		if (residual.norm() <= tolerance) {
		    return;
		}

		// To store the Jacobi matrix
		Eigen::Matrix4d DF;
		computeDF(DF, vNext, k1, k2);


		Eigen::Matrix4d DG = 0.5 * dt * DF - Eigen::Matrix4d::Identity();
	
		// Solves DG deltaX = G(u) for deltaX
		Eigen::Vector4d deltaX = DG.lu().solve(-residual);

		vNext = vNext + deltaX;
	    }
	//// CMEA_END_TEMPLATE

	// If we reach this point, something wrong happened.
	throw std::runtime_error("Did not reach tolerance in Newton iteration");
}
//----------------NewtonEnd----------------

//----------------mainBegin----------------
int main() {
    double T = 20.0;
    size_t N = 1000;
    double dt = T / N;

    // This makes saving the results a lot easier. u1[n] will hold u1 at timestep n
    std::vector<double> u1(N+1), u2(N+1), u3(N+1), u4(N+1), time(N+1);

    // Set initial start values
    u1[0] = 0.2; 
    u2[0] = 0.3;
    u3[0] = 0.1;
    u4[0] = 0.4;
    time[0] = 0;

    // From the exercise
    double k1 = 0.4; 
    double k2 = 0.1;

    for (int n = 0; n < N; ++n) {
	Eigen::Vector4d vPrevious;
	vPrevious << u1[n], u2[n], u3[n], u4[n];

	Eigen::Vector4d vNext;

	//// CMEA_START_TEMPLATE
	newtonSolve(vNext, vPrevious, vPrevious, k1, k2, dt, 1e-15, 100);
	//// CMEA_END_TEMPLATE

	u1[n + 1] = vNext[0];
	u2[n + 1] = vNext[1];
	u3[n + 1] = vNext[2];
	u4[n + 1] = vNext[3];

	time[n + 1] = time[n] + dt;
    }
    writeToFile("u1.txt", u1);
    writeToFile("u2.txt", u2);
    writeToFile("u3.txt", u3);
    writeToFile("u4.txt", u4);
    writeToFile("time.txt", time);
    return 0;
}
//----------------mainEnd----------------
