#pragma once
#include <Eigen/Core>
#include <vector>
#include "coronaoutbreak.hpp"

class DIRKSolver {
public:
	DIRKSolver(const CoronaOutbreak& coronaOutbreak_)
		: coronaOutbreak(coronaOutbreak_) {

	}

	///
	/// Evaluates function G1(y) from task b)
	/// @param[out] G evaluation of function
	/// @param[in] y input to G
	/// @param[in] tn current time
	/// @param[in] Un computed value of U = [S,E,I,R] at time tn
	/// @param[in] dt timestep
	///
	//----------------G1Start----------------
	void computeG1(Eigen::VectorXd& G, Eigen::VectorXd y, double tn,
		Eigen::VectorXd Un, double dt) {
		//// CMEA_START_TEMPLATE
		int dim = y.size();
		Eigen::VectorXd Fy(dim);
		coronaOutbreak.computeF(Fy, tn + mu * dt, y);
		G = Un + dt * mu * Fy - y;
		//// CMEA_END_TEMPLATE
	}
	//----------------G1End----------------

	///
	/// Evaluates function G2(y1, y) from task b)
	/// @param[out] G evaluation of function
	/// @param[in] y input to G
	/// @param[in] tn current time
	/// @param[in] Un computed value of U = [S,E,I,R] at time tn
	/// @param[in] dt timestep
	/// @param[in] y1 computed value for first intermediate RK value
	///
	void computeG2(Eigen::VectorXd& G, Eigen::VectorXd y, double tn,
		Eigen::VectorXd Un, double dt, Eigen::VectorXd y1) {
		//// CMEA_START_TEMPLATE
		int dim = y.size();
		Eigen::VectorXd Fy1(dim), Fy(dim);
		coronaOutbreak.computeF(Fy1, tn + mu * dt, y1);
		coronaOutbreak.computeF(Fy, tn + (mu - nu)*dt, y);
		G = Un - nu * dt * Fy1 + mu * dt * Fy - y;
		//// CMEA_END_TEMPLATE
	}

	///
	/// Find a solution to JG1*x = -G1 with Newton's method
	/// @param[out] x solution to the system
	/// @param[in] Un computed value of U = [S,E,I,R] at time tn
	/// @param[in] dt timestep
	/// @param[in] tn current time
	/// @param[in] tolerance if Newton increment smaller, successfully converged
	/// @param[in] maxIterations  max Newton iterations to try before failing
	///
	void newtonSolveY1(Eigen::VectorXd& u, Eigen::VectorXd Un, double dt,
		double tn, double tolerance, int maxIterations) {
		int dim = Un.size();
		Eigen::VectorXd RHSG1(dim), Fu(dim), x(dim);
		Eigen::MatrixXd JG1(dim,dim), JFu(dim,dim);
		u = Un;
	
	   // Write your loop for the Newton solver here.
	   // You will need to use coronaOutbreak.computeF(...) 
	   // and coronaOutbreak.computeJF(...)
		//// CMEA_START_TEMPLATE

		for (int iteration = 0; iteration < maxIterations; ++iteration) {
			coronaOutbreak.computeJF(JFu, tn + mu * dt, u);
			coronaOutbreak.computeF(Fu, tn + mu * dt, u);
			Eigen::MatrixXd JG1 = dt * mu * JFu - Eigen::MatrixXd::Identity(dim,dim);
			Eigen::VectorXd RHSG1(dim);
			computeG1(RHSG1,u,tn,Un,dt);
			x = JG1.lu().solve(-RHSG1);

			if ( x.norm() <= tolerance ) {
				return;
			}

			u = u + x;
		}

		//// CMEA_END_TEMPLATE

		// If we reach this point, something wrong happened.
		throw std::runtime_error("Did not reach tolerance in Newton iteration in Y1");
	}

	///
	/// Find a solution to JG2*x = -G2 with Newton's method
	/// @param[out] x solution to the system
	/// @param[in] Un computed value of U = [S, E, I, R] at time tn
	/// @param[in] y1 previous intermediate value for RK method
	/// @param[in] dt timestep
	/// @param[in] tn current time
	/// @param[in] tolerance if Newton increment smaller, successfully converged
	/// @param[in] maxIterations  max Newton iterations to try before failing
	///
	//----------------NewtonG2Start----------------
	void newtonSolveY2(Eigen::VectorXd& v, Eigen::VectorXd Un,
		Eigen::VectorXd y1, double dt, double tn, double tolerance, int maxIterations) {

	   // Use newtonSolveY1 as a model for this
		//// CMEA_START_TEMPLATE
		int dim = Un.size();
		Eigen::VectorXd RHSG2(dim), Fv(dim), Fy(dim), x(dim);
		Eigen::MatrixXd JG2(dim,dim), JFv(dim,dim);
		coronaOutbreak.computeF(Fy,tn + mu * dt,y1);
		v = y1;

		for (int iteration = 0; iteration < maxIterations; ++iteration) {
			coronaOutbreak.computeJF(JFv, tn + (mu - nu)*dt, v);
			coronaOutbreak.computeF(Fy, tn + mu * dt, y1);
			coronaOutbreak.computeF(Fv, tn + (mu - nu)*dt, v);
			Eigen::MatrixXd JG2 = dt * mu * JFv - Eigen::MatrixXd::Identity(dim,dim);
			Eigen::VectorXd RHSG2;
			computeG2(RHSG2,v,tn,Un,dt,y1);

			x = JG2.lu().solve(-RHSG2);

			if ( x.norm() <= tolerance ) {
				return;
			}

			v = v + x;
		}

		// If we reach this point, something wrong happened.
		throw std::runtime_error("Did not reach tolerance in Newton iteration in Y2");
		//// CMEA_END_TEMPLATE
	}
	//----------------NewtonG2End----------------



	///
	/// Compute N timesteps of DIRK(2,3)
	/// @param[in/out] u should be a vector of size 5, where each
	///                component is a vector of size N+1. u[i][0]
	///                should have the initial value stored before
	///                calling the funtion
	///
	/// @param[out] time should be of length N+1
	///
	/// @param[in] T the final time
	/// @param[in] N the number of timesteps
	///
	//----------------DirkStart----------------
	void solve(std::vector<std::vector<double> >& u, std::vector<double>& time,
		double T, int N) {

		const double dt = T / N;

		// Your main loop goes here. At iteration n,
		// 1) Find Y_1 with newtonSolveY1 (resp. Y2)
		// 2) Compute U^{n+1} with F(Y1), F(Y2)
		// 3) Write the values at u[...][n] 
		// 4) Compute D and write time[n]

		//// CMEA_START_TEMPLATE
		int dim = u.size()-1;
		Eigen::VectorXd Fy1(dim), Fy2(dim);

		for (int i = 1; i < N + 1; ++i) {
			double tn = time[i - 1];
			Eigen::VectorXd uPrevious(dim);

			for(int k=0; k<dim; ++k){
				uPrevious(k) = u[k][i - 1];
			}

			Eigen::VectorXd y1(dim), y2(dim);
			newtonSolveY1(y1, uPrevious, dt, time[i - 1], 1e-10, 100);
			newtonSolveY2(y2, uPrevious, y1, dt, time[i - 1], 1e-10, 100);
			coronaOutbreak.computeF(Fy1, tn + mu * dt, y1);
			coronaOutbreak.computeF(Fy2, tn + (mu - nu)*dt, y2);
			Eigen::VectorXd uNext = uPrevious + (dt / 2.) * (Fy1 + Fy2);

			u[0][i] = uNext[0];
			u[1][i] = uNext[1];
			u[2][i] = uNext[2];
			u[3][i] = uNext[3];

			coronaOutbreak.computeD(u,dt,i);

			time[i] = time[i - 1] + dt;
		}

		//// CMEA_END_TEMPLATE
	}
	//----------------DirkEnd----------------

private:

	const double mu = 0.5 + 0.5 / sqrt(3);
	const double nu = 1. / sqrt(3);

	CoronaOutbreak coronaOutbreak;

}; // end class DIRKSolver
