#pragma once
#include <Eigen/Core>
#include "coronaoutbreak.hpp"
#include <vector>

class ForwardEuler {
public:
	ForwardEuler(const CoronaOutbreak& coronaOutbreak_) :
		coronaOutbreak(coronaOutbreak_) {
		// empty
	}

	///
	/// Compute N timesteps of Forward-Euler
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
	void solve(std::vector<std::vector<double> >& u,
		std::vector<double>& time,
		double T, int N) {

		const double dt = T / N;
		Eigen::VectorXd Fu(u.size()-1);
		Eigen::VectorXd uPrevious(u.size()-1);

		for (int i = 1; i < N + 1; ++i) {
			for (int j=0; j<u.size()-1; ++j){
				uPrevious[j] = u[j][i-1];
			}

			coronaOutbreak.computeF(Fu, time[i - 1], uPrevious);
			Eigen::VectorXd uNext = uPrevious + dt * Fu;

			for (int j=0; j<u.size()-1; ++j){
				u[j][i] = uNext[j];
			}

			coronaOutbreak.computeD(u,dt,i);
			
			time[i] = time[i - 1] + dt;
		}
	}

private:
	CoronaOutbreak coronaOutbreak;
}; // end ForwardEuler class
