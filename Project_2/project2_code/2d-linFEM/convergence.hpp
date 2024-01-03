#pragma once
#include <functional>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "fem_solve.hpp"
#include "L2_norm.hpp"
#include "H1_norm.hpp"
#include "writer.hpp"
#include <igl/readMESH.h>

//----------------convBegin----------------
//! Compute the L2 and H1-error for different meshes.
//! 
//! @param baseMeshName the basename for the mesh (eg. Square).
//!                     the method will automatically append _<n>.stl
//!
//! @param maxLevel the max mesh level to use. Will read up to and including
//!                 baseMeshName_<maxLevel>.stl
//!
//! @param f the function f as in the exercise (LHS)
//! @param g the function g as in the exercise (boundary value)
//! @param[in] sigma the function sigma as in the exercise
//! @param[in] r the parameter r from the lecture notes
//! @param exactSol the exact solution (found by eg. hand computation or googling).
//! @param exactSol_grad the gradient of the exact solution.
void convergenceAnalysis(const std::string& baseMeshName, int maxLevel,
			 const std::function<double(double, double)> f,
			 const std::function<double(double, double)>& sigma,
			 const std::function<double(double, double)> g, double r,
			 const std::function<double(double, double)> exactSol,
			 const std::function<Eigen::Vector2d(double,double)> exactSol_grad)
{
    std::vector<double> differences_L2;
    std::vector<double> differences_H1;
    std::vector<int> numberOfDegreesOfFreedom;
    for (int i = 0; i <= maxLevel; i++) {
        Vector u;
        Eigen::MatrixXd vertices;
        Eigen::MatrixXi triangles;
        Eigen::MatrixXi tetrahedra;

        std::stringstream basenameSS;
        basenameSS << baseMeshName << "_" << i;
        std::string basename = basenameSS.str();

        igl::readMESH(std::string(CMEA_DATA_PATH)
            + basename + ".mesh", vertices, tetrahedra, triangles);
        std::cout << "Computing convergence. At: " << basename << std::endl;
	// solve finite element system
	//// CMEA_START_TEMPLATE
        const int degreeOfFreedom = solveFiniteElement(u, vertices, triangles, f,
						       sigma, g , r);
	//// CMEA_END_TEMPLATE

	//compute L2-error and save it in the corresponding vector
	//// CMEA_START_TEMPLATE
        differences_L2.push_back(computeL2Difference(vertices, triangles, u, exactSol));
	//// CMEA_END_TEMPLATE

	//compute H1-error and save it in the corresponding vector
	//// CMEA_START_TEMPLATE
	differences_H1.push_back(computeH1Difference(vertices, triangles, u, exactSol_grad));
	//// CMEA_END_TEMPLATE

	// store number of dofs in vector
	//// CMEA_START_TEMPLATE
        numberOfDegreesOfFreedom.push_back(degreeOfFreedom);
	//// CMEA_END_TEMPLATE

        writeToFile(basename + "_values.txt", u);
        writeMatrixToFile(basename + "_vertices.txt", vertices);
        writeMatrixToFile(basename + "_triangles.txt", triangles);
    }

    std::stringstream errorL2Filename;
    errorL2Filename << baseMeshName << "_errors.txt";

    std::stringstream errorH1Filename;
    errorH1Filename << baseMeshName << "_errorsH1.txt";

    writeToFile(errorL2Filename.str(), differences_L2);
    writeToFile(errorH1Filename.str(), differences_H1);
    writeToFile(baseMeshName + "_resolutions.txt", numberOfDegreesOfFreedom);
}
//----------------convEnd----------------
