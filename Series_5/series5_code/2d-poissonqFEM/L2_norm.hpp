#pragma once
#include <functional>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "coordinate_transform.hpp"
#include "integrate.hpp"
#include "shape.hpp"

//----------------L2Begin----------------
//! Computes the L^2 differences between 
//! u1 (considered as coefficients for Quadratic FEM) and u2.
double computeL2Difference(const Eigen::MatrixXd& vertices,
			   const Eigen::MatrixXi& dofs,
			   const Eigen::VectorXd& u1,
			   const std::function<double(double, double)>& u2) 
{
  const int numberOfElements = dofs.rows();

  double error = 0;
  for (int i = 0; i < numberOfElements; ++i) {
    auto& idSet = dofs.row(i);

    const auto& a = vertices.row(idSet(0));
    const auto& b = vertices.row(idSet(1));
    const auto& c = vertices.row(idSet(2));

    auto coordinateTransform = makeCoordinateTransform(b - a, c - a);
    auto volumeFactor = std::abs(coordinateTransform.determinant());

    //// CMEA_START_TEMPLATE
    error += integrate([&](double x, double y) {
	Eigen::Vector2d z = coordinateTransform * Eigen::Vector2d(x, y) + Eigen::Vector2d(a(0), a(1));

	double approximateValue = 0;
	for(int j=0; j<6; j++){
	  approximateValue+= u1(idSet(j))*shapefun(j,x,y);
	}

	return std::pow(std::abs(u2(z(0), z(1)) - approximateValue), 2) * volumeFactor;
      });
    //// CMEA_END_TEMPLATE
  }

  return std::sqrt(error);
}
//----------------L2End----------------
