#pragma once
#include <functional>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "coordinate_transform.hpp"
#include "integrate.hpp"
#include "shape.hpp"

//----------------H1Begin----------------
//! Computes the H^1 differences between 
//! u1 (considered as coefficients for Quadratic FEM) and u2grad.
double computeH1Difference(const Eigen::MatrixXd& vertices,
			   const Eigen::MatrixXi& dofs,
			   const Eigen::VectorXd& u1,
			   const std::function<Eigen::Vector2d(double, double)>& u2grad) 
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
    Eigen::Matrix2d elementMap = coordinateTransform.inverse().transpose();

    //// CMEA_START_TEMPLATE
    error += integrate([&](double x, double y) {
	Eigen::Vector2d z = coordinateTransform * Eigen::Vector2d(x, y) + Eigen::Vector2d(a(0), a(1));

        Eigen::Vector2d approximate_grad(0.,0.);
	for(int j=0; j<6; j++){
	  approximate_grad+= u1(idSet(j))* elementMap* gradientShapefun(j,x,y);
	}
	
	Eigen::Vector2d difference_grad = u2grad(z(0), z(1)) - approximate_grad;
	return difference_grad.dot(difference_grad) * volumeFactor;
      });
    //// CMEA_END_TEMPLATE
  }

  return std::sqrt(error);
}
//----------------H1End----------------
