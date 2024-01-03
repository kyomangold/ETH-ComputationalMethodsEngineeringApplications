#pragma once
#include <Eigen/Core>
#include "load_vector.hpp"

//----------------AssembleVectorBegin----------------
//! Assemble the load vector into the full right hand side
//! for the linear system
//!
//! @param[out] F will at the end contain the RHS values for each vertex.
//! @param[in] vertices a list of triangle vertices
//! @param[in] dofs a list of the dofs' indices in each triangle
//! @param[in] f the RHS function f.
void assembleLoadVector(Eigen::VectorXd& F,
			const Eigen::MatrixXd& vertices,
			const Eigen::MatrixXi& dofs, const int& N,
			const std::function<double(double, double)>& f)
{
     const int numberOfElements = dofs.rows();

     F.resize(N);
     F.setZero();
     //// CMEA_START_TEMPLATE
     for (int i = 0; i < numberOfElements; ++i) {
         const auto& indexSet = dofs.row(i);

         const auto& a = vertices.row(indexSet(0));
         const auto& b = vertices.row(indexSet(1));
         const auto& c = vertices.row(indexSet(2));

         Eigen::VectorXd elementVector(6);
         computeLoadVector(elementVector, a, b, c, f);

         for (int i = 0; i < 6; ++i) {
            F(indexSet(i)) += elementVector(i);
         }
     }
     //// CMEA_END_TEMPLATE
}
//----------------AssembleVectorEnd----------------
