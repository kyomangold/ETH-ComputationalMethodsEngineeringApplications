#pragma once
#include <Eigen/Core>


//! The gradient of the shape function (on the reference element)
//! 
//! We have three shape functions
//!
//! @param i integer between 0 and 2 (inclusive). Decides which shape function to return.
//! @param x x coordinate in the reference element.
//! @param y y coordinate in the reference element.
inline Eigen::Vector2d gradientLambda(const int i, double x, double y) {
    //// CMEA_START_TEMPLATE
    return Eigen::Vector2d(-1 + (i > 0) + (i==1),
                           -1 + (i > 0) + (i==2));
    //// CMEA_RETURN_TEMPLATE
    //// CMEA_END_TEMPLATE
}
