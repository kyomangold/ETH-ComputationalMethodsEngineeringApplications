#pragma once
#include <Eigen/Core>
#include "shape.hpp"


//! The gradient of the shape function (on the reference element) for LINEAR FEM
//! 
//! We have three shape functions
//!
//! @param i integer between 0 and 2 (inclusive). Decides which shape function to return.
//! @param x x coordinate in the reference element.
//! @param y y coordinate in the reference element.
inline Eigen::Vector2d gradientLambda(const int i, double x, double y) {
    return Eigen::Vector2d(-1 + (i > 0) + (i==1),
                           -1 + (i > 0) + (i==2));
    return Eigen::Vector2d(0,0); //remove when implemented
}

//----------------gradshapefunBegin----------------
//! The gradient of the shape function (on the reference element) for QUADRATIC FEM
//! 
//! We have six shape functions
//!
//! @param i integer between 0 and 5 (inclusive). Decides which shape function to return.
//! @param x x coordinate in the reference element.
//! @param y y coordinate in the reference element.
inline Eigen::Vector2d gradientShapefun(const int i, double x, double y) {
  Eigen::Vector2d value;
  //// CMEA_START_TEMPLATE
  switch (i){ 
  case 0: value =(4.*lambda(0,x,y)-1.)*gradientLambda(0,x,y); break;
  case 1: value =(4.*lambda(1,x,y)-1.)*gradientLambda(1,x,y); break;
  case 2: value =(4.*lambda(2,x,y)-1.)*gradientLambda(2,x,y); break;
  case 3: {
    value = 4.*(gradientLambda(0,x,y)*lambda(1,x,y)
	       +lambda(0,x,y)*gradientLambda(1,x,y));
    break;
  }
  case 4: {
    value = 4.*(gradientLambda(1,x,y)*lambda(2,x,y)
	       +lambda(1,x,y)*gradientLambda(2,x,y));
    break;
  }
  case 5: {
    value = 4.*(gradientLambda(0,x,y)*lambda(2,x,y)
	       +lambda(0,x,y)*gradientLambda(2,x,y));
    break;
  }
  }
  //// CMEA_END_TEMPLATE
  return value;
}
//----------------gradshapefunEnd----------------
