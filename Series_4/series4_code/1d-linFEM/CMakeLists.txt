# Mandatory information that CMAKE needs
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)

# Project name
PROJECT (1d-linFEM)
INCLUDE("cmake/CommonHeader.cmake")
INCLUDE("cmake/Compiler.cmake")
# So that Visual Studio will define M_PI
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)

INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES("..")
INCLUDE("cmake/GetEigen.cmake")

SET(SRC_FEM "fem.cpp") # Source
ADD_EXECUTABLE(fem ${SRC_FEM})
ADD_DEPENDENCIES(fem Eigen)
