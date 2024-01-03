# Mandatory information that CMAKE needs
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)

# Project name
PROJECT (CMAEIIS1EX1)


# Prepare for Hunter (currently unused)
include("cmake/CommonHeader.cmake")

# Set some compiler flags
include ("cmake/Compiler.cmake")

# get Eigen
include("cmake/GetEigen.cmake")

SET(SRC "trapezoidal.cpp") # Source
ADD_EXECUTABLE(trapezoidal ${SRC}) # Compile as library
add_dependencies(trapezoidal Eigen)

