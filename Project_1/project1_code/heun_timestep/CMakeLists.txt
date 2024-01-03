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

SET(SRC "heun.cpp") # Source
ADD_EXECUTABLE(heun ${SRC})
add_dependencies(heun Eigen)

SET(SRC "heunconv.cpp") # Source
ADD_EXECUTABLE(heunconv ${SRC})
add_dependencies(heunconv Eigen)