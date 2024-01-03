# Mandatory information that CMAKE needs
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)

# Project name
PROJECT (CMAEIIS2MULTISTEP)

# Prepare for Hunter (currently unused)
include("cmake/CommonHeader.cmake")

# Set some compiler flags
include ("cmake/Compiler.cmake")

# get Eigen
include("cmake/GetEigen.cmake")


INCLUDE_DIRECTORIES(".")
SET(SRC_FE "forward_euler.cpp") # Source
ADD_EXECUTABLE(forward_euler ${SRC_FE}) 
ADD_DEPENDENCIES(forward_euler Eigen)

SET(SRC_MULTISTEP "multistep.cpp") # Source
ADD_EXECUTABLE(multistep ${SRC_MULTISTEP}) 
ADD_DEPENDENCIES(multistep Eigen)
