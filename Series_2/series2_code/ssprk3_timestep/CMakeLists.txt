# Mandatory information that CMAKE needs
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)

# Project name
PROJECT (CMAEIIS2EX1)


# Prepare for Hunter (currently unused)
include("cmake/CommonHeader.cmake")

# Set some compiler flags
include ("cmake/Compiler.cmake")

# get Eigen
include("cmake/GetEigen.cmake")

SET(SRC "ssprk3.cpp") # Source
ADD_EXECUTABLE(ssprk3 ${SRC}) # Compile as library
ADD_DEPENDENCIES(ssprk3 Eigen)

SET(SRCb "ssprk3conv.cpp")
ADD_EXECUTABLE(ssprk3conv ${SRCb})
ADD_DEPENDENCIES(ssprk3conv Eigen)
