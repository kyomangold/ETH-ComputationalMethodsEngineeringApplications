# Mandatory information that CMAKE needs
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)

# Project name
PROJECT (heat_1dfd)
INCLUDE("cmake/CommonHeader.cmake")
INCLUDE("cmake/Compiler.cmake")
# So that Visual Studio will define M_PI
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)

INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES("..")
INCLUDE("cmake/GetEigen.cmake")


INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES("..")

SET(SRC_FD "heat_1dfd.cpp") # Source
ADD_EXECUTABLE(heat_1dfd ${SRC_FD}) 
ADD_DEPENDENCIES(heat_1dfd Eigen)
