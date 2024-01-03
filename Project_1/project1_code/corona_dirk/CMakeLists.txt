# Mandatory information that CMAKE needs
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)

# Project name
PROJECT (CMEA2_Corona_project)

# Prepare for Hunter (currently unused)
include("cmake/CommonHeader.cmake")

# Set some compiler flags
include ("cmake/Compiler.cmake")

# get Eigen
include("cmake/GetEigen.cmake")

SET(SRC "forward_euler.cpp" "forwardeulersolver.hpp" "coronaoutbreak.hpp")
ADD_EXECUTABLE(forward_euler ${SRC})
add_dependencies(forward_euler Eigen)

SET(SRCb "dirk.cpp" "dirksolver.hpp")
ADD_EXECUTABLE(dirk ${SRCb})
add_dependencies(dirk Eigen)

SET(SRCc "dirkconv.cpp")
ADD_EXECUTABLE(dirkconv ${SRCc})
add_dependencies(dirkconv Eigen)

add_subdirectory(unittest)
