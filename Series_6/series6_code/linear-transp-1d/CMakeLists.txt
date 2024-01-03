# Mandatory information that CMAKE needs
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)

# Project name
PROJECT (CMAEIIS5FD1D)
if(CMAKE_VERSION VERSION_GREATER 3.1.0 OR CMAKE_VERSION VERSION_EQUAL 3.1.0)
    # use c++11
    # only valid for new platforms
    set(CMAKE_CXX_STANDARD 11)
else()
    if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        # for older cmake versions
        # (note, this CXX flag is only valid for clang and gcc, for MSVC it is not needed)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()
endif()
# So that Visual Studio will define M_PI
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)

INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES("..")
INCLUDE("cmake/GetEigen.cmake")

SET(SRC_LT "linear_transport.cpp") # Source
ADD_EXECUTABLE(linear_transport ${SRC_LT}) 

ADD_DEPENDENCIES(linear_transport Eigen)
