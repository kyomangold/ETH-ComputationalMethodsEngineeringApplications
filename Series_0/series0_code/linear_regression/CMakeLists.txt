cmake_minimum_required (VERSION 3.1.0)

include("cmake/CommonHeader.cmake")
project("numpde_regression")
include ("cmake/Compiler.cmake")

include("cmake/GetEigen.cmake")

FILE(GLOB SRC *.cpp *.hpp)

ADD_EXECUTABLE(regression ${SRC})
add_dependencies(regression Eigen)

TARGET_LINK_LIBRARIES(regression ${LIBRARY_NAME})

