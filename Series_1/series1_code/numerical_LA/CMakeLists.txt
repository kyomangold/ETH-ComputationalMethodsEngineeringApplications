cmake_minimum_required (VERSION 3.1.0)

include("cmake/CommonHeader.cmake")
project("cmea_linalg")
include ("cmake/Compiler.cmake")

include("cmake/GetEigen.cmake")

FILE(GLOB SRC *.cpp *.hpp)

ADD_EXECUTABLE(linalg ${SRC})
add_dependencies(linalg Eigen)

TARGET_LINK_LIBRARIES(linalg ${LIBRARY_NAME})

