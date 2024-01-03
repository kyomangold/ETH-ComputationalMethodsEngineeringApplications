cmake_minimum_required (VERSION 3.1.0)

include("cmake/CommonHeader.cmake")
project("numpde_midpoint")



SET(LIBRARY_NAME "midpoint")
FILE(GLOB LIBRARY_SRC ${LIBRARY_NAME}/*.cpp ${LIBRARY_NAME}/*.hpp)
ADD_LIBRARY(${LIBRARY_NAME} ${LIBRARY_SRC})

INCLUDE_DIRECTORIES("." ${LIBRARY_NAME})

FOREACH(DIRECTORY "test_single" "test_convergence")
  FILE(GLOB SRC ${DIRECTORY}/*.cpp ${DIRECTORY}/*.hpp)

  ADD_EXECUTABLE(${DIRECTORY} ${SRC})
  TARGET_LINK_LIBRARIES(${DIRECTORY} ${LIBRARY_NAME})
ENDFOREACH()
