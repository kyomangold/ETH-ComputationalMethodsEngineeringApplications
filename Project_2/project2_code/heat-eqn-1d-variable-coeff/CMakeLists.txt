# Mandatory information that CMAKE needs
CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)

# Project name
PROJECT (heat_1dfd)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

INCLUDE("cmake/CommonHeader.cmake")
INCLUDE("cmake/Compiler.cmake")
# So that Visual Studio will define M_PI
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
INCLUDE("cmake/GetEigen.cmake")

add_subdirectory(base_template)

add_library(create_poisson_matrix "create_poisson_matrix.cpp")
target_include_directories(create_poisson_matrix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(create_poisson_matrix Eigen::headers)

add_library(forward_euler "forward_euler.cpp")
target_include_directories(forward_euler PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(forward_euler PUBLIC Eigen::headers create_poisson_matrix)

add_library(crank_nicolson "crank_nicolson.cpp")
target_include_directories(crank_nicolson PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(crank_nicolson PUBLIC Eigen::headers create_poisson_matrix)

add_subdirectory(run_files)

add_subdirectory(unittest)
