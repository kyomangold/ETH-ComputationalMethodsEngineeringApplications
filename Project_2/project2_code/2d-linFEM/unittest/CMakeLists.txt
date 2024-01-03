cmake_minimum_required (VERSION 2.8.8)



ExternalProject_Add(
    GTest
    URL https://github.com/google/googletest/archive/release-1.7.0.zip
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gtest_source
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/gtest_binary
    CMAKE_ARGS "-Dgtest_force_shared_crt=ON"
    INSTALL_COMMAND  ""
    )

# We need the thread package for GTest
FIND_PACKAGE(Threads REQUIRED)



FILE(GLOB_RECURSE SRC *.cpp)

INCLUDE_DIRECTORIES(  ${CMAKE_CURRENT_BINARY_DIR}/gtest_source/include "..")



LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/gtest_binary
                 ${CMAKE_CURRENT_BINARY_DIR}/gtest_binary/${CMAKE_BUILD_TYPE})
ADD_EXECUTABLE(unittest ${SRC})
ADD_DEPENDENCIES(unittest GTest IGL Eigen)
TARGET_LINK_LIBRARIES(unittest gtest ${THREADS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
