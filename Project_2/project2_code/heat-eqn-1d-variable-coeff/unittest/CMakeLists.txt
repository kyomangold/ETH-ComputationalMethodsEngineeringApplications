cmake_minimum_required (VERSION 3.1.0)

ExternalProject_Add(
    GTestDownload
    URL https://github.com/google/googletest/archive/release-1.8.1.zip
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gtest_source
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/gtest_binary
    CMAKE_ARGS "-Dgtest_force_shared_crt=ON" ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/gtest_binary
    INSTALL_COMMAND  ""
    )

add_library(GTest::gtest INTERFACE IMPORTED)
add_dependencies(GTest::gtest GTestDownload)
set_property(TARGET GTest::gtest PROPERTY INTERFACE_LINK_LIBRARIES ${THREADS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} gtest$<$<AND:$<OR:$<PLATFORM_ID:Windows>,$<AND:$<PLATFORM_ID:Darwin>,$<CXX_COMPILER_ID:AppleClang>>>,$<CONFIG:Debug>>:d>)

# We need the thread package for GTest
FIND_PACKAGE(Threads REQUIRED)



LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/gtest_binary
  ${CMAKE_CURRENT_BINARY_DIR}/gtest_binary/googlemock/gtest/${CMAKE_BUILD_TYPE}
  ${CMAKE_CURRENT_BINARY_DIR}/gtest_binary/googlemock/gtest/)


ADD_EXECUTABLE(unittest
    main.cpp
    test_create_poisson_matrix.cpp
    test_forward_euler.cpp
    test_crank_nicolson.cpp)

TARGET_INCLUDE_DIRECTORIES(unittest PRIVATE SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/gtest_source/googletest/include)
TARGET_LINK_LIBRARIES(unittest GTest::gtest Threads::Threads create_poisson_matrix forward_euler crank_nicolson base_template)
