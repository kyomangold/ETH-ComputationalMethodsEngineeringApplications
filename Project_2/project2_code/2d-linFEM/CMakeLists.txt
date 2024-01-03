# Mandatory information that CMAKE needs
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)

# Project name
PROJECT (2d-linFEM)


IF(NOT WIN32)
    SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11")
    SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Either: '', 'Release' or 'Debug'")
ENDIF()


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

find_package (Threads REQUIRED)



if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    # (this is only applicable on clang)
    # ignore some mathgl warnings
    add_definitions(-Wno-deprecated-register -Wno-return-type-c-linkage)
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    add_definitions(-Wno-deprecated-declarations -Wno-ignored-attributes)
endif()


INCLUDE(ExternalProject)

ExternalProject_Add(
    Eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.zip
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

ExternalProject_Add(
    IGL
    URL https://github.com/libigl/libigl/archive/master.zip
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/IGL
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

    
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/Eigen)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/IGL/include)
INCLUDE_DIRECTORIES(".")

ADD_SUBDIRECTORY(unittest)

ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
ADD_DEFINITIONS(-DCMEA_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/")

FILE(GLOB SRC *.cpp *.hpp)
ADD_EXECUTABLE(fem2d ${SRC})
TARGET_LINK_LIBRARIES(fem2d ${THREADS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

ADD_DEPENDENCIES(fem2d IGL Eigen )


