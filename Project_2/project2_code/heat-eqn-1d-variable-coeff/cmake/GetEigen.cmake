#
# include, download external libraries
#
include(ExternalProject)
#  all external projects should use the same compiler
set(EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}" "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}" )

# eigen
find_package(Eigen3)

if (EIGEN3_FOUND)
        add_library(Eigen::headers INTERFACE IMPORTED)
        set_property(TARGET Eigen::headers PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${EIGEN3_INCLUDE_DIR})
else()
	SET(DOWNLOADING_EIGEN ON)
	#  if not found system wide download
	message("-- Downloading Eigen3")
	ExternalProject_Add(
            EigenDownload
	    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.zip
	    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen
	    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen_install
	    CMAKE_ARGS ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/Eigen_install)

        add_library(Eigen::headers INTERFACE IMPORTED)
        add_dependencies(Eigen::headers EigenDownload)
        INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/Eigen_install/include/eigen3)
endif()
