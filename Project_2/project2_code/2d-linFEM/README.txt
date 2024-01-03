NumPDE

# Overview


# Requirements
  * A C++ compiler WITH C++11 support (GCC, Clang and MSVC tested).
  * cmake [optional]
  * EITHER MATLAB or Python (with numpy + matplotlib) for plotting.


# Building with CMake

We strongly recommend building the source using cmake.

## Building with cmake from the command line

In the template_code folder, first make a new build folder:

    mkdir build
    cd build

then issue cmake:

    cmake ..

You can build the whole program using

    make

## Running the executables

If you used cmake, you can now run the executables using

    ./fem2d

## Running the unittests

We have made some simple unittests for the code that checks if you have
implemented the functions correctly. You can run them by typing

     ./unittest/unittest

in the build folder.

NOTE: The unittests are no guarantee that your code is correct, just a small
indication!

