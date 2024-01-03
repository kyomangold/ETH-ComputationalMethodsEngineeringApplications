#include <Eigen/Core>
#include <sstream>
#include "fem_solve.hpp"
#include "writer.hpp"

double f_square(double x, double y) {
    return 2 * M_PI * M_PI * sin(M_PI *x) * sin(M_PI * y);
}



int main(int, char**) {
  try {
    Vector u;

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi triangles;
    Eigen::MatrixXi tetrahedra;

    igl::readMESH(CMEA_DATA_PATH "square_5.mesh", vertices, tetrahedra, triangles);

    solveFiniteElement(u, vertices, triangles, f_square);

    writeToFile("square_5_values.txt", u);
    writeMatrixToFile("square_5_vertices.txt", vertices);
    writeMatrixToFile("square_5_triangles.txt", triangles);

  }
  catch (std::runtime_error& e) {
    std::cerr << "An error occurred. Error message: " << std::endl;
    std::cerr << "    \"" << e.what() << "\"" << std::endl;
    return EXIT_FAILURE;
  } 
  catch (...) {
    std::cerr << "An unknown error occurred." << std::endl;
    throw ;
  }

  return EXIT_SUCCESS;

}
