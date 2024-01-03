#include <Eigen/Core>
#include <sstream>
#include "writer.hpp"
#include <igl/readMESH.h>
#include <igl/readSTL.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include "dofs.hpp"
#include "fem_solve.hpp"
#include "convergence.hpp"


typedef Eigen::VectorXd Vector;

double f_square(double x, double y) {
    return 2 * M_PI * M_PI * sin(M_PI *x) * sin(M_PI * y);
}

double uex_square(double x, double y) {
    return sin(M_PI *x) * sin(M_PI * y);
}

Eigen::Vector2d uex_grad_square(double x, double y) {
  Eigen::Vector2d grad;
  grad << M_PI*cos(M_PI *x)*sin(M_PI*y),
          M_PI*sin(M_PI *x)*cos(M_PI*y);
  return grad;
}

int main(int, char**) {
  try {
    Vector u;

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi triangles;
    Eigen::MatrixXi tetrahedra;

    igl::readMESH(CMEA_DATA_PATH "square_1.mesh", vertices, tetrahedra, triangles);

    QDofs quadraticDofs(vertices, triangles);
    solveFiniteElement(u, quadraticDofs, f_square);
    
    writeToFile("square_5_values.txt", u.segment(0,vertices.rows()));
    writeMatrixToFile("square_5_vertices.txt", vertices);
    writeMatrixToFile("square_5_triangles.txt", triangles);

    convergenceAnalysis("square", 7, f_square, uex_square, uex_grad_square);
    
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
