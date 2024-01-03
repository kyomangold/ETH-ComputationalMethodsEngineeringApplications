#include <Eigen/Core>
#include "stiffness_matrix.hpp"
#include "fem_solve.hpp"
#include "writer.hpp"
#include <sstream>
#include "Lshape.hpp"
#include "convergence.hpp"


double sigma(double x, double y){
  // return 1.0;
   return 0.01*(x+2)*(x+2);
}

/* 
// Data for square
double f_square(double x, double y) {
    return 2 * M_PI * M_PI * sin(M_PI *x) * sin(M_PI * y);
}

double g_square(double x, double y) {
    return 0;
}

double uex_square(double x, double y) {
    return sin(M_PI *x) * sin(M_PI * y);
}
*/

Eigen::Vector2d uex_grad_square(double x, double y) {
  Eigen::Vector2d grad;
  grad << M_PI*cos(M_PI *x)*sin(M_PI*y),
          M_PI*sin(M_PI *x)*cos(M_PI*y);
  return grad;
}

Eigen::Vector2d u_square_ex_grad(double x, double y) {
  Eigen::Vector2d gradient;
  gradient << M_PI*cos(M_PI*x)*sin(M_PI*y), M_PI*sin(M_PI*x)*cos(M_PI*y);
  return gradient;
}

int main(int, char**) {
  try {

    solveL(0.5, sigma);

    /*
    std::cout << "Convergence Analysis" << std::endl;
    convergenceAnalysis("square", 7, f_square, sigma, g_square, 0,
			uex_square, uex_grad_square);

    convergenceAnalysis("Lshape", 7, f_lshape, sigma, g_lshape, 0,
			g_lshape, g_grad_lshape);
    */

    
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
