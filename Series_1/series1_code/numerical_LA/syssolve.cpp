#include <iostream>
#include <Eigen/Dense>
#include <chrono>


///
/// Solve linear system Ax=b by computing inv(A)*b
/// @param[in] A is an NxN Eigen::MatrixXd
/// @param[in] b is an N Eigen::VectorXd
/// @return    x which approximately solves Ax=b
///
//----------------SolveInverseBegin----------------
Eigen::VectorXd solve_inverse(Eigen::MatrixXd A, Eigen::VectorXd b) {
  Eigen::VectorXd solution;
  auto begin = std::chrono::high_resolution_clock::now();
  
  //// CMEA_START_TEMPLATE
  solution = A.inverse()*b; // Solve the system
  //// CMEA_END_TEMPLATE

  auto end = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout << "Time elapsed, inverse solver (ms): " << time/1e6 << std::endl;
  return solution;
}
//----------------SolveInverseEnd----------------


///
/// Solve linear system Ax=b using an LU solver
/// @param[in] A is an NxN Eigen::MatrixXd
/// @param[in] b is an N Eigen::VectorXd
/// @return    x which approximately solves Ax=b
///
//----------------SolveLUBegin----------------
Eigen::VectorXd solve_lu(Eigen::MatrixXd A, Eigen::VectorXd b) {
  Eigen::VectorXd solution;
  auto begin = std::chrono::high_resolution_clock::now();

  //// CMEA_START_TEMPLATE
  solution = A.lu().solve(b); // Solve the system
  //// CMEA_END_TEMPLATE

  auto end = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout << "Time elapsed,      LU solver (ms): " << time/1e6 << std::endl;
  return solution;
}
//----------------SolveLUEnd----------------


int main(int argc, char **argv){
  int N = 10;
  if (argc > 1)  {
    // Read N from command line
    // We use atof because we want to allow things like 1e7
    N = int(atof(argv[1]));
  }

  // Construct Hilbert's matrix
  Eigen::MatrixXd H(N,N);
  for(int i=0; i<N; i++) {
    for(int j=0; j<N; j++) {
      H(i,j) = 1./(1 + i + j);
    }
  }

  // Declare Eigen vector type for doubles
  using vector_t =  Eigen::VectorXd;

  // Create right hand side for the equation
  vector_t rhs(N);
  rhs.setOnes();

  //----------------UseInverseBegin----------------
  // Use solve_inverse to compute the solution, and find the residual
  //// CMEA_START_TEMPLATE
  auto xInv = solve_inverse(H,rhs);
  auto norm_rhs = rhs.norm();
  auto errorInv = (H*xInv - rhs).norm()/norm_rhs;
  std::cout << "Relative norm of Ax-b, inverse solver: " << errorInv << std::endl;
  //// CMEA_END_TEMPLATE
  //----------------UseInverseEnd----------------

  //----------------UseLUBegin----------------
  // Use solve_lu to compute the solution, and find the residual
  //// CMEA_START_TEMPLATE
  auto xLU = solve_lu(H, rhs);
  auto errorLU = (H*xLU - rhs).norm()/norm_rhs;
  std::cout << "Relative norm of Ax-b,      LU solver: " << errorLU << std::endl;
  //// CMEA_END_TEMPLATE
  //----------------UseLUEnd----------------

  return 0;
}
