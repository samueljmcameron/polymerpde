
#include "PolyConfig.h"
#include "polymer.hpp"
#include <Eigen/SparseLU>

#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{

  int N = 10;
  std::cout << argv[0] << " version " << poly_VERSION_MAJOR << "."
	    << poly_VERSION_MINOR << std::endl;

  Eigen::Vector3d x0(1,1,1);

  Eigen::Vector3d F(1000.0,0.0,0.0);

  Polymer pmer(N,2.0,0.001,100.0,x0);


  pmer.init_R();
  
  std::ofstream file("initialconditions.txt");
  if (file.is_open()) {
    Eigen::MatrixXd m(N+5,3);

    m << pmer.R_x, pmer.R_y, pmer.R_z;

    file << m << std::endl;
  }


  pmer.set_M();
  pmer.update_C(F);
  
  pmer.M.makeCompressed();

  Eigen::SparseLU<SpMat> tensionsolver;
  tensionsolver.analyzePattern(pmer.M);
  
  tensionsolver.factorize(pmer.M);
  if (tensionsolver.info() != Eigen::Success ) {
    std::cout << "Failed at factorize step!" << std::endl;
    return -1;
  }
  pmer.tension = tensionsolver.solve(pmer.C);
  if (tensionsolver.info() != Eigen::Success ) {
    std::cout << "Failed at solve step! " << std::endl;
    return -1;
  }
  
  std::ofstream tensionfile("initialtension.txt");
  if (tensionfile.is_open()) {
    tensionfile << pmer.tension << std::endl;
  }

  std::cout << pmer.tension(1) << std::endl;
  pmer.set_A();
  std::cout << pmer.A << std::endl;
  pmer.update_B(F);
  //std::cout << pmer.B_x << std::endl;
  //std::cout << pmer.B_y << std::endl;
  Eigen::SparseLU<SpMat> positionsolver;

  pmer.A.makeCompressed();
  positionsolver.analyzePattern(pmer.A);
  positionsolver.factorize(pmer.A);
  if (positionsolver.info() != Eigen::Success ) {
    std::cout << "Failed at factorize step!" << std::endl;
    return -1;
  }
  pmer.R_x = positionsolver.solve(pmer.B_x);
  if (positionsolver.info() != Eigen::Success ) {
    std::cout << "Failed at R_x solve step! " << std::endl;
    return -1;
  }

  pmer.R_y = positionsolver.solve(pmer.B_y);
  if (positionsolver.info() != Eigen::Success ) {
    std::cout << "Failed at R_y solve step! " << std::endl;
    return -1;
  }

  std::ofstream secondfile("singletimestep.txt");
  if (file.is_open()) {
    Eigen::MatrixXd m(N+5,3);

    m << pmer.R_x, pmer.R_y, pmer.R_z;

    secondfile << m << std::endl;
  }

  
  return 0;

}
