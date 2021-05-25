
#include "PolyConfig.h"
#include "zerokappapolymer.hpp"
#include <Eigen/SparseLU>

#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[])
{ 
  std::cout << argv[0] << " version " << poly_VERSION_MAJOR << "."
	    << poly_VERSION_MINOR << std::endl;
  int N = std::stoi(argv[1]);
  double L = std::stod(argv[2]);
  double kappa = std::stod(argv[3]);
  double dt = std::stod(argv[4]);
  double Fx = std::stod(argv[5]);
  double Fy = std::stod(argv[6]);
  double Fz = std::stod(argv[7]);

		     

  std::cout << dt << " " << N << std::endl;
  
  Eigen::Vector3d x0(1,1,1);

  Eigen::Vector3d F(1000.0,0.0,0.0);

  ZeroKappaPolymer pmer(N,L,dt,x0);


  pmer.init_R(kappa,F);
  
  std::ofstream file("initialconditions.txt");
  if (file.is_open()) {
    Eigen::MatrixXd m(N+2,3);

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


  pmer.set_A();

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

  pmer.R_z = positionsolver.solve(pmer.B_z);
  if (positionsolver.info() != Eigen::Success ) {
    std::cout << "Failed at R_z solve step! " << std::endl;
    return -1;
  }

  std::ofstream secondfile("singletimestep.txt");
  if (file.is_open()) {
    Eigen::MatrixXd m(N+2,3);

    m << pmer.R_x, pmer.R_y, pmer.R_z;

    secondfile << m << std::endl;
  }
  
  return 0;

}
