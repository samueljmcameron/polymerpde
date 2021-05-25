
#include "PolyConfig.h"
#include "polymer.hpp"
#include <Eigen/SparseLU>

#include <iostream>
#include <fstream>
#include <string>

#include <chrono>

int update_polymer(Polymer &pmer,Eigen::Vector3d &F,
		   Eigen::SparseLU<SpMat> &tensionsolver,
		   Eigen::SparseLU<SpMat> &positionsolver);

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

  Eigen::Vector3d x0(1,1,1);

  Eigen::Vector3d F(Fx,Fy,Fz);

  Polymer pmer(N,L,dt,kappa,x0);


  pmer.init_R(F);
  
  std::ofstream file("initialconditions.txt");
  if (file.is_open()) {
    Eigen::MatrixXd m(pmer.R_x.size(),4);

    m << pmer.ArcLength, pmer.R_x, pmer.R_y, pmer.R_z;

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

  for (int i = 0; i < 1000; i++) {
  
    update_polymer(pmer,F,tensionsolver,positionsolver);
  }


  std::ofstream firstfile("firsttimestep.txt");
  if (firstfile.is_open()) {
    Eigen::MatrixXd m(pmer.R_x.size(),3);

    m << pmer.R_x, pmer.R_y, pmer.R_z;

    firstfile << m << std::endl;
  }
  auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < 9000; i++) {
  
    update_polymer(pmer,F,tensionsolver,positionsolver);
  }

  auto stop = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);

  std::cout << duration.count() << std::endl;
  std::ofstream secondfile("secondtimestep.txt");
  if (secondfile.is_open()) {
    Eigen::MatrixXd m(pmer.R_x.size(),3);

    m << pmer.R_x, pmer.R_y, pmer.R_z;

    secondfile << m << std::endl;
  }


  
  return 0;

}


/* -------------------------------------------------------------------------- */
/* update polymer one step */
/* -------------------------------------------------------------------------- */
int update_polymer(Polymer &pmer,Eigen::Vector3d &F,
		   Eigen::SparseLU<SpMat> &tensionsolver,
		   Eigen::SparseLU<SpMat> &positionsolver)
{

  pmer.update_M();

  pmer.update_C(F);
  
  //  pmer.M.makeCompressed();

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
  

  pmer.update_A();


  pmer.update_B(F);


  //  pmer.A.makeCompressed();
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
  return 0;
}
