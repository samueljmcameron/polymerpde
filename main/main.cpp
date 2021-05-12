
#include "PolyConfig.h"
#include "polymer.hpp"

#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{

  int N = 1000;
  std::cout << argv[0] << " version " << poly_VERSION_MAJOR << "."
	    << poly_VERSION_MINOR << std::endl;

  Eigen::Vector3d x0(1,1,1);

  Eigen::Vector3d F(1.0,0.0,0.0);

  Polymer pmer(N,250.0,0.01,100.0,x0);


  pmer.init_R();
  
  std::ofstream file("initialconditions.txt");
  if (file.is_open()) {
    Eigen::MatrixXd m(N+5,3);

    m << pmer.R_x, pmer.R_y, pmer.R_z;

    file << m << std::endl;
  }

  pmer.update_C(F);

  
  std::ofstream tensionfile("initialtension.txt");
  if (tensionfile.is_open()) {
    tensionfile << pmer.tension << std::endl;
  }
  
  return 0;

}
