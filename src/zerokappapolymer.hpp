
#ifndef ZEROKAPPAPOLYMER_HPP
#define ZEROKAPPAPOLYMER_HPP

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "base.hpp"

typedef Eigen::SparseMatrix<double> SpMat; // declares column-major
typedef Eigen::Triplet<double> T;

class ZeroKappaPolymer : public Base {
public:
  

  SpMat A; // matrix for setting R_new, i.e. A R_new = B
  Eigen::VectorXd B_x,B_y,B_z; // vector for setting R_new, i.e. A R_new = B


  Eigen::VectorXd C; // vector for setting tension_new, i.e. Mx tension_new = C
  

  // constructor
  ZeroKappaPolymer(int,double,double,Eigen::Vector3d);
  ~ZeroKappaPolymer();  
  Eigen::VectorXd s_i_plus_vec();
  Eigen::VectorXd s_i_minus_vec();
  
  void set_A();
  void update_A();
  void update_B(Eigen::Vector3d& );

  void update_C(Eigen::Vector3d& );

  void init_R(double,Eigen::Vector3d& );
  
private:

  double set_gamma_0_plus();
  double set_gamma_0_minus();
  double set_d_i(int );    
  double set_s_i_plus(int );
  double set_s_i_minus(int );


  
  std::vector<T> init_A_coeffsmatrix();



  

};

#endif
