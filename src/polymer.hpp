
#ifndef POLYMER_HPP
#define POLYMER_HPP

#include <Eigen/Core>
#include <Eigen/Sparse>


typedef Eigen::SparseMatrix<double> SpMat; // declares column-major
typedef Eigen::Triplet<double> T;

class Polymer {
public:

  Eigen::VectorXd R_x,R_y,R_z; // x,y,z component of polymer direction at time step n+1


  Eigen::VectorXd tension; // tension of the polymer at time step n+1
  

  SpMat A; // matrix for setting R_new, i.e. A R_new = B
  Eigen::VectorXd B_x,B_y,B_z; // vector for setting R_new, i.e. A R_new = B


  SpMat M;           // matrix for setting tension_new, i.e. Mx tension_new = C
  Eigen::VectorXd C; // vector for setting tension_new, i.e. Mx tension_new = C
  

  // constructor
  Polymer(int,double,double,double,Eigen::Vector3d);
  ~Polymer();  
  Eigen::VectorXd s_i_plus_vec();
  Eigen::VectorXd s_i_minus_vec();
  
  void set_A();
  void update_A();
  void update_B(Eigen::Vector3d& );

  void set_M();
  void update_M();
  void update_C(Eigen::Vector3d& );

  void init_R();
  double kappa;         // bending modulus
  
private:
  const int N;                // number of grid spaces (N+1 grid points)
  const double L;             // arc length
  const double Delta_t;       // time discretisation
  const double Delta_s;       // arc length discretisation
  
  double alpha;         // ratio Delta_s^4/Delta_t
  double set_s_i_plus(int ); // 
  double set_s_i_minus(int );
  double set_d_i(int );
  std::vector<T> init_A_coeffsmatrix();

  Eigen::Vector3d x0;
  double set_q_i(int );

  std::vector<T> init_M_coeffsmatrix();


  double set_c_i(int );

  
  double FDiff_first(int ,Eigen::VectorXd&);
  double FDiff_second(int ,Eigen::VectorXd&);
  double FDiff_third(int ,Eigen::VectorXd&);
  double FDiff_fourth(int ,Eigen::VectorXd&);


  

};

#endif
