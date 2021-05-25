
#ifndef BASE_HPP
#define BASE_HPP

#include <Eigen/Core>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat; // declares column-major
typedef Eigen::Triplet<double> T;

class Base {
public:

  Eigen::VectorXd R_x,R_y,R_z; // x,y,z component of polymer direction at time step n+1

  Eigen::VectorXd ArcLength;


  Eigen::VectorXd tension; // tension of the polymer at time step n+1
  



  SpMat M;           // matrix for setting tension_new, i.e. Mx tension_new = C
  

  // constructor
  Base(int,double,double,Eigen::Vector3d);
  ~Base();  

  void set_M();
  void update_M();

  virtual void init_R(double, Eigen::Vector3d& , int);
  
protected:
  const int N;                // number of grid spaces (N+1 grid points)
  const double L;             // arc length
  const double Delta_t;       // time discretisation
  const double Delta_s;       // arc length discretisation
  double alpha;         // ratio Delta_s^4/Delta_t
  
  Eigen::Vector3d x0;

  double FDiff_first(int ,Eigen::VectorXd&);
  double FDiff_second(int ,Eigen::VectorXd&);
  double FDiff_third(int ,Eigen::VectorXd&);
  double FDiff_fourth(int ,Eigen::VectorXd&);
  double BackwardDiff_first(int , Eigen::VectorXd&);
  double BackwardDiff_third(int ,Eigen::VectorXd&);
  double BackwardDiff_fourth(int ,Eigen::VectorXd&);



  virtual double set_gamma_0_plus(double );
  virtual double set_gamma_0_minus(double );
  virtual double set_s_i_plus(double ,int );
  virtual double set_s_i_minus(double ,int );
  virtual double set_d_i(double ,int );

private:
  double set_q_i(int );

  std::vector<T> init_M_coeffsmatrix();
  


  

};

#endif
