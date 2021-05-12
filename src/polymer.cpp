#include "polymer.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>
#include <cmath>
#include <random>

#define TENSION_OFFSET 1
#define POSITION_OFFSET 2

/* -------------------------------------------------------------------------- */
/* testing function! */
/* -------------------------------------------------------------------------- */

Eigen::VectorXd Polymer::s_i_plus_vec()
{

  Eigen::VectorXd vec(N);
  
  for (int i = 0; i < N; i++) {
    vec(i) = set_s_i_plus(i+TENSION_OFFSET);

  }
  return vec;

}


/* -------------------------------------------------------------------------- */
/* testing function! */
/* -------------------------------------------------------------------------- */
Eigen::VectorXd Polymer::s_i_minus_vec()
{

  Eigen::VectorXd vec(N);
  
  for (int i = 0; i < N; i++) {
    vec(i) = set_s_i_minus(i+TENSION_OFFSET);

  }
  return vec;

}


/* -------------------------------------------------------------------------- */
/* Setting the initial condition as 2d curve a la Seifert et al (PRL, 1996). */
/* -------------------------------------------------------------------------- */
void Polymer::init_R()
{

  Eigen::VectorXd qs(N);

  qs.setLinSpaced(M_PI/L,M_PI/Delta_s);

  std::mt19937 gen(2480);

  std::normal_distribution<double> distribution(0,1);



  Eigen::VectorXd a_qs(N);

  for (int j = 0; j < N; j++) {
    a_qs(j) = distribution(gen)/sqrt(kappa*L*qs(j)*qs(j)*qs(j)*qs(j));
  }


  R_y.setZero();
  R_x.setZero();
  R_z.setZero();
  
  for (int i = 1; i< N+1; i++) {
    R_y(i+POSITION_OFFSET) = (a_qs.array() *Eigen::sin((Delta_s*i*qs).array())).sum();
  }

  double s0 = 0.0;
  double tmp;
  for (int i = 1; i < N+1; i++) {
    tmp = FDiff_first(i+POSITION_OFFSET,R_y);
    s0 += sqrt(1-tmp*tmp/Delta_s/Delta_s);
    R_x(i+POSITION_OFFSET) = s0*Delta_s;
  }

  // zero second derivative BC at s = 0
  R_x(POSITION_OFFSET-1) = 2*R_x(POSITION_OFFSET)-R_x(POSITION_OFFSET+1);
  R_y(POSITION_OFFSET-1) = 2*R_y(POSITION_OFFSET)-R_y(POSITION_OFFSET+1);
  R_z(POSITION_OFFSET-1) = 2*R_z(POSITION_OFFSET)-R_z(POSITION_OFFSET+1);
  
  // zero second derivative BC at s = L
  R_x(N+POSITION_OFFSET+1) = 2*R_x(N+POSITION_OFFSET)-R_x(N+POSITION_OFFSET-1);
  R_y(N+POSITION_OFFSET+1) = 2*R_y(N+POSITION_OFFSET)-R_y(N+POSITION_OFFSET-1);
  R_z(N+POSITION_OFFSET+1) = 2*R_z(N+POSITION_OFFSET)-R_z(N+POSITION_OFFSET-1);


  // solve for free boundary conditions at s = 0
  
  Eigen::Vector3d sol;
  Eigen::Matrix3d mat;
  Eigen::Vector3d vec;

  double d1x = FDiff_first(POSITION_OFFSET,R_x)/Delta_s;
  double d1y = FDiff_first(POSITION_OFFSET,R_y)/Delta_s;
  double d1z = FDiff_first(POSITION_OFFSET,R_z)/Delta_s;

  tmp = kappa/(2*Delta_s*Delta_s*Delta_s);
  
  mat(0,0) = -tmp*(d1x*d1x-1); mat(0,1) = -tmp*d1x*d1y; mat(0,2) = -tmp*d1x*d1z;
  mat(1,0) = -tmp*d1x*d1y; mat(1,1) = -tmp*(d1y*d1y-1); mat(1,2) = -tmp*d1y*d1z;
  mat(2,0) = -tmp*d1x*d1z; mat(2,1) = -tmp*d1y*d1z; mat(2,2) = -tmp*(d1z*d1z-1);

  vec(0) = 4*Delta_s*d1x - R_x(POSITION_OFFSET+2);
  vec(1) = 4*Delta_s*d1y - R_y(POSITION_OFFSET+2);
  vec(2) = 4*Delta_s*d1z - R_z(POSITION_OFFSET+2);
  vec = -1*mat*vec;

  sol = mat.colPivHouseholderQr().solve(vec);
  
  R_x(POSITION_OFFSET-2) = sol(0);
  R_y(POSITION_OFFSET-2) = sol(1);
  R_z(POSITION_OFFSET-2) = sol(2);

  // solve for free boundary conditions at s = L

  d1x = FDiff_first(POSITION_OFFSET+N,R_x)/Delta_s;
  d1y = FDiff_first(POSITION_OFFSET+N,R_y)/Delta_s;
  d1z = FDiff_first(POSITION_OFFSET+N,R_z)/Delta_s;

  
  mat(0,0) = tmp*(d1x*d1x-1); mat(0,1) = tmp*d1x*d1y; mat(0,2) = tmp*d1x*d1z;
  mat(1,0) = tmp*d1x*d1y; mat(1,1) = tmp*(d1y*d1y-1); mat(1,2) = tmp*d1y*d1z;
  mat(2,0) = tmp*d1x*d1z; mat(2,1) = tmp*d1y*d1z; mat(2,2) = tmp*(d1z*d1z-1);

  vec(0) = 4*Delta_s*d1x + R_x(POSITION_OFFSET+N-2);
  vec(1) = 4*Delta_s*d1y + R_y(POSITION_OFFSET+N-2);
  vec(2) = 4*Delta_s*d1z + R_z(POSITION_OFFSET+N-2);
  vec = mat*vec;
  
  sol = mat.colPivHouseholderQr().solve(vec);

  R_x(POSITION_OFFSET+N+2) = sol(0);
  R_y(POSITION_OFFSET+N+2) = sol(1);
  R_z(POSITION_OFFSET+N+2) = sol(2);

  R_x.array() += x0(0);
  R_y.array() += x0(1);
  R_z.array() += x0(2);
  
  return;

}


/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
Polymer::Polymer(int Nin, double Lin, double Delta_tin, double kappain,Eigen::Vector3d x0in)
  : N(Nin), L(Lin), Delta_s(Lin/Nin), Delta_t(Delta_tin), kappa(kappain)
{
  

  alpha = Delta_s*Delta_s*Delta_s*Delta_s/Delta_t;


  // domain has N+1 grid points, but need an
  // extra 4 ghost points since 4th order PDE
  R_x.resize(N+5);
  R_y.resize(N+5);
  R_z.resize(N+5);

  A.resize(N+5,N+5);
  B_x.resize(N+5);
  B_y.resize(N+5);
  B_z.resize(N+5);

  x0 = x0in;

  
  M.resize(N+3,N+3);
  C.resize(N+3);
  // only need extra 2 ghost points because
  // tension is only second order PDE
  tension.resize(N+3);

  tension.setZero();

  R_x.setLinSpaced(-2*Delta_s,(L+2*Delta_s)/2.0);
  R_x = R_x.cwiseProduct(R_x);
  
  R_y.setZero();
  R_z.setZero();

  R_x.array() += x0(0);
  R_y.array() += x0(1);
  R_z.array() += x0(2);

}

/* -------------------------------------------------------------------------- */
/* Destructor */
/* -------------------------------------------------------------------------- */
Polymer::~Polymer()
{
}


/* -------------------------------------------------------------------------- */
/* Initialise matrix A in A * R_new = B. Only call at start of simulation. */
/* -------------------------------------------------------------------------- */
void Polymer::set_A()
{

  std::vector<T> coefficients = init_A_coeffsmatrix();

  A.setFromTriplets(coefficients.begin(),coefficients.end());
  
}

/* -------------------------------------------------------------------------- */
/* Helper function to set up matrix A. */
/* -------------------------------------------------------------------------- */
std::vector<T>  Polymer::init_A_coeffsmatrix()
{
  std::vector<T> coeffs;
  
  // set boundary conditions
  coeffs.push_back(T(0,0,kappa));
  coeffs.push_back(T(0,1,-2*kappa-Delta_s*Delta_s*tension(TENSION_OFFSET)));
  coeffs.push_back(T(0,3,2*kappa+Delta_s*Delta_s*tension(TENSION_OFFSET)));
  coeffs.push_back(T(0,4,-kappa));
  coeffs.push_back(T(1,1,kappa));
  coeffs.push_back(T(1,2,-2*kappa));
  coeffs.push_back(T(1,3,kappa));
  coeffs.push_back(T(N+2,N+2,1));
  coeffs.push_back(T(N+3,N+1,kappa));
  coeffs.push_back(T(N+3,N+2,-2*kappa));
  coeffs.push_back(T(N+3,N+3,kappa));
  coeffs.push_back(T(N+4,N,kappa));
  coeffs.push_back(T(N+4,N+1,set_s_i_minus(N+TENSION_OFFSET)));
  coeffs.push_back(T(N+4,N+2,set_d_i(N+TENSION_OFFSET)));
  coeffs.push_back(T(N+4,N+3,set_s_i_plus(N+TENSION_OFFSET)));
  coeffs.push_back(T(N+4,N+4,kappa));

  
  for (int i = POSITION_OFFSET; i < N+POSITION_OFFSET; i++) {

    coeffs.push_back(T(i,i-2,kappa));
    coeffs.push_back(T(i,i-1,set_s_i_minus(i-POSITION_OFFSET+TENSION_OFFSET)));
    coeffs.push_back(T(i,i,set_d_i(i-POSITION_OFFSET+TENSION_OFFSET)));
    coeffs.push_back(T(i,i+1,set_s_i_plus(i-POSITION_OFFSET+TENSION_OFFSET)));
    coeffs.push_back(T(i,i+2,kappa));
  }

  return coeffs;
	
}

/* -------------------------------------------------------------------------- */
/* Update matrix A. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void Polymer::update_A()
{
  // update first 5 cols using slow method
  A.coeffRef(0,1) = -2*kappa-Delta_s*Delta_s*tension(TENSION_OFFSET);
  A.coeffRef(0,3) = -2*kappa+Delta_s*Delta_s*tension(TENSION_OFFSET);

  for (int i = 2; i < 6; i++) {

    A.coeffRef(i,i-1) = set_s_i_minus(i-POSITION_OFFSET+TENSION_OFFSET);
    if (i < 5) 
      A.coeffRef(i,i) = set_d_i(i-POSITION_OFFSET+TENSION_OFFSET);
    if (i < 4) 
      A.coeffRef(i,i+1) = set_s_i_plus(i-POSITION_OFFSET+TENSION_OFFSET);
  }


  // update last 4 cols using slow method
  int i = N;
  A.coeffRef(i,i+1) = set_s_i_plus(i-POSITION_OFFSET+TENSION_OFFSET);

  i = N+1;
  A.coeffRef(i,i) = set_d_i(i-POSITION_OFFSET+TENSION_OFFSET);
  A.coeffRef(i,i+1) = set_s_i_plus(i-POSITION_OFFSET+TENSION_OFFSET);
  
  A.coeffRef(N+4,N+1) = set_s_i_minus(N+TENSION_OFFSET);
  A.coeffRef(N+4,N+2) = set_d_i(N+TENSION_OFFSET);
  A.coeffRef(N+4,N+3) = set_s_i_plus(N+TENSION_OFFSET);

  // update middle columns
  for (int k = 5; k < N+1; k++) {
    int count = 0;
    for (SpMat::InnerIterator it(A,k); it; ++it) {


      if (count == 1) {

      	it.valueRef() = set_s_i_plus(k-1-POSITION_OFFSET+TENSION_OFFSET);
      } else if (count == 2) {
      	it.valueRef() = set_d_i(k-POSITION_OFFSET+TENSION_OFFSET);
      } else if (count == 3) {
      	it.valueRef() = set_s_i_minus(k+1-POSITION_OFFSET+TENSION_OFFSET);
      }

      count += 1;
    }

  }
  return;
}





/* -------------------------------------------------------------------------- */
/* Update the three vectors B_x, B_y, B_z in A * R_new_i = B_i.
/* -------------------------------------------------------------------------- */
void Polymer::update_B(Eigen::Vector3d& F)
{

  B_x = alpha*R_x;
  B_x(0) = 2*Delta_s*Delta_s*Delta_s*F(0);
  B_x(1) = 0.0;
  B_x(N+4) = B_x(N+2);
  B_x(N+3) = 0.0;
  B_x(N+2) = x0(0);


  B_y = alpha*R_y;
  B_y(0) = 2*Delta_s*Delta_s*Delta_s*F(1);
  B_y(1) = 0.0;
  B_y(N+4) = B_y(N+2);
  B_y(N+3) = 0.0;
  B_y(N+2) = x0(1);

  B_z = alpha*R_z;
  B_z(0) = 2*Delta_s*Delta_s*Delta_s*F(2);
  B_z(1) = 0.0;
  B_z(N+4) = B_z(N+2);
  B_z(N+3) = 0.0;
  B_z(N+2) = x0(2);

  return;
  
  
}



/* -------------------------------------------------------------------------- */
/* Initialise matrix M in M * tension = C. Only call at start of simulation. */
/* -------------------------------------------------------------------------- */
void Polymer::set_M()
{

  std::vector<T> coefficients = init_M_coeffsmatrix();

  M.setFromTriplets(coefficients.begin(),coefficients.end());

  return;

}

/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix M. */
/* -------------------------------------------------------------------------- */
std::vector<T> Polymer::init_M_coeffsmatrix()
{
  std::vector<T> coeffs;

  double ds4 = Delta_s*Delta_s*Delta_s*Delta_s;



  // set boundary conditions
  coeffs.push_back(T(0,0,ds4));
  coeffs.push_back(T(0,1,set_q_i(TENSION_OFFSET)));
  coeffs.push_back(T(0,2,ds4));
  coeffs.push_back(T(1,1,ds4));

  coeffs.push_back(T(N+2,N,-ds4));
  coeffs.push_back(T(N+2,N+2,ds4));

  

  
  for (int i = TENSION_OFFSET+1; i < N+TENSION_OFFSET+1; i++) {

    coeffs.push_back(T(i,i-1,ds4));
    coeffs.push_back(T(i,i,set_q_i(i-TENSION_OFFSET)));
    coeffs.push_back(T(i,i+1,ds4));

  }

  return coeffs;

  
}


/* -------------------------------------------------------------------------- */
/* Update matrix M. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void Polymer::update_M()
{
  // update first three cols by hand
  M.coeffRef(0,1) = set_q_i(TENSION_OFFSET);

  M.coeffRef(2,2) = set_q_i(TENSION_OFFSET+1);

  // update middle columns
  for (int k = 3; k < N+2; k++) {
    int count = 0;
    for (SpMat::InnerIterator it(M,k); it; ++it) {
      if (count == 1) {
      	it.valueRef() = set_q_i(k-TENSION_OFFSET);
      }
      count += 1;
    }
  }
}

/* -------------------------------------------------------------------------- */
/* Update the vectors C in M * tension = C.
/* -------------------------------------------------------------------------- */
void Polymer::update_C(Eigen::Vector3d& F)
{



  double d1x = FDiff_first(POSITION_OFFSET,R_x);
  double d1y = FDiff_first(POSITION_OFFSET,R_y);
  double d1z = FDiff_first(POSITION_OFFSET,R_z);

  double dd2x = FDiff_second(POSITION_OFFSET,R_x);
  double dd2y = FDiff_second(POSITION_OFFSET,R_y);
  double dd2z = FDiff_second(POSITION_OFFSET,R_z);

  double ddd3x = FDiff_third(POSITION_OFFSET,R_x);
  double ddd3y = FDiff_third(POSITION_OFFSET,R_y);
  double ddd3z = FDiff_third(POSITION_OFFSET,R_z);

  double dddd4x = FDiff_fourth(POSITION_OFFSET,R_x);
  double dddd4y = FDiff_fourth(POSITION_OFFSET,R_y);
  double dddd4z = FDiff_fourth(POSITION_OFFSET,R_z);

  C(0) = -kappa*(4*(dd2x*dddd4x + dd2y*dddd4y + dd2z*dddd4z)
		 + 3*(ddd3x*ddd3x+ddd3y*ddd3y+ddd3z*ddd3z));

  C(1) = (kappa*(d1x*ddd3x + d1y*ddd3y + d1z*ddd3z)
	  + Delta_s*Delta_s*Delta_s*(d1x*F(0) + d1y*F(1) + d1z*F(2)));

  for (int i = TENSION_OFFSET+1; i< N+TENSION_OFFSET+1; i++) {
    C(i) = set_c_i(i-TENSION_OFFSET+POSITION_OFFSET);
  }

  d1x = FDiff_first(N+POSITION_OFFSET,R_x);
  d1y = FDiff_first(N+POSITION_OFFSET,R_y);
  d1z = FDiff_first(N+POSITION_OFFSET,R_z);
  dddd4x = FDiff_fourth(N+POSITION_OFFSET,R_x);
  dddd4y = FDiff_fourth(N+POSITION_OFFSET,R_y);
  dddd4z = FDiff_fourth(N+POSITION_OFFSET,R_z);

  
  C(N+2) = kappa*(d1x*dddd4x + d1y*dddd4y + d1z*dddd4z);

  return;
  
}


double Polymer::set_c_i(int i)
{

  double dd2x = FDiff_second(i,R_x);
  double dd2y = FDiff_second(i,R_y);
  double dd2z = FDiff_second(i,R_z);

  double ddd3x = FDiff_third(i,R_x);
  double ddd3y = FDiff_third(i,R_y);
  double ddd3z = FDiff_third(i,R_z);

  double dddd4x = FDiff_fourth(i,R_x);
  double dddd4y = FDiff_fourth(i,R_y);
  double dddd4z = FDiff_fourth(i,R_z);


  return -kappa*(4*(dd2x*dddd4x + dd2y*dddd4y + dd2z*dddd4z)
		 + 3*(ddd3x*ddd3x+ddd3y*ddd3y+ddd3z*ddd3z));
  
}

/* -------------------------------------------------------------------------- */
double Polymer::set_s_i_plus(int i)
{
  return -4*kappa - Delta_s*Delta_s/4*(4*tension(i)
				       + tension(i+1) - tension(i-1));
}

/* -------------------------------------------------------------------------- */
double Polymer::set_s_i_minus(int i)
{
  return -4*kappa - Delta_s*Delta_s/4*(4*tension(i)
				       - tension(i+1) + tension(i-1));
}


/* -------------------------------------------------------------------------- */
double Polymer::set_d_i(int i)
{
  return alpha + 6*kappa + 2*Delta_s*Delta_s*tension(i);
}

/* -------------------------------------------------------------------------- */
double Polymer::set_q_i(int i)
{
  double drx = FDiff_second(i,R_x);
  double dry = FDiff_second(i,R_y);
  double drz = FDiff_second(i,R_z);
  return -Delta_s*Delta_s*(drx*drx + dry*dry + drz*drz
			   + 2*Delta_s*Delta_s);
}



/* -------------------------------------------------------------------------- */
/* First order finite difference */
/* -------------------------------------------------------------------------- */
double Polymer::FDiff_first(int i,Eigen::VectorXd& vec)
{
  return (vec(i+1)-vec(i-1))/2.0;
  
}
/* -------------------------------------------------------------------------- */
/* Second order finite difference */
/* -------------------------------------------------------------------------- */
double Polymer::FDiff_second(int i,Eigen::VectorXd& vec)
{
  return (vec(i+1)+vec(i-1)-2*vec(i));
}
/* -------------------------------------------------------------------------- */
/* Third order finite difference */
/* -------------------------------------------------------------------------- */
double Polymer::FDiff_third(int i,Eigen::VectorXd& vec)
{
  return (vec(i+2)-vec(i-2)-2*vec(i+1)+2*vec(i-1))/2.0;
}

/* -------------------------------------------------------------------------- */
/* Fourth order finite difference */
/* -------------------------------------------------------------------------- */
double Polymer::FDiff_fourth(int i,Eigen::VectorXd& vec)
{
  return vec(i+2)+vec(i-2)-4*vec(i+1)-4*vec(i-1)+6*vec(i);
}
