#include "zerokappapolymer.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>
#include <cmath>

#define POSITION_OFFSET 1
#define GHOST_POINTS 1

/* -------------------------------------------------------------------------- */
/* testing function! */
/* -------------------------------------------------------------------------- */

Eigen::VectorXd ZeroKappaPolymer::s_i_plus_vec()
{

  Eigen::VectorXd vec(N);
  
  for (int i = 0; i < N; i++) {
    vec(i) = set_s_i_plus(i);

  }
  return vec;

}


/* -------------------------------------------------------------------------- */
/* testing function! */
/* -------------------------------------------------------------------------- */
Eigen::VectorXd ZeroKappaPolymer::s_i_minus_vec()
{

  Eigen::VectorXd vec(N);
  
  for (int i = 0; i < N; i++) {
    vec(i) = set_s_i_minus(i);

  }
  return vec;

}


/* -------------------------------------------------------------------------- */
/* Setting the initial condition as 2d curve a la Seifert et al (PRL, 1996). */
/* -------------------------------------------------------------------------- */
void ZeroKappaPolymer::init_R(double kappa,Eigen::Vector3d& F)
{


  Base::init_R(kappa,F,POSITION_OFFSET);
  
  // solve for free boundary conditions at s = 0

  double fmag = F.norm();

  R_x(POSITION_OFFSET-1) = R_x(POSITION_OFFSET+1) - 2*Delta_s*F(0)/fmag;
  R_y(POSITION_OFFSET-1) = R_y(POSITION_OFFSET+1) - 2*Delta_s*F(1)/fmag;
  R_z(POSITION_OFFSET-1) = R_z(POSITION_OFFSET+1) - 2*Delta_s*F(2)/fmag;
  
  return;

}


/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
ZeroKappaPolymer::ZeroKappaPolymer(int Nin, double Lin, double Delta_tin,
				   Eigen::Vector3d x0in)
  : Base(Nin,Lin,Delta_tin,x0in)
{

  // domain has N+1 grid points, but need an
  // extra 4 ghost points since 4th order PDE
  R_x.resize(N+1+GHOST_POINTS);
  R_y.resize(N+1+GHOST_POINTS);
  R_z.resize(N+1+GHOST_POINTS);

  A.resize(N+1+GHOST_POINTS,N+1+GHOST_POINTS);
  B_x.resize(N+1+GHOST_POINTS);
  B_y.resize(N+1+GHOST_POINTS);
  B_z.resize(N+1+GHOST_POINTS);

  x0 = x0in;

  
  C.resize(N+2);
  // only need extra 2 ghost points because
  // tension is only second order PDE

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
ZeroKappaPolymer::~ZeroKappaPolymer()
{
}


/* -------------------------------------------------------------------------- */
/* Initialise matrix A in A * R_new = B. Only call at start of simulation. */
/* -------------------------------------------------------------------------- */
void ZeroKappaPolymer::set_A()
{

  std::vector<T> coefficients = init_A_coeffsmatrix();

  A.setFromTriplets(coefficients.begin(),coefficients.end());
  
}

/* -------------------------------------------------------------------------- */
/* Helper function to set up matrix A. */
/* -------------------------------------------------------------------------- */
std::vector<T>  ZeroKappaPolymer::init_A_coeffsmatrix()
{
  std::vector<T> coeffs;
  
  // set boundary conditions
  coeffs.push_back(T(-1+POSITION_OFFSET,-1+POSITION_OFFSET,-tension(0)));
  coeffs.push_back(T(-1+POSITION_OFFSET,1+POSITION_OFFSET,tension(0)));
  coeffs.push_back(T(N+POSITION_OFFSET,N+POSITION_OFFSET,1));

  int k = POSITION_OFFSET;

  coeffs.push_back(T(k,k-1,set_gamma_0_minus()));
  coeffs.push_back(T(k,k,set_d_i(0)));
  coeffs.push_back(T(k,k+1, set_gamma_0_plus()));

  
  for (int i = 1; i < N; i++) {
    k = i + POSITION_OFFSET;
    
    coeffs.push_back(T(k,k-1,
		       set_s_i_minus(i)));
    coeffs.push_back(T(k,k,
		       set_d_i(i)));
    coeffs.push_back(T(k,k+1,
		       set_s_i_plus(i)));
  }

  return coeffs;
	
}

/* -------------------------------------------------------------------------- */
/* Update matrix A. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void ZeroKappaPolymer::update_A()
{

  int k;
  
  // update first 3 cols using slow method
  A.coeffRef(-1+POSITION_OFFSET,-1+POSITION_OFFSET) = -tension(0);
  A.coeffRef(-1+POSITION_OFFSET,1+POSITION_OFFSET) = tension(0);

  k = POSITION_OFFSET;
  
  A.coeffRef(k,k) = set_d_i(0);
  A.coeffRef(k,k-1) = set_gamma_0_minus();  
  A.coeffRef(k,k+1) = set_gamma_0_plus();


  
  for (int i = 1; i < 3; i++) {
    k = i + POSITION_OFFSET;
    A.coeffRef(k,k-1) = set_s_i_minus(i);
    if (i < 2) 
      A.coeffRef(k,k) = set_d_i(i);

  }
  
  // update last two cols using slow method

  for (int i = N-2; i < N; i++) {
    k = i + POSITION_OFFSET;
    A.coeffRef(k,k+1) = set_s_i_plus(i);
    if (i == N-1) {
      A.coeffRef(k,k) = set_d_i(i);
      A.coeffRef(k,k+1) = set_s_i_plus(i);
    }
  }



  // update middle columns
  for (int i = 2; i < N-1; i++) {
    int count = 0;
    for (SpMat::InnerIterator it(A,i+POSITION_OFFSET); it; ++it) {

      if (count == 0) {

      	it.valueRef() = set_s_i_plus(i-1);
      } else if (count == 1) {
      	it.valueRef() = set_d_i(i);
      } else if (count == 2) {
      	it.valueRef() = set_s_i_minus(i+1);
      }

      count += 1;
    }

  }
  return;
}





/* -------------------------------------------------------------------------- */
/* Update the three vectors B_x, B_y, B_z in A * R_new_i = B_i.
/* -------------------------------------------------------------------------- */
void ZeroKappaPolymer::update_B(Eigen::Vector3d& F)
{

  double ds3 = Delta_s*Delta_s*Delta_s;
  B_x = alpha*R_x;
  B_x(-1+POSITION_OFFSET) = 2*ds3*F(0);

  B_x(N+POSITION_OFFSET) = x0(0);


  B_y = alpha*R_y;
  B_y(-1+POSITION_OFFSET) = 2*ds3*F(1);
  B_y(N+POSITION_OFFSET) = x0(1);

  B_z = alpha*R_z;
  B_z(-1+POSITION_OFFSET) = 2*ds3*F(2);
  B_z(N+POSITION_OFFSET) = x0(2);

  return;
  
  
}



/* -------------------------------------------------------------------------- */
/* Update the vectors C in M * tension = C.
/* -------------------------------------------------------------------------- */
void ZeroKappaPolymer::update_C(Eigen::Vector3d& F)
{

  C.setZero();
  double d1x = FDiff_first(POSITION_OFFSET,R_x);
  double d1y = FDiff_first(POSITION_OFFSET,R_y);
  double d1z = FDiff_first(POSITION_OFFSET,R_z);
  C(0) = (Delta_s*Delta_s*Delta_s*(d1x*F(0) + d1y*F(1) + d1z*F(2)));


  return;
}


/* -------------------------------------------------------------------------- */
double ZeroKappaPolymer::set_gamma_0_plus()
{
  return Base::set_gamma_0_plus(0.0);
}

/* -------------------------------------------------------------------------- */
double ZeroKappaPolymer::set_gamma_0_minus()
{
  return Base::set_gamma_0_minus(0.0);
}



/* -------------------------------------------------------------------------- */
double ZeroKappaPolymer::set_s_i_plus(int i)
{
  return Base::set_s_i_plus(0.0,i);
}

/* -------------------------------------------------------------------------- */
double ZeroKappaPolymer::set_s_i_minus(int i)
{
  return Base::set_s_i_minus(0.0,i);
}


/* -------------------------------------------------------------------------- */
double ZeroKappaPolymer::set_d_i(int i)
{
  return Base::set_d_i(0.0,i);
}
