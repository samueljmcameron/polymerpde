#include "polymer.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>
#include <cmath>

#define POSITION_OFFSET 2
#define GHOST_POINTS 3

/* -------------------------------------------------------------------------- */
/* testing function! */
/* -------------------------------------------------------------------------- */

Eigen::VectorXd Polymer::s_i_plus_vec()
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
Eigen::VectorXd Polymer::s_i_minus_vec()
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
void Polymer::init_R(Eigen::Vector3d& F)
{


  ArcLength.setLinSpaced(-2*Delta_s,(L+Delta_s));

  R_x = -0.2*ArcLength;
  R_y = -0.3*ArcLength;
  

  R_z = sqrt(1-0.2*0.2-0.3*0.3)*ArcLength;
  
  R_x.array() -= R_x(N+POSITION_OFFSET) - x0(0);
  R_y.array() -= R_y(N+POSITION_OFFSET) - x0(1);
  R_z.array() -= R_z(N+POSITION_OFFSET) - x0(2);

  
  //  Base::init_R(kappa,F,POSITION_OFFSET);
  // zero second derivative BC at s = 0
  R_x(POSITION_OFFSET-1) = 2*R_x(POSITION_OFFSET)-R_x(POSITION_OFFSET+1);
  R_y(POSITION_OFFSET-1) = 2*R_y(POSITION_OFFSET)-R_y(POSITION_OFFSET+1);
  R_z(POSITION_OFFSET-1) = 2*R_z(POSITION_OFFSET)-R_z(POSITION_OFFSET+1);
  
  // zero second derivative BC at s = L
  R_x(N+POSITION_OFFSET+1) = 2*R_x(N+POSITION_OFFSET)-R_x(N+POSITION_OFFSET-1);
  R_y(N+POSITION_OFFSET+1) = 2*R_y(N+POSITION_OFFSET)-R_y(N+POSITION_OFFSET-1);
  R_z(N+POSITION_OFFSET+1) = 2*R_z(N+POSITION_OFFSET)-R_z(N+POSITION_OFFSET-1);


  // solve for free boundary conditions at s = 0

  double tmp = 2*Delta_s*Delta_s*Delta_s/kappa;

  Eigen::Vector3d R_0;
  R_0(0) = R_x(POSITION_OFFSET);
  R_0(1) = R_y(POSITION_OFFSET);
  R_0(2) = R_z(POSITION_OFFSET);

  Eigen::Vector3d R_2;
  R_2(0) = R_x(POSITION_OFFSET+2);
  R_2(1) = R_y(POSITION_OFFSET+2);
  R_2(2) = R_z(POSITION_OFFSET+2);

  Eigen::Vector3d ds_R;
  ds_R(0) = FDiff_first(POSITION_OFFSET,R_x)/Delta_s;
  ds_R(1) = FDiff_first(POSITION_OFFSET,R_y)/Delta_s;
  ds_R(2) = FDiff_first(POSITION_OFFSET,R_z)/Delta_s;


  Eigen::Vector3d d;

  d = R_2 - 4*Delta_s*ds_R + tmp*F;

  

  double c = R_0.dot(R_0)-2*R_0.dot(d)+d.dot(d) - 4*(Delta_s*Delta_s);
  double b = 2*tmp*ds_R.dot(R_0-d);
  double a = tmp*tmp*ds_R.dot(ds_R);


  double u = (-b + sqrt(b*b - 4*a*c))/(2*a);

  std::cout << " u = " << u << std::endl;

  R_x(POSITION_OFFSET-2) = d(0) - tmp*u*ds_R(0);
  R_y(POSITION_OFFSET-2) = d(1) - tmp*u*ds_R(1);
  R_z(POSITION_OFFSET-2) = d(2) - tmp*u*ds_R(2);

  tmp = (kappa*FDiff_third(POSITION_OFFSET,R_x)/(Delta_s*Delta_s*Delta_s)
	 +F(0) - u * FDiff_first(POSITION_OFFSET,R_x)/Delta_s);

  std:: cout << "x component of output with + root = " << tmp << std::endl;

  tmp = (kappa*FDiff_third(POSITION_OFFSET,R_y)/(Delta_s*Delta_s*Delta_s)
	 +F(1) - u * FDiff_first(POSITION_OFFSET,R_y)/Delta_s);

  std:: cout << "y component of output with + root = " << tmp << std::endl;

  tmp = (kappa*FDiff_third(POSITION_OFFSET,R_z)/(Delta_s*Delta_s*Delta_s)
	 +F(2) - u * FDiff_first(POSITION_OFFSET,R_z)/Delta_s);

  std:: cout << "z component of output with + root = " << tmp << std::endl;
  

  
  return;

}


/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
Polymer::Polymer(int Nin, double Lin, double Delta_tin, double kappain,Eigen::Vector3d x0in)
  : Base(Nin,Lin,Delta_tin,x0in), kappa(kappain)
{
  


  // domain has N+1 grid points, but need an
  // extra 4 ghost points since 4th order PDE
  R_x.resize(N+1+GHOST_POINTS);
  R_y.resize(N+1+GHOST_POINTS);
  R_z.resize(N+1+GHOST_POINTS);

  ArcLength.resize(N+1+GHOST_POINTS);
  A.resize(N+1+GHOST_POINTS,N+1+GHOST_POINTS);
  B_x.resize(N+1+GHOST_POINTS);
  B_y.resize(N+1+GHOST_POINTS);
  B_z.resize(N+1+GHOST_POINTS);

  x0 = x0in;

 
  C.resize(N+2);
  // only need extra 2 ghost points because
  // tension is only second order PDE

  tension.setZero();


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
  coeffs.push_back(T(POSITION_OFFSET-2,POSITION_OFFSET-2,kappa));
  coeffs.push_back(T(POSITION_OFFSET-2,POSITION_OFFSET-1,
		     -2*kappa-Delta_s*Delta_s*tension(0)));
  coeffs.push_back(T(POSITION_OFFSET-2,POSITION_OFFSET+1,2*kappa+Delta_s*Delta_s*tension(0)));
  coeffs.push_back(T(POSITION_OFFSET-2,POSITION_OFFSET+2,-kappa));
  coeffs.push_back(T(POSITION_OFFSET-1,POSITION_OFFSET-1,kappa));
  coeffs.push_back(T(POSITION_OFFSET-1,POSITION_OFFSET,-2*kappa));
  coeffs.push_back(T(POSITION_OFFSET-1,POSITION_OFFSET+1,kappa));
  coeffs.push_back(T(POSITION_OFFSET+N,POSITION_OFFSET+N,1));
  coeffs.push_back(T(POSITION_OFFSET+N+1,POSITION_OFFSET+N-1,kappa));
  coeffs.push_back(T(POSITION_OFFSET+N+1,POSITION_OFFSET+N,-2*kappa));
  coeffs.push_back(T(POSITION_OFFSET+N+1,POSITION_OFFSET+N+1,kappa));

  int k = POSITION_OFFSET;

  coeffs.push_back(T(k,k-2,kappa));
  coeffs.push_back(T(k,k-1,set_gamma_0_minus()));
  coeffs.push_back(T(k,k,set_d_i(0)));
  coeffs.push_back(T(k,k+1,set_gamma_0_plus()));
  coeffs.push_back(T(k,k+2,kappa));

  
  for (int i = 1; i < N; i++) {
    k = i + POSITION_OFFSET;
    coeffs.push_back(T(k,k-2,kappa));
    coeffs.push_back(T(k,k-1,set_s_i_minus(i)));
    coeffs.push_back(T(k,k,set_d_i(i)));
    coeffs.push_back(T(k,k+1,set_s_i_plus(i)));
    coeffs.push_back(T(k,k+2,kappa));
  }

  return coeffs;
	
}

/* -------------------------------------------------------------------------- */
/* Update matrix A. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void Polymer::update_A()
{
  // update block matrix with rows spanning from 0 to 3+POSITION_OFFSET
  //   and columns spanning from 0 to 2+POSITION_OFFSET using slow method
  
  A.coeffRef(POSITION_OFFSET-2,POSITION_OFFSET-1)
    = -2*kappa-Delta_s*Delta_s*tension(0);
  A.coeffRef(POSITION_OFFSET-2,POSITION_OFFSET+1)
    = 2*kappa+Delta_s*Delta_s*tension(0);

  int k=POSITION_OFFSET;      // k is row number

  A.coeffRef(k,k-1) = set_gamma_0_minus();
  A.coeffRef(k,k) = set_d_i(0);
  A.coeffRef(k,k+1) = set_gamma_0_plus();
  
  for (int i = 1; i < 4; i++) { //  row = i +POSIITON_OFFSET

    k = i+POSITION_OFFSET;

    A.coeffRef(k,k-1) = set_s_i_minus(i);
    if (i < 3) 
      A.coeffRef(k,k) = set_d_i(i);
    if (i < 2) 
      A.coeffRef(k,k+1) = set_s_i_plus(i);
  }


  // update block matrix with rows spanning from POSITION_OFFSET+N-2 to
  //  POSITION_OFFSET+N-1 and columns spanning from POSITION_OFFSET+N-1
  //  to POSITION_OFFSET+N+1

  for (int i = N-2; i < N; i++) {
    k = i + POSITION_OFFSET;

    A.coeffRef(k,k+1) = set_s_i_plus(i);
    
    if (i > N- 2)
      A.coeffRef(k,k) = set_d_i(i);

  }
  

  // update middle columns
  for (int i = 3; i < N-1; i++) { // columns = i + POSITION_OFFSET
    k = i + POSITION_OFFSET;
    int count = 0;
    for (SpMat::InnerIterator it(A,k); it; ++it) {


      if (count == 1) {

      	it.valueRef() = set_s_i_plus(i-1);
      } else if (count == 2) {
      	it.valueRef() = set_d_i(i);
      } else if (count == 3) {
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
void Polymer::update_B(Eigen::Vector3d& F)
{

  double ds3 = Delta_s*Delta_s*Delta_s;
  B_x = alpha*R_x;
  B_x(0) = 2*ds3*F(0);
  B_x(POSITION_OFFSET-1) = 0.0;
  B_x(POSITION_OFFSET+N) = x0(0);
  B_x(POSITION_OFFSET+N+1) = 0.0;


  B_y = alpha*R_y;
  B_y(0) = 2*ds3*F(1);
  B_y(POSITION_OFFSET-1) = 0.0;
  B_y(POSITION_OFFSET+N) = x0(1);
  B_y(POSITION_OFFSET+N+1) = 0.0;

  B_z = alpha*R_z;
  B_z(0) = 2*ds3*F(2);
  B_z(POSITION_OFFSET-1) = 0.0;
  B_z(POSITION_OFFSET+N) = x0(2);
  B_z(POSITION_OFFSET+N+1) = 0.0;

  return;
  
  
}




/* -------------------------------------------------------------------------- */
/* Update the vectors C in M * tension = C.
/* -------------------------------------------------------------------------- */
void Polymer::update_C(Eigen::Vector3d& F)
{



  double d1x = FDiff_first(POSITION_OFFSET,R_x);
  double d1y = FDiff_first(POSITION_OFFSET,R_y);
  double d1z = FDiff_first(POSITION_OFFSET,R_z);

  double ddd3x = FDiff_third(POSITION_OFFSET,R_x);
  double ddd3y = FDiff_third(POSITION_OFFSET,R_y);
  double ddd3z = FDiff_third(POSITION_OFFSET,R_z);

  C(0) = (kappa*(d1x*ddd3x + d1y*ddd3y + d1z*ddd3z)
	  + Delta_s*Delta_s*Delta_s*(d1x*F(0) + d1y*F(1) + d1z*F(2)));

  for (int i = 1; i< N; i++) {
    C(i) = set_c_i(i+POSITION_OFFSET);
  }

  C(N) = set_c_i_last(N+POSITION_OFFSET);
  
  d1x = FDiff_first(N+POSITION_OFFSET,R_x);
  d1y = FDiff_first(N+POSITION_OFFSET,R_y);
  d1z = FDiff_first(N+POSITION_OFFSET,R_z);

  double dddd4x = BackwardDiff_fourth(N+POSITION_OFFSET,R_x);
  double dddd4y = BackwardDiff_fourth(N+POSITION_OFFSET,R_y);
  double dddd4z = BackwardDiff_fourth(N+POSITION_OFFSET,R_z);
  
  C(N+1) = kappa*(d1x*dddd4x + d1y*dddd4y + d1z*dddd4z);

  return;
  
}


double Polymer::set_c_i_last(int i)
{
  
  double dd2x = FDiff_second(i,R_x);
  double dd2y = FDiff_second(i,R_y);
  double dd2z = FDiff_second(i,R_z);

  double ddd3x = BackwardDiff_third(i,R_x);
  double ddd3y = BackwardDiff_third(i,R_y);
  double ddd3z = BackwardDiff_third(i,R_z);

  double dddd4x = BackwardDiff_fourth(i,R_x);
  double dddd4y = BackwardDiff_fourth(i,R_y);
  double dddd4z = BackwardDiff_fourth(i,R_z);

  return -kappa*(4*(dd2x*dddd4x + dd2y*dddd4y + dd2z*dddd4z)
		 + 3*(ddd3x*ddd3x+ddd3y*ddd3y+ddd3z*ddd3z));
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
double Polymer::set_gamma_0_plus()
{
  return Base::set_gamma_0_plus(kappa);
}

/* -------------------------------------------------------------------------- */
double Polymer::set_gamma_0_minus()
{
  return Base::set_gamma_0_minus(kappa);
}



/* -------------------------------------------------------------------------- */
double Polymer::set_s_i_plus(int i)
{
  return Base::set_s_i_plus(kappa,i);
}

/* -------------------------------------------------------------------------- */
double Polymer::set_s_i_minus(int i)
{
  return Base::set_s_i_minus(kappa,i);
}


/* -------------------------------------------------------------------------- */
double Polymer::set_d_i(int i)
{
  return Base::set_d_i(kappa,i);
}

