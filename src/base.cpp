#include "base.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>
#include <cmath>
#include <random>
#include <stdexcept>


/* -------------------------------------------------------------------------- */
/* Constructor */
/* -------------------------------------------------------------------------- */
Base::Base(int Nin, double Lin, double Delta_tin,
				   Eigen::Vector3d x0in)
  : N(Nin), L(Lin), Delta_s(Lin/Nin), Delta_t(Delta_tin)
{

  alpha = Delta_s*Delta_s*Delta_s*Delta_s/Delta_t;

  x0 = x0in;
  M.resize(N+2,N+2);
  tension.resize(N+2);
  

}

/* -------------------------------------------------------------------------- */
/* Destructor */
/* -------------------------------------------------------------------------- */
Base::~Base()
{
}


/* -------------------------------------------------------------------------- */
/* Setting the internal points as 2d curve a la Seifert et al (PRL, 1996). */
/* -------------------------------------------------------------------------- */
void Base::init_R(double kappa,Eigen::Vector3d& F,int position_offset)
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
  std::cout << R_x(N+1) << std::endl;  
  for (int i = 1; i< N+1; i++) {
    R_y(i+position_offset) = (a_qs.array() *Eigen::sin((Delta_s*i*qs).array())).sum();
  }

  double s0 = 0.0;
  double tmp;
  for (int i = 1; i < N; i++) {
    tmp = FDiff_first(i+position_offset,R_y);
    s0 += sqrt(1-tmp*tmp/Delta_s/Delta_s);

    R_x(i+position_offset) = s0*Delta_s;

  }

  if (!R_x.allFinite()) {
    throw std::invalid_argument("Nan or inf in R_x when setting initial condition.");
  }
  tmp = BackwardDiff_first(N+position_offset,R_y);
  s0 += sqrt(1-tmp*tmp/Delta_s/Delta_s);
  R_x(N+position_offset) = s0*Delta_s;

  R_x.array() -= R_x(N+position_offset) - x0(0);
  R_y.array() -= R_y(N+position_offset) - x0(1);
  R_z.array() -= R_z(N+position_offset) - x0(2);
  
  return;

}





/* -------------------------------------------------------------------------- */
/* Initialise matrix M in M * tension = C. Only call at start of simulation. */
/* -------------------------------------------------------------------------- */
void Base::set_M()
{

  std::vector<T> coefficients = init_M_coeffsmatrix();

  M.setFromTriplets(coefficients.begin(),coefficients.end());

  return;

}

/* -------------------------------------------------------------------------- */
/* Helper function to initialise matrix M. */
/* -------------------------------------------------------------------------- */
std::vector<T> Base::init_M_coeffsmatrix()
{
  std::vector<T> coeffs;

  double ds4 = Delta_s*Delta_s*Delta_s*Delta_s;

  // set boundary conditions
  coeffs.push_back(T(0,0,ds4));

  coeffs.push_back(T(N+1,N-1,-ds4));
  coeffs.push_back(T(N+1,N+1,ds4));

  

  
  for (int i = 1; i < N+1; i++) {

    coeffs.push_back(T(i,i-1,ds4));
    coeffs.push_back(T(i,i,set_q_i(i)));
    coeffs.push_back(T(i,i+1,ds4));

  }

  return coeffs;

  
}


/* -------------------------------------------------------------------------- */
/* Update matrix M. Call at every time step aside from initial step. */
/* -------------------------------------------------------------------------- */
void Base::update_M()
{
  // update first two cols by hand

  M.coeffRef(1,1) = set_q_i(1);

  // update middle columns
  for (int k = 2; k < N+1; k++) {
    int count = 0;
    for (SpMat::InnerIterator it(M,k); it; ++it) {
      if (count == 1) {
      	it.valueRef() = set_q_i(k);
      }
      count += 1;
    }
  }
}



/* -------------------------------------------------------------------------- */
double Base::set_q_i(int i)
{
  double drx = FDiff_second(i,R_x);
  double dry = FDiff_second(i,R_y);
  double drz = FDiff_second(i,R_z);
  return -Delta_s*Delta_s*(drx*drx + dry*dry + drz*drz
			   + 2*Delta_s*Delta_s);
}




/* -------------------------------------------------------------------------- */
double Base::set_gamma_0_plus(double kappa)
{
  return -4*kappa-Delta_s*Delta_s/4*(4*tension(0) - tension(2)
				     +4*tension(1) - 3*tension(0));
}

/* -------------------------------------------------------------------------- */
double Base::set_gamma_0_minus(double kappa)
{
  return -4*kappa-Delta_s*Delta_s/4*(4*tension(0) + tension(2)
				     - 4*tension(1) + 3*tension(0));
}



/* -------------------------------------------------------------------------- */
double Base::set_s_i_plus(double kappa,int i)
{
  return -4*kappa-Delta_s*Delta_s/4*(4*tension(i) + tension(i+1)
				     - tension(i-1));
}

/* -------------------------------------------------------------------------- */
double Base::set_s_i_minus(double kappa,int i)
{
  return -4*kappa-Delta_s*Delta_s/4*(4*tension(i) - tension(i+1)
				     + tension(i-1));
}

/* -------------------------------------------------------------------------- */
double Base::set_d_i(double kappa,int i)
{
  return alpha + 6*kappa +  2*Delta_s*Delta_s*tension(i);
}





/* -------------------------------------------------------------------------- */
/* First order centred finite difference */
/* -------------------------------------------------------------------------- */
double Base::FDiff_first(int i,Eigen::VectorXd& vec)
{
  return (vec(i+1)-vec(i-1))/2.0;
  
}
/* -------------------------------------------------------------------------- */
/* Second order centred finite difference */
/* -------------------------------------------------------------------------- */
double Base::FDiff_second(int i,Eigen::VectorXd& vec)
{
  return (vec(i+1)+vec(i-1)-2*vec(i));
}
/* -------------------------------------------------------------------------- */
/* Third order centred finite difference */
/* -------------------------------------------------------------------------- */
double Base::FDiff_third(int i,Eigen::VectorXd& vec)
{
  return (vec(i+2)-vec(i-2)-2*vec(i+1)+2*vec(i-1))/2.0;
}

/* -------------------------------------------------------------------------- */
/* Fourth order centred finite difference */
/* -------------------------------------------------------------------------- */
double Base::FDiff_fourth(int i,Eigen::VectorXd& vec)
{
  return vec(i+2)+vec(i-2)-4*vec(i+1)-4*vec(i-1)+6*vec(i);
}

/* -------------------------------------------------------------------------- */
/* First order backward finite difference */
/* -------------------------------------------------------------------------- */
double Base::BackwardDiff_first(int i, Eigen::VectorXd& vec)
{
  return (vec(i-2) - 4 * vec(i-1) + 3* vec(i))/2.0;
}

/* -------------------------------------------------------------------------- */
/* Third order semi-backward finite difference */
/* -------------------------------------------------------------------------- */
double Base::BackwardDiff_third(int i,Eigen::VectorXd& vec)
{
  return (3*vec(i+1)-10*vec(i)+12*vec(i-1)-6*vec(i-2)+vec(i-3))/2.0;
}



/* -------------------------------------------------------------------------- */
/* Fourth order semi-backward finite difference */
/* -------------------------------------------------------------------------- */
double Base::BackwardDiff_fourth(int i,Eigen::VectorXd& vec)
{
  return 2*vec(i+1)-9*vec(i)+16*vec(i-1)-14*vec(i-2)+6*vec(i-3)-vec(i-4);
}
