#include "integrator.hpp"


using namespace std::complex_literals;

Integrator::Integrator(MPI_Comm comm,const GridData& fourier, const int seed,
		       const SolutionParams& solparams,const double dt)
  : ft_phi(comm,"ft_phi",fourier), ft_nonlinear(comm,"ft_nl",fourier),
    fourier(fourier),seed(seed),gen(seed),
    real_dist(-0.5,0.5),mobility(solparams.mobility),gamma(solparams.gamma),
    temp(solparams.temp),chi(solparams.chi), volFH(solparams.volFH),
    dt(dt),normalization(sqrt(1.0/(fourier.get_Nx()*fourier.get_Ny()*fourier.get_Nz())))
{

  double tmp = fourier.get_dx()*fourier.get_dy()*fourier.get_dz()/(2*M_PI*2*M_PI*2*M_PI);
  local0start = ft_phi.get_local0start();
  complexprefactor = sqrt(12*temp*mobility*tmp);
  realprefactor = sqrt(24*temp*mobility*tmp);
  sqrtdt = sqrt(dt);

  std::cout << (fourier.get_dz()*fourier.get_dz()*fourier.get_Nz()/2*fourier.get_Nz()/2+fourier.get_dy()*fourier.get_dy()*fourier.get_Ny()/2*fourier.get_Ny()/2+fourier.get_dx()*fourier.get_dx()*fourier.get_Nx()/2*fourier.get_Nx()/2) << std::endl;
  
}

Integrator::~Integrator()
{
}

void Integrator::initialize(fftw_MPI_3Darray<double> & phi,
			  const double initial_volfrac, const double variance)
{
  
  for (int i = 0; i < phi.axis_size(0); i++) {
    for (int j = 0; j < phi.axis_size(1); j++) {
      for (int k = 0; k < phi.axis_size(2); k++) {
	phi(i,j,k) = initial_volfrac + variance*real_dist(gen);
      }
    }
  }
  return;
}


void Integrator::nonlinear(fftw_MPI_3Darray<double>& nonlinear,
			 const fftw_MPI_3Darray<double>& phi)
{
  for (int i = 0; i < nonlinear.axis_size(0); i++) {
    for (int j = 0; j < nonlinear.axis_size(1); j++) {
      for (int k = 0; k < nonlinear.axis_size(2); k++) {
	nonlinear(i,j,k)
	  = temp/volFH*(log(phi(i,j,k)/(1-phi(i,j,k)))+chi*(1-2*phi(i,j,k)));

	  //-quad*phi(i,j,k)+quartic*phi(i,j,k)*phi(i,j,k)*phi(i,j,k);
      }
    }
  }
  return;

}

void Integrator::integrate(int i, int j, int k)
{

  if (i + local0start > fourier.get_Nz()/2) {
    qz = fourier.get_dz()*(fourier.get_Nz()- i - local0start);
  } else {
    qz = fourier.get_dz()*(i + local0start);
  }

  if (j > fourier.get_Ny()/2) {
    qy = fourier.get_dy()*(fourier.get_Ny()-j);
  } else {
    qy = fourier.get_dy()*j;
  }  

  qx = fourier.get_dx()*k;
  q2 = qx*qx + qy*qy + qz*qz;
  noise = complexprefactor*sqrt(q2)*sqrtdt*(real_dist(gen)+1i*real_dist(gen));

  
  ode(ft_phi(i,j,k),ft_nonlinear(i,j,k),noise,q2);
  
  
  //  ft_phi(i,j,k) = 1.0*(i+local0start) + 1i*(1.0*j);
  return;
}


void Integrator::integrate_real(int i, int j, int k)
{

  if (i + local0start > fourier.get_Nz()/2) {
    qz = fourier.get_dz()*(fourier.get_Nz()- i - local0start);
  } else {
    qz = fourier.get_dz()*(i + local0start);
  }

  if (j > fourier.get_Ny()/2) {
    qy = fourier.get_dy()*(fourier.get_Ny()-j);
  } else {
    qy = fourier.get_dy()*j;
  }
  
  qx = fourier.get_dx()*k;
  q2 = qx*qx + qy*qy + qz*qz;
  
  noise = realprefactor*sqrt(q2)*sqrtdt*real_dist(gen);

  
  ode(ft_phi(i,j,k),ft_nonlinear(i,j,k),noise,q2);
  
  //  ft_phi(i,j,k) = (1.0*(i+local0start) + 1i*(1.0*j))/2.0;
  return;
}

void Integrator::ode(std::complex<double> & y, std::complex<double> ynl,
		   std::complex<double> rnd, double q2)
{
  y = (y-mobility*q2*dt*(ynl+gamma*q2*y))*normalization*normalization +  rnd;
  //  y = y*(1-mobility*q2*dt)*normalization*normalization + rnd;
  return;
}


/*
void Integrator::ode(std::complex<double> & y, std::complex<double> ynl,
		   std::complex<double> rnd, double q2)
{
  y = (((y-mobility*q2*ynl*dt)*normalization*normalization + rnd )/(1+mobility*gamma*q2*q2*dt));
  return;
}
*/
