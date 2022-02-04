#include "timestep.hpp"


using namespace std::complex_literals;

TimeStep::TimeStep(MPI_Comm comm,const int Nz, const int Ny, const int Nx,
		   const int Lz, const int Ly, const int Lx, const int seed,
		   const double mobility, const double gamma, const double temp,
		   const double chi, const double volFH,double dt)
  : ft_phi(comm,"ft_phi",Nz,Ny,Nx), ft_nonlinear(comm,"ft_nl",Nz,Ny,Nx),
    dqz(2*M_PI/Lz),dqy(2*M_PI/Ly),dqx(2*M_PI/Lx),seed(seed),gen(seed),
    real_dist(-0.5,0.5),mobility(mobility),gamma(gamma),temp(temp),chi(chi),
    volFH(volFH),dt(dt),normalization(sqrt(1.0/(Nx*Ny*Nz))),Nz(Nz),Ny(Ny),Nx(Nx)
{

  local0start = ft_phi.get_local0start();
  complexprefactor = sqrt(12*temp*mobility/(Lx*Ly*Lz));
  realprefactor = sqrt(24*temp*mobility/(Lx*Ly*Lz));
  sqrtdt = sqrt(dt);

  std::cout << (dqz*dqz*Nz/2*Nz/2+dqy*dqy*Ny/2*Ny/2+dqx*dqx*Nx/2*Nx/2) << std::endl;
  
}

TimeStep::~TimeStep()
{
}

void TimeStep::initialize(fftw_MPI_3Darray<double> & phi,
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


void TimeStep::nonlinear(fftw_MPI_3Darray<double>& nonlinear,
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

void TimeStep::update(int i, int j, int k)
{

  if (i + local0start > Nz/2) {
    qz = dqz*(Nz- i - local0start);
  } else {
    qz = dqz*(i + local0start);
  }

  if (j > Ny/2) {
    qy = dqy*(Ny-j);
  } else {
    qy = dqy*j;
  }  

  qx = dqx*k;
  q2 = qx*qx + qy*qy + qz*qz;
  noise = complexprefactor*sqrt(q2)*sqrtdt*(real_dist(gen)+1i*real_dist(gen));

  
  ode(ft_phi(i,j,k),ft_nonlinear(i,j,k),noise,q2);
  
  
  //  ft_phi(i,j,k) = 1.0*(i+local0start) + 1i*(1.0*j);
  return;
}


void TimeStep::update_real(int i, int j, int k)
{

  if (i + local0start > Nz/2) {
    qz = dqz*(Nz- i - local0start);
  } else {
    qz = dqz*(i + local0start);
  }

  if (j > Ny/2) {
    qy = dqy*(Ny-j);
  } else {
    qy = dqy*j;
  }
  
  qx = dqx*k;
  q2 = qx*qx + qy*qy + qz*qz;
  
  noise = realprefactor*sqrt(q2)*sqrtdt*real_dist(gen);

  
  ode(ft_phi(i,j,k),ft_nonlinear(i,j,k),noise,q2);
  
  //  ft_phi(i,j,k) = (1.0*(i+local0start) + 1i*(1.0*j))/2.0;
  return;
}

void TimeStep::ode(std::complex<double> & y, std::complex<double> ynl,
		   std::complex<double> rnd, double q2)
{
  y = (y-mobility*q2*dt*(ynl+gamma*q2*y))*normalization*normalization +  rnd;
  //  y = y*(1-mobility*q2*dt)*normalization*normalization + rnd;
  return;
}


/*
void TimeStep::ode(std::complex<double> & y, std::complex<double> ynl,
		   std::complex<double> rnd, double q2)
{
  y = (((y-mobility*q2*ynl*dt)*normalization*normalization + rnd )/(1+mobility*gamma*q2*q2*dt));
  return;
}
*/
