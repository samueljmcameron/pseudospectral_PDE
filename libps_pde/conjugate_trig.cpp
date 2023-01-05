

#include "conjugate_trig.hpp"

#include <cmath>


using namespace psPDE;

ConjugateTrig::ConjugateTrig(const MPI_Comm comm,const int mpi_size, const int id,
					     fftw_MPI_3Darray<std::complex<double>> * ft_tomod)
  : Conjugate(comm,mpi_size,id,ft_tomod)
{

  complexprefactor = 1.0;
  realprefactor = sqrt(2)*complexprefactor;
  sqrtdt = 1.0;


  int Nx = ft_tomod->boxgrid[0];
  int Ny = ft_tomod->boxgrid[1];
  int Nz = ft_tomod->boxgrid[2];
  
  double dqx = ft_tomod->dx;
  double dqy = ft_tomod->dy;
  double dqz = ft_tomod->dz;

  

  double dx = 2*M_PI/dqx/Nx;
  double dy = 2*M_PI/dqy/Nx;
  double dz = 2*M_PI/dqz/Nx;

  x0 = dx*4;
  y0 = -2*dy;
  z0 = 13*dz;

}


void ConjugateTrig::complex_update(int i , int j, int k)
{


  qz = qzs[i];
  qy = qys[j];
  qx = ft_array->dx*k;

  
  (*ft_array)(i,j,k).real(complexprefactor*cos(qx*x0+qy*y0+qz*z0));// *real_dist(gen));
  (*ft_array)(i,j,k).imag(complexprefactor*sin(qx*x0+qy*y0+qz*z0));// *real_dist(gen));

  return;
  
}


void ConjugateTrig::real_update(int i, int j, int k)
{

  qz = qzs[i];
  qy = qys[j];
  qx = ft_array->dx*k;

  (*ft_array)(i,j,k) = realprefactor*cos(qx*x0+qy*y0+qz*z0);


  return;
}

void ConjugateTrig::origin_update()
{


  (*ft_array)(0,0,0) = realprefactor;


  return;
}

