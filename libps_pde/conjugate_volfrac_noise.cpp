

#include "conjugate_volfrac_noise.hpp"

#include <cmath>


using namespace psPDE;

ConjugateVolFracNoise::ConjugateVolFracNoise(const MPI_Comm comm,const int mpi_size, const int id,
					     fftw_MPI_3Darray<std::complex<double>> * ft_tomod)
  : Conjugate(comm,mpi_size,id,ft_tomod)
{

  complexprefactor = 1.0;
  realprefactor = sqrt(2)*complexprefactor;
  sqrtdt = 1.0;
}


void ConjugateVolFracNoise::complex_update(int i , int j, int k)
{

  
  qz = qzs[i];
  qy = qys[j];
  qx = ft_array->dx*k;

  q2 = qx*qx + qy*qy + qz*qz;
  (*ft_array)(i,j,k).real(complexprefactor*sqrt(q2)*sqrtdt);// *real_dist(gen));
  (*ft_array)(i,j,k).imag(complexprefactor*sqrt(q2)*sqrtdt);// *real_dist(gen));

  return;
  
}


void ConjugateVolFracNoise::real_update(int i, int j, int k)
{

  qz = qzs[i];
  qy = qys[j];
  qx = ft_array->dx*k;

  q2 = qx*qx + qy*qy + qz*qz;
  (*ft_array)(i,j,k) = realprefactor*sqrt(q2)*sqrtdt;// *real_dist(gen);

  return;
}

void ConjugateVolFracNoise::origin_update()
{

  (*ft_array)(0,0,0) = 0.0;

  return;
}

