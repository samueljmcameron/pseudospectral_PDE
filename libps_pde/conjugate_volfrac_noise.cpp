#include <stdexcept>
#include <cmath>

#include "randompll.hpp"
#include "conjugate_volfrac_noise.hpp"




using namespace psPDE;

ConjugateVolFracNoise::ConjugateVolFracNoise(Domain & domain, Grid &grid)
  : Conjugate(domain,grid)
{

  if (!grid.ft_noise)
    throw std::runtime_error("Calling conjugate/volfrac/noise on grid that "
			     "doesn't have noise grid_style.");
  
  ft_array = grid.ft_noise.get();
  
  complexprefactor = 1.0;
  realprefactor = sqrt(2)*complexprefactor;
  sqrtdt = 1.0;

  setup();
}


void ConjugateVolFracNoise::complex_update(int i , int j, int k)
{


  double qx,qy,qz,q2;
  qz = qzs[i];
  qy = qys[j];
  qx = domain.dqx()*k;

  q2 = qx*qx + qy*qy + qz*qz;
  (*ft_array)(i,j,k).real(complexprefactor*sqrt(q2)*sqrtdt);// *real_dist(gen));
  (*ft_array)(i,j,k).imag(complexprefactor*sqrt(q2)*sqrtdt);// *real_dist(gen));

  return;
  
}


void ConjugateVolFracNoise::real_update(int i, int j, int k)
{

  double qx,qy,qz,q2;
  qz = qzs[i];
  qy = qys[j];
  qx = domain.dqx()*k;

  q2 = qx*qx + qy*qy + qz*qz;
  (*ft_array)(i,j,k) = realprefactor*sqrt(q2)*sqrtdt;// *real_dist(gen);

  return;
}

void ConjugateVolFracNoise::origin_update()
{

  (*ft_array)(0,0,0) = 0.0;

  return;
}

