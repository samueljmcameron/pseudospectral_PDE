

#include "conjugate_trig.hpp"

#include <cmath>


using namespace psPDE;

ConjugateTrig::ConjugateTrig(Grid *grid)
  : Conjugate(grid)
{

  if (!grid->ft_phi)
    throw std::runtime_error("Calling conjugate/trig on grid that "
			     "doesn't have concentration grid_style.");
  
  ft_array = grid->ft_phi.get();
  
  complexprefactor = 1.0;
  realprefactor = sqrt(2)*complexprefactor;
  sqrtdt = 1.0;


  int Nx = grid->ft_boxgrid[0];
  int Ny = grid->ft_boxgrid[1];
  int Nz = grid->ft_boxgrid[2];
  
  double dqx = grid->dqx();
  double dqy = grid->dqy();
  double dqz = grid->dqz();

  

  double dx = 2*M_PI/dqx/Nx;
  double dy = 2*M_PI/dqy/Ny;
  double dz = 2*M_PI/dqz/Nz;

  x0 = dx*4;
  y0 = -2*dy;
  z0 = 13*dz;

  setup();
}


void ConjugateTrig::complex_update(int i , int j, int k)
{


  double qx,qy,qz;
  qz = qzs[i];
  qy = qys[j];
  qx = grid->dqx()*k;

  
  (*ft_array)(i,j,k).real(complexprefactor*cos(qx*x0+qy*y0+qz*z0));// *real_dist(gen));
  (*ft_array)(i,j,k).imag(complexprefactor*sin(qx*x0+qy*y0+qz*z0));// *real_dist(gen));

  return;
  
}


void ConjugateTrig::real_update(int i, int j, int k)
{

  double qx,qy,qz;
  qz = qzs[i];
  qy = qys[j];
  qx = grid->dqx()*k;

  (*ft_array)(i,j,k) = realprefactor*cos(qx*x0+qy*y0+qz*z0);


  return;
}

void ConjugateTrig::origin_update()
{


  (*ft_array)(0,0,0) = realprefactor;


  return;
}

