#include <stdexcept>
#include <cmath>

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


void ConjugateVolFracNoise::readCoeffs(const std::vector<std::string> &v_line)
{


  mobility = temp = 1;

  int iarg = 0;

  while (iarg < v_line.size()) {

    if (v_line[iarg] == "mobility") {
      mobility = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else if (v_line[iarg] == "temp") {
      temp = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else if (v_line[iarg] == "seed") {
      seed = std::stod(v_line[iarg+1]);
      iarg += 2;
      seed_flag = true;
      
    } else {
      throw std::runtime_error("Error: invalid conjugate/volfrac command");
    }
  }

  if (!seed_flag)
    throw std::runtime_error("Error: need seed specified in conjugate/volfrac command");

  gen.seed(seed);

}

void ConjugateVolFracNoise::reset_dt(double timestep)
{
  dt = timestep;
  double invLcubed = 1.0/(domain.period[0]*domain.period[1]*domain.period[2]);


  complexprefactor = sqrt(12*temp*mobility*invLcubed);
  realprefactor = sqrt(24*temp*mobility*invLcubed);

  sqrtdt = sqrt(dt);

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

