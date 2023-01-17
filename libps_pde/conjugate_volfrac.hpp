#ifndef PSPDE_CONJUGATE_VOLFRAC_HPP
#define PSPDE_CONJUGATE_VOLFRAC_HPP

#include <random>

#include "conjugate.hpp"

#include <complex>
#include <vector>
#include <string>

namespace psPDE {


class ConjugateVolFrac : public Conjugate {
public:
  ConjugateVolFrac(Domain &, Grid &);

  virtual void readCoeffs(const std::vector<std::string> &) override;
  

  virtual void reset_dt(double) override;
private:


  std::mt19937 gen;

  std::uniform_real_distribution<double> real_dist;

  double normalization;

  double mobility,temp,volFH,gamma,dt;

  int seed;

  bool seed_flag;

  std::complex<double> noise;

  fftw_MPI_3Darray<std::complex<double>> &ft_phi;
  fftw_MPI_3Darray<std::complex<double>> &ft_nonlinear;


  
  
  virtual void complex_update(int,int,int) override;
  virtual void real_update(int,int,int) override;
  virtual void origin_update() override;

};

}

#endif
