#ifndef PSPDE_CONJUGATE_VOLFRAC_NOISE_HPP
#define PSPDE_CONJUGATE_VOLFRAC_NOISE_HPP

#include <random>

#include "conjugate.hpp"

namespace psPDE {


class ConjugateVolFracNoise : public Conjugate {
public:
  ConjugateVolFracNoise(Domain &,Grid &);
  
private:


  std::mt19937 gen;

  std::uniform_real_distribution<double> real_dist;

  double normalization;

  double mobility,temp,volFH,gamma;

  bool seed_flag;
  
  
  virtual void complex_update(int,int,int) override;
  virtual void real_update(int,int,int) override;
  virtual void origin_update() override;

};

}

#endif
