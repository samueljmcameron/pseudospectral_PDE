#ifndef PSPDE_CONJUGATE_VOLFRAC_NOISE_HPP
#define PSPDE_CONJUGATE_VOLFRAC_NOISE_HPP

#include "conjugate.hpp"

namespace psPDE {


class ConjugateVolFracNoise : public Conjugate {
public:
  ConjugateVolFracNoise(const MPI_Comm,const int , const int,
			fftw_MPI_3Darray<std::complex<double>> *);
  
private:

  virtual void complex_update(int,int,int) override;
  virtual void real_update(int,int,int) override;
  virtual void origin_update() override;

};

}

#endif