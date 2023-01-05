#ifndef PSPDE_CONJUGATE_TRIG_HPP
#define PSPDE_CONJUGATE_TRIG_HPP

#include "conjugate.hpp"

namespace psPDE {


class ConjugateTrig : public Conjugate {
public:
  ConjugateTrig(const MPI_Comm,const int , const int,
		    fftw_MPI_3Darray<std::complex<double>> *);

  double x0,y0,z0;
private:



  virtual void complex_update(int,int,int) override;
  virtual void real_update(int,int,int) override;
  virtual void origin_update() override;

};

}

#endif
